/*
 * cDuctSegmentStriated.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

//#include <thread>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <numeric>
#include <optional>
#include <omp.h>
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include "h5pp/h5pp.h"

#include "cCell.hpp"
#include "cDuctSegment.hpp"
#include "cDuctSegmentStriated.hpp"
#include "cCellStriated.hpp"
#include "cCVode.hpp"
#include "utils.hpp"

//#define DEBUGFODE
//#define DEBUGFODELOADX
//#define DEBUGWRITEXDOT

using namespace dss;

// the function that will be called by the SUNDIALS solver
static int ode_func(realtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  // pointer to DuctSegmentStriated object
  cDuctSegmentStriated* pt_dss = static_cast<cDuctSegmentStriated*>(user_data);

  // create input and output arrays for calling cell flow function
  int nvars = pt_dss->get_nvars();
  Array1Nd ymat(1, nvars);
  Array1Nd ydotmat(1, nvars);

  // copy input from sundials data structure to array for calling secretion function
  #pragma omp parallel for
  for (int i = 0; i < nvars; i++) { 
    ymat(i) = NV_Ith_S(y, i);
  }

  // call secretion function
  pt_dss->f_ODE(ymat, ydotmat);

  // copy result back into sundials data structure
  #pragma omp parallel for
  for (int i = 0; i < nvars; i++) {
    NV_Ith_S(ydot, i) = ydotmat(i);
  }

  return (0);
}

cDuctSegmentStriated::cDuctSegmentStriated(cMiniGlandDuct* _parent, int _seg_number) : cDuctSegment(_parent, _seg_number), stepnum(0), outputnum(0) {
  out << "<DuctSegmentStriated> initialiser" << std::endl;

  // Model input setup (TODO: should these be read from parameter file? or from another class...)

  double L_int = 1;  // um length of lumen discretisation interval
  PSflow = 100 / 10;  // um3/s volumetric primary saliva flow rate

  Conc.Int(Na) = 140.2;  // concentration of interstitium
  Conc.Int(K) = 5.3;
  Conc.Int(Cl) = 102.6;
  Conc.Int(HCO) = 24.7;
  Conc.Int(H) = 1000 * pow(10, -7.35);
  Conc.Int(CO) = 1.28;

  Conc.PS(Na) = 143.5;  // concentration of primary saliva (TODO: this to come from the Acinus??)
  Conc.PS(K) = 5.2;
  Conc.PS(Cl) = 114.5;
  Conc.PS(HCO) = 34.2;
  Conc.PS(H) = 1000 * pow(10, -7.35);
  Conc.PS(CO) = 1.28;


  // create parameters structure
  get_parameters();

  // mesh stuff (each cell has it's own mesh data, lumen segments here)
  process_mesh_info(L_int);

  // setup initial conditions
  setup_IC();

  // setup the cells too
  int ncells = cells.size();
  Eigen::VectorXf cellz(ncells);
  for (int i = 0; i < ncells; i++) {
    // have to cast to cCellStriated to get methods defined only on that class
    cCellStriated *cell_striated = static_cast<cCellStriated*>(cells[i]);

    // initialise the cell
    cell_striated->init(P);
    
    // process mesh info
    cell_striated->process_mesh_info(lumen_prop.segment);

    // store centroid z coordinate for postprocessing
    cellz(i) = static_cast<float>(cell_striated->get_mean_z());
  }

  // allocate solver vectors
  int num_var = LUMENALCOUNT*lumen_prop.n_int + CELLULARCOUNT*cells.size();
  x.resize(1, num_var);
  dxdt.resize(1, num_var);
  gather_x(x);

  // setting up the solver
  solver = new cCVode(out, p.at("odeSolverAbsTol"), p.at("odeSolverRelTol"));
  solver->init(ode_func, x, static_cast<void*>(this));

  // create hdf5 dataset
  int num_steps = std::ceil(p.at("totalT") / (p.at("delT") * p.at("Tstride"))) + 1;
  out << "<DuctSegmentStriated> output data size: " << num_steps << " x " << num_var << std::endl;
  resultsh5_filename = id + "_results.h5";
  resultsh5_dataset = id + "/x";

  // initialise the file
  h5pp::File resultsh5(resultsh5_filename, h5pp::FilePermission::REPLACE);

  // create the dataset for x (and xdot for debugging)
  Eigen::VectorXf xf(num_var);
  resultsh5.createDataset(xf, resultsh5_dataset, {num_steps, num_var});
#ifdef DEBUGWRITEXDOT
  resultsh5.createDataset(xf, id + "/xdot", {num_steps, num_var});
#endif

  // create dataset for electroneutrality
  Eigen::VectorXf ef(ncells);
  resultsh5.createDataset(ef, id + "/electroneutrality", {num_steps, ncells});

  // store some attributes (output time interval, lumen vars, etc)
  resultsh5.writeAttribute(LUMENALCOUNT, "number of lumenal variables", id);
  resultsh5.writeAttribute(lumen_prop.n_int, "number of lumen segments", id);
  resultsh5.writeAttribute(CELLULARCOUNT, "number of cellular variables", id);
  resultsh5.writeAttribute(static_cast<int>(cells.size()), "number of cells", id);
  double outputdt = p.at("delT") * p.at("Tstride");
  resultsh5.writeAttribute(outputdt, "output time interval", "/");

  // store cell centroid z components for postprocessing
  resultsh5.writeDataset(cellz, id + "/zcells");

  // store lumen segments
  Eigen::VectorXf segf(lumen_prop.n_int + 1);
  for (int i = 0; i < lumen_prop.n_int + 1; i++) {
    segf(i) = static_cast<float>(lumen_prop.segment[i]);
  }
  resultsh5.writeDataset(segf, id + "/segment");

  // store t=0
  save_results();
}

cDuctSegmentStriated::~cDuctSegmentStriated() {
  delete solver;
}

void cDuctSegmentStriated::process_mesh_info(double L) {
  out << "cDuctSegmentStriated> process mesh info..." << std::endl;
  lumen_prop.L = L;  // discretisation interval
  // min/max z coord (was read from mesh, changed to match matlab)
//  double lumen_start = std::min(vertex_out(2), vertex_in(2));
//  double lumen_end = std::max(vertex_out(2), vertex_in(2));
  double lumen_start = std::numeric_limits<double>::max();
  double lumen_end = std::numeric_limits<double>::lowest();
  int ncells = cells.size();
  for (int i = 0; i < ncells; i++) {
    // have to cast to cCellStriated to get methods defined only on that class
    cCellStriated *cell_striated = static_cast<cCellStriated*>(cells[i]);
    double cellminz = cell_striated->get_min_z();
    double cellmaxz = cell_striated->get_max_z();
    lumen_start = std::min(cellminz, lumen_start);
    lumen_end = std::max(cellmaxz, lumen_end);
  }
  double lumen_length = lumen_end - lumen_start;
  lumen_prop.n_int = ceil(lumen_length / L);
  double lumen_radius = inner_diameter / 2.0;
  lumen_prop.X_area = M_PI * pow(lumen_radius, 2);
  lumen_prop.volume = lumen_prop.X_area * L;
  out << "<cDuctSegmentStriated> lumen start, end, length: " << lumen_start << ", " << lumen_end << ", " << lumen_length << std::endl;
  out << "                       lumen radius, volume: " << lumen_radius << ", " << lumen_prop.volume << std::endl;
  out << "                       n_int: " << lumen_prop.n_int << std::endl;
  lumen_prop.segment.resize(lumen_prop.n_int + 1);
  out << "                       segment:";
  for (int i = 0; i <= lumen_prop.n_int; i++) {
    lumen_prop.segment[i] = i * L + lumen_start;
    out << "  " << lumen_prop.segment[i];
  }
  out << std::endl;
}

void cDuctSegmentStriated::get_parameters() {
  P.ConI = Conc.Int;
  P.ConP = Conc.PS;

  P.PSflow = PSflow;

  // apical channel conductances
  P.G_ENaC = p.at("G_ENaC");
  P.G_CFTR = p.at("G_CFTR");
  P.G_BK = p.at("G_BK");

  // basolateral channel conductances
  P.G_K_B = p.at("G_K_B");

  // apical or basolateral transporter rates
  P.NBC.alpha = p.at("NBC_alpha");
  P.NBC.k5_p = p.at("NBC_k5_p"); // 1/s
  P.NBC.k5_m = p.at("NBC_k5_m"); // 1/s
  P.NBC.k6_p = p.at("NBC_k6_p"); // 1/s
  P.NBC.k6_m = p.at("NBC_k6_m"); // 1/s

  P.AE2.alpha_A = p.at("AE2_alpha_A");
  P.AE2.alpha_B = p.at("AE2_alpha_B");
  P.AE2.k3_p = p.at("AE2_k3_p"); // 1/s
  P.AE2.k3_m = p.at("AE2_k3_m"); // 1/s
  P.AE2.k4_p = p.at("AE2_k4_p"); // 1/s
  P.AE2.k4_m = p.at("AE2_k4_m"); // 1/s

  P.NHE.alpha_A = p.at("NHE_alpha_A");
  P.NHE.alpha_B = p.at("NHE_alpha_B");
  P.NHE.k1_p = p.at("NHE_k1_p"); // 1/s
  P.NHE.k1_m = p.at("NHE_k1_m"); // 1/s
  P.NHE.k2_p = p.at("NHE_k2_p"); // 1/s
  P.NHE.k2_m = p.at("NHE_k2_m"); // 1/s

  // CO2 permeability
  P.p_CO = p.at("p_CO"); // 1/s

  // CO2 bicarbonate buffering
  P.buf.k_p = p.at("buf_k_p"); // /s
  P.buf.k_m = p.at("buf_k_m"); // /mMs

  // sodium potassium pump rates
  P.NKA.alpha_A = p.at("NKA_alpha_A"); // mol/m2
  P.NKA.alpha_B = p.at("NKA_alpha_B"); // mol/m2
  P.NKA.r = p.at("NKA_r"); // mM-3s-1
  P.NKA.beta = p.at("NKA_beta"); // mM-1

	// paracellular conductances
  P.G_P_Na = p.at("G_P_Na"); // S/m2
  P.G_P_K = p.at("G_P_K"); // S/m2
  P.G_P_Cl = p.at("G_P_Cl"); // S/m2

  // water permeability across membranes
  P.L_A = p.at("L_A"); //0.6e1; // um/s
  P.L_B = p.at("L_B"); // um/s

  // universal physical constants
  P.V = PSflow; // um^3 volumetric flow rate of primary saliva

  // osmolarity adjusting constants
  P.chi_C = p.at("chi_C"); // mol (40 mM * 1000 um3  = xxx e-18 mol)
  P.phi_A = p.at("phi_A"); // mM (fong 2016)
  P.phi_B = p.at("phi_B"); // mM (Mangos 1972)
}

void cDuctSegmentStriated::setup_IC() {
   x_l.resize(Eigen::NoChange, lumen_prop.n_int);
   dxldt.resize(Eigen::NoChange, lumen_prop.n_int);

  for (int i = 0; i < lumen_prop.n_int; i++) {  // looping over lumen segments
    // lumenal initial concentration
    x_l(Na, i) = 143.5;  // TODO: move these to parameter file
    x_l(K, i) = 5.2;
    x_l(Cl, i) = 114.5;
    x_l(HCO, i) = 34.2;
    x_l(H, i) = 1000 * pow(10, -7.35);
    x_l(CO, i) = 1.28;
  }
}

void cDuctSegmentStriated::distribute_x(const Array1Nd &x_in) {
  int n_c = cells.size();
  int n_l = lumen_prop.n_int;
  for (int i = 0; i < n_c; i++) {
    cCellStriated *cell_striated = static_cast<cCellStriated*>(cells[i]);
    int start = i * CELLULARCOUNT;
    cell_striated->x_c.row(0) = x_in(0, Eigen::seq(start, start+CELLULARCOUNT-1));
  }
  int s_l = CELLULARCOUNT * n_c;
  for (int i = 0; i < n_l; i++) {
    x_l.col(i) = x_in(Eigen::seq(s_l+i*LUMENALCOUNT, s_l+(i+1)*LUMENALCOUNT-1));
  }
}

void cDuctSegmentStriated::gather_x(Array1Nd &x_out) {
  int n_c = cells.size();
  int n_l = lumen_prop.n_int;
  for (int i = 0; i < n_c; i++) {
    cCellStriated *cell_striated = static_cast<cCellStriated*>(cells[i]);
    x_out(0, Eigen::seq(i*CELLULARCOUNT, (i+1)*CELLULARCOUNT-1)) = cell_striated->x_c.row(0);
  }
  int s_l = CELLULARCOUNT * n_c;
  for (int i = 0; i < n_l; i++) {
    x_out(0, Eigen::seq(s_l+i*LUMENALCOUNT, s_l+(i+1)*LUMENALCOUNT-1)) = x_l.col(i);
  }

}

void cDuctSegmentStriated::f_ODE(const Array1Nd &x_in, Array1Nd &dxdt) {
  // populate x_l and x_c from x_in
  distribute_x(x_in);

  int n_c = cells.size();
  int n_l = lumen_prop.n_int;

  // constant parameters
  double L = lumen_prop.L;
  double A_L = lumen_prop.X_area;

  // loop through the cells to populate the rate of change for each cell/variable
  #pragma omp parallel for
  for (int i = 0; i < n_c; i++) {
    static_cast<cCellStriated*>(cells[i])->f_ODE(x_l, lumen_prop);
  }

  // setup a vector to record the rate of change of lumen fluid flow
  Array1Nd dwAdt(1, n_l);
  dwAdt.setZero();

  // setup the ode rate of change matrices
  dxldt.setZero();

  // sum values from cells
  for (int i = 0; i < n_c; i++) {
    cCellStriated* cell_striated = static_cast<cCellStriated*>(cells[i]);
    dwAdt += cell_striated->dwAdt;
    dxldt += cell_striated->dxldt;
  }

  // % compute the fluid flow rate in the lumen
  // v = ones(1,n_l) * P.PSflow; % um^3/s volume flow rate of fluid out of each lumen segment
  // v_up = ones(1,n_l) * P.PSflow; % um^3/s volume flow rate of fluid into each lumen segment
  Array1Nd v(1, n_l);
  v = P.PSflow;
  Array1Nd v_up(1, n_l);
  v_up = P.PSflow;

  // accumulate the fluid as secreted from cells along the lumen
  std::partial_sum(dwAdt.begin(), dwAdt.end(), v.begin());
  v += P.PSflow;

  // construct a matrix to represent the upstream variable value for each lumen segment
  ArrayNFC x_up;
  x_up.resize(Eigen::NoChange, n_l);
  x_up.setZero();
  x_up.col(0) = P.ConP;
  
  // % fill up the upstream flow rate(v_up)/variable(x_up) with downstream flow rate/variable
  // if n_l>1 % if there are more than one lumen segment
  //   v_up(2:n_l) = v(1:n_l-1);
  //   x_up(:,2:n_l) = x_l(:,1:n_l-1);
  // end
  if (n_l > 1) {
    v_up(0, Eigen::seq(1, Eigen::last)) = v(0, Eigen::seq(0, Eigen::last - 1));
    x_up(Eigen::all, Eigen::seq(1, Eigen::last)) = x_l(Eigen::all, Eigen::seq(0, Eigen::last - 1));
  }

  // % convert volume flow rate to linear flow speed
  // v = v./A_L; % um/s 
  // v_up = v_up./A_L; % um/s
  v = v / A_L;
  v_up = v_up / A_L;

  // % 1D finite difference discretisation of the lumen, backward differences scheme
  // for i = 1:6
  //   dxldt(i,:) = dxldt(i,:) + (v_up.*x_up(i,:) - v.*x_l(i,:))./L;
  // end
  for (int i = 0; i < LUMENALCOUNT; i++) {
    dxldt.row(i) += (v_up * x_up.row(i) - v * x_l.row(i)) / L;
  }

  // % flatten the matrix to a column vector
  // dxdt = [dxcdt(:); dxldt(:)];
  for (int i = 0; i < n_c; i++) {
    cCellStriated *cell_striated = static_cast<cCellStriated*>(cells[i]);
    int idx = i * CELLULARCOUNT;
    dxdt(0, Eigen::seq(idx, idx+CELLULARCOUNT-1)) = cell_striated->dxcdt.row(0);
  }
  int s_l = CELLULARCOUNT * n_c;
  for (int i = 0; i < n_l; i++) {
    int idx = s_l + i * LUMENALCOUNT;
    dxdt(0, Eigen::seq(idx, idx+LUMENALCOUNT-1)) = dxldt.col(i);
  }
  
#ifdef DEBUGFODE
    out << std::scientific << std::setprecision(8);
    out << "================ DEBUG =================" << std::endl;
    out << "initial P.S. flow rate: %2.2f  um3 " << (v_up(0)*A_L) << std::endl;
    out << "final P.S. flow rate:   %2.2f  um3 " << (v(Eigen::last)*A_L) << std::endl;
    out << "percentage:             %2.2f  " << (v(Eigen::last)-v_up(0))/v_up(0)*100 << std::endl;
    out << "================ END DEBUG =================" << std::endl;
#endif
}

int cDuctSegmentStriated::get_nvars() {
  int n_c = cells.size();
  int n_l = lumen_prop.n_int;
  int nvars = n_c * CELLULARCOUNT + n_l * LUMENALCOUNT;
  return nvars;
}

void cDuctSegmentStriated::step(double current_time, double timestep) {
  out << "<DuctSegmentStriated> step - threads in use: " << omp_get_num_threads() << std::endl;

  // Testing: call f_ODE once
  if (stepnum == 0) {
    out << "writing x and xdot for debugging" << std::endl;
    Array1Nd testx(1, get_nvars());
    Array1Nd testxdot(1, get_nvars());
    gather_x(testx);

#ifdef DEBUGFODELOADX
    // load x from HDF5 file for debugging
    h5pp::File hxfile("testx.h5", h5pp::FilePermission::READONLY);
    Eigen::VectorXf testxf = hxfile.readDataset<Eigen::VectorXf>("x");
    testx = testxf.cast<double>();
    distribute_x(testx);
    gather_x(testx);
#endif

    f_ODE(testx, testxdot);
    std::ofstream xfile("xdump.txt");
    xfile << std::fixed << std::setprecision(15);
    xfile << testx.transpose();
    xfile.close();
    std::ofstream xdotfile("xdotdump.txt");
    xdotfile << std::fixed << std::setprecision(15);
    xdotfile << testxdot.transpose();
    xdotfile.close();

#ifdef DEBUGFODE
    utils::fatal_error("stopping for debugging", out);
#endif
  }
  // End testing

  // call the solver
  gather_x(x);
  solver->run(current_time, current_time + timestep, x);
  solver->PrintFinalStatsBrief();

  // store results
  stepnum++;
  if (stepnum % int(p.at("Tstride")) == 0) {
    save_results();
  }
}

void cDuctSegmentStriated::save_results() {
  // append to variable in HDF5 file...
  h5pp::File resultsh5(resultsh5_filename, h5pp::FilePermission::READWRITE);
  int nv = get_nvars();
  Eigen::VectorXf xf(nv);
  xf = x.cast<float>();
  resultsh5.writeHyperslab(xf, resultsh5_dataset, h5pp::Hyperslab({outputnum, 0}, {1, nv}));

#ifdef DEBUGWRITEXDOT
  // temporarily outputting xdot for debugging
  Array1Nd xdot(1, nv);
  f_ODE(x, xdot);
  Eigen::VectorXf xdotf(nv);
  xdotf = xdot.cast<float>();
  resultsh5.writeHyperslab(xdotf, id + "/xdot", h5pp::Hyperslab({outputnum, 0}, {1, nv}));
#endif

  // outputting electroneutrality check
  int nc = cells.size();
  Eigen::VectorXf ef(nc);
  for (int i = 0; i < nc; i++) {
    double cell_e = static_cast<cCellStriated*>(cells[i])->compute_electroneutrality_check();
    ef(i) = static_cast<float>(cell_e);
  }
  resultsh5.writeHyperslab(ef, id + "/electroneutrality", h5pp::Hyperslab({outputnum, 0}, {1, nc}));

  outputnum++;
}
