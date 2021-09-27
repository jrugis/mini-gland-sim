/*
 * cDuct.cpp
 *
 *  Created on: 26/07/2021
 *      Author: jrugis, cscott
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include "h5pp/h5pp.h"

#include "global_defs.hpp"
#include "utils.hpp"
#include "cMiniGland.hpp"
#include "cLTree.hpp"
#include "cSICell.hpp"
#include "cCVode.hpp"
#include "cDuct.hpp"

using namespace duct;

// the function that will be called by the SUNDIALS solver
static int ode_func(realtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  // pointer to Duct object
  cDuct* pt_duct = static_cast<cDuct*>(user_data);

  // create input and output arrays for calling cell flow function
  int nvars = pt_duct->get_nvars();
  Array1Nd ymat(1, nvars);
  Array1Nd ydotmat(1, nvars);

  // copy input from sundials data structure to array for calling secretion function
  #pragma omp parallel for
  for (int i = 0; i < nvars; i++) { 
    ymat(i) = NV_Ith_S(y, i);
  }

  // call secretion function
  pt_duct->f_ODE(ymat, ydotmat);

  // copy result back into sundials data structure
  #pragma omp parallel for
  for (int i = 0; i < nvars; i++) {
    NV_Ith_S(ydot, i) = ydotmat(i);
  }

  return (0);
}

cDuct::cDuct(cMiniGland* _parent) : parent(_parent), stepnum(0), outputnum(0)
{
  id = "_duct";    // there's only one duct object
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);
  p = parent->p; // the parameters INIReader object

  // search the mesh directory for striated and intercalated cell meshes
  for (const auto &file : std::filesystem::directory_iterator(MESH_FILE_DIR)){
	std::string fpath = std::filesystem::path(file.path());
    if(fpath.find("Cell_I") != std::string::npos) icells.push_back(new cSICell(this, fpath, INTERCALATED));
    if(fpath.find("Cell_S") != std::string::npos) scells.push_back(new cSICell(this, fpath, STRIATED));
  }
  
  //out << "<Duct> Duct segment count: " << nlsegs << std::endl; 
  out << "<Duct> Intercalated cell count: " << icells.size() << std::endl;
  out << "<Duct> Striated cell count: " << scells.size() << std::endl;

  // create parameters structure
  p = parent->p;  // pointer to ini reader object on parent
  get_parameters();

/*
  // mesh stuff (each cell has it's own mesh data, lumen segments here)
  process_mesh_info();

  // setup initial conditions
  setup_IC();

  // setup arrays for ODE calculation
  setup_arrays();

  // setup the cells too
  Eigen::VectorXf cellz(nscells);
  for (int i = 0; i < nscells; i++) {
    cSICell *cell_striated = scells[i];

    // initialise the cell
    cell_striated->init(P);
    
    // process mesh info
    cell_striated->process_mesh_info(lumen_prop.segment);

    // store centroid z coordinate for postprocessing
    cellz(i) = static_cast<float>(cell_striated->get_mean_z());
  }

  // allocate solver vectors
  int num_var = get_nvars();
  x.resize(1, num_var);
  dxdt.resize(1, num_var);
  gather_x(x);

  // setting up the solver
  double abstol = utils::get_parameter_real(p, "odeSolver", "odeSolverAbsTol", out);
  double reltol = utils::get_parameter_real(p, "odeSolver", "odeSolverRelTol", out);
  solver = new cCVode(out, abstol, reltol);
  solver->init(ode_func, x, static_cast<void*>(this));

  // simulation time parameters
  double totalT = utils::get_parameter_real(p, "time", "totalT", out);
  double delT = utils::get_parameter_real(p, "time", "delT", out);
  Tstride = utils::get_parameter_integer(p, "time", "Tstride", out);

  // create hdf5 dataset
  int num_steps = std::ceil(totalT / (delT * Tstride)) + 1;
  out << "<Duct> output data size: " << num_steps << " x " << num_var << std::endl;
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
  Eigen::VectorXf ef(nscells);
  resultsh5.createDataset(ef, id + "/electroneutrality", {num_steps, nscells});

  // store some attributes (output time interval, lumen vars, etc)
  resultsh5.writeAttribute(LUMENALCOUNT, "number of lumenal variables", id);
  resultsh5.writeAttribute(lumen_prop.n_int, "number of lumen segments", id);
  resultsh5.writeAttribute(CELLULARCOUNT, "number of cellular variables", id);
  resultsh5.writeAttribute(nscells, "number of cells", id);
  double outputdt = delT * Tstride;
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
  */
}

cDuct::~cDuct()
{
  for (unsigned int i = 0; i < icells.size(); i++) delete icells[i]; // delete the intercalated cells
  for (unsigned int i = 0; i < scells.size(); i++) delete scells[i]; // delete the striated cells
  out.close();
  //delete solver;
}

void cDuct::distribute_x(const Array1Nd &x_in) {
  int n_c = scells.size();
  int n_l = lumen_prop.n_int;
  for (int i = 0; i < n_c; i++) {
    int start = i * CELLULARCOUNT;
    scells[i]->x_c.row(0) = x_in(0, Eigen::seq(start, start+CELLULARCOUNT-1));
  }
  int s_l = CELLULARCOUNT * n_c;
  for (int i = 0; i < n_l; i++) {
    x_l.col(i) = x_in(Eigen::seq(s_l+i*LUMENALCOUNT, s_l+(i+1)*LUMENALCOUNT-1));
  }
}

void cDuct::gather_x(Array1Nd &x_out) {
  int n_c = scells.size();
  int n_l = lumen_prop.n_int;
  for (int i = 0; i < n_c; i++) {
    x_out(0, Eigen::seq(i*CELLULARCOUNT, (i+1)*CELLULARCOUNT-1)) = scells[i]->x_c.row(0);
  }
  int s_l = CELLULARCOUNT * n_c;
  for (int i = 0; i < n_l; i++) {
    x_out(0, Eigen::seq(s_l+i*LUMENALCOUNT, s_l+(i+1)*LUMENALCOUNT-1)) = x_l.col(i);
  }

}
int cDuct::get_nvars() {
  int n_sc = scells.size();
  int n_l = lumen_prop.n_int;
  int nvars = n_sc * CELLULARCOUNT + n_l * LUMENALCOUNT;
  return nvars;
}

void cDuct::process_mesh_info() {
  out << "<Duct> process mesh info..." << std::endl;
  lumen_prop.L = L_int;  // discretisation interval
  // min/max z coord (was read from mesh, changed to match matlab)
//  double lumen_start = std::min(vertex_out(2), vertex_in(2));
//  double lumen_end = std::max(vertex_out(2), vertex_in(2));
  double lumen_start = std::numeric_limits<double>::max();
  double lumen_end = std::numeric_limits<double>::lowest();
  int ncells = scells.size();
  for (int i = 0; i < ncells; i++) {
    cSICell* cell_striated = scells[i];
    double cellminz = cell_striated->get_min_z();
    double cellmaxz = cell_striated->get_max_z();
    lumen_start = std::min(cellminz, lumen_start);
    lumen_end = std::max(cellmaxz, lumen_end);
  }
  double lumen_length = lumen_end - lumen_start;
  lumen_prop.n_int = ceil(lumen_length / L_int);
  double inner_diameter = 8.0; // TODO must come from somewhere else
  out << ">>>>>>>>>> Remember to get inner_diameter from somewhere" << std::endl;
  double lumen_radius = inner_diameter / 2.0;
  lumen_prop.X_area = M_PI * pow(lumen_radius, 2);
  lumen_prop.volume = lumen_prop.X_area * L_int;
  out << "<Duct> lumen start, end, length: " << lumen_start << ", " << lumen_end << ", " << lumen_length << std::endl;
  out << "       lumen radius, volume: " << lumen_radius << ", " << lumen_prop.volume << std::endl;
  out << "       n_int: " << lumen_prop.n_int << std::endl;
  lumen_prop.segment.resize(lumen_prop.n_int + 1);
  out << "       segment:";
  for (int i = 0; i <= lumen_prop.n_int; i++) {
    lumen_prop.segment[i] = i * L_int + lumen_start;
    out << "  " << lumen_prop.segment[i];
  }
  out << std::endl;
}

void cDuct::setup_IC() {
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

void cDuct::get_parameters() {
  L_int = 1.0;  // lumen discretisation interval

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

  P.ConI = Conc.Int;
  P.ConP = Conc.PS;

  P.PSflow = PSflow;

  // load common parameters from the ini file

  // apical channel conductances
//  P.G_ENaC = utils::get_parameter_real(p, "duct_common", "G_ENaC", out);
  P.G_CFTR = utils::get_parameter_real(p, "duct_common", "G_CFTR", out);
  P.G_BK = utils::get_parameter_real(p, "duct_common", "G_BK", out);

  // basolateral channel conductances
  P.G_K_B = utils::get_parameter_real(p, "duct_common", "G_K_B", out);

  // apical or basolateral transporter rates
  P.NBC.alpha = utils::get_parameter_real(p, "duct_common", "NBC_alpha", out);
  P.NBC.k5_p = utils::get_parameter_real(p, "duct_common", "NBC_k5_p", out);
  P.NBC.k5_m = utils::get_parameter_real(p, "duct_common", "NBC_k5_m", out);
  P.NBC.k6_p = utils::get_parameter_real(p, "duct_common", "NBC_k6_p", out);
  P.NBC.k6_m = utils::get_parameter_real(p, "duct_common", "NBC_k6_m", out);

  P.AE2.alpha_A = utils::get_parameter_real(p, "duct_common", "AE2_alpha_A", out);
  P.AE2.alpha_B = utils::get_parameter_real(p, "duct_common", "AE2_alpha_B", out);
  P.AE2.k3_p = utils::get_parameter_real(p, "duct_common", "AE2_k3_p", out);
  P.AE2.k3_m = utils::get_parameter_real(p, "duct_common", "AE2_k3_m", out);
  P.AE2.k4_p = utils::get_parameter_real(p, "duct_common", "AE2_k4_p", out);
  P.AE2.k4_m = utils::get_parameter_real(p, "duct_common", "AE2_k4_m", out);

  P.NHE.alpha_A = utils::get_parameter_real(p, "duct_common", "NHE_alpha_A", out);
  P.NHE.alpha_B = utils::get_parameter_real(p, "duct_common", "NHE_alpha_B", out);
  P.NHE.k1_p = utils::get_parameter_real(p, "duct_common", "NHE_k1_p", out);
  P.NHE.k1_m = utils::get_parameter_real(p, "duct_common", "NHE_k1_m", out);
  P.NHE.k2_p = utils::get_parameter_real(p, "duct_common", "NHE_k2_p", out);
  P.NHE.k2_m = utils::get_parameter_real(p, "duct_common", "NHE_k2_m", out);

  // CO2 permeability
  P.p_CO = utils::get_parameter_real(p, "duct_common", "p_CO", out);

  // CO2 bicarbonate buffering
  P.buf.k_p = utils::get_parameter_real(p, "duct_common", "buf_k_p", out);
  P.buf.k_m = utils::get_parameter_real(p, "duct_common", "buf_k_m", out);

  // sodium potassium pump rates
//  P.NKA.alpha_A = utils::get_parameter_real(p, "duct_common", "NKA_alpha_A", out);
//  P.NKA.alpha_B = utils::get_parameter_real(p, "duct_common", "NKA_alpha_B", out);
  P.NKA.r = utils::get_parameter_real(p, "duct_common", "NKA_r", out);
  P.NKA.beta = utils::get_parameter_real(p, "duct_common", "NKA_beta", out);

	// paracellular conductances
//  P.G_P_Na = utils::get_parameter_real(p, "duct_common", "G_P_Na", out);
//  P.G_P_K = utils::get_parameter_real(p, "duct_common", "G_P_K", out);
//  P.G_P_Cl = utils::get_parameter_real(p, "duct_common", "G_P_Cl", out);

  // water permeability across membranes
//  P.L_A = utils::get_parameter_real(p, "duct_common", "L_A", out);
//  P.L_B = utils::get_parameter_real(p, "duct_common", "L_B", out);

  // universal physical constants
  P.V = PSflow; // um^3 volumetric flow rate of primary saliva

  // osmolarity adjusting constants
  P.chi_C = utils::get_parameter_real(p, "duct_common", "chi_C", out);
  P.phi_A = utils::get_parameter_real(p, "duct_common", "phi_A", out);
  P.phi_B = utils::get_parameter_real(p, "duct_common", "phi_B", out);
}

void cDuct::step(double t, double dt)
{
/*	
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
  solver->run(t, t + dt, x);
  solver->PrintFinalStatsBrief();

  // store results
  stepnum++;
  if (stepnum % Tstride == 0) {
    save_results();
  }
*/
}

void cDuct::setup_arrays() {
  // setting up arrays once at the beginning for performance
  int n_l = lumen_prop.n_int;

  dwAdt.resize(Eigen::NoChange, n_l);
  v.resize(Eigen::NoChange, n_l);
  v_up.resize(Eigen::NoChange, n_l);
  x_up.resize(Eigen::NoChange, n_l);
}

void cDuct::f_ODE(const Array1Nd &x_in, Array1Nd &dxdt) {
  // populate x_l and x_c from x_in
  distribute_x(x_in);

  int n_c = scells.size();
  int n_l = lumen_prop.n_int;

  // constant parameters
  double L = lumen_prop.L;
  double A_L = lumen_prop.X_area;

  // loop through the cells to populate the rate of change for each cell/variable
  #pragma omp parallel for
  for (int i = 0; i < n_c; i++) {
    scells[i]->f_ODE(x_l, lumen_prop);
  }

  // setup a vector to record the rate of change of lumen fluid flow
  dwAdt.setZero();

  // setup the ode rate of change matrices
  dxldt.setZero();

  // sum values from cells
  for (int i = 0; i < n_c; i++) {
    cSICell* cell_striated = scells[i];
    dwAdt += cell_striated->dwAdt;
    dxldt += cell_striated->dxldt;
  }

  // % compute the fluid flow rate in the lumen
  // v = ones(1,n_l) * P.PSflow; % um^3/s volume flow rate of fluid out of each lumen segment
  // v_up = ones(1,n_l) * P.PSflow; % um^3/s volume flow rate of fluid into each lumen segment
  v = P.PSflow;
  v_up = P.PSflow;

  // accumulate the fluid as secreted from cells along the lumen
  std::partial_sum(dwAdt.begin(), dwAdt.end(), v.begin());
  v += P.PSflow;

  // construct a matrix to represent the upstream variable value for each lumen segment
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
    int idx = i * CELLULARCOUNT;
    dxdt(0, Eigen::seq(idx, idx+CELLULARCOUNT-1)) = scells[i]->dxcdt.row(0);
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

void cDuct::save_results() {
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
  int nc = scells.size();
  Eigen::VectorXf ef(nc);
  for (int i = 0; i < nc; i++) {
    double cell_e = scells[i]->compute_electroneutrality_check();
    ef(i) = static_cast<float>(cell_e);
  }
  resultsh5.writeHyperslab(ef, id + "/electroneutrality", h5pp::Hyperslab({outputnum, 0}, {1, nc}));

  outputnum++;
}
