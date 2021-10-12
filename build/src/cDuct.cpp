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

  // setup cells
  for (cSICell* scell : scells) {
    scell->setup(P);
  }

  // mesh stuff (each cell has its own mesh data, lumen segments here)
  process_mesh_info();

  // setup initial conditions
  setup_IC();

  // setup arrays for ODE calculation
  setup_arrays();

  // setup the cells too
  int nscells = scells.size();
  Eigen::VectorXf cellz(nscells);
  for (int i = 0; i < nscells; i++) {
    // store centroid z coordinate for postprocessing
    cellz(i) = static_cast<float>(scells[i]->get_mean_z());
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
  resultsh5.writeAttribute(n_disc, "number of lumen discs", id);
  resultsh5.writeAttribute(CELLULARCOUNT, "number of cellular variables", id);
  resultsh5.writeAttribute(nscells, "number of cells", id);
  double outputdt = delT * Tstride;
  resultsh5.writeAttribute(outputdt, "output time interval", "/");

  // store cell centroid z components for postprocessing
  resultsh5.writeDataset(cellz, id + "/zcells");

  // store lumen segments
//  Eigen::VectorXf segf(n_disc + 1);
//  for (int i = 0; i < lumen_prop.n_int + 1; i++) {
//    segf(i) = static_cast<float>(lumen_prop.segment[i]);
//  }
//  resultsh5.writeDataset(segf, id + "/segment");

  // store t=0
  save_results();
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
  int n_l = n_disc;
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
  int n_l = n_disc;
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
  int n_l = n_disc;
  int nvars = n_sc * CELLULARCOUNT + n_l * LUMENALCOUNT;
  return nvars;
}

void cDuct::process_mesh_info() {
  out << "<Duct> process mesh info..." << std::endl;

  // segment is the sections of duct provided in the mesh file
  int n_seg = parent->ltree->segs.rows();
  out << "<Duct> number of segments = " << n_seg << std::endl;

  // segments are indexed from node 0 or duct outlet
  Array1Nd seg_length(n_seg);
  seg_length = (parent->ltree->nodes(parent->ltree->segs.col(0), Eigen::all) -
                parent->ltree->nodes(parent->ltree->segs.col(1), Eigen::all)).rowwise().norm();
  out << "<Duct> segment lengths = " << seg_length << std::endl;

  // discs are further discretisation of the duct segments
  n_disc = 0;

  // which segment the disc belongs to
  Array1Ni d_s_Vec;

  // keep track of the output segment/disc of each segment/disc, in terms of water flow
  Array1Ni seg_out_Vec = Array1Ni::Constant(n_seg, -1);

  for (int i = 0; i < n_seg; i++) {
    // number of discs in this segment
    int n = ceil(seg_length(i) / L_int);

    // resize arrays
    disc_length.conservativeResize(n_disc + n);
    d_s_Vec.conservativeResize(n_disc + n);
    disc_X_area.conservativeResize(n_disc + n);
    disc_out_Vec.conservativeResize(n_disc + n);  // TODO: should this be initialised to something...

    // the first n-1 discs have length L_int
    disc_length(Eigen::seq(n_disc, Eigen::last-1)) = L_int;

    // the last disc has the remainder as length
    if (fmod(seg_length(i), L_int) != 0) {
      disc_length(n_disc+n-1) = fmod(seg_length(i), L_int);
    } else {
      disc_length(n_disc+n-1) = L_int;
    }

    // record the duct segment the discs belongs to
    d_s_Vec(Eigen::seq(n_disc, Eigen::last)) = i;  // NOTE: we're storing zero-based index, whereas matlab stores 1-based

    // disc_mid_point is used to interpolate disc radius
    Array1Nd disc_mid_point = Array1Nd::Zero(n);
    double disc_length_sum = 0.0;
    for (int j = 0; j < n; j++) {
      disc_mid_point(j) = disc_length(j) / 2.0 + disc_length_sum;
      disc_length_sum += disc_length(j);
    }
    out << "<Duct> disc_mid_point = " << disc_mid_point << std::endl;
    double radii_this = parent->ltree->radii(i);
    double radii_next = parent->ltree->radii(i+1);
    disc_X_area(Eigen::seq(n_disc, Eigen::last)) = M_PI *
        (radii_this + (radii_next - radii_this) / seg_length(i) * disc_mid_point).pow(2);

    // seg_out(i) is the output segment of duct segment i
    // i.e. seg_out(2) = 1, the output segment of segment 2 is 1.
    // seg_out = find(segments(:,1) == segments(i,2)); 
    // find the input node that equals the output node of the current segment.
    // TODO: can the be more than one seg_out or is it always a single value
    int seg_out = -1;
    for (int j = 0; j < n_seg; j++) {
      if (parent->ltree->segs(j, 0) == parent->ltree->segs(i, 1)) {
        seg_out = j;
        break;
      }
    }
    out << "<Duct> seg_out for segment " << i << " is " << seg_out << std::endl;

    // if there is an output segment, or it is not the last segment
    if (seg_out != -1) {
      seg_out_Vec(i) = seg_out;

      // find the disc index of the corresponding output segment
      // the output disc of the first disc is the last disc of the output segment
      // (first - close to node 0, last - far from node 0)
      int seg_out_disc = -1;
      for (int j = 0; j < n_disc + n; j++) {
        if (d_s_Vec(j) == seg_out) {
          seg_out_disc = j;
        }
      }
      disc_out_Vec(n_disc) = seg_out_disc;
    } else {
      disc_out_Vec(n_disc) = -1;  // TODO: not sure if required
    }
    // the other discs in the segment are consecutive
    for (int j = 1; j < n; j++) {
      disc_out_Vec(n_disc + j) = n_disc + j - 1;
    }

    // increment disc counter
    n_disc += n;
  }
  disc_volume = disc_X_area * disc_length;

  out << "<Duct> total number of discs = " << n_disc << std::endl;
  out << "<Duct> disc lengths = " << disc_length << std::endl;
  out << "<Duct> d_s_Vec = " << d_s_Vec << std::endl;
  out << "<Duct> disc_X_area = " << disc_X_area << std::endl;
  out << "<Duct> disc_out_Vec = " << disc_out_Vec << std::endl;
  out << "<Duct> seg_out_Vec = " << seg_out_Vec << std::endl;
  out << "<Duct> disc_volume = " << disc_volume << std::endl;

  // now for the cells
  for (cSICell* scell : scells) {
    scell->process_mesh_info(seg_out_Vec, seg_length, d_s_Vec);
  }
}

void cDuct::setup_IC() {
   x_l.resize(Eigen::NoChange, n_disc);
   dxldt.resize(Eigen::NoChange, n_disc);

  for (int i = 0; i < n_disc; i++) {  // looping over lumen segments
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

  PSflow = 100.0 / 10.0;  // um3/s volumetric primary saliva flow rate

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
  P.NBC.alpha_A = utils::get_parameter_real(p, "duct_common", "NBC_alpha_A", out);
  P.NBC.alpha_B = utils::get_parameter_real(p, "duct_common", "NBC_alpha_B", out);
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

/*	
  // call the solver
  gather_x(x);
  solver->run(t, t + dt, x);
  solver->PrintFinalStatsBrief();
*/

  // store results
  stepnum++;
  if (stepnum % Tstride == 0) {
    save_results();
  }
}

void cDuct::setup_arrays() {
  // setting up arrays once at the beginning for performance
  int n_l = n_disc;

  dwAdt.resize(Eigen::NoChange, n_l);
  v.resize(Eigen::NoChange, n_l);
  v_up.resize(Eigen::NoChange, n_l);
  x_up.resize(Eigen::NoChange, n_l);
}

double cDuct::accum_fluid(const int duct_idx) {
  // need to accumulate fluid flow from all those that feed into this one
  double accum = dwAdt(duct_idx);
  for (int i = 0; i < n_disc; i++) {
    if (disc_out_Vec(i) == duct_idx) {
      accum += accum_fluid(i);
    }
  }
  return accum;
}

void cDuct::f_ODE(const Array1Nd &x_in, Array1Nd &dxdt) {
  // populate x_l and x_c from x_in
  distribute_x(x_in);

  int n_c = scells.size();
  int n_l = n_disc;

  // loop through the cells to populate the rate of change for each cell/variable
  #pragma omp parallel for
  for (int i = 0; i < n_c; i++) {
    scells[i]->f_ODE(x_l);
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
  for (int i = 0; i < n_l; i++) {
    // starting from disc i, tracing back its input disc and adding up the disc water secretion recursively
    v(i) += accum_fluid(i);
  }
  out << "DEBUG: v: " << v << std::endl;

  // construct a matrix to represent the upstream variable value for each lumen segment
  x_up.setZero();
  x_up(Eigen::all, Eigen::last) = P.ConP;
  
  // % fill up the upstream water flow v_up and variables x_up
  if (n_l > 1) {
    for (int j = 0; j < n_l; j++) {
      // find the index of the upstream disc(s)
      std::vector<int> disc_up;
      for (int k = 0; k < n_disc; k++) {
        if (disc_out_Vec(k) == j) {
          disc_up.push_back(k);
        }
      }

      if (disc_up.size() > 0) {  // there exists upstream discs
        v_up(j) = v(disc_up).sum();  // summing up upstream discs in case of a joint branch
        x_up.col(j) = x_l(Eigen::all, disc_up).rowwise().mean();  // take inflow variable as mean of upstream variable
        // !!! might need to correct to a weighted average in future versions!!!
      }
      else {  // for the most upstream discs, use primary saliva and input flow rate
        v_up(j) = P.PSflow;
        x_up(Eigen::all, j) = P.ConP;
      }
    }
  }

  // % convert volume flow rate to linear flow speed
  // v = v./A_L; % um/s 
  // v_up = v_up./A_L; % um/s
  v = v / disc_X_area.array();
  v_up = v_up / disc_X_area;

  // % 1D finite difference discretisation of the lumen, backward differences scheme
  // for i = 1:6
  //   dxldt(i,:) = dxldt(i,:) + (v_up.*x_up(i,:) - v.*x_l(i,:))./lumen_prop.disc_length;
  // end
  for (int i = 0; i < LUMENALCOUNT; i++) {
    dxldt.row(i) += (v_up * x_up.row(i) - v * x_l.row(i)) / disc_length;
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
