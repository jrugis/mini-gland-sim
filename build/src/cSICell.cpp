/*
 * cSICell.cpp
 *
 *  Created on: 26/09/2021
 *      Author: jrugis, cscott
 */

#include <cmath>
#include <iostream>
#include <string>
#include <set>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cDuct.hpp"
#include "cSIMesh.hpp"
#include "cSICell.hpp"

cSICell::cSICell(cDuct* _parent, std::string fname, cell_groups _cell_group) : parent(_parent), cell_group(_cell_group)
{
  id = "cell_" + fname.substr(fname.find("Cell")+5,4); // NOTE: one based cell id
  cell_number = std::stoi(id.substr(6,3)); 
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);
  out << std::fixed << std::setprecision(16);
  p = parent->p; // the parameters map
  mesh = new cSIMesh(fname, out);

  // calc areas of different surface regions
  api_area = 0.0;
  baslat_area = 0.0;
  napical = 0;
  mean_z = 0.0;
  for (int i = 0; i < mesh->nfaces; i++) {
    if (mesh->face_types(i) == APICAL) {
      api_area += mesh->face_areas(i);
    }
    else {  // basal or basolateral
      baslat_area += mesh->face_areas(i);
    }
    for (int j = 0; j < 3; j++) {
      mean_z += mesh->verts.row(mesh->faces(i, j))(2);
    }
  }
  mean_z = mean_z / static_cast<double>(mesh->nfaces * 3);
  out << "<SICell> num apical triangles = " << napical << std::endl;
  out << "<SICell> api_area = " << api_area << std::endl;
  out << "<SICell> baslat_area = " << baslat_area << std::endl;
  out << "<SICell> mean_z = " << mean_z << std::endl;
}

cSICell::~cSICell() {
  delete mesh;
  out.close();
}

void cSICell::setup(duct::parameters_t &parent_P) {
  // init in separate function so parent initialiser has completed
  setup_parameters(parent_P);
  setup_IC();
}

void cSICell::setup_arrays() {
  // setting up arrays once at the beginning for performance
  int n_loc_disc = loc_disc.size();
  Na_A.resize(Eigen::NoChange, n_loc_disc);
  K_A.resize(Eigen::NoChange, n_loc_disc);
  Cl_A.resize(Eigen::NoChange, n_loc_disc);
  HCO_A.resize(Eigen::NoChange, n_loc_disc);
  H_A.resize(Eigen::NoChange, n_loc_disc);
  CO_A.resize(Eigen::NoChange, n_loc_disc);
  J_A.resize(Eigen::NoChange, n_loc_disc);
  J_CDF_A.resize(Eigen::NoChange, n_loc_disc);
  J_buf_A.resize(Eigen::NoChange, n_loc_disc);
  J_NHE_A.resize(Eigen::NoChange, n_loc_disc);
  J_AE2_A.resize(Eigen::NoChange, n_loc_disc);
  J_NBC_A.resize(Eigen::NoChange, n_loc_disc);
  V_A_Cl.resize(Eigen::NoChange, n_loc_disc);
  I_CFTR.resize(Eigen::NoChange, n_loc_disc);
  V_A_HCO.resize(Eigen::NoChange, n_loc_disc);
  I_CFTR_B.resize(Eigen::NoChange, n_loc_disc);
  V_A_K.resize(Eigen::NoChange, n_loc_disc);
  I_BK.resize(Eigen::NoChange, n_loc_disc);
  V_A_Na.resize(Eigen::NoChange, n_loc_disc);
  I_ENaC.resize(Eigen::NoChange, n_loc_disc);
  J_NKA_A.resize(Eigen::NoChange, n_loc_disc);
  V_P_Na.resize(Eigen::NoChange, n_loc_disc);
  V_P_K.resize(Eigen::NoChange, n_loc_disc);
  V_P_Cl.resize(Eigen::NoChange, n_loc_disc);
  I_P_Na.resize(Eigen::NoChange, n_loc_disc);
  I_P_K.resize(Eigen::NoChange, n_loc_disc);
  I_P_Cl.resize(Eigen::NoChange, n_loc_disc);

  // constants
  w_A.resize(Eigen::NoChange, n_loc_disc);
  w_A = parent->disc_volume(loc_disc);
  A_A_disc.resize(Eigen::NoChange, n_loc_disc);
  A_A_disc = api_area_discs(loc_disc);
}

void cSICell::setup_parameters(duct::parameters_t &parent_P) {
  // copy of parents local parameters object
  P = parent_P;

  // read in specific striated parameters
  //
  // apical channel conductances
  P.G_ENaC = utils::get_parameter_real(p, "striated", "G_ENaC", out);

  // sodium potassium pump rates
  P.NKA.alpha_A = utils::get_parameter_real(p, "striated", "NKA_alpha_A", out);
  P.NKA.alpha_B = utils::get_parameter_real(p, "striated", "NKA_alpha_B", out);

	// paracellular conductances
  P.G_P_Na = utils::get_parameter_real(p, "striated", "G_P_Na", out);
  P.G_P_K = utils::get_parameter_real(p, "striated", "G_P_K", out);
  P.G_P_Cl = utils::get_parameter_real(p, "striated", "G_P_Cl", out);

  // water permeability across membranes
  P.L_A = utils::get_parameter_real(p, "striated", "L_A", out);
  P.L_B = utils::get_parameter_real(p, "striated", "L_B", out);
}

void cSICell::process_mesh_info(const Array1Ni &seg_out_Vec, const Array1Nd &seg_length, const Array1Ni &d_s_Vec) {
  // duct indices
  std::set<int> duct_seg, duct_seg_tot;
  for (int i = 0; i < mesh->nfaces; i++) {
    // the duct index covered by the cell
    duct_seg_tot.insert(mesh->duct_idx(i));

    // the duct index corresponding to the apical triangles
    if (mesh->face_types(i) == APICAL) {
      duct_seg.insert(mesh->duct_idx(i));
    }
  }
  out << "<SICell> duct_seg_tot:";
  for (int idx : duct_seg_tot) out << " " << idx;
  out << std::endl;
  out << "<SICell> duct_seg:";
  for (int idx : duct_seg) out << " " << idx;
  out << std::endl;

  // convert distance along duct segment to distance from node 0
  Array1Nd total_dist_along_duct(mesh->nfaces);

  // Since triangle dist_along_duct is segment based, we need to consider segment by segment
  for (int s : duct_seg_tot) {
    // MATLAB: tri_seg_idx = find(duct_idx == s);
    std::vector<int> tri_seg_idx;
    for (int i = 0; i < mesh->nfaces; i++) {
      if (mesh->duct_idx(i) == s) {
        tri_seg_idx.push_back(i);
      }
    }
    double dist_start_seg = calc_dist_start_seg(s, seg_out_Vec, seg_length);

		// calculate the total distance
    //   - discs are indexed from node 0 of the whole duct
    //   - in contrast, the triangle's distance along duct is measured from acinus
    //   - thus dist along duct needs to be reversed.
    // total_dist_along_duct(tri_seg_idx) = seg_length(s) - dist_along_duct(tri_seg_idx) + dist_start_seg;
    total_dist_along_duct(tri_seg_idx) = seg_length(s) - mesh->dist_along_duct.array()(tri_seg_idx) + dist_start_seg;
  }

  // the average distance along the duct for a cell
  mean_dist = total_dist_along_duct.mean();

  // array for disc apical areas
  api_area_discs.resize(parent->n_disc);
  api_area_discs.setZero();

  // For loop to put apical triangles into the corresponding discs, based on
  // their total distance along duct.
  // This needs to be done segment by segment, in case of duct branches.
  // (triangles on two branches have similar distance but belong to diff faces)
  for (int s : duct_seg) {
    out << "DEBUG duct_seg loop " << s << std::endl;
    // apical indices corresponding to segment s
    std::vector<int> tri_seg_idx;
    for (int i = 0; i < mesh->nfaces; i++) {
      if (mesh->face_types(i) == APICAL && mesh->duct_idx(i) == s) {
        tri_seg_idx.push_back(i);
      }
    }
    Array1Nd total_dist_along_duct_seg(tri_seg_idx.size());
    total_dist_along_duct_seg = total_dist_along_duct(tri_seg_idx);
    out << "  DEBUG tri_seg_idx.size() = " << tri_seg_idx.size() << std::endl;

    // find all discs in segment s
    std::vector<int> all_discs_in_seg;
    for (int i = 0; i < parent->n_disc; i++) {
      if (d_s_Vec(i) == s) {
        all_discs_in_seg.push_back(i);
      }
    }
    int distal_disc = *std::min_element(all_discs_in_seg.begin(), all_discs_in_seg.end());  // close to node 0 disc
    int proxim_disc = *std::max_element(all_discs_in_seg.begin(), all_discs_in_seg.end());  // far from node 0 disc
    out << "  DEBUG distal and proxim discs: " << distal_disc << ", " << proxim_disc << std::endl;

    // make an array of bins for apical triangle distances, using the discs
    // consider all the discs from node 0 to segment s
    double dist_start_seg = calc_dist_start_seg(s, seg_out_Vec, seg_length);
    Array1Nd disc_edges(proxim_disc + 2);  // +2 because indexing from 0
    disc_edges.setZero();
    disc_edges(distal_disc) = dist_start_seg;
    for (int k = distal_disc + 1; k < proxim_disc + 2; k++) {
      disc_edges(k) = disc_edges(k - 1) + parent->disc_length(k - 1);
    }
    disc_edges(proxim_disc + 1) += 1e-12;  // for the last bin, allow equal on upper edge too
    out << "  DEBUG disc_edges (len=" << proxim_disc+1 << "):";
    for (auto val : disc_edges) out << " " << val;
    out << std::endl;

    // NOTE: using "distal_disc - 1" due to some faces belonging below distal_disc (not sure whether this should be happening)
    distal_disc = std::max(0, distal_disc - 1);
    // drop the distance along the duct into the disc edge bins
    MatrixN1i api_disc_conn(tri_seg_idx.size());
    api_disc_conn.setConstant(-1);
    for (std::vector<int>::size_type i = 0; i < tri_seg_idx.size(); i++) {
      bool binned = false;
      for (int j = distal_disc; j < proxim_disc + 1; j++) {
        if ((disc_edges(j) <= total_dist_along_duct_seg(i)) && (total_dist_along_duct_seg(i) < disc_edges(j+1))) {
          api_disc_conn(i) = j;
          binned = true;
          break;
        }
      }
      if (not binned) {
        out << "ERROR could not discretise " << i << " of " << tri_seg_idx.size() << ": " << total_dist_along_duct_seg(i) << std::endl;
        utils::fatal_error("failed to discretise discs", out);
      }
    }
    out << "  DEBUG api_disc_conn:";
    for (auto val : api_disc_conn) out << " " << val;
    out << std::endl;

    for (std::vector<int>::size_type i = 0; i < tri_seg_idx.size(); i++) {
      // disc index for this apical triangle
      int disc_idx = api_disc_conn(i);
      int face_idx = tri_seg_idx[i];
      api_area_discs(disc_idx) += mesh->face_areas(face_idx);
    }
  }
  out << "  DEBUG api_area_discs:";
  for (auto v : api_area_discs) out << " " << v;
  out << std::endl;

  // discs with faces
  for (int i = 0; i < parent->n_disc; i++) {
    if (api_area_discs(i) > 0.0) {
      loc_disc.push_back(i);
    }
  }
  out << "  DEBUG loc_disc:";
  for (auto v : loc_disc) out << " " << v;
  out << std::endl;

  // scale the rates based on cell surface areas
  double A = 104.719755;  // the apical area used to scale the conductances G
  double A_A = api_area;
  double A_B = baslat_area;
  scaled_rates.L_A = P.L_A * A / A_A;
  scaled_rates.L_B = P.L_B * A / A_B;
  scaled_rates.G_ENaC = P.G_ENaC * A / A_A;
  scaled_rates.G_CFTR = P.G_CFTR * A / A_A;
  scaled_rates.G_BK   = P.G_BK * A / A_A;
  scaled_rates.G_K_B  = P.G_K_B * A / A_B;
  scaled_rates.G_P_Na = P.G_P_Na * A / A_A;
  scaled_rates.G_P_K  = P.G_P_K * A / A_A;
  scaled_rates.G_P_Cl = P.G_P_Cl * A / A_A;
  scaled_rates.NKA    = P.NKA;
  scaled_rates.NKA.alpha_A = P.NKA.alpha_A * A / A_A;
  scaled_rates.NKA.alpha_B = P.NKA.alpha_B * A / A_B;

  // setup ode function arrays
  dwAdt.resize(Eigen::NoChange, parent->n_disc);
  dwAdt.setZero();
  dxldt.resize(Eigen::NoChange, parent->n_disc);
  dxldt.setZero();

  // set up arrays for f_ODE calculation
  setup_arrays();
}

double cSICell::calc_dist_start_seg(const int seg_idx, const Array1Ni &seg_out_Vec, const Array1Nd &seg_length) {
  // distance_start = distance from node 0 at the start of seg_idx
  double dist_start_seg = 0.0;

  // starting from the current segment, tracing its output segment until the root
  int k = seg_idx;
  while (seg_out_Vec(k) != -1) {
    k = seg_out_Vec(k);

    // adding up the output segment lengths
    dist_start_seg += seg_length(k);
  }

  return dist_start_seg;
}

void cSICell::setup_IC(){
  // cellular initial concentration
  x_c(0) = -26.1257;  // TODO: where should these 3 be read from?
  x_c(1) = -52.2513;
  x_c(2) = 1000;
  x_c(3) = 17;  // TODO: move these to parameter file
  x_c(4) = 140;
  x_c(5) = 22;
  x_c(6) = 75;
  x_c(7) = 1000 * pow(10, -7.35);
  x_c(8) = 1.28;
}

void cSICell::f_ODE(const duct::ArrayNFC &x_l) {
  // constant parameters
  double Na_B = P.ConI(duct::Na);
  double K_B = P.ConI(duct::K);
  double Cl_B = P.ConI(duct::Cl);
  double HCO_B = P.ConI(duct::HCO);
  double H_B = P.ConI(duct::H);
  double CO_B = P.ConI(duct::CO);

  double phi_A = P.phi_A; // mol per lumen interval volume
  double phi_B = P.phi_B; // mM

  double alpha_NHE_A = P.NHE.alpha_A;
  double alpha_NHE_B = P.NHE.alpha_B;
  double k1_p = P.NHE.k1_p; // 1/s
  double k1_m = P.NHE.k1_m; // 1/s
  double k2_p = P.NHE.k2_p; // 1/s
  double k2_m = P.NHE.k2_m; // 1/s

  double alpha_AE2_A = P.AE2.alpha_A;
  double alpha_AE2_B = P.AE2.alpha_B;
  double k3_p = P.AE2.k3_p; // 1/s
  double k3_m = P.AE2.k3_m; // 1/s
  double k4_p = P.AE2.k4_p; // 1/s
  double k4_m = P.AE2.k4_m; // 1/s

  double alpha_NBC_A = P.NBC.alpha_A;
  double alpha_NBC_B = P.NBC.alpha_B;
  double k5_p = P.NBC.k5_p; // 1/s
  double k5_m = P.NBC.k5_m; // 1/s
  double k6_p = P.NBC.k6_p; // 1/s
  double k6_m = P.NBC.k6_m; // 1/s

  double r_NKA = P.NKA.r; // mM-3s-1
  double beta_NKA = P.NKA.beta; // mM-1

  double p_CO = P.p_CO; // 1/s 
  double k_buf_p = P.buf.k_p; // /s
  double k_buf_m = P.buf.k_m; // /mMs
  // setup the ODE rate of change matrix
  dxcdt.setZero();

  // cell specific parameters
  double G_ENaC = scaled_rates.G_ENaC;
  double G_CFTR = scaled_rates.G_CFTR;
  double G_BK   = scaled_rates.G_BK;
  double G_K_B  = scaled_rates.G_K_B;
  double G_P_Na = scaled_rates.G_P_Na;
  double G_P_K  = scaled_rates.G_P_K;
  double G_P_Cl = scaled_rates.G_P_Cl;

  double alpha_NKA_A = scaled_rates.NKA.alpha_A;
  double alpha_NKA_B = scaled_rates.NKA.alpha_B;

  double L_B = scaled_rates.L_B; // um/s 
  double L_A = scaled_rates.L_A; // um/s 

  // cellular variables
  double V_A = x_c(0);
  double V_B = x_c(1);
  double w_C = x_c(2);
  double Na_C = x_c(3);
  double K_C = x_c(4);
  double Cl_C = x_c(5);
  double HCO_C = x_c(6);
  double H_C = x_c(7);
  double CO_C = x_c(8);

  // lumenal variables of all the lumen discs this cell interfaces with
  double A_A = api_area;
  double A_B = baslat_area;

  Na_A = x_l(0, loc_disc);
  K_A = x_l(1, loc_disc);
  Cl_A = x_l(2, loc_disc);
  HCO_A = x_l(3, loc_disc);
  H_A = x_l(4, loc_disc);
  CO_A = x_l(5, loc_disc);

  // water transport
  // osm_c = chi_C./w_C*1e18; % osmolarity of cell due to proteins (chi)
  //
  // J_B = 1e-18*L_B.*V_w.*(Na_C + K_C + Cl_C + HCO_C + osm_c - Na_B - K_B - Cl_B - HCO_B - phi_B); % um/s 
  // J_A = 1e-18*L_A.*V_w.*(Na_A + K_A + Cl_A + HCO_A + phi_A - Na_C - K_C - Cl_C - HCO_C - osm_c); % um/s [1, n_loc_int]
  // dwdt = A_B * J_B - sum(A_A_int .* J_A); % um^3/s
  // dwAdt(1,loc_int) = dwAdt(1,loc_int) + A_A_int .* J_A; % um^3/s [1, n_loc_int]
  double osm_c = P.chi_C / w_C * 1e18;
  double J_B = 1e-18*L_B*V_w*(Na_C + K_C + Cl_C + HCO_C + osm_c - Na_B - K_B - Cl_B - HCO_B - phi_B);
  J_A = 1e-18*L_A*V_w*(Na_A + K_A + Cl_A + HCO_A + phi_A - Na_C - K_C - Cl_C - HCO_C - osm_c);
  double dwdt = A_B * J_B - (A_A_disc * J_A).sum();
  dwAdt(0, loc_disc) = A_A_disc * J_A;

  // % CDF C02 Diffusion 
  // J_CDF_A = p_CO * (CO_C - CO_A).* w_C .* A_A_int./A_A; % e-18 mol/s [1, n_loc_int]
  // J_CDF_B = p_CO * (CO_C - CO_B).* w_C; % e-18 mol/s
  J_CDF_A = p_CO * (CO_C - CO_A) * w_C * A_A_disc / A_A;
  double J_CDF_B = p_CO * (CO_C - CO_B) * w_C;

  // % buf CO2 buffering
  // J_buf_C = (k_buf_p*CO_C - k_buf_m.*HCO_C.*H_C).* w_C; % e-18 mol/s 
  // J_buf_A =  k_buf_p*CO_A - k_buf_m.*HCO_A.*H_A; %mM/s [1,n_loc_int]
  double J_buf_C = (k_buf_p * CO_C - k_buf_m * HCO_C * H_C) * w_C;
  J_buf_A = k_buf_p * CO_A - k_buf_m * HCO_A * H_A;

  // % NHE
  // J_NHE_A = alpha_NHE_A*(k1_p*k2_p*Na_A.*H_C-k1_m*k2_m*Na_C.*H_A)./(k1_p*Na_A+k2_p*H_C+k1_m*H_A+k2_m*Na_C).*w_C.* A_A_int./A_A; % e-18 mol/s [1,n_loc_int]
  // J_NHE_B = alpha_NHE_B*(k1_p*k2_p*Na_B *H_C-k1_m*k2_m*Na_C *H_B)./(k1_p*Na_B+k2_p*H_C+k1_m*H_B+k2_m*Na_C).*w_C; % e-18 mol/s 
  J_NHE_A = alpha_NHE_A * (k1_p * k2_p * Na_A * H_C - k1_m * k2_m * Na_C * H_A) / (k1_p * Na_A + k2_p * H_C + k1_m * H_A + k2_m * Na_C) * w_C * A_A_disc / A_A;
  double J_NHE_B = alpha_NHE_B * (k1_p * k2_p * Na_B * H_C - k1_m * k2_m * Na_C * H_B) / (k1_p * Na_B + k2_p * H_C + k1_m * H_B + k2_m * Na_C) * w_C;

  // % AE2
  // J_AE2_A = alpha_AE2_A*(k3_p*k4_p*Cl_A.*HCO_C - k3_m*k4_m*Cl_C.*HCO_A)./(k3_p*Cl_A+k4_p*HCO_C+k3_m*HCO_A+k4_m*Cl_C).*w_C.* A_A_int./A_A; % e-18 mol/s [1,n_loc_int]
  // J_AE2_B = alpha_AE2_B*(k3_p*k4_p*Cl_B.*HCO_C - k3_m*k4_m*Cl_C.*HCO_B)./(k3_p*Cl_B+k4_p*HCO_C+k3_m*HCO_B+k4_m*Cl_C).*w_C; % e-18 mol/s 
  J_AE2_A = alpha_AE2_A * (k3_p * k4_p * Cl_A * HCO_C - k3_m * k4_m * Cl_C * HCO_A) / (k3_p * Cl_A + k4_p * HCO_C + k3_m * HCO_A + k4_m * Cl_C) * w_C * A_A_disc / A_A;
  double J_AE2_B = alpha_AE2_B*(k3_p*k4_p*Cl_B*HCO_C - k3_m*k4_m*Cl_C*HCO_B)/(k3_p*Cl_B+k4_p*HCO_C+k3_m*HCO_B+k4_m*Cl_C)*w_C;

  // % NBC
  // J_NBC_B = alpha_NBC_B.*(k5_p*k6_p*Na_C.*HCO_C-k5_m*k6_m*Na_B*HCO_B)./(k5_p.*Na_C.*HCO_C+k6_p*k5_m+k6_m*Na_B*HCO_B).*w_C; % e-18 mol/s
  // J_NBC_A = alpha_NBC_A.*(k5_p*k6_p*Na_C.*HCO_C-k5_m*k6_m.*Na_A.*HCO_A)./(k5_p.*Na_C.*HCO_C+k6_p*k5_m+k6_m.*Na_A.*HCO_A).*w_C.* A_A_disc./A_A; % e-18 mol/s [1,n_loc_disc]
  double J_NBC_B = alpha_NBC_B * (k5_p * k6_p * Na_C * HCO_C - k5_m * k6_m * Na_B * HCO_B) / (k5_p * Na_C * HCO_C + k6_p * k5_m + k6_m * Na_B * HCO_B) * w_C;
  J_NBC_A = alpha_NBC_A * (k5_p * k6_p * Na_C * HCO_C - k5_m * k6_m * Na_A * HCO_A) / (k5_p * Na_C * HCO_C + k6_p * k5_m + k6_m * Na_A * HCO_A) * w_C * A_A_disc / A_A;

  // % CFTR Apical 
  // V_A_Cl = 1e3*R*T/(-1*F).*log(Cl_A./Cl_C); % mV [1,n_loc_int]
  // I_CFTR = G_CFTR .* A_A_int .* (V_A - V_A_Cl); % e-6 nA [1,n_loc_int]
  V_A_Cl = 1e3*R*T/(-1*F_const)*(Cl_A/Cl_C).log();
  I_CFTR = G_CFTR * A_A_disc * (V_A - V_A_Cl);

  // % CFTR_B Apical
  // V_A_HCO = 1e3*R*T/((-1)*F).*log(HCO_A./HCO_C); % mV [1,n_loc_int]
  // I_CFTR_B = 0.25 * G_CFTR .* A_A_int .* (V_A - V_A_HCO); % e-6 nA [1,n_loc_int]
  V_A_HCO = 1e3 * R * T / ((-1) * F_const) * (HCO_A / HCO_C).log();
  I_CFTR_B = 0.25 * G_CFTR * A_A_disc * (V_A - V_A_HCO);

  // % I_BK Apical
  // V_A_K = 1e3*R*T/F.*log(K_A./K_C); % mV  [1,n_loc_int]
  // I_BK = G_BK .* A_A_int .* (V_A - V_A_K); % e-6 nA 
  V_A_K = 1e3*R*T/F_const*(K_A/K_C).log();
  I_BK = G_BK * A_A_disc * (V_A - V_A_K);

  // % I_K_B Basolateral
  // V_B_K = 1e3*R*T/F.*log(K_B./K_C); % mV
  // I_K_B = G_K_B * A_B .* (V_B - V_B_K); % e-6 nA 
  double V_B_K = 1e3*R*T/F_const*std::log(K_B/K_C);
  double I_K_B = G_K_B * A_B * (V_B - V_B_K);

  // % ENaC Apical
  // V_A_Na = 1e3*R*T/F*log(Na_A./Na_C); % mV  [1,n_loc_int]
  // I_ENaC = G_ENaC .* A_A_int .* (V_A - V_A_Na); % e-6 nA
  V_A_Na = 1e3*R*T/F_const*(Na_A/Na_C).log();
  I_ENaC = G_ENaC * A_A_disc * (V_A - V_A_Na);

  // % NaKATPase, NKA 
  // J_NKA_A = A_A_int .* alpha_NKA_A * r_NKA .*(K_A.^2.*Na_C.^3)./(K_A.^2+beta_NKA*Na_C.^3); % 10^-12 mol/s [1,n_loc_int]
  // J_NKA_B = A_B .* alpha_NKA_B * r_NKA .*(K_B.^2.*Na_C.^3)./(K_B.^2+beta_NKA*Na_C.^3); % 10^-12 mol/s
  J_NKA_A = A_A_disc * alpha_NKA_A * r_NKA * (K_A.pow(2)*std::pow(Na_C,3))/(K_A.pow(2)+beta_NKA*std::pow(Na_C,3));
  double J_NKA_B = A_B * alpha_NKA_B * r_NKA *(std::pow(K_B,2)*std::pow(Na_C,3))/(std::pow(K_B,2)+beta_NKA*std::pow(Na_C,3));

  // % Paracellular currents
  // V_T      = V_A - V_B; % mV
  // V_P_Na   = 1e3*R*T/F.*log(Na_A/Na_B); % mV [1,n_loc_int]
  // V_P_K    = 1e3*R*T/F.*log(K_A/K_B); % mV [1,n_loc_int]
  // V_P_Cl   = 1e3*R*T/(-F).*log(Cl_A/Cl_B); % mV [1,n_loc_int]
  // I_P_Na   = G_P_Na .* A_A_int .* (V_T - V_P_Na); % e-6 nA [1,n_loc_int]
  // I_P_K    = G_P_K .* A_A_int .* (V_T - V_P_K); % e-6 nA [1,n_loc_int]
  // I_P_Cl   = G_P_Cl .* A_A_int .* (V_T - V_P_Cl); % e-6 nA [1,n_loc_int]
  double V_T = V_A - V_B;
  V_P_Na   = 1e3*R*T/F_const*(Na_A/Na_B).log();
  V_P_K    = 1e3*R*T/F_const*(K_A/K_B).log();
  V_P_Cl   = 1e3*R*T/(-F_const)*(Cl_A/Cl_B).log();
  I_P_Na   = G_P_Na * A_A_disc * (V_T - V_P_Na);
  I_P_K    = G_P_K * A_A_disc * (V_T - V_P_K);
  I_P_Cl   = G_P_Cl * A_A_disc * (V_T - V_P_Cl);

  // % V_A e-15 c/s
  // dxcdt(1,i) = -(sum(F*J_NKA_A*1e3 + I_ENaC + I_BK + I_CFTR + I_CFTR_B + I_P_Na + I_P_K + I_P_Cl));
  dxcdt(0) = -(F_const*J_NKA_A*1e3 + I_ENaC + I_BK + I_CFTR + I_CFTR_B + I_P_Na + I_P_K + I_P_Cl).sum();
  // % V_B e-15 c/s
  // dxcdt(2,i) = -(F*J_NKA_B*1e3 + I_K_B - sum(I_P_Na + I_P_K + I_P_Cl));
  dxcdt(1) = -(F_const*J_NKA_B*1e3 + I_K_B - (I_P_Na + I_P_K + I_P_Cl).sum());
  // % w_C um^3
  // dxcdt(3,i) = dwdt;
  dxcdt(2) = dwdt;
  // % Na_C mM/s
  // dxcdt(4,i) = -dwdt*Na_C/w_C + 1e3*(-sum(I_ENaC)./(F*w_C)) - 1e6*(3*(J_NKA_B+sum(J_NKA_A))/w_C) + sum(J_NBC_A)/w_C + J_NBC_B/w_C + sum(J_NHE_A)/w_C + J_NHE_B/w_C;
  dxcdt(3) = -dwdt*Na_C/w_C + 1e3*(-I_ENaC.sum()/(F_const*w_C)) - 1e6*(3*(J_NKA_B+J_NKA_A.sum())/w_C) + J_NBC_A.sum()/w_C + J_NBC_B/w_C + J_NHE_A.sum()/w_C + J_NHE_B/w_C;
  // % K_C mM/s
  // dxcdt(5,i) = -dwdt*K_C/w_C + 1e3*(-sum(I_BK)./(F*w_C) - I_K_B./(F*w_C)) + 1e6*(2*(J_NKA_B+sum(J_NKA_A))/w_C);
  dxcdt(4) = -dwdt*K_C/w_C + 1e3*(-I_BK.sum()/(F_const*w_C) - I_K_B/(F_const*w_C)) + 1e6*(2*(J_NKA_B+J_NKA_A.sum())/w_C);
  // % Cl_C mM/s
  // dxcdt(6,i) = -dwdt*Cl_C/w_C + 1e3*(sum(I_CFTR)./(F*w_C)) + sum(J_AE2_A)/w_C + J_AE2_B/w_C;
  dxcdt(5) = -dwdt*Cl_C/w_C + 1e3*(I_CFTR.sum()/(F_const*w_C)) + J_AE2_A.sum()/w_C + J_AE2_B/w_C;
  // % HCO_C mM/s
  // dxcdt(7,i) = -dwdt*HCO_C/w_C + 1e3*(sum(I_CFTR_B)./(F*w_C)) + sum(J_NBC_A)/w_C + J_NBC_B/w_C - sum(J_AE2_A)/w_C - J_AE2_B/w_C + J_buf_C/w_C;
  dxcdt(6) = -dwdt*HCO_C/w_C + 1e3*(I_CFTR_B.sum()/(F_const*w_C)) + J_NBC_A.sum()/w_C + J_NBC_B/w_C - J_AE2_A.sum()/w_C - J_AE2_B/w_C + J_buf_C/w_C;
  // % H_C mM/s
  // dxcdt(8,i) = -dwdt*H_C/w_C - sum(J_NHE_A)/w_C - J_NHE_B/w_C + J_buf_C/w_C;
  dxcdt(7) = -dwdt*H_C/w_C - J_NHE_A.sum()/w_C - J_NHE_B/w_C + J_buf_C/w_C;
  // % CO_C mM/s
  // dxcdt(9,i) = -dwdt*CO_C/w_C - sum(J_CDF_A)/w_C - J_CDF_B/w_C - J_buf_C/w_C;
  dxcdt(8) = -dwdt*CO_C/w_C - J_CDF_A.sum()/w_C - J_CDF_B/w_C - J_buf_C/w_C;

  // % Na_A mM/s
  // dxldt(1,loc_disc) = dxldt(1,loc_disc) + 1e6*(3*J_NKA_A./w_A) + 1e3*(I_ENaC./(F*w_A)) + 1e3*(I_P_Na./(F*w_A)) - J_NHE_A./w_A - J_NBC_A./w_A;
  dxldt(0,loc_disc) = 1e6*(3*J_NKA_A/w_A) + 1e3*(I_ENaC/(F_const*w_A)) + 1e3*(I_P_Na/(F_const*w_A)) - J_NHE_A/w_A - J_NBC_A/w_A;
  // % K_A mM/s
  // dxldt(2,loc_disc) = dxldt(2,loc_disc) - 1e6*(2*J_NKA_A./w_A) + 1e3*(I_BK./(F*w_A)) + 1e3*(I_P_K./(F*w_A));
  dxldt(1,loc_disc) = - 1e6*(2*J_NKA_A/w_A) + 1e3*(I_BK/(F_const*w_A)) + 1e3*(I_P_K/(F_const*w_A));
  // % Cl_A mM/s
  // dxldt(3,loc_disc) = dxldt(3,loc_disc) + 1e3*(-I_CFTR./(F*w_A)) + 1e3*(-I_P_Cl./(F*w_A)) - J_AE2_A./w_A;
  dxldt(2,loc_disc) = 1e3*(-I_CFTR/(F_const*w_A)) + 1e3*(-I_P_Cl/(F_const*w_A)) - J_AE2_A/w_A;
  // % HCO_A mM/s
  // dxldt(4,loc_disc) = dxldt(4,loc_disc) + 1e3*(-I_CFTR_B./(F*w_A)) - J_NBC_A./w_A + J_AE2_A./w_A + J_buf_A;
  dxldt(3,loc_disc) = 1e3*(-I_CFTR_B/(F_const*w_A)) - J_NBC_A/w_A + J_AE2_A/w_A + J_buf_A;
  // % H_A mM/s
  // dxldt(5,loc_disc) = dxldt(5,loc_disc) + J_NHE_A./w_A + J_buf_A;
  dxldt(4,loc_disc) = J_NHE_A/w_A + J_buf_A;
  // % CO_A mM/s
  // dxldt(6,loc_disc) = dxldt(6,loc_disc) + J_CDF_A./w_A - J_buf_A;
  dxldt(5,loc_disc) = J_CDF_A/w_A - J_buf_A;

#ifdef DEBUGFODE
  out << std::scientific << std::setprecision(8);
  out << "================ DEBUG =================" << std::endl;
  out << "I_ENaC:       %.8d nA  " << I_ENaC*1e-6 << std::endl;
  out << "I_BK:         %.8d nA  " << I_BK*1e-6 << std::endl;
  out << "I_K_B:        %.8d nA  " << I_K_B*1e-6 << std::endl;
  out << "J_NKA_A:      %.8d nA  " << J_NKA_A*F_const*1e-3 << std::endl;
  out << "J_NKA_B:      %.8d nA  " << J_NKA_B*F_const*1e-3 << std::endl;
  out << "I_CFTR:       %.8d nA  " << I_CFTR*1e-6 << std::endl;
  out << "I_CFTR_B:     %.8d nA  " << I_CFTR_B*1e-6 << std::endl;
  out << "J_AE2_A:      %.8d nA  " << J_AE2_A*F_const*1e-9 << std::endl;
  out << "J_AE2_B:      %.8d nA  " << J_AE2_B*F_const*1e-9 << std::endl;
  out << "J_NBC:        %.8d nA  " << J_NBC*F_const*1e-9 << std::endl;
  out << "J_NHE_A:      %.8d nA  " << J_NHE_A*F_const*1e-9 << std::endl;
  out << "J_NHE_B:      %.8d nA  " << J_NHE_B*F_const*1e-9 << std::endl;
  out << "J_CDF_A:      %.8d nA  " << J_CDF_A*F_const*1e-9 << std::endl;
  out << "J_CDF_B:      %.8d nA  " << J_CDF_B*F_const*1e-9 << std::endl;
  out << "J_buf_A:      %.8d nA  " << J_buf_A*F_const*w_A*1e-9 << std::endl;
  out << "J_buf_C:      %.8d nA  " << J_buf_C*F_const*1e-9 << std::endl;
  out << "J_A:          %.8d nA  " << J_A << std::endl;
  out << "J_B:          %.8d nA  " << J_B << std::endl;
  out << " " << std::endl;
  out << "V_A_K:     %.8d mV " <<  V_A_K << std::endl;
  out << "V_B_K:     %.8d mV " <<  V_B_K << std::endl;
  out << "V_A_Cl:    %.8d mV " <<  V_A_Cl << std::endl;
  out << "V_A_Na:    %.8d mV " <<  V_A_Na << std::endl;
  out << "V_A:       %.8d mV " <<  x_c(0,0) << std::endl;
  out << "V_B:       %.8d mV " <<  x_c(0,1) << std::endl;
  out << "V_T:       %.8d mV " <<  V_T << std::endl;
  out << " " << std::endl;
  out << "V_P_Na:    %.8d mV " <<  V_P_Na << std::endl;
  out << "V_P_K:     %.8d mV " <<  V_P_K << std::endl;
  out << "V_P_Cl:    %.8d mV " <<  V_P_Cl << std::endl;
  out << "I_P_Na:    %.8d nA " << I_P_Na*1e-6 << std::endl;
  out << "I_P_K:     %.8d nA " << I_P_K*1e-6 << std::endl;
  out << "I_P_Cl:    %.8d nA " << I_P_Cl*1e-6 << std::endl;
  out << " " << std::endl;
  out << "Na flux A: %.8d nA " << I_ENaC*1e-6 - J_NHE_A*F_const*1e-9 + I_P_Na*1e-6 + 3*J_NKA_A*1e-3*F_const << std::endl;
  out << "K flux A:  %.8d nA " << I_BK*1e-6 + I_P_K*1e-6 - 2*J_NKA_A*1e-3*F_const << std::endl;
  out << "Cl flux A: %.8d nA " <<  I_P_Cl*1e-6 + I_CFTR*1e-6 + J_AE2_A*F_const*1e-9 << std::endl;
  out << "HC flux A: %.8d nA " <<  I_CFTR_B*1e-6 - J_AE2_A*F_const*1e-9 - J_buf_A*F_const*w_A*1e-9 << std::endl;
  out << "================ END DEBUG =================" << std::endl;
#endif
}

const double cSICell::compute_electroneutrality_check() {
  // in mol
  double e = (x_c(3) + x_c(4) - x_c(5) - x_c(6) + x_c(7)) * x_c(2) * 1e-18 - 1.5 * P.chi_C;

  return e;
}
