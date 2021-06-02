/*
 * cCellStriated.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#include <chrono> // std::chrono::seconds
#include <thread> // std::this_thread::sleep_for
#include <unordered_set>
#include <algorithm>

#include "cCell.hpp"
#include "cCellStriated.hpp"

using namespace S;
using namespace dss;

cCellStriated::cCellStriated(cDuctSegment* parent, int cell_number) : cCell(parent, cell_number)
{
  out << "<CellStriated> @constructor" << std::endl;
}

void cCellStriated::init(dss::parameters_t &parent_P) {
  // init in separate function so parent initialiser has completed
  setup_parameters(parent_P);
  init_const();
  setup_IC();
}

void cCellStriated::setup_parameters(dss::parameters_t &parent_P) {
  // copy of parents local parameters object
  P = parent_P;
}

void cCellStriated::process_mesh_info(std::vector<double>& lumen_segment) {
  // mean z coordinate of each apical triangle vertex
  std::vector<double> api_coord_z_mean(napical);
  int api_count = 0;
  for (int i = 0; i < nfaces; i++) {
    if (face_types(i) == 0) {  // apical
      double zsum = 0.0;
      for (int j = 0; j < 3; j++) {  // for each vertex
        int vertidx = faces(i, j);
        zsum += verts(vertidx, 2);
      }
      api_coord_z_mean[api_count++] = zsum / 3.0;
    }
  }

  // sort apical triangles into corresponding lumen segments
  api_lumen_conn.resize(napical, -1);
  for (int i = 0; i < napical; i++) {
    int nsegments = lumen_segment.size() - 1;
    for (int j = 0; j < nsegments; j++) {
      if ((lumen_segment[j] <= api_coord_z_mean[i]) && (api_coord_z_mean[i] < lumen_segment[j+1])) {
        api_lumen_conn[i] = j;
        break;
      }
    }
  }

  // find the unique lumen segments for this cell
  std::unordered_set<int> loc_int_set(api_lumen_conn.begin(), api_lumen_conn.end());
  loc_int = std::vector<int>(loc_int_set.begin(), loc_int_set.end());
  std::sort(loc_int.begin(), loc_int.end());
  out << loc_int.size() <<  " unique lumen segments for this cell:";
  for (auto const &e: loc_int) {
    out << " " << e;
  }
  out << std::endl;

  // loop through the lumen segments to calculate areas per segment
  int n_loc_int = loc_int.size();
  api_area_int.resize(Eigen::NoChange, n_loc_int);
  api_area_int.setZero();
  for (int i = 0; i < napical; i++) {
    auto it = std::find(loc_int.begin(), loc_int.end(), api_lumen_conn[i]);
    int idx = std::distance(loc_int.begin(), it);
    api_area_int(0, idx) += api_face_area[i];
  }
  out << "lumen segment areas:";
  for (auto const &e: api_area_int.reshaped()) {
    out << "  " << e;
  }
  out << std::endl;

  // scale the rates based on cell surface areas
  double A = 104.719755;  // the apical area used to scale the conductances G
  double A_A = api_area;
  double A_B = baslat_area;
  scaled_rates.L_A = P.L_A * A / A_A;
  scaled_rates.L_A = P.L_B * A / A_B;
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
}

void cCellStriated::init_const(){
  // calc areas of different surface regions
  api_area = 0.0;
  baslat_area = 0.0;
  napical = 0;
  for (int i = 0; i < nfaces; i++) {
    if (face_types(i) == 0) {  // apical
      api_area += face_areas(i);
      api_face_area.push_back(face_areas(i));
      napical++;
    }
    else {  // basal or basolateral
      baslat_area += face_areas(i);
    }
  }
  out << "<CellStriated> num apical triangles = " << napical << std::endl;
  out << "<CellStriated> api_area = " << api_area << std::endl;
  out << "<CellStriated> baslat_area = " << baslat_area << std::endl;

  // initial cell volume
  //s.V0 = parent->element_data.col(VOL_e).sum();  

  // apical surface area
  //s.Sa = 0.0;   
  //for (int n = 0; n < parent->mesh->apical_triangles_count; n++) {
  //  int apical_tri = parent->mesh->apical_triangles(n);
  //  s.Sa += parent->surface_data(apical_tri, AREA_s);
  //}

  // basal surface area
  //s.Sb = 0.0;   
  //for (int n = 0; n < parent->mesh->basal_triangles_count; n++) {
  //  int basal_tri = parent->mesh->basal_triangles(n);
  //  s.Sb += parent->surface_data(basal_tri, AREA_s);
  //}

  //s.St = s.Sa / 0.943551157250391;
  //s.aNkcc1 = ( 0.0063812 ) * s.Sb;
}

void cCellStriated::setup_IC(){
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

void cCellStriated::f_ODE(const dss::ArrayNFC &x_l, const dss::lumen_prop_t &lumen_prop, Array1Nd &dwAdt) {
  // constant parameters
  double Na_B = P.ConI(dss::Na);
  double K_B = P.ConI(dss::K);
  double Cl_B = P.ConI(dss::Cl);
  double HCO_B = P.ConI(dss::HCO);
  double H_B = P.ConI(dss::H);
  double CO_B = P.ConI(dss::CO);

  double w_A = lumen_prop.volume;
  double L = lumen_prop.L;
  double A_L = lumen_prop.X_area;
  double chi_C = P.chi_C; // mol
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

  double alpha_NBC = P.NBC.alpha;
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
  double A_A = api_area;
  double A_B = baslat_area;

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

  // lumenal variables of all the lumen segments this cell interfaces with
  Array1Nd Na_A(1, loc_int.size());
  Na_A = x_l(0, loc_int);
  Array1Nd K_A(1, loc_int.size());
  K_A = x_l(1, loc_int);
  Array1Nd Cl_A(1, loc_int.size());
  Cl_A = x_l(2, loc_int);
  Array1Nd HCO_A(1, loc_int.size());
  HCO_A = x_l(3, loc_int);
  Array1Nd H_A(1, loc_int.size());
  H_A = x_l(4, loc_int);
  Array1Nd CO_A(1, loc_int.size());
  CO_A = x_l(5, loc_int);

  // water transport
  // osm_c = chi_C./w_C*1e18; % osmolarity of cell due to proteins (chi)
  //
  // J_B = 1e-18*L_B.*V_w.*(Na_C + K_C + Cl_C + HCO_C + osm_c - Na_B - K_B - Cl_B - HCO_B - phi_B); % um/s 
  // J_A = 1e-18*L_A.*V_w.*(Na_A + K_A + Cl_A + HCO_A + phi_A - Na_C - K_C - Cl_C - HCO_C - osm_c); % um/s [1, n_loc_int]
  // dwdt = A_B * J_B - sum(A_A_int .* J_A); % um^3/s
  // dwAdt(1,loc_int) = dwAdt(1,loc_int) + A_A_int .* J_A; % um^3/s [1, n_loc_int]
  double osm_c = P.chi_C / w_C;
  double J_B = 1e-18 * L_B * V_w * (Na_C + K_C + Cl_C + HCO_C + osm_c - Na_B - K_B - Cl_B - HCO_B - phi_B);
  Array1Nd J_A(1, loc_int.size());
  J_A = 1e-18 * L_A * V_w * (Na_A + K_A + Cl_A + HCO_A + phi_A - Na_C - K_C - Cl_C - HCO_C - osm_c);
  double dwdt = A_B * J_B * (api_area_int * J_A).sum();
  dwAdt(0, loc_int) += api_area_int * J_A;




}
