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

cCellStriated::cCellStriated(cDuctSegment* parent, int cell_number) : cCell(parent, cell_number)
{
  out << "<CellStriated> @constructor" << std::endl;

  get_parameters();

  init_const();

  setup_IC();
}

void cCellStriated::get_parameters() {
  // apical channel conductances
  P.G_ENaC = p.at("G_ENaC");
  P.G_CFTR = p.at("G_CFTR");
  P.G_BK = p.at("G_BK");

  // basolateral channel conductances
  P.G_K_B = p.at("G_K_B");

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
  api_int_area.resize(n_loc_int, 0.0);
  for (int i = 0; i < napical; i++) {
    auto it = std::find(loc_int.begin(), loc_int.end(), api_lumen_conn[i]);
    int idx = std::distance(loc_int.begin(), it);
    api_int_area[idx] += api_face_area[i];
  }
  out << "lumen segment areas:";
  for (auto const &e: api_int_area) {
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
  xc(V_A) = -26.1257;  // TODO: where should these 3 be read from?
  xc(V_B) = -52.2513;
  xc(w_C) = 1000;
  xc(Na_C) = 17;  // TODO: move these to parameter file
  xc(K_C) = 140;
  xc(Cl_C) = 22;
  xc(HCO_C) = 75;
  xc(H_C) = 1000 * pow(10, -7.35);
  xc(CO_C) = 1.28;
}
