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

cCellStriated::cCellStriated(cDuctSegment* parent, int cell_number) : cCell(parent, cell_number)
{
  out << "<CellStriated> @constructor" << std::endl;

  // calc areas of different regions
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
}

void cCellStriated::process_geom(std::vector<double>& lumen_segment) {
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
  std::vector<int> api_lumen_conn(napical);
//  out << "discretize:" << std::endl;
  for (int i = 0; i < napical; i++) {
    int nsegments = lumen_segment.size() - 1;
    for (int j = 0; j < nsegments; j++) {
      if ((lumen_segment[j] <= api_coord_z_mean[i]) && (api_coord_z_mean[i] < lumen_segment[j+1])) {
        api_lumen_conn[i] = j;
//        out << "  " << i << ": " << api_coord_z_mean[i] << " -> " << j << std::endl;
        break;
      }
    }
  }

  // find the unique lumen segments for this cell
  std::unordered_set<int> loc_int_set(api_lumen_conn.begin(), api_lumen_conn.end());
  std::vector<int> loc_int(loc_int_set.begin(), loc_int_set.end());
  std::sort(loc_int.begin(), loc_int.end());
  out << loc_int.size() <<  " unique lumen segments for this cell:";
  for (auto const &e: loc_int) {
    out << " " << e;
  }
  out << std::endl;

  // loop through the lumen segments to calculate areas per segment
  int n_loc_int = loc_int.size();
  std::vector<double> api_int_area(n_loc_int, 0.0);
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


}

void cCellStriated::init_const(){
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

void cCellStriated::init_solvec(){
  //solvec(Nal) = p.at("Nal0");
  //solvec(Kl) = p.at("Kl0");
  //solvec(Cll) = p.at("Nal0") + p.at("Kl0");
  //solvec(VOL) = s.V0;
  //solvec(Na) = p.at("Na0");
  //solvec(K) = p.at("K0");
  //solvec(Cl) = p.at("Cl0");
  //solvec(HCO3) = p.at("HCO30");
  //solvec(H) = 1e3 * pow(10, -(p.at("pHi")));
  //solvec(Va) = p.at("Va0");
  //solvec(Vb) = p.at("Vb0");
}
