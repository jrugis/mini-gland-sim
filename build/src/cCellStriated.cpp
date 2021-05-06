/*
 * cCellStriated.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#include <chrono> // std::chrono::seconds
#include <thread> // std::this_thread::sleep_for

#include "cCell.hpp"
#include "cCellStriated.hpp"

cCellStriated::cCellStriated(cDuctSegment* parent, int cell_number) : cCell(parent, cell_number)
{
  out << "<CellStriated> @constructor" << std::endl;
}

void cCellStriated::step()
{
  // TO DO
  // ...
  std::this_thread::sleep_for(std::chrono::milliseconds(200));

  out << "<CellStriated> step - threads in use: " << omp_get_num_threads() << std::endl;
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
