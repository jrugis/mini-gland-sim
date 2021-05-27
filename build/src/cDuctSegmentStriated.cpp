/*
 * cDuctSegmentStriated.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

//#include <thread>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <omp.h>

#include "cCell.hpp"
#include "cDuctSegment.hpp"
#include "cDuctSegmentStriated.hpp"
#include "cCellStriated.hpp"

using namespace dss;

cDuctSegmentStriated::cDuctSegmentStriated(cMiniGlandDuct* parent, int seg_number) : cDuctSegment(parent, seg_number) {
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

  Conc.PS(Na) = 143.5;  // concentration of primary saliva
  Conc.PS(K) = 5.2;
  Conc.PS(Cl) = 114.5;
  Conc.PS(HCO) = 34.2;
  Conc.PS(H) = 1000 * pow(10, -7.35);
  Conc.PS(CO) = 1.28;

  Conc.CIC(Na) = 17;  // cellular initial concentration
  Conc.CIC(K) = 140;
  Conc.CIC(Cl) = 22;
  Conc.CIC(HCO) = 75;
  Conc.CIC(H) = 1000 * pow(10, -7.35);
  Conc.CIC(CO) = 1.28;

  Conc.LIC(Na) = 143.5;  // lumenal initial concentration
  Conc.LIC(K) = 5.2;
  Conc.LIC(Cl) = 114.5;
  Conc.LIC(HCO) = 34.2;
  Conc.LIC(H) = 1000 * pow(10, -7.35);
  Conc.LIC(CO) = 1.28;

  // create parameters structure
  get_parameters();

  // mesh stuff (each cell has it's own mesh data, lumen segments here)
  get_lumen_geom(L_int);
  get_cell_geom();


}

void cDuctSegmentStriated::get_cell_geom() {
  int ncells = cells.size();
  for (int i = 0; i < ncells; i++) {
    static_cast<cCellStriated*>(cells[i])->process_geom(lumen_prop.segment);
  }

}

void cDuctSegmentStriated::get_lumen_geom(double L) {
  out << "<cDuctSegmentStriated> Lumen min/max in z: " << vertex_in(2) << " - " << vertex_out(2) << std::endl;

  lumen_prop.L = L;  // discretisation interval

  double lumen_start = std::min(vertex_out(2), vertex_in(2));
  double lumen_end = std::max(vertex_out(2), vertex_in(2));
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
  // TODO: should be read from parameter file?

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

void cDuctSegmentStriated::step()
{
  // combine cells fluid flow  --  TO DO
  // ....

  out << "<DuctSegmentStriated> step - threads in use: " << omp_get_num_threads() << std::endl;
}

