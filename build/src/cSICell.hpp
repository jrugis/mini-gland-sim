/*
 * cSICell.hpp
 *
 *  Created on: 26/09/2021
 *      Author: jrugis, cscott
 */

#ifndef CSICELL_H_
#define CSICELL_H_

#include <fstream>
#include <string>
#include <vector>

#include "cDuct.hpp"
#include "inih/cpp/INIReader.h"
#include "global_defs.hpp"

#define CELLULARCOUNT 9

namespace scell {
  // solution vector
  typedef Eigen::Array<double, 1, CELLULARCOUNT> Array1CC;

  struct scaled_rates_t {
    double L_A;
    double L_B;
    double G_ENaC;
    double G_CFTR;
    double G_BK;
    double G_K_B;
    double G_P_Na;
    double G_P_K;
    double G_P_Cl;
    duct::nka_t NKA;
  };
}

class cSIMesh;

class cSICell {

public:
  cSICell(cDuct* _parent, std::string _fname, cell_groups _cell_type);
  ~cSICell();
  scell::Array1CC x_c, dxcdt;   // solution vectors for ions
  duct::ArrayNFC dxldt;
  Array1Nd dwAdt;
  void step();
  const double get_min_z() { return min_z; }
  const double get_max_z() { return max_z; }
  const double get_mean_z() { return mean_z; }
  void process_mesh_info();
  void f_ODE(const duct::ArrayNFC &xl);
  void setup(duct::parameters_t &parent_P);
  const double compute_electroneutrality_check();

protected:
	  
private:
  cDuct* parent;
  std::string id;
  std::ofstream out;   // runtime diagnostic file for this object
  INIReader *p;        // model parameters
  duct::parameters_t P;

  int cell_number;              // this cell number
  cell_groups cell_group;       // STRIATED or INTERCALATED
  cSIMesh *mesh;                // this cell mesh
  double api_area, baslat_area; // surface areas for different regions
  int napical;                  // number of apical triangles
  double min_z, max_z, mean_z;
  std::vector<int> api_lumen_conn;
  std::vector<double> api_face_area;
  std::vector<int> loc_int;
  Array1Nd api_area_int;
  scell::scaled_rates_t scaled_rates;
  void setup_parameters(duct::parameters_t &parent_P);
  void setup_IC();
  void setup_arrays();
  Array1Nd Na_A, K_A, Cl_A, HCO_A, H_A, CO_A;  // intermediate calculation arrays as members for performance
  Array1Nd J_A, J_CDF_A, J_buf_A, J_NHE_A, J_AE2_A;
  Array1Nd V_A_Cl, I_CFTR, V_A_HCO, I_CFTR_B, V_A_K;
  Array1Nd I_BK, V_A_Na, I_ENaC, J_NKA_A, V_P_Na, V_P_K, V_P_Cl;
  Array1Nd I_P_Na, I_P_K, I_P_Cl;
};


#endif /* CSICELL_H_ */
