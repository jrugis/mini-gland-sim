/*
 * cSCell.hpp
 *
 *  Created on: 21/4/2021
 *      Author: jrugis, cscott
 */

#ifndef CSCELL_H_
#define CSCELL_H_

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

class cCMesh;

class cSCell {

public:
  cSCell(cDuct* _parent, int _cell_number);
  ~cSCell();
  scell::Array1CC x_c, dxcdt;   // solution vectors for ions
  duct::ArrayNFC dxldt;
  Array1Nd dwAdt;
  void step();
  const double get_min_z() { return min_z; }
  const double get_max_z() { return max_z; }
  const double get_mean_z() { return mean_z; }
  void process_mesh_info(std::vector<double>& lumen_segment);
  void f_ODE(const duct::ArrayNFC &xl, const duct::lumen_prop_t &lumen_prop);
  void init(duct::parameters_t &parent_P);
  const double compute_electroneutrality_check();

protected:
	  
private:
  cDuct* parent;
  std::string id;
  std::ofstream out;   // runtime diagnostic file for this object
  INIReader *p;        // model parameters
  duct::parameters_t P;

  int cell_number;     // this cell number
  cCMesh *mesh;        // this cell mesh
  double api_area, baslat_area;      // surface areas for different regions
  int napical;                       // number of apical triangles
  double min_z, max_z, mean_z;
  std::vector<int> api_lumen_conn;
  std::vector<double> api_face_area;
  std::vector<int> loc_int;
  Array1Nd api_area_int;
  scell::scaled_rates_t scaled_rates;
  void setup_parameters(duct::parameters_t &parent_P);
  void setup_IC();
  void init_const();
};


#endif /* CSCELL_H_ */
