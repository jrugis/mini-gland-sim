/*
 * cCell.hpp
 *
 *  Created on: 6/5/2021
 *      Author: jrugis
 */

#ifndef CCELLSTRIATED_H_
#define CCELLSTRIATED_H_

#include <vector>
#include "global_defs.hpp"
#include "cCell.hpp"
#include "cDuctSegmentStriated.hpp"

class cDuctSegment;

#define CELLULARCOUNT 9

namespace S{
  // solution vector
//  enum solution_values { Nal, Kl, Cll, VOL, Na, K, Cl, HCO3, H, Va, Vb, IONCOUNT };
//  enum cellular_values { V_A, V_B, w_C, Na_C, K_C, Cl_C, HCO_C, H_C, CO_C, CELLULARCOUNT };
  typedef Eigen::Array<double, 1, CELLULARCOUNT> Array1CC;

  // invariant cell properties
//  struct constant_values {
//    double aNaK;
//    double aaNkcc1;
//    double aGtNa;
//    double aGtK;
//    double aGCl;
//    double aGK;
//    double aG1;
//    double aG4;
//    double aGB;
//    double aSt;
//    double aSb;
//    double aSa;
//    double aV0;
//  };

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
    dss::nka_t NKA;
  };
}

class cCellStriated : public cCell {
  public:
  cCellStriated(cDuctSegment* parent, int cell_number);
  S::Array1CC x_c, dxcdt;   // solution vectors for ions
  void process_mesh_info(std::vector<double>& lumen_segment);
  void f_ODE(const dss::ArrayNFC &xl, const dss::lumen_prop_t &lumen_prop, dss::ArrayNFC &dxldt, Array1Nd &dwAdt);
  void init(dss::parameters_t &parent_P);
  const double get_min_z() { return min_z; }
  const double get_max_z() { return max_z; }
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW   // required when using fixed-size vectorizable Eigen object(s)

  private:
  dss::parameters_t P;
  double api_area, baslat_area;      // surface areas for different regions
  int napical;                       // number of apical triangles
  std::vector<int> api_lumen_conn;
  std::vector<double> api_face_area;
  std::vector<int> loc_int;
  //std::vector<double> api_area_int;
  Array1Nd api_area_int;
  S::scaled_rates_t scaled_rates;
  double min_z, max_z;
  void setup_parameters(dss::parameters_t &parent_P);
  void setup_IC();
  void init_const();
};

#endif /* CCELLSTRIATED_H_ */
