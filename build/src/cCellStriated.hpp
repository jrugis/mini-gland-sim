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

class cDuctSegment;

namespace S{
  // solution vector
  enum solution_values { Nal, Kl, Cll, VOL, Na, K, Cl, HCO3, H, Va, Vb, IONCOUNT };
  enum cellular_values { V_A, V_B, w_C, Na_C, K_C, Cl_C, HCO_C, H_C, CO_C, CELLULARCOUNT };
  typedef Eigen::Array<double, 1, CELLULARCOUNT> Array1CC;

  // invariant cell properties
  struct constant_values {
    double aNaK;
    double aaNkcc1;
    double aGtNa;
    double aGtK;
    double aGCl;
    double aGK;
    double aG1;
    double aG4;
    double aGB;
    double aSt;
    double aSb;
    double aSa;
    double aV0;
  };

  struct nka_t {
    double alpha_A;
    double alpha_B;
    double r;
    double beta;
  };

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
    nka_t NKA;
  };

  struct parameters_t {
    // apical channel conductances
    double G_ENaC;
    double G_CFTR;
    double G_BK;

    // basolateral channels conductances
    double G_K_B;
    
    // paracellular conductances
    double G_P_Na;
    double G_P_K;
    double G_P_Cl;

    // water permeability across membranes
    double L_A;
    double L_B;

    // sodium potassium pump rates
    nka_t NKA;
  };
}

class cCellStriated : public cCell {
  public:
  cCellStriated(cDuctSegment* parent, int cell_number);
  void process_mesh_info(std::vector<double>& lumen_segment);
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW   // required when using fixed-size vectorizable Eigen object(s)

  private:
  S::parameters_t P;
  S::Array1CC xc, dxcdt;   // solution vectors for ions
  S::constant_values c;              // secretion constants vector
  double api_area, baslat_area;      // surface areas for different regions
  int napical;                       // number of apical triangles
  std::vector<int> api_lumen_conn;
  std::vector<double> api_face_area;
  std::vector<int> loc_int;
  std::vector<double> api_int_area;
  S::scaled_rates_t scaled_rates;
  void setup_IC();
  void init_const();
  void get_parameters();
};

#endif /* CCELLSTRIATED_H_ */
