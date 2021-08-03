/*
 * cDuct.hpp
 *
 *  Created on: 21/4/2021
 *      Author: jrugis
 */

#ifndef CDUCT_H_
#define CDUCT_H_

#include <fstream>
#include <string>
#include <vector>

#include "inih/cpp/INIReader.h"
#include "cCVode.hpp"
#include "global_defs.hpp"

// number of lumenal variables
#define LUMENALCOUNT 6

namespace duct {
  // solution values enums
//  enum lumenal_values { Na_A, K_A, Cl_A, HCO_A, H_A, CO_A, LUMENALCOUNT };
  enum field_values { Na, K, Cl, HCO, H, CO, FIELDCOUNT };
  typedef Eigen::Array<double, 1, LUMENALCOUNT> Array1LC;
  typedef Eigen::Array<double, 1, FIELDCOUNT> Array1FC;
  typedef Eigen::Array<double, FIELDCOUNT, Eigen::Dynamic> ArrayNFC;

  // parameters
  struct nbc_t {
    double alpha;
    double k5_p;
    double k5_m;
    double k6_p;
    double k6_m;
  };

  struct ae2_t {
    double alpha_A;
    double alpha_B;
    double k3_p;
    double k3_m;
    double k4_p;
    double k4_m;
  };

  struct nhe_t {
    double alpha_A;
    double alpha_B;
    double k1_p;
    double k1_m;
    double k2_p;
    double k2_m;
  };

  struct buf_t {
    double k_p;
    double k_m;
  };

  struct nka_t {
    double alpha_A;
    double alpha_B;
    double r;
    double beta;
  };

  struct parameters_t {
    Array1FC ConI;
    Array1FC ConP;
    double PSflow;

    // apical channel conductances
    double G_ENaC;
    double G_CFTR;
    double G_BK;

    // basolateral channels conductances
    double G_K_B;

    // apical or basolateral transporter rates
    nbc_t NBC;
    ae2_t AE2;
    nhe_t NHE;

    // CO2 permeability
    double p_CO;

    // CO2 bicarbonate buffering
    buf_t buf;

    // sodium potassium pump rates
    nka_t NKA;

    // paracellular conductances
    double G_P_Na;
    double G_P_K;
    double G_P_Cl;

    // water permeability across membranes
    double L_A;
    double L_B;

    // universal physical constants
    double V;

    // osmolarity adjusting constants
    double chi_C;
    double phi_A;
    double phi_B;
  };

  // concentrations
  struct conc_t {
    Array1FC Int;  // interstium
    Array1FC PS;  // primary saliva
    Array1FC CIC;  // cellular initial concentration
    Array1FC LIC;  // lumenal initial concentration
  };

  // lumen properties
  struct lumen_prop_t {
    std::vector<double> segment;
    double volume;
    double X_area;
    int n_int;
    double L;
  };
}

class cMiniGland;
class cSCell;

class cDuct {
  friend cSCell; // child object

public:
  cDuct(cMiniGland* parent);
  ~cDuct();
  void step(double t, double dt);
  int get_nvars();
  void f_ODE(const Array1Nd &x_in, Array1Nd &dxdt);

protected:
  std::string id;
  std::ofstream out;    // runtime diagnostic file for this object
  INIReader* p;         // model parameters
  int nlsegs;           // number of duct lumen segments
  int ilsegs;           // index of first lumen segment

private:
  cMiniGland* parent;
  std::vector<cSCell*> scells;    // the striated cells associated with this duct
  double L_int;  // Lumen segment discretisation interval
  double PSflow;  // um3/s volumetric primary saliva flow rate
  duct::conc_t Conc;
  duct::parameters_t P;
  duct::lumen_prop_t lumen_prop;
  duct::ArrayNFC x_l, dxldt;  // solution vector and derivative
  Array1Nd x, dxdt;
  cCVode *solver;
  int stepnum, outputnum, Tstride;
  std::string resultsh5_dataset, resultsh5_filename;
  void get_parameters();
  void process_mesh_info();
  void setup_IC();
  void distribute_x(const Array1Nd &x_in);
  void gather_x(Array1Nd &x_out);
  void save_results();
};

#endif /* CDUCT_H_ */
