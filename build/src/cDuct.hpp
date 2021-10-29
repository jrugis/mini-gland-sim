/*
 * cDuct.hpp
 *
 *  Created on: 21/4/2021
 *      Author: jrugis, cscott
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
    double alpha_A;
    double alpha_B;
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
}

class cMiniGland;
class cSICell;

class cDuct {
  friend cSICell; // child object

public:
  cDuct(cMiniGland* parent);
  ~cDuct();
  void step(double t, double dt);
  int get_nvars();
  void f_ODE(const double t, const Array1Nd &x_in, Array1Nd &dxdt);

protected:
  std::string id;
  std::ofstream out;    // runtime diagnostic file for this object
  INIReader* p;         // model parameters
  int nlsegs;           // number of duct lumen segments
  int ilsegs;           // index of first lumen segment

  int n_disc;
  Array1Nd disc_length;
  Array1Nd disc_volume;
  Array1Nd disc_X_area;
  Array1Ni disc_out_Vec;

private:
  cMiniGland* parent;
  std::vector<cSICell*> icells;    // the intercalated cells associated with this duct
  std::vector<cSICell*> scells;    //     striated 
  double L_int;   // Lumen segment discretisation interval
  double PSflow;  // um3/s volumetric primary saliva flow rate
  duct::conc_t Conc;
  duct::parameters_t P;
  duct::ArrayNFC x_l, dxldt;  // solution vector and derivative
  Array1Nd x, dxdt;
  Array1Nd dwAdt, v, v_up;
  duct::ArrayNFC x_up;
  cCVode *solver;
  int stepnum, outputnum, Tstride;
  std::string resultsh5_dataset, resultsh5_filename;
  void get_parameters();
  void process_mesh_info();
  void setup_IC();
  void distribute_x(const Array1Nd &x_in);
  void gather_x(Array1Nd &x_out);
  void save_results();
  void setup_arrays();
  double accum_fluid(const int duct_idx);
};

#endif /* CDUCT_H_ */
