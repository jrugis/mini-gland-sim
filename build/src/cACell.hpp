/*
 * cACell.hpp
 *
 *  Created on: 2/8/2021
 *      Author: jrugis
 */

#ifndef CACELL_H_
#define CACELL_H_

#include <fstream>
#include <string>

#include "inih/cpp/INIReader.h"
#include "global_defs.hpp"

class cAcinus;
class cCMesh;

#define DIFVARS 3                        // number of diffusing node variables - c, ip
//#define NONDIFVARS 2                     // number of non-diffusing variables - g, h
//#define VARIABLES (DIFVARS + NONDIFVARS) // total number of node variables
//#define REF_MASS_SIZE 4                  // reference mass dimension

// some convenience typedefs
//typedef Eigen::Array<double, Eigen::Dynamic, 1> ArrayX1C;
typedef Eigen::Array<double, 1, DIFVARS> Array1VC;
//typedef Eigen::Array<double, REF_MASS_SIZE, REF_MASS_SIZE> ArrayRefMass;
//typedef Eigen::Triplet<double> Triplet;

class cACell {
public:
  cACell(cAcinus* _parent, int _cell_number);
  ~cACell();
  void step();

protected:
	  
private:
  cAcinus* parent;
  std::string id;
  std::ofstream out;   // runtime diagnostic file for this object
  INIReader *p;        // model parameters

  int cell_number;     // this cell number
  cCMesh *mesh;        // this cell mesh
  
  void make_matrices();
  //void initialise();
  
  //MatrixN1d get_load(double delta_time, bool plc);
  //Array1VC get_rhs_reactions(double c, double ip, double ce, double g, double ryr_f, double plc_f);
  //Array1VC get_bndy_reactions(double c, double ip, double ce, double g, double ryr_f, double plc_f);
  //Array1VC get_conductances(double c, double ip, double ce, double g, double ryr_f, double plc_f);
  //void secretion();
};


#endif /* CACELL_H_ */
