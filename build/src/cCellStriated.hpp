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
  typedef Eigen::Array<double, 1, IONCOUNT> Array1IC;

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
}

class cCellStriated : public cCell {
  public:
  cCellStriated(cDuctSegment* parent, int cell_number);
  void process_geom(std::vector<double>& lumen_segment);
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW   // required when using fixed-size vectorizable Eigen object(s)

  private:
  S::Array1IC solvec, prev_solvec;   // solution vectors for ions
  S::constant_values c;              // secretion constants vector
  double api_area, baslat_area;      // surface areas for different regions
  int napical;                       // number of apical triangles
  MatrixN1d api_face_areas;
  std::vector<double> api_face_area;
  void init_solvec();
  void init_const();
};

#endif /* CCELLSTRIATED_H_ */
