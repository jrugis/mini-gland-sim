/*
 * cCell.hpp
 *
 *  Created on: 6/5/2021
 *      Author: jrugis
 */

#ifndef CCELLSTRIATED_H_
#define CCELLSTRIATED_H_

#include "global_defs.hpp"

class cDuctSegment;

namespace S{
  // solution vector components
  enum solution_values { Nal, Kl, Cll, VOL, Na, K, Cl, HCO3, H, Va, Vb, IONCOUNT };
  typedef Eigen::Array<double, 1, IONCOUNT> Array1IC;

  // invariant cell properties
  class constant_values {
    public:
    double aNaK, aNkcc1, GtNa, GtK, GCl, GK, G1, G4, GB, St, Sb, Sa, V0;
  };
}

class cCellStriated : public cCell {
  public:
  cCellStriated(cDuctSegment* parent, int cell_number);
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW   // required when using fixed-size vectorizable Eigen object(s)
  virtual void step();

  private:
  S::Array1IC solvec, prev_solvec;   // solution vectors for ions
  S::constant_values s;              // secretion constants vector
  void init_solvec();
  void init_const();
};

#endif /* CCELLSTRIATED_H_ */
