/*
 * cDuctSegmentAcinus.hpp
 *
 *  Created on: 15/7/2021
 *      Author: jrugis
 */

#ifndef CDUCTSEGMENTACINUS_H_
#define CDUCTSEGMENTACINUS_H_

#include <vector>
#include <fstream>

#include "global_defs.hpp"
#include "cDuctSegment.hpp"

class cMiniGlandDuct;

class cDuctSegmentAcinus : public cDuctSegment {
  public:
  cDuctSegmentAcinus(cMiniGlandDuct* parent, int seg_number);
  ~cDuctSegmentAcinus();
  virtual void step(double t, double dt);
  //EIGEN_MAKE_ALIGNED_OPERATOR_NEW   // required when using fixed-size vectorizable Eigen object(s)

  protected:
  //void get_parameters();
  //void save_results();

  //int stepnum, outputnum, Tstride;
  //std::string resultsh5_dataset, resultsh5_filename;
};

#endif /* CDUCTSEGMENTACINUS_H_ */
