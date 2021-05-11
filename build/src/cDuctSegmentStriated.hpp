/*
 * cDuctSegmentStriated.hpp
 *
 *  Created on: 11/5/2021
 *      Author: jrugis
 */

#ifndef CDUCTSEGMENTSTRIATED_H_
#define CDUCTSEGMENTSTRIATED_H_

class cMiniGlandDuct;

class cDuctSegmentStriated : public cDuctSegment {
  public:
  cDuctSegmentStriated(cMiniGlandDuct* parent, int seg_number);
  virtual void step();
};

#endif /* CDUCTSEGMENTSTRIATED_H_ */
