/*
 * cDuctSegment.hpp
 *
 *  Created on: 21/4/2021
 *      Author: jrugis
 */

#ifndef CDUCTSEGMENT_H_
#define CDUCTSEGMENT_H_

#include <fstream>
#include <string>
#include <unordered_map>
//#include <vector>

class cMiniGlandDuct;

class cDuctSegment {
  public:
  cDuctSegment(cMiniGlandDuct* parent, int seg_number);
  ~cDuctSegment();
  void run();

  private:
  cMiniGlandDuct* parent;
  std::string id;
  int seg_number;       // this segment number                 
  std::ofstream out;    // runtime diagnostic file for this object
  std::unordered_map<std::string, double> p; // the model parameters
};

#endif /* CDUCTSEGMENT_H_ */
