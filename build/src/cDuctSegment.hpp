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
#include <vector>

#include "global_defs.hpp"

class cMiniGlandDuct;
class cCell;

//enum eSeg_type {acinus, intercalated, striated};

class cDuctSegment {
  friend cCell; // child object
	
  public:
  cDuctSegment(cMiniGlandDuct* parent, int seg_number);
  ~cDuctSegment();
  void run();

  private:
  cMiniGlandDuct* parent;
  std::string id;
  std::ofstream out;                         // runtime diagnostic file for this object
  std::unordered_map<std::string, double> p; // the model parameters

  int seg_number;                          // this segment number
  Vector3d vertex_in, vertex_out;          // input and output center points for this duct segment
  double inner_diameter, outer_diameter;   // inner and outer diameter of this duct segment
  int seg_type;                            // type of this duct segment 
  std::vector<cCell*> cells;               // the cells associated with this duct segment
  
  void get_segment_data();
};

#endif /* CDUCTSEGMENT_H_ */
