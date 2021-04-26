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

class cDuctSegment {
  friend cCell;                  // child object

  public:
  cDuctSegment(cMiniGlandDuct* parent, int seg_number);
  virtual ~cDuctSegment();
  virtual void step(){};  // redefined in the derived classes

  // the actual in/out fluid flow data variables should go here for the duct object to read & write

  protected:
  cMiniGlandDuct* parent;
  std::string id;
  std::ofstream out;                         // runtime diagnostic file for this object
  std::unordered_map<std::string, double> p; // the model parameters

  int seg_number;                          // this segment number
  int seg_type;                            // type of this duct segment 
  Vector3d vertex_in, vertex_out;          // input and output center points for this duct segment
  double inner_diameter, outer_diameter;   // inner and outer diameter of this duct segment
  std::vector<cCell*> cells;               // the cells associated with this duct segment
};

class cDuctSegmentAcinus : public cDuctSegment {
public:
  cDuctSegmentAcinus(cMiniGlandDuct* parent, int seg_number); 
  virtual void step();	
};

class cDuctSegmentIntercalated : public cDuctSegment {
public:
  cDuctSegmentIntercalated(cMiniGlandDuct* parent, int seg_number);
  virtual void step();	
};

class cDuctSegmentStriated : public cDuctSegment {
public:
  cDuctSegmentStriated(cMiniGlandDuct* parent, int seg_number);
  virtual void step();	
};

#endif /* CDUCTSEGMENT_H_ */
