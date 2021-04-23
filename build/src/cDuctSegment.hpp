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
  virtual ~cDuctSegment();

  private:
  cMiniGlandDuct* parent;
  std::string id;
  std::ofstream out;                         // runtime diagnostic file for this object
  std::unordered_map<std::string, double> p; // the model parameters
  int seg_number;                          // this segment number
  int seg_type;                            // type of this duct segment 

  Vector3d vertex_in, vertex_out;          // input and output center points for this duct segment
  double inner_diameter, outer_diameter;   // inner and outer diameter of this duct segment
  std::vector<cCell*> cells;               // the cells associated with this duct segment
  
  std::tuple<int,int> get_segment_data();
  virtual void run();
  virtual void create_cells(std::tuple<int,int> cell_info);
};

class cDuctSegmentAcinus : public cDuctSegment {
  private:
  virtual void run();
  virtual void create_cells(std::tuple<int,int> cell_info);
};

class cDuctSegmentIntercalated : public cDuctSegment {
  private:
  virtual void run();
  virtual void create_cells(std::tuple<int,int> cell_info);
};

class cDuctSegmentStriated : public cDuctSegment {
  private:
  virtual void run();
  virtual void create_cells(std::tuple<int,int> cell_info);
};

#endif /* CDUCTSEGMENT_H_ */
