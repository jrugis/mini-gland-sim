/*
 * cMiniGlandDuct.hpp
 *
 *  Created on: 21/4/2021
 *      Author: jrugis
 */

#ifndef CMINIGLANDDUCT_H_
#define CMINIGLANDDUCT_H_

#include <fstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

class cDuctSegment;

class cMiniGlandDuct {
  friend cDuctSegment; // child object

  public:
  cMiniGlandDuct();
  ~cMiniGlandDuct();
  void run();

  private:
  std::string id;
  std::ofstream out;                         // the runtime diagnostic file for this object
  std::unordered_map<std::string, double> p; // the model parameters

  std::vector<cDuctSegment*> segments;           // the duct segment objects
  std::vector<std::tuple<int,int,int>> seg_data; // vertex_in, vertex_out, seg_type for each duct segment
};

#endif /* CMINIGLANDDUCT_H_ */
