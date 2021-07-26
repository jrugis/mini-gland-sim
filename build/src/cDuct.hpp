/*
 * cDuct.hpp
 *
 *  Created on: 21/4/2021
 *      Author: jrugis
 */

#ifndef CDUCT_H_
#define CDUCT_H_

#include <fstream>
#include <string>
#include <vector>

#include "inih/cpp/INIReader.h"
#include "global_defs.hpp"

class cMiniGland;
class cSCell;

class cDuct {
  friend cSCell; // child object

public:
  cDuct(cMiniGland* parent);
  ~cDuct();
  void step(double t, double dt);

protected:
  std::string id;
  std::ofstream out;    // runtime diagnostic file for this object
  INIReader* p;         // model parameters
  int nlsegs;           // number of duct lumen segments
  int ilsegs;           // index of first lumen segment

private:
  cMiniGland* parent;
  std::vector<cSCell*> scells;    // the striated cells associated with this duct
};

#endif /* CDUCT_H_ */
