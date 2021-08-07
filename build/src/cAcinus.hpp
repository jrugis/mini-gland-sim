/*
 * cAcinus.hpp
 *
 *  Created on: 2/8/2021
 *      Author: jrugis
 */

#ifndef CACINUS_H_
#define CACINUS_H_

#include <fstream>
#include <string>
#include <vector>

#include "inih/cpp/INIReader.h"
#include "global_defs.hpp"

class cMiniGland;
class cACell;

class cAcinus {
  friend cACell; // child object

public:
  cAcinus(cMiniGland* parent, int _acinus_number);
  ~cAcinus();
  void step(double t, double dt);

protected:
  std::string id;
  std::ofstream out;    // runtime diagnostic file for this object
  INIReader* p;         // model parameters
  int nlsegs;           // number of duct lumen segments
  int ilsegs;           // index of first lumen segment

private:
  cMiniGland* parent;
  int acinus_number;     // this acinus number
  std::vector<cACell*> cells;    // the cells associated with this acinus
};

#endif /* CACINUS_H_ */
