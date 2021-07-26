/*
 * cMiniGland.hpp
 *
 *  Created on: 21/4/2021
 *      Author: jrugis
 */

#ifndef CMINIGLAND_H_
#define CMINIGLAND_H_

#include <fstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "inih/cpp/INIReader.h"
#include "global_defs.hpp"

class cDuct;
class cLTree;

class cMiniGland {
  friend cDuct;  // child object
  friend cLTree; //

public:
  cMiniGland();
  ~cMiniGland();
  void run();

protected:
  cLTree *ltree;                  // a single lumen tree object
  
private:
  std::string id;
  std::ofstream out;        // the runtime diagnostic file for this object
  INIReader *p;
  
  //std::vector<cAcinus*> acinii; // the acinii objects
  cDuct *duct;                    // a single duct object
};

#endif /* CMINIGLAND_H_ */
