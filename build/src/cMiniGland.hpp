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

class cAcinus;
class cDuct;
class cLTree;

class cMiniGland {
    friend cAcinus;  // child object
    friend cDuct;    //
    friend cLTree;   //

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
  
  std::vector<cAcinus*> acinii; // the acinus objects
  cDuct *duct;                    // a single duct object
};

#endif /* CMINIGLAND_H_ */
