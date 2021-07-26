/*
 * cSCell.hpp
 *
 *  Created on: 21/4/2021
 *      Author: jrugis
 */

#ifndef CSCELL_H_
#define CSCELL_H_

#include <fstream>
#include <string>

#include "inih/cpp/INIReader.h"
#include "global_defs.hpp"

class cDuct;
class cCMesh;

class cSCell {

public:
  cSCell(cDuct* _parent, int _cell_number);
  ~cSCell();
  void step();

protected:
	  
private:
  cDuct* parent;
  std::string id;
  std::ofstream out;   // runtime diagnostic file for this object
  INIReader *p;        // model parameters

  int cell_number;     // this cell number
  cCMesh *mesh;        // this cell mesh
};


#endif /* CSCELL_H_ */
