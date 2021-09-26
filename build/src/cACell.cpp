/*
 * cACell.cpp
 *
 *  Created on: 2/8/2021
 *      Author: jrugis
 */

#include <cmath>
#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cAcinus.hpp"
//#include "cCMesh.hpp"
#include "cACell.hpp"

cACell::cACell(cAcinus* _parent, int _cell_number) : parent(_parent), cell_number(_cell_number)
{
  id = parent->id + "c" + std::to_string(cell_number + 1); // NOTE: one based cell id
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);
  p = parent->p; // the parameters map
  //mesh = new cCMesh(cell_number, out);

  make_matrices(); // create the constant matrices

}

cACell::~cACell() {
  //delete mesh;
  out.close();
}

void cACell::make_matrices() {
	
}
