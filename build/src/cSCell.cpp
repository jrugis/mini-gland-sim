/*
 * cSCell.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#include <cmath>
#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cDuct.hpp"
#include "cCMesh.hpp"
#include "cSCell.hpp"

cSCell::cSCell(cDuct* _parent, int _cell_number) : parent(_parent), cell_number(_cell_number)
{
  id = parent->id + "c" + std::to_string(cell_number + 1); // NOTE: one based cell id
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);
  p = parent->p; // the parameters map
  mesh = new cCMesh(cell_number, out);
  
// ***********************************************  
// - TEMPORARY - to be replaced later  
  int a = 0;
  int bl = 0;
  int b = 0;
  for (int i = 0; i < mesh->nfaces; i++) {
	double x = mesh->face_centers(i,0);  
	double y = mesh->face_centers(i,1);  
	double dist = sqrt(pow(x,2) + pow(y,2));
	if ((dist - 4.0) < 1.5) { a++; mesh->face_types(i) = APICAL;}
    else if ((23.8 - dist) < 3.0) { b++; mesh->face_types(i) = BASAL;}
    else { bl++; mesh->face_types(i) = BASOLATERAL;}
  }
  out << "<SCell> Apical face count: " << a << std::endl;
  out << "<SCell> Basolateral face count: " << bl << std::endl;
  out << "<SCell> Basal face count: " << b << std::endl;
// ***********************************************  

}

cSCell::~cSCell() {
  delete mesh;
  out.close();
}
