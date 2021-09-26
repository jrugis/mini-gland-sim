/*
 * cSIMesh.hpp
 *
 *  Created on: 26/09/2021
 *      Author: jrugis
 */

#ifndef CSIMESH_H_
#define CSIMESH_H_

#include "global_defs.hpp"

class cSIMesh {
public:
  cSIMesh(std::string fname, std::ofstream& out);
  ~cSIMesh();
  int nverts, nfaces; // number of vertices, faces for this cell
  MatrixN3d verts;
  MatrixN3i faces;
  MatrixN1i face_types;      // face types: 0=apical, 1=basolateral, 3=basal
  MatrixN1d face_areas;      //      areas
  MatrixN1i duct_idx;        //      index of nearest duct segment
  MatrixN1d dist_from_duct;  //      distance from duct
  MatrixN1d dist_along_duct; //      distance along duct
};

#endif /* CSIMESH_H_ */
