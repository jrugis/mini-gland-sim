/*
 * cCMesh.hpp
 *
 *  Created on: 26/07/2021
 *      Author: jrugis
 */

#ifndef CCMESH_H_
#define CCMESH_H_

#include "global_defs.hpp"

class cCMesh {
public:
  cCMesh(int cell_number, std::ofstream& out);
  ~cCMesh();
  int nverts, nfaces, ntets; // number of vertices, faces, tets for this cell
  MatrixN3d verts;
  MatrixN3i faces;
  MatrixN4i tets;
  MatrixN3d face_centers; //      centers
  MatrixN1d face_areas;   //      areas
  MatrixN1i face_types;   // face types for this cell: 0=apical, 1=basolateral, 3=basal

protected:

private:

};

#endif /* CCMESH_H_ */
