/*
 * cCell.hpp
 *
 *  Created on: 21/4/2021
 *      Author: jrugis
 */

#ifndef CCELL_H_
#define CCELL_H_

#include <fstream>
#include <string>
#include <unordered_map>

#include "global_defs.hpp"

class cDuctSegment;
// class cCellAcinus;
// class cCellIntercalated;
// class cCellStriated;

class cCell {
  // friend cCellAcinus;        // derived classes
  // friend cCellIntercalated;
  // friend cCellStriated;

  public:
  cCell(cDuctSegment* parent, int cell_number);
  virtual ~cCell();
  virtual void step(){}; // defined in dervied class

  protected:
  cDuctSegment* parent;
  std::string id;
  std::ofstream out;                         // runtime diagnostic file for this object
  std::unordered_map<std::string, double> p; // the model parameters

  int cell_number;           // this cell number
  int nverts, nfaces, ntets; // number of vertices, faces, tets for this cell
  MatrixN3d verts;
  MatrixN3i faces;
  MatrixN4i tets;
  MatrixN1i face_types;   // face types for this cell: 0=apical, 1=basolateral, 3=basal
  MatrixN3d face_centers; //      centers
  MatrixN1d face_areas;   //      areas
};

class cCellAcinus : public cCell {
  public:
  cCellAcinus(cDuctSegment* parent, int cell_number) : cCell(parent, cell_number){};
  virtual void step();
};

class cCellIntercalated : public cCell {
  public:
  cCellIntercalated(cDuctSegment* parent, int cell_number) : cCell(parent, cell_number){};
  virtual void step();
};

#endif /* CCELL_H_ */
