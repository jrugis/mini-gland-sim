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

class cDuctSegment;

class cCell {
  public:
  cCell(cDuctSegment* parent, int cell_number);
  ~cCell();
  void run();

  private:
  cDuctSegment* parent;
  std::string id;
  std::ofstream out;                         // runtime diagnostic file for this object
  std::unordered_map<std::string, double> p; // the model parameters

  int cell_number;           // this cell number
  int nverts, nfaces, ntets; // number of vertices, faces, tets for this cell
  MatrixN3d verts;           
  MatrixN3i faces;
  MatrixN4i tets;

  void get_cell_data();
};

#endif /* CCELL_H_ */
