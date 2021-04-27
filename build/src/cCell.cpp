/*
 * cCell.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cDuctSegment.hpp"
#include "cCell.hpp"

cCell::cCell(cDuctSegment* _parent, int _cell_number) : parent(_parent), cell_number(_cell_number)
{
  id = parent->id + "c" + std::to_string(cell_number+1); // NOTE: one based cell id
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);
  p = parent->p; // the parameters map

  // get the duct cell mesh data
  std::ifstream mesh_file;
  utils::mesh_open(mesh_file, out);                                           // open the mesh file
  int tnverts = utils::mesh_get_count(mesh_file, std::string("vertex"));      // get the total mesh vertex count
  int tnfaces = utils::mesh_get_count(mesh_file, std::string("face"));        // get the total mesh face count
  int tntets = utils::mesh_get_count(mesh_file, std::string("tetrahedron"));  // get the total mesh tetrahedron count
  int tncells = utils::mesh_get_count(mesh_file, std::string("cell"));        // get the total mesh cell count
  utils::mesh_end_header(mesh_file);                                          // skip over the rest of the header

  // TO DO next.................
  // skip to this cell
  // get n&i for verts, faces, tets
  // get vertices for this cell
  // get faces for this cell
  // get tetrahedra for this cell
  
  mesh_file.close();
  
}

cCell::~cCell() { out.close(); }
