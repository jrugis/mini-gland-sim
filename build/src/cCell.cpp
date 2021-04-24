/*
 * cCell.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/tokenizer.hpp>

//#include <iomanip>
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

  // get the cell data summary over all cells
  std::string line;                    // file line buffer
  std::vector<std::string> tokens;     // tokenized line
  // open the mesh file
  std::ifstream mesh_file(MESH_FILE_NAME); 
  if (not mesh_file.is_open()) { utils::fatal_error("cell file " + std::string(MESH_FILE_NAME) + " could not be opened", out); }
  // get the total mesh vertex count
  while (getline(mesh_file, line)) {
  	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  	if (tokens[1]=="vertex" ) break;
  }
  int tnverts = atoi(tokens[2].c_str());
  // get the total face count
  while (getline(mesh_file, line)) {
  	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  	if (tokens[1]=="face" ) break;
  }
  int tnfaces = atoi(tokens[2].c_str());
  // get the total tetrahedron count
  while (getline(mesh_file, line)) {
  	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  	if (tokens[1]=="tetrahedron" ) break;
  }
  int tntets = atoi(tokens[2].c_str());
  while (getline(mesh_file, line)) {
  	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  	if (tokens[1]=="cell" ) break;
  }
  int tncells = atoi(tokens[2].c_str());
  // skip over the rest of the header
  while (getline(mesh_file, line)) {
  	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  	if (tokens[0]=="end_header" ) break;
  }

  // TO DO next.................
  // skip to this cell
  // get n&i for verts, faces, tets
  // get vertices for this cell
  // get faces for this cell
  // get tetrahedra for this cell
  
  mesh_file.close();
  
}

cCell::~cCell() { out.close(); }

void cCell::step()
{
}
