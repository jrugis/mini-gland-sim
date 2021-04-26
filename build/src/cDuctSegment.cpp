/*
 * cDuctSegment.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/tokenizer.hpp>

#include <iostream>
#include <string>
#include <tuple>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cCell.hpp"
#include "cMiniGlandDuct.hpp"
#include "cDuctSegment.hpp"

cDuctSegment::cDuctSegment(cMiniGlandDuct* _parent, int _seg_number) : parent(_parent), seg_number(_seg_number)
{
  id = parent->id + "s" + std::to_string(seg_number+1); // NOTE: one based segment id
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);
  p = parent->p; // the parameters map

  // get the duct segment data
  std::string line;                    // file line buffer
  std::vector<std::string> tokens;     // tokenized line
  std::vector<Vector3d> verts;
  // open the mesh file
  std::ifstream mesh_file(MESH_FILE_NAME); 
  if (not mesh_file.is_open()) { utils::fatal_error("mesh file " + std::string(MESH_FILE_NAME) + " could not be opened", out); }
  // get the duct node vertex count
  while (getline(mesh_file, line)) {
  	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  	if (tokens[1]=="duct_node" ) break;
  }
  int nnodes = atoi(tokens[2].c_str());
  // skip over the rest of the header
  while (getline(mesh_file, line)) {
  	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  	if (tokens[0]=="end_header" ) break;
  }
  // get the node data
  for(int i=0; i<nnodes; i++){
  	getline(mesh_file, line);
	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
	verts.push_back(Vector3d(atof(tokens[0].c_str()), atof(tokens[1].c_str()), atof(tokens[2].c_str())));
  }

  // skip to the correct segment and get the data
  for(int i=0; i<seg_number+1; i++) getline(mesh_file, line);
  boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  vertex_in = verts[atoi(tokens[0].c_str())];
  vertex_out = verts[atoi(tokens[1].c_str())];
  inner_diameter = atof(tokens[2].c_str());
  outer_diameter = atof(tokens[3].c_str());
  seg_type = atoi(tokens[4].c_str());
  int ncells = atoi(tokens[5].c_str());
  int icells = atoi(tokens[6].c_str());
  mesh_file.close();

  // create the cells
  out << "<DuctSegment> Cell count: " << ncells << std::endl;
  for(int i=0; i<ncells; i++){
	if(seg_type == ACINUS) cells.push_back(new cCellAcinus(this, i+icells));
	else if(seg_type == INTERCALATED) cells.push_back(new cCellIntercalated(this, i+icells));
	else if(seg_type == STRIATED) cells.push_back(new cCellStriated(this, i+icells));
  }
}

cDuctSegment::~cDuctSegment() { 
  for(unsigned int i = 0; i<cells.size(); i++) delete cells[i]; // delete the cells
  out.close(); 
}
