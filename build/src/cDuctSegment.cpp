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

  utils::get_parameters(PARAMETER_FILE_NAME, p, out); // NOTE: all the parameters are in this file
  create_cells(get_segment_data());
  out << "<DuctSegment> Cell count: " << cells.size() << std::endl;
}

cDuctSegment::~cDuctSegment() { 
  for(int i = 0; i<cells.size(); i++) delete cells[i]; // delete the cells
  out.close(); 
}

std::tuple<int,int> cDuctSegment::get_segment_data(){
  std::string line;                    // file line buffer
  std::vector<std::string> tokens;     // tokenized line
  std::vector<Vector3d> verts;

  std::ifstream mesh_file(MESH_FILE_NAME); // open the mesh file
  if (not mesh_file.is_open()) { utils::fatal_error("cell file " + std::string(MESH_FILE_NAME) + " could not be opened", out); }

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

  mesh_file.close();

  // return tuple of cell count and start index
  return(std::make_tuple( atoi(tokens[5].c_str()), atoi(tokens[6].c_str())));
}

// TO DO 
// create cells of same type as segment given cell count and start index
void cDuctSegment::create_cells(std::tuple<int,int> cell_info){}

// TO DO change "run" to "step", in cells too...
void cDuctSegment::run()
{
}
