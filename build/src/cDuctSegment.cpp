/*
 * cDuctSegment.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cCell.hpp"
#include "cCellStriated.hpp"
#include "cDuctSegment.hpp"
#include "cMiniGlandDuct.hpp"

cDuctSegment::cDuctSegment(cMiniGlandDuct* _parent, int _seg_number) : parent(_parent), seg_number(_seg_number)
{
  id = parent->id + "s" + std::to_string(seg_number + 1); // NOTE: one based segment id
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);
  p = parent->p; // the parameters map

  // get the duct segment mesh data
  std::ifstream mesh_file;
  utils::mesh_open(mesh_file, out);                                        // open the mesh file
  int nnodes = utils::mesh_get_count(mesh_file, std::string("duct_node")); // get the duct node count
  utils::mesh_end_header(mesh_file);                                       // skip over the rest of the header

  // get the duct node data
  std::vector<Vector3d> verts;
  for (int i = 0; i < nnodes; i++) {
    std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
    verts.push_back(Vector3d(std::stof(tokens[0]), std::stof(tokens[1]), std::stof(tokens[2])));
  }

  // skip to the correct segment and get the segment data
  utils::mesh_skip_lines(mesh_file, seg_number); // skip to the correct segment
  std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
  vertex_in = verts[std::stoi(tokens[0])];
  vertex_out = verts[std::stoi(tokens[1])];
  inner_diameter = std::stof(tokens[2]);
  outer_diameter = std::stof(tokens[3]);
  seg_type = std::stoi(tokens[4]);
  int ncells = std::stoi(tokens[5]);
  int icells = std::stoi(tokens[6]);
  mesh_file.close();

  // create the cells
  out << "<DuctSegment> Cell count: " << ncells << std::endl;
  for (int i = 0; i < ncells; i++) {
    if (seg_type == ACINUS)
      cells.push_back(new cCellAcinus(this, i + icells));
    else if (seg_type == INTERCALATED)
      cells.push_back(new cCellIntercalated(this, i + icells));
    else if (seg_type == STRIATED)
      cells.push_back(new cCellStriated(this, i + icells));
  }
}

cDuctSegment::~cDuctSegment()
{
  for (unsigned int i = 0; i < cells.size(); i++) delete cells[i]; // delete the cells
  out.close();
}
