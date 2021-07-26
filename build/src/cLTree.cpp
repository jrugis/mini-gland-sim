/*
 * cLTree.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cMiniGland.hpp"
#include "cLTree.hpp"

cLTree::cLTree(cMiniGland* _parent) : parent(_parent)
{
  id = parent->id + "lt1";  // there's only a single lumen tree
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);

  // get the lumen tree data
  std::ifstream mesh_file;
  utils::mesh_open(mesh_file, out);                                           // open the mesh file
  int nacinii = utils::mesh_get_count(mesh_file, std::string("acinii"));      // get acinii count
  int nduct = utils::mesh_get_count(mesh_file, std::string("duct"));          //     duct 
  int nnodes = utils::mesh_get_count(mesh_file, std::string("lumen_node"));   //   lumen node
  int nsegs = utils::mesh_get_count(mesh_file, std::string("lumen_segment")); //   lumen segment
  utils::mesh_end_header(mesh_file);                                          // skip over the rest of the header
  
  utils::mesh_skip_lines(mesh_file, nacinii + nduct);  // skip to the lumen data
  nodes.resize(nnodes, Eigen::NoChange);               // resize the nodes matrix
  for (int i = 0; i < nnodes; i++) {                   // get nodes
	std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
	nodes.row(i) = Vector3d(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]));
  }
  segs.resize(nsegs, Eigen::NoChange);              // resize the segments matrix
  diams.resize(nsegs, Eigen::NoChange);             // resize the diameters matrix
  for (int i = 0; i < nsegs; i++) {                 // get nodes
    std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
	segs.row(i) = Vector2i(std::stoi(tokens[0]), std::stoi(tokens[1]));
    diams(i) = std::stod(tokens[2]);
  }
  mesh_file.close();
  out << "<LTree> Lumen tree node count: " << nnodes << std::endl; 
  out << "<LTree> Lumen tree segment count: " << nsegs << std::endl; 


  // get the duct segment mesh data
  //std::ifstream mesh_file;
  //utils::mesh_open(mesh_file, out);                                        // open the mesh file
  //int nnodes = utils::mesh_get_count(mesh_file, std::string("duct_node")); // get the duct node count
  //utils::mesh_end_header(mesh_file);                                       // skip over the rest of the header

  // get the duct node data
  //std::vector<Vector3d> verts;
  //for (int i = 0; i < nnodes; i++) {
  //  std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
  //  verts.push_back(Vector3d(std::stof(tokens[0]), std::stof(tokens[1]), std::stof(tokens[2])));
  //}

  // skip to the correct segment and get the segment data
  //utils::mesh_skip_lines(mesh_file, seg_number); // skip to the correct segment
  //std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
  //vertex_in = verts[std::stoi(tokens[0])];
  //vertex_out = verts[std::stoi(tokens[1])];
  //inner_diameter = std::stof(tokens[2]);
  //outer_diameter = std::stof(tokens[3]);
  //seg_type = std::stoi(tokens[4]);
  //int ncells = std::stoi(tokens[5]);
  //int icells = std::stoi(tokens[6]);
  //mesh_file.close();

  // create the cells
  //out << "<DuctSegment> Cell count: " << ncells << std::endl;
  //for (int i = 0; i < ncells; i++) {
  //  if (seg_type == ACINUS)
  //    cells.push_back(new cCellAcinus(this, i + icells));
  //  else if (seg_type == INTERCALATED)
  //    cells.push_back(new cCellIntercalated(this, i + icells));
  //  else if (seg_type == STRIATED)
  //    cells.push_back(new cCellStriated(this, i + icells));
  //}
  
}

cLTree::~cLTree()
{
  out.close();
}
