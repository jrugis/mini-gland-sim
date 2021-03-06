/*
 * cLTree.cpp
 *
 *  Created on: 26/09/2021
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
  id = "_lumen_tree";  // there's only a single lumen tree
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);

  // get the lumen tree data
  std::ifstream mesh_file;
  utils::mesh_open(mesh_file, DUCT_FILE_NAME, out);                                           // open the mesh file
  int nnodes = utils::mesh_get_count(mesh_file, std::string("duct_node"));   // duct nodes count
  int nsegs = utils::mesh_get_count(mesh_file, std::string("duct_segment")); // duct segments count
  utils::mesh_end_header(mesh_file);                                         // skip over the rest of the header
  out << "<LTree> Lumen tree node count: " << nnodes << std::endl; 
  out << "<LTree> Lumen tree segment count: " << nsegs << std::endl; 

  nodes.resize(nnodes, Eigen::NoChange);     // resize the nodes matrix
  radii.resize(nnodes, Eigen::NoChange);     // resize the diameters matrix
  for (int i = 0; i < nnodes; i++) {         // get nodes
	std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
	nodes.row(i) = Vector3d(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]));
    radii(i) = std::stod(tokens[3]);
  }
  segs.resize(nsegs, Eigen::NoChange);       // resize the segments matrix
  types.resize(nsegs, Eigen::NoChange);        // resize the segment ids matrix
  for (int i = 0; i < nsegs; i++) {          // get segments
    std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
	segs.row(i) = Vector2i(std::stoi(tokens[0]), std::stoi(tokens[1]));
    types(i) = std::stoi(tokens[2]);
  }
  mesh_file.close();
}

cLTree::~cLTree()
{
  out.close();
}
