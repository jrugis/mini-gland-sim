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

  out << "<LTree> Lumen tree node count: " << nnodes << std::endl; 
  out << "<LTree> Lumen tree segment count: " << nsegs << std::endl; 
  utils::mesh_skip_lines(mesh_file, nacinii + nduct);  // skip to the lumen data

  nodes.resize(nnodes, Eigen::NoChange);               // resize the nodes matrix
  radii.resize(nnodes, Eigen::NoChange);                // resize the diameters matrix
  for (int i = 0; i < nnodes; i++) {                   // get nodes
	std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
	nodes.row(i) = Vector3d(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]));
    radii(i) = std::stod(tokens[3]);
  }

  segs.resize(nsegs, Eigen::NoChange);              // resize the segments matrix
  for (int i = 0; i < nsegs; i++) {                 // get nodes
    std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
	segs.row(i) = Vector2i(std::stoi(tokens[0]), std::stoi(tokens[1]));
  }
  mesh_file.close();
}

cLTree::~cLTree()
{
  out.close();
}
