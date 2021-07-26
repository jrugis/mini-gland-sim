/*
 * cDuct.cpp
 *
 *  Created on: 26/07/2021
 *      Author: jrugis
 */

#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cMiniGland.hpp"
#include "cLTree.hpp"
#include "cSCell.hpp"
#include "cDuct.hpp"

cDuct::cDuct(cMiniGland* _parent) : parent(_parent)
{
  id = parent->id + "d1";
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);
  p = parent->p; // the parameters INIReader object

  // get the duct lumen tree data
  std::ifstream mesh_file;
  utils::mesh_open(mesh_file, out);                                      // open the mesh file
  int nacinii = utils::mesh_get_count(mesh_file, std::string("acinii")); // number of acinii
  utils::mesh_end_header(mesh_file);                                     // skip over the rest of the header

  utils::mesh_skip_lines(mesh_file, nacinii);                            // skip to the duct data
  std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
  int nicells = std::stoi(tokens[0]);                                     // number of intercalated cell                                   
  int iicells = std::stoi(tokens[1]);                                     // first intercalated cell index
  int nscells = std::stoi(tokens[2]);                                     // number of striated cell                                   
  int iscells = std::stoi(tokens[3]);                                     // first striated cell index
  nlsegs = std::stoi(tokens[4]);                                      // number of duct lumen segments                                   
  ilsegs = std::stoi(tokens[5]);                                      // first duct lumen segment index
  mesh_file.close();

  out << "<Duct> Lumen segment count: " << nlsegs << std::endl; 
  out << "<Duct> Intercalated cell count: " << nicells << std::endl;
  out << "<Duct> Striated cell count: " << nscells << std::endl;

  // create the cells
  for (int i = 0; i < nscells; i++) {
	scells.push_back(new cSCell(this, i + iscells));
  }
}

cDuct::~cDuct()
{
  for (unsigned int i = 0; i < scells.size(); i++) delete scells[i]; // delete the striated cells
  out.close();
}

void cDuct::step(double t, double dt)
{
//********************************
  // TEMPORARY EXAMPLE: to get first duct lumen tree segment endpoints...
  cLTree *lt = parent->ltree;
  Vector3d in = lt->nodes.row(lt->segs(ilsegs + 0, 0));	
  Vector3d out = lt->nodes.row(lt->segs(ilsegs + 0, 1));	
//********************************

 // Much more goes here.,.
	
	 	
}