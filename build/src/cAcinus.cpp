/*
 * cAcinus.cpp
 *
 *  Created on: 2/8/2021
 *      Author: jrugis
 */

#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cMiniGland.hpp"
#include "cLTree.hpp"
#include "cACell.hpp"
#include "cAcinus.hpp"

cAcinus::cAcinus(cMiniGland* _parent, int _acinus_number) : parent(_parent), acinus_number(_acinus_number)
{
  id = parent->id + "a" + std::to_string(acinus_number + 1); // NOTE: one based acinus id;
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);
  p = parent->p; // the parameters INIReader object
  /*
  // get the acinus data
  std::ifstream mesh_file;
  utils::mesh_open(mesh_file, out);                                      // open the mesh file
  //int nacinii = utils::mesh_get_count(mesh_file, std::string("acinii")); // number of acinii
  utils::mesh_end_header(mesh_file);                                     // skip over the rest of the header

  utils::mesh_skip_lines(mesh_file, acinus_number-1);                    // skip to the correct acinus
  std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
  int ncells = std::stoi(tokens[0]);                                     // number of cells                                   
  int icells = std::stoi(tokens[1]);                                     // first cell index
  nlsegs = std::stoi(tokens[2]);                                      // number of duct lumen segments                                   
  ilsegs = std::stoi(tokens[3]);                                      // first duct lumen segment index
  mesh_file.close();

  out << "<Acinus> Lumen segment count: " << nlsegs << std::endl; 
  out << "<Acinus> Cell count: " << ncells << std::endl;

  // create the cells
  for (int i = 0; i < ncells; i++) {
	cells.push_back(new cACell(this, i + icells));
  }
*/
}

cAcinus::~cAcinus()
{
  //for (unsigned int i = 0; i < cells.size(); i++) delete cells[i]; // delete the striated cells
  out.close();
}

void cAcinus::step(double t, double dt)
{
//********************************
  // TEMPORARY EXAMPLE: to get first duct lumen tree segment endpoints...
  //cLTree *lt = parent->ltree;
  //Vector3d in = lt->nodes.row(lt->segs(ilsegs + 0, 0));	
  //Vector3d out = lt->nodes.row(lt->segs(ilsegs + 0, 1));	
//********************************

 // Much more goes here.,.
	
	 	
}