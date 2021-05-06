/*
 * cDuctSegmentStriated.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

//#include <thread>

#include "cCell.hpp"
#include "cDuctSegment.hpp"

cDuctSegmentStriated::cDuctSegmentStriated(cMiniGlandDuct* parent, int seg_number) : cDuctSegment(parent, seg_number) {}

void cDuctSegmentStriated::step()
{
#pragma omp parallel for
  for (auto cell : cells) { cell->step(); }

  // combine cells fluid flow  --  TO DO
  // ....

  // NOTE - do this to the cells vector in the derived class access the the cell data
  //  Base* basepointer = new Derived;
  //  static_cast<Derived*>(basepointer)->derived_int; // Can now access

  out << "<DuctSegmentStriated> step - threads in use: " << omp_get_num_threads() << std::endl;
}

