/*
 * cDuctSegmentStriated.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

//#include <thread>

#include "cCell.hpp"
#include "cDuctSegment.hpp"
#include "cDuctSegmentStriated.hpp"

cDuctSegmentStriated::cDuctSegmentStriated(cMiniGlandDuct* parent, int seg_number) : cDuctSegment(parent, seg_number) {}

void cDuctSegmentStriated::step()
{
  // combine cells fluid flow  --  TO DO
  // ....

  out << "<DuctSegmentStriated> step - threads in use: " << omp_get_num_threads() << std::endl;
}

