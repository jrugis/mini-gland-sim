/*
 * cCellStriated.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#include <chrono> // std::chrono::seconds
#include <thread> // std::this_thread::sleep_for

#include "cCell.hpp"

cCellStriated::cCellStriated(cDuctSegment* parent, int cell_number) : cCell(parent, cell_number)
{
  out << "<CellStriated> @constructor" << std::endl;
}

void cCellStriated::step()
{
  // TO DO
  // ...
  std::this_thread::sleep_for(std::chrono::milliseconds(200));

  // out << "<CellStriated> threads in use: " << omp_get_num_threads() << std::endl;
  out << "<CellStriated> step" << std::endl;
}
