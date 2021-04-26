/*
 * cCellStriated.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds

#include "cCell.hpp"

void cCellStriated::step(){
  // TO DO
  // ...
  std::this_thread::sleep_for (std::chrono::milliseconds(20));

  //out << "<CellStriated> threads in use: " << omp_get_num_threads() << std::endl;
  out << "<CellStriated> step" << std::endl;
}
