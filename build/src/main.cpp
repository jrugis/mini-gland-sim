/*
 * main.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 *
 *  NOTE:
 *    environment variables:
 *      OPM_NUM_THREADS = 30
 *      OMP_NESTED = TRUE 
 */

#include <iostream>
#include <omp.h>
#include <time.h>
#include <unistd.h>

#include "cMiniGlandDuct.hpp"

#define TEMP_SIZE 40

// the main program function
int main(int argc, char** args)
{
  std::string host_name;
  struct timespec start, end;
  int duration;

  clock_gettime(CLOCK_REALTIME, &start);

  // get the hostname
  char temp[TEMP_SIZE];
  gethostname(temp, TEMP_SIZE);
  host_name = temp;
  std::cout << "<main> running on host: " << host_name << std::endl;
  std::cout << "<main> maximum threads: " << omp_get_max_threads() << std::endl;
  
  // run the simulation
  cMiniGlandDuct* duct = new cMiniGlandDuct();
  duct->run();
  delete duct;

  // output total running time
  clock_gettime(CLOCK_REALTIME, &end);
  duration = end.tv_sec - start.tv_sec;
  std::cout << "<main> execution time: " << duration << "s" << std::endl;

  return 0;
}
