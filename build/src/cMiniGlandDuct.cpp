/*
 * cMiniGlandDuct.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

//#include <iomanip>
#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"
#include "cMiniGlandDuct.hpp"

cMiniGlandDuct::cMiniGlandDuct()
{
  id = "d1";       // there's only a single mini-gland duct
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);
  utils::get_parameters(PARAMETER_FILE_NAME, p, out); // NOTE: all the parameters are in this file

  // TO DO NEXT...
  // populate the segments vector


}

cMiniGlandDuct::~cMiniGlandDuct()
{
  out.close();
}

void cMiniGlandDuct::run()
{
/*  double t = 0.0;
  double solver_dt = p.at("delT");
  double error;
  struct timespec start, end;
  double elapsed;

  // simulation time stepping and synchronization
  clock_gettime(CLOCK_REALTIME, &start);
  while ((p.at("totalT") - t) > 0.000001) { // HARD CODED: assumes solver_dt always > 1us
    error = snd_recv(t, solver_dt);         // invoke the calcium solver
    if (error != 0.0) {
      // ...   change time step???
    }
    // invoke lumen step ???

    clock_gettime(CLOCK_REALTIME, &end);
    elapsed = (end.tv_sec - start.tv_sec) + ((end.tv_nsec - start.tv_nsec) / 1000000000.0);
    out << std::fixed << std::setprecision(3);
    out << "<Acinus> step duration: " << elapsed << "s" << std::endl;
    start = end;
    t += solver_dt;
  }

  // instruct cells to finish
  snd_recv(t, 0.0);
	*/
}
