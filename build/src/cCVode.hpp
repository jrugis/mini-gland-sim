#ifndef _CCVODE_H_
#define _CCVODE_H_

#include <fstream>
#include <string>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */

#include "global_defs.hpp"

class cCVode {
  public:
  cCVode(std::ofstream& out_, realtype abstol_, realtype reltol_);
  ~cCVode();
  void init(CVRhsFn rhs_func, Array1Nd& yini, void* user_data);
  void run(realtype t, realtype tend, Array1Nd& yout);
  void PrintFinalStatsBrief();
  void PrintFinalStatsDetailed();

  private:
  std::ofstream& out;
  bool initialised;
  sunindextype nvars;
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
  void* cvode_mem;
  realtype abstol, reltol;

  void check_retval(void* returnvalue, std::string funcname, int opt);
};

#endif
