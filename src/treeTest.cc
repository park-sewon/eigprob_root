#include "../include/linsys.h"
#include "../include/intervals.h"
#include "../include/polynomial.h"
#include "../include/polynomialroot.h"
#include "../include/utilities.h"

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include <ctime>
#include <unistd.h>

#include <math.h>
#include <functional>
#include <iRRAM/limit_templates.h>
// #include <tuple>
using namespace iRRAM;
	
int iRRAM_compute(const int & argc, char ** const & argv)
{

  int sep = atoi(argv[1]);
  int mode = atoi(argv[2]);
  int d = atoi(argv[3]);
  POLYNOMIAL P = Wilkinson(d, sep);
  print(P);  
  std::vector<std::pair<RATIONALINTERVAL, int> > Q;


  std::clock_t begin_time = std::clock();
  if (mode == 4)
  {
      Q = RealDistinctRootSeparationNaive(P,d);        
  }
  else if (mode == 5)
  {
      Q = RealDistinctRootSeparationBounded(P,d);        
  }
  float duration = float( std::clock () - begin_time ) /  CLOCKS_PER_SEC;


  std::cout <<"["<<::getpid()<<"]"<< ", "<< sep << ", " << mode <<", " << MAXDEPTH <<", " << duration <<"\n" << std::flush<<"\n\n";

  for(int i=0; i<d; i++)
  {
    cout<<Q[i].second << " in (" << REAL(Q[i].first.center - Q[i].first.radius) <<", " <<REAL(Q[i].first.center + Q[i].first.radius) <<")\n";
  }


  return 0;
}

int main(int argc, char ** argv)
{
	iRRAM_initialize(argc, argv);
	return iRRAM::exec(iRRAM_compute, argc, argv);
}
