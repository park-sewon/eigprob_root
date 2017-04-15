#include "../include/linsys.h"
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
int first_run = 1;
int iRRAM_compute(const int & argc, char ** const & argv)
{

  int sep = 1;
  int model = 0;
  int d = atoi(argv[1]);
  int k = atoi(argv[2]);
  int mode = atoi(argv[3]);

  if(first_run == 1)
    std::cout <<"the process (algorithm: "<< mode <<")runs infinite re-iteration if there exists no error. The pid is "<< ::getpid()<<"\n";
  first_run = 0;


  REALMATRIX m = REALMATRIX(d,d);
  REALMATRIX O = REALMATRIX(d,d);
	std::vector<std::pair< REALMATRIX, REAL> >  Q;


	int t = 0;

// building diag matrix
  for(int i = 0; i<d;i++)
    m(i,i) = REAL(argv[i+4]);

// building seed matrix for generating orthogonal matrix in Haar measure
	for(int i = 0; i < d; i++)
  {
 		for(int j = 0; j < d; j++)
		{
  		O(i,j) = REAL(argv[t+4+d]);
  		t += 1;
		}
	}

	O = QRDecomposition(O);
	m = O*m*transpose(O);

// m is real symmetrix matrix

	std::clock_t begin_time = std::clock();
	if (mode == 1)
	{
  		Q = diagonalizeNaiveEig(m,k);
	}
	else if (mode == 2)
	{
  		Q = diagonalizeTrisectionEig(m,k);
	}
	else if (mode == 3)
	{
  		Q = diagonalizeNewtonEig(m,k);	
  }
  else if (mode == 4)
  {
      Q = diagonalizeCombinedEig(m,k);        
  }
	float duration = float( std::clock () - begin_time ) /  CLOCKS_PER_SEC;

  REAL max = 0;
  sizetype error;   
  int max_errorexp; 

  REAL eigen;
  REALMATRIX R, S;
  for(int i=0; i<(int)Q.size(); i++)
  {
    R = Q[i].first;
    for(unsigned int j =0; j<R.maxrow; j++)
      for(unsigned int k = 0; k<R.maxcolumn; k++)
        R(j,k) = R(j,k) * Q[i].second;
    S = m*Q[i].first - R;
    S.geterror(error);
    if(i==0)
      max_errorexp =  (int)std::floor(std::log2(error.mantissa)) + error.exponent;
    else
      max_errorexp = std::max(max_errorexp, (int)std::floor(std::log2(error.mantissa)) + error.exponent);
    max = maximum(maxnorm(S), max);
  }
  std::cout << "Testing validity of mode "<< mode << " upto precision " << max_errorexp << "\n";
  if (max > 0 )
  {
    std::cout << "invalid computation.. aborting...\n";
    return 1;
  }
  return 0;
}

int main(int argc, char ** argv)
{
	iRRAM_initialize(argc, argv);
	return iRRAM::exec(iRRAM_compute, argc, argv);
}
