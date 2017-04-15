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
	
sizetype geterror( const std::vector<std::pair< REALMATRIX, REAL> > & Q)
{
  sizetype error, lerror;
  Q[0].first.geterror(error);
  for (unsigned int i=0; i< Q.size(); i++)
  {
    lerror = geterror(Q[i].first);
    sizetype_max(error,error,lerror);
  }
  return error;
}


int iRRAM_compute(const int & argc, char ** const & argv)
{

  int sep = atoi(argv[1]);

  int mode = atoi(argv[2]);
  int d = 100;
  REALMATRIX m = REALMATRIX(d,d);
  std::vector<std::pair< REALMATRIX, REAL> >  Q;
  REALMATRIX O;
  int k = d;
  int t = 0;


  O = REALMATRIX(d,d);

  for (int i = 0; i < d; i ++)
    m(i,i) = power(sep,-i);
  
  for(int i = 0; i < d; i++)
  {
    for(int j = 0; j < d; j++)
    {
      O(i,j) = REAL(argv[t+3]);
      t += 1;
    }
  }
  O = QRDecomposition(O);
  m = O*m*transpose(O);



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
  else if (mode == 5)
  {
      Q = diagonalizeBoundedEig(m,k);        
  }

  float duration = float( std::clock () - begin_time ) /  CLOCKS_PER_SEC;


  // printing statistics of the iteration. 
  sizetype lerror = geterror(Q);
  int errorexp = (int) std::floor(std::log2(lerror.mantissa)) + lerror.exponent;

  std::cout <<"["<<::getpid()<<"]"<< ", "<< sep << ", " << mode <<", "<<std::floor(std::log2(lerror.mantissa)) + lerror.exponent << ", " 
  << duration <<", "<< TIME_CHAR_POLYNOMIAL << ", " << TIME_ROOT_FINDING << ", " <<TIME_ROOT_SEPARATION <<", "<< TIME_ROOT_APPROXIMATION << ", " <<TIME_GAUSSIAN_ELIMINATION 
  << ", " <<  state.ACTUAL_STACK.actual_prec << ", " 
  << MAXDEPTH << ", " << MAX_TRISECTION << ", " << MAX_NEWTON << ", "<< NewTonTimesAll << ", " << NewTonTimes << ", " << NewTon  << ", " 
  << float( std::clock ()) /  CLOCKS_PER_SEC  <<", "<<swrite(SEPARATION,20) << "\n" << std::flush;

  if (errorexp < - 10000)
    return 0; 


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




	// int model = atoi(argv[1]);	
	// int multi =  atoi(argv[2]);
	// int mode = atoi(argv[3]);

	// REALMATRIX m = REALMATRIX(10 * multi,10 * multi);
	// REALMATRIX Q;
 //  REALMATRIX O;
 //  int d = 10 * multi;
	// int k;
	// int t = 0;



 // 	if (model == 1)
 // 	{
 //    k = 10;
	// 	O = REALMATRIX(10 * multi,10 * multi);

 //    for (int j = 0; j < multi; j ++)
 //  	 for (int i = 0; i < 10; i++)
 //        m(i+10*j,i+10*j) = power(1000,-i);
    
	// 	for(int i = 0; i < 10 * multi; i++)
 //    {
 //   		for(int j = 0; j < 10 * multi; j++)
 //  		{
 //    		O(i,j) = REAL(argv[t+4]);
 //    		t += 1;
 //  		}
 //  	}

 //  	O = QRDecomposition(O);
 //  	m = O*m*transpose(O);
	// }



	// std::clock_t begin_time = std::clock();
	// if (mode == 1)
	// {
 //  		Q = diagonalizeNaive(m,k);
	// }
	// else if (mode == 2)
	// {
 //  		Q = diagonalizeTrisection(m,k);
	// }
	// else if (mode == 3)
	// {
 //  		Q = diagonalizeNewton(m,k);	
 //  }
 //  else if (mode == 4)
 //  {
 //      Q = diagonalizeCombined(m,k);        
 //  }
	// float duration = float( std::clock () - begin_time ) /  CLOCKS_PER_SEC;


 //  // printing statistics of the iteration. 
 //  sizetype error;   Q.geterror(error);
 //  int errorexp = std::floor(std::log2(error.mantissa)) + error.exponent;

 //  std::cout <<"["<<::getpid()<<"]"<< ", "<< model << ", " << d << ", " << mode <<", "<<std::floor(std::log2(error.mantissa)) + error.exponent << ", " 
 //  << duration <<", "<< TIME_CHAR_POLYNOMIAL << ", " << TIME_ROOT_FINDING << ", " <<TIME_ROOT_SEPARATION <<", "<< TIME_ROOT_APPROXIMATION << ", " <<TIME_GAUSSIAN_ELIMINATION 
 //  << ", " <<  state.ACTUAL_STACK.actual_prec << ", " 
 //  << MAXDEPTH << ", " << MAX_TRISECTION << ", " << MAX_NEWTON << ", "<< NewTonTimesAll << ", " << NewTonTimes << ", " << NewTon  << ", " 
 //  << float( std::clock ()) /  CLOCKS_PER_SEC  <<", "<<swrite(SEPARATION,20) << "\n" << std::flush;

 //    // infinite iteration
 //    REAL x, y;
 //    x =1; y=1;
 //    if(x >y)
 //        return 1;
 //    return 0;
}

int main(int argc, char ** argv)
{
	iRRAM_initialize(argc, argv);
	return iRRAM::exec(iRRAM_compute, argc, argv);
}
