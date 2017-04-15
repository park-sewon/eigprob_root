#ifndef RPOLYNOMIALROOT_H
#define RPOLYNOMIALROOT_H

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "polynomial.h"
#include "intervals.h"

#include <utility>
using namespace iRRAM;



//GLOBAL VRIABLES FOR MEASUREMENT PURPOSE:
extern float TIME_ROOT_SEPARATION;
extern int FLAG_ROOT_SEPARATION;
extern float TIME_ROOT_FINDING;
extern float TIME_ROOT_APPROXIMATION;
extern float TIME_CHAR_POLYNOMIAL;
extern float TIME_GAUSSIAN_ELIMINATION;
extern int MAXDEPTH;
extern int NewTonTimesAll;
extern int NewTonTimes;
extern int NewTon;
extern int MAXDEPTHSINGLE;
extern int NewTonTimesAllSINGLE;
extern int NewTonTimesSINGLE;
extern int NewTonSINGLE;
extern DYADIC SEPARATION;
extern int PHASE;

extern int MAX_TRISECTION;
extern int MAX_NEWTON;



// iRRAM namespace functions that is needed to define
// a vector of real numbers continuous.
namespace iRRAM{template <> struct is_continuous<std::vector< std::pair<REAL, int> > > : public std::true_type{};}
namespace iRRAM{
inline sizetype geterror( const std::vector<std::pair<REAL, int> > & l){	sizetype error, lerror;
	l[0].first.geterror(error);
	for (unsigned int i=0; i< l.size(); i++)
	{
		lerror = geterror(l[i].first);
		sizetype_max(error,error,lerror);
	}
	return error;}}
namespace iRRAM{
inline void seterror(std::vector<std::pair<REAL, int> > & l, sizetype & error)
{
	for (unsigned int i=0; i< l.size(); i++)
		l[i].first.seterror(error);
}}
POLYNOMIAL Giterate(POLYNOMIAL , int );
REAL trisect(int , POLYNOMIAL , REAL , REAL );



std::vector<std::pair<RATIONALINTERVAL, int> >
RealDistinctRootSeparationNaive(POLYNOMIAL , int );
std::vector<std::pair<RATIONALINTERVAL, int> >
RealDistinctRootSeparationBounded(POLYNOMIAL , int );

// below functions f(P, k) returns a 
// vector of pair of roots and their multiplicities
// where k is the number of distinct roots and 
// P is promised to have all of its roots on the real axis. 
std::vector<std::pair<REAL, int> >
RealRootFindingNaive(POLYNOMIAL , int );
std::vector<std::pair<REAL, int> >
RealRootFindingTrisection(POLYNOMIAL , int );
std::vector<std::pair<REAL, int> >
RealRootFindingNewton(POLYNOMIAL , int );
std::vector<std::pair<REAL, int> >
RealRootFindingCombined(POLYNOMIAL , int );
std::vector<std::pair<REAL, int> >
RealRootFindingBounded(POLYNOMIAL , int );
#endif
