#ifndef INTERVALS_H
#define INTERVALS_H

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include <utility>
using namespace iRRAM;



//********************************************//
//Open interval with real end points.         //
//The data is composed of center and radius   //
//********************************************//
class OPENINTERVAL
{
	public:
		REAL center;
		REAL radius;
		int depth;
		OPENINTERVAL();
		OPENINTERVAL(REAL, REAL);
		~OPENINTERVAL();

		OPENINTERVAL multiply(REAL);
		OPENINTERVAL subdivide(int );
};

//********************************************//
//interval with rational end points.          //
//The data is composed of center and radius   //
//********************************************//
class RATIONALINTERVAL
{
	public:
		RATIONAL center;
		RATIONAL radius;
		int depth;
		RATIONALINTERVAL();
		RATIONALINTERVAL(RATIONAL, RATIONAL);

		~RATIONALINTERVAL();

		RATIONALINTERVAL multiply(RATIONAL);

		RATIONALINTERVAL subdivide(int );
		RATIONALINTERVAL quaddivide(int );
};


//********************************************//
//A set of adjoint equal sized rational       //
//intervals (closed). The data is managed by  //
//enclosing interval and the number of        //
//subdivision.								  //
//********************************************//
class INTERVALCOMPONENT
{
	public:
		RATIONAL upper;
		RATIONAL lower;
		INTEGER Nc;
		int kc;
		int depth;
		INTEGER subdivision;
		INTERVALCOMPONENT();
		INTERVALCOMPONENT(RATIONALINTERVAL I);
		INTERVALCOMPONENT(RATIONAL a, RATIONAL b);

		~INTERVALCOMPONENT();

		void add(RATIONALINTERVAL I);
		
		INTEGER size();
		
		RATIONAL Wc();		
		RATIONAL wc();		
		RATIONAL min();		
		RATIONAL max();		
		RATIONAL Mc();		
		RATIONAL Rc();		
		bool isempty();		
		void split(INTEGER k);
		void split(int k);
		RATIONALINTERVAL operator [](int );
		RATIONALINTERVAL operator [](INTEGER );

};



// OPENINTERVAL FUNCTIONS
LAZY_BOOLEAN intersect(OPENINTERVAL, OPENINTERVAL);
OPENINTERVAL intersection(OPENINTERVAL , OPENINTERVAL );



// RATIONALINTERVAL FUNCTIONS
bool containedIn(RATIONALINTERVAL , RATIONALINTERVAL );
RATIONALINTERVAL intersection(RATIONALINTERVAL , RATIONALINTERVAL );

RATIONALINTERVAL operator + (const RATIONALINTERVAL&, const RATIONAL&);
RATIONALINTERVAL operator - (const RATIONALINTERVAL&);
RATIONALINTERVAL operator - (const RATIONALINTERVAL&, const RATIONAL&);
RATIONALINTERVAL operator * (const RATIONALINTERVAL&, const RATIONAL&);


bool adj(RATIONALINTERVAL , RATIONALINTERVAL );
bool intersect(RATIONALINTERVAL , RATIONALINTERVAL );



bool intersect(RATIONALINTERVAL , OPENINTERVAL );



// INTERVALCOMPONENT FUNCTIONS
bool adj(INTERVALCOMPONENT , RATIONALINTERVAL );
bool intersect(INTERVALCOMPONENT , RATIONALINTERVAL );
INTERVALCOMPONENT wraper(INTERVALCOMPONENT , OPENINTERVAL );
RATIONALINTERVAL COMPONENTtoINTERVAL(INTERVALCOMPONENT );








#endif
