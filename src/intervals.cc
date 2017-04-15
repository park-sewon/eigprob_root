#include "../include/intervals.h"
#include "../include/utilities.h"
#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include <stdexcept>

#include <stack>
#include <queue>
#include <utility>
#include <vector>

using namespace iRRAM;




// OPENINTERVAL MEMBER FUNCTIONS
OPENINTERVAL::OPENINTERVAL()
{
	center = 0;
	radius = 1/2;
	depth = 1;
}

OPENINTERVAL::OPENINTERVAL(REAL c, REAL r)
{
	center = c;
	radius = r;
	depth = 1;
}
OPENINTERVAL::~OPENINTERVAL()
{
}



OPENINTERVAL  OPENINTERVAL::multiply(REAL f)
{

	return OPENINTERVAL(center, radius*f); 
}

OPENINTERVAL  OPENINTERVAL::subdivide(int n)
{
	OPENINTERVAL S;
	if (n == 1)
		S =  OPENINTERVAL(center - radius/2, radius / 2);
	else if (n == 2)
		S =  OPENINTERVAL(center + radius/2, radius / 2);

	S.depth = depth + 1;
	return S;
}



// OPENINTERVAL FUNCTIONS

LAZY_BOOLEAN intersect(OPENINTERVAL A, OPENINTERVAL B)
{
	return abs(A.center - B.center) < A.radius + B.radius;
}

LAZY_BOOLEAN in(OPENINTERVAL A, REAL z)
{
	return abs(A.center - z) < A.radius;
}

OPENINTERVAL intersection(OPENINTERVAL A, OPENINTERVAL B)
{
	REAL a, b, c, d;
	a = A.center - A.radius;
	b = A.center + A.radius;
	c = B.center - B.radius;
	d = B.center + B.radius;
	return OPENINTERVAL((minimum(b,d) + maximum(a, c))/2, (minimum(b,d) - maximum(a, c))/2);
}





// RATIONALINTERVAL MEMBER FUNCTIONS

RATIONALINTERVAL::RATIONALINTERVAL()
{
	center = 0;
	radius = 1/2;
	depth = 1;
}
RATIONALINTERVAL::RATIONALINTERVAL(RATIONAL c, RATIONAL r)
{
	center = c;
	radius = r;
	depth = 1;
}
RATIONALINTERVAL::~RATIONALINTERVAL()
{
}

RATIONALINTERVAL  RATIONALINTERVAL::multiply(RATIONAL f)
{	
	RATIONALINTERVAL II = RATIONALINTERVAL(center, radius*f);
	II.depth = depth;

	return II;
}

RATIONALINTERVAL  RATIONALINTERVAL::subdivide(int n)
{
	RATIONALINTERVAL S;
	RATIONAL nr = radius / 2;
	if (n == 1)
		S =  RATIONALINTERVAL(center - nr, nr);
	else if (n == 2)
		S =  RATIONALINTERVAL(center + nr, nr);

	S.depth = depth + 1;
	return S;
}
RATIONALINTERVAL  RATIONALINTERVAL::quaddivide(int n)
{
	RATIONALINTERVAL S;
	if (n == 1)
		S =  RATIONALINTERVAL(center - radius*3/4, radius / 4);
	else if (n == 2)
		S =  RATIONALINTERVAL(center - radius*1/4, radius / 4);
	else if (n == 2)
		S =  RATIONALINTERVAL(center + radius*1/4, radius / 4);
	else if (n == 2)
		S =  RATIONALINTERVAL(center + radius*3/4, radius / 4);
	S.depth = depth + 1;
	return S;
}

RATIONALINTERVAL operator + (const RATIONALINTERVAL& I, const RATIONAL& Q)
{
	return RATIONALINTERVAL(I.center + Q, I.radius);
}
RATIONALINTERVAL operator - (const RATIONALINTERVAL& I)
{
	return RATIONALINTERVAL(-I.center, I.radius);

}
RATIONALINTERVAL operator - (const RATIONALINTERVAL& I, const RATIONAL& Q)
{
	return RATIONALINTERVAL(I.center - Q, I.radius);
}
RATIONALINTERVAL operator * (const RATIONALINTERVAL& I, const RATIONAL& Q)
{
	return RATIONALINTERVAL(I.center, I.radius*Q); 
}

// RATIONALINTERVAL FUNCTIONS


bool containedIn(RATIONALINTERVAL A, RATIONALINTERVAL B)
{
	return (B.center - B.radius <= A.center - A.radius && B.center + B.radius <= A.center +A.radius);
}


RATIONALINTERVAL intersection(RATIONALINTERVAL A, RATIONALINTERVAL B)
{
	RATIONAL a, b, c, d;
	a = A.center - A.radius;
	b = A.center + A.radius;
	c = B.center - B.radius;
	d = B.center + B.radius;
	return RATIONALINTERVAL((minimum(b,d) + maximum(a, c))/2, (minimum(b,d) - maximum(a, c))/2);
}

bool adj(RATIONALINTERVAL A, RATIONALINTERVAL B)
{
	RATIONAL d = A.center - B.center ;
	if (d < 0)
		d = -d;
	if(d == A.radius + B.radius)
		return true;
	return false;
}

bool intersect(RATIONALINTERVAL A, RATIONALINTERVAL B)
{
	RATIONAL d = A.center - B.center;
	if (d < 0)
		d = -d;
	if(d <= A.radius + B.radius)
		return true;
	return false;
}

bool intersect(RATIONALINTERVAL C, OPENINTERVAL I)
{
	REAL d = abs(C.center - I.center);
	REAL r = abs(C.radius + I.radius);
	return choose(r<d, C.radius * 3 / 2 +I.radius > d) == 2;
}



// INTERVALCOMPONENT MEMBER FUNCTIONS

INTERVALCOMPONENT::INTERVALCOMPONENT()
{
	subdivision = 0;
	depth = 0;

}
INTERVALCOMPONENT::INTERVALCOMPONENT(RATIONALINTERVAL I)
{
	depth = 0;
	lower = I.center - I.radius;
	upper = I.center + I.radius;
	subdivision = 1;
}
INTERVALCOMPONENT::~INTERVALCOMPONENT()
{}
void INTERVALCOMPONENT::add(RATIONALINTERVAL I)
{
	if(subdivision == 0)
	{
		upper = I.center + I.radius;
		lower = I.center - I.radius;
	}
	else
	{
		if(lower > I.center)
			lower = I.center - I.radius;
		else
			upper = I.center + I.radius;
	}
	subdivision += 1;
}
INTEGER INTERVALCOMPONENT::size()
{
	return subdivision;
}
RATIONAL INTERVALCOMPONENT::Wc()
{
	return upper - lower;
}
RATIONAL INTERVALCOMPONENT::wc()
{
	if (0 == subdivision)
		throw std::invalid_argument( "sub-interval out of index" );


	return (upper - lower) / subdivision;
}
RATIONAL INTERVALCOMPONENT::min()
{
	return lower;
}
RATIONAL INTERVALCOMPONENT::max()
{	
	return upper;
}
RATIONAL INTERVALCOMPONENT::Mc()
{
	return (upper + lower)/2;
}
RATIONAL INTERVALCOMPONENT::Rc()
{
	return (upper - lower) * 3 / 4;
}
bool INTERVALCOMPONENT::isempty()
{
	return subdivision == 0;
}
void INTERVALCOMPONENT::split(INTEGER k)
{
	depth += 1;
	subdivision = subdivision *  k;
}

void INTERVALCOMPONENT::split(int k)
{
	depth += 1;
	subdivision = subdivision *  k;
}

RATIONALINTERVAL toINTERVAL(INTERVALCOMPONENT C)
{
	return RATIONALINTERVAL((C.min()+C.max())/2, (C.max()-C.min())/2);
}

RATIONALINTERVAL INTERVALCOMPONENT::operator [](int i)
{
	if (i>= subdivision)
		throw std::invalid_argument( "sub-interval out of index" );

	RATIONAL center, radius;

	radius = (upper - lower)/2;
	assert (subdivision != 0);
	return RATIONALINTERVAL(lower + radius / subdivision + i*2*radius / subdivision , radius / subdivision);
};

RATIONALINTERVAL INTERVALCOMPONENT::operator [](INTEGER i)
{
	if (i>= subdivision)
		throw std::invalid_argument( "sub-interval out of index" );

	RATIONAL center, radius;

	radius = (upper - lower)/2;
	assert (subdivision != 0);
	return RATIONALINTERVAL(lower + radius / subdivision + i*2*radius / subdivision , radius / subdivision);
};
bool adj(INTERVALCOMPONENT C, RATIONALINTERVAL I)
{
	return adj(toINTERVAL(C), I);
}


bool intersect(INTERVALCOMPONENT C, RATIONALINTERVAL I)
{
	return intersect(toINTERVAL(C), I);
}


INTERVALCOMPONENT wraper(INTERVALCOMPONENT C, OPENINTERVAL I)
{
	INTERVALCOMPONENT w =	INTERVALCOMPONENT();
	for (INTEGER i=0 ; i< C.size(); i= i + 1)
	{
		if (intersect(C[i], I))
		{
			w.add(C[i]);
		}
	}
	w.depth = C.depth;
	w.Nc = C.Nc;
	w.kc = C.kc;
	return w; 
}


