#include "../include/polynomial.h"
#include "../include/utilities.h"

#include "iRRAM/lib.h"
#include "iRRAM/core.h"

#include <stack>
#include <queue>
#include <utility>
#include <vector>
#include <cmath>
using namespace iRRAM;


int geterror_exp(const POLYNOMIAL & P)
{
	sizetype error; 
  	P.coef[0].geterror(error);
  	int ex = error.exponent;
  	int tmp;

	for(int i=1; i< P.degree; i++)
	{
		P.coef[0].geterror(error);
		tmp = error.exponent;
	  	if (tmp >ex)
	  		ex = tmp;
	}

	return tmp;
}



POLYNOMIAL::POLYNOMIAL(int d, REAL* c)
{
	std::vector<REAL> co;
	co.reserve(d+1);
	for (int i = 0; i<d+1; i++)
	{	
		co.push_back(c[i]);
	}
	degree = d;
	coef = co;
}


POLYNOMIAL::POLYNOMIAL(int d, std::vector<REAL> v)
{
	std::vector<REAL> co;
	co.reserve(d+1);
	for (int i = 0; i<d+1; i++)
	{	
		co.push_back(v[i]);
	}
	degree = d;
	coef = co;
}

POLYNOMIAL::POLYNOMIAL(REAL x)
{
	std::vector<REAL> co;
	co.reserve(1);
	co.push_back(x);
	degree = 0;
	coef = co;
}


POLYNOMIAL::POLYNOMIAL()
{
}

POLYNOMIAL::~POLYNOMIAL()
{
}


REAL POLYNOMIAL::operator () (const REAL& z)
{
	REAL fx = 0;
	for(int i=0;i<degree + 1;i++)
		fx = fx + power(z,i) * coef[i];
	return fx;

}


POLYNOMIAL operator + (const POLYNOMIAL& P , const POLYNOMIAL& Q)
{
	int d;
	if ( P.degree > Q.degree)
		d =  P.degree;
	else
		d = Q.degree;

	REAL c[d+1];
	for(int i=0; i< d+1; i++)
	{
		if(i < P.degree +1 && i < Q.degree+1)
			c[i] = P.coef[i] + Q.coef[i];
		else if (i > P.degree)
			c[i] = Q.coef[i];
		else if (i > Q.degree)
			c[i] = P.coef[i];
	}

	POLYNOMIAL R = POLYNOMIAL(d, c);
	return R;
}
POLYNOMIAL operator + (const POLYNOMIAL& P , const REAL& a)
{
	REAL c[P.degree+1];
	for(int i=0; i< P.degree+1; i++)
	{
		c[i] = P.coef[i];
	}
	c[0] += a;
	return POLYNOMIAL(P.degree, c);
}
POLYNOMIAL operator + (const REAL& a , const POLYNOMIAL& P)
{
	return P + a;
}

POLYNOMIAL operator - (const POLYNOMIAL& P)
{

	REAL c[P.degree+1];
	for(int i=0; i< P.degree+1; i++)
		c[i] =  - P.coef[i];

	POLYNOMIAL R = POLYNOMIAL(P.degree, c);
	return R;
}
POLYNOMIAL operator - (const POLYNOMIAL& P, const POLYNOMIAL& Q)
{
	int d;
	int p = P.degree;
	int q = Q.degree;

	if (p > q)
		d = p;
	else
		d = q;

	REAL c[d+1];

	for(int i=0; i< d+1; i++)
	{
		if(i < p +1 && i < q+1)
			c[i] = P.coef[i] - Q.coef[i];
		else if (i > p)
			c[i] = 0 - Q.coef[i];
		else if (i > q)
			c[i] = P.coef[i];
	}

	POLYNOMIAL R = POLYNOMIAL(d, c);
	return R;
}
POLYNOMIAL operator - (const POLYNOMIAL& P, const REAL& c)
{
	return P + REAL(-c);
}
POLYNOMIAL operator - (const REAL& c, const POLYNOMIAL& P)
{
	return (-P) + c;
}

	

POLYNOMIAL operator * (const POLYNOMIAL& P, const POLYNOMIAL& Q)
{
	int d = P.degree + Q.degree;
	int p = P.degree;
	int q = Q.degree;
	REAL c[d+1];

	for (int i=0; i < p + 1; i++)
	{
		for(int j=0; j < q +1; j++)
		{
			c[i+j] = c[i+j] +  P.coef[i] * Q.coef[j];
		}
	}
	POLYNOMIAL R = POLYNOMIAL(d, c);
	
	return R;
}

POLYNOMIAL operator * (const POLYNOMIAL& P , const REAL& a)
{
	REAL c[P.degree+1];
	for(int i=0; i< P.degree+1; i++)
	{
		c[i] = P.coef[i] * a;
	}
	return POLYNOMIAL(P.degree, c);
}
POLYNOMIAL operator * (const REAL& c, const POLYNOMIAL& P)
{
	return P*c;
}
POLYNOMIAL operator / (const POLYNOMIAL& P, const REAL& c)
{
	return P*REAL(1/c);
}



REAL evaluate(POLYNOMIAL P, REAL z)
{
	REAL fx = 0;
	for(int i=0;i<P.degree + 1;i++)
		fx = fx + power(z,i) * P.coef[i];
	return fx;

}

POLYNOMIAL deriv(POLYNOMIAL P, int k)
{
	int d;
	if (k == 0)
		return POLYNOMIAL(P.degree, P.coef);

	if (k > P.degree)
	{
		d = 0;
		REAL C[1];
		C[0] = 0;
		return POLYNOMIAL(d, C);

	}
	else
	{
		d = P.degree - k;
		REAL tmp;
		REAL C[d+1];
		for(int i=0; i< d+1; i++)
		{
			tmp = 1;
			for (int j=0; j< k; j++)
			{
				tmp *= (j+i+1);
			}
			C[i] = P.coef[k+i] * tmp;
		}
	
		return POLYNOMIAL(d, C);
	}
}

// kth taylor coefficient of P at z
REAL CoefAt(POLYNOMIAL P, int k, REAL z)
{
	REAL tmp[P.degree+1];
	for(int i=0; i<P.degree + 1; i++)
		tmp[i] = P.coef[i];

	for(int j=0; j<k; j++)
	{
		for(int i=0;i<P.degree-j;i++)
		{
			tmp[i] = tmp[i+1] * (i+1);
		}
	}
	REAL fx = 0;
	for(int i=0;i < P.degree + 1 - k; i++)
	{
		fx = fx + power(z,i) * tmp[i];
	}
	return  fx / factorial(k);
}	

// translation and dilation by g(x) = f(a + bx)
POLYNOMIAL translation(POLYNOMIAL P, REAL a, REAL b)
{
	REAL C[P.degree +1];
	for(int i=0; i<P.degree + 1; i++)
	{
		C[i] = power(a, i) * deriv(P,i)(b) / factorial(i);
	}

	return POLYNOMIAL(P.degree, C);
}

// Chebyshev polynomial of degree n
POLYNOMIAL Chebyshev(int n)
{
	REAL co[1];
	co[0] = 1;
	POLYNOMIAL Q = POLYNOMIAL(0,co);

	if(n == 0) return Q;
	REAL c[2];
	c[0] = 0;
	c[1] = 1;
	POLYNOMIAL P = POLYNOMIAL(1, c);
	if(n == 1) return P;
	POLYNOMIAL R;
	REAL coo[2];
	coo[0] = 0;
	coo[1] = 2;
	POLYNOMIAL tx = POLYNOMIAL(1, coo);
	for(int k=2; k<n+1; k++)
	{
		R = (tx * P) - Q;
		Q = P;
		P = R;
	}	
	return R;
}



// monomial Polynomial having roots uniformly in [-1,1]
// n >= 2 
POLYNOMIAL Uniform(int n)
{
	REAL co[2];
	co[1] = 1;
	co[0] = 1;
	POLYNOMIAL Q = POLYNOMIAL(1, co);
	for(int i=1; i<n; i++)
	{
		co[0] = 1 - REAL(2)*i / (n-1);
		Q = Q * POLYNOMIAL(1, co);
	}
	return Q;
}

// monomial Polynomial having roots uniformly in [-1,1]
// n >= 2 
POLYNOMIAL Wilkinson(int n, REAL b)
{
	REAL co[2];
	co[1] = 1;
	co[0] = -1;
	POLYNOMIAL Q = POLYNOMIAL(1, co);
	for(int i=1; i<n; i++)
	{
		co[0] = 0 - REAL(1) / power(b,i);
		Q = Q * POLYNOMIAL(1, co);
	}
	return Q;
}
