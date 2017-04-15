#ifndef RPOLYNOMIAL_H
#define RPOLYNOMIAL_H

#include "iRRAM/lib.h"
#include "iRRAM/core.h"

#include <utility>
using namespace iRRAM;


class POLYNOMIAL
{
	public:
		std::vector<REAL> coef;
		int degree;

		POLYNOMIAL(int, REAL *);
		POLYNOMIAL(int, std::vector<REAL>);
		POLYNOMIAL(REAL);
		POLYNOMIAL();
		~POLYNOMIAL();
		REAL operator () (const REAL&);
};

POLYNOMIAL operator + (const POLYNOMIAL&, const POLYNOMIAL&);
POLYNOMIAL operator + (const POLYNOMIAL&, const REAL&);
POLYNOMIAL operator + (const REAL&, const POLYNOMIAL&);
POLYNOMIAL operator - (const POLYNOMIAL&);
POLYNOMIAL operator - (const POLYNOMIAL&, const POLYNOMIAL&);
POLYNOMIAL operator - (const POLYNOMIAL&, const REAL&);
POLYNOMIAL operator - (const REAL&, const POLYNOMIAL&);
POLYNOMIAL operator * (const POLYNOMIAL&, const POLYNOMIAL&);
POLYNOMIAL operator * (const POLYNOMIAL&, const REAL&);
POLYNOMIAL operator * (const REAL&, const POLYNOMIAL&);
POLYNOMIAL operator / (const POLYNOMIAL&, const REAL&);

namespace iRRAM{template <> struct is_continuous<POLYNOMIAL> : public std::true_type{};}
// void geterror(const POLYNOMIAL &);
int geterror_exp(const POLYNOMIAL & );

// Q = deriv(P,k) := Q^{(k)}
POLYNOMIAL deriv(POLYNOMIAL, int);
// y = evaluate(P, x) := P(x)
REAL evaluate(POLYNOMIAL, REAL);
// Taylor coefficients: a_k = CoefAt(P, k, z) := P^{(k)}(z)/k!
REAL CoefAt(POLYNOMIAL, int, REAL);

// Example polynomials
POLYNOMIAL Chebyshev(int );
POLYNOMIAL Uniform(int );
POLYNOMIAL Wilkinson(int , REAL );

// Translation and dilation by g(x) = f(ax + b)
POLYNOMIAL translation(POLYNOMIAL , REAL , REAL );





#endif
