#include "../include/utilities.h"
#include "iRRAM/lib.h"
#include "iRRAM/core.h"

#include <stack>
#include <queue>
#include <utility>
#include <math.h>  
#include <vector>

using namespace iRRAM;


INTEGER factorial(int n)
{
	INTEGER ans = 1;
	for(int i=1;i<n+1;i++)
		ans*=i;
	return ans;
}

void print(COMPLEX c)
{
	cout << real(c) << " + i" << imag(c);
}

COMPLEX power(COMPLEX C, int n)
{
	COMPLEX D = COMPLEX(1,0);
	
	for(int i=0; i<n; i++)
		D = D*C;
	if(n<0)
	{
		for(int i=0; i>n; i--)
			D =D/ C;
	}
	return D;

}

DYADIC minimum(DYADIC a, DYADIC b)
{
	if ( a > b)
		return b;
	return a;
}

DYADIC maximum(DYADIC a, DYADIC b)
{
	if ( a > b)
		return a;
	return b;
}


RATIONAL minimum(RATIONAL a, RATIONAL b)
{
	if ( a > b)
		return b;
	return a;
}

RATIONAL maximum(RATIONAL a, RATIONAL b)
{
	if ( a > b)
		return a;
	return b;
}


void print(OPENINTERVAL D)
{
	cout <<"("<<D.center - D.radius<<", "<< D.center + D.radius<<")\n";
}
void print(RATIONALINTERVAL D)
{
	cout <<"("<<D.center - D.radius<<", "<< D.center + D.radius<<")\n";
}

void print(INTERVALCOMPONENT C)
{
	cout <<"--- components ---\n";
	for(int i=0; i < C.size(); i++)
		print(C[i]);
	cout <<"------------------\n";
}

void print(POLYNOMIAL P)
{
	for(int i=0;i<P.degree+1;i++)
	{
		if(i!=0)
			cout<<" + ";
		cout<< P.coef[P.degree-i] << "x^" << P.degree -i<< " ";
	}
	cout<<"\n";
}
