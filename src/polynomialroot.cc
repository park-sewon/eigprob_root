#include "../include/polynomialroot.h"

#include "../include/intervals.h"
#include "../include/utilities.h"

#include "iRRAM/lib.h"
#include "iRRAM/core.h"

#include <stack>
#include <queue>
#include <utility>
#include <vector>
#include <cmath>
#include <time.h>

using namespace iRRAM;

// global variable for splitting computation paths throughout re-iterations:
float TIME_ROOT_SEPARATION = 0;
float TIME_ROOT_FINDING = 0;
float TIME_ROOT_APPROXIMATION = 0;
float TIME_GAUSSIAN_ELIMINATION = 0;
float TIME_CHAR_POLYNOMIAL =0;

int MAXDEPTH = 0;
int NewTonTimesAll = 0;
int NewTonTimes = 0;
int NewTon = 0;
DYADIC SEPARATION = 0;
int MAX_TRISECTION=0;
int MAX_NEWTON=0;



// < relation to use pre-impelented sorting function and prioirty queue.
struct less_than_key
{
    inline bool operator() (const std::pair<RATIONALINTERVAL, int>& C1, 
    	const std::pair<RATIONALINTERVAL, int>& C2)
	{
		return (C1.first.center < C2.first.center);
	}
};
struct IC_LESS
{
    inline bool operator() (const INTERVALCOMPONENT& C1, 
    	const INTERVALCOMPONENT& C2)
	{
		return (C1.lower + C1.upper < C2.lower + C2.upper);
	}
};

struct wider
{
    inline bool operator() (const std::pair<RATIONALINTERVAL, int>& C1, 
    	const std::pair<RATIONALINTERVAL, int>& C2)
	{
		return (C1.first.radius < C2.first.radius);
	}
};



std::vector<std::pair<RATIONALINTERVAL, int> > 
casting(std::vector<INTERVALCOMPONENT> P)
{

	std::vector<std::pair<RATIONALINTERVAL, int> > Q;

	for (int i =0; i<(int) P.size(); i++)
	{	
		Q.push_back(std::pair<RATIONALINTERVAL, int> (RATIONALINTERVAL(P[i].Mc(), P[i].Wc()/2), P[i].kc));
	}
	return Q;
} 




// Graeffe iteration 
POLYNOMIAL polyE(POLYNOMIAL P)
{
	double td = P.degree / 2;
	int d = std::floor(td);

	REAL C[d+1];

	for(int i = 0; i<d+1; i++)
	{
		C[i] = P.coef[i*2];
	}
	return POLYNOMIAL(td, C);
}

POLYNOMIAL polyO(POLYNOMIAL P)
{
	double td = (P.degree-1) / 2;
	int d = std::floor(td);

	REAL C[d+1];

	for(int i = 0; i<d+1; i++)
	{
		C[i] = P.coef[1+i*2];
	}
	return POLYNOMIAL(td, C);
}

POLYNOMIAL Giterate(POLYNOMIAL P, int N)
{

	POLYNOMIAL Q = POLYNOMIAL(P.degree, P.coef);
	REAL c[2];
	c[0] = 0; c[1] = 1;
	int k = pow(-1, P.degree);
	POLYNOMIAL T = POLYNOMIAL(1, c);
	if(k > 0)
		for(int i=1; i<N+1; i++)
		{
			Q = (polyE(Q) * polyE(Q)) - (T * (polyO(Q) * polyO(Q)));
		}
	else
		for(int i=0; i<N; i++)
		{
			Q = (T*(polyO(Q) * polyO(Q))) - (polyE(Q) * polyE(Q));
		}
	return Q;
}

// Soft pallet test: \tilde{\mathcal{T}}_k(P,D)
bool softTTest(POLYNOMIAL P, int k, RATIONALINTERVAL D)
{
	REAL m = D.center;
	REAL r = D.radius;
	REAL LHS = abs(CoefAt(P,k,m)) * power(r,k);
	REAL RHS = 0;
	bool ans;
	for(int i=0;i<P.degree+1;i++)
	{
		if(i!=k)
		{
			RHS += abs(CoefAt(P,i,m)) * power(r,i);
		}
	}
	ans = choose(LHS>RHS, LHS<RHS, LHS*2 < RHS*3 && 2*RHS < 3*LHS) == 1;
	return ans;
}

// Soft pallet test on Graeffe iteration: \tilde{\mathcal{T}}_k(G_{\log(1+\log d) + 5}(P_D),D(0,1))
// on a interval with rational endpoints
bool softGTest(POLYNOMIAL P, int k, RATIONALINTERVAL D)
{
	int N = std::ceil(std::log(1+(std::log(P.degree) / std::log(2)))/std::log(2)) + 5;
	POLYNOMIAL Q = Giterate(translation(P, D.radius, D.center), N);
	RATIONALINTERVAL U(0, 1);
	return softTTest(Q, k, U);
}


// Soft pallet test on Graeffe iteration: \tilde{\mathcal{T}}_k(G_{\log(1+\log d) + 5}(P_D),D(0,1))
// on a interval with real endpoints
bool softGTest(POLYNOMIAL P, int k, OPENINTERVAL D)
{
	int N = std::ceil(std::log(1+(std::log(P.degree) / std::log(2)))/std::log(2)) + 5;
	POLYNOMIAL Q = Giterate(translation(P, D.radius, D.center), N);
	RATIONALINTERVAL U(0, 1);
	return softTTest(Q, k, U);
}


// First i such that T_i(G_N(f_{m,r}), 0, 1) /\ T_i(G_N(f_{m,3r}), 0, 1)holds 
//    0 if T_0(G_N(f_{m,r}), 0, 1) 
//   -1 if - T_i(G_N(f_{m,r}), 0, 1) for all i
// -i-1 if T_i(G_N(f_{m,r}), 0, 1) /\  -T_i(G_N(f_{m,3r}), 0, 1)
int softGStarThree(POLYNOMIAL P, RATIONALINTERVAL B, int k)
{
	int N = std::ceil(std::log(1+(std::log(P.degree) / std::log(2)))/std::log(2)) + 5;
	POLYNOMIAL Q = Giterate(translation(P, B.radius, B.center), N);
	POLYNOMIAL R = Giterate(translation(P, 3*B.radius, B.center), N);

	RATIONALINTERVAL U(0, 1);

	for(int i =0;i<k+1;i++)
	{
		if(softTTest(Q, i, U))
		{
			if (i == 0)
				return 0;
			if(softTTest(R, i, U))
				return i;
			else
				return -i-1;
		}
	}
	return - 1;
}

// First i such that \tilde{T}_i(G_N(f_{m,r}), 0, 1) holds 
//  -1 if all fail
int softGStar(POLYNOMIAL P, RATIONALINTERVAL B, int k)
{
	int N = std::ceil(std::log(1+(std::log(P.degree) / std::log(2)))/std::log(2)) + 5;
	POLYNOMIAL Q = Giterate(translation(P, B.radius, B.center), N);
	RATIONALINTERVAL U(0, 1);

	for(int i =0;i<k+1;i++)
		if(softTTest(Q, i, U))
			return i;
	return - 1;
}



// Return a list of (D_i, d_i) such that each D_i is disjoint (1,3)--isolating for positive number of roots.
// Implementation based on simple naive subdivision algorithm with Graeffe's iteration
std::vector<std::pair<RATIONALINTERVAL, int> >
RealDistinctRootSeparationNaive(POLYNOMIAL P, int distinct)
{

	int NODEINDEPTH[4096] = {0};

	int n = P.degree;
	RATIONAL min = 1;
	std::pair<RATIONALINTERVAL, int> tmp;

	RATIONALINTERVAL D = RATIONALINTERVAL(0,1);
	while(!softGTest(P, n, D))
	{
		D = D.multiply(2);
	}

	RATIONAL epsilon = D.radius / 2 / distinct;


	std::queue<RATIONALINTERVAL > Q1;
	Q1.push(D);

    std::priority_queue<std::pair<RATIONALINTERVAL, int> , std::vector<std::pair<RATIONALINTERVAL, int> >, wider> Q0;

	int k;
	int maxDepth = 0;
	int f = 0;
	int ff = 0;

	std::vector<std::pair<RATIONALINTERVAL, int> > solutionDiscs;
	solutionDiscs.reserve(P.degree);

	// while(!Q1.empty())
	// {
	// 	RATIONALINTERVAL a = Q1.front();
	// 	Q1.pop();
	// 	if (softGTest(P, 0, a) == false)
	// 	{
	// 		if (a.radius < epsilon)
	// 			Q0.push(std::pair<RATIONALINTERVAL, int>(a, n));
	// 		else
	// 			for(int si=1;si<3;si++)
	// 				Q1.push(a.subdivide(si));
	// 	}
	// }
	
	RATIONALINTERVAL a = Q1.front();
	Q0.push(std::pair<RATIONALINTERVAL, int>(a, n));



	while(!Q0.empty())
	{
		std::pair<RATIONALINTERVAL, int> fp = Q0.top();	Q0.pop();
		RATIONALINTERVAL a = fp.first;
		int nroot = fp.second;

		ff = 0;
		for (int i = 0; i < solutionDiscs.size(); i++)
		{
			if (intersect(a, solutionDiscs[i].first))
			{
				ff = 1; break;
			}
		}
		if (ff == 0) 
		{
			NODEINDEPTH[a.depth]+=1;


			maxDepth = std::max(maxDepth, a.depth);
			k = softGStarThree(P,a.multiply(RATIONAL(3)/2),nroot);

			if (k==-1)   // failed
			{	
				for(int si=1;si<3;si++)
					Q0.push(std::pair<RATIONALINTERVAL, int>(a.subdivide(si), nroot));

			}
			else if (k < 0)	// T(3D) failed but T(D) succeeded with i = -k - 1
			{
				for(int si=1;si<3;si++)
					Q0.push(std::pair<RATIONALINTERVAL, int>(a.subdivide(si), -k - 1));

			}
			else if (k != 0) // D, 3D has k > 0 roots
			{
				if(a.radius * 3 / 2 > epsilon)
				{
					for(int si=1;si<3;si++)
						Q0.push(std::pair<RATIONALINTERVAL, int>(a.subdivide(si), k));
				}
				else
				{
					solutionDiscs.push_back(std::pair<RATIONALINTERVAL, int> (a.multiply(RATIONAL(3)/2), k));
				}
			}
		}

		if (Q0.empty() && solutionDiscs.size() != distinct)
		{
			for(int i=0; i< (int) solutionDiscs.size(); i++)
			{
				NODEINDEPTH[solutionDiscs[i].first.depth] -=1;
				Q0.push(solutionDiscs[i]);
				min = minimum(min, solutionDiscs[i].first.radius);
			}
			solutionDiscs.clear();
			epsilon =  min / 2;		
		}

	}

	
	std::sort(solutionDiscs.begin(), solutionDiscs.end(), less_than_key());

	MAXDEPTH =  maxDepth;

	for(int i=0; i<MAXDEPTH+1; i++)
		cout << NODEINDEPTH[i]<<", ";
	cout<<"\n";


	return solutionDiscs;
}



// Return a list of (D_i, d_i) such that each D_i is disjoint (1,3)--isolating for positive number of roots.
// Implementation based on simple naive subdivision algorithm with Graeffe's iteration
std::vector<std::pair<RATIONALINTERVAL, int> >
RealDistinctRootSeparationBounded(POLYNOMIAL P, int distinct)
{
	int NODEINDEPTH[4096] = {0};

	int n = P.degree;
	RATIONAL min = 1;
	std::pair<RATIONALINTERVAL, int> tmp;

	RATIONALINTERVAL D = RATIONALINTERVAL(0,1);
	while(!softGTest(P, n, D))
	{
		D = D.multiply(2);
	}

	RATIONAL epsilon = D.radius / 2 / distinct;

    std::priority_queue<std::pair<RATIONALINTERVAL, int> , std::vector<std::pair<RATIONALINTERVAL, int> >, wider> Q0;

	int k;
	int maxDepth = 0;
	int f = 0;
	int ff = 0;

	std::vector<std::pair<RATIONALINTERVAL, int> > solutionDiscs;
	solutionDiscs.reserve(P.degree);

	Q0.push(std::pair<RATIONALINTERVAL, int>(D, n));

	// main loop
	int intersectionChecker = 0;
	int solnum=0;
	int curruntDepth = 0;
	while(1)
	{


		std::pair<RATIONALINTERVAL, int> fp = Q0.top();	Q0.pop();
		RATIONALINTERVAL a = fp.first;
		int nroot = fp.second;

		if (curruntDepth != a.depth)
		{
			curruntDepth = a.depth;
			if(solutionDiscs.size() == distinct)
				break;
			else
				solutionDiscs.clear();
		}


		NODEINDEPTH[a.depth]+= 1;


		maxDepth = std::max(maxDepth, a.depth);
		
		k = softGStarThree(P,a.multiply(RATIONAL(3)/2), nroot);

		if (k==-1)   // failed
		{	
			for(int si=1;si<3;si++)
				Q0.push(std::pair<RATIONALINTERVAL, int>(a.subdivide(si), nroot));

		}
		else if (k < 0)	// T(3D) failed but T(D) succeeded with i = -k - 1
		{
			for(int si=1;si<3;si++)
				Q0.push(std::pair<RATIONALINTERVAL, int>(a.subdivide(si), -k - 1));
		}
		else if (k != 0) // D, 3D has k > 0 roots
		{
			// merge push
			intersectionChecker = 0;
			for (unsigned int j=0; j< solutionDiscs.size(); j++)
			{
				if (intersect(solutionDiscs[j].first, a.multiply(RATIONAL(3)/2)))
				{
					solutionDiscs[j].first = intersection(a.multiply(RATIONAL(3)/2), solutionDiscs[j].first);
					intersectionChecker = 1;
					break;
				}

			}
			if (intersectionChecker == 0)
				solutionDiscs.push_back(std::pair<RATIONALINTERVAL, int> (a.multiply(RATIONAL(3)/2), k));
		

			// push it back to the main queue:
			for(int si=1;si<3;si++)
				Q0.push(std::pair<RATIONALINTERVAL, int>(a.subdivide(si), nroot));

		}

	}

	
	std::sort(solutionDiscs.begin(), solutionDiscs.end(), less_than_key());

	MAXDEPTH =  maxDepth;
	for(int i=0; i<MAXDEPTH+1; i++)
		cout << NODEINDEPTH[i]<<", ";
	cout<<"\n";

	return solutionDiscs;
}


std::vector<std::pair<REAL, int> >
RealDistinctRootApproximationNaive(int prec, POLYNOMIAL P, int distinct)
{
	int d = P.degree;
	RATIONAL min = 1;
	std::pair<RATIONALINTERVAL, int> tmp;

	RATIONALINTERVAL D = RATIONALINTERVAL(0,1);
	while(!softGTest(P, d, D))
	{
		D = D.multiply(2);
	}

	REAL epsilon = power(2, prec);


	std::queue<RATIONALINTERVAL > Q1;
	Q1.push(D);

    std::priority_queue<std::pair<RATIONALINTERVAL, int> , std::vector<std::pair<RATIONALINTERVAL, int> >, wider> Q0;

	int k;
	int maxDepth = 0;
	int f = 0;
	int ff = 0;
	std::vector<std::pair<REAL, int> > distinctRoots;

	std::vector<std::pair<RATIONALINTERVAL, int> > solutionDiscs;
	solutionDiscs.reserve(P.degree);
	Q0.push(std::pair<RATIONALINTERVAL, int>(D, d));

	while(!Q0.empty())
	{
		std::pair<RATIONALINTERVAL, int> fp = Q0.top();	Q0.pop();
		RATIONALINTERVAL a = fp.first;
		int nroot = fp.second;

		ff = 0;
		for (int i = 0; i < solutionDiscs.size(); i++)
		{
			if (intersect(a, solutionDiscs[i].first))
			{
				ff = 1; break;
			}
		}
		if (ff == 0) 
		{
			maxDepth = std::max(maxDepth, a.depth);
			k = softGStarThree(P,a.multiply(RATIONAL(3)/2),nroot);

			if (k==-1)   // failed
			{	
				for(int si=1;si<3;si++)
					Q0.push(std::pair<RATIONALINTERVAL, int>(a.subdivide(si), nroot));

			}
			else if (k < 0)	// T(3D) failed but T(D) succeeded with i = -k - 1
			{
				for(int si=1;si<3;si++)
					Q0.push(std::pair<RATIONALINTERVAL, int>(a.subdivide(si), -k - 1));

			}
			else if (k != 0) // D, 3D has k > 0 roots
			{
				if(a.radius * 3 / 2 > epsilon)
				{
					for(int si=1;si<3;si++)
						Q0.push(std::pair<RATIONALINTERVAL, int>(a.subdivide(si), k));
				}
				else
				{
					solutionDiscs.push_back(std::pair<RATIONALINTERVAL, int> (a.multiply(RATIONAL(3)/2), k));
				}
			}
		}

		if (Q0.empty() && solutionDiscs.size() != distinct)
		{
			for(int i=0; i< (int) solutionDiscs.size(); i++)
			{
				Q0.push(solutionDiscs[i]);
				// min = minimum(min, solutionDiscs[i].first.radius);
			}
			solutionDiscs.clear();
			epsilon =  epsilon / 2;		
		}

	}




	
	std::sort(solutionDiscs.begin(), solutionDiscs.end(), less_than_key());
	for(int i=0; i< distinct; i++)
		distinctRoots.push_back(std::pair< REAL, int >(REAL(solutionDiscs[i].first.center), solutionDiscs[i].second));

	MAXDEPTH =  maxDepth;

	return distinctRoots;
}


REAL trisection(int prec, POLYNOMIAL P, REAL a, REAL b)
{	
	int depth = 0;
	REAL returnvalue;
	REAL epsilon = power(2,prec);
	while (1)
	{
		if (choose(P((2*a + b)/3) *  P(b)  <0, P((a + 2*b)/3) *  P(a)  <0) == 1)
		{
			a = (2*a  + b) / 3;
		}
		else
		{
			b = (a + 2*b ) / 3;
		}
		depth += 1;
		if (choose(abs(a-b) /2< epsilon, abs(a-b) > epsilon) == 1)
		{	


			MAX_TRISECTION = std::max(MAX_TRISECTION, depth);
			returnvalue = (b+a)/2;
			return a;
		}
	}

}
REAL RootApprox(int p, POLYNOMIAL P, REAL z1, REAL z2)
{ 
	REAL epsilon = power(2, p);
	REAL delta = abs(z2 - z1);
	if (choose(delta < epsilon, delta > epsilon /2) == 1)
		return z1;
	else
		if (choose(P(z1) != 0, P(z2) != 0) == 1 )
			return z2;
		else
			return z1;
}

INTERVALCOMPONENT
Newton(POLYNOMIAL P, INTERVALCOMPONENT C, REAL epsilon)
{

	RATIONAL Xc = C.max() + C.wc()/2;
	RATIONAL r = C.Wc() / 2;
	REAL L, R;
	REAL Xprime;
	INTERVALCOMPONENT D;

	L = 4 * r * abs(deriv(P, 1)(Xc));
	R = abs(P(Xc));

	NewTonTimesAll += 1;

	if( choose(L>R, L<R, 2*L< 3*R && 2*R < 3*L) != 2) // when f'(Xc) != 0
	{
		NewTonTimes += 1;

		Xprime = Xc - C.kc * P(Xc)/ deriv(P,1)(Xc);
		OPENINTERVAL II = OPENINTERVAL(Xprime,REAL( C.wc()) / 8 / C.Nc);

		if( softGTest(P, C.kc, II))
		{
			NewTon += 1;

			C.split(2*C.Nc);   //C.depth ++
			INTEGER index = ((Xprime -C.lower)/C.wc()).as_INTEGER();
			if (index == 0)
				index = index + 1;
			if (index == C.size())
				index = index - 1;
			D = INTERVALCOMPONENT();
			D.add(C[index -1]);
			D.add(C[index]);
			D.Nc = C.Nc * C.Nc;
			D.depth = C.depth;
		}
	}
	return D;
}



REAL
RecursiveNewton(int prec, POLYNOMIAL P, OPENINTERVAL D, int nroot)
{
	REAL xc, xprime;
	REAL epsilon = power(2, prec);
	REAL wc; 
	INTEGER nc;
	wc = 2*D.radius;
	nc = 4;

	int depth = 0;
	while(choose(D.radius > epsilon / 2, D.radius < epsilon) == 1)
	{
		xc = D.center + D.radius + wc / 2;
		xprime = xc - nroot * P(xc) / deriv(P,1)(xc);
		D = OPENINTERVAL(xprime, REAL(wc) / 8 / nc);
		wc = wc / (2*nc);
		nc = nc * nc;
		depth += 1;
	}

	MAX_NEWTON = std::max(MAX_NEWTON, depth);
	return D.center;  
}


std::vector< INTERVALCOMPONENT > 
Bisect(POLYNOMIAL P, INTERVALCOMPONENT C)
{
	std::vector< INTERVALCOMPONENT > U;
	RATIONALINTERVAL I;
	int flg;
	bool specialflg = true;

	C.split(2);  //C.depth ++

	for (int i =0; i< C.size(); i++)
	{
		I = C[i];

		if ( softGTest(P,0, I.multiply(RATIONAL(3)/2)) == false)
		{
			flg = 0;
			for(int j=0; j < (int)U.size(); j++)
			{	
				if(adj(U[j], I))
				{	
					U[j].add(I);
					flg = 1;
					break;
				}
			}
			if (flg == 0)	// If there does not exists a adjoint component in U
			{
				INTERVALCOMPONENT T = INTERVALCOMPONENT(I);
				T.depth = C.depth;
				U.push_back(T);
			}
		}
		
	}

	if ((int)U.size() ==1)
		specialflg = false;

	for (int i=0; i<(int)U.size(); i++)
	{
		if (specialflg)
			U[i].Nc = 4;
		else
			if (sqrt(C.Nc) >4)
				U[i].Nc = sqrt(C.Nc);
			else
				U[i].Nc = 4;
	}
	return U;
}

std::vector< std::pair<REAL, int> > 
RealDistinctRootApproximation(int prec, POLYNOMIAL P, int distinct)
{

	int NewTonTimesAll = 0;
	int NewTonTimes = 0;
	int NewTon = 0;
	

	int n = P.degree;

	std::pair<RATIONALINTERVAL, int> tmp;

	RATIONALINTERVAL D = RATIONALINTERVAL(0,1);
	while(!softGTest(P, n, D))
	{
		D = D.multiply(2);
	}
	std::vector< INTERVALCOMPONENT > Q1, Qout;
	
	std::vector< INTERVALCOMPONENT > bisectresult;
	
	INTERVALCOMPONENT C;
	C = INTERVALCOMPONENT(D);
	C.Nc = 4;
	Q1.push_back(C);
	INTERVALCOMPONENT tc, Cprime;
	REAL epsilon = power(2,prec);
	

	int fll;
	while(1)
	{
		if(Q1.empty() && (int) Qout.size() < distinct)
		{
			for(int i=0; i<(int)Qout.size(); i++)
				Q1.push_back(Qout[i]);

			epsilon /= 2;
			Qout.clear();
			fll = 1; 
		}

		if(Q1.empty())
			break;
		
// 		Swap Q1[end -1] <-> Q1[max wc]
		int pi = 0;
		RATIONAL qtmp = 0;
		RATIONAL tmp = 0;
		for(int i=0; i< (int) Q1.size(); i++)
		{
			qtmp = Q1[i].Wc();
			if (qtmp > tmp)
			{
				pi = i;
				tmp = qtmp;
			}
		}
		tc = Q1[pi];
		Q1[pi] = Q1[(int) Q1.size() - 1];
		Q1[(int) Q1.size() -1] = tc;
// 		swap done


		C = Q1.back();
		Q1.pop_back();

		fll = 0;
		RATIONALINTERVAL II = RATIONALINTERVAL(C.Mc(),C.Rc());

		for(int i=0; i< (int) Q1.size(); i++)
		{
			if(intersect(Q1[i], II.multiply(4)))
			{
				fll = 1;
				break;
			}
		}

		if (fll == 0)
		{
			C.kc = softGStar(P, II, P.degree);
			
			if(C.kc >0)
			{
				if(choose(C.Wc() > epsilon, C.Wc() < epsilon*2) == 1)
				{
					
					Cprime= Newton(P,C,epsilon);

					if(Cprime.isempty() == false)
					{
						NewTon += 1;
						Q1.push_back(Cprime);
						continue;
					}
				}
				else if (C.size() < 4)
				{
					Qout.push_back(C);
					continue;
				}
			}
		}

		bisectresult = Bisect(P,C);
		for(int i=0; i<(int)bisectresult.size(); i++)
			Q1.push_back(bisectresult[i]);
	}
	std::sort(Qout.begin(), Qout.end(), IC_LESS());
	
	int m = 0;
	for (int i =0; i<(int) Qout.size(); i++)
		m = std::max(Qout[i].depth, m);
	MAXDEPTH = m;

	std::vector< std::pair<REAL, int> > roots;
	for (int i =0; i<(int) Qout.size(); i++)
		roots.push_back(std::pair<REAL, int> (REAL(Qout[i].Mc()), Qout[i].kc));



	return roots;
}


std::vector<std::pair<REAL, int> >
RealRootFindingNaive(POLYNOMIAL P, int distinct)
{
	std::vector<std::pair<REAL, int> > roots;
	roots.reserve(distinct);

	roots = limit(RealDistinctRootApproximationNaive, P, distinct);

	// computing root separation
	REAL min=100;
	for (int i =0; i<(int)roots.size(); i++)
	{
		if (i != (int)roots.size() -1)
			min = minimum(min, abs(roots[i+1].first - roots[i].first));
	}

	SEPARATION = approx(min, -50);
	return roots;
}

std::vector<std::pair<REAL, int> >
RealRootFindingTrisection(POLYNOMIAL P, int distinct)
{
	std::vector<std::pair<RATIONALINTERVAL, int> > cluster;
	std::vector<std::pair<REAL, int> > roots;
	roots.reserve(distinct);

	REAL f;

	clock_t begin_time = clock();
	cluster = RealDistinctRootSeparationNaive(P,distinct);
	TIME_ROOT_SEPARATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
	

	REAL tmp, tmp1, tmp2;
	POLYNOMIAL Q;

	begin_time = clock();
	for (int i =0; i<(int)cluster.size(); i++)
	{

		if(cluster[i].second % 2 == 0)
		{

			int N = 1;
			while(1)
			{
				if (std::pow(3, std::pow(2, N)) >= RATIONAL(P.degree) / cluster[i].second)
					break;
				N += 1;
			}

			Q = deriv(Giterate(translation(P, cluster[i].first.radius, cluster[i].first.center), N), 1);
			tmp = limit(trisection, Q, -REAL(1),REAL(1));

			for(int j=0; j< N; j++)
				tmp = sqrt(tmp);
			
			tmp1 = tmp*cluster[i].first.radius + cluster[i].first.center;
			tmp2 = 0 - tmp*cluster[i].first.radius + cluster[i].first.center; 
			
			roots.push_back(std::pair<REAL, int>(limit(RootApprox, P, tmp1, tmp2), cluster[i].second));
		}

		else
		{
			sizetype error;
			REAL rt;

			rt = limit(trisection,	P, REAL(cluster[i].first.center - cluster[i].first.radius)
				,REAL(cluster[i].first.center + cluster[i].first.radius));
			roots.push_back(std::pair<REAL, int> (rt, cluster[i].second));
		}

	}
	TIME_ROOT_APPROXIMATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;



	// computing root separation
	REAL min=100;
	for (int i =0; i<(int)roots.size(); i++)
	{
		if (i != (int)roots.size() -1)
			min = minimum(min, abs(roots[i+1].first - roots[i].first));
	}

	SEPARATION = approx(min, -50);
	return roots;
}



std::vector<std::pair<REAL, int> >
RealRootFindingNewton(POLYNOMIAL P, int distinct)
{

	std::vector<std::pair<REAL, int> > roots;
	roots.reserve(distinct);

	roots = limit(RealDistinctRootApproximation, P, distinct);

	// computing root separation
	REAL min=100;
	for (int i =0; i<(int)roots.size(); i++)
	{
		if (i != (int)roots.size() -1)
			min = minimum(min, abs(roots[i+1].first - roots[i].first));
	}

	SEPARATION = approx(min, -50);
	return roots;
}

std::vector<std::pair<REAL, int> >
RealRootFindingCombined(POLYNOMIAL P, int distinct)
{		
	std::vector<std::pair<RATIONALINTERVAL, int> > cluster;
	std::vector<std::pair<REAL, int> > roots;
	roots.reserve(distinct);

	clock_t begin_time = clock();
	cluster = RealDistinctRootSeparationNaive(P,distinct);
	TIME_ROOT_SEPARATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	REAL tmp, tmp1, tmp2;
	REAL cent, rad;
	REAL rapprox;
	int approxnum;
	POLYNOMIAL Q;


	begin_time = clock();
	for (int i =0; i<(int)cluster.size(); i++)
	{
		if(cluster[i].second % 2 == 0)
		{
			int N = 1;
			while(1)
			{
				if (std::pow(3, std::pow(2, N)) >= RATIONAL(3*P.degree) / cluster[i].second)
					break;
				N += 1;
			}

			Q = deriv(Giterate(translation(P, cluster[i].first.radius, cluster[i].first.center), N), 1);
			
			rapprox =  REAL(1) / (power(2,19)*P.degree*P.degree*4);

			approxnum = (int)((log(rapprox)/log(2)).as_INTEGER()) - 2;
			rad = power(2, approxnum);
			cent = trisection(approxnum, Q, REAL(-1), REAL(1));
			

			tmp = limit(RecursiveNewton, Q, OPENINTERVAL(cent, rad), cluster[i].second-1);

			for(int j=0; j< N; j++)
				tmp = sqrt(tmp);
			
			tmp1 = tmp*cluster[i].first.radius + cluster[i].first.center;
			tmp2 = 0 - tmp*cluster[i].first.radius + cluster[i].first.center; 
			roots.push_back(std::pair<REAL, int>(limit(RootApprox, P, tmp1, tmp2), cluster[i].second));
		}
		else
		{
			rapprox =  cluster[i].first.radius / (power(2,19)*P.degree*P.degree*4);
			approxnum = (int)((log(rapprox)/log(2)).as_INTEGER()) - 2;
			rad = power(2, approxnum);
			cent = trisection(approxnum, P, REAL(cluster[i].first.center - cluster[i].first.radius), REAL(cluster[i].first.center + cluster[i].first.radius));
			
			roots.push_back(std::pair<REAL, int> (limit(RecursiveNewton, P, OPENINTERVAL(cent, rad), cluster[i].second), cluster[i].second));
		}
	}
	TIME_ROOT_APPROXIMATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;


	REAL min=100;
	for (int i =0; i<(int)roots.size(); i++)
	{
		if (i != (int)roots.size() -1)
			min = minimum(min, abs(roots[i+1].first - roots[i].first));
	}
	SEPARATION = approx(min,-50);
	return roots;

}		

std::vector<std::pair<REAL, int> >
RealRootFindingBounded(POLYNOMIAL P, int distinct)
{		
	std::vector<std::pair<RATIONALINTERVAL, int> > cluster;
	std::vector<std::pair<REAL, int> > roots;
	roots.reserve(distinct);

	clock_t begin_time = clock();
	cluster = RealDistinctRootSeparationBounded(P,distinct);
	TIME_ROOT_SEPARATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	REAL tmp, tmp1, tmp2;
	REAL cent, rad;
	REAL rapprox;
	int approxnum;
	POLYNOMIAL Q;


	begin_time = clock();
	for (int i =0; i<(int)cluster.size(); i++)
	{
		if(cluster[i].second % 2 == 0)
		{
			int N = 1;
			while(1)
			{
				if (std::pow(3, std::pow(2, N)) >= RATIONAL(3*P.degree) / cluster[i].second)
					break;
				N += 1;
			}

			Q = deriv(Giterate(translation(P, cluster[i].first.radius, cluster[i].first.center), N), 1);
			
			rapprox =  REAL(1) / (power(2,19)*P.degree*P.degree*4);

			approxnum = (int)((log(rapprox)/log(2)).as_INTEGER()) - 2;
			rad = power(2, approxnum);
			cent = trisection(approxnum, Q, REAL(-1), REAL(1));
			

			tmp = limit(RecursiveNewton, Q, OPENINTERVAL(cent, rad), cluster[i].second-1);

			for(int j=0; j< N; j++)
				tmp = sqrt(tmp);
			
			tmp1 = tmp*cluster[i].first.radius + cluster[i].first.center;
			tmp2 = 0 - tmp*cluster[i].first.radius + cluster[i].first.center; 
			roots.push_back(std::pair<REAL, int>(limit(RootApprox, P, tmp1, tmp2), cluster[i].second));
		}
		else
		{
			rapprox =  cluster[i].first.radius / (power(2,19)*P.degree*P.degree*4);
			approxnum = (int)((log(rapprox)/log(2)).as_INTEGER()) - 2;
			rad = power(2, approxnum);
			cent = trisection(approxnum, P, REAL(cluster[i].first.center - cluster[i].first.radius), REAL(cluster[i].first.center + cluster[i].first.radius));
			
			roots.push_back(std::pair<REAL, int> (limit(RecursiveNewton, P, OPENINTERVAL(cent, rad), cluster[i].second), cluster[i].second));
		}
	}
	TIME_ROOT_APPROXIMATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;


	REAL min=100;
	for (int i =0; i<(int)roots.size(); i++)
	{
		if (i != (int)roots.size() -1)
			min = minimum(min, abs(roots[i+1].first - roots[i].first));
	}
	SEPARATION = approx(min,-50);
	return roots;

}		
