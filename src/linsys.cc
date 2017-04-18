#include "../include/linsys.h"
#include "../include/polynomialroot.h"
#include "../include/complexmatrix.h"
#include "../include/utilities.h"
#include "iRRAM/lib.h"
#include "iRRAM/core.h"

#include <iostream>
#include <utility>
#include <cstdlib>
#include <sstream>
#include <time.h>
#include <random>
using namespace iRRAM;



REAL trace(REALMATRIX mat)
{
	REAL tr = 0;
	for (int i=0; i<(int)mat.maxrow; i++)
		tr += mat(i,i);
	return tr;
}


REAL trace(COMPLEXMATRIX mat)
{
	COMPLEX tr = 0;
	for (int i=0; i<(int)mat.maxrow; i++)
		tr += mat(i,i);
	return tr;
}

REALMATRIX power(REALMATRIX mat, int k)
{
	REALMATRIX tmp = mat;
	for(int i=0; i<k; i++)
		tmp = tmp * mat;
	return tmp;

}



void print(REALMATRIX M)
{
	for (int i=0; i<(int)M.maxrow; i++)
	{
		cout<< "[ ";
		for (int j=0; j<(int)M.maxcolumn; j++)
		{
			cout <<  M(i,j) << " ";
		}
		cout <<"] \n";
	}
	cout<<"\n";
}



REALMATRIX concat(REALMATRIX A, REALMATRIX B)
{
	if(A.maxrow == 0)
		return B;
	if(B.maxrow == 0)
		return A;
	
	REALMATRIX M(A.maxrow, A.maxcolumn + B.maxcolumn);
	for(int r = 0; r<(int)A.maxrow; r++)
	{
		for(int c = 0; c<(int)A.maxcolumn+(int)B.maxcolumn; c++)
		{
			if (c<(int)A.maxcolumn)
			{
				M(r,c) = A(r,c);
			}
			else
			{
				M(r,c) = B(r,c-A.maxcolumn);
			}
		}
	}
	return M;
}

COMPLEXMATRIX concat(COMPLEXMATRIX A, COMPLEXMATRIX B)
{
	if(A.maxrow == 0)
		return B;
	if(B.maxrow == 0)
		return A;
	
	COMPLEXMATRIX M(A.maxrow, A.maxcolumn + B.maxcolumn);
	for(int r = 0; r<(int)A.maxrow; r++)
	{
		for(int c = 0; c<(int)A.maxcolumn+(int)B.maxcolumn; c++)
		{
			if (c<(int)A.maxcolumn)
			{
				M(r,c) = A(r,c);
			}
			else
			{
				M(r,c) = B(r,c-A.maxcolumn);
			}
		}
	}
	return M;
}


REAL inner(REALMATRIX u, REALMATRIX v)
{
	assert(u.maxrow == v.maxrow);
	assert(u.maxcolumn == 1);
	assert(v.maxcolumn == 1);

	REAL sum = 0;
	for(int i=0; i<(int)u.maxrow; i++)
	{
		sum += u(i,0) * v(i,0);
	}
	return sum;
}

REALMATRIX colVector(REALMATRIX u, int idx)
{
	REALMATRIX v = REALMATRIX(u.maxrow, 1);
	for (int i=0; i<(int) u.maxrow; i++)
	{
		v(i,0) = u(i,idx);
	}
	return v;

}

REALMATRIX colVector(REALMATRIX u, int r, int c)
{
	REALMATRIX v = REALMATRIX(u.maxrow - r, 1);
	for (int i=r; i<(int) u.maxrow; i++)
	{
		v(i,0) = u(i,c);
	}
	return v;

}

REALMATRIX projection(REALMATRIX u, REALMATRIX v)
{
	return u * (inner(v,u)/inner(u,u));
}

// Orthonormalize columnvectors of M
REALMATRIX GramSchmidt(REALMATRIX M)
{
	REALMATRIX Q;
	REALMATRIX v;
	for(int i=0; i<(int) M.maxcolumn; i++)
	{
		v = colVector(M,i);
		for(int j=0; j < i; j++)
		{
			v = v - projection(colVector(Q, j), colVector(M, i)); 
		}
		v = v /sqrt(inner(v,v));
		Q = concat(Q, v);
	}
	return Q;
}


COMPLEXMATRIX inner(COMPLEXMATRIX u, COMPLEXMATRIX v)
{
	assert(u.maxrow == v.maxrow);
	assert(u.maxcolumn == 1);
	assert(v.maxcolumn == 1);

	REAL sum = 0;
	for(int i=0; i<(int)u.maxrow; i++)
	{
		sum += u(i,0) * v(i,0);
	}
	return sum;
}

COMPLEXMATRIX colVector(COMPLEXMATRIX u, int idx)
{
	COMPLEXMATRIX v = COMPLEXMATRIX(u.maxrow, 1);
	for (int i=0; i<(int) u.maxrow; i++)
	{
		v(i,0) = u(i,idx);
	}
	return v;

}

REALMATRIX colVector(COMPLEXMATRIX u, int r, int c)
{
	COMPLEXMATRIX v = COMPLEXMATRIX(u.maxrow - r, 1);
	for (int i=r; i<(int) u.maxrow; i++)
	{
		v(i,0) = u(i,c);
	}
	return v;

}

COMPLEXMATRIX projection(COMPLEXMATRIX u, COMPLEXMATRIX v)
{
	return u * (inner(v,u)/inner(u,u));
}

// Orthonormalize columnvectors of M
COMPLEXMATRIX GramSchmidt(COMPLEXMATRIX M)
{
	COMPLEXMATRIX Q;
	COMPLEXMATRIX v;
	for(int i=0; i<(int) M.maxcolumn; i++)
	{
		v = colVector(M,i);
		for(int j=0; j < i; j++)
		{
			v = v - projection(colVector(Q, j), colVector(M, i)); 
		}
		v = v /sqrt(inner(v,v));
		Q = concat(Q, v);
	}
	return Q;
}

REALMATRIX linearSys(REALMATRIX M, REALMATRIX b)
{
	REALMATRIX W = M;
	REALMATRIX w = b;
	REAL tmpfactor;
	REAL tmp = 0;
	REAL tmpV;
	for(int i=0; i<(int)W.maxrow; i++)
	{
		// finding half maximum index
		tmp = 0;
		for(int t = i; t<(int)W.maxrow; t++)
		{
			tmp = maximum(abs(W(t,i)), tmp);
		}

		for(int t = i; t<(int)W.maxrow; t++)
		{
			if(choose(abs(W(t,i)) > tmp / 2, abs(W(t,i)) < tmp) == 1)
			{
				for(int tc=0; tc<(int)W.maxcolumn; tc++)
				{
					tmpV = W(i,tc);
					W(i,tc) = W(t,tc);
					W(t,tc) = tmpV;
				}
				for(int tc=0; tc<(int)w.maxcolumn; tc++)
				{
					tmpV = w(i,tc);
					w(i,tc) = w(t,tc);
					w(t,tc) = tmpV;
				}
				
				tmpfactor = 1/W(i,i);
				
				for(int tc=0; tc<(int)W.maxcolumn; tc++)
				{
					W(i,tc) = W(i,tc)*tmpfactor;
				}
				for(int tc=0; tc<(int)w.maxcolumn; tc++)
				{
					w(i,tc) = w(i,tc)*tmpfactor;
				}
				
				for(int j=0; j<(int)W.maxrow; j++)
				{	
					if(i!= j)
					{
						tmpfactor = 0-W(j,i)/W(i,i);

						for(int tc=0; tc<(int)W.maxcolumn; tc++)
						{
							W(j,tc) = W(j,tc)+ W(i,tc)*tmpfactor;
						}
						for(int tc=0; tc<(int)w.maxcolumn; tc++)
						{
							w(j,tc) = w(j,tc)+ w(i,tc)*tmpfactor;
						}
						
					}
				}
			break;
			}
		}
	}
	return w;
}



POLYNOMIAL traceFormulae(REALMATRIX M)
{
	REAL coef[M.maxcolumn +1];
	int n = (int)M.maxcolumn;
	coef[n] = 1;
	REAL T[n];
	REALMATRIX tmp = M;
	for(int k=0; k< n; k++)
	{
		T[k] = trace(tmp);
		tmp =  tmp * M;
	}
	for(int k=0; k<n; k++)
	{
		for(int i=0; i< k+1; i++)
		{
			coef[n-k-1] += T[k-i]*coef[n-i];
		}
		coef[n-k-1] *= -1;
		coef[n-k-1] /= (k+1);
	}

	POLYNOMIAL P = POLYNOMIAL(n, coef);
	return P;
}

POLYNOMIAL HermitianTraceFormulae(COMPLEXMATRIX M)
{
	REAL coef[M.maxcolumn +1];
	int n = (int)M.maxcolumn;
	coef[n] = 1;
	REAL T[n];
	COMPLEXMATRIX tmp = M;
	for(int k=0; k< n; k++)
	{
		T[k] = real(trace(tmp));
		tmp =  tmp * M;
	}
	for(int k=0; k<n; k++)
	{
		for(int i=0; i< k+1; i++)
		{
			coef[n-k-1] += T[k-i]*coef[n-i];
		}
		coef[n-k-1] *= -1;
		coef[n-k-1] /= (k+1);
	}

	POLYNOMIAL P = POLYNOMIAL(n, coef);
	return P;
}



REAL determinant(REALMATRIX M)
{
	REALMATRIX W = M;
	REAL det = 1;
	REAL tmp = 0;
	REAL tmpV;
	for(int i=0; i<(int)W.maxrow; i++)
	{
		tmp = 0;
		// finding half maximum index
		for(int t = i; t<(int)W.maxrow; t++)
			tmp = maximum(abs(W(t,i)), tmp);
		for(int t = i; t<(int)W.maxrow; t++)
		{
			if(choose(abs(W(t,i)) > tmp / 2, abs(W(t,i)) < tmp) == 1)
			{
				for(int tc=0; tc<(int)W.maxcolumn; tc++)
				{
					tmpV = W(i,tc);
					W(i,tc) = W(t,tc);
					W(t,tc) = tmpV;
				}				

				det *= -1;
				
				for(int j=i+1; j<(int)W.maxrow; j++)
				{
					tmpV = W(j,i)/W(i,i);
					for(int tc=0; tc<(int)W.maxcolumn; tc++)
					{
						W(j,tc) = W(j,tc) - W(i,tc)*tmpV;
					}
				}
			break;
			}
		}
	}
	for(int i=0; i<(int)M.maxrow; i++)
		det *= W(i,i);
	return det;
}

REAL rootBound(REALMATRIX M)
{
	REAL tmp = 0;
	for(int r=0; r<(int)M.maxrow; r++)
		for(int c=0; c<(int)M.maxcolumn; c++)
			tmp = maximum(abs(M(r,c)), tmp);
	tmp = maximum(tmp,2);
	REAL maxCoef = power(2, ((int)M.maxrow / 2)*(iRRAM::log((int)M.maxrow)/iRRAM::log(REAL(2)) + iRRAM::log(tmp*tmp)/iRRAM::log(REAL(2)) + REAL(0.21163175)));
	REAL raucheBound = 1 + maxCoef;
	return raucheBound;
}

POLYNOMIAL interpolatingMethod(REALMATRIX M)
{
	std::vector<std::pair<REAL, REAL> > interpole((int)M.maxcolumn);
	REAL rb = rootBound(M);
	REAL xi;
	REAL yi;
	REALMATRIX tmpM = M;

	REALMATRIX A = REALMATRIX((int)M.maxrow+1, (int)M.maxrow+1);
	REALMATRIX B = REALMATRIX((int)M.maxrow+1, 1);
	for(int i=0; i< (int)M.maxcolumn+1; i++)
	{
		xi = rb + i;
		for(int k=0; k<(int)M.maxcolumn+1; k++)
			A(i,k) = power(xi,k);

		tmpM = M;
		for(int j=0; j<(int)M.maxcolumn; j++)
			tmpM(j,j) = M(j,j) - xi;
		yi = determinant(tmpM);

		B(i,0) = yi;
	}

	REALMATRIX COEF = linearSys(A,B);

	REAL c[(int)M.maxrow+1];
	for(int i=0; i<(int)M.maxrow+1; i++)
		c[i] = COEF(i,0);
	POLYNOMIAL Q = POLYNOMIAL((int)M.maxrow, c);

	return Q;
}

POLYNOMIAL charPoly(REALMATRIX M)
{
	return traceFormulae(M);
}


POLYNOMIAL HermitianCharPoly(COMPLEXMATRIX M)
{
	return HermitianTraceFormulae(M);
}


POLYNOMIAL charPoly(REALMATRIX M, bool t)
{
	if(t)
		return traceFormulae(M);
	else
		return interpolatingMethod(M);
}





REALMATRIX inv(REALMATRIX M)
{

	REALMATRIX INV = REALMATRIX(M.maxcolumn,M.maxcolumn);
	for (int i = 0 ; i < (int)M.maxcolumn; i ++)
		INV(i,i) = 1;
	REALMATRIX ins = linearSys(M,INV);
	return ins;
}


REALMATRIX VSimilar(REALMATRIX M)
{
	REALMATRIX V = REALMATRIX(M.maxrow, M.maxcolumn);
	for(int r = 0; r<(int)M.maxrow; r ++)
	{
		for(int c = 0; c<(int)M.maxcolumn; c++)
		{
			V(r,c) = power(REAL(r+1), c);
		}
	}
	return (inv(V)*M)*V;
}


REALMATRIX eigenVector(REALMATRIX A, REAL eigenValue, int nulldim)
{

	REAL max, pivotvalue, scale, temp;
	int i,j,k, temp_index, flg, ROW, COL;
	int pivoti = 0;
	int pivotj = 0;
	ROW =  (int)A.maxrow;
	COL = (int)A.maxcolumn;
	int rank = ROW - nulldim;
	
	REALMATRIX solutionMatrix = REALMATRIX(rank, COL - rank);
	REAL tmpV;

	for (int i=0; i<COL; i++)
	{
		A(i,i) = A(i,i)-eigenValue;
	}

	int record[COL];
    for (int i=0; i<COL; i++) record[i] = i;


	for (i=0; i < rank; i++)
	{


		max = 0;
		for(int ti = i; ti < ROW; ti ++)
		{
			for(int tj = i; tj <COL; tj ++)
			{
				max = maximum(abs(A(ti,tj)), max);
			}
		}

		flg = 0;
    	for (k=i; k<COL; k++)
    	{
    		for (j=i; j<ROW; j++)
    		{

    			if(1==choose(abs(A(j,k)) > max / 2, abs(A(j,k)) < max))
    			{
					pivotvalue = A(j,k);
    				pivoti = j;
    				pivotj = k;	
    				flg = 1;
    				break;
    			}
    		}
    		if (flg==1) break;
    	}

    	if (pivoti != i)
    	{
			for(int tc=0; tc<(int)A.maxcolumn; tc++)
			{
				tmpV = A(i,tc);
				A(i,tc) = A(pivoti,tc);
				A(pivoti,tc) = tmpV;
			}
    	}

		if (pivotj != i)
		{
			
			for(int tc=0; tc<(int)A.maxrow; tc++)
			{
				tmpV = A(tc,i);
				A(tc,i) = A(tc,pivotj);
				A(tc,pivotj) = tmpV;
			}
            temp_index = record[i];
            record[i] = record[pivotj];
            record[pivotj] = temp_index;
		}

		for(int tc=0; tc<(int)A.maxcolumn; tc++)
		{
			A(i,tc) = A(i,tc) / pivotvalue;
		}

		for (j=0; j<ROW; j++) 
		{
			if (j == i)
				continue;
			scale = A(j,i);
			

			for(int tc=0; tc<(int)A.maxcolumn; tc++)
			{
				A(j,tc) = A(j,tc) -scale*A(i,tc);
			}
		}
	}


	// copy ** to solution matrix. 
	for(i = 0; i<rank; i++)
	{
		for(j=0; j<COL - rank; j++)
		{
			solutionMatrix(i,j) = A(i,j + rank);
		}
	}


	REALMATRIX B = REALMATRIX(ROW,nulldim);
	for(int i=0; i<COL - rank; i++)
	{
		for(int j=0; j<rank; j++)
		{
			B(record[j],i) =  solutionMatrix(j,i);
		}
		for (int k=0; k<COL - rank; k++)
		{
			if(i == k)
				B(record[rank + k], i) =  REAL(-1);
			else
				B(record[rank + k], i) =  REAL(0);
		}
	}
	B = GramSchmidt(B);

	return B;
}


COMPLEXMATRIX eigenVector(COMPLEXMATRIX A, COMPLEX eigenValue, int nulldim)
{

	REAL max;
	COMPLEX pivotvalue, scale, temp;
	int i,j,k, temp_index, flg, ROW, COL;
	int pivoti = 0;
	int pivotj = 0;
	ROW =  (int)A.maxrow;
	COL = (int)A.maxcolumn;
	int rank = ROW - nulldim;
	
	COMPLEXMATRIX solutionMatrix = COMPLEXMATRIX(rank, COL - rank);
	COMPLEX tmpV;

	for (int i=0; i<COL; i++)
	{
		A(i,i) = A(i,i) - eigenValue;
	}

	int record[COL];
    for (int i=0; i<COL; i++) record[i] = i;


	for (i=0; i < rank; i++)
	{

		max = 0;
		for(int ti = i; ti < ROW; ti ++)
		{
			for(int tj = i; tj <COL; tj ++)
			{
				max = maximum(abs(A(ti,tj)), max);
			}
		}

		flg = 0;
    	for (k=i; k<COL; k++)
    	{
    		for (j=i; j<ROW; j++)
    		{

    			if(1==choose(abs(A(j,k)) > max / 2, abs(A(j,k)) < max))
    			{
					pivotvalue = A(j,k);
    				pivoti = j;
    				pivotj = k;	
    				flg = 1;
    				break;
    			}
    		}
    		if (flg==1) break;
    	}

    	if (pivoti != i)
    	{
			for(int tc=0; tc<(int)A.maxcolumn; tc++)
			{
				tmpV = A(i,tc);
				A(i,tc) = A(pivoti,tc);
				A(pivoti,tc) = tmpV;
			}
    	}

		if (pivotj != i)
		{
			
			for(int tc=0; tc<(int)A.maxrow; tc++)
			{
				tmpV = A(tc,i);
				A(tc,i) = A(tc,pivotj);
				A(tc,pivotj) = tmpV;
			}
            temp_index = record[i];
            record[i] = record[pivotj];
            record[pivotj] = temp_index;
		}

		for(int tc=0; tc<(int)A.maxcolumn; tc++)
		{
			A(i,tc) = A(i,tc) / pivotvalue;
		}

		for (j=0; j<ROW; j++) 
		{
			if (j == i)
				continue;
			scale = A(j,i);
			

			for(int tc=0; tc<(int)A.maxcolumn; tc++)
			{
				A(j,tc) = A(j,tc) - scale*A(i,tc);
			}
		}
	}


	// copy ** to solution matrix. 
	for(i = 0; i<rank; i++)
	{
		for(j=0; j<COL - rank; j++)
		{
			solutionMatrix(i,j) = A(i,j + rank);
		}
	}


	COMPLEXMATRIX B = COMPLEXMATRIX(ROW,nulldim);
	for(int i=0; i<COL - rank; i++)
	{
		for(int j=0; j<rank; j++)
		{
			B(record[j],i) =  solutionMatrix(j,i);
		}
		for (int k=0; k<COL - rank; k++)
		{
			if(i == k)
				B(record[rank + k], i) =  REAL(-1);
			else
				B(record[rank + k], i) =  REAL(0);
		}
	}
	B = GramSchmidt(B);

	return B;
}


REALMATRIX diagonalizeNaive(REALMATRIX M, int d)
{
	assert(M.maxcolumn == M.maxrow);
	
	clock_t begin_time = clock();
	POLYNOMIAL P = traceFormulae(M);
	TIME_CHAR_POLYNOMIAL = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	begin_time = clock();
	std::vector<std::pair<REAL,int> > roots = RealRootFindingNaive(P, d);
 	TIME_ROOT_FINDING =  float( clock () - begin_time ) /  CLOCKS_PER_SEC;

 	
  	REALMATRIX m;
	REALMATRIX Q;
  	REAL eigenVal;

	if (roots[0].second == M.maxcolumn)
	{
		Q =  REALMATRIX(M.maxcolumn, M.maxcolumn);
		for (int i=0; i<M.maxcolumn; i++)
			Q(i,i) = 1;
		return Q;
  	}

 	begin_time = clock();	
	for (int i =0; i<(int)roots.size(); i++)
	{	
		

		m = eigenVector(M, roots[i].first, roots[i].second);
		

		Q = concat(m,Q);
	}
	TIME_GAUSSIAN_ELIMINATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	return Q;
}


REALMATRIX diagonalizeTrisection(REALMATRIX M, int d)
{
	assert(M.maxcolumn == M.maxrow);
	
	clock_t begin_time = clock();
	POLYNOMIAL P = traceFormulae(M);
	TIME_CHAR_POLYNOMIAL = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	begin_time = clock();
	std::vector<std::pair<REAL,int> > roots = RealRootFindingTrisection(P, d);
 	TIME_ROOT_FINDING =  float( clock () - begin_time ) /  CLOCKS_PER_SEC;

 	
  	REALMATRIX m;
	REALMATRIX Q;
  	REAL eigenVal;

	if (roots[0].second == M.maxcolumn)
	{
		Q =  REALMATRIX(M.maxcolumn, M.maxcolumn);
		for (int i=0; i<M.maxcolumn; i++)
			Q(i,i) = 1;
		return Q;
  	}

 	begin_time = clock();	
	for (int i =0; i<(int)roots.size(); i++)
	{	
		

		m = eigenVector(M, roots[i].first, roots[i].second);
		

		Q = concat(m,Q);
	}
	TIME_GAUSSIAN_ELIMINATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	return Q;
}


REALMATRIX diagonalizeNewton(REALMATRIX M, int d)
{
	assert(M.maxcolumn == M.maxrow);

	clock_t begin_time = clock();
	POLYNOMIAL P = traceFormulae(M);
	TIME_CHAR_POLYNOMIAL = float( clock () - begin_time ) /  CLOCKS_PER_SEC;



	begin_time = clock();
	std::vector<std::pair<REAL,int> > roots = RealRootFindingNewton(P, d);
	TIME_ROOT_FINDING =  float( clock () - begin_time ) /  CLOCKS_PER_SEC;
  	

  	REALMATRIX m;
	REALMATRIX Q;
  	REAL eigenVal;

	if (roots[0].second == M.maxcolumn)
	{
		Q =  REALMATRIX(M.maxcolumn, M.maxcolumn);
		for (int i=0; i<M.maxcolumn; i++)
			Q(i,i) = 1;
		return Q;
  	}

	begin_time = clock();
	for (int i =0; i<(int)roots.size(); i++)
	{	
		m = eigenVector(M, roots[i].first, roots[i].second);
		Q = concat(m,Q);
	}
	TIME_GAUSSIAN_ELIMINATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	return Q;
}


REALMATRIX diagonalizeCombined(REALMATRIX M, int d)
{
	assert(M.maxcolumn == M.maxrow);

	clock_t begin_time = clock();
	POLYNOMIAL P = traceFormulae(M);
	TIME_CHAR_POLYNOMIAL = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	begin_time = clock();
	std::vector<std::pair<REAL,int> > roots = RealRootFindingCombined(P, d);
	TIME_ROOT_FINDING =  float( clock () - begin_time ) /  CLOCKS_PER_SEC;
  	

  	REALMATRIX m;
	REALMATRIX Q;
  	REAL eigenVal;

	if (roots[0].second == M.maxcolumn)
	{
		Q =  REALMATRIX(M.maxcolumn, M.maxcolumn);
		for (int i=0; i<M.maxcolumn; i++)
			Q(i,i) = 1;
		return Q;
  	}

 	begin_time = clock();	
	for (int i =0; i<(int)roots.size(); i++)
	{	
		m = eigenVector(M, roots[i].first, roots[i].second);
		Q = concat(m,Q);
	}
	TIME_GAUSSIAN_ELIMINATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	return Q;
}

REALMATRIX diagonalizeBounded(REALMATRIX M, int d)
{
	assert(M.maxcolumn == M.maxrow);

	clock_t begin_time = clock();
	POLYNOMIAL P = traceFormulae(M);
	TIME_CHAR_POLYNOMIAL = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	begin_time = clock();
	std::vector<std::pair<REAL,int> > roots = RealRootFindingBounded(P, d);
	TIME_ROOT_FINDING =  float( clock () - begin_time ) /  CLOCKS_PER_SEC;
  	

  	REALMATRIX m;
	REALMATRIX Q;
  	REAL eigenVal;

	if (roots[0].second == M.maxcolumn)
	{
		Q =  REALMATRIX(M.maxcolumn, M.maxcolumn);
		for (int i=0; i<M.maxcolumn; i++)
			Q(i,i) = 1;
		return Q;
  	}

 	begin_time = clock();	
	for (int i =0; i<(int)roots.size(); i++)
	{	
		m = eigenVector(M, roots[i].first, roots[i].second);
		Q = concat(m,Q);
	}
	TIME_GAUSSIAN_ELIMINATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	return Q;
}

COMPLEXMATRIX HermitianDiagonalizeBounded(COMPLEXMATRIX M, int d)
{
	assert(M.maxcolumn == M.maxrow);

	clock_t begin_time = clock();
	POLYNOMIAL P = HermitiantraceFormulae(M);
	TIME_CHAR_POLYNOMIAL = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	begin_time = clock();
	std::vector<std::pair<REAL,int> > roots = RealRootFindingBounded(P, d);
	TIME_ROOT_FINDING =  float( clock () - begin_time ) /  CLOCKS_PER_SEC;
  	

  	COMPLEXMATRIX m;
	COMPLEXMATRIX Q;
  	REAL eigenVal;

	if (roots[0].second == M.maxcolumn)
	{
		Q =  COMPLEXMATRIX(M.maxcolumn, M.maxcolumn);
		for (int i=0; i<M.maxcolumn; i++)
			Q(i,i) = 1;
		return Q;
  	}

 	begin_time = clock();	
	for (int i =0; i<(int)roots.size(); i++)
	{	
		m = eigenVector(M, roots[i].first, roots[i].second);
		Q = concat(m,Q);
	}
	TIME_GAUSSIAN_ELIMINATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	return Q;
}



REALMATRIX transpose(REALMATRIX M)
{
	REALMATRIX N =REALMATRIX(M.maxcolumn, M.maxrow);
	for(int i=0; i < (int)M.maxcolumn; i++)
		for(int j=0; j < (int)M.maxrow; j++)
			N(i,j) = M(j,i);
	return N;
}

COMPLEXMATRIX transpose(COMPLEXMATRIX M)
{
	COMPLEXMATRIX N =COMPLEXMATRIX(M.maxcolumn, M.maxrow);
	for(int i=0; i < (int)M.maxcolumn; i++)
		for(int j=0; j < (int)M.maxrow; j++)
			N(i,j) = M(j,i);
	return N;
}


REALMATRIX basisVector(int i, int n)
{
	REALMATRIX M = REALMATRIX(n, 1);
	M(i,0) = 1;
	return M;
}

REALMATRIX subMatrix(REALMATRIX M, int r, int c)
{
	REALMATRIX N = REALMATRIX(M.maxrow-r, M.maxcolumn-c);
	for(int i=r; i <(int)M.maxrow; i++)
		for(int j=c; j<(int)M.maxcolumn; j++)
			N(i-r, j-c) = M(i,j);
	return N;

}


// QR decomposition for a regular matrix X
// return a normalized Q, orthogonal matrix
REALMATRIX QRDecomposition(REALMATRIX X)
{

	int n = (int) X.maxrow;
	REALMATRIX Y = X;
	REAL alpha;
	REALMATRIX u;
	REALMATRIX v;
	REALMATRIX Q[n];
	REALMATRIX x;
	REAL rr[n];
	for (int i=0; i < n-1; i ++)
	{
		x = colVector(Y, 0);
		alpha = sqrt(inner(x, x));
		u = x - alpha * basisVector(0, n - i);
		v = u / sqrt(inner(u,u));
		Q[i] = eye(n-i) - 2 * v * transpose(v);
		Y = Q[i]*Y;
		rr[i] = Y(0,0);
		Y = subMatrix(Y, 1,1);
	}
	rr[n-1] = Y(0,0);

	REALMATRIX T;
	for (int i=0; i< n-1; i++)
	{
		T = eye(n);
		for(int j=i; j<n;j++)
		{
			for(int k=i;k<n;k++)
			{
				T(j,k) = Q[i](k-i,j-i);
			}
		}
		Q[i] = T;
	}
	T = Q[0];
	for (int i=1;i<n-1;i++)
		T = T * Q[i];
	

	REALMATRIX lambda = REALMATRIX(n,n);
	for(int i=0;i<n;i++)
		if(rr[i] > 0)
			lambda(i,i) = 1;
		else
			lambda(i,i) = -1;

	return T*lambda;

}

REALMATRIX HaarMeasure(int n)
{
	
	REALMATRIX X = REALMATRIX(n,n);
	std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,1.0);
     
	for(int i=0; i<n;i++)
		for(int j=0; j<n;j++)
			X(i,j) = distribution(generator);



	REALMATRIX Y = X;
	REAL alpha;
	REALMATRIX u;
	REALMATRIX v;
	REALMATRIX Q[n];
	REALMATRIX x;
	REAL rr[n];
	for (int i=0; i < n-1; i ++)
	{
		x = colVector(Y, 0);
		alpha = sqrt(inner(x, x));
		u = x - alpha * basisVector(0, n - i);
		v = u / sqrt(inner(u,u));
		Q[i] = eye(n-i) - 2 * v * transpose(v);
		Y = Q[i]*Y;
		rr[i] = Y(0,0);
		Y = subMatrix(Y, 1,1);
	}
	rr[n-1] = Y(0,0);

	REALMATRIX T;
	for (int i=0; i< n-1; i++)
	{
		T = eye(n);
		for(int j=i; j<n;j++)
		{
			for(int k=i;k<n;k++)
			{
				T(j,k) = Q[i](k-i,j-i);
			}
		}
		Q[i] = T;
	}
	T = Q[0];
	for (int i=1;i<n-1;i++)
		T = T * Q[i];
	

	REALMATRIX lambda = REALMATRIX(n,n);
	for(int i=0;i<n;i++)
		if(rr[i] > 0)
			lambda(i,i) = 1;
		else
			lambda(i,i) = -1;

	return T*lambda;
}


std::vector<std::pair< REALMATRIX, REAL> > 
diagonalizeNaiveEig(REALMATRIX M, int d)
{
	assert(M.maxcolumn == M.maxrow);
	
	clock_t begin_time = clock();
	POLYNOMIAL P = traceFormulae(M);
	TIME_CHAR_POLYNOMIAL = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	begin_time = clock();
	std::vector<std::pair<REAL,int> > roots = RealRootFindingNaive(P, d);
 	TIME_ROOT_FINDING =  float( clock () - begin_time ) /  CLOCKS_PER_SEC;

 	std::vector<std::pair< REALMATRIX, REAL> > ans; 
  	REALMATRIX m;
	REALMATRIX Q;
  	REAL eigenVal;

	if (roots[0].second == M.maxcolumn)
	{
		Q =  REALMATRIX(M.maxcolumn, M.maxcolumn);
		for (int i=0; i<M.maxcolumn; i++)
			Q(i,i) = 1;
		ans.push_back(std::pair<REALMATRIX, REAL>(Q,roots[0].first));


		return ans;
  	}

 	begin_time = clock();	
	for (int i =0; i<(int)roots.size(); i++)
	{	
		 ans.push_back(std::pair<REALMATRIX, REAL>(eigenVector(M, roots[i].first, roots[i].second),roots[i].first));
	}
	TIME_GAUSSIAN_ELIMINATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	return ans;
}


std::vector<std::pair< REALMATRIX, REAL> > 
diagonalizeTrisectionEig(REALMATRIX M, int d)
{
	assert(M.maxcolumn == M.maxrow);
	
	clock_t begin_time = clock();
	POLYNOMIAL P = traceFormulae(M);
	TIME_CHAR_POLYNOMIAL = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	begin_time = clock();
	std::vector<std::pair<REAL,int> > roots = RealRootFindingTrisection(P, d);
 	TIME_ROOT_FINDING =  float( clock () - begin_time ) /  CLOCKS_PER_SEC;

  	std::vector<std::pair< REALMATRIX, REAL> > ans; 	
  	REALMATRIX m;
	REALMATRIX Q;
  	REAL eigenVal;

	if (roots[0].second == M.maxcolumn)
	{
		Q =  REALMATRIX(M.maxcolumn, M.maxcolumn);
		for (int i=0; i<M.maxcolumn; i++)
			Q(i,i) = 1;
		ans.push_back(std::pair<REALMATRIX, REAL>(Q,roots[0].first));


		return ans;
  	}

 	begin_time = clock();	
	for (int i =0; i<(int)roots.size(); i++)
	{	
		ans.push_back(std::pair<REALMATRIX, REAL>(eigenVector(M, roots[i].first, roots[i].second),roots[i].first));
	}
	TIME_GAUSSIAN_ELIMINATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	return ans;
}


std::vector<std::pair< REALMATRIX, REAL> > 
diagonalizeNewtonEig(REALMATRIX M, int d)
{
	assert(M.maxcolumn == M.maxrow);

	clock_t begin_time = clock();
	POLYNOMIAL P = traceFormulae(M);
	TIME_CHAR_POLYNOMIAL = float( clock () - begin_time ) /  CLOCKS_PER_SEC;



	begin_time = clock();
	std::vector<std::pair<REAL,int> > roots = RealRootFindingNewton(P, d);
	TIME_ROOT_FINDING =  float( clock () - begin_time ) /  CLOCKS_PER_SEC;
  	
  	std::vector<std::pair< REALMATRIX, REAL> > ans; 	
  	REALMATRIX m;
	REALMATRIX Q;
  	REAL eigenVal;

	if (roots[0].second == M.maxcolumn)
	{
		Q =  REALMATRIX(M.maxcolumn, M.maxcolumn);
		for (int i=0; i<M.maxcolumn; i++)
			Q(i,i) = 1;
		ans.push_back(std::pair<REALMATRIX, REAL>(Q,roots[0].first));


		return ans;
  	}

	begin_time = clock();
	for (int i =0; i<(int)roots.size(); i++)
	{	
		ans.push_back(std::pair<REALMATRIX, REAL>(eigenVector(M, roots[i].first, roots[i].second),roots[i].first));
	}
	TIME_GAUSSIAN_ELIMINATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	return ans;
}


std::vector<std::pair< REALMATRIX, REAL> > 
diagonalizeCombinedEig(REALMATRIX M, int d)
{
	assert(M.maxcolumn == M.maxrow);

	clock_t begin_time = clock();
	POLYNOMIAL P = traceFormulae(M);
	TIME_CHAR_POLYNOMIAL = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	begin_time = clock();
	std::vector<std::pair<REAL,int> > roots = RealRootFindingCombined(P, d);
	TIME_ROOT_FINDING =  float( clock () - begin_time ) /  CLOCKS_PER_SEC;
  	
  	std::vector<std::pair< REALMATRIX, REAL> > ans; 	
  	REALMATRIX m;
	REALMATRIX Q;
  	REAL eigenVal;

	if (roots[0].second == M.maxcolumn)
	{
		Q =  REALMATRIX(M.maxcolumn, M.maxcolumn);
		for (int i=0; i<M.maxcolumn; i++)
			Q(i,i) = 1;
	
		ans.push_back(std::pair<REALMATRIX, REAL>(Q,roots[0].first));


		return ans;
  	}

 	begin_time = clock();	
	for (int i =0; i<(int)roots.size(); i++)
	{	
		ans.push_back(std::pair<REALMATRIX, REAL>(eigenVector(M, roots[i].first, roots[i].second),roots[i].first));
	}
	TIME_GAUSSIAN_ELIMINATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	return ans;
}

std::vector<std::pair< REALMATRIX, REAL> > 
diagonalizeBoundedEig(REALMATRIX M, int d)
{
	assert(M.maxcolumn == M.maxrow);

	clock_t begin_time = clock();
	POLYNOMIAL P = traceFormulae(M);
	TIME_CHAR_POLYNOMIAL = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	begin_time = clock();
	std::vector<std::pair<REAL,int> > roots = RealRootFindingBounded(P, d);
	TIME_ROOT_FINDING =  float( clock () - begin_time ) /  CLOCKS_PER_SEC;
  	
  	std::vector<std::pair< REALMATRIX, REAL> > ans; 	
  	REALMATRIX m;
	REALMATRIX Q;
  	REAL eigenVal;

	if (roots[0].second == M.maxcolumn)
	{
		Q =  REALMATRIX(M.maxcolumn, M.maxcolumn);
		for (int i=0; i<M.maxcolumn; i++)
			Q(i,i) = 1;
	
		ans.push_back(std::pair<REALMATRIX, REAL>(Q,roots[0].first));


		return ans;
  	}

 	begin_time = clock();	
	for (int i =0; i<(int)roots.size(); i++)
	{	
		ans.push_back(std::pair<REALMATRIX, REAL>(eigenVector(M, roots[i].first, roots[i].second),roots[i].first));
	}
	TIME_GAUSSIAN_ELIMINATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	return ans;
}


std::vector<std::pair< COMPLEXMATRIX, REAL> > 
HermitianDiagonalizeBoundedEig(COMPLEXMATRIX M, int d)
{
	assert(M.maxcolumn == M.maxrow);

	clock_t begin_time = clock();
	POLYNOMIAL P = HermitianTraceFormulae(M);
	TIME_CHAR_POLYNOMIAL = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	begin_time = clock();
	std::vector<std::pair<REAL,int> > roots = RealRootFindingBounded(P, d);
	TIME_ROOT_FINDING =  float( clock () - begin_time ) /  CLOCKS_PER_SEC;
  	
  	std::vector<std::pair< REALMATRIX, REAL> > ans; 	
  	COMPLEXMATRIX m;
	COMPLEXMATRIX Q;
  	REAL eigenVal;

	if (roots[0].second == M.maxcolumn)
	{
		Q =  COMPLEXMATRIX(M.maxcolumn, M.maxcolumn);
		for (int i=0; i<M.maxcolumn; i++)
			Q(i,i) = 1;
	
		ans.push_back(std::pair<COMPLEXMATRIX, REAL>(Q,roots[0].first));


		return ans;
  	}

 	begin_time = clock();	
	for (int i =0; i<(int)roots.size(); i++)
	{	
		ans.push_back(std::pair<COMPLEXMATRIX, REAL>(eigenVector(M, roots[i].first, roots[i].second),roots[i].first));
	}
	TIME_GAUSSIAN_ELIMINATION = float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	return ans;
}