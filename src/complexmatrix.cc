#include "../include/intervals.h"
#include "../include/utilities.h"
#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include <stdexcept>

#include <stack>
#include <queue>
#include <utility>
#include <vector>
#include <cstdlib>

using namespace iRRAM;

/*
REALMATRIX.cc -- routines for REALMATRIX class
 
Copyright (C) 2001-2004 Norbert Mueller
 
This file is part of the iRRAM Library.
 
The iRRAM Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Library General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at your
option) any later version.
 
The iRRAM Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
License for more details.
 
You should have received a copy of the GNU Library General Public License
along with the iRRAM Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. 
*/

// namespace iRRAM {
/**************************************************************************/
#define ELEMENT(x,i,j) x.values[(i)*x.maxcolumn+j]

COMPLEXMATRIX::COMPLEXMATRIX(unsigned int rows,unsigned int columns)
{
//fprintf(stderr," +- Kreiere 1m: %x\n",this);
  if ( (rows == 0)||  (columns == 0 ) ) { 
     fprintf(stderr,"Error in dimensioning real matrix of size [%d,%d]\n",
         rows,columns);
     exit(1);
     }
  maxrow=rows;maxcolumn=columns;
  values=new COMPLEX[rows*columns];
}

COMPLEXMATRIX eye (unsigned int rows) {
  COMPLEXMATRIX prod(rows,rows); REAL one=1;
  for (unsigned int i=0;i<rows;i++)ELEMENT(prod,i,i)=one;
  return prod;
}

COMPLEXMATRIX zeroes  (unsigned int rows,
                    unsigned int columns) {
  REAL zero=0;
  COMPLEXMATRIX result(rows,columns);
  for (unsigned int i=0;i<rows;i++) 
  for (unsigned int j=0;j<columns;j++) 
	ELEMENT(result,i,j)=zero;
  return result;
}

COMPLEXMATRIX ones    (unsigned int rows,
                    unsigned int columns) {
  COMPLEXMATRIX result(rows,columns);
  REAL one=1;
  for (unsigned int i=0;i<rows;i++) 
  for (unsigned int j=0;j<columns;j++) 
	ELEMENT(result,i,j)=one;
  return result;
}

COMPLEXMATRIX& COMPLEXMATRIX::operator = (const COMPLEXMATRIX& y) {
  if ( maxrow * maxcolumn != y.maxrow*y.maxcolumn) {
     delete []values;
     values=new COMPLEX[y.maxrow*y.maxcolumn];
     }
  maxrow = y.maxrow;
  maxcolumn = y.maxcolumn;
  unsigned int size=maxrow*maxcolumn;
  for (unsigned int i=0;i<size;i++) values[i]=y.values[i];
  return(*this);
}

COMPLEXMATRIX& COMPLEXMATRIX::operator = (const REALMATRIX& y) {
  if ( maxrow * maxcolumn != y.maxrow*y.maxcolumn) {
     delete []values;
     values=new COMPLEX[y.maxrow*y.maxcolumn];
     }
  maxrow = y.maxrow;
  maxcolumn = y.maxcolumn;
  unsigned int size=maxrow*maxcolumn;
  for (unsigned int i=0;i<size;i++) values[i]=y.values[i];
  return(*this);
}



COMPLEXMATRIX operator + (const COMPLEXMATRIX& x, const COMPLEXMATRIX& y) {
  COMPLEXMATRIX sum(x.maxrow,x.maxcolumn);
  
  if ( x.maxrow != y.maxrow || x.maxcolumn != y.maxcolumn) { 
     fprintf(stderr,"Error in adding real matrices of different sizes \n");
     exit(1);
     }
  unsigned int size=x.maxrow*x.maxcolumn;
  for (unsigned int i=0; i<size;i++) sum.values[i]=x.values[i]+y.values[i];
  return sum;
}

COMPLEXMATRIX operator - (const COMPLEXMATRIX& x, const COMPLEXMATRIX& y) {
  COMPLEXMATRIX sum(x.maxrow,x.maxcolumn);
  
  if ( x.maxrow != y.maxrow || x.maxcolumn != y.maxcolumn) { 
     fprintf(stderr,"Error in adding real matrices of different sizes \n");
     exit(1);
     }
  unsigned int size=x.maxrow*x.maxcolumn;
  for (unsigned int i=0; i<size;i++) sum.values[i]=x.values[i]-y.values[i];
  return sum;
}

COMPLEXMATRIX operator * (const COMPLEXMATRIX& x ,const COMPLEXMATRIX& y) {
  COMPLEXMATRIX prod(x.maxrow,y.maxcolumn);
  
  if ( x.maxcolumn != y.maxrow) { 
     fprintf(stderr,"Error in multiplying real matrices of different sizes \n");
     exit(1);
     }
  for (unsigned int i=0; i<x.maxrow;i++) {
  for (unsigned int j=0; j<y.maxcolumn;j++) {
    ELEMENT(prod,i,j)=ELEMENT(x,i,0)*ELEMENT(y,0,j);
    for (unsigned int k=1; k<x.maxcolumn;k++) {
      ELEMENT(prod,i,j)=ELEMENT(prod,i,j)+ELEMENT(x,i,k)*ELEMENT(y,k,j);
    } } }
  return prod;
}


COMPLEXMATRIX::~COMPLEXMATRIX() {
//fprintf(stderr," +- Destructor 1m: %x\n",this);
  delete []values;
}

COMPLEXMATRIX::COMPLEXMATRIX(){
  maxrow=0;maxcolumn=0;values=NULL;
}

COMPLEXMATRIX::COMPLEXMATRIX(const COMPLEXMATRIX& y){
  maxrow=y.maxrow;maxcolumn=y.maxcolumn;
  unsigned int size=maxrow*maxcolumn;
  values = new COMPLEX[size];
  for (unsigned int i=0; i < size; i++) values[i]=y.values[i];
}

COMPLEX&  COMPLEXMATRIX::element (unsigned int i, unsigned int j) const {
  if ((i>= maxrow) || (j>= maxcolumn) ) {
     fprintf(stderr,"Illegal indices [%d,%d] for real matrix [%d,%d]\n",
           i,j,maxrow,maxcolumn);
     exit(1);
     }
  return values[i*maxcolumn+j];
}

COMPLEX&  COMPLEXMATRIX::operator () (unsigned int i, unsigned int j) const {
  if ((i>= maxrow) || (j>= maxcolumn) ) {
     fprintf(stderr,"Illegal indices [%d,%d] for real matrix [%d,%d]\n",
           i,j,maxrow,maxcolumn);
     exit(1);
     }
  return values[i*maxcolumn+j];
}


COMPLEXMATRIX operator * (const COMPLEXMATRIX& x, const COMPLEX& y) {
  COMPLEXMATRIX prod(x.maxrow,x.maxcolumn); 
  for (unsigned int i=0; i<x.maxcolumn;i++) {
  for (unsigned int j=0; j<x.maxrow;j++) {
      ELEMENT(prod,i,j)=ELEMENT(x,i,j)*y;
    } }
  return prod;
}

COMPLEXMATRIX operator / (const COMPLEXMATRIX& x, const COMPLEX& y) {
  COMPLEXMATRIX prod(x.maxrow,x.maxcolumn); 
  for (unsigned int i=0; i<x.maxcolumn;i++) {
  for (unsigned int j=0; j<x.maxrow;j++) {
      ELEMENT(prod,i,j)=ELEMENT(x,i,j)/y;
    } }
  return prod;
}

// int bound (const REALMATRIX& x, const int k){
//   for (unsigned int i=0;i<x.maxcolumn;i++) 
//   for (unsigned int j=0;j<x.maxrow;j++) 
// 	if ( ! bound(ELEMENT(x,i,j),k) ) return 0;
//   return 1;
// }

// void REALMATRIX::adderror (sizetype error)
// { 
//   for (unsigned int i=0;i<(*this).maxcolumn;i++) 
//   for (unsigned int j=0;j<(*this).maxrow;j++) 
//       ELEMENT((*this),i,j).adderror(error);
// }

// void REALMATRIX::seterror (sizetype error)
// { 
//   for (unsigned int i=0;i<(*this).maxcolumn;i++) 
//   for (unsigned int j=0;j<(*this).maxrow;j++) 
//       ELEMENT((*this),i,j).seterror(error);
// }

// void REALMATRIX::geterror (sizetype& error) const
// {
//   unsigned int i,j;
//   sizetype lerror; 
//   ELEMENT((*this),0,0).geterror(error); 
//   for (i=0;i<(*this).maxcolumn;i++) 
//   for (j=0;j<(*this).maxrow;j++) {
//       ELEMENT((*this),i,j).geterror(lerror);
//       sizetype_max(error,error,lerror);
//   }
// }


// REAL       maxnorm (const REALMATRIX& x){
//   REAL m=abs(x.values[0]);
//   unsigned int size=x.maxrow*x.maxcolumn;
//   for (unsigned int i=1; i<size;i++) m=maximum(m,abs(x.values[i]));
//   return m;
// }

// REAL       rowsumnorm (const REALMATRIX& x){
//   REAL m=0;
//   REAL sum;
//   for (unsigned int i=0;i<x.maxrow;i++) {
//     sum=0;
//     for (unsigned int j=0;j<x.maxcolumn;j++) sum=sum+abs(x(i,j));
//     m=maximum(m,sum);
//   }
//   return m;
// }

// REAL       colsumnorm (const REALMATRIX& x){
//   REAL m=0;
//   REAL sum;
//   for (unsigned int j=0;j<x.maxcolumn;j++) {
//     sum=0;
//     for (unsigned int i=0;i<x.maxrow;i++) sum=sum+abs(x(i,j));
//     m=maximum(m,sum);
//   }
//   return m;
// }
// } // namespace iRRAM