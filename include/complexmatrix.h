#ifndef INTERVALS_H
#define INTERVALS_H

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include <utility>
using namespace iRRAM;

#ifndef COMPLEXMATRIX_H
#define COMPLEXMATRIX_H

class COMPLEXMATRIX
{
public:
// Constructors: -------------------------------

COMPLEXMATRIX(unsigned int rows,unsigned int columns);
COMPLEXMATRIX();
COMPLEXMATRIX(const COMPLEXMATRIX& y);

// Copy Constructor: ---------------------------

COMPLEXMATRIX&   operator = (const COMPLEXMATRIX& y);
COMPLEXMATRIX&   operator = (const REALMATRIX& y);


// Destructor: ---------------------------------

~REALMATRIX();

// Access to matrix elements: ------------------

REAL&   operator ()           (unsigned int row,
                               unsigned int column) const;

REAL&   element               (unsigned int row,
                               unsigned int column) const;

// Standard Arithmetic: ------------------------

friend COMPLEXMATRIX  operator + (const COMPLEXMATRIX& x,
                               const COMPLEXMATRIX& y);
friend COMPLEXMATRIX  operator - (const COMPLEXMATRIX& x,
                               const COMPLEXMATRIX& y);
friend COMPLEXMATRIX  operator * (const COMPLEXMATRIX& x,
                               const COMPLEXMATRIX& y);
// friend COMPLEXMATRIX  operator / (const COMPLEXMATRIX& x,
//                                const COMPLEXMATRIX& y);

// Arithmetic with Scalar: ---------------------
friend COMPLEXMATRIX  operator * (const COMPLEXMATRIX& x,
                               const COMPLEX& y);
inline friend COMPLEXMATRIX operator * (const COMPLEX& y,
                               const COMPLEXMATRIX& x)
       {return x*y;};

friend COMPLEXMATRIX  operator / (const COMPLEXMATRIX& x,
                               const COMPLEX& y);

// Information on Dimensions: ------------------

friend inline unsigned int rows (const COMPLEXMATRIX& x)
                {return x.maxrow;};
friend inline unsigned int columns (const COMPLEXMATRIX& x)
                {return x.maxcolumn;};

// friend int bound (const COMPLEXMATRIX& x, 
//                   const int k);

// Linear Algebra: -----------------------------

friend COMPLEXMATRIX eye     (unsigned int rows);
friend COMPLEXMATRIX zeroes  (unsigned int rows,
                           unsigned int columns);
friend COMPLEXMATRIX ones    (unsigned int rows,
                           unsigned int columns);


// friend REAL       maxnorm (const REALMATRIX& x);
// friend REAL       rowsumnorm (const REALMATRIX& x);
// friend REAL       colsumnorm (const REALMATRIX& x);

// friend REALMATRIX solve (
//        REALMATRIX& lside,
//        REALMATRIX& rside,
//        int use_pivot);

// // Limits: --------------------------

// friend REALMATRIX limit_lip (REALMATRIX f(int, const REALMATRIX&),
//                            int lip,
//                            bool on_domain(const REALMATRIX&),
//                            const REALMATRIX& x);


public:
COMPLEX*  values; 
unsigned int    maxrow,maxcolumn;
void adderror (sizetype error);
void seterror (sizetype error);
void geterror (sizetype& error) const;
};

COMPLEXMATRIX eye     (unsigned int rows);
COMPLEXMATRIX zeroes  (unsigned int rows,
                           unsigned int columns);
COMPLEXMATRIX ones    (unsigned int rows,
                           unsigned int columns);
// REAL       maxnorm (const REALMATRIX& x);
// REAL       rowsumnorm (const REALMATRIX& x);
// REAL       colsumnorm (const REALMATRIX& x);
// REALMATRIX solve (
//        REALMATRIX& lside,
//        REALMATRIX& rside,
//        int use_pivot);
// REALMATRIX limit_lip (REALMATRIX f(int, const REALMATRIX&),
//                            int lip,
//                            bool on_domain(const REALMATRIX&),
//                            const REALMATRIX& x);

} // namespace iRRAM

#endif