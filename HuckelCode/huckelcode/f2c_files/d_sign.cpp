#include "f2c.h"
/****************************************************************************
*
*                   Procedure d_sign
*
* Arguments:  a,b: doubles
*
* Returns: double
*
* Action:  This mimics the fortran 77 function sign.
*
*
*        b\a|  <0  |   >=0  |
*        ---|------|--------|
*         <0|  a   |  -a    |
*        >=0| -a   |   a    |
*
*
*****************************************************************************/
double d_sign(double a,double b)
{
    double x;

    if( a >= 0 ) x = a;
    else x = -a;

    if( b >= 0 ) return x;
    else return -x;
}