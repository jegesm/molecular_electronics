/* diag.f -- translated by f2c (version 19950602).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

#include <omp.h>
#include <iostream>
#include <vector>
#include "f2c.h"

extern double beginMP;

extern int cchol(int *n, int *nd, real *a, int *fail);
extern int ctred2(int *n, int *nd, real *a, real *b,
                  real *d,real *e,real *f);

extern int ctql2(int *n, int *nd, real *d, real *e,
                 real *f,real *a,real *b,int *fail);
extern double d_sign(double a,double b);

/* SUBROUTINE CCHOL COMPUTES CHOLESKI */
/* DECOMPOSITION OF GIVEN COMPLEX POSITIVE DEFINITE */
/* MATRIX A. */

/* INPUT DATA */

/*  N --- ORDER OF MATRIX */
/*  ND -- DIMENSION OF ARRAY A (IT CAN BE */
/*        GREATER THAN OR EQUAL TO N) */
/*  A --- GIVEN MATRIX */
/*        IT IS SUPPOSED TO BE STORED IN THE */
/*        FOLLOWING WAY.  DIAGONAL ELEMENTS, */
/*        BEING REAL, ARE STORED ON THE DIAGONAL, */
/*        REAL PARTS OF OFFDIAGONAL ELEMENTS */
/*        ARE STORED IN THE LOWER TRIANGLE OF A, */
/*        IMAG PARTS OF THE LOWER TRIANGLE ARE */
/*        STORED IN THE UPPER TRIANGLE OF A. */

/*     EXIT INFORMATION */

/*  A --- COMPLEX ELEMENTS OF MATRIX L, DEFINED BY */
/*     A=L*L(H) */
/*        ARE STORED IN THE SAME WAY AS ORIGINAL */
/*        ELEMENTS OF A, THAT IS REAL PARTS OF THE */
/*        LOWER TRIANGLE OF L IN THE LOWER TRIANGLE */
/*        OF A AND THE CORRESPONDING IMAG PARTS IN */
/*        THE UPPER TRIANGLE OF A. */
/*  FAIL --- IS SET TO ZERO IF THE DECOMPOSITION WAS */
/*           SUCCESFUL AND TO NONZERO IF */
/*           THE MATRIX WAS NOT POSITIVE DEFINITE. */

/*    PROGRAMMED BY E. ZAKRAJSEK */
/*     JUNE 20, 1974 */



//Refactoring was carried out, original contained GOTO language constructs
int cchol(int *n, int *nd, real *a, int *fail)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */

    int ia, ka;

/*     SUPPOSE DECOMPOSITION WILL FAIL */

    /* Parameter adjustments */
    a_dim1 = *nd;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    *fail = 1;

    i__1 = *n;

    /*     TEST FOR POSITIVE DEFINITNESS */
    bool notposdef = false;
    for (int i = 1; i <= i__1; ++i) {
        if (a[i + i * a_dim1] <= 0.) {
            notposdef=true;
        }
    }
    if(notposdef)
        return 0;



        /*      COMPUTE COLUMN I */
        for (int i = 1; i <= i__1; ++i) {
            a[i + i * a_dim1] = (double) std::sqrt((double) a[i + i * a_dim1]);
        }


        //i cannot be equal to i__1, previously achieved with goto
        for (int i = 1; i < i__1; ++i) {
            for (int j = i + 1; j <= i__1; ++j) { //2-4596
                a[j + i * a_dim1] /= a[i + i * a_dim1];
                a[i + j * a_dim1] /= a[i + i * a_dim1];
            }
        }

        /*     REDUCE REMAINING COLUMNS */
        for (int i = 1; i < i__1; ++i) {
            //#pragma omp for
            for (int k = ia; k <= i__1; ++k) {
                a[k + k * a_dim1] = a[k + k * a_dim1] - a[k + i * a_dim1] *
                                                        a[k + i * a_dim1] - a[i + k * a_dim1] * a[i + k * a_dim1];
            }
        }

        for (int i = 1; i < i__1; ++i) {
            int ia = i + 1;
            for (int k = ia; k < i__1; ++k) { //k cannot be equal to i__1, previously achieved with goto
                int ka = k + 1;
                for (int j = ka; j <= i__1; ++j) {
                    a[j + k * a_dim1] = a[j + k * a_dim1] - a[j + i * a_dim1] * a[
                            k + i * a_dim1] - a[i + j * a_dim1] * a[i + k *
                                                                        a_dim1];
                }
            }
        }
        for (int i = 1; i < i__1; ++i) {
            int ia = i + 1;
            for (int k = ia; k < i__1; ++k) { //k cannot be equal to i__1, previously achieved with goto
                int ka = k + 1;
                for (int j = ka; j <= i__1; ++j) {
                    a[k + j * a_dim1] = a[k + j * a_dim1] - a[i + j * a_dim1] * a[
                            k + i * a_dim1] + a[j + i * a_dim1] * a[i + k *
                                                                        a_dim1];
                }
            }
        }

    *fail = 0;
    return 0;
}

int ctred2(int *n, int *nd, real *a, real *b, real *d, real *e, real *f)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    real d__1;

    /* Builtin functions */

    /* Local variables */
    static real chep, c, g;
    static int i, j, k, l;
    static real r, s, t;
    static int ia, kk;
    static real sm, ali, all, alr;



/*     SUBROUTINE CTRED2 REDUCES GIVEN COMPLEX */
/*     HERMITIAN MATRIX TO A TRIDIAGONAL FORM */

/*     PARAMETERS */

/*     N    --- ORDER OF THE MATRIX */
/*     ND   --- DIMENSION OF ARRAYS A AND B */
/*     A    --- GIVEN MATRIX, REPLACED BY REAL PART */
/*              OF THE TRANSFORMATION MATRIX */
/*     B    --- IMAG PART OF TRANSFORMATION MATRIX */
/*     D    --- DIAGONAL PART OF THE TRIADIAGONAL MATRIX */
/*     E    --- REAL PART OF THE CODIAGONAL OF THE */
/*              TRIDIAGONAL MATRIX */
/*              (LAST N-1 LOCATIONS) */
/*     F    --- IMAG PARTS OF THE LOWER CODIAGONAL. */

/*     THE GIVEN MATRIX SHOULD BE STORED IN THE */
/*     FOLLOWING WAY */

/*          --- DIAGONAL ELEMENTS IN THE DIAGONAL */
/*          --- REAL PART OF THE LOWER TRIANGLE IN THE */
/*              LOWER TRIANGLE */
/*          --- IMAG PARTS OF THE LOWER TRIANGLE */
/*              IN THE UPPER TRIANGLE */


/*     PROGRAMMED BY E. ZAKRAJSEK */
/*     JUNE 20,1974 */



    /* Parameter adjustments */
    --f;
    --e;
    --d;
    b_dim1 = *nd;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *nd;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    chep = 1.3877787807814457e-17;
    d[1] = a[a_dim1 + 1];
    if (*n == 1) {
        goto L31;
    }

/* MAIN K LOOP */

    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
        l = k - 1;

/*     COMPUTE NORM */

        all = 0.;
        for (i = k; i <= i__1; ++i) {
/* L10: */
            all = all + a[i + l * a_dim1] * a[i + l * a_dim1] + a[l + i *
                                                                      a_dim1] * a[l + i * a_dim1];
        }
        all = (double)std::sqrt((double)all);

/*     COMPUTE CONSTANTS */

        c = 1.;
        s = 0.;
        r = (double)std::sqrt((double)(a[k + l * a_dim1] * a[k + l * a_dim1] + a[l + k * a_dim1] *
                                                                               a[l + k * a_dim1]));
        if (std::fabs(r) < 1e-50) {
            r = 0.;
        }
        if (r == 0.) {
            goto L11;
        }
        c = a[k + l * a_dim1] / r;
        s = a[l + k * a_dim1] / r;
        L11:
        alr = all * c;
        ali = all * s;
        a[l + l * a_dim1] = 0.;

/*     TEST FOR SUPERFLUOUS TRANSFORMATION */

        sm = all * (all + r);
        if (std::fabs(sm) < 1e-50) {
            sm = 0.;
        }
        if (sm == 0.) {
            goto L20;
        }
        g = 1. / sm;
        a[l + l * a_dim1] = g;

        a[k + l * a_dim1] += alr;
        a[l + k * a_dim1] += ali;

/*     NOW COMPUTE U=A*W */
/*     AND STORE INTO (E,F) */

        t = 0.;
        for (i = k; i <= i__1; ++i) {
            c = a[i + i * a_dim1] * a[i + l * a_dim1];
            s = a[i + i * a_dim1] * a[l + i * a_dim1];
            if (i == k) {
                goto L13;
            }
            ia = i - 1;
            i__3 = ia;
            for (j = k; j <= i__3; ++j) {
                c = c + a[i + j * a_dim1] * a[j + l * a_dim1] - a[j + i *
                                                                      a_dim1] * a[l + j * a_dim1];
/* L12: */
                s = s + a[i + j * a_dim1] * a[l + j * a_dim1] + a[j + i *
                                                                      a_dim1] * a[j + l * a_dim1];
            }
            L13:
            if (i == *n) {
                goto L15;
            }
            ia = i + 1;
            for (j = ia; j <= i__1; ++j) {
                c = c + a[j + i * a_dim1] * a[j + l * a_dim1] + a[i + j *
                                                                      a_dim1] * a[l + j * a_dim1];
/* L14: */
                s = s + a[j + i * a_dim1] * a[l + j * a_dim1] - a[i + j *
                                                                      a_dim1] * a[j + l * a_dim1];
            }
            L15:
            e[i] = g * c;
            f[i] = g * s;
/* L16: */
            t = t + a[i + l * a_dim1] * c + a[l + i * a_dim1] * s;
        }
        t *= g * g;

/*    TRANSFORM  MATRIX */
        for (i = k; i <= i__1; ++i) {
            a[i + i * a_dim1] = a[i + i * a_dim1] - (a[i + l * a_dim1] * e[i]
                                                     + a[l + i * a_dim1] * f[i]) * 2. + t * (a[i + l * a_dim1]
                                                                                             * a[i + l * a_dim1] + a[l + i * a_dim1] * a[l + i *
                                                                                                                                             a_dim1]);
            if (i == k) {
                goto L18;
            }
            ia = i - 1;
            i__3 = ia;
            for (j = k; j <= i__3; ++j) {
                a[i + j * a_dim1] = a[i + j * a_dim1] - a[i + l * a_dim1] * e[
                        j] - a[l + i * a_dim1] * f[j] - a[j + l * a_dim1] * e[
                        i] - a[l + j * a_dim1] * f[i] + t * (a[i + l * a_dim1]
                                                             * a[j + l * a_dim1] + a[l + i * a_dim1] * a[l + j *
                                                                                                             a_dim1]);
                a[j + i * a_dim1] = a[j + i * a_dim1] - a[l + i * a_dim1] * e[
                        j] + a[i + l * a_dim1] * f[j] + a[l + j * a_dim1] * e[
                        i] - a[j + l * a_dim1] * f[i] + t * (a[l + i * a_dim1]
                                                             * a[j + l * a_dim1] - a[i + l * a_dim1] * a[l + j *
                                                                                                             a_dim1]);
/* L17: */
            }
            L18:
            ;
        }

/*     STORE DIAGONAL AND CODIAGONAL ELEMENTS */

        L20:
        d[k] = a[k + k * a_dim1];
        e[k] = -alr;
        f[k] = -ali;
/* L30: */
    }

/*     NOW ACCUMULATE TRANSFORMATIONS */

    L31:
    a[*n + *n * a_dim1] = 1.;
    b[*n + *n * b_dim1] = 0.;
    if (*n == 1) {
        return 0;
    }
    i__1 = *n;
    for (kk = 2; kk <= i__1; ++kk) {
        k = *n - kk + 2;
        l = k - 1;

/*     SKIP TRANSFORMATION IF UNIT */

        if ((d__1 = a[l + l * a_dim1], std::fabs(d__1)) < 1e-50) {
            a[l + l * a_dim1] = 0.;
        }
        if (a[l + l * a_dim1] == 0.) {
            goto L36;
        }

/*     COMPUTE PRODUCT */

        for (j = k; j <= i__1; ++j) {
            c = 0.;
            s = 0.;
            for (i = k; i <= i__1; ++i) {
                c = c + a[i + l * a_dim1] * a[i + j * a_dim1] + a[l + i *
                                                                      a_dim1] * b[i + j * b_dim1];
/* L33: */
                s = s + a[i + l * a_dim1] * b[i + j * b_dim1] - a[l + i *
                                                                      a_dim1] * a[i + j * a_dim1];
            }
            c *= a[l + l * a_dim1];
            s *= a[l + l * a_dim1];
            for (i = k; i <= i__1; ++i) {
                a[i + j * a_dim1] = a[i + j * a_dim1] - c * a[i + l * a_dim1]
                                    + s * a[l + i * a_dim1];
/* L34: */
                b[i + j * b_dim1] = b[i + j * b_dim1] - c * a[l + i * a_dim1]
                                    - s * a[i + l * a_dim1];
            }
/* L35: */
        }

/*     MAKE NEW LINE */

        L36:
        for (i = k; i <= i__1; ++i) {
            a[i + l * a_dim1] = 0.;
            a[l + i * a_dim1] = 0.;
            b[i + l * b_dim1] = 0.;
/* L37: */
            b[l + i * b_dim1] = 0.;
        }
        a[l + l * a_dim1] = 1.;
        b[l + l * b_dim1] = 0.;
/* L40: */
    }
    return 0;
} /* ctred2 */

/* Subroutine */ int ctql2(int *n, int *nd, real *d, real *e, real *f, real *a, real *b, int *fail)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    real d__1, d__2;

    /* Local variables */
    real chep, c, g, h;
    int i, j, k, l, m;
    real p, r, s;
    int i1;
    real bb;
    int ia;
    real ff;
    int ma;
    real hi, hr;



/*     SUBROUTINE CTQL2 COMPUTES THE EIGENVALUES AND */
/*     EIGENVECTORS OF A COMPLEX HERMITIAN TRIDIAGONAL */
/*     MATRIX */

/*     PARAMETERS */

/*     N    --- ORDER OF MATRIX */
/*     ND   --- DIMENSION OF A AND B */
/*     D    --- DIAGONAL GIVEN */
/*     E    --- REAL PART OF CODIAGONAL GIVEN */
/*              (LAST N-1 LOCATIONS) */
/*     F    --- IMAG PART OF THE LOWER CODIAGONAL */
/*     A    --- REAL PART OF EIGENVECTORS */
/*     B    --- IMAG PART OF EIGENVECTORS */
/*     FAIL --- RECEIVES VALUE OF 1 INSTEAD OF ZERO */
/*              IF SOME EIGENVALUE TAKES MORE THAN 30 */
/*              ITERATIONS. */


/*     EIGENVALUES ARE OBTAINED IN INCREASING OF */
/*     MAGNITUDE IN VECTOR D, EIGENVECTORS ARE STORED */
/*     BY COLUMNS.  ARRAYS A AND B SHOULD BE PRESET TO */
/*     SOME UNITARY MATRIX SUCH AS THE IDENTITY MATRIX */
/*     OR THE MATRIX PRODUCED BY CTRED2. */


/*     PROGRAMMED BY E.  ZAKRAJSEK */
/*     JUNE 21, 1974 */




/*     *************************************** */
/*     *                                     * */
/*     * NEXT LINE OF PROGRAM DEFINES        * */
/*     * MACHINE DEPENDENT CONSTANT CHEP     * */
/*     * DEFINED AS THE SMALLEST REAL        * */
/*     * NUMBER FOR WHICH                    * */
/*     *                                     * */
/*     *        1.0+CHEP .GT. 1.0            * */
/*     *                                     * */
/*     *************************************** */

    /* Parameter adjustments */
    b_dim1 = *nd;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *nd;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --f;
    --e;
    --d;

    /* Function Body */
    chep = 1.3877787807814457e-17;

/*     FIRST MAKE REAL CODIAGONAL MOVED DOWN */
/*     TO FIRST LOCATION */

    if (*n == 1) {
        goto L12;
    }
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
        r = (double)std::sqrt((double)(e[k] * e[k] + f[k] * f[k]));
        if (std::fabs(r) < 1e-50) {
            r = 0.;
        }
        if (r == 0.) {
            goto L11;
        }
        c = e[k] / r;
        s = f[k] / r;

/*     ACCUMULATE ROTATION */
        for (i = 1; i <= i__1; ++i) {
            p = a[i + k * a_dim1] * c - b[i + k * b_dim1] * s;
            b[i + k * b_dim1] = a[i + k * a_dim1] * s + b[i + k * b_dim1] * c;
/* L10: */
            a[i + k * a_dim1] = p;
        }

/*     TRANSFORM NEXT E */

        if (k == *n) {
            goto L11;
        }
        l = k + 1;
        p = e[l] * c - f[l] * s;
        f[l] = e[l] * s + f[l] * c;
        e[l] = p;
        L11:
        e[k - 1] = r;
    }
    L12:
    e[*n] = 0.;

/*     INITIALIZE */

    bb = 0.;
    ff = 0.;
    *fail = 1;

/*     MAIN LOOP */

    i__1 = *n;
    for (l = 1; l <= i__1; ++l) {
        j = 0;
        h = chep * ((d__1 = d[l], std::fabs(d__1)) + (d__2 = e[l], std::fabs(d__2)));
        if (bb < h) {
            bb = h;
        }

/*     LOOK FOR SMALL SUBDIAGONAL ELEMENT */
        for (m = l; m <= i__1; ++m) {
            if ((d__1 = e[m], std::fabs(d__1)) <= bb) {
                goto L21;
            }
/* L20: */
        }
        L21:
        if (m == l) {
            goto L31;
        }

/*     NEXT ITERATION */

        L24:
        if (j == 30) {
            return 0;
        }
        ++j;

/*     FORM SHIFT */

        p = (d[l + 1] - d[l]) / (e[l] * 2.);
        r = std::sqrt((double)(p * p + 1.));
        h = d[l] - e[l] / (p + d_sign(r, p));

        for (i = l; i <= i__1; ++i) {
/* L25: */
            d[i] -= h;
        }
        ff += h;

/*     QL TRANSFORMATION */

        p = d[m];
        c = 1.;
        s = 0.;
        ma = m - 1;

        i__2 = ma;
        for (ia = l; ia <= i__2; ++ia) {
            i = ma - ia + l;
            i1 = i + 1;
            g = c * e[i];
            h = c * p;
            if (std::fabs(p) < (d__1 = e[i], std::fabs(d__1))) {
                goto L26;
            }

            c = e[i] / p;
            r = (double)std::sqrt((double)(c * c + 1.));
            e[i1] = s * p * r;
            s = c / r;
            c = 1. / r;
            goto L27;

            L26:
            c = p / e[i];
            r = (double)std::sqrt((double)(c * c + 1.));
            e[i1] = s * e[i] * r;
            s = 1. / r;
            c /= r;

            L27:
            p = c * d[i] - s * g;
            d[i1] = h + s * (c * g + s * d[i]);

/*     FORM VECTOR */
            for (k = 1; k <= i__1; ++k) {
                hr = a[k + i1 * a_dim1];
                hi = b[k + i1 * b_dim1];
                a[k + i1 * a_dim1] = s * a[k + i * a_dim1] + c * hr;
                b[k + i1 * b_dim1] = s * b[k + i * b_dim1] + c * hi;
                a[k + i * a_dim1] = c * a[k + i * a_dim1] - s * hr;
/* L28: */
                b[k + i * b_dim1] = c * b[k + i * b_dim1] - s * hi;
            }

/* L30: */
        }

        e[l] = s * p;
        d[l] = c * p;
        if ((d__1 = e[l], std::fabs(d__1)) > bb) {
            goto L24;
        }

/*     ROOT FOUND */

        L31:
        d[l] += ff;
/* L32: */
    }

/*     ORDER EIGENVALUES AND EIGENVECTORS */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
        k = i;
        for (j = i; j <= i__1; ++j) {
            if (d[j] < d[k]) {
                k = j;
            }
/* L40: */
        }


        if (k == i) {
            goto L42;
        }
        p = d[i];
        d[i] = d[k];
        d[k] = p;

        for (j = 1; j <= i__1; ++j) {
            p = a[j + i * a_dim1];
            a[j + i * a_dim1] = a[j + k * a_dim1];
            a[j + k * a_dim1] = p;

            p = b[j + i * b_dim1];
            b[j + i * b_dim1] = b[j + k * b_dim1];
/* L41: */
            b[j + k * b_dim1] = p;
        }
        L42:
        ;
    }
    *fail = 0;
    return 0;
} /* ctql2 */


