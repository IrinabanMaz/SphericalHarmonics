#pragma once

#include<vector>


void legendre_compute_glr(int n, double x[], double w[]);
void legendre_compute_glr0(int n, double* p, double* pp);
void legendre_compute_glr1(int n, double* roots, double* ders);
void legendre_compute_glr2(double p, int n, double* roots, double* ders);


double rk2_leg(double t, double tn, double x, int n);
double ts_mult(double* u, double h, int n);



template<class T>
std::vector<T> vectorFromPtr(T* objs, size_t size)
{
    std::vector<T> temp;
    temp.resize(size);

    for (size_t i = 0; i < size; i++)
        temp[i] = objs[i];

    return temp;
}


class MathConstants
{
public:
    double PI = 4.0 * atan(1.0);
    double NUMGLNODES;
    double NUMTRAPNODES;
    std::vector<double> GLnodes;
    std::vector<double> GLweights;


    MathConstants(const int& ntrap = 8, const int& ngl = 8)
    {
        NUMTRAPNODES = ntrap;
        NUMGLNODES = ngl;
        double* glweights = new double[ngl];
        double* glnodes = new double[ngl];

        legendre_compute_glr(ngl, glnodes, glweights);

        GLweights = vectorFromPtr(glweights, ngl);
        GLnodes = vectorFromPtr(glnodes, ngl);

        delete[] glweights;
        delete[] glnodes;


    }

    MathConstants& operator =(const MathConstants& consts)
    {
        NUMGLNODES = consts.NUMGLNODES;
        NUMTRAPNODES = consts.NUMTRAPNODES;
        GLweights = consts.GLweights;
        GLnodes = consts.GLnodes;

        return *this;
    }

    MathConstants(const MathConstants& consts)
    {
        NUMGLNODES = consts.NUMGLNODES;
        NUMTRAPNODES = consts.NUMTRAPNODES;
        GLweights = consts.GLweights;
        GLnodes = consts.GLnodes;
    }
};





/****************************************************************************
*
*  Purpose:
*
*    LEGENDRE_COMPUTE_GLR: Legendre quadrature by the Glaser-Liu-Rokhlin method.
*
*  Licensing:
*
*    This code is distributed under the GNU LGPL license.
*
*  Modified:
*
*    20 October 2009
*
*  Author:
*
*    Original C++ version by Nick Hale.
*    This C++ version by John Burkardt.
*
*  Reference:
*
*    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin,
*    A fast algorithm for the calculation of the roots of special functions,
*    SIAM Journal on Scientific Computing,
*    Volume 29, Number 4, pages 1420-1438, 2007.
*
*  Parameters:
*
*    Input, int N, the order.
*
*    Output, double X[N], the abscissas.
*
*    Output, double W[N], the weights.
*/
void legendre_compute_glr(int n, double x[], double w[])


{
    int i;
    double p;
    double pp;
    double w_sum;
    //
    //  Get the value and derivative of the N-th Legendre polynomial at 0.
    //
    legendre_compute_glr0(n, &p, &pp);
    //
    //  If N is odd, then zero is a root.
    //  
    if (n % 2 == 1)
    {
        x[(n - 1) / 2] = p;
        w[(n - 1) / 2] = pp;
    }
    //
    //  If N is even, we have to call a function to find the first root.
    //
    else
    {
        legendre_compute_glr2(p, n, &x[n / 2], &w[n / 2]);
    }
    //
    //  Get the complete set of roots and derivatives.
    //
    legendre_compute_glr1(n, x, w);
    //
    //  Compute the W.
    //
    for (i = 0; i < n; i++)
    {
        w[i] = 2.0 / (1.0 - x[i]) / (1.0 + x[i]) / w[i] / w[i];
    }
    w_sum = 0.0;
    for (i = 0; i < n; i++)
    {
        w_sum = w_sum + w[i];
    }
    for (i = 0; i < n; i++)
    {
        w[i] = 2.0 * w[i] / w_sum;
    }
    return;
}
//****************************************************************************
/****************************************************************************
*
*  Purpose:
*
*    LEGENDRE_COMPUTE_GLR0 gets a starting value for the fast algorithm.
*
*  Licensing:
*
*    This code is distributed under the GNU LGPL license.
*
*  Modified:
*
*    19 October 2009
*
*  Author:
*
*    Original C++ version by Nick Hale.
*    This C++ version by John Burkardt.
*
*  Reference:
*
*    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin,
*    A fast algorithm for the calculation of the roots of special functions,
*    SIAM Journal on Scientific Computing,
*    Volume 29, Number 4, pages 1420-1438, 2007.
*
*  Parameters:
*
*    Input, int N, the order of the Legendre polynomial.
*
*    Output, double *P, *PP, the value of the N-th Legendre polynomial
*    and its derivative at 0.
*/
void legendre_compute_glr0(int n, double* p, double* pp)


{
    double dk;
    int k;
    double pm1;
    double pm2;
    double ppm1;
    double ppm2;

    pm2 = 0.0;
    pm1 = 1.0;
    ppm2 = 0.0;
    ppm1 = 0.0;

    for (k = 0; k < n; k++)
    {
        dk = (double)k;
        *p = -dk * pm2 / (dk + 1.0);
        *pp = ((2.0 * dk + 1.0) * pm1 - dk * ppm2) / (dk + 1.0);
        pm2 = pm1;
        pm1 = *p;
        ppm2 = ppm1;
        ppm1 = *pp;
    }
    return;
}
//****************************************************************************


/****************************************************************************
*
*  Purpose:
*
*    LEGENDRE_COMPUTE_GLR1 gets the complete set of Legendre points and weights.
*
*  Discussion:
*
*    This routine requires that a starting estimate be provided for one
*    root and its derivative.  This information will be stored in entry
*    (N+1)/2 if N is odd, or N/2 if N is even, of X and W.
*
*  Licensing:
*
*    This code is distributed under the GNU LGPL license.
*
*  Modified:
*
*    19 October 2009
*
*  Author:
*
*    Original C++ version by Nick Hale.
*    This C++ version by John Burkardt.
*
*  Reference:
*
*    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin,
*    A fast algorithm for the calculation of the roots of special functions,
*    SIAM Journal on Scientific Computing,
*    Volume 29, Number 4, pages 1420-1438, 2007.
*
*  Parameters:
*
*    Input, int N, the order of the Legendre polynomial.
*
*    Input/output, double X[N].  On input, a starting value
*    has been set in one entry.  On output, the roots of the Legendre
*    polynomial.
*
*    Input/output, double W[N].  On input, a starting value
*    has been set in one entry.  On output, the derivatives of the Legendre
*    polynomial at the zeros.
*
*  Local Parameters:
*
*    Local, int M, the number of terms in the Taylor expansion.
*/
void legendre_compute_glr1(int n, double* x, double* w)


{
    double dk;
    double dn;
    double h;
    int j;
    int k;
    int l;
    int m = 30;
    int n2;
    static double pi = 3.141592653589793;
    int s;
    double* u;
    double* up;
    double xp;

    if (n % 2 == 1)
    {
        n2 = (n - 1) / 2 - 1;
        s = 1;
    }
    else
    {
        n2 = n / 2 - 1;
        s = 0;
    }

    u = new double[m + 2];
    up = new double[m + 1];

    dn = (double)n;

    for (j = n2 + 1; j < n - 1; j++)
    {
        xp = x[j];

        h = rk2_leg(pi / 2.0, -pi / 2.0, xp, n) - xp;

        u[0] = 0.0;
        u[1] = 0.0;
        u[2] = w[j];

        up[0] = 0.0;
        up[1] = u[2];

        for (k = 0; k <= m - 2; k++)
        {
            dk = (double)k;

            u[k + 3] =
                (
                    2.0 * xp * (dk + 1.0) * u[k + 2]
                    + (dk * (dk + 1.0) - dn * (dn + 1.0)) * u[k + 1] / (dk + 1.0)
                    ) / (1.0 - xp) / (1.0 + xp) / (dk + 2.0);

            up[k + 2] = (dk + 2.0) * u[k + 3];
        }

        for (l = 0; l < 5; l++)
        {
            h = h - ts_mult(u, h, m) / ts_mult(up, h, m - 1);
        }

        x[j + 1] = xp + h;
        w[j + 1] = ts_mult(up, h, m - 1);
    }

    for (k = 0; k <= n2 + s; k++)
    {
        x[k] = -x[n - 1 - k];
        w[k] = w[n - 1 - k];
    }

    delete[] u;
    delete[] up;
    return;
}
//****************************************************************************80

void legendre_compute_glr2(double pn0, int n, double* x1, double* d1)

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_COMPUTE_GLR2 finds the first real root.
//
//  Discussion:
//
//    This function is only called if N is even.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2009
//
//  Author:
//
//    Original C++ version by Nick Hale.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
//    A fast algorithm for the calculation of the roots of special functions, 
//    SIAM Journal on Scientific Computing,
//    Volume 29, Number 4, pages 1420-1438, 2007.
//
//  Parameters:
//
//    Input, double PN0, the value of the N-th Legendre polynomial
//    at 0.
//
//    Input, int N, the order of the Legendre polynomial.
//
//    Output, double *X1, the first real root.
//
//    Output, double *D1, the derivative at X1.
//
//  Local Parameters:
//
//    Local, int M, the number of terms in the Taylor expansion.
//
{
    double dk;
    double dn;
    int k;
    int l;
    int m = 30;
    static double pi = 3.141592653589793;
    double t;
    double* u;
    double* up;

    t = 0.0;
    *x1 = rk2_leg(t, -pi / 2.0, 0.0, n);

    u = new double[m + 2];
    up = new double[m + 1];

    dn = (double)n;
    //
    //  U[0] and UP[0] are never used.
    //  U[M+1] is set, but not used, and UP[M] is set and not used.
    //  What gives?
    //
    u[0] = 0.0;
    u[1] = pn0;

    up[0] = 0.0;

    for (k = 0; k <= m - 2; k = k + 2)
    {
        dk = (double)k;

        u[k + 2] = 0.0;
        u[k + 3] = (dk * (dk + 1.0) - dn * (dn + 1.0)) * u[k + 1]
            / (dk + 1.0) / (dk + 2.0);

        up[k + 1] = 0.0;
        up[k + 2] = (dk + 2.0) * u[k + 3];
    }

    for (l = 0; l < 5; l++)
    {
        *x1 = *x1 - ts_mult(u, *x1, m) / ts_mult(up, *x1, m - 1);
    }
    *d1 = ts_mult(up, *x1, m - 1);


    delete[] u;
    delete[] up;
    return;
}
/****************************************************************************
*
*  Purpose:
*
*    RK2_LEG advances the value of X(T) using a Runge-Kutta method.
*
*  Licensing:
*
*    This code is distributed under the GNU LGPL license.
*
*  Modified:
*
*    22 October 2009
*
*  Author:
*
*    Original C++ version by Nick Hale.
*    This C++ version by John Burkardt.
*
*  Parameters:
*
*    Input, double T1, T2, the range of the integration interval.
*
*    Input, double X, the value of X at T1.
*
*    Input, int N, the number of steps to take.
*
*    Output, double RK2_LEG, the value of X at T2.
*/
double rk2_leg(double t1, double t2, double x, int n)


{
    double f;
    double h;
    int j;
    double k1;
    double k2;
    int m = 10;
    double snn1;
    double t;

    h = (t2 - t1) / (double)m;
    snn1 = sqrt((double)(n * (n + 1)));
    t = t1;

    for (j = 0; j < m; j++)
    {
        f = (1.0 - x) * (1.0 + x);
        k1 = -h * f / (snn1 * sqrt(f) - 0.5 * x * sin(2.0 * t));
        x = x + k1;

        t = t + h;

        f = (1.0 - x) * (1.0 + x);
        k2 = -h * f / (snn1 * sqrt(f) - 0.5 * x * sin(2.0 * t));
        x = x + 0.5 * (k2 - k1);
    }
    return x;
}
/****************************************************************************
*
*  Purpose:
*
*    TS_MULT evaluates a polynomial.
*
*  Licensing:
*
*    This code is distributed under the GNU LGPL license.
*
*  Modified:
*
*    17 May 2013
*
*  Author:
*
*    Original C++ version by Nick Hale.
*    This C++ version by John Burkardt.
*
*  Parameters:
*
*    Input, double U[N+1], the polynomial coefficients.
*    U[0] is ignored.
*
*    Input, double H, the polynomial argument.
*
*    Input, int N, the number of terms to compute.
*
*    Output, double TS_MULT, the value of the polynomial.
*/
double ts_mult(double* u, double h, int n)


{
    double hk;
    int k;
    double ts;

    ts = 0.0;
    hk = 1.0;
    for (k = 1; k <= n; k++)
    {
        ts = ts + u[k] * hk;
        hk = hk * h;
    }
    return ts;
}
