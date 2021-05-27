// SphericalHarmonics.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>
#include "SphericalHarmonics.h"
#include "Legendre.h"
#include "SphericalCalc.h"











//Encapsulation of the Scalar Harmonic Function Y.
//How to use:
//1) declare and initialize a variable Y y(m , n , kappa);
//2) call the function as y(x); where x is a spherical coordinate point.
class Y : public SurfaceScalarFunction
{
private:
    int m;
    int n;
    //kappa value for computing gradient. See SurfaceGrad for more.
    double kappa;
    Legendre poly;

public:
    Y(int a, int b , double k)
    {
        m = a;
        n = b;
        kappa = k;
        

    }


    //Implements the function Y 
    std::complex<double> operator()(SurfaceCoord s)
    {
        
        poly.populate(cos(s.theta.real()) , m, n);
        const std::complex<double> I(0.0, 1.0);
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * MATHPI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;
        
        coef *= poly.getValue(m, n);

        //std::cout << coef << std::endl;
        
        return coef * exp( I * (double)m * s.phi);


    }
};

//Encapsulation of gradient of Scalar Harmonic function. Not used.
class G : public SphericalVectorField
{
private:
    int m;
    int n;
    double h;
    double k;
public:
    G(int a, int b , double steps , double kap)
    {
        m = a;
        n = b;
        h = steps;
        k = kap;
    }

    PolarCoord operator()(PolarCoord x)
    {
        Y y(m, n , k);
        return surfaceGrad(&y, x.s, h , k);

    }


};


//Spherical Harmonic Basis function encapsulations.
class V : public SphericalVectorField
{
    int m;
    int n;
    double h;
    double k;
public:

    V(){}
    V(int a, int b, double steps , double kap)
    {
        m = a;
        n = b;
        h = steps;
        k = kap;
    }

    PolarCoord operator()(PolarCoord x)
    {
        PolarCoord temp;
        Y y(m, n , k);

        temp = -((double)n+1) * y(x.s) * PolarCoord(1, x.s);
        temp = temp + surfaceGrad(&y, x.s, h , k);

        return temp;

    }

    void reset(int a, int b)
    {
        m = a;
        n = b;
    }

};
class W : public SphericalVectorField
{
    int m;
    int n;
    double h;
    double k;
public:
    W(){}

    W(int a, int b, double steps , double kap)
    {
        m = a;
        n = b;
        h = steps;
        k = kap;
    }

    PolarCoord operator()(PolarCoord x)
    {
        PolarCoord temp;
        Y y(m, n , k);

        temp = (double)n * y(x.s) * PolarCoord(1 , x.s);
        temp = surfaceGrad(&y, x.s, h, k) + temp;

        return temp;

    }
    void reset(int a, int b)
    {
        m = a;
        n = b;
    }
};
class X : public SphericalVectorField
{
    int m;
    int n;
    double h;
    double k;

public:

    X(){}

    X(int a, int b, double steps, double kap)
    {
        m = a;
        n = b;
        h = steps;
        k = kap;
    }

    PolarCoord operator()(PolarCoord x)
    {
        PolarCoord temp;
        PolarCoord er = PolarCoord(1.0, x.s);
        Y y(m, n , k);

        
        temp = surfaceGrad(&y, x.s, h , k);

        return RectToSphere(cross(SphereToRect(x), SphereToRect(temp)));

    }

    void reset(int a, int b)
    {
        m = a;
        n = b;
    }

};

//Encapsulation of a single term in the expansion of a function.
class SphericalHarmonic : public SphericalVectorField
{
private:

    std::complex<double> coef;
    std::complex<double> denom;
    SphericalVectorField* f;

public:
    SphericalHarmonic() { coef = 1; f = nullptr; }

    SphericalHarmonic(std::complex<double> c, SphericalVectorField& f0)
    {
        coef = c;
        f = &f0;
        denom = L2InnerProduct(f, f, 100);
    }

    PolarCoord operator()(PolarCoord x)
    {
        PolarCoord fx = (*f)(x);
        PolarCoord Fnormalized = (norm(denom) < 10e-9) ? PolarCoord(0.0 , SurfaceCoord(0.0 , 0.0)) : (1.0 / norm(denom) )* fx;
        return coef * Fnormalized;
    }

};

//The function to approximate.
class rho : public SphericalVectorField
{

private:

    double step;
    double k;
public:
    rho(){}
    rho(double steps, double kappa) 
    {
        step = steps;
        k = kappa;
    }

    PolarCoord operator()(PolarCoord x)
    {
        //identity function
        //return x;

        //scale function.
        //return 2 * x;


        //rotation of x by pi/4.
        //return PolarCoord(x.rho , x.s + SurfaceCoord(MATHPI / 4.0 , 0));

        V z(2, 2, step, k);

        return z(x);
    }
};

int main()
{
    std::setprecision(12);
    
    
    //numerical parameters.
    
    //series truncation
    const int N = 5;
    
    //step size for gradient.
    const double GRADSTEP = 10e-6;

    
    //kappa for gradient computation.
    const double KAPPA = 10e-6;

    //number of panel points per dimension in integrals.
    const int NUMGRIDS = 100;

    //Basis functions. one for each term in the series.
    static V v[2*N + 1][N];
    static W w[2*N + 1][N];
    static X x[2*N + 1][N];

    rho r(GRADSTEP, KAPPA);
   
    //set up the basis functions. an index of [m+N][n] corresponds to negative values of m.
    for(int m = 0; m < N; m++)
        for (int n = 0; n < N; n++)
        {
            v[m][n] = V(m, n, GRADSTEP , KAPPA);
            v[m + N][n] = V(-m, n, GRADSTEP , KAPPA);
            w[m][n] = W(m, n, GRADSTEP, KAPPA);
            w[m + N][n] = W(-m, n, GRADSTEP, KAPPA);
            x[m][n] = X(m, n, GRADSTEP, KAPPA);
            x[m + N][n] = X(-m, n, GRADSTEP , KAPPA);
        }

    //calculate the inner products;

    //for V
    static std::complex<double> rhohatV[2 * N + 1][N];
    std::cout << "Computing inner product coefficiencts for V... " << std::endl;
    for (int n = 0; n < N; n++)
        for (int m = 0; m <= n; m++)
        {

        
             rhohatV[m][n] = L2InnerProduct(&r, &v[m][n], NUMGRIDS);

             rhohatV[m+N][n] = L2InnerProduct(&r, &v[m + N][n], NUMGRIDS);

             std::cout << rhohatV[m][n] <<" "<<  rhohatV[m + N][n] << std::endl;
       
        }
    std::cout << "Done!" << std::endl;

    //for W.
    static std::complex<double> rhohatW[2 * N + 1][N];
    std::cout << "Computing inner product coefficiencts for W... " << std::endl;
    for (int n = 0; n < N; n++)
        for (int m = 0; m <= n; m++)
        {

            rhohatW[m][n] = L2InnerProduct(&r, &w[m][n], NUMGRIDS);

            rhohatW[m + N][n] = L2InnerProduct(&r, &w[m + N][n], NUMGRIDS);

            std::cout << rhohatW[m][n] << " " << rhohatW[m + N][n] << std::endl;

        }
    std::cout << "Done!" << std::endl;

    //and for X.
    static std::complex<double> rhohatX[2 * N + 1][N];
    std::cout << "Computing inner product coefficiencts for X... " << std::endl;
    for (int n = 0; n < N; n++)
        for (int m = 0; m <= n; m++)
        {

            rhohatX[m][n] = L2InnerProduct(&r, &x[m][n], NUMGRIDS);

            rhohatX[m + N][n] = L2InnerProduct(&r, &x[m + N][n], NUMGRIDS);

            std::cout << rhohatX[m][n] << " " << rhohatX[m + N][n] << std::endl;

        }
    std::cout << "Done!" << std::endl <<std::endl;

    std::cout << "Gathering summation... ";
  

    //declare terms in the series.
    SphericalHarmonic Vterm[2 * N + 1][N];
    SphericalHarmonic Wterm[2 * N + 1][N];
    SphericalHarmonic Xterm[2 * N + 1][N];


    //declare the approximation.
    VectorFieldSum rhoapprox;

    //construct the summation. Note we must compute the L2 norm of the basis functions to do so, so this takes a few seconds..
        for (int n = 0; n < N; n++)
            for (int m = 0; m <= n; m++)
            {
                
                Vterm[m][n] = SphericalHarmonic(rhohatV[m][n], v[m][n]);
                if(m > 0)
                Vterm[m+N][n] = SphericalHarmonic(rhohatV[m+N][n], v[m+N][n]);
                

                Wterm[m][n] = SphericalHarmonic(rhohatW[m][n], w[m][n]);
                if(m > 0)
                Wterm[m + N][n] = SphericalHarmonic(rhohatW[m+N][n], w[m+N][n]);

                

                Xterm[m][n] = SphericalHarmonic(rhohatX[m][n], x[m][n]);
                if (m > 0)
                Xterm[m + N][n] = SphericalHarmonic(rhohatX[m + N][n], x[m + N][n]);

                rhoapprox.append(Vterm[m][n]);
                if (m > 0)
                    rhoapprox.append(Vterm[m + N][n]);

                rhoapprox.append(Wterm[m][n]);
                if (m > 0)
                    rhoapprox.append(Wterm[m + N][n]);
                
                rhoapprox.append(Xterm[m][n]);
                if (m > 0)
                    rhoapprox.append(Xterm[m + N][n]);

                
                std::cout << "|";

            }
    std::cout << "Done!" << std::endl << std::endl;


    //Compute the error in approximating rho by it's expansion in Spherical Harmonics.
    std::cout << "Error in approximation(L2): ";
    std::cout << L2Difference(&r, &rhoapprox, NUMGRIDS) << std::endl;

    std::cout << "Error in approximation(LInf): ";
    std::cout << LInfDifference(&r, &rhoapprox, NUMGRIDS) << std::endl;

    return 0;
}


