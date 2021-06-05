// SphericalHarmonics.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>
#include "SphericalHarmonics.h"
#include "Legendre.h"
#include "SphericalCalc.h"



typedef SphereCoord(*svf)(SphereCoord);







//Encapsulations of the Scalar Harmonic Function Y.
//How to use:
//1) declare and initialize a variable Y y(m , n , kappa);
//2) call the function as y(x); where x is a spherical coordinate point.
class YReal : public SurfaceScalarFunction
{
private:
    int m;
    int n;
    Legendre poly;

public:
    YReal(int a, int b )
    {
        m = a;
        n = b;
       
    }


    //Implements the function Y 
    double operator()(SurfaceCoord s)
    {
    
            
        poly.populate(cos(s.theta) , m, n);
       
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;
        
        coef *= poly.getValue(m, n);

        //std::cout << coef << std::endl;
        
        return coef * cos((double)m * s.phi);


    }
};

class YImag : public SurfaceScalarFunction
{
private:
    int m;
    int n;
    Legendre poly;

public:
    YImag(int a, int b)
    {
        m = a;
        n = b;

    }


    //Implements the function Y 
    double operator()(SurfaceCoord s)
    {

        poly.populate(cos(s.theta), m, n);
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        coef *= poly.getValue(m, n);

        //std::cout << coef << std::endl;

        return coef * sin((double)m * s.phi);


    }
};



//Spherical Harmonic Basis function encapsulations.
class VReal : public SphericalVectorField
{
    int m;
    int n;
    double h;
public:

    VReal(){}
    VReal(int a, int b, double steps , double kap)
    {
        m = a;
        n = b;
        h = steps;
    }

    SphereCoord operator()(SphereCoord x)
    {
        SphereCoord temp;
        YReal y(m, n );

        temp = -((double)n+1) * y(x.s) * SphereCoord(1, x.s);
        temp = temp + surfaceGrad(&y, x.s, h);

        return temp;

    }

    void reset(int a, int b)
    {
        m = a;
        n = b;
    }

};

class VImag : public SphericalVectorField
{
    int m;
    int n;
    double h;
public:

    VImag() {}
    VImag(int a, int b, double steps, double kap)
    {
        m = a;
        n = b;
        h = steps;
    }

    SphereCoord operator()(SphereCoord x)
    {
        SphereCoord temp;
        YImag y(m, n);

        temp = -((double)n + 1) * y(x.s) * SphereCoord(1, x.s);
        temp = temp + surfaceGrad(&y, x.s, h);

        return temp;

    }

    void reset(int a, int b)
    {
        m = a;
        n = b;
    }

};



class WReal : public SphericalVectorField
{
    int m;
    int n;
    double h;
public:
    WReal(){}

    WReal(int a, int b, double steps , double kap)
    {
        m = a;
        n = b;
        h = steps;
    }

    SphereCoord operator()(SphereCoord x)
    {
        SphereCoord temp;
        YReal y(m, n);

        temp = (double)n * y(x.s) * SphereCoord(1 , x.s);
        temp = surfaceGrad(&y, x.s, h) + temp;

        return temp;

    }
    void reset(int a, int b)
    {
        m = a;
        n = b;
    }
};

class WImag : public SphericalVectorField
{
    int m;
    int n;
    double h;
public:
    WImag() {}

    WImag(int a, int b, double steps, double kap)
    {
        m = a;
        n = b;
        h = steps;
    }

    SphereCoord operator()(SphereCoord x)
    {
        SphereCoord temp;
        YImag y(m, n);

        temp = (double)n * y(x.s) * SphereCoord(1, x.s);
        temp = surfaceGrad(&y, x.s, h) + temp;

        return temp;

    }
    void reset(int a, int b)
    {
        m = a;
        n = b;
    }
};




class XReal : public SphericalVectorField
{
    int m;
    int n;
    double h;

public:

    XReal(){}

    XReal(int a, int b, double steps, double kap)
    {
        m = a;
        n = b;
        h = steps;
    }

    SphereCoord operator()(SphereCoord x)
    {
        SphereCoord temp;
        SphereCoord er = SphereCoord(1.0, x.s);
        YReal y(m, n);

        
        temp = surfaceGrad(&y, x.s, h);

        return cross(er , temp);

    }

    void reset(int a, int b)
    {
        m = a;
        n = b;
    }

};

class XImag : public SphericalVectorField
{
    int m;
    int n;
    double h;

public:

    XImag() {}

    XImag(int a, int b, double steps, double kap)
    {
        m = a;
        n = b;
        h = steps;
    }

    SphereCoord operator()(SphereCoord x)
    {
        SphereCoord temp;
        SphereCoord er = SphereCoord(1.0, x.s);
        YImag y(m, n);


        temp = surfaceGrad(&y, x.s, h);

        return cross(er , temp);

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

    double coefr;
    double coefi;
    double denom;
    SphericalVectorField* fr;
    SphericalVectorField* fi;


public:
    SphericalHarmonic() { coefr = 1; fr = nullptr; coefi = 1; fi = nullptr;
    }

    SphericalHarmonic(double cr, double ci ,  SphericalVectorField& f0r, SphericalVectorField& f0i)
    {
        coefr = cr;
        coefi = ci;
        fr = &f0r;
        fi = &f0i;

        denom = L2InnerProduct(fr, fr, 20) + L2InnerProduct(fi , fi , 20);
    }

    SphereCoord operator()(SphereCoord x)
    {
        SphereCoord frx = (*fr)(x);
        SphereCoord Frnormalized; 
        if (denom < 1e-18)
            Frnormalized = SphereCoord(0.0, SurfaceCoord(0.0, 0.0));
        else
            Frnormalized = (1.0 /denom)* frx;

        SphereCoord fix = (*fi)(x);
        SphereCoord Finormalized;
        if (denom < 1e-18)
            Finormalized = SphereCoord(0.0, SurfaceCoord(0.0, 0.0));
        else
            Finormalized = (1.0 / denom) * fix;

        return coefr * Frnormalized + coefi * Finormalized;
    }

};

//The function to approximate.


int main()
{

    std::cout << std::fixed << std::setprecision(5);

    
   
    
    //numerical parameters.
    
    //series truncation
    const int N = 5;
    
    //step size for gradient.
    const double GRADSTEP = 1e-5;

    //number of panel points per dimension in integrals.
    const int NUMGRIDS = 20;



    //Basis functions. one for each term in the series.
    static VReal vr[2*N + 1][N];
    static WReal wr[2*N + 1][N];
    static XReal xr[2*N + 1][N];

    static VImag vi[2 * N + 1][N];
    static WImag wi[2 * N + 1][N];
    static XImag xi[2 * N + 1][N];
    
    auto  r = SphericalVectorField([](SphereCoord a) {   return VReal(2,2,1e-5,0)(a); });


        //set up the basis functions. an index of [m+N][n] corresponds to negative values of m.
        for (int n = 0; n < N; n++)
            for (int m = 0; m <= n; m++)
            {
                vr[m][n] = VReal(m, n, GRADSTEP, 0.0);
                vr[m + N][n] = VReal(-m, n, GRADSTEP, 0.0);
                wr[m][n] = WReal(m, n, GRADSTEP, 0.0);
                wr[m + N][n] = WReal(-m, n, GRADSTEP, 0.0);
                xr[m][n] = XReal(m, n, GRADSTEP, 0.0);
                xr[m + N][n] = XReal(-m, n, GRADSTEP, 0.0);

                vi[m][n] = VImag(m, n, GRADSTEP, 0.0);
                vi[m + N][n] = VImag(-m, n, GRADSTEP, 0.0);
                wi[m][n] = WImag(m, n, GRADSTEP, 0.0);
                wi[m + N][n] = WImag(-m, n, GRADSTEP, 0.0);
                xi[m][n] = XImag(m, n, GRADSTEP, 0.0);
                xi[m + N][n] = XImag(-m, n, GRADSTEP, 0.0);
            }

        //calculate the inner products;

        //for V
        static double rhohatVr[2 * N + 1][N];
        static double rhohatVi[2 * N + 1][N];
        std::cout << "Computing inner product coefficiencts for V... " << std::endl;
        for (int n = 0; n < N; n++)
            for (int m = 0; m <= n; m++)
            {


                rhohatVr[m][n] = L2InnerProduct(&r, &vr[m][n], NUMGRIDS);

                rhohatVr[m + N][n] = L2InnerProduct(&r, &vr[m + N][n], NUMGRIDS);

                rhohatVi[m][n] = L2InnerProduct(&r, &vi[m][n], NUMGRIDS);

                rhohatVi[m + N][n] = L2InnerProduct(&r, &vi[m + N][n], NUMGRIDS);

                //std::cout << "|";
                std::cout << "m = " << m << ", " << "n = " << n << "   " << rhohatVr[m][n] << " " << rhohatVr[m + N][n] << std::endl
                    << std::setw(22) << rhohatVi[m][n] << " " << rhohatVi[m + N][n] << "\n\n";

            }
        std::cout << "Done!" << std::endl;

        //for W.
        static double rhohatWr[2 * N + 1][N];
        static double rhohatWi[2 * N + 1][N];
        std::cout << "Computing inner product coefficiencts for W... " << std::endl;
        for (int n = 0; n < N; n++)
            for (int m = 0; m <= n; m++)
            {


                rhohatWr[m][n] = L2InnerProduct(&r, &wr[m][n], NUMGRIDS);

                rhohatWr[m + N][n] = L2InnerProduct(&r, &wr[m + N][n], NUMGRIDS);

                rhohatWi[m][n] = L2InnerProduct(&r, &wi[m][n], NUMGRIDS);

                rhohatWi[m + N][n] = L2InnerProduct(&r, &wi[m + N][n], NUMGRIDS);
                //std::cout << "|";
                std::cout << "m = " << m << ", " << "n = " << n << "   " << rhohatWr[m][n] << " " << rhohatWr[m + N][n] << std::endl
                    << std::setw(22) << rhohatWi[m][n] << " " << rhohatWi[m + N][n] << "\n\n";

            }
        std::cout << "Done!" << std::endl;

        //and for X.
        static double rhohatXr[2 * N + 1][N];
        static double rhohatXi[2 * N + 1][N];
        std::cout << "Computing inner product coefficiencts for X... " << std::endl;
        for (int n = 0; n < N; n++)
            for (int m = 0; m <= n; m++)
            {


                rhohatXr[m][n] = L2InnerProduct(&r, &xr[m][n], NUMGRIDS);

                rhohatXr[m + N][n] = L2InnerProduct(&r, &xr[m + N][n], NUMGRIDS);

                rhohatXi[m][n] = L2InnerProduct(&r, &xi[m][n], NUMGRIDS);

                rhohatXi[m + N][n] = L2InnerProduct(&r, &xi[m + N][n], NUMGRIDS);
                //std::cout << "|";
                std::cout << "m = " << m << ", " << "n = " << n << "   " << rhohatXr[m][n] << " " << rhohatXr[m + N][n] << std::endl
                    << std::setw(22) << rhohatXi[m][n] << " " << rhohatXi[m + N][n] << "\n\n";

            }
        std::cout << "Done!" << std::endl;
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

                Vterm[m][n] = SphericalHarmonic(rhohatVr[m][n] ,rhohatVi[m][n], vr[m][n] , vi[m][n]);
                if (m > 0)
                    Vterm[m + N][n] = SphericalHarmonic(rhohatVr[m+N][n], rhohatVi[m+N][n], vr[m+N][n], vi[m+N][n]);


                Wterm[m][n] = SphericalHarmonic(rhohatWr[m][n], rhohatWi[m][n], wr[m][n], wi[m][n]);
                if (m > 0)
                    Wterm[m + N][n] = SphericalHarmonic(rhohatWr[m + N][n], rhohatWi[m + N][n], wr[m + N][n], wi[m + N][n]);



                Xterm[m][n] = SphericalHarmonic(rhohatXr[m][n], rhohatXi[m][n], xr[m][n], xi[m][n]);
                if (m > 0)
                    Xterm[m + N][n] = SphericalHarmonic(rhohatXr[m + N][n], rhohatXi[m + N][n], xr[m + N][n], xi[m + N][n]);

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

        double abserr = L2Difference(&r, &rhoapprox, NUMGRIDS);
        //Compute the error in approximating rho by it's expansion in Spherical Harmonics.
        std::cout << "Absolute Error in approximation(L2): ";
        std::cout << abserr<< std::endl;
        std::cout << std::endl;
    
        std::cout << "Relative Error in approximation(L2): ";
        std::cout << abserr / sqrt(L2InnerProduct(&r , &r , NUMGRIDS)) << std::endl;
        std::cout << std::endl;
           

    return 0;

}


