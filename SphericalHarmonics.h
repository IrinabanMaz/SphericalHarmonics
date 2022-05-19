#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <random>
#include <string>
#include <fstream>

#include "Legendre.h"
#include "SphericalCalc.h"

#define min(a , b) a < b? a : b
#define max(a , b) a > b? a : b






//Encapsulations of the Scalar Harmonic Function Y.
//How to use:
//1) declare and initialize a variable Y y(m , n );
//2) call the function as y(x); where x is a spherical coordinate point.

class YReal : public SurfaceScalarFunction
{
private:
    int m;
    int n;
    Legendre poly;

public:

    YReal() { m = 0; n = 0; }
    YReal(int a, int b)
    {
        m = a;
        n = b;
    }


    //Implements the function Y 
    double operator()(SurfaceCoord s)
    {



        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        coef *= poly(cos(s.theta), m, n);

        //std::cout << coef << std::endl;

        return coef * cos((double)m * s.phi);


    }

    double dphi(SurfaceCoord s)
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        coef *= poly(cos(s.theta), m, n);

        //std::cout << coef << std::endl;

        return - coef *m * sin((double)m * s.phi);
    }

    double dphidphi(SurfaceCoord s)
    {
        return -m * m * operator()(s);
    }


    double dtheta(SurfaceCoord s)
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly.populate(cos(s.theta), m + 1, n + 1);

        coef *= poly.dtheta(m, n, s.theta);


        return coef * cos((double)m * s.phi);
    }

    double dthetadtheta(SurfaceCoord s)
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly.populate(cos(s.theta), m + 1, n + 1);

        coef *= poly.dtheta(m, n, s.theta , 2);


        return coef * cos((double)m * s.phi);
    }

    double dthetadphi(SurfaceCoord s)
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly.populate(cos(s.theta), m + 1, n + 1);

        coef *= poly.dtheta(m, n, s.theta);


        return -   coef * m* sin((double)m * s.phi);
    }

};






class YImag : public SurfaceScalarFunction
{
private:
    int m;
    int n;
    Legendre poly;

public:

    YImag() { m = 0; n = 0; }
    YImag(int a, int b)
    {
        m = a;
        n = b;

    }


    //Implements the function Y 
    double operator()(SurfaceCoord s)
    {


        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        coef *= poly(cos(s.theta), m, n);

        //std::cout << coef << std::endl;

        return coef * sin((double)m * s.phi);


    }

    double dphi(SurfaceCoord s)
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        coef *= poly(cos(s.theta), m, n);

        //std::cout << coef << std::endl;

        return coef * m * cos((double)m * s.phi);
    }

    double dphidphi(SurfaceCoord s)
    {
        return -m * m * operator()(s);
    }

    double dtheta(SurfaceCoord s)
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly.populate(cos(s.theta), m + 1, n + 1);

        coef *= poly.dtheta(m, n, s.theta);


        return coef * sin((double)m * s.phi);
    }

    double dthetadtheta(SurfaceCoord s)
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly.populate(cos(s.theta), m + 1, n + 1);

        coef *= poly.dtheta(m, n, s.theta , 2);


        return coef * sin((double)m * s.phi);
    }

    double dthetadphi(SurfaceCoord s)
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly.populate(cos(s.theta), m + 1, n + 1);

        coef *= poly.dtheta(m, n, s.theta);


        return coef * m * cos((double)m * s.phi);
    }

};

class GReal : public SphericalVectorField
{
private:
    int m;
    int n;
    Legendre p;
    YReal yr;

public:
    GReal() {}
    GReal(int m0, int n0)
    {
        m = m0;
        n = n0;

        yr = YReal(m, n);
    }

    SphereCoord operator()(SphereCoord x)
    {

        //std::cout << coef << std::endl;

        double dtheta = yr.dtheta(x.s);

        double dphi = yr.dphi(x.s);

        return dtheta * e_theta(x) + dphi / sin(x.s.theta) * e_phi(x);


    }

    double eTheta(SphereCoord x)
    {
        return yr.dtheta(x.s);
    }

    double ePhi(SphereCoord x)
    {
        return yr.dphi(x.s) / sin(x.s.theta);
    }

    double eTheta_dTheta(SphereCoord x)
    {
        return yr.dthetadtheta(x.s);
    }

    double eTheta_dPhi(SphereCoord x)
    {
        return yr.dthetadphi(x.s);
    }

    double ePhi_dPhi(SphereCoord x)
    {
        return yr.dphidphi(x.s) / sin(x.s.theta);
    }

    double ePhi_dTheta(SphereCoord x)
    {
        return (yr.dthetadphi(x.s) * sin(x.s.theta) - cos(x.s.theta) * yr.dphi(x.s)) / (sin(x.s.theta) * sin(x.s.theta));
    }





};

class GImag : public SphericalVectorField
{
private:
    int m;
    int n;
    Legendre p;
    YImag yi;

public:
    GImag() {}
    GImag(int m0, int n0)
    {
        m = m0;
        n = n0;

        yi = YImag(m, n);
    }

    double eTheta(SphereCoord x)
    {
        return yi.dtheta(x.s);
    }

    double ePhi(SphereCoord x)
    {
        return yi.dphi(x.s) / sin(x.s.theta);
    }

    SphereCoord operator()(SphereCoord x)
    {



        //std::cout << coef << std::endl;

        double dtheta = yi.dtheta(x.s);

        double dphi = yi.dphi(x.s);

        return dtheta * e_theta(x) + dphi / sin(x.s.theta) * e_phi(x);


    }


    double eTheta_dTheta(SphereCoord x)
    {
        return yi.dthetadtheta(x.s);
    }

    double eTheta_dPhi(SphereCoord x)
    {
        return yi.dthetadphi(x.s);
    }

    double ePhi_dPhi(SphereCoord x)
    {
        return yi.dphidphi(x.s) / sin(x.s.theta);
    }

    double ePhi_dTheta(SphereCoord x)
    {
        return (yi.dthetadphi(x.s) * sin(x.s.theta) - cos(x.s.theta) * yi.dphi(x.s)) / (sin(x.s.theta) * sin(x.s.theta));
    }


};

//Spherical Harmonic Basis function encapsulations.
class VReal : public SphericalVectorField
{
    int m;
    int n;
    YReal yr;
    GReal gr;

public:

    VReal() {}
    VReal(int a, int b)
    {
        m = a;
        n = b;
        yr = YReal(m, n);
        gr = GReal(m, n);
    }

    SphereCoord operator()(SphereCoord x)
    {
        return -(double)(n + 1) * yr(x.s) * e_r(x) + gr(x);
    }

    void reset(int a, int b)
    {
        m = a;
        n = b;
        yr = YReal(m, n);
        gr = GReal(m, n);
    }


    double eR(SphereCoord x)
    {
        return -(double)(n + 1) * yr(x.s);
    }

    double eTheta(SphereCoord x)
    {
        return yr.dtheta(x.s);
    }

    double ePhi(SphereCoord x)
    {
        return yr.dphi(x.s) / sin(x.s.theta);
    }

    double eTheta_dTheta(SphereCoord x)
    {
        return yr.dthetadtheta(x.s);
    }

    double eTheta_dPhi(SphereCoord x)
    {
        return yr.dthetadphi(x.s);
    }

    double ePhi_dPhi(SphereCoord x)
    {
        return yr.dphidphi(x.s) / sin(x.s.theta);
    }

    double ePhi_dTheta(SphereCoord x)
    {
        return (yr.dthetadphi(x.s) * sin(x.s.theta) - cos(x.s.theta) * yr.dphi(x.s)) / (sin(x.s.theta) * sin(x.s.theta));
    }

    double eR_dTheta(SphereCoord x)
    {
        return yr.dtheta(x.s);
    }

    double eR_dPhi(SphereCoord x)
    {
        return yr.dphi(x.s);
    }


};

class VImag : public SphericalVectorField
{
    int m;
    int n;
    YImag yi;
    GImag gi;

public:

    VImag() {}
    VImag(int a, int b)
    {
        m = a;
        n = b;
        yi = YImag(m, n);
        gi = GImag(m, n);
    }

    SphereCoord operator()(SphereCoord x)
    {
        return -(double)(n + 1) * yi(x.s) * e_r(x) + gi(x);
    }



    double eR(SphereCoord x)
    {
        return -(double)(n + 1) * yi(x.s);
    }

    double eTheta(SphereCoord x)
    {
        return yi.dtheta(x.s);
    }

    double ePhi(SphereCoord x)
    {
        return yi.dphi(x.s) / sin(x.s.theta);
    }

    double eTheta_dTheta(SphereCoord x)
    {
        return yi.dthetadtheta(x.s);
    }

    double eTheta_dPhi(SphereCoord x)
    {
        return yi.dthetadphi(x.s);
    }

    double ePhi_dPhi(SphereCoord x)
    {
        return yi.dphidphi(x.s) / sin(x.s.theta);
    }

    double ePhi_dTheta(SphereCoord x)
    {
        return (yi.dthetadphi(x.s) * sin(x.s.theta) - cos(x.s.theta) * yi.dphi(x.s)) / (sin(x.s.theta) * sin(x.s.theta));
    }

    double eR_dTheta(SphereCoord x)
    {
        return yi.dtheta(x.s);
    }

    double eR_dPhi(SphereCoord x)
    {
        return yi.dphi(x.s);
    }


    void reset(int a, int b)
    {
        m = a;
        n = b;
        yi = YImag(m, n);
        gi = GImag(m, n);
    }

};



class WReal : public SphericalVectorField
{
    int m;
    int n;
    YReal yr;
    GReal gr;

public:

    WReal() {}
    WReal(int a, int b)
    {
        m = a;
        n = b;
        yr = YReal(m, n);
        gr = GReal(m, n);
    }

    SphereCoord operator()(SphereCoord x)
    {
        return (double)n * yr(x.s) * e_r(x) + gr(x);
    }

    void reset(int a, int b)
    {
        m = a;
        n = b;
        yr = YReal(m, n);
        gr = GReal(m, n);
    }

    double eR(SphereCoord x)
    {
        return n * yr(x.s);
    }

    double eTheta(SphereCoord x)
    {
        return yr.dtheta(x.s);
    }

    double ePhi(SphereCoord x)
    {
        return yr.dphi(x.s) / sin(x.s.theta);
    }

    double eTheta_dTheta(SphereCoord x)
    {
        return yr.dthetadtheta(x.s);
    }

    double eTheta_dPhi(SphereCoord x)
    {
        return yr.dthetadphi(x.s);
    }

    double ePhi_dPhi(SphereCoord x)
    {
        return yr.dphidphi(x.s) / sin(x.s.theta);
    }

    double ePhi_dTheta(SphereCoord x)
    {
        return (yr.dthetadphi(x.s) * sin(x.s.theta) - cos(x.s.theta) * yr.dphi(x.s)) / (sin(x.s.theta) * sin(x.s.theta));
    }
    double eR_dTheta(SphereCoord x)
    {
        return yr.dtheta(x.s);
    }

    double eR_dPhi(SphereCoord x)
    {
        return yr.dphi(x.s);
    }
};

class WImag : public SphericalVectorField
{
    int m;
    int n;
    YImag yi;
    GImag gi;

public:

    WImag() {}
    WImag(int a, int b)
    {
        m = a;
        n = b;
        yi = YImag(m, n);
        gi = GImag(m, n);
    }

    SphereCoord operator()(SphereCoord x)
    {
        return (double)n * yi(x.s) * e_r(x) + gi(x);
    }

    void reset(int a, int b)
    {
        m = a;
        n = b;
        yi = YImag(m, n);
        gi = GImag(m, n);
    }

    double eR(SphereCoord x)
    {
        return -(double)(n + 1) * yi(x.s);
    }

    double eTheta(SphereCoord x)
    {
        return yi.dtheta(x.s);
    }

    double ePhi(SphereCoord x)
    {
        return yi.dphi(x.s) / sin(x.s.theta);
    }

    double eTheta_dTheta(SphereCoord x)
    {
        return yi.dthetadtheta(x.s);
    }

    double eTheta_dPhi(SphereCoord x)
    {
        return yi.dthetadphi(x.s);
    }

    double ePhi_dPhi(SphereCoord x)
    {
        return yi.dphidphi(x.s) / sin(x.s.theta);
    }

    double ePhi_dTheta(SphereCoord x)
    {
        return (yi.dthetadphi(x.s) * sin(x.s.theta) - cos(x.s.theta) * yi.dphi(x.s)) / (sin(x.s.theta) * sin(x.s.theta));
    }

    double eR_dTheta(SphereCoord x)
    {
        return yi.dtheta(x.s);
    }

    double eR_dPhi(SphereCoord x)
    {
        return yi.dphi(x.s);
    }

};




class XReal : public SphericalVectorField
{
    int m;
    int n;
    GReal gr;
    YReal yr;

public:

    XReal() {}
    XReal(int a, int b)
    {
        m = a;
        n = b;

        yr = YReal(m, n);
        gr = GReal(m, n);
    }

    SphereCoord operator()(SphereCoord x)
    {
        return cross(e_r(x), gr(x));
    }

    void reset(int a, int b)
    {
        m = a;
        n = b;
        gr = GReal(m, n);
    }

    double eR(SphereCoord x)
    {
        return yr.dphi(x.s) / sin(x.s.theta);
    }

    double eTheta(SphereCoord x)
    {
        return 0;
    }

    double ePhi(SphereCoord x)
    {
        return yr.dtheta(x.s);
    }

    double eTheta_dTheta(SphereCoord x)
    {
        return 0;
    }

    double eTheta_dPhi(SphereCoord x)
    {
        return 0;
    }

    double ePhi_dPhi(SphereCoord x)
    {
        return yr.dtheta(x.s);
    }
    
    double ePhi_dTheta(SphereCoord x)
    {
        return yr.dthetadtheta(x.s);
    }

    double eR_dTheta(SphereCoord x)
    {
        return (yr.dthetadphi(x.s) * sin(x.s.theta) - cos(x.s.theta) * yr.dphi(x.s)) / (sin(x.s.theta) * sin(x.s.theta));
    }

    double eR_dPhi(SphereCoord x)
    {
        return yr.dthetadphi(x.s);
    }
};

class XImag : public SphericalVectorField
{
    int m;
    int n;
    YImag yi;
    GImag gi;

public:

    XImag() {}
    XImag(int a, int b)
    {
        m = a;
        n = b;

        yi = YImag(m, n);
        gi = GImag(m, n);
    }

    SphereCoord operator()(SphereCoord x)
    {
        return cross(e_r(x), gi(x));
    }

    void reset(int a, int b)
    {
        m = a;
        n = b;
        gi = GImag(m, n);
    }
    double eR(SphereCoord x)
    {
        return yi.dphi(x.s) / sin(x.s.theta);
    }

    double eTheta(SphereCoord x)
    {
        return 0;
    }

    double ePhi(SphereCoord x)
    {
        return yi.dtheta(x.s);
    }

    double eTheta_dTheta(SphereCoord x)
    {
        return 0;
    }

    double eTheta_dPhi(SphereCoord x)
    {
        return 0;
    }

    double ePhi_dPhi(SphereCoord x)
    {
        return yi.dtheta(x.s);
    }

    double ePhi_dTheta(SphereCoord x)
    {
        return yi.dthetadtheta(x.s);
    }

    double eR_dTheta(SphereCoord x)
    {
        return (yi.dthetadphi(x.s) * sin(x.s.theta) - cos(x.s.theta) * yi.dphi(x.s)) / (sin(x.s.theta) * sin(x.s.theta));
    }

    double eR_dPhi(SphereCoord x)
    {
        return yi.dthetadphi(x.s);
    }

};

//Encapsulation of a single term in the expansion of a function.
class SphericalHarmonic : public SphericalVectorField
{
private:


    int m;
    int n;
    VReal vr;
    VImag vi;
    WReal wr;
    WImag wi;
    XReal xr;
    XImag xi;
public:
    double rhohatVr;
    double rhohatVi;
    double rhohatWr;
    double rhohatWi;
    double rhohatXr;
    double rhohatXi;



    SphericalHarmonic() {}

    SphericalHarmonic(int m0, int n0, std::array<double, 6> rhohats)
    {
        m = m0;
        n = n0;

        vr = VReal(m, n);
        vi = VImag(m, n);
        wr = WReal(m, n);
        wi = WImag(m, n);
        xr = XReal(m, n);
        xi = XImag(m, n);

        rhohatVr = rhohats[0];
        rhohatVi = rhohats[1];
        rhohatWr = rhohats[2];
        rhohatWi = rhohats[3];
        rhohatXr = rhohats[4];
        rhohatXi = rhohats[5];

    }

    SphereCoord operator()(SphereCoord x)
    {
        if (n == 0)
            return rhohatVr * vr(x) + rhohatVi * vi(x);

        SphereCoord temp = rhohatVr * vr(x) + rhohatVi * vi(x);
        temp = temp + rhohatWr * wr(x) + rhohatWi * wi(x);
        return temp + rhohatXr * xr(x) + rhohatXi * xi(x);
    }

};

class VSHSeries : public VectorFieldSum
{

private:



    Legendre P;



public:
    int N;
    int NUMGRIDS;
    std::vector<std::vector<SphericalHarmonic>> terms;




    VSHSeries(int n, int numgrids)
    {
        N = n;
        NUMGRIDS = numgrids;

        terms.resize(n + 1);

        for (int i = 0; i <= n; i++)
            terms[i].resize(2 * i + 1);



    }

    void approximate(SphericalVectorField* r)
    {

        std::array<double, 6> rhohats;
        VReal vr;
        VImag vi;
        WReal wr;
        WImag wi;
        XReal xr;
        XImag xi;

        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                if (m == 0)
                {
                    vr.reset(m, n);
                    rhohats[0] = L2InnerProduct(r, &vr, NUMGRIDS) / (double)((2 * n + 1) * (n + 1));

                    vi.reset(m, n);
                    rhohats[1] = L2InnerProduct(r, &vi, NUMGRIDS) / (double)((2 * n + 1) * (n + 1));

                    wr.reset(m, n);
                    rhohats[2] = L2InnerProduct(r, &wr, NUMGRIDS) / (double)((2 * n + 1) * n);

                    wi.reset(m, n);
                    rhohats[3] = L2InnerProduct(r, &wi, NUMGRIDS) / (double)((2 * n + 1) * n);

                    xr.reset(m, n);
                    rhohats[4] = L2InnerProduct(r, &xr, NUMGRIDS) / (double)((n + 1) * n);

                    xi.reset(m, n);
                    rhohats[5] = L2InnerProduct(r, &xi, NUMGRIDS) / (double)((n + 1) * n);

                    terms[n][m] = SphericalHarmonic(m, n, rhohats);
                }
                else if (m > 0)
                {
                    vr.reset(m, n);
                    rhohats[0] = 2.0 * L2InnerProduct(r, &vr, NUMGRIDS) / (double)((2 * n + 1) * (n + 1));

                    vi.reset(m, n);
                    rhohats[1] = 2.0 * L2InnerProduct(r, &vi, NUMGRIDS) / (double)((2 * n + 1) * (n + 1));

                    wr.reset(m, n);
                    rhohats[2] = 2.0 * L2InnerProduct(r, &wr, NUMGRIDS) / (double)((2 * n + 1) * n);

                    wi.reset(m, n);
                    rhohats[3] = 2.0 * L2InnerProduct(r, &wi, NUMGRIDS) / (double)((2 * n + 1) * n);

                    xr.reset(m, n);
                    rhohats[4] = 2.0 * L2InnerProduct(r, &xr, NUMGRIDS) / (double)((n + 1) * n);

                    xi.reset(m, n);
                    rhohats[5] = 2.0 * L2InnerProduct(r, &xi, NUMGRIDS) / (double)((n + 1) * n);

                    terms[n][m] = SphericalHarmonic(m, n, rhohats);
                }
                /*
                   std::cout << "Coefficients for m = " << m << ", n = " << n << "\n";
                   std::cout << "Real part of V " << rhohats[0] << "\n";
                   std::cout << "Imaginary part of V " << rhohats[1] << "\n\n";
                   std::cout << "Real part of W " << rhohats[2] << "\n";
                   std::cout << "Imaginary part of W " << rhohats[3] << "\n\n";
                   std::cout << "Real part of X " << rhohats[4] << "\n";
                   std::cout << "Imaginary part of X " << rhohats[5] << "\n\n\n";
                   */
                append(terms[n][m]);
                /*
                if (m > 0)
                {
                    vr.reset(-m, n);
                    rhohats[0] = 2.0*L2InnerProduct(r, &vr, NUMGRIDS) / (double)((2 * n + 1)*(n+1));

                    vi.reset(-m, n);
                    rhohats[1] = 2.0*L2InnerProduct(r, &vi, NUMGRIDS) / (double)((2 * n + 1) * (n + 1));

                    wr.reset(-m, n);
                    rhohats[2] = 2.0*L2InnerProduct(r, &wr, NUMGRIDS) / (double)((2 * n + 1) *n);

                    wi.reset(-m, n);
                    rhohats[3] = 2.0*L2InnerProduct(r, &wi, NUMGRIDS) / (double)((2 * n + 1) * n);

                    xr.reset(-m, n);
                    rhohats[4] = 2.0*L2InnerProduct(r, &xr, NUMGRIDS) / (double)((n + 1) * n);

                    xi.reset(-m, n);
                    rhohats[5] = 2.0*L2InnerProduct(r, &xi, NUMGRIDS) / (double)((n + 1) * n);

                    terms[n][m + n] = SphericalHarmonic(-m, n, rhohats);

                    append(terms[n][m + n]);
                }*/
            }


        //calculate the inner products;

    //for V

     //   std::cout << "Computing inner product coefficiencts for V... " << std::endl;
 /*       std::cout << "\\begin{center}\n";
        std::cout << "\\begin{tabular}{ | c ||";
            for (int i = 1; i <= N; i++)
                std::cout << "c |";
        std::cout << "}\n";
        std::cout << "\\hline \n";
        std::cout << "$V_n^m$ & m = 0 &";
        for (int i = 1; i < N - 1; i++)
            std::cout <<"m = " << i << " & ";
        std::cout << N - 1 << "\\\\ \\hline \\hline \n"; */
        /*
        for (int n = 0; n <= N; n++)
        {
            //std::cout << "n = " << n;
            for (int m = 0; m <= n; m++)
            {


                rhohatVr[m][n] = L2InnerProduct(&r, &vr[m][n], NUMGRIDS);

                rhohatVr[m + N][n] = L2InnerProduct(&r, &vr[m + N][n], NUMGRIDS);

                rhohatVi[m][n] = L2InnerProduct(&r, &vi[m][n], NUMGRIDS);

                rhohatVi[m + N][n] = L2InnerProduct(&r, &vi[m + N][n], NUMGRIDS);

                //std::cout << "|";
                std::cout << "m = " << m << ", " << "n = " << n << "  Inner product of rho with V_n^m real part:  " << rhohatVr[m][n] << "  || V_n^-m" << rhohatVr[m + N][n] << std::endl
                    << "              Inner product of rho with V_n^m imaginary part   :" << rhohatVi[m][n] << "  ||  V_n^-m " << rhohatVi[m + N][n] << "\n\n";

                //std::cout << " & $" << rhohatVr[m][n] << " + i " << rhohatVi[m][n] << "$";

            }
           // std::cout << "\\\\ \\hline \n";
        }
       // std::cout << "\\end{tabular}\n";
       // std::cout << "\\end{center}\n";
        std::cout << "Done!" << std::endl;

        //for W.
        double** rhohatWr = new double* [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            rhohatWr[i] = new double[N + 1];

        double** rhohatWi = new double* [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            rhohatWi[i] = new double[N + 1];
        std::cout << "Computing inner product coefficiencts for W... " << std::endl;
        for (int n = 0; n <= N; n++)
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
        double** rhohatXr = new double* [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            rhohatXr[i] = new double[N + 1];

        double** rhohatXi = new double* [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            rhohatXi[i] = new double[N + 1];
        std::cout << "Computing inner product coefficiencts for X... " << std::endl;
        for (int n = 0; n <= N; n++)
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

        //construct the summation. Note we must compute the L2 norm of the basis functions to do so, so this takes a few seconds..
        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {

                Vterm[m][n] = SphericalHarmonic(rhohatVr[m][n], rhohatVi[m][n], vr[m][n], vi[m][n]);
                Vterm[m][n].denom = (double)((2 * n + 1) * (n + 1));
                if (m > 0)
                {
                    Vterm[m + N][n] = SphericalHarmonic(rhohatVr[m + N][n], rhohatVi[m + N][n], vr[m + N][n], vi[m + N][n]);
                    Vterm[m + N][n].denom = (double)((2 * n + 1) * (n + 1));
                }

                Wterm[m][n] = SphericalHarmonic(rhohatWr[m][n], rhohatWi[m][n], wr[m][n], wi[m][n]);
                Wterm[m][n].denom = (double)(n * (2 * n + 1));
                if (m > 0)
                {
                    Wterm[m + N][n] = SphericalHarmonic(rhohatWr[m + N][n], rhohatWi[m + N][n], wr[m + N][n], wi[m + N][n]);
                    Wterm[m + N][n].denom = (double)(n * (2 * n + 1));
                }


                Xterm[m][n] = SphericalHarmonic(rhohatXr[m][n], rhohatXi[m][n], xr[m][n], xi[m][n]);
                Xterm[m][n].denom = (double)(n * (n + 1));
                if (m > 0)
                {
                    Xterm[m + N][n] = SphericalHarmonic(rhohatXr[m + N][n], rhohatXi[m + N][n], xr[m + N][n], xi[m + N][n]);
                    Xterm[m + N][n].denom = (double)(n * (n + 1));
                }
               append(Vterm[m][n]);
                if (m > 0)
                   append(Vterm[m + N][n]);

                append(Wterm[m][n]);
                if (m > 0)
                    append(Wterm[m + N][n]);

                append(Xterm[m][n]);
                if (m > 0)
                    append(Xterm[m + N][n]);


                std::cout << "|";

            }*/

            //std::cout << "Done!" << std::endl << std::endl;



    }




    SphereCoord operator()(SphereCoord x)
    {
        P.populate(cos(x.s.theta), N, N);

        return VectorFieldSum::operator()(x);
    }
};