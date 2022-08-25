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
    YReal() { m = 0; n = 0; name = "YReal"; }
    YReal(int a, int b)
    {
        m = a;
        n = b;
        name = "YReal";
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

    YImag() { m = 0; n = 0; name = "YImag"; }
    YImag(int a, int b)
    {
        m = a;
        n = b;
        name = "YImag";

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

class Yrtest : public SphericalScalarFunction
{
private:
    YReal yr;
public:
    double operator()(SphereCoord x)
    {
        return yr(x.s);
    }


    void testhelper(int m ,int n)
    {
            
            yr = YReal(m, n);
            NdPhi dphi(this);
            NdPhi dphidphi(&dphi);
            NdTheta dtheta(this);
            NdTheta dthetadtheta(&dtheta);
            NdTheta dphidtheta(&dphi);


            double dphierr = 0.0;
            double dphidphierr = 0.0;
            double dthetaerr = 0.0;
            double dthetadthetaerr = 0.0;
            double dphidthetaerr = 0.0;

            for (int p = 0; p < NUMTRAPNODES; p++)
                for (int i = 0; i < NUMGLNODES; i++)
                {
                    SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
                    SphereCoord x(1, s);
                    double diff = dphi(x) - yr.dphi(x.s);
                    dphierr +=sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                    diff = dphidphi(x) - yr.dphidphi(x.s);
                    dphidphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                    diff = dtheta(x) - yr.dtheta(x.s);
                    dthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                    diff = dthetadtheta(x) - yr.dthetadtheta(x.s);
                    dthetadthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                    diff = dphidtheta(x) - yr.dthetadphi(x.s);
                    dphidthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                    /*if (dot(diff, diff) > 1e-12)
                    {
                        std::cout << " error of " << dot(diff, diff) << " at " << x << std::endl;
                        std::cout << "value of f: " << (*f1)(x) << std::endl;
                        std::cout << "value of g: " << (*f2)(x) << std::endl;

                    }
                    */
                }

            std::cout << "Error in computing dphi: " << dphierr << "\n";
            std::cout << "Error in computing dphidphi: " << dphidphierr << "\n";
            std::cout << "Error in computing dtheta: " << dthetaerr << "\n";
            std::cout << "Error in computing dthetadtheta: " << dthetadthetaerr << "\n";
            std::cout << " Error in computing dthetadphi: " << dphidthetaerr << "\n";
        }
    

};

class Yitest : public SphericalScalarFunction
{
private:
    YImag yi;
public:
    double operator()(SphereCoord x)
    {
        return yi(x.s);
    }


    void testhelper(int m, int n)
    {
       
        yi = YImag(m, n);
        NdPhi dphi(this);
        NdPhi dphidphi(&dphi);
        NdTheta dtheta(this);
        NdTheta dthetadtheta(&dtheta);
        NdTheta dphidtheta(&dphi);


        double dphierr = 0.0;
        double dphidphierr = 0.0;
        double dthetaerr = 0.0;
        double dthetadthetaerr = 0.0;
        double dphidthetaerr = 0.0;

        for (int p = 0; p < NUMTRAPNODES; p++)
            for (int i = 0; i < NUMTRAPNODES; i++)
            {
                SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
                SphereCoord x(1, s);
                double diff = dphi(x) - yi.dphi(x.s);
                dphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = dphidphi(x) - yi.dphidphi(x.s);
                dphidphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = dtheta(x) - yi.dtheta(x.s);
                dthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = dthetadtheta(x) - yi.dthetadtheta(x.s);
                dthetadthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = dphidtheta(x) - yi.dthetadphi(x.s);
                dphidthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                /*if (dot(diff, diff) > 1e-12)
                {
                    std::cout << " error of " << dot(diff, diff) << " at " << x << std::endl;
                    std::cout << "value of f: " << (*f1)(x) << std::endl;
                    std::cout << "value of g: " << (*f2)(x) << std::endl;

                }
                */
            }

        std::cout << "Error in computing dphi: " << dphierr << "\n";
        std::cout << "Error in computing dphidphi: " << dphidphierr << "\n";
        std::cout << "Error in computing dtheta: " << dthetaerr << "\n";
        std::cout << "Error in computing dthetadtheta: " << dthetadthetaerr << "\n";
        std::cout << " Error in computing dthetadphi: " << dphidthetaerr << "\n";
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

    RectCoord operator()(SphereCoord x)
    {
        const double coef = sqrt(0.375 / PI);
       /*
        if ((n == 1) && (m == 1))
        {
            return coef * (cos(x.s.theta) * cos(x.s.phi) * e_theta(x) + sin(x.s.phi) * e_phi(x));
        }
        */
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


    void testhelper()
    {
        dot_ePhi dep(this);
        dot_eTheta det(this);

        NdPhi epdp(&dep);
        NdPhi etdp(&det);

        NdTheta epdt(&dep);
        NdTheta etdt(&det);


        double ephierr = 0.0;
        double ethetaerr = 0.0;
        double ephidphierr = 0.0;
        double ephidthetaerr = 0.0;
        double ethetadthetaerr = 0.0;
        double ethetadphierr = 0.0;

        for (int p = 0; p < NUMTRAPNODES; p++)
            for (int i = 0; i < NUMTRAPNODES; i++)
            {
                SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
                SphereCoord x(1, s);

                double diff = ePhi(x) - dep(x);
                ephierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta(x) - det(x);
                ethetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dPhi(x) - epdp(x);
                ephidphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dTheta(x) - etdt(x);
                ethetadthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                /*if (dot(diff, diff) > 1e-12)
                {
                    std::cout << " error of " << dot(diff, diff) << " at " << x << std::endl;
                    std::cout << "value of f: " << (*f1)(x) << std::endl;
                    std::cout << "value of g: " << (*f2)(x) << std::endl;

                }
                */
            }

        std::cout << "Error in computing ephi: " << ephierr << "\n";
        std::cout << "Error in computing ephidphi: " << ephidphierr << "\n";
        std::cout << "Error in computing etheta: " << ethetaerr << "\n";
        std::cout << "Error in computing ethetadtheta: " << ethetadthetaerr << "\n";
        std::cout << " Error in computing ethetadphi: " << ethetadphierr << "\n";
        std::cout << " Error in computing ephidtheta: " << ephidthetaerr << "\n";
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

    

    RectCoord operator()(SphereCoord x)
    {


        const double coef = sqrt(0.375 / PI);
        /*
        if ((n == 1) && (m == 1))
        {
            return coef * (cos(x.s.theta) * sin(x.s.phi) * e_theta(x) - cos(x.s.phi) * e_phi(x));
        }
        */
        double dtheta = yi.dtheta(x.s);

        double dphi = yi.dphi(x.s);

        return dtheta * e_theta(x) + dphi / sin(x.s.theta) * e_phi(x);


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


    void testhelper()
    {
        dot_ePhi dep(this);
        dot_eTheta det(this);

        NdPhi epdp(&dep);
        NdPhi etdp(&det);

        NdTheta epdt(&dep);
        NdTheta etdt(&det);


        double ephierr = 0.0;
        double ethetaerr = 0.0;
        double ephidphierr = 0.0;
        double ephidthetaerr = 0.0;
        double ethetadthetaerr = 0.0;
        double ethetadphierr = 0.0;

        for (int p = 0; p < NUMTRAPNODES; p++)
            for (int i = 0; i < NUMTRAPNODES; i++)
            {
                SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
                SphereCoord x(1, s);

                double diff = ePhi(x) - dep(x);
                ephierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta(x) - det(x);
                ethetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dPhi(x) - epdp(x);
                ephidphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dTheta(x) - etdt(x);
                ethetadthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                /*if (dot(diff, diff) > 1e-12)
                {
                    std::cout << " error of " << dot(diff, diff) << " at " << x << std::endl;
                    std::cout << "value of f: " << (*f1)(x) << std::endl;
                    std::cout << "value of g: " << (*f2)(x) << std::endl;

                }
                */
            }

        std::cout << "Error in computing ephi: " << ephierr << "\n";
        std::cout << "Error in computing ephidphi: " << ephidphierr << "\n";
        std::cout << "Error in computing etheta: " << ethetaerr << "\n";
        std::cout << "Error in computing ethetadtheta: " << ethetadthetaerr << "\n";
        std::cout << " Error in computing ethetadphi: " << ethetadphierr << "\n";
        std::cout << " Error in computing ephidtheta: " << ephidthetaerr << "\n";
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

    RectCoord operator()(SphereCoord x)
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
        return gr.eTheta(x);
    }

    double ePhi(SphereCoord x)
    {
        return gr.ePhi(x);
    }

    double eTheta_dTheta(SphereCoord x)
    {
        return gr.eTheta_dTheta(x);
    }

    double eTheta_dPhi(SphereCoord x)
    {
        return gr.eTheta_dPhi(x);
    }

    double ePhi_dPhi(SphereCoord x)
    {
        return gr.ePhi_dPhi(x);
    }

    double ePhi_dTheta(SphereCoord x)
    {
        return gr.ePhi_dTheta(x);
    }

    double eR_dTheta(SphereCoord x)
    {
        return -(double)(n+1)*yr.dtheta(x.s);
    }

    double eR_dPhi(SphereCoord x)
    {
        return -(double)(n + 1) * yr.dphi(x.s);
    }

    void testhelper()
    {
        dot_ePhi dep(this);
        dot_eTheta det(this);
        dot_eR der(this);

        NdPhi epdp(&dep);
        NdPhi etdp(&det);
        NdPhi erdp(&der);

        NdTheta epdt(&dep);
        NdTheta etdt(&det);
        NdTheta erdt(&der);



        double ephierr = 0.0;
        double ethetaerr = 0.0;
        double ererr = 0.0;
        double ephidphierr = 0.0;
        double ephidthetaerr = 0.0;
        double erdphierr = 0.0;
        double ethetadthetaerr = 0.0;
        double ethetadphierr = 0.0;
        double erdthetaerr = 0.0;


        for (int p = 0; p < NUMTRAPNODES; p++)
            for (int i = 0; i < NUMTRAPNODES; i++)
            {
                SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
                SphereCoord x(1, s);

                double diff = ePhi(x) - dep(x);
                ephierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta(x) - det(x);
                ethetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                
                diff = eR(x) - der(x);
                ethetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dPhi(x) - epdp(x);
                ephidphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                

                diff = eTheta_dTheta(x) - etdt(x);
                ethetadthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                
                diff = eR_dTheta(x) - erdt(x);
                erdthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dPhi(x) - erdp(x);
                erdphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                /*if (dot(diff, diff) > 1e-12)
                {
                    std::cout << " error of " << dot(diff, diff) << " at " << x << std::endl;
                    std::cout << "value of f: " << (*f1)(x) << std::endl;
                    std::cout << "value of g: " << (*f2)(x) << std::endl;

                }
                */
            }

        std::cout << "Error in computing ephi: " << ephierr << "\n";
        std::cout << "Error in computing etheta: " << ethetaerr << "\n";
        std::cout << "Error in computing er: " << ererr << "\n";

        std::cout << "Error in computing ephidphi: " << ephidphierr << "\n";
        std::cout << "Error in computing ephidtheta: " << ephidthetaerr << "\n";

        std::cout << "Error in computing ethetadtheta: " << ethetadthetaerr << "\n";
        std::cout << "Error in computing ethetadphi: " << ethetadphierr << "\n";

        std::cout << "Error in computing erdtheta: " << erdthetaerr << "\n";
        std::cout << "Error in computing erdphi: " << erdphierr << "\n";
        
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

    RectCoord operator()(SphereCoord x)
    {
        return -(double)(n + 1) * yi(x.s) * e_r(x) + gi(x);
    }



    double eR(SphereCoord x)
    {
        return -(double)(n + 1) * yi(x.s);
    }

    double eTheta(SphereCoord x)
    {
        return gi.eTheta(x);
    }

    double ePhi(SphereCoord x)
    {
        return gi.ePhi(x);
    }

    double eTheta_dTheta(SphereCoord x)
    {
        return gi.eTheta_dTheta(x);
    }

    double eTheta_dPhi(SphereCoord x)
    {
        return gi.eTheta_dPhi(x);
    }

    double ePhi_dPhi(SphereCoord x)
    {
        return gi.ePhi_dPhi(x);
    }

    double ePhi_dTheta(SphereCoord x)
    {
        return gi.ePhi_dTheta(x);
    }

    double eR_dTheta(SphereCoord x)
    {
        return -(double)(n + 1) * yi.dtheta(x.s);
    }

    double eR_dPhi(SphereCoord x)
    {
        return -(double)(n + 1) * yi.dphi(x.s);
    }


    void reset(int a, int b)
    {
        m = a;
        n = b;
        yi = YImag(m, n);
        gi = GImag(m, n);
    }

    void testhelper()
    {
        dot_ePhi dep(this);
        dot_eTheta det(this);
        dot_eR der(this);

        NdPhi epdp(&dep);
        NdPhi etdp(&det);
        NdPhi erdp(&der);

        NdTheta epdt(&dep);
        NdTheta etdt(&det);
        NdTheta erdt(&der);



        double ephierr = 0.0;
        double ethetaerr = 0.0;
        double ererr = 0.0;
        double ephidphierr = 0.0;
        double ephidthetaerr = 0.0;
        double erdphierr = 0.0;
        double ethetadthetaerr = 0.0;
        double ethetadphierr = 0.0;
        double erdthetaerr = 0.0;


        for (int p = 0; p < NUMTRAPNODES; p++)
            for (int i = 0; i < NUMTRAPNODES; i++)
            {
                SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
                SphereCoord x(1, s);

                double diff = ePhi(x) - dep(x);
                ephierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta(x) - det(x);
                ethetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR(x) - der(x);
                ethetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dPhi(x) - epdp(x);
                ephidphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;



                diff = eTheta_dTheta(x) - etdt(x);
                ethetadthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dTheta(x) - erdt(x);
                erdthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dPhi(x) - erdp(x);
                erdphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                /*if (dot(diff, diff) > 1e-12)
                {
                    std::cout << " error of " << dot(diff, diff) << " at " << x << std::endl;
                    std::cout << "value of f: " << (*f1)(x) << std::endl;
                    std::cout << "value of g: " << (*f2)(x) << std::endl;

                }
                */
            }

        std::cout << "Error in computing ephi: " << ephierr << "\n";
        std::cout << "Error in computing etheta: " << ethetaerr << "\n";
        std::cout << "Error in computing er: " << ererr << "\n";

        std::cout << "Error in computing ephidphi: " << ephidphierr << "\n";
        std::cout << "Error in computing ephidtheta: " << ephidthetaerr << "\n";

        std::cout << "Error in computing ethetadtheta: " << ethetadthetaerr << "\n";
        std::cout << "Error in computing ethetadphi: " << ethetadphierr << "\n";

        std::cout << "Error in computing erdtheta: " << erdthetaerr << "\n";
        std::cout << "Error in computing erdphi: " << erdphierr << "\n";

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

    RectCoord operator()(SphereCoord x)
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
        return gr.eTheta(x.s);
    }

    double ePhi(SphereCoord x)
    {
        return gr.ePhi(x);
    }

    double eTheta_dTheta(SphereCoord x)
    {
        return gr.eTheta_dTheta(x);
    }

    double eTheta_dPhi(SphereCoord x)
    {
        return gr.eTheta_dPhi(x);
    }

    double ePhi_dPhi(SphereCoord x)
    {
        return gr.ePhi_dPhi(x);
    }

    double ePhi_dTheta(SphereCoord x)
    {
        return gr.ePhi_dTheta(x);
    }
    double eR_dTheta(SphereCoord x)
    {
        return (double)n * yr.dtheta(x.s);
    }

    double eR_dPhi(SphereCoord x)
    {
        return (double)n * yr.dphi(x.s);
    }

    void testhelper()
    {
        dot_ePhi dep(this);
        dot_eTheta det(this);
        dot_eR der(this);

        NdPhi epdp(&dep);
        NdPhi etdp(&det);
        NdPhi erdp(&der);

        NdTheta epdt(&dep);
        NdTheta etdt(&det);
        NdTheta erdt(&der);



        double ephierr = 0.0;
        double ethetaerr = 0.0;
        double ererr = 0.0;
        double ephidphierr = 0.0;
        double ephidthetaerr = 0.0;
        double erdphierr = 0.0;
        double ethetadthetaerr = 0.0;
        double ethetadphierr = 0.0;
        double erdthetaerr = 0.0;


        for (int p = 0; p < NUMTRAPNODES; p++)
            for (int i = 0; i < NUMTRAPNODES; i++)
            {
                SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
                SphereCoord x(1, s);

                double diff = ePhi(x) - dep(x);
                ephierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta(x) - det(x);
                ethetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR(x) - der(x);
                ethetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dPhi(x) - epdp(x);
                ephidphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;



                diff = eTheta_dTheta(x) - etdt(x);
                ethetadthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dTheta(x) - erdt(x);
                erdthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dPhi(x) - erdp(x);
                erdphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                /*if (dot(diff, diff) > 1e-12)
                {
                    std::cout << " error of " << dot(diff, diff) << " at " << x << std::endl;
                    std::cout << "value of f: " << (*f1)(x) << std::endl;
                    std::cout << "value of g: " << (*f2)(x) << std::endl;

                }
                */
            }

        std::cout << "Error in computing ephi: " << ephierr << "\n";
        std::cout << "Error in computing etheta: " << ethetaerr << "\n";
        std::cout << "Error in computing er: " << ererr << "\n";

        std::cout << "Error in computing ephidphi: " << ephidphierr << "\n";
        std::cout << "Error in computing ephidtheta: " << ephidthetaerr << "\n";

        std::cout << "Error in computing ethetadtheta: " << ethetadthetaerr << "\n";
        std::cout << "Error in computing ethetadphi: " << ethetadphierr << "\n";

        std::cout << "Error in computing erdtheta: " << erdthetaerr << "\n";
        std::cout << "Error in computing erdphi: " << erdphierr << "\n";

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

    RectCoord operator()(SphereCoord x)
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
        return n * yi(x.s);
    }

    double eTheta(SphereCoord x)
    {
        return gi.eTheta(x.s);
    }

    double ePhi(SphereCoord x)
    {
        return gi.ePhi(x);
    }

    double eTheta_dTheta(SphereCoord x)
    {
        return gi.eTheta_dTheta(x);
    }

    double eTheta_dPhi(SphereCoord x)
    {
        return gi.eTheta_dPhi(x);
    }

    double ePhi_dPhi(SphereCoord x)
    {
        return gi.ePhi_dPhi(x);
    }

    double ePhi_dTheta(SphereCoord x)
    {
        return gi.ePhi_dTheta(x);
    }
    double eR_dTheta(SphereCoord x)
    {
        return (double)n * yi.dtheta(x.s);
    }

    double eR_dPhi(SphereCoord x)
    {
        return (double)n * yi.dphi(x.s);
    }


    void testhelper()
    {
        dot_ePhi dep(this);
        dot_eTheta det(this);
        dot_eR der(this);

        NdPhi epdp(&dep);
        NdPhi etdp(&det);
        NdPhi erdp(&der);

        NdTheta epdt(&dep);
        NdTheta etdt(&det);
        NdTheta erdt(&der);



        double ephierr = 0.0;
        double ethetaerr = 0.0;
        double ererr = 0.0;
        double ephidphierr = 0.0;
        double ephidthetaerr = 0.0;
        double erdphierr = 0.0;
        double ethetadthetaerr = 0.0;
        double ethetadphierr = 0.0;
        double erdthetaerr = 0.0;


        for (int p = 0; p < NUMTRAPNODES; p++)
            for (int i = 0; i < NUMTRAPNODES; i++)
            {
                SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
                SphereCoord x(1, s);

                double diff = ePhi(x) - dep(x);
                ephierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta(x) - det(x);
                ethetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR(x) - der(x);
                ethetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dPhi(x) - epdp(x);
                ephidphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;



                diff = eTheta_dTheta(x) - etdt(x);
                ethetadthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dTheta(x) - erdt(x);
                erdthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dPhi(x) - erdp(x);
                erdphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                /*if (dot(diff, diff) > 1e-12)
                {
                    std::cout << " error of " << dot(diff, diff) << " at " << x << std::endl;
                    std::cout << "value of f: " << (*f1)(x) << std::endl;
                    std::cout << "value of g: " << (*f2)(x) << std::endl;

                }
                */
            }

        std::cout << "Error in computing ephi: " << ephierr << "\n";
        std::cout << "Error in computing etheta: " << ethetaerr << "\n";
        std::cout << "Error in computing er: " << ererr << "\n";

        std::cout << "Error in computing ephidphi: " << ephidphierr << "\n";
        std::cout << "Error in computing ephidtheta: " << ephidthetaerr << "\n";

        std::cout << "Error in computing ethetadtheta: " << ethetadthetaerr << "\n";
        std::cout << "Error in computing ethetadphi: " << ethetadphierr << "\n";

        std::cout << "Error in computing erdtheta: " << erdthetaerr << "\n";
        std::cout << "Error in computing erdphi: " << erdphierr << "\n";

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

    RectCoord operator()(SphereCoord x)
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
        return 0;
    }

    double eTheta(SphereCoord x)
    {
        return -gr.ePhi(x);
    }

    double ePhi(SphereCoord x)
    {
        return gr.eTheta(x);
    }

    double eTheta_dTheta(SphereCoord x)
    {
        return -gr.ePhi_dTheta(x);
    }

    double eTheta_dPhi(SphereCoord x)
    {
        return -gr.ePhi_dPhi(x);
    }

    double ePhi_dPhi(SphereCoord x)
    {
        return gr.eTheta_dPhi(x);
    }

    double ePhi_dTheta(SphereCoord x)
    {
        return gr.eTheta_dTheta(x);
    }

    double eR_dTheta(SphereCoord x)
    {
        return 0;
    }

    double eR_dPhi(SphereCoord x)
    {
        return 0;
    }

    void testhelper()
    {
        dot_ePhi dep(this);
        dot_eTheta det(this);
        dot_eR der(this);

        NdPhi epdp(&dep);
        NdPhi etdp(&det);
        NdPhi erdp(&der);

        NdTheta epdt(&dep);
        NdTheta etdt(&det);
        NdTheta erdt(&der);



        double ephierr = 0.0;
        double ethetaerr = 0.0;
        double ererr = 0.0;
        double ephidphierr = 0.0;
        double ephidthetaerr = 0.0;
        double erdphierr = 0.0;
        double ethetadthetaerr = 0.0;
        double ethetadphierr = 0.0;
        double erdthetaerr = 0.0;


        for (int p = 0; p < NUMTRAPNODES; p++)
            for (int i = 0; i < NUMTRAPNODES; i++)
            {
                SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
                SphereCoord x(1, s);

                double diff = ePhi(x) - dep(x);
                ephierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta(x) - det(x);
                ethetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR(x) - der(x);
                ethetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dPhi(x) - epdp(x);
                ephidphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;



                diff = eTheta_dTheta(x) - etdt(x);
                ethetadthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dTheta(x) - erdt(x);
                erdthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dPhi(x) - erdp(x);
                erdphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                /*if (dot(diff, diff) > 1e-12)
                {
                    std::cout << " error of " << dot(diff, diff) << " at " << x << std::endl;
                    std::cout << "value of f: " << (*f1)(x) << std::endl;
                    std::cout << "value of g: " << (*f2)(x) << std::endl;

                }
                */
            }

        std::cout << "Error in computing ephi: " << ephierr << "\n";
        std::cout << "Error in computing etheta: " << ethetaerr << "\n";
        std::cout << "Error in computing er: " << ererr << "\n";

        std::cout << "Error in computing ephidphi: " << ephidphierr << "\n";
        std::cout << "Error in computing ephidtheta: " << ephidthetaerr << "\n";

        std::cout << "Error in computing ethetadtheta: " << ethetadthetaerr << "\n";
        std::cout << "Error in computing ethetadphi: " << ethetadphierr << "\n";

        std::cout << "Error in computing erdtheta: " << erdthetaerr << "\n";
        std::cout << "Error in computing erdphi: " << erdphierr << "\n";

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

    RectCoord operator()(SphereCoord x)
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
        return 0;
    }

    double eTheta(SphereCoord x)
    {
        return -gi.ePhi(x);
    }

    double ePhi(SphereCoord x)
    {
        return gi.eTheta(x);
    }

    double eTheta_dTheta(SphereCoord x)
    {
        return -gi.ePhi_dTheta(x);
    }

    double eTheta_dPhi(SphereCoord x)
    {
        return -gi.ePhi_dPhi(x);
    }

    double ePhi_dPhi(SphereCoord x)
    {
        return gi.eTheta_dPhi(x);
    }

    double ePhi_dTheta(SphereCoord x)
    {
        return gi.eTheta_dTheta(x);
    }

    double eR_dTheta(SphereCoord x)
    {
        return 0;
    }

    double eR_dPhi(SphereCoord x)
    {
        return 0;
    }

    void testhelper()
    {
        dot_ePhi dep(this);
        dot_eTheta det(this);
        dot_eR der(this);

        NdPhi epdp(&dep);
        NdPhi etdp(&det);
        NdPhi erdp(&der);

        NdTheta epdt(&dep);
        NdTheta etdt(&det);
        NdTheta erdt(&der);



        double ephierr = 0.0;
        double ethetaerr = 0.0;
        double ererr = 0.0;
        double ephidphierr = 0.0;
        double ephidthetaerr = 0.0;
        double erdphierr = 0.0;
        double ethetadthetaerr = 0.0;
        double ethetadphierr = 0.0;
        double erdthetaerr = 0.0;


        for (int p = 0; p < NUMTRAPNODES; p++)
            for (int i = 0; i < NUMTRAPNODES; i++)
            {
                SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
                SphereCoord x(1, s);

                double diff = ePhi(x) - dep(x);
                ephierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta(x) - det(x);
                ethetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR(x) - der(x);
                ethetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dPhi(x) - epdp(x);
                ephidphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;



                diff = eTheta_dTheta(x) - etdt(x);
                ethetadthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dTheta(x) - erdt(x);
                erdthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dPhi(x) - erdp(x);
                erdphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                /*if (dot(diff, diff) > 1e-12)
                {
                    std::cout << " error of " << dot(diff, diff) << " at " << x << std::endl;
                    std::cout << "value of f: " << (*f1)(x) << std::endl;
                    std::cout << "value of g: " << (*f2)(x) << std::endl;

                }
                */
            }

        std::cout << "Error in computing ephi: " << ephierr << "\n";
        std::cout << "Error in computing etheta: " << ethetaerr << "\n";
        std::cout << "Error in computing er: " << ererr << "\n";

        std::cout << "Error in computing ephidphi: " << ephidphierr << "\n";
        std::cout << "Error in computing ephidtheta: " << ephidthetaerr << "\n";

        std::cout << "Error in computing ethetadtheta: " << ethetadthetaerr << "\n";
        std::cout << "Error in computing ethetadphi: " << ethetadphierr << "\n";

        std::cout << "Error in computing erdtheta: " << erdthetaerr << "\n";
        std::cout << "Error in computing erdphi: " << erdphierr << "\n";

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

    SphericalHarmonic operator +(SphericalHarmonic t)
    {
        std::array<double, 6> temp;

        temp[0] = rhohatVr + t.rhohatVr;
        temp[1] = rhohatVi + t.rhohatVi;
        temp[2] = rhohatWr + t.rhohatWr;
        temp[3] = rhohatWi + t.rhohatWi;
        temp[4] = rhohatXr + t.rhohatXr;
        temp[5] = rhohatXi + t.rhohatXi;

        return SphericalHarmonic(m, n, temp);
    }

    SphericalHarmonic operator -(SphericalHarmonic t)
    {
        std::array<double, 6> temp;

        temp[0] = rhohatVr - t.rhohatVr;
        temp[1] = rhohatVi - t.rhohatVi;
        temp[2] = rhohatWr - t.rhohatWr;
        temp[3] = rhohatWi - t.rhohatWi;
        temp[4] = rhohatXr - t.rhohatXr;
        temp[5] = rhohatXi - t.rhohatXi;

        return SphericalHarmonic(m, n, temp);
    }

    double dot(SphericalHarmonic s)
    {
        if (m != s.m || n != s.n)
            return 0.0;
        else
            return rhohatVr * s.rhohatVr + rhohatVi * s.rhohatVi + rhohatWr * s.rhohatWr + rhohatWi * s.rhohatWi + rhohatXr * s.rhohatXr + rhohatXi * s.rhohatXi;
    }


    RectCoord operator()(SphereCoord x)
    {
        if (n == 0)
            return rhohatVr * vr(x) + rhohatVi * vi(x);

        RectCoord temp = rhohatVr * vr(x) + rhohatVi * vi(x);
        temp = temp + (rhohatWr * wr(x) + rhohatWi * wi(x));
        return temp + (rhohatXr * xr(x) + rhohatXi * xi(x));
    }

};

class VSHSeries : public VectorFieldSum
{

private:



    Legendre P;



public:
    int N;
    int NUMGRIDS;
    RectCoord center;
    std::vector<std::vector<SphericalHarmonic>> terms;




    VSHSeries(int n, int numgrids , RectCoord c = RectCoord())
    {
        N = n;
        NUMGRIDS = numgrids;

        terms.resize(n + 1);

        for (int i = 0; i <= n; i++)
            terms[i].resize(i + 1);

        center = c;



    }

    void approximate(SphericalVectorField* r, int n = 0, int ng = 0)
    {

        if (n > 0)
        {
            terms.resize(n + 1);

            for (int i = 0; i <= n; i++)
                terms[i].resize( i + 1);

        }

        if (ng > 0)
            NUMGRIDS = ng;



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
                    rhohats[0] = L2InnerProduct(r, &vr, NUMGRIDS, center) / (double)((2 * n + 1) * (n + 1));

                    vi.reset(m, n);
                    rhohats[1] = L2InnerProduct(r, &vi, NUMGRIDS, center) / (double)((2 * n + 1) * (n + 1));

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

            }
    }

    void approximate(SphereData & data, int n = 0, int ng = 0)
    {

        if (n > 0)
        {
            terms.resize(n + 1);

            for (int i = 0; i <= n; i++)
                terms[i].resize( i + 1);

        }

        if (ng > 0)
            NUMGRIDS = ng;



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
                    rhohats[0] = L2InnerProductDiscrete(data, &vr, NUMGRIDS, center) / (double)((2 * n + 1) * (n + 1));

                    vi.reset(m, n);
                    rhohats[1] = L2InnerProductDiscrete(data, &vi, NUMGRIDS, center) / (double)((2 * n + 1) * (n + 1));

                    wr.reset(m, n);
                    rhohats[2] = L2InnerProductDiscrete(data, &wr, NUMGRIDS) / (double)((2 * n + 1) * n);

                    wi.reset(m, n);
                    rhohats[3] = L2InnerProductDiscrete(data, &wi, NUMGRIDS) / (double)((2 * n + 1) * n);

                    xr.reset(m, n);
                    rhohats[4] = L2InnerProductDiscrete(data, &xr, NUMGRIDS) / (double)((n + 1) * n);

                    xi.reset(m, n);
                    rhohats[5] = L2InnerProductDiscrete(data, &xi, NUMGRIDS) / (double)((n + 1) * n);

                    terms[n][m] = SphericalHarmonic(m, n, rhohats);
                }
                else if (m > 0)
                {
                    vr.reset(m, n);
                    rhohats[0] = 2.0 * L2InnerProductDiscrete(data, &vr, NUMGRIDS) / (double)((2 * n + 1) * (n + 1));

                    vi.reset(m, n);
                    rhohats[1] = 2.0 * L2InnerProductDiscrete(data, &vi, NUMGRIDS) / (double)((2 * n + 1) * (n + 1));

                    wr.reset(m, n);
                    rhohats[2] = 2.0 * L2InnerProductDiscrete(data, &wr, NUMGRIDS) / (double)((2 * n + 1) * n);

                    wi.reset(m, n);
                    rhohats[3] = 2.0 * L2InnerProductDiscrete(data, &wi, NUMGRIDS) / (double)((2 * n + 1) * n);

                    xr.reset(m, n);
                    rhohats[4] = 2.0 * L2InnerProductDiscrete(data, &xr, NUMGRIDS) / (double)((n + 1) * n);

                    xi.reset(m, n);
                    rhohats[5] = 2.0 * L2InnerProductDiscrete(data, &xi, NUMGRIDS) / (double)((n + 1) * n);

                    terms[n][m] = SphericalHarmonic(m, n, rhohats);
                }

            }
    }

    VSHSeries operator +(const VSHSeries & s)
    {
        VSHSeries temp(N, NUMGRIDS, center);

        temp.terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            temp.terms[i].resize(i + 1);



        for(int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                temp.terms[n][m] = terms[n][m] + s.terms[n][m];

        return temp;
            
    }

    VSHSeries operator -(VSHSeries s)
    {
        VSHSeries temp(N, NUMGRIDS, center);

        temp.terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            temp.terms[i].resize(i + 1);



        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                temp.terms[n][m] = terms[n][m] - s.terms[n][m];

        return temp;

    }

    double dot(VSHSeries & series)
    {
        double total = 0.0;

        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                total += terms[n][m].dot(series.terms[n][m]);

        return total;
    }

    RectCoord operator()(SphereCoord x)
    {

        SphereCoord temp = recenter(x, center);

        P.populate(cos(temp.s.theta), N, N);

        return VectorFieldSum::operator()(temp);
    }
};