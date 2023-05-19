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

/**
* \file "SphericalHarmonics.h"
* \brief defninition of the Spherical Harmonic basis functions, and series representation.
* 
* Here we intruduce the Spherical Harmonic basis functions, The eigenfunctions for the Laplacian on a sphere.
* The complex valued scalar spherical harmonics are given in spherical coordinates by 
* @f{equation}{
   Y_n^m(\theta , \phi) = \alpha_{n,m} P_n^m(\cos(\theta)) e^{i m \phi}. \quad |m| \leq n.
 @f}
* Throughout this work, we do not deal in complex numbers, and define the real and imaginary parts of @f$ Y_n^m @f$ separately.
* From the above, we define the Vector Spherical Harmonics to be 
* @f{eqnarray}{
\boldsymbol{V}_n^m &=& -(n+1) Y_n^m \boldsymbol{e}_r + \nabla_\Gamma Y_n^m \\
\boldsymbol{W}_n^m &=& n Y_n^m 
\boldsymbol{e}_r + \nabla_\Gamma Y_n^m \\
\boldsymbol{X}_n^m &=& \boldsymbol{e}_r \times \nabla_\Gamma Y_n^m.
 @f}
 Given a density @f$ \boldsymbol{\rho}: \Gamma \rightarrow R^3 @f$ periodic in the azimuth angle, we can write it in the form  
 @f{equation}{
   \boldsymbol{\rho}(\theta , \phi) =\sum_{n = 0}^\infty\sum_{-n \leq m \leq n} \frac{\hat{\rho}_{n,m}^V}{\|\boldsymbol{V}_n^m\|_2^2}\boldsymbol{V}_n^m(\theta, \phi) +  \frac{\hat{\rho}_{n,m}^W}{\|\boldsymbol{W}_n^m\|_2^2}\boldsymbol{W}_n^m(\theta, \phi)
    +\frac{\hat{\rho}_{n,m}^X}{\|\boldsymbol{X}_n^m\|_2^2}\boldsymbol{X}_n^m(\theta, \phi).
 @f}
 With the @f$ \hat{\rho} @f$'s are found by forming the inner product
  @f{equation}{
    \hat{\rho}_{n,m}^U = \int_{\Gamma} \boldsymbol{\rho}(\theta , \phi) \cdot \bar{\boldsymbol{U}}_n^m(\theta , \phi) dS_x
 @f}
*
*/



/**
* @brief minimum macro.
*/
#define min(a , b) a < b? a : b

/**
* @brief maximum macro.
*/
#define max(a , b) a > b? a : b


Legendre legendrePoly;


/**
 * @brief The Real part of @f$ Y_n^m @f$.
 * 
 * Computes the real part of the scalar spherical harmonic and derivatives up to second order.
*/
class YReal : public SurfaceScalarFunction
{
private:
    int m;
    int n;
    Legendre* poly;

public:
    /**
     * @brief Constructor, sets the values of m and n.
     * @param a the order(m).
     * @param b the degree (n).
     * @note we require @f$ 0 \leq m \leq n @f$.
    */
    YReal(int a = 0, int b = 0) : poly(&legendrePoly)
    {
        m = a;
        n = b;
        name = "YReal";

        if (m < 0)
            std::cout << "Index error in YReal, m < 0" << std::endl;

        if(m > n)
            std::cout << "Index error in YReal, m > n" << std::endl;
    }

    YReal(const YReal& yr) : poly(&legendrePoly)
    {
        m = yr.m;
        n = yr.n;

        name = "YReal";
    }

    const YReal & operator =(const YReal& yr)
    {
        m = yr.m;
        n = yr.n;

        name = "YReal";

        return *this;
    }

    void operator ~()
    {
    }
    /**
     * @brief Evaluates @f$ \mathfrak{R}(Y_n^m) @f$ at s.
     * @return @f{equation}{
      \alpha_{n,m} P_n^m(\cos(\theta)) \cos{m \phi}
       @f}
    */
    double operator()(const SurfaceCoord & s) const 
    {



        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        coef *= (*poly)(cos(s.theta), m, n);

        //std::cout << coef << std::endl;

        return coef * cos((double)m * s.phi);


    }

    /**
     * @brief Computes the partial derivative of @f$ Y_n^m @f$ with respect to @f$ \phi @f$.
     * @return  @f{equation}{
      \frac{\partial \mathfrak{R}(Y_n^m)}{\partial \phi} = -m \alpha_{n,m} P_n^m(\cos(\theta)) \sin{m \phi}
       @f}
    */
    double dPhi(const SurfaceCoord & s) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        coef *= (*poly)(cos(s.theta), m, n);

        //std::cout << coef << std::endl;

        return - coef *m * sin((double)m * s.phi);
    }

    /**
     * @brief Computes the second order partial derivative of @f$ Y_n^m @f$ with respect to @f$ \phi @f$ twice.
     * @return  @f{equation}{
      \frac{\partial^2 \mathfrak{R}(Y_n^m)}{\partial \phi^2} = -m^2 \mathfrak{R}(Y_n^m)
       @f}
    */
    double dPhidPhi(const SurfaceCoord & s) const
    {
        double m_(m);
        return -m_ * m_ * operator()(s);
    }

    /**
     * @brief Computes the partial derivative of @f$ Y_n^m @f$ with respect to @f$ \theta @f$.
     * @return @f{equation}{
      \alpha_{n,m} \frac{\partial P_n^m(\cos(\theta))}{\partial \theta} \cos{m \phi}
       @f}
       See also: @ref  Legendre::dTheta()
    */
    double dTheta(const SurfaceCoord & s) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly->populate(cos(s.theta), m + 1, n + 1);

        coef *= poly->dTheta(m, n, s.theta);


        return coef * cos((double)m * s.phi);
    }

    /**
     * @brief Computes the second partial derivative of @f$ Y_n^m @f$ with respect to @f$ \theta @f$ twice.
     * @return @f{equation}{
      \alpha_{n,m} \frac{\partial^2 P_n^m(\cos(\theta))}{\partial \theta^2} \cos{m \phi}
       @f}
       See also: @ref Legendre::dTheta()
    */
    double dThetadTheta(const SurfaceCoord & s) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly->populate(cos(s.theta), m + 2, n + 2);

        coef *= poly->dTheta(m, n, s.theta , 2);


        return coef * cos((double)m * s.phi);
    }

    /**
     * @brief Computes the second mixed partial derivative of @f$ Y_n^m @f$ with respect to @f$ \theta @f$ and @f$ \phi @f$.
     * @return @f{equation}{
      -m \alpha_{n,m} \frac{\partial P_n^m(\cos(\theta))}{\partial \theta} \sin{m \phi}
       @f}
       See also: @ref Legendre::dTheta() @ref dPhi()
    */
    double dThetadPhi(const SurfaceCoord & s) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly->populate(cos(s.theta), m + 1, n + 1);

        coef *= poly->dTheta(m, n, s.theta);


        return -   coef * m* sin((double)m * s.phi);
    }

    

};





/**
 * @brief The Imaginary part of @f$ Y_n^m @f$.
 *
 * Computes the imaginary part of the scalar spherical harmonic and derivatives up to second order.
*/
class YImag : public SurfaceScalarFunction
{
private:
    int m;
    int n;
    Legendre* poly;

public:

    /**
     * @brief Constructor, sets the values of m and n.
     * @param a the order(m).
     * @param b the degree (n).
     * @note we require @f$ 0 \leq m \leq n @f$.
    */
    YImag(int a = 0, int b = 0) : poly(&legendrePoly)
    {
        m = a;
        n = b;
        name = "YImag";

        if (m < 0)
            std::cout << "Index error in YImag, m < 0" << std::endl;

        if (m > n)
            std::cout << "Index error in YImag, m > n" << std::endl;

    }

    YImag(const YImag& yi) : poly(&legendrePoly)
    {
        m = yi.m;
        n = yi.n;

        name = "YReal";
    }

    const YImag& operator =(const YImag& yi)
    {
        m = yi.m;
        n = yi.n;

        name = "YReal";

        return *this;
    }

    void operator ~()
    {
    }

    /**
     * @brief Evaluates @f$ \mathfrak{I}(Y_n^m) @f$ at s.
     * @return @f{equation}{
      \alpha_{n,m} P_n^m(\cos(\theta)) \sin{m \phi}
       @f}
    */
    double operator()(const SurfaceCoord & s) const 
    {


        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        coef *= (*poly)(cos(s.theta), m, n);

        //std::cout << coef << std::endl;

        return coef * sin((double)m * s.phi);


    }

    /**
     * @brief Computes the partial derivative of @f$ Y_n^m @f$ with respect to @f$ \phi @f$.
     * @return  @f{equation}{
      \frac{\partial \mathfrak{R}(Y_n^m)}{\partial \phi} = m \alpha_{n,m} P_n^m(\cos(\theta)) \cos{m \phi}
       @f}
    */
    double dPhi(const SurfaceCoord & s) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        coef *= (*poly)(cos(s.theta), m, n);

        //std::cout << coef << std::endl;

        return coef * m * cos((double)m * s.phi);
    }

    /**
     * @brief Computes the second order partial derivative of @f$ Y_n^m @f$ with respect to @f$ \phi @f$ twice.
     * @return  @f{equation}{
      \frac{\partial^2 \mathfrak{I}(Y_n^m)}{\partial \phi^2} = -m^2 \mathfrak{I}(Y_n^m)
       @f}
    */
    double dPhidPhi(const SurfaceCoord & s) const
    {
        double m_(m);
        return -m_ * m_ * operator()(s);
    }

    /**
     * @brief Computes the partial derivative of @f$ Y_n^m @f$ with respect to @f$ \theta @f$.
     * @return @f{equation}{
      \alpha_{n,m} \frac{\partial P_n^m(\cos(\theta))}{\partial \theta} \cos{m \phi}
       @f}
       See also: @ref  Legendre::dTheta()
    */
    double dTheta(const SurfaceCoord & s) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly->populate(cos(s.theta), m + 1, n + 1);

        coef *= poly->dTheta(m, n, s.theta);


        return coef * sin((double)m * s.phi);
    }

    /**
     * @brief Computes the second partial derivative of @f$ Y_n^m @f$ with respect to @f$ \theta @f$ twice.
     * @return @f{equation}{
      \alpha_{n,m} \frac{\partial^2 P_n^m(\cos(\theta))}{\partial \theta^2} \cos{m \phi}
       @f}
       See also: @ref Legendre::dTheta()
    */
    double dThetadTheta(const SurfaceCoord & s) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly->populate(cos(s.theta), m + 1, n + 1);

        coef *= poly->dTheta(m, n, s.theta , 2);


        return coef * sin((double)m * s.phi);
    }

    /**
     * @brief Computes the second mixed partial derivative of @f$ Y_n^m @f$ with respect to @f$ \theta @f$ and @f$ \phi @f$.
     * @return @f{equation}{
      m \alpha_{n,m} \frac{\partial P_n^m(\cos(\theta))}{\partial \theta} \cos{m \phi}
       @f}
       See also: @ref Legendre::dTheta() @ref dPhi()
    */
    double dThetadPhi(const SurfaceCoord & s) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly->populate(cos(s.theta), m + 1, n + 1);

        coef *= poly->dTheta(m, n, s.theta);


        return coef * m * cos((double)m * s.phi);
    }

};

/**
 * @brief Class testing the derivatives of the real part of @f$ Y_n^m @f$. 
*/
class Yrtest : public SphericalScalarFunction
{
private:
    YReal yr;
public:
    /**
     * @brief Evaluates to YReal.
    */
    double operator()(const SphereCoord & x) const 
    {
        return yr(x.s);
    }

    /**
     * @brief Compares the values of the derivatives given by member function to those given by generic numerical methods.
     *The numerical methods are given by the NdPhi, NdTheta classes, possibly with repeated applications.
    */
    void testhelper(int m ,int n)
    {
            
            yr = YReal(m, n);
            NdPhi dPhi(this);
            NdPhi dPhidPhi(&dPhi);
            NdTheta dTheta(this);
            NdTheta dThetadTheta(&dTheta);
            NdTheta dPhidTheta(&dPhi);


            double dPhierr = 0.0;
            double dPhidPhierr = 0.0;
            double dThetaerr = 0.0;
            double dThetadThetaerr = 0.0;
            double dPhidThetaerr = 0.0;

            for (int p = 0; p < NUMTRAPNODES; p++)
                for (int i = 0; i < NUMGLNODES; i++)
                {
                    SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
                    SphereCoord x(1, s);
                    double diff = dPhi(x) - yr.dPhi(x.s);
                    dPhierr +=sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                    diff = dPhidPhi(x) - yr.dPhidPhi(x.s);
                    dPhidPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                    diff = dTheta(x) - yr.dTheta(x.s);
                    dThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                    diff = dThetadTheta(x) - yr.dThetadTheta(x.s);
                    dThetadThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                    diff = dPhidTheta(x) - yr.dThetadPhi(x.s);
                    dPhidThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                    /*if (dot(diff, diff) > 1e-12)
                    {
                        std::cout << " error of " << dot(diff, diff) << " at " << x << std::endl;
                        std::cout << "value of f: " << (*f1)(x) << std::endl;
                        std::cout << "value of g: " << (*f2)(x) << std::endl;

                    }
                    */
                }

            std::cout << "Error in computing dPhi: " << dPhierr << "\n";
            std::cout << "Error in computing dPhidPhi: " << dPhidPhierr << "\n";
            std::cout << "Error in computing dTheta: " << dThetaerr << "\n";
            std::cout << "Error in computing dThetadTheta: " << dThetadThetaerr << "\n";
            std::cout << " Error in computing dThetadPhi: " << dPhidThetaerr << "\n";
        }
    

};

/**
 * @brief Class testing the derivatives of the imaginary part of @f$ Y_n^m @f$.
*/
class Yitest : public SphericalScalarFunction
{
private:
    YImag yi;
public:

    /**
     * @brief Evaluates to YImag.
    */
    double  operator()(const SphereCoord & x) const
    {
        return yi(x.s);
    }

    /**
     * @brief Compares the values of the derivatives given by member function to those given by generic numerical methods.
     *The numerical methods are given by the NdPhi, NdTheta classes, possibly with repeated applications.
    */
    void testhelper(int m, int n)
    {
       
        yi = YImag(m, n);
        NdPhi dPhi(this);
        NdPhi dPhidPhi(&dPhi);
        NdTheta dTheta(this);
        NdTheta dThetadTheta(&dTheta);
        NdTheta dPhidTheta(&dPhi);


        double dPhierr = 0.0;
        double dPhidPhierr = 0.0;
        double dThetaerr = 0.0;
        double dThetadThetaerr = 0.0;
        double dPhidThetaerr = 0.0;

        for (int p = 0; p < NUMTRAPNODES; p++)
            for (int i = 0; i < NUMTRAPNODES; i++)
            {
                SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
                SphereCoord x(1, s);
                double diff = dPhi(x) - yi.dPhi(x.s);
                dPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = dPhidPhi(x) - yi.dPhidPhi(x.s);
                dPhidPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = dTheta(x) - yi.dTheta(x.s);
                dThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = dThetadTheta(x) - yi.dThetadTheta(x.s);
                dThetadThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = dPhidTheta(x) - yi.dThetadPhi(x.s);
                dPhidThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                /*if (dot(diff, diff) > 1e-12)
                {
                    std::cout << " error of " << dot(diff, diff) << " at " << x << std::endl;
                    std::cout << "value of f: " << (*f1)(x) << std::endl;
                    std::cout << "value of g: " << (*f2)(x) << std::endl;

                }
                */
            }

        std::cout << "Error in computing dPhi: " << dPhierr << "\n";
        std::cout << "Error in computing dPhidPhi: " << dPhidPhierr << "\n";
        std::cout << "Error in computing dTheta: " << dThetaerr << "\n";
        std::cout << "Error in computing dThetadTheta: " << dThetadThetaerr << "\n";
        std::cout << " Error in computing dThetadPhi: " << dPhidThetaerr << "\n";
    }


};
/**
 * @brief Functor for surface gradient of @f$ \mathfrak{R}(Y_n^m) @f$. Includes derivatives and spherical basis component extractions. 
*/
class GReal : public SphericalVectorField
{
private:
    int m;
    int n;
    YReal yr;

public:
    GReal() {}
    /**
     * @brief Assigns the order and degree to the object.
    */
    GReal(int m0, int n0)
    {
        m = m0;
        n = n0;

        yr = YReal(m, n);
    }

    /**
     * @brief The surface gradient of @f$ \mathfrak{R}(Y_n^m) @f$.
     * @return @f{equation}{
     * \nabla_\Gamma \mathfrak{R}(Y_n^m) = \frac{\partial \mathfrak{R}(Y_n^m)}{\partial \theta} \boldsymbol{e}_\theta + \frac{1}{\sin\theta} \frac{\partial \mathfrak{R}(Y_n^m)}{\partial\phi}\boldsymbol{e}_\phi
     * @f} evaluated at x. The partial derivatives are YReal::dPhi() and YReal::dTheta().
    */
    RectCoord  operator()(const SphereCoord & x) const
    {
        double dTheta = yr.dTheta(x.s);

        double dPhi = yr.dPhi(x.s);

        return dTheta * e_theta(x) + dPhi / sin(x.s.theta) * e_phi(x);


    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \nabla_\Gamma \mathfrak{R}(Y_n^m) @f$.
     * @return @f$ \frac{\partial \mathfrak{R}(Y_n^m)}{\partial \theta} @f$ evaluated at x.
     * 
     * See also: YReal::dTheta().
    */
    double eTheta(const SphereCoord & x) const
    {
        return yr.dTheta(x.s);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \nabla_\Gamma \mathfrak{R}(Y_n^m) @f$.
     * @return @f$ \frac{1}{\sin\theta} \frac{\partial \mathfrak{R}(Y_n^m)}{\partial \phi} @f$ evaluated at x.
     *
     * See also: YReal::dPhi().
    */
    double ePhi(const SphereCoord & x) const
    {
        return yr.dPhi(x.s) / sin(x.s.theta);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \nabla_\Gamma \mathfrak{R}(Y_n^m)@f$, differentiated with respect to @f$\theta @f$.
     * @return @f$ \frac{\partial^2 \mathfrak{R}(Y_n^m)}{\partial \theta^2} @f$ evaluated at x.
     *
     * See also: YReal::dThetadTheta().
    */
    double eTheta_dTheta(const SphereCoord & x) const
    {
        return yr.dThetadTheta(x.s);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \nabla_\Gamma \mathfrak{R}(Y_n^m)@f$, differentiated with respect to @f$\phi @f$.
     * @return @f$ \frac{\partial^2 \mathfrak{R}(Y_n^m)}{\partial \theta \partial \phi} @f$ evaluated at x.
     *
     * See also: YReal::dThetadPhi().
    */
    double eTheta_dPhi(const SphereCoord & x) const
    {
        return yr.dThetadPhi(x.s);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \nabla_\Gamma \mathfrak{R}(Y_n^m)@f$, differentiated with respect to @f$\phi @f$.
     * @return @f$ \frac{1}{\sin \theta}\frac{\partial^2 \mathfrak{R}(Y_n^m)}{\partial \phi^2} @f$ evaluated at x.
     *
     * See also: YReal::dPhidPhi().
    */
    double ePhi_dPhi(const SphereCoord & x) const
    {
        return yr.dPhidPhi(x.s) / sin(x.s.theta);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \nabla_\Gamma \mathfrak{R}(Y_n^m)@f$, differentiated with respect to @f$\theta @f$.
     * @return @f{equation}{ 
     \frac{\sin \theta \ \frac{\partial^2 \mathfrak{R}(Y_n^m)}{\partial \theta \partial \phi} 
            - \cos \theta \frac{\partial \mathfrak{R}(Y_n^m)}{\partial \phi} }{\sin^2 \theta}
     @f} evaluated at x.
     *
     * See also: YReal::dThetadPhi() , YReal::dPhi().
    */
    double ePhi_dTheta(const SphereCoord & x) const
    {
        return (yr.dThetadPhi(x.s) * sin(x.s.theta) - cos(x.s.theta) * yr.dPhi(x.s)) / (sin(x.s.theta) * sin(x.s.theta));
    }

    /**
     * @brief Compares the values of the derivatives given by member function to those given by generic numerical methods.
     *The numerical methods are given by the NdPhi, NdTheta classes, possibly with repeated applications.
    */
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
        double ephidPhierr = 0.0;
        double ephidThetaerr = 0.0;
        double ethetadThetaerr = 0.0;
        double ethetadPhierr = 0.0;

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
                ephidPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dTheta(x) - etdt(x);
                ethetadThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
               
            }

        std::cout << "Error in computing ephi: " << ephierr << "\n";
        std::cout << "Error in computing ephidPhi: " << ephidPhierr << "\n";
        std::cout << "Error in computing etheta: " << ethetaerr << "\n";
        std::cout << "Error in computing ethetadTheta: " << ethetadThetaerr << "\n";
        std::cout << " Error in computing ethetadPhi: " << ethetadPhierr << "\n";
        std::cout << " Error in computing ephidTheta: " << ephidThetaerr << "\n";
    }
};



/**
 * @brief Functor for surface gradient of @f$ \mathfrak{I}(Y_n^m) @f$. Includes derivatives and spherical basis component extractions.
*/
class GImag : public SphericalVectorField
{
private:
    int m;
    int n;
    YImag yi;

public:
    GImag() {}


    /**
     * @brief Assigns the order and degree to the object.
    */
    GImag(int m0, int n0)
    {
        m = m0;
        n = n0;

        yi = YImag(m, n);
    }

    
    /**
     * @brief The surface gradient of @f$ \mathfrak{I}(Y_n^m) @f$.
     * @return @f{equation}{
     * \nabla_\Gamma \mathfrak{I}(Y_n^m) = \frac{\partial \mathfrak{I}(Y_n^m)}{\partial \theta} \boldsymbol{e}_\theta + \frac{1}{\sin\theta} \frac{\partial \mathfrak{I}(Y_n^m)}{\partial\phi}\boldsymbol{e}_\phi
     * @f} evaluated at x. The partial derivatives are YImag::dPhi() and YImag::dTheta().
    */
    RectCoord  operator()(const SphereCoord & x) const
    {


        double dTheta = yi.dTheta(x.s);

        double dPhi = yi.dPhi(x.s);

        return dTheta * e_theta(x) + dPhi / sin(x.s.theta) * e_phi(x);


    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \nabla_\Gamma \mathfrak{I}(Y_n^m) @f$.
     * @return @f$ \frac{\partial \mathfrak{I}(Y_n^m)}{\partial \theta} @f$ evaluated at x.
     *
     * See also: YImag::dTheta().
    */
    double eTheta(const SphereCoord & x) const
    {
        return yi.dTheta(x.s);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \nabla_\Gamma \mathfrak{I}(Y_n^m) @f$.
     * @return @f$ \frac{1}{\sin\theta} \frac{\partial \mathfrak{I}(Y_n^m)}{\partial \phi} @f$ evaluated at x.
     *
     * See also: YImag::dPhi().
    */
    double ePhi(const SphereCoord & x) const
    {
        return yi.dPhi(x.s) / sin(x.s.theta);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \nabla_\Gamma \mathfrak{I}(Y_n^m)@f$, differentiated with respect to @f$\theta @f$.
     * @return @f$ \frac{\partial^2 \mathfrak{I}(Y_n^m)}{\partial \theta^2} @f$ evaluated at x.
     *
     * See also: YImag::dThetadTheta().
    */
    double eTheta_dTheta(const SphereCoord & x) const
    {
        return yi.dThetadTheta(x.s);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \nabla_\Gamma \mathfrak{I}(Y_n^m)@f$, differentiated with respect to @f$\phi @f$.
     * @return @f$ \frac{\partial^2 \mathfrak{I}(Y_n^m)}{\partial \theta \partial \phi} @f$ evaluated at x.
     *
     * See also: YImag::dThetadPhi().
    */
    double eTheta_dPhi(const SphereCoord & x) const
    {
        return yi.dThetadPhi(x.s);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \nabla_\Gamma \mathfrak{I}(Y_n^m)@f$, differentiated with respect to @f$\phi @f$.
     * @return @f$ \frac{1}{\sin \theta}\frac{\partial^2 \mathfrak{I}(Y_n^m)}{\partial \phi^2} @f$ evaluated at x.
     *
     * See also: YImag::dPhidPhi().
    */
    double ePhi_dPhi(const SphereCoord & x) const
    {
        return yi.dPhidPhi(x.s) / sin(x.s.theta);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \nabla_\Gamma \mathfrak{I}(Y_n^m)@f$, differentiated with respect to @f$\theta @f$.
     * @return @f{equation}{
     \frac{\sin \theta \ \frac{\partial^2 \mathfrak{I}(Y_n^m)}{\partial \theta \partial \phi}
            - \cos \theta \frac{\partial \mathfrak{I}(Y_n^m)}{\partial \phi} }{\sin^2 \theta}
     @f} evaluated at x.
     *
     * See also: YImag::dThetadPhi() , YImag::dPhi().
    */
    double ePhi_dTheta(const SphereCoord & x) const
    {
        return (yi.dThetadPhi(x.s) * sin(x.s.theta) - cos(x.s.theta) * yi.dPhi(x.s)) / (sin(x.s.theta) * sin(x.s.theta));
    }

    /**
     * @brief Compares the values of the derivatives given by member function to those given by generic numerical methods.
     *The numerical methods are given by the NdPhi, NdTheta classes, possibly with repeated applications.
    */
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
        double ephidPhierr = 0.0;
        double ephidThetaerr = 0.0;
        double ethetadThetaerr = 0.0;
        double ethetadPhierr = 0.0;

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
                ephidPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dTheta(x) - etdt(x);
                ethetadThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                
            }

        std::cout << "Error in computing ephi: " << ephierr << "\n";
        std::cout << "Error in computing ephidPhi: " << ephidPhierr << "\n";
        std::cout << "Error in computing etheta: " << ethetaerr << "\n";
        std::cout << "Error in computing ethetadTheta: " << ethetadThetaerr << "\n";
        std::cout << " Error in computing ethetadPhi: " << ethetadPhierr << "\n";
        std::cout << " Error in computing ephidTheta: " << ephidThetaerr << "\n";
    }

};

/**
 * @brief Functor for @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$. Also has functions for component extraction and differentiation.
*/
class VReal : public SphericalVectorField
{
    int m;
    int n;
    YReal yr;
    GReal gr;

public:


    /**
     * @brief Initializes the order and degree.
    */
    VReal(int a = 0, int b = 0)
    {
        m = a;
        n = b;
        yr = YReal(m, n);
        gr = GReal(m, n);
    }

    /**
     * @brief Evaluates  @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$.
     * @return @f$ -(n+1) \mathfrak{R}(\boldsymbol{Y}_n^m)(x) \boldsymbol{e}_rs + \nabla_\Gamma \mathfrak{R}(\boldsymbol{Y}_n^m)(x) @f$.
     * 
     * See YReal::operator()() , GReal::operator()().
    */
    RectCoord operator()(const SphereCoord & x) const 
    {
        return -((double)n + 1.0) * yr(x.s) * e_r(x) + gr(x);
    }

    /**
     * @brief resets the indices.
    */
    void reset(int a, int b)
    {
        m = a;
        n = b;
        yr = YReal(m, n);
        gr = GReal(m, n);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$.
     * @return  @f$ -(n+1) \mathfrak{R}(\boldsymbol{Y}_n^m)(x) @f$.
     * 
     * See: YReal.
    */
    double eR(const SphereCoord & x) const
    {
        return -((double)n + 1.0) * yr(x.s);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$.
     * @return  @f$ GReal::eTheta() @f$.
    */
    double eTheta(const SphereCoord & x) const
    {
        return gr.eTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$.
     * @return GReal::ePhi().
    */
    double ePhi(const SphereCoord & x) const
    {
        return gr.ePhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return GReal::eTheta_dTheta().
    */
    double eTheta_dTheta(const SphereCoord & x) const
    {
        return gr.eTheta_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GReal::eTheta_dPhi().
    */
    double eTheta_dPhi(const SphereCoord & x) const
    {
        return gr.eTheta_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GReal::ePhi_dPhi().
    */
    double ePhi_dPhi(const SphereCoord & x) const
    {
        return gr.ePhi_dPhi(x);
    }
    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return GReal::ePhi_dTheta().
    */
    double ePhi_dTheta(const SphereCoord & x) const
    {
        return gr.ePhi_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return  @f$ -(n+1)\frac{\partial \mathfrak{R}(\boldsymbol{Y}_n^m)}{\partial \theta}  @f$.
     * 
     * See: YReal::dTheta()
    */
    double eR_dTheta(const SphereCoord & x) const
    {
        return -((double) n+1.0)*yr.dTheta(x.s);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return  @f$ -(n+1)\frac{\partial \mathfrak{R}(\boldsymbol{Y}_n^m)}{\partial \phi}  @f$.
     *
     * See: YReal::dPhi()
    */
    double eR_dPhi(const SphereCoord & x) const
    {
        return -((double)n + 1.0) * yr.dPhi(x.s);
    }

    /**
     * @brief Compares the values of the derivatives given by member function to those given by generic numerical methods.
     *The numerical methods are given by the NdPhi, NdTheta classes, possibly with repeated applications.
    */
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
        double ephidPhierr = 0.0;
        double ephidThetaerr = 0.0;
        double erdPhierr = 0.0;
        double ethetadThetaerr = 0.0;
        double ethetadPhierr = 0.0;
        double erdThetaerr = 0.0;


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
                ephidPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                

                diff = eTheta_dTheta(x) - etdt(x);
                ethetadThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                
                diff = eR_dTheta(x) - erdt(x);
                erdThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dPhi(x) - erdp(x);
                erdPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
                
            }

        std::cout << "Error in computing ephi: " << ephierr << "\n";
        std::cout << "Error in computing etheta: " << ethetaerr << "\n";
        std::cout << "Error in computing er: " << ererr << "\n";

        std::cout << "Error in computing ephidPhi: " << ephidPhierr << "\n";
        std::cout << "Error in computing ephidTheta: " << ephidThetaerr << "\n";

        std::cout << "Error in computing ethetadTheta: " << ethetadThetaerr << "\n";
        std::cout << "Error in computing ethetadPhi: " << ethetadPhierr << "\n";

        std::cout << "Error in computing erdTheta: " << erdThetaerr << "\n";
        std::cout << "Error in computing erdPhi: " << erdPhierr << "\n";
        
    }


};

/**
 * @brief Functor for @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$. Also has functions for component extraction and differentiation.
*/
class VImag : public SphericalVectorField
{
    int m;
    int n;
    YImag yi;
    GImag gi;

public:

    /**
     * @brief Initializes the order and degree.
    */
    VImag(int a = 0, int b=0)
    {
        m = a;
        n = b;
        yi = YImag(m, n);
        gi = GImag(m, n);
    }

    /**
     * @brief Evaluates  @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$.
     * @return @f$ -(n+1) \mathfrak{I}(\boldsymbol{Y}_n^m)(x) \boldsymbol{e}_rs + \nabla_\Gamma \mathfrak{I}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See YImag::operator()() , GImag::operator()().
    */
    RectCoord operator()(const SphereCoord & x) const 
    {
        return -((double)n + 1.0) * yi(x.s) * e_r(x) + gi(x);
    }


    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$.
     * @return  @f$ -(n+1) \mathfrak{I}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See: YImag::operator()().
    */
    double eR(const SphereCoord & x) const
    {
        return -((double)n + 1.0) * yi(x.s);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$.
     * @return GImag::eTheta().
    */
    double eTheta(const SphereCoord & x) const
    {
        return gi.eTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$.
     * @return GImag::ePhi().
    */
    double ePhi(const SphereCoord & x) const
    {
        return gi.ePhi(x);
    }

    /**
    * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
    * @return GImag::eTheta_dTheta().
    */
    double eTheta_dTheta(const SphereCoord & x) const
    {
        return gi.eTheta_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GImag::eTheta_dPhi().
    */
    double eTheta_dPhi(const SphereCoord & x) const
    {
        return gi.eTheta_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GImag::ePhi_dPhi().
    */
    double ePhi_dPhi(const SphereCoord & x) const
    {
        return gi.ePhi_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return GImag::ePhi_dTheta().
    */
    double ePhi_dTheta(const SphereCoord& x) const
    {
        return gi.ePhi_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return  @f$ -(n+1)\frac{\partial \mathfrak{R}(\boldsymbol{Y}_n^m)}{\partial \theta}  @f$.
     *
     * See: YImag::dTheta()
    */
    double eR_dTheta(const SphereCoord& x) const
    {
        return -((double)n + 1.0) * yi.dTheta(x.s);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return  @f$ -(n+1)\frac{\partial \mathfrak{I}(\boldsymbol{Y}_n^m)}{\partial \phi}  @f$.
     *
     * See: YImag::dPhi()
    */
    double eR_dPhi(const SphereCoord & x) const
    {
        return -(double)(n + 1) * yi.dPhi(x.s);
    }

    /**
     * @brief resets the indices.
    */
    void reset(int a, int b)
    {
        m = a;
        n = b;
        yi = YImag(m, n);
        gi = GImag(m, n);
    }

    /**
     * @brief Compares the values of the derivatives given by member function to those given by generic numerical methods.
     *The numerical methods are given by the NdPhi, NdTheta classes, possibly with repeated applications.
    */
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
        double ephidPhierr = 0.0;
        double ephidThetaerr = 0.0;
        double erdPhierr = 0.0;
        double ethetadThetaerr = 0.0;
        double ethetadPhierr = 0.0;
        double erdThetaerr = 0.0;


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
                ephidPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;



                diff = eTheta_dTheta(x) - etdt(x);
                ethetadThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dTheta(x) - erdt(x);
                erdThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dPhi(x) - erdp(x);
                erdPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
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

        std::cout << "Error in computing ephidPhi: " << ephidPhierr << "\n";
        std::cout << "Error in computing ephidTheta: " << ephidThetaerr << "\n";

        std::cout << "Error in computing ethetadTheta: " << ethetadThetaerr << "\n";
        std::cout << "Error in computing ethetadPhi: " << ethetadPhierr << "\n";

        std::cout << "Error in computing erdTheta: " << erdThetaerr << "\n";
        std::cout << "Error in computing erdPhi: " << erdPhierr << "\n";

    }

};


/**
 * @brief Functor for @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$. Also has functions for component extraction and differentiation.
*/
class WReal : public SphericalVectorField
{
    int m;
    int n;
    YReal yr;
    GReal gr;

public:

    /**
     * @brief Initializes the order and degree.
    */
    WReal(int a = 0, int b = 1)
    {
        m = a;
        n = b;
        yr = YReal(m, n);
        gr = GReal(m, n);
    }

    /**
     * @brief Evaluates  @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$.
     * @return @f$ n \mathfrak{R}(\boldsymbol{Y}_n^m)(x) \boldsymbol{e}_rs + \nabla_\Gamma \mathfrak{R}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See YReal::operator()() , GReal::operator()().
    */
    RectCoord  operator()(const SphereCoord & x) const
    {
        return (double)n * yr(x.s) * e_r(x) + gr(x);
    }

    /**
     * @brief resets the indices.
    */
    void reset(int a, int b)
    {
        m = a;
        n = b;
        yr = YReal(m, n);
        gr = GReal(m, n);
    }


    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$.
     * @return  @f$ n \mathfrak{R}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See: YReal::operator()().
    */
    double eR(const SphereCoord & x) const
    {
        return n * yr(x.s);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$.
     * @return GReal::eTheta().
    */
    double eTheta(const SphereCoord & x) const
    {
        return gr.eTheta(x.s);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$.
     * @return GReal::ePhi().
    */
    double ePhi(const SphereCoord & x) const
    {
        return gr.ePhi(x);
    }

    /**
    * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
    * @return GReal::eTheta_dTheta().
    */
    double eTheta_dTheta(const SphereCoord & x) const
    {
        return gr.eTheta_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GReal::eTheta_dPhi().
    */
    double eTheta_dPhi(const SphereCoord & x) const
    {
        return gr.eTheta_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GReal::ePhi_dPhi().
    */
    double ePhi_dPhi(const SphereCoord & x) const
    {
        return gr.ePhi_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return GReal::ePhi_dTheta().
    */
    double ePhi_dTheta(const SphereCoord& x) const
    {
        return gr.ePhi_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return  @f$ n \frac{\partial \mathfrak{R}(\boldsymbol{Y}_n^m)}{\partial \theta}  @f$.
     *
     * See: YReal::dTheta()
    */
    double eR_dTheta(const SphereCoord& x) const
    {
        return (double)n * yr.dTheta(x.s);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return  @f$ n \frac{\partial \mathfrak{R}(\boldsymbol{Y}_n^m)}{\partial \phi}  @f$.
     *
     * See: YReal::dPhi()
    */
    double eR_dPhi(const SphereCoord & x) const
    {
        return (double)n * yr.dPhi(x.s);
    }

    /**
     * @brief Compares the values of the derivatives given by member function to those given by generic numerical methods.
     *The numerical methods are given by the NdPhi, NdTheta classes, possibly with repeated applications.
    */
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
        double ephidPhierr = 0.0;
        double ephidThetaerr = 0.0;
        double erdPhierr = 0.0;
        double ethetadThetaerr = 0.0;
        double ethetadPhierr = 0.0;
        double erdThetaerr = 0.0;


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
                ephidPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;



                diff = eTheta_dTheta(x) - etdt(x);
                ethetadThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dTheta(x) - erdt(x);
                erdThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dPhi(x) - erdp(x);
                erdPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
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

        std::cout << "Error in computing ephidPhi: " << ephidPhierr << "\n";
        std::cout << "Error in computing ephidTheta: " << ephidThetaerr << "\n";

        std::cout << "Error in computing ethetadTheta: " << ethetadThetaerr << "\n";
        std::cout << "Error in computing ethetadPhi: " << ethetadPhierr << "\n";

        std::cout << "Error in computing erdTheta: " << erdThetaerr << "\n";
        std::cout << "Error in computing erdPhi: " << erdPhierr << "\n";

    }

};

/**
 * @brief Functor for @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$. Also has functions for component extraction and differentiation.
*/
class WImag : public SphericalVectorField
{
    int m;
    int n;
    YImag yi;
    GImag gi;

public:


    /**
     * @brief Initializes the order and degree.
    */
    WImag(int a = 0, int b = 1)
    {
        m = a;
        n = b;
        yi = YImag(m, n);
        gi = GImag(m, n);
    }

    /**
     * @brief Evaluates  @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$.
     * @return @f$ n \mathfrak{I}(\boldsymbol{Y}_n^m)(x) \boldsymbol{e}_rs + \nabla_\Gamma \mathfrak{I}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See YImag::operator()() , GImag::operator()().
    */
    RectCoord operator()(const SphereCoord & x) const
    {
        return (double)n * yi(x.s) * e_r(x) + gi(x);
    }

    /**
     * @brief resets the indices.
    */
    void reset(int a, int b)
    {
        m = a;
        n = b;
        yi = YImag(m, n);
        gi = GImag(m, n);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$.
     * @return  @f$ n \mathfrak{I}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See: YImag::operator()().
    */
    double eR(const SphereCoord & x) const
    {
        return n * yi(x.s);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$.
     * @return GImag::eTheta().
    */
    double eTheta(const SphereCoord & x) const
    {
        return gi.eTheta(x.s);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$.
     * @return GImag::ePhi().
    */
    double ePhi(const SphereCoord & x) const
    {
        return gi.ePhi(x);
    }

    /**
    * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
    * @return GImag::eTheta_dTheta().
    */
    double eTheta_dTheta(const SphereCoord & x) const
    {
        return gi.eTheta_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GImag::eTheta_dPhi().
    */
    double eTheta_dPhi(const SphereCoord & x) const
    {
        return gi.eTheta_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GImag::ePhi_dPhi().
    */
    double ePhi_dPhi(const SphereCoord & x) const
    {
        return gi.ePhi_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return GImag::ePhi_dTheta().
    */
    double ePhi_dTheta(const SphereCoord& x) const
    {
        return gi.ePhi_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return  @f$ n \frac{\partial \mathfrak{R}(\boldsymbol{Y}_n^m)}{\partial \theta}  @f$.
     *
     * See: YImag::dTheta()
    */
    double eR_dTheta(const SphereCoord& x) const
    {
        return (double)n * yi.dTheta(x.s);
    }


    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return  @f$ n \frac{\partial \mathfrak{I}(\boldsymbol{Y}_n^m)}{\partial \phi}  @f$.
     *
     * See: YImag::dPhi()
    */
    double eR_dPhi(const SphereCoord& x) const
    {
        return (double)n * yi.dPhi(x.s);
    }

    /**
     * @brief Compares the values of the derivatives given by member function to those given by generic numerical methods.
     *The numerical methods are given by the NdPhi, NdTheta classes, possibly with repeated applications.
    */
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
        double ephidPhierr = 0.0;
        double ephidThetaerr = 0.0;
        double erdPhierr = 0.0;
        double ethetadThetaerr = 0.0;
        double ethetadPhierr = 0.0;
        double erdThetaerr = 0.0;


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
                ephidPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;



                diff = eTheta_dTheta(x) - etdt(x);
                ethetadThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dTheta(x) - erdt(x);
                erdThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dPhi(x) - erdp(x);
                erdPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
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

        std::cout << "Error in computing ephidPhi: " << ephidPhierr << "\n";
        std::cout << "Error in computing ephidTheta: " << ephidThetaerr << "\n";

        std::cout << "Error in computing ethetadTheta: " << ethetadThetaerr << "\n";
        std::cout << "Error in computing ethetadPhi: " << ethetadPhierr << "\n";

        std::cout << "Error in computing erdTheta: " << erdThetaerr << "\n";
        std::cout << "Error in computing erdPhi: " << erdPhierr << "\n";

    }
};



/**
 * @brief Functor for @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$. Also has functions for component extraction and differentiation.
*/
class XReal : public SphericalVectorField
{
    int m;
    int n;
    GReal gr;
    YReal yr;

public:

    /**
     * @brief Initializes the order and degree.
    */
    XReal(int a = 0, int b = 1)
    {
        m = a;
        n = b;

        yr = YReal(m, n);
        gr = GReal(m, n);
    }

    /**
     * @brief Evaluates  @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$.
     * @return @f$  \boldsymbol{e}_rs \times \nabla_\Gamma \mathfrak{R}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See GReal::operator()().
    */
    RectCoord operator()(const SphereCoord & x) const
    {
        return cross(e_r(x), gr(x));
    }

    /**
     * @brief resets the indices.
    */
    void reset(int a, int b)
    {
        m = a;
        n = b;
        gr = GReal(m, n);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$.
     * @return 0.
    */
    double eR(const SphereCoord & x) const
    {
        return 0;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$.
     * @return the negative of GReal::ePhi().
    */
    double eTheta(const SphereCoord & x) const
    {
        return -gr.ePhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$.
     * @return GReal::eTheta().
    */
    double ePhi(const SphereCoord & x) const
    {
        return gr.eTheta(x);
    }

    /**
    * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
    * @return The negative of GReal::ePhi_dTheta().
    */
    double eTheta_dTheta(const SphereCoord & x) const
    {
        return -gr.ePhi_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return The negative of GReal::ePhi_dPhi().
    */
    double eTheta_dPhi(const SphereCoord & x) const
    {
        return -gr.ePhi_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GReal::eTheta_dPhi().
    */
    double ePhi_dPhi(const SphereCoord & x) const
    {
        return gr.eTheta_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return GReal::eTheta_dTheta().
    */
    double ePhi_dTheta(const SphereCoord& x) const
    {
        return gr.eTheta_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return 0.
    */
    double eR_dTheta(const SphereCoord& x) const
    {
        return 0;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return 0.
    */
    double eR_dPhi(SphereCoord x) const
    {
        return 0;
    }

    /**
     * @brief Compares the values of the derivatives given by member function to those given by generic numerical methods.
     *The numerical methods are given by the NdPhi, NdTheta classes, possibly with repeated applications.
    */
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
        double ephidPhierr = 0.0;
        double ephidThetaerr = 0.0;
        double erdPhierr = 0.0;
        double ethetadThetaerr = 0.0;
        double ethetadPhierr = 0.0;
        double erdThetaerr = 0.0;


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
                ephidPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;



                diff = eTheta_dTheta(x) - etdt(x);
                ethetadThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dTheta(x) - erdt(x);
                erdThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dPhi(x) - erdp(x);
                erdPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
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

        std::cout << "Error in computing ephidPhi: " << ephidPhierr << "\n";
        std::cout << "Error in computing ephidTheta: " << ephidThetaerr << "\n";

        std::cout << "Error in computing ethetadTheta: " << ethetadThetaerr << "\n";
        std::cout << "Error in computing ethetadPhi: " << ethetadPhierr << "\n";

        std::cout << "Error in computing erdTheta: " << erdThetaerr << "\n";
        std::cout << "Error in computing erdPhi: " << erdPhierr << "\n";

    }
};

/**
 * @brief Functor for @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$. Also has functions for component extraction and differentiation.
*/
class XImag : public SphericalVectorField
{
    int m;
    int n;
    YImag yi;
    GImag gi;

public:

    /**
     * @brief Initializes the order and degree.
    */
    XImag(int a = 0, int b = 1)
    {
        m = a;
        n = b;

        yi = YImag(m, n);
        gi = GImag(m, n);
    }

    /**
     * @brief Evaluates  @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$.
     * @return @f$  \boldsymbol{e}_rs \times \nabla_\Gamma \mathfrak{I}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See GImag::operator()().
    */
    RectCoord operator()(const SphereCoord & x) const
    {
        return cross(e_r(x), gi(x));
    }

    /**
     * @brief resets the indices.
    */
    void reset(int a, int b)
    {
        m = a;
        n = b;
        gi = GImag(m, n);
    }
    
    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$.
     * @return 0.
    */
    double eR(const SphereCoord & x) const
    {
        return 0;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$.
     * @return The negative GImag::ePhi().
    */
    double eTheta(const SphereCoord & x) const
    {
        return -gi.ePhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$.
     * @return GImag::eTheta().
    */
    double ePhi(const SphereCoord & x) const
    {
        return gi.eTheta(x);
    }

    /**
    * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
    * @return The negative of GImag::ePhi_dTheta().
    */
    double eTheta_dTheta(const SphereCoord & x) const
    {
        return -gi.ePhi_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return The negative of GImag::ePhi_dPhi().
    */
    double eTheta_dPhi(const SphereCoord & x) const
    {
        return -gi.ePhi_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GImag::eTheta_dPhi().
    */
    double ePhi_dPhi(const SphereCoord & x) const
    {
        return gi.eTheta_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return GImag::eTheta_dTheta().
    */
    double ePhi_dTheta(const SphereCoord& x) const
    {
        return gi.eTheta_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return 0.
    */
    double eR_dTheta(const SphereCoord& x) const
    {
        return 0;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return 0.
    */
    double eR_dPhi(SphereCoord x) const
    {
        return 0;
    }

    /**
     * @brief Compares the values of the derivatives given by member function to those given by generic numerical methods.
     *The numerical methods are given by the NdPhi, NdTheta classes, possibly with repeated applications.
    */
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
        double ephidPhierr = 0.0;
        double ephidThetaerr = 0.0;
        double erdPhierr = 0.0;
        double ethetadThetaerr = 0.0;
        double ethetadPhierr = 0.0;
        double erdThetaerr = 0.0;


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
                ephidPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = ePhi_dTheta(x) - epdt(x);
                ephidThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;



                diff = eTheta_dTheta(x) - etdt(x);
                ethetadThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dTheta(x) - erdt(x);
                erdThetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dPhi(x) - erdp(x);
                erdPhierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
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

        std::cout << "Error in computing ephidPhi: " << ephidPhierr << "\n";
        std::cout << "Error in computing ephidTheta: " << ephidThetaerr << "\n";

        std::cout << "Error in computing ethetadTheta: " << ethetadThetaerr << "\n";
        std::cout << "Error in computing ethetadPhi: " << ethetadPhierr << "\n";

        std::cout << "Error in computing erdTheta: " << erdThetaerr << "\n";
        std::cout << "Error in computing erdPhi: " << erdPhierr << "\n";

    }

};

/**
 * @brief Single term in the VSHSeries class. Contains each Vector Spherical Harmonic Basis function for as given index. 
*/
class VSHTerm : public SphericalVectorField
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



    VSHTerm() {}

    /**
     * @brief Constructs a term in the series given the normalized coefficients.
     * @param rhohats in order, the VReal, VImag, WReal, WImag, XReal, XImag coefficients to be assigned.
    */
    VSHTerm(int m0, int n0, std::array<double, 6> rhohats)
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

    VSHTerm(const VSHTerm& term)
    {
        m = term.m;
        n = term.n;

        vr = VReal(m, n);
        vi = VImag(m, n);
        wr = WReal(m, n);
        wi = WImag(m, n);
        xr = XReal(m, n);
        xi = XImag(m, n);

        rhohatVr = term.rhohatVr;
        rhohatVi = term.rhohatVi;
        rhohatWr = term.rhohatWr;
        rhohatWi = term.rhohatWi;
        rhohatXr = term.rhohatXr;
        rhohatXi = term.rhohatXi;
    }

    VSHTerm operator =(const VSHTerm& term)
    {
        m = term.m;
        n = term.n;

        vr = VReal(m, n);
        vi = VImag(m, n);
        wr = WReal(m, n);
        wi = WImag(m, n);
        xr = XReal(m, n);
        xi = XImag(m, n);

        rhohatVr = term.rhohatVr;
        rhohatVi = term.rhohatVi;
        rhohatWr = term.rhohatWr;
        rhohatWi = term.rhohatWi;
        rhohatXr = term.rhohatXr;
        rhohatXi = term.rhohatXi;

        return *this;
    }

    /**
     * @brief Addition operator between terms.
     * @param t The term to add.
     * @return A term whose coefficients are sums of the two terms.
    */
    VSHTerm operator +(const VSHTerm & t) const
    {
        std::array<double, 6> temp;

        temp[0] = rhohatVr + t.rhohatVr;
        temp[1] = rhohatVi + t.rhohatVi;
        temp[2] = rhohatWr + t.rhohatWr;
        temp[3] = rhohatWi + t.rhohatWi;
        temp[4] = rhohatXr + t.rhohatXr;
        temp[5] = rhohatXi + t.rhohatXi;

        return VSHTerm(m, n, temp);
    }

    /**
     * @brief Subtraction operator between terms.
     * @param t The term to add.
     * @return A term whose coefficients are the difference of the two terms.
    */
    VSHTerm operator -(const VSHTerm & t) const
    {
        std::array<double, 6> temp;

        temp[0] = rhohatVr - t.rhohatVr;
        temp[1] = rhohatVi - t.rhohatVi;
        temp[2] = rhohatWr - t.rhohatWr;
        temp[3] = rhohatWi - t.rhohatWi;
        temp[4] = rhohatXr - t.rhohatXr;
        temp[5] = rhohatXi - t.rhohatXi;

        return VSHTerm(m, n, temp);
    }

    /**
     * @brief Unry minus operator.
     * @return A term whose coefficients are the negative of the calling term.
    */
    VSHTerm operator -() const
    {
        std::array<double, 6> temp;

        temp[0] = -rhohatVr;
        temp[1] = -rhohatVi;
        temp[2] = -rhohatWr;
        temp[3] = -rhohatWi;
        temp[4] = -rhohatXr;
        temp[5] = -rhohatXi;

        return VSHTerm(m, n, temp);
    }


    VSHTerm operator *(const double& a) const
    {
        std::array<double, 6> temp;

        temp[0] = rhohatVr *a;
        temp[1] = rhohatVi * a;
        temp[2] = rhohatWr * a;
        temp[3] = rhohatWi * a;
        temp[4] = rhohatXr * a;
        temp[5] = rhohatXi * a;

        return VSHTerm(m, n, temp);
    }

    /**
     * @brief Computes the inner product between terms.
    */
    double dot(const VSHTerm & s) const
    {
        if (m != s.m || n != s.n)
            return 0.0;
        else
            return rhohatVr * s.rhohatVr + rhohatVi * s.rhohatVi + rhohatWr * s.rhohatWr + rhohatWi * s.rhohatWi + rhohatXr * s.rhohatXr + rhohatXi * s.rhohatXi;
    }

    /**
     * @brief Evaluates the term at the point x.
     * @return In terms of complex numbers, @f$   \hat{\rho}_{n,m}^{V*}\boldsymbol{V}_n^m(\theta, \phi) +  \hat{\rho}_{n,m}^{W*}\boldsymbol{W}_n^m(\theta, \phi)
    +\hat{\rho}_{n,m}^{X*}\boldsymbol{X}_n^m(\theta, \phi)  @f$. The evaluation uses real functions.
    */
    RectCoord operator()(const SphereCoord & x) const
    {
        if (n == 0)
            return 2 * rhohatVr * vr(x);

        RectCoord temp = rhohatVr * vr(x) + rhohatVi * vi(x);
        temp = temp + (rhohatWr * wr(x) + rhohatWi * wi(x));
        return temp + (rhohatXr * xr(x) + rhohatXi * xi(x));
    }

};

/**
 * @brief The Vector Spherical Harmonic Series Functor. Computes and evaluates the fourier series representation of a SphericalVectorField.
*/
class VSHSeries : public VectorFieldSum
{

private:



    Legendre* P;

    

public:
    int N;
    RectCoord center;
    std::vector<std::vector<VSHTerm>> terms;



    /**
     * @brief Constructs the series. n is the maximum value of the degree. there will be @f$ n(n+1)/2 @f$ terms total.
     * @param c The center of the coordinate system of the series.
    */
    VSHSeries(int n , RectCoord c = RectCoord()): P(&legendrePoly)
    {
        N = n;

        terms.resize(n + 1);

        for (int i = 0; i <= n; i++)
            terms[i].resize(i + 1);

        center = c;



    }

    VSHSeries(const VSHSeries & vsh) : P(&legendrePoly)
    {
        N = vsh.N;

        terms = vsh.terms;

        center = vsh.center;

    }

    VSHSeries(VSHSeries&& vsh) noexcept : P(&legendrePoly)
    {
        N = vsh.N;

        terms = std::move(vsh.terms);

        center = vsh.center;

    }

    VSHSeries & operator =(const VSHSeries& vsh)
    {
        N = vsh.N;

        terms = vsh.terms;

        center = vsh.center;

        return *this;

    }

    VSHSeries& operator =(VSHSeries&& vsh) noexcept
    {
        N = vsh.N;

        terms = std::move(vsh.terms);

        center = vsh.center;

        return *this;

    }

    void operator ~()
    {
    }

    /**
     * @brief Approximates the passed SphericalVectorField.
     * @param r the function to approximate.
     * @param n number of terms in the series.
    */
    void approximate(SphericalVectorField* r, int n = 0)
    {

        if (n > 0)
        {
            terms.resize(n + 1);

            for (int i = 0; i <= n; i++)
                terms[i].resize( i + 1);

        }




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
                    double n_ = double(n);
                    vr.reset(m, n);
                    rhohats[0] =  L2InnerProduct(r, &vr, center) / (double)((2 * n_ + 1) * (n_ + 1));

                    vi.reset(m, n);
                    rhohats[1] = 0.0;

                    wr.reset(m, n);
                    rhohats[2] =  L2InnerProduct(r, &wr , center) / (double)((2 * n_ + 1) * n_);

                    wi.reset(m, n);
                    rhohats[3] = 0.0;

                    xr.reset(m, n);
                    rhohats[4] =  L2InnerProduct(r, &xr, center) / (double)((n_ + 1) * n_);

                    xi.reset(m, n);
                    rhohats[5] = 0.0;

                    terms[n][m] = VSHTerm(m, n, rhohats);
                }
                else if (m > 0)
                {
                    double n_ = n;
                    vr.reset(m, n);
                    rhohats[0] = 2.0 * L2InnerProduct(r, &vr, center) / (double)((2 * n_ + 1) * (n_ + 1));

                    vi.reset(m, n);
                    rhohats[1] = 2.0 * L2InnerProduct(r, &vi, center) / (double)((2 * n_ + 1) * (n_ + 1));

                    wr.reset(m, n);
                    rhohats[2] = 2.0 * L2InnerProduct(r, &wr, center) / (double)((2 * n_ + 1) * n_);

                    wi.reset(m, n);
                    rhohats[3] = 2.0 * L2InnerProduct(r, &wi, center) / (double)((2 * n_ + 1) * n_);

                    xr.reset(m, n);
                    rhohats[4] = 2.0 * L2InnerProduct(r, &xr, center) / (double)((n_ + 1) * n_);

                    xi.reset(m, n);
                    rhohats[5] = 2.0 * L2InnerProduct(r, &xi, center) / (double)((n_ + 1) * n_);

                    terms[n][m] = VSHTerm(m, n, rhohats);
                }
                
                append(&terms[n][m]);
            }
    }

    /**
     * @brief Approximates A function defined on the surface using data in quadrature nodes.
     * @param data boundary data.
     * @param n maximum degree of the series.
    */
    void approximate(SphereData & data, int n = 0)
    {

        if (n > 0)
        {
            terms.resize(n + 1);

            for (int i = 0; i <= n; i++)
                terms[i].resize( i + 1);

        }




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
                    double n_ = double(n);
                    vr.reset(m, n);
                    rhohats[0] = L2InnerProductDiscrete(data, &vr,center) / ((2 * n_ + 1) * (n_ + 1));

                    vi.reset(m, n);
                    rhohats[1] = L2InnerProductDiscrete(data, &vi,  center) / ((2 * n_ + 1) * (n_ + 1));

                    wr.reset(m, n);
                    rhohats[2] = L2InnerProductDiscrete(data, &wr, center) / ((2 * n_ + 1) * n_);

                    wi.reset(m, n);
                    rhohats[3] = L2InnerProductDiscrete(data, &wi, center) / (double)((2 * n_ + 1) * n_);

                    xr.reset(m, n);
                    rhohats[4] = L2InnerProductDiscrete(data, &xr, center) / (double)((n_ + 1) * n_);

                    xi.reset(m, n);
                    rhohats[5] = L2InnerProductDiscrete(data, &xi, center) / (double)((n_ + 1) * n_);

                    terms[n][m] = VSHTerm(m, n, rhohats);
                }
                else if (m > 0)
                {
                    double n_ = double(n);
                    vr.reset(m, n);
                    rhohats[0] = 2.0 * L2InnerProductDiscrete(data, &vr, center) / (double)((2 * n_ + 1) * (n_ + 1));
                    vi.reset(m, n);
                    rhohats[1] = 2.0 * L2InnerProductDiscrete(data, &vi, center) / (double)((2 * n_ + 1) * (n_ + 1));

                    wr.reset(m, n);
                    rhohats[2] = 2.0 * L2InnerProductDiscrete(data, &wr, center) / (double)((2 * n_ + 1) * n_);

                    wi.reset(m, n);
                    rhohats[3] = 2.0 * L2InnerProductDiscrete(data, &wi, center) / (double)((2 * n_ + 1) * n_);
                    xr.reset(m, n);
                    rhohats[4] = 2.0 * L2InnerProductDiscrete(data, &xr, center) / (double)((n_ + 1) * n_);

                    xi.reset(m, n);
                    rhohats[5] = 2.0 * L2InnerProductDiscrete(data, &xi, center) / (double)((n_ + 1) * n_);

                    terms[n][m] = VSHTerm(m, n, rhohats);
                }
                append(&terms[n][m]);
            }
    }

    /**
     * @brief Adds two VSHSeries.
     * @param s the second series to add.
    */
    VSHSeries operator +(const VSHSeries & s) const
    {
        VSHSeries temp(N, center);

        temp.terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            temp.terms[i].resize(i + 1);



        for(int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                temp.terms[n][m] = terms[n][m] + s.terms[n][m];

        return temp;
            
    }

    /**
     * @brief subtracts two VSHSeries.
     * @param s the second series to add.
    */
    VSHSeries operator -(const VSHSeries & s) const
    {
        VSHSeries temp(N,  center);

        temp.terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            temp.terms[i].resize(i + 1);



        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                temp.terms[n][m] = terms[n][m] - s.terms[n][m];

        return temp;

    }

    VSHSeries operator -() const
    {
        VSHSeries temp(N, center);

        temp.terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            temp.terms[i].resize(i + 1);



        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                temp.terms[n][m] = -terms[n][m];

        return temp;
    }

    VSHSeries operator *(const double& a) const
    {
        VSHSeries temp(N, center);

        temp.terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            temp.terms[i].resize(i + 1);



        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                temp.terms[n][m] = terms[n][m] * a;

        return temp;
    }

    /**
     * @brief Computes the @f$ L^2 @f$ inner product of two series.
     * @param series the second series in the product.
    */
    double dot(const VSHSeries & series) const
    {
        double total = 0.0;

        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                total += terms[n][m].dot(series.terms[n][m]);

        return total;
    }

    /**
     * @brief evaluates the series. Recasts x in the spherical coordinate system of the series.
     * @return 
    */
    RectCoord operator()(const SphereCoord & x) const
    {

        SphereCoord temp = recenter(x, center);

        P->populate(cos(temp.s.theta), N, N);

        return VectorFieldSum::operator()(temp);
    }
};