#pragma once


#include <cmath>
#include <complex>
#include <vector>

#define MATHPI 4.0 * atan(1)

//structures for the different coordinate systems. Surface coordinates are for the unit sphere. Polar coordinates are actually Spherical Coordinates.
//all coordinates are complex.
struct SurfaceCoord
{
public:

    //zenith angle.
    std::complex<double> phi = 0;
    //azimuthal angle.
    std::complex<double> theta = 0;

    //addition for spherical coordinates.(note, addition here is essentially rotation).
    SurfaceCoord operator +(SurfaceCoord s2)
    {
        SurfaceCoord temp;
        temp.phi = phi + s2.phi;
        temp.theta = theta + s2.theta;

        return temp;
    }

    SurfaceCoord()
    {
        theta = 0;
        phi = 0;
    }

    SurfaceCoord(std::complex<double> a, std::complex<double> b)
    {
        theta = a;
        phi = b;
    }
};
struct SphereCoord
{
    std::complex<double> rho;
    SurfaceCoord s;

    SphereCoord()
    {
        rho = 1;
    }
    SphereCoord(std::complex<double> r, SurfaceCoord s1)
    {
        rho = r;
        s = s1;
    }
    SphereCoord(SurfaceCoord s1)
    {
        rho = 1;
        s = s1;
    }







};
struct RectCoord
{
    std::complex<double> x;
    std::complex<double> y;
    std::complex<double> z;

    RectCoord() : x(0), y(0), z(0) {}
    RectCoord(std::complex<double>x0, std::complex<double> y0, std::complex<double> z0) :x(x0), y(y0), z(z0) {}
    RectCoord(SphereCoord r)
    {
        x = r.rho * sin(r.s.theta) * cos(r.s.phi);
        y = r.rho * sin(r.s.theta) * sin(r.s.phi);
        z = r.rho * cos(r.s.theta);
    }
};


//Arithmetic operators for rectangular coordinates.
RectCoord operator +(RectCoord a, RectCoord b)
{
    return RectCoord(a.x + b.x, a.y + b.y, a.z + b.z);
}
RectCoord operator -(RectCoord a, RectCoord b)
{
    return RectCoord(a.x - b.x, a.y - b.y, a.z - b.z);
}
RectCoord operator *(std::complex<double> c, RectCoord a)
{
    return RectCoord(c * a.x, c * a.y, c * a.z);
}
RectCoord operator *(RectCoord a , std::complex<double> c)
{
    return RectCoord(c * a.x, c * a.y, c * a.z);
}

//dot product for rectangular coordinates.
std::complex<double> dot(RectCoord a, RectCoord b)
{
    return a.x * conj(b.x) + a.y * conj(b.y) + a.z * conj(b.z);
}

//vector cross product.
RectCoord cross(RectCoord a, RectCoord b)
{
    RectCoord temp;
    temp.x = conj(a.y * b.z - a.z * b.y);
    temp.y = conj(a.z * b.x - a.x * b.z);
    temp.z = conj(a.x * b.y - a.y * b.x);

    return temp;
}

//Spherical and rectangular conversion functions.
SphereCoord RectToSphere(RectCoord x)
{
    SphereCoord temp;
    
    temp.rho = sqrt(x.x * conj(x.x) + x.y * conj(x.y) + x.z * conj(x.z));

    temp.s.phi = atan2(x.y.real(), x.x.real());
   
    std::complex<double> planenorm = sqrt(x.x * conj(x.x) + x.y * conj(x.y));

    temp.s.theta = atan2(planenorm.real(), x.z.real());
    

    return temp;

}
RectCoord SphereToRect(SphereCoord r)
{
    return RectCoord(r);
}
SurfaceCoord RectToSurface(RectCoord x)
{
    return RectToSphere(x).s;
}
RectCoord SurfaceToRect(SurfaceCoord s)
{
    SphereCoord temp(1, s);
    return SphereToRect(temp);
}

//Polar coordinate arithmetic operators.
SphereCoord operator +(SphereCoord r, SphereCoord q)
{
    return RectToSphere(SphereToRect(r) + SphereToRect(q));
}
SphereCoord operator -(SphereCoord r, SphereCoord q)
{
    return RectToSphere(SphereToRect(r) - SphereToRect(q));
}
SphereCoord operator *(std::complex<double> a, SphereCoord x)
{
    return SphereCoord(a * x.rho, SurfaceCoord(x.s));
}

//Polar Coordinate dot product operator.
std::complex<double> dot(SphereCoord x, SphereCoord y)
{
    return dot(SphereToRect(x), SphereToRect(y));
}

//returns the norm squared of the given polar coordinates.
double norm(SphereCoord x)
{
    return norm(x.rho);
}

//Encapulation of a Surface to scalar function. In all such encapsulations, to create a function, create a derived class and implement the function operator().
//Then create a variable of the derived type and "call it". See code in main file for examples. In particular the rho function(immediately before main function)
// is the most simple.
class SurfaceScalarFunction
{
public:
    virtual std::complex<double> operator()(SurfaceCoord s) = 0;

};

//Computes the numerical gradient of a Surface to Scalar function using a symmetric difference with specified stepsize.
//kappa exists to prevent division by zero and should be set small for accuracy.
SphereCoord surfaceGrad(SurfaceScalarFunction* u, SurfaceCoord s, double stepsize)
{ 
    
    //Construct unit vectors.
    SphereCoord ephi , etheta;
    ephi  = RectToSphere(RectCoord(-sin(s.phi),cos(s.phi) ,0));
    etheta = RectToSphere(RectCoord(cos(s.theta) * cos(s.phi) , cos(s.theta) * sin(s.phi) , -sin(s.theta)));

    //evaluate partial derivative with respect to theta.
    SurfaceCoord phiPlus = SurfaceCoord(s.theta, s.phi + stepsize / 2);
    SurfaceCoord phiMinus = SurfaceCoord(s.theta, s.phi - stepsize / 2);

    std::complex<double> dphi = 0;
    dphi = (*u)(phiPlus) - (*u)(phiMinus);
    dphi /= sin(s.theta) * stepsize;

    //add contribution to gradient.
    SphereCoord temp;
    temp = dphi * ephi;

    //evaulate partial derivative with respect to phi
    SurfaceCoord thetaPlus = SurfaceCoord(s.theta + stepsize / 2, s.phi);
    SurfaceCoord thetaMinus = SurfaceCoord(s.theta - stepsize / 2, s.phi);

    std::complex<double> dtheta = 0;
    dtheta = (*u)(phiPlus) - (*u)(phiMinus);
    dtheta /= stepsize;

    //add contribution to gradient and return.
    temp = temp + dtheta * etheta;
    return temp;

    

}

//Encapsulates a Spherical coordinate to Scalar function.
class SphericalScalarFunction
{
public:
    virtual std::complex<double> operator()(SphereCoord) = 0;
};

//Encapsulates a Spherical coordinate to Spherical Coordinate function.
class SphericalVectorField
{
public:
    virtual SphereCoord operator()(SphereCoord) = 0;
};

//Serves as a container allowing a Sum of vector fields to be called as a single function. to add a new function to the sum,
// call the append function and pass an encapsulated function to it.
class VectorFieldSum : public SphericalVectorField
{
private:

    int numterms;
    std::vector<SphericalVectorField* > fs;
public:

    VectorFieldSum() { numterms = 0; }
    VectorFieldSum(SphericalVectorField& f, SphericalVectorField& g)
    {
        numterms = 2;
        fs.push_back(&f);
        fs.push_back(&g);
    }

    VectorFieldSum(VectorFieldSum f, SphericalVectorField& g)
    {
        numterms = f.numterms + 1;
        fs = f.fs;
        fs.push_back(&g);
    }

    VectorFieldSum(SphericalVectorField& f, VectorFieldSum g)
    {
        numterms = g.numterms + 1;
        fs.resize(numterms);
        fs.push_back(&f);
        for (int i = 0; i < g.numterms; i++)
            fs[i + 1] = g.fs[i];

    }

    VectorFieldSum(VectorFieldSum f, VectorFieldSum g)
    {
        numterms = g.numterms + f.numterms;

        fs = f.fs;
        fs.resize(numterms);

        for (int i = 0; i < g.numterms; i++)
            fs[i + f.numterms] = g.fs[i];


    }

    void append(SphericalVectorField& f)
    {
        fs.push_back(&f);
        numterms++;
    }


    SphereCoord operator()(SphereCoord x)
    {
        SphereCoord temp = (0.0, SurfaceCoord(0.0, 0.0));
        for (int i = 0; i < numterms; i++)
        {
            temp = temp + (*(fs[i]))(x);
        }

        return temp;
    }

};


//Computes the L2 inner product of two vector fields on the Surface of a unit sphere using a trapezoidal rule in both coordinates.
std::complex<double> L2InnerProduct(SphericalVectorField* f1, SphericalVectorField* f2, int n)
{

    double GLnodes[16] = { -0.0950125098376374  , 0.0950125098376374  , -0.2816035507792589 , 0.2816035507792589 ,  -0.4580167776572274  , 0.4580167776572274 
                          - 0.6178762444026438  , 0.6178762444026438 , -0.7554044083550030 , 0.7554044083550030 , -0.8656312023878318 , 0.8656312023878318  ,
                          -0.9445750230732326 , 0.9445750230732326  , -0.9894009349916499  , 0.9894009349916499 };

    double GLweights[16] = { 0.1894506104550685 , 0.1894506104550685 , 0.1826034150449236  , 0.1826034150449236  , 0.1691565193950025  , 0.1691565193950025 ,
                            0.1495959888165767  , 0.1495959888165767 , 0.1246289712555339  , 0.1246289712555339  , 0.0951585116824928  , 0.0951585116824928 ,
                            0.0622535239386479  , 0.0622535239386479 , 0.0271524594117541  , 0.0271524594117541 };
    std::complex<double> total = 0;

    for (double p = 0.0; p < 2 * MATHPI; p += 2.0 * MATHPI / n)
        for (int i = 0; i < 16; i++)
        {
            SurfaceCoord s( MATHPI / 2.0 * (GLnodes[i] + 1), p);
            SphereCoord x(s);
            total += sin(s.theta) * GLweights[i] * dot((*f1)(x), (*f2)(x)) *  (MATHPI / n) * MATHPI;
        }

    return total;
}

//computes the L2 difference of two vector fields on the Surface of a unit sphere using a trapezoidal rule in both coordinates.
double L2Difference(SphericalVectorField* f1, SphericalVectorField* f2, int n)
{
    double GLnodes[16] = { -0.0950125098376374  , 0.0950125098376374  , -0.2816035507792589 , 0.2816035507792589 ,  -0.4580167776572274  , 0.4580167776572274
                          - 0.6178762444026438  , 0.6178762444026438 , -0.7554044083550030 , 0.7554044083550030 , -0.8656312023878318 , 0.8656312023878318  ,
                          -0.9445750230732326 , 0.9445750230732326  , -0.9894009349916499  , 0.9894009349916499 };

    double GLweights[16] = { 0.1894506104550685 , 0.1894506104550685 , 0.1826034150449236  , 0.1826034150449236  , 0.1691565193950025  , 0.1691565193950025 ,
                            0.1495959888165767  , 0.1495959888165767 , 0.1246289712555339  , 0.1246289712555339  , 0.0951585116824928  , 0.0951585116824928 ,
                            0.0622535239386479  , 0.0622535239386479 , 0.0271524594117541  , 0.0271524594117541 };
    std::complex<double> total = 0;

    for (double p = 0.0; p < 2 * MATHPI; p += 2.0 * MATHPI / n)
        for (int i = 0; i < 16; i++)
        {
            SurfaceCoord s(MATHPI  / 2.0 * (GLnodes[i] + 1), p);
            SphereCoord x(s);
            SphereCoord diff = (*f1)(x) - (*f2)(x);
            total += sin(s.theta) * GLweights[i] * dot(diff, diff) * (MATHPI / n) * MATHPI / 2.0;
        }

    return sqrt(total.real());
}

//computes the sup norm of two vector fields on the surface of a unit sphere.
double LInfDifference(SphericalVectorField* f1, SphericalVectorField* f2, int n)
{
    double sup = 0;
    double curr = 0;

    for (double t = 0.0; t < MATHPI; t += MATHPI / n)
        for (double p = 0.0; p < 2.0 * MATHPI; p += 2.0 * MATHPI / n)
        {
            SurfaceCoord s(t, p);
            SphereCoord x(s);
            RectCoord diff = RectCoord((*f1)(x)) - RectCoord((*f2)(x));
            curr = sqrt(dot(diff , diff)).real();
            if (curr > sup)
                sup = curr;

        }

    return sup;
}