#pragma once


#include <cmath>
#include <complex>
#include <vector>

#define MATHPI 4 * atan(1)

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
        phi = a;
        theta = b;
    }
};
struct PolarCoord
{
    std::complex<double> rho;
    SurfaceCoord s;

    PolarCoord()
    {
        rho = 1;
    }
    PolarCoord(std::complex<double> r, SurfaceCoord s1)
    {
        rho = r;
        s = s1;
    }
    PolarCoord(SurfaceCoord s1)
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
    RectCoord(PolarCoord r)
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
    temp.x = a.y * b.z - a.z * b.y;
    temp.y = a.z * b.x - a.x * b.z;
    temp.z = a.x * b.y - a.y * b.x;

    return temp;
}

//Spherical and rectangular conversion functions.
PolarCoord RectToSphere(RectCoord x)
{
    PolarCoord temp;
    temp.rho = sqrt(x.x * x.x + x.y * x.y + x.z * x.z);

    std::complex<double> slope = x.y / x.x;
    if (norm(x.x) > 10e-9)
        temp.s.phi = atan(slope);
    else if (abs(x.y) > 10e-9)
        temp.s.phi = MATHPI / 2 * x.y / abs(x.y);
    else
        temp.s.phi = 0;
    if (norm(temp.rho) > 10e-9)

        temp.s.theta = acos(x.z / temp.rho);
    else
        temp.s.theta = 0;

    return temp;

}
RectCoord SphereToRect(PolarCoord r)
{
    return RectCoord(r);
}
SurfaceCoord RectToSurface(RectCoord x)
{
    return RectToSphere(x).s;
}
RectCoord SurfaceToRect(SurfaceCoord s)
{
    PolarCoord temp(1, s);
    return SphereToRect(temp);
}

//Polar coordinate arithmetic operators.
PolarCoord operator +(PolarCoord r, PolarCoord q)
{
    return RectToSphere(SphereToRect(r) + SphereToRect(q));
}
PolarCoord operator -(PolarCoord r, PolarCoord q)
{
    return RectToSphere(SphereToRect(r) - SphereToRect(q));
}
PolarCoord operator *(std::complex<double> a, PolarCoord x)
{
    return PolarCoord(a * x.rho, SurfaceCoord(x.s));
}

//Polar Coordinate dot product operator.
std::complex<double> dot(PolarCoord x, PolarCoord y)
{
    return dot(SphereToRect(x), SphereToRect(y));
}

//returns the norm squared of the given polar coordinates.
double norm(PolarCoord x)
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
PolarCoord surfaceGrad(SurfaceScalarFunction* u, SurfaceCoord s, double stepsize, double kappa)
{ 
    //Construct unit vectors.
    PolarCoord ephi , etheta;
    ephi  = RectToSphere(RectCoord(-sin(s.phi),cos(s.phi) ,0));
    etheta = RectToSphere(RectCoord(cos(s.theta) * cos(s.phi) , cos(s.theta) * sin(s.phi) , -sin(s.theta)));

    //evaluate partial derivative with respect to theta.
    SurfaceCoord thetaPlus = SurfaceCoord(s.phi, s.theta + stepsize / 2);
    SurfaceCoord thetaMinus = SurfaceCoord(s.phi, s.theta - stepsize / 2);

    std::complex<double> dtheta = 0;
    dtheta = (*u)(thetaPlus) - (*u)(thetaMinus);
    dtheta /= stepsize;

    //add contribution to gradient.
    PolarCoord temp;
    temp = temp + dtheta * etheta;

    //evaulate partial derivative with respect to phi
    SurfaceCoord phiPlus = SurfaceCoord(s.phi + stepsize / 2, s.theta);
    SurfaceCoord phiMinus = SurfaceCoord(s.phi - stepsize / 2, s.theta);

    std::complex<double> dphi = 0;
    dphi += (*u)(phiPlus) - (*u)(phiMinus);
    dphi /= sin(s.theta) * stepsize + kappa * stepsize;

    //add contribution to gradient and return.
    temp = temp + dphi * etheta;
    return temp;
}

//Encapsulates a Spherical coordinate to Scalar function.
class SphericalScalarFunction
{
public:
    virtual std::complex<double> operator()(PolarCoord) = 0;
};

//Encapsulates a Spherical coordinate to Spherical Coordinate function.
class SphericalVectorField
{
public:
    virtual PolarCoord operator()(PolarCoord) = 0;
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


    PolarCoord operator()(PolarCoord x)
    {
        PolarCoord temp = (0.0, SurfaceCoord(0.0, 0.0));
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
    std::complex<double> total = 0;

    for (double p = 0.0; p < MATHPI; p += MATHPI / n)
        for (double t = 0.0; t < 2.0 * MATHPI; t += 2.0 * MATHPI / n)
        {
            SurfaceCoord s(p, t);
            PolarCoord x(s);
            total += sin(p) * dot((*f1)(x), (*f2)(x)) * 2.0 * (MATHPI / n) * (MATHPI / n);
        }

    return total;
}

//computes the L2 difference of two vector fields on the Surface of a unit sphere using a trapezoidal rule in both coordinates.
double L2Difference(SphericalVectorField* f1, SphericalVectorField* f2, int n)
{
    std::complex<double> total = 0;

    for (double p = 0.0; p < MATHPI; p += MATHPI / n)
        for (double t = 0.0; t < 2.0 * MATHPI; t += 2.0 * MATHPI / n)
        {
            SurfaceCoord s(p, t);
            PolarCoord x(s);
            RectCoord diff = RectCoord((*f1)(x)) - RectCoord((*f2)(x));
            total += sin(p) * dot(diff, diff) * 2.0 * (MATHPI / n) * (MATHPI / n);
        }

    return sqrt(total.real());
}

//computes the sup norm of two vector fields on the surface of a unit sphere.
double LInfDifference(SphericalVectorField* f1, SphericalVectorField* f2, int n)
{
    double sup = 0;
    double curr = 0;

    for (double p = 0.0; p < MATHPI; p += MATHPI / n)
        for (double t = 0.0; t < 2.0 * MATHPI; t += 2.0 * MATHPI / n)
        {
            SurfaceCoord s(p, t);
            PolarCoord x(s);
            RectCoord diff = RectCoord((*f1)(x)) - RectCoord((*f2)(x));
            curr = sqrt(dot(diff , diff)).real();
            if (curr > sup)
                sup = curr;

        }

    return sup;
}