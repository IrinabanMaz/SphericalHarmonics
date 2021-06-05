#pragma once


#include <cmath>
#include <complex>
#include <vector>

const double PI = 4 * atan(1);


//structures for the different coordinate systems. Surface coordinates are for the unit sphere. Polar coordinates are actually Spherical Coordinates.
//all coordinates are complex.
struct SurfaceCoord
{
public:

    //zenith angle.
    double phi = 0;
    //azimuthal angle.
    double theta = 0;

    //addition for spherical coordinates.(note, addition here is essentially rotation).
    SurfaceCoord operator +(SurfaceCoord s2)
    {
        SurfaceCoord temp;
        temp.phi = phi + s2.phi;
        temp.theta = theta + s2.theta;

        return temp;
    }

    SurfaceCoord():
        theta(0),
        phi(0){}
    

    SurfaceCoord(double a, double b):
       theta(a),
        phi(b){}
   
};
struct SphereCoord
{
    double rho;
    SurfaceCoord s;

    SphereCoord()
    {
        rho = 1;
    }
    SphereCoord(double r, SurfaceCoord s1)
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
    double x;
    double y;
    double z;

    RectCoord():
        x(0), 
        y(0),
        z(0){}
    RectCoord(double x0, double y0, double z0):
        x(x0),
        y(y0),
        z(z0){}
    RectCoord(SphereCoord r):
        x(r.rho * sin(r.s.theta) * cos(r.s.phi)),
        y(r.rho * sin(r.s.theta) * sin(r.s.phi)),
        z(r.rho * cos(r.s.theta)){}
    
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
RectCoord operator *(double c, RectCoord a)
{
    return RectCoord(c * a.x, c * a.y, c * a.z);
}
RectCoord operator *(RectCoord a , double c)
{
    return RectCoord(c * a.x, c * a.y, c * a.z);
}

//dot product for rectangular coordinates.
double dot(RectCoord a, RectCoord b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
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
SphereCoord RectToSphere(RectCoord x)
{
    SphereCoord temp;
    
    temp.rho = sqrt(dot(x , x));
    
    temp.s.phi = atan2(x.y , x.x);
    
    double planenorm = sqrt(x.x * x.x + x.y * x.y);

    temp.s.theta = atan2(planenorm, x.z);
    

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
SphereCoord operator *(double a, SphereCoord x)
{
    if (a >= 0)
        return SphereCoord(a * x.rho, SurfaceCoord(x.s));
    else if (x.s.phi < PI)
        return SphereCoord(-a * x.rho, SurfaceCoord(x.s.theta, x.s.phi + PI));
    else
        return SphereCoord(-a * x.rho, SurfaceCoord(x.s.theta, x.s.phi - PI));
}

//Polar Coordinate dot product operator.
double dot(SphereCoord x, SphereCoord y)
{
    return dot(SphereToRect(x), SphereToRect(y));
}

SphereCoord cross(SphereCoord x, SphereCoord y)
{
    return RectToSphere(cross(SphereToRect(x), SphereToRect(y)));
}

//returns the norm squared of the given polar coordinates.
double norm(SphereCoord x)
{
    return x.rho;
}

//Encapulation of a Surface to scalar function. In all such encapsulations, to create a function, create a derived class and implement the function operator().
//Then create a variable of the derived type and "call it". See code in main file for examples. In particular the rho function(immediately before main function)
// is the most simple.
class SurfaceScalarFunction
{
private:
    double (*f)(SurfaceCoord s);
public:
    SurfaceScalarFunction(): f(nullptr){}
    SurfaceScalarFunction(double (*f0)(SurfaceCoord s)) : f(f0){}
    virtual double operator()(SurfaceCoord s)
    {
        if (f != nullptr)
            return (*f)(s);
        else
            return 0.0;
    }

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
    SurfaceCoord phiPlus = SurfaceCoord(s.theta, s.phi + stepsize / 2.0);
    SurfaceCoord phiMinus = SurfaceCoord(s.theta, s.phi - stepsize / 2.0);

    double dphi = 0;
    dphi = (*u)(phiPlus) - (*u)(phiMinus);
    dphi /= sin(s.theta) * stepsize;

    //add contribution to gradient.
    SphereCoord temp;
    temp = dphi * ephi;

    //evaulate partial derivative with respect to phi
    SurfaceCoord thetaPlus = SurfaceCoord(s.theta + stepsize / 2.0, s.phi);
    SurfaceCoord thetaMinus = SurfaceCoord(s.theta - stepsize / 2.0, s.phi);

    double dtheta = 0;
    dtheta = (*u)(thetaPlus) - (*u)(thetaMinus);
    dtheta /= stepsize;

    //add contribution to gradient and return.
    temp = temp + dtheta * etheta;
    return temp;

    

}

//Encapsulates a Spherical coordinate to Scalar function.
class SphericalScalarFunction
{
private:
    double (*f)(SphereCoord s);
public:
    SphericalScalarFunction() : f(nullptr) {}
    SphericalScalarFunction(double (*f0)(SphereCoord s)) : f(f0) {}
    virtual double operator()(SurfaceCoord s)
    {
        if (f != nullptr)
            return (*f)(s);
        else
            return 0.0;
    }

};


//Encapsulates a Spherical coordinate to Spherical Coordinate function.
class SphericalVectorField
{
private:
    SphereCoord (*f)(SphereCoord s);
public:
    SphericalVectorField() : f(nullptr) {}
    SphericalVectorField(SphereCoord (*f0)(SphereCoord s)) : f(f0) {}
    virtual SphereCoord operator()(SphereCoord s)
    {
        if (f != nullptr)
            return (*f)(s);
        else
            return s;
    }

};

class SurfaceGrad : public SphericalVectorField
{
private:
    SurfaceScalarFunction* arg;
    double step;
public:
    SurfaceGrad(SurfaceScalarFunction* t , double st = 1e-5): arg(t) , step(st){}

    SphereCoord operator()(SphereCoord r)
    {
        return surfaceGrad(arg, r.s, step);
    }
    
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

/*
template<class T>
concept SVF = std::is_base_of < SphericalVectorField, T >::value;

template<SVF T,SVF V>
class VFPLUS : public SphericalVectorField
{
private:
    T left;
    V right;
public:
    VFPLUS(T t, V v)
    {
        left = t;
        right = v;
    }

    SphereCoord operator()(SphereCoord r)
    {
        return left(r) + right(r);
    }
};

template <SVF T, SVF V>
VFPLUS<T, V> operator +(T t, V v)
{
    return VFPLUS<T, V>(t, v);
}

*/

//Computes the L2 inner product of two vector fields on the Surface of a unit sphere using a trapezoidal rule in both coordinates.
double L2InnerProduct(SphericalVectorField* f1, SphericalVectorField* f2, int n)
{

    const double GLnodes[16] = { -0.0950125098376374  , 0.0950125098376374  , -0.2816035507792589 , 0.2816035507792589 ,  -0.4580167776572274  , 0.4580167776572274 
                          - 0.6178762444026438  , 0.6178762444026438 , -0.7554044083550030 , 0.7554044083550030 , -0.8656312023878318 , 0.8656312023878318  ,
                          -0.9445750230732326 , 0.9445750230732326  , -0.9894009349916499  , 0.9894009349916499 };

    const double GLweights[16] = { 0.1894506104550685 , 0.1894506104550685 , 0.1826034150449236  , 0.1826034150449236  , 0.1691565193950025  , 0.1691565193950025 ,
                            0.1495959888165767  , 0.1495959888165767 , 0.1246289712555339  , 0.1246289712555339  , 0.0951585116824928  , 0.0951585116824928 ,
                            0.0622535239386479  , 0.0622535239386479 , 0.0271524594117541  , 0.0271524594117541 };
    double total = 0;

    for (double p = 0.0; p < 2 * PI; p += 2.0 * PI / n)
        for (int i = 0; i < 16; i++)
        {
            SurfaceCoord s( PI / 2.0 * (GLnodes[i] + 1), p);
            SphereCoord x(s);
            total += sin(s.theta) * GLweights[i] * dot((*f1)(x), (*f2)(x)) *  (PI / n) * PI;
        }

    return total;
}



//computes the L2 difference of two vector fields on the Surface of a unit sphere using a trapezoidal rule in both coordinates.
double L2Difference(SphericalVectorField* f1, SphericalVectorField* f2, int n)
{
    const double GLnodes[16] = { -0.0950125098376374  , 0.0950125098376374  , -0.2816035507792589 , 0.2816035507792589 ,  -0.4580167776572274  , 0.4580167776572274
                          - 0.6178762444026438  , 0.6178762444026438 , -0.7554044083550030 , 0.7554044083550030 , -0.8656312023878318 , 0.8656312023878318  ,
                          -0.9445750230732326 , 0.9445750230732326  , -0.9894009349916499  , 0.9894009349916499 };

    const double GLweights[16] = { 0.1894506104550685 , 0.1894506104550685 , 0.1826034150449236  , 0.1826034150449236  , 0.1691565193950025  , 0.1691565193950025 ,
                            0.1495959888165767  , 0.1495959888165767 , 0.1246289712555339  , 0.1246289712555339  , 0.0951585116824928  , 0.0951585116824928 ,
                            0.0622535239386479  , 0.0622535239386479 , 0.0271524594117541  , 0.0271524594117541 };
    double total = 0;

    for (double p = 0.0; p < 2 * PI; p += 2.0 * PI / n)
        for (int i = 0; i < 16; i++)
        {
            SurfaceCoord s(PI  / 2.0 * (GLnodes[i] + 1), p);
            SphereCoord x(s);
            SphereCoord diff = (*f1)(x) - (*f2)(x);
            total += sin(s.theta) * GLweights[i] * dot(diff, diff) * (PI / n) * PI;
        }

    return sqrt(total);
}

//computes the sup norm of two vector fields on the surface of a unit sphere.
double LInfDifference(SphericalVectorField* f1, SphericalVectorField* f2, int n)
{
    double sup = 0;
    double curr = 0;

    for (double t = 0.0; t < PI; t += PI / n)
        for (double p = 0.0; p < 2.0 * PI; p += 2.0 * PI / n)
        {
            SurfaceCoord s(t, p);
            SphereCoord x(s);
            RectCoord diff = RectCoord((*f1)(x)) - RectCoord((*f2)(x));
            curr = sqrt(dot(diff , diff));
            if (curr > sup)
                sup = curr;

        }

    return sup;
}