#pragma once


#include <cmath>
#include <complex>
#include <vector>
#include <memory>

const double PI = 4.0 * atan(1.0);

//structures for the different coordinate systems. Surface coordinates are for the unit sphere.
struct SurfaceCoord
{

    //azimuthal angle.
    double theta = 0;
    //zenith angle.
    double phi = 0;
    

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
        rho = 0;
    }

    SphereCoord(double r)
    {
        rho = r;
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

std::ostream & operator <<(std::ostream &  o, SphereCoord s)
{
    o << "( " << s.rho << ", " << s.s.phi << ", " << s.s.theta <<")";
    return o;
}
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
    return RectToSphere(a * SphereToRect(x));
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

//returns the norm squared of the given spherical coordinates.
double norm(SphereCoord x)
{
    return x.rho;
}


SphereCoord e_r(SphereCoord x)
{
    return SphereCoord(1 , x.s);
}

SphereCoord e_phi(SphereCoord x)
{
    return RectToSphere(RectCoord(-sin(x.s.phi), cos(x.s.phi), 0));
}

SphereCoord e_theta(SphereCoord x)
{
    return RectToSphere(RectCoord(cos(x.s.theta) * cos(x.s.phi), cos(x.s.theta) * sin(x.s.phi), -sin(x.s.theta)));
}

template<class DOM , class RAN>
class Functor
{
private:
    RAN (*f)(DOM s);
public:
    Functor() : f(nullptr) {}
    Functor(RAN (*f0)(DOM s)) : f(f0) {}
    virtual RAN operator()(DOM s)
    {
        if (f != nullptr)
            return (*f)(s);
        else
            return RAN();
    }

};



//Encapulation of a Surface to scalar function. In all such encapsulations, to create a function, create a derived class and implement the function operator().
//Then create a variable of the derived type and "call it". See code in main file for examples. In particular the rho function(immediately before main function)
// is the most simple.
/*
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
*/
typedef Functor<SurfaceCoord, double> SurfaceScalarFunction;
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
typedef Functor<SphereCoord, double> SphericalScalarFunction;

//Encapsulates a Spherical coordinate to Spherical Coordinate function.
typedef Functor<SphereCoord, SphereCoord> SphericalVectorField;

template<class DOM, class RAN>
class FunctorSum: public Functor<DOM , RAN>
{
private:

    int numterms;
    std::vector<Functor<DOM , RAN >*> fs;
public:

    FunctorSum()
    {
        numterms = 0;
        fs.clear();
    }
    void append(Functor<DOM , RAN>& f)
    {
        fs.push_back(&f);
        numterms++;
    }


    virtual RAN operator()(DOM x)
    {
        RAN temp = RAN(0.0);
        for (int i = 0; i < numterms; i++)
        {
            temp = temp + (*(fs[i]))(x);
        }

        return temp;
    }
};

typedef FunctorSum<SphereCoord, double>       SphereScalFunctionSum;

typedef FunctorSum<SphereCoord, SphereCoord>  VectorFieldSum;

//Serves as a container allowing a Sum of vector fields to be called as a single function. to add a new function to the sum,
// call the append function and pass an encapsulated function to it.


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
{   /*const double GLnodes[16] = {-0.0950125098376374  , 0.0950125098376374  , -0.2816035507792589 , 0.2816035507792589 ,  -0.4580167776572274  , 0.4580167776572274,
                          -0.6178762444026438  , 0.6178762444026438 , -0.7554044083550030 , 0.7554044083550030 , -0.8656312023878318 , 0.8656312023878318  ,
                          -0.9445750230732326 , 0.9445750230732326  , -0.9894009349916499  , 0.9894009349916499 };

    */
    const double GLnodes[32] = { -0.0483076656877383, 0.0483076656877383,
                                 -0.1444719615827965, 0.1444719615827965,
                                 -0.2392873622521371, 0.2392873622521371,
                                 -0.3318686022821277, 0.3318686022821277,
                                 -0.4213512761306353, 0.4213512761306353,
                                 -0.5068999089322294, 0.5068999089322294,
                                 -0.5877157572407623, 0.5877157572407623,
                                 -0.6630442669302152, 0.6630442669302152,
                                 -0.7321821187402897, 0.7321821187402897,
                                 -0.7944837959679424, 0.7944837959679424,
                                 -0.8493676137325700, 0.8493676137325700,
                                 -0.8963211557660521, 0.8963211557660521,
                                 -0.9349060759377397, 0.9349060759377397,
                                 -0.9647622555875064, 0.9647622555875064,
                                 -0.9856115115452684, 0.9856115115452684,
                                 -0.9972638618494816, 0.9972638618494816 };

    /*
    const double GLweights[16] = { 0.1894506104550685 , 0.1894506104550685 , 0.1826034150449236  , 0.1826034150449236  , 0.1691565193950025  , 0.1691565193950025 ,
                            0.1495959888165767  , 0.1495959888165767 , 0.1246289712555339  , 0.1246289712555339  , 0.0951585116824928  , 0.0951585116824928 ,
                            0.0622535239386479  , 0.0622535239386479 , 0.0271524594117541  , 0.0271524594117541 };
    */
    const double GLweights[32] = { 0.0965400885147278 , 0.0965400885147278,
                                   0.0956387200792749 , 0.0956387200792749,
                                   0.0938443990808046 , 0.0938443990808046,
                                   0.0911738786957639 , 0.0911738786957639,
                                   0.0876520930044038 , 0.0876520930044038,
                                   0.0833119242269467 , 0.0833119242269467,
                                   0.0781938957870703 , 0.0781938957870703,
                                   0.0723457941088485 , 0.0723457941088485,
                                   0.0658222227763618 , 0.0658222227763618,
                                   0.0586840934785355 , 0.0586840934785355,
                                   0.0509980592623762 , 0.0509980592623762,
                                   0.0428358980222267 , 0.0428358980222267,
                                   0.0342738629130214 , 0.0342738629130214,
                                   0.0253920653092621 , 0.0253920653092621,
                                   0.0162743947309057 , 0.0162743947309057,
                                   0.0070186100094701 , 0.0070186100094701 };

    double total = 0;

    for (int p = 0; p < n; p++)
        for (int i = 0; i < 32; i++)
        {
            SurfaceCoord s( PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)n);
            SphereCoord x(s);
            total += sin(s.theta) * GLweights[i] * dot((*f1)(x), (*f2)(x)) *  (PI / n) * PI;
        }

    return total;
}



//computes the L2 difference of two vector fields on the Surface of a unit sphere using a trapezoidal rule in both coordinates.
double L2Difference(SphericalVectorField* f1, SphericalVectorField* f2, int n , double radius = 1)
{
    /*const double GLnodes[16] = {-0.0950125098376374  , 0.0950125098376374  , -0.2816035507792589 , 0.2816035507792589 ,  -0.4580167776572274  , 0.4580167776572274,
                          -0.6178762444026438  , 0.6178762444026438 , -0.7554044083550030 , 0.7554044083550030 , -0.8656312023878318 , 0.8656312023878318  ,
                          -0.9445750230732326 , 0.9445750230732326  , -0.9894009349916499  , 0.9894009349916499 };

    */
    const double GLnodes[32] = { -0.0483076656877383, 0.0483076656877383,
                                 -0.1444719615827965, 0.1444719615827965,
                                 -0.2392873622521371, 0.2392873622521371,
                                 -0.3318686022821277, 0.3318686022821277,
                                 -0.4213512761306353, 0.4213512761306353,
                                 -0.5068999089322294, 0.5068999089322294,
                                 -0.5877157572407623, 0.5877157572407623,
                                 -0.6630442669302152, 0.6630442669302152,
                                 -0.7321821187402897, 0.7321821187402897,
                                 -0.7944837959679424, 0.7944837959679424,
                                 -0.8493676137325700, 0.8493676137325700,
                                 -0.8963211557660521, 0.8963211557660521,
                                 -0.9349060759377397, 0.9349060759377397,
                                 -0.9647622555875064, 0.9647622555875064,
                                 -0.9856115115452684, 0.9856115115452684,
                                 -0.9972638618494816, 0.9972638618494816 };

    /*
    const double GLweights[16] = { 0.1894506104550685 , 0.1894506104550685 , 0.1826034150449236  , 0.1826034150449236  , 0.1691565193950025  , 0.1691565193950025 ,
                            0.1495959888165767  , 0.1495959888165767 , 0.1246289712555339  , 0.1246289712555339  , 0.0951585116824928  , 0.0951585116824928 ,
                            0.0622535239386479  , 0.0622535239386479 , 0.0271524594117541  , 0.0271524594117541 };
    */
    const double GLweights[32] = { 0.0965400885147278 , 0.0965400885147278,
                                   0.0956387200792749 , 0.0956387200792749,
                                   0.0938443990808046 , 0.0938443990808046,
                                   0.0911738786957639 , 0.0911738786957639,
                                   0.0876520930044038 , 0.0876520930044038,
                                   0.0833119242269467 , 0.0833119242269467,
                                   0.0781938957870703 , 0.0781938957870703,
                                   0.0723457941088485 , 0.0723457941088485,
                                   0.0658222227763618 , 0.0658222227763618,
                                   0.0586840934785355 , 0.0586840934785355,
                                   0.0509980592623762 , 0.0509980592623762,
                                   0.0428358980222267 , 0.0428358980222267,
                                   0.0342738629130214 , 0.0342738629130214,
                                   0.0253920653092621 , 0.0253920653092621,
                                   0.0162743947309057 , 0.0162743947309057,
                                   0.0070186100094701 , 0.0070186100094701 };

    double total = 0;

    for (int p = 0; p < n; p++)
        for (int i = 0; i < 32; i++)
        {
            SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)n);
            SphereCoord x(radius , s);
            SphereCoord diff = (*f1)(x) - (*f2)(x);
            total += radius * radius * sin(s.theta) * GLweights[i] * dot(diff, diff) * (PI / n) * PI;
            /*if (dot(diff, diff) > 1e-12)
            {
                std::cout << " error of " << dot(diff, diff) << " at " << x << std::endl;
                std::cout << "value of f: " << (*f1)(x) << std::endl;
                std::cout << "value of g: " << (*f2)(x) << std::endl;

            }
            */
        }

    return sqrt(total);
}

double max(RectCoord x)
{
    double temp;

    if (x.x > x.y)
        temp = x.x;
    else
        temp = x.y;

    if (temp > x.z)
        return temp;
    else
        return x.z;
}

//computes the sup norm of two vector fields on the surface of a unit sphere.
double LInfDifference(SphericalVectorField* f1, SphericalVectorField* f2, int n)
{
    const double GLnodes[32] = { -0.0483076656877383, 0.0483076656877383,
                                -0.1444719615827965, 0.1444719615827965,
                                -0.2392873622521371, 0.2392873622521371,
                                -0.3318686022821277, 0.3318686022821277,
                                -0.4213512761306353, 0.4213512761306353,
                                -0.5068999089322294, 0.5068999089322294,
                                -0.5877157572407623, 0.5877157572407623,
                                -0.6630442669302152, 0.6630442669302152,
                                -0.7321821187402897, 0.7321821187402897,
                                -0.7944837959679424, 0.7944837959679424,
                                -0.8493676137325700, 0.8493676137325700,
                                -0.8963211557660521, 0.8963211557660521,
                                -0.9349060759377397, 0.9349060759377397,
                                -0.9647622555875064, 0.9647622555875064,
                                -0.9856115115452684, 0.9856115115452684,
                                -0.9972638618494816, 0.9972638618494816 };

    double sup = 0;
    double curr = 0;

    for (int i = 0; i < n; i++)
        for (double j = 0; j < 32; j++)
        {
            SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)j / (double)n);
            SphereCoord x(s);
            SphereCoord diff = (*f1)(x) - (*f2)(x);
            curr = max(RectCoord(diff));
            if (curr > sup)
                sup = curr;
            
            if (curr > 1e-10)
            {
                std::cout << "Error of" << curr << " at " << x << "\n";
            }
        }

    return sup;
}
