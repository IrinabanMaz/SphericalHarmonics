/**
*  @file SphericalCalc.h
* 
* @brief Provides basic coordinate systems and a frame for using functions on those coordinate systems as data. 
* 
* We introduce here the Rectangular and Spherical domains and their arithmetic operations. With these domains we define Functors acting on them.
* General numerical methods for integrating and differentiating functions are provided. In addition, functions can be used in a discrete
* or continuous manner. We state at this point that all quadratures are to use the same method. To compute an integral of the form 
*  @f{equation}{
* \int_\Gamma f(x) dS_x = \int_0^\pi \int_0^{2 \pi} f(\theta , \phi) r^2 \sin\theta \ d\phi d\theta  
@f}
*  We use a Gauss Legendre Quadrature for the @f$\theta @f$ integral and trapezoidal in the @f$ \phi @f$ integral.
*/



#pragma once


#include <cmath>
#include <complex>
#include <vector>
#include <memory>
#include <string>



const int MAXTRAPNODES = 64;
const int MAXGLNODES = 64;

double GLnodes[MAXGLNODES] = {-0.0950125098376374  , 0.0950125098376374  , -0.2816035507792589 , 0.2816035507792589 ,  -0.4580167776572274  , 0.4580167776572274,
                          -0.6178762444026438  , 0.6178762444026438 , -0.7554044083550030 , 0.7554044083550030 , -0.8656312023878318 , 0.8656312023878318  ,
                          -0.9445750230732326 , 0.9445750230732326  , -0.9894009349916499  , 0.9894009349916499 };

    

/*double GLnodes[MAXGLNODES] = {-0.0483076656877383, 0.0483076656877383,
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
 */                            

double GLweights[MAXGLNODES] = { 0.1894506104550685 , 0.1894506104550685 , 0.1826034150449236  , 0.1826034150449236  , 0.1691565193950025  , 0.1691565193950025 ,
                        0.1495959888165767  , 0.1495959888165767 , 0.1246289712555339  , 0.1246289712555339  , 0.0951585116824928  , 0.0951585116824928 ,
                        0.0622535239386479  , 0.0622535239386479 , 0.0271524594117541  , 0.0271524594117541 };
                        

/*
double GLweights[MAXGLNODES] = { 0.0965400885147278 , 0.0965400885147278,
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
  */                             

double NUMGLNODES = 16;
double NUMTRAPNODES = 16;

const double PI = 4.0 * atan(1.0);

/**
 * @brief Struct representing a point on the surface of a sphere in spherical coordinates.
 * 
*/
struct SurfaceCoord
{

    double theta = 0; /** < The zenith angle of the sphere. Takes values @f 0 \leq \theta \leq \pi @f */
    double phi = 0; /** < The Azimuthal angle of the sphere. Takes values @f 0 \leq \theta \leq  2 \pi @f */
    

    /**
     * @brief Addition for surface coordinte(unused).
     * @param s2 The second point.
     * @return A SurfaceCoord whose zenith and azimuth are the sums of the passed arguments, modulo the sphere's surface. 
    */
    SurfaceCoord operator +(const SurfaceCoord & s2) const
    {
        SurfaceCoord temp;
        temp.phi = phi + s2.phi;
        temp.theta = theta + s2.theta;


        while (temp.theta > 2.0 * PI)
            temp.theta -= 2.0 * PI;

        if (temp.theta > PI)
        {
            temp.theta = 2.0 * PI - temp.theta;
            temp.phi += PI;
        }

        while (temp.phi > 2.0 * PI)
            temp.phi -= 2.0 * PI;
            

        return temp;
    }

    /**
     * @brief Subtraction for surface coordinte(unused).
     * @param s2 The second point.
     * @return A SurfaceCoord whose zenith and azimuth are the sums of the passed arguments, modulo the sphere's surface.
    */
    SurfaceCoord operator -(const SurfaceCoord & s2) const
    {
        SurfaceCoord temp;
        temp.phi = phi - s2.phi;
        temp.theta = theta - s2.theta;


        while (temp.theta < - PI)
            temp.theta += 2.0 * PI;

        if (temp.theta < 0.0)
        {
            temp.theta = - temp.theta;
            temp.phi += PI;
        }

        while (temp.phi > 2.0 * PI)
            temp.phi -= 2.0 * PI;

        while (temp.phi < 0.0)
            temp.phi += 2.0 * PI;

        return temp;
    }

    SurfaceCoord():
        theta(0),
        phi(0){}
    

    SurfaceCoord(double t, double p):
       theta(t),
        phi(p){}

    
   
};



/**
 * @brief A point in  @f$ R^3 @f$ in the eulerian description.
*/
struct RectCoord
{
    double x; /**< the first coordinate of the point. */
    double y; /**< the second coordinate of the point. */
    double z; /**< the third coordinate of the point. */

    /**
     * @brief Constructor for RectoCoord.
     * @param x0 Value to be assigned to the first coordinate.
     * @param y0 Value to be assigned to the second coordinate.
     * @param z0 Value to be assigned to the third coordinate.
    */
    RectCoord(double x0 = 0.0, double y0 = 0.0, double z0 = 0.0):
        x(x0),
        y(y0),
        z(z0){}

    
    RectCoord& operator =(const double& a)
    {
        x = a;
        y = a;
        z = a;
        return *this;
    }
    
};

/**
 * @brief Addition for two RectCoord points.
*/
RectCoord operator +(const RectCoord& a,const RectCoord & b)
{
    return RectCoord(a.x + b.x, a.y + b.y, a.z + b.z);
}

/**
 * @brief Subtraction for two RectCoord points.
*/
RectCoord operator -(const RectCoord & a,const RectCoord& b)
{
    return RectCoord(a.x - b.x, a.y - b.y, a.z - b.z);
}

RectCoord operator -(const RectCoord& a)
{
    return RectCoord(-a.x, -a.y, -a.z);
}

/**
 * @brief Scalar multiplication for RectCoord.
*/
RectCoord operator *(double c, RectCoord a)
{
    return RectCoord(c * a.x, c * a.y, c * a.z);
}
/**
 * @brief Scalar multiplication for RectCoord.
*/
RectCoord operator *(RectCoord a , double c)
{
    return RectCoord(c * a.x, c * a.y, c * a.z);
}

RectCoord operator / (RectCoord x, double a)
{
    return RectCoord(x.x / a, x.y / a, x.z / a);
}
RectCoord operator +=(RectCoord & x, const RectCoord& y)
{
    x = x + y;
    return x;
}
RectCoord operator -=(RectCoord& x,const RectCoord& y)
{
    x = x - y;
    return x;
}

RectCoord operator *=(RectCoord & x,const double &  a)
{
    x = x * a;
    return x;
}

RectCoord operator /=(RectCoord & x, const double & a)
{
    x = x / a;
    return x;
}



/**
 * @brief dot product for rectangular coordinates.
*/
double dot(RectCoord a, RectCoord b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/**
 * @brief Vector cross product.
*/
RectCoord cross(RectCoord a, RectCoord b)
{
    RectCoord temp;
    temp.x = a.y * b.z - a.z * b.y;
    temp.y = a.z * b.x - a.x * b.z;
    temp.z = a.x * b.y - a.y * b.x;

    return temp;
}

/**
 * @brief Length of a RectCoord.
*/
double norm(const RectCoord & x)
{
    return sqrt(dot(x, x));
}

/**
 * @brief string output for a RectCoord.
 * 
 * Output will be of the form [x, y, z].
*/
std::ostream& operator <<(std::ostream& o, RectCoord x)
{
    o << "[" << x.x << ", " << x.y << ", " << x.z << "]";
    return o;
}
bool operator ==(RectCoord x, RectCoord y)
{
    return norm(x - y) <= 1e-15;
}
/**
 * @brief Rotates the vector v about the given axis with given angle.
 * @param v The vector to rotate.
 * @param axis The axis to revolve about.
 * @param angle The angle between input and resulting vector.
 * @return The rotated angle.
 * 
 * The rotation is computed via Rodriguez formula, where @f$ k @f$ is the axis and @f$\theta @f$ is the angle.
 * @f$ v_\mbox{rot} = v \cos\theta + (k \times v)\sin\theta + k(k\cdot v)(1 - \cos\theta)@f$.
*/
RectCoord rotate(RectCoord v, RectCoord axis, double angle)
{
    return v * cos(angle) + cross(axis, v) * sin(angle) + axis * dot(axis, v)*(1.0 - cos(angle));
}

/**
 * @brief A point in a spherical coordinate system, centered at a point with a given radius.
*/
struct SphereCoord
{
    double rho;/**< The radius of the point in the spherical coordinate system. */
    double radius; /**< the radius of the coordinate system itself.*/
    RectCoord center;/**< The center of the coordinate system. */
    SurfaceCoord s;/**<The surface data for the point.*/
    SurfaceCoord northPole; /** The point on the spherical surface corresponding to @f$\theta = \phi = 0 @f$/

    /**
     * @brief Contructor for the spherical coordinate. 
     * @param rh Radius of the point.
     * @param s1 surface data for the point.
     * @param c center of the coordinate system, optional and defaults to origin.
     * @param r radius of the coordinate system, optional and defaults to 1.
    */
    SphereCoord(double rh = 0.0, SurfaceCoord s1  = SurfaceCoord(), RectCoord c = RectCoord(),SurfaceCoord axis = SurfaceCoord(), double r = 1.0)
    {
        rho = rh;
        s = s1;
        center = c;
        radius = r;
        northPole = axis;
    }
    /**
     * @brief Constructor with unit point radius.
    */
    SphereCoord(SurfaceCoord s1, RectCoord c = RectCoord(), SurfaceCoord axis = SurfaceCoord(), double r = 1.0)
    {
        rho = 1;
        s = s1;
        center = c;
        radius = r;
        northPole = axis;
    }
};


/**
 * @brief String output for the SphereCoord.
 * 
 * Output will be of the form @f$(\rho, \theta, \phi)@f$. Information about the coordinate system is not printed.
*/
std::ostream & operator <<(std::ostream &  o, SphereCoord s)
{
    o << "( " << s.rho << ", " << s.s.theta << ", " << s.s.phi <<")";
    return o;
}
/**
 * @brief Converts a RectCoord into a SphereCoord with specified system parameters.
 * @param x the point to convert.
 * @param c the center of the spherical coordinate system.
 * @param axis The north pole of the coordinate system.
 * @param radius the radius of the coordinate system.
 * @return the SphereCoord corresponging to the RectCoord in the given spherical coordinate system.
 * 
 * We use the conversion formulas:
 * 
 * @f{eqnarray}{
  \rho &=&  \frac{\|x - c\|}{\mbox{radius}}\\
  \phi &=& \mbox{atan2}(x_2 - c_2 , x_1 - c_1) \\
  \phi &=& \mbox{atan2}(\sqrt{(x_1 - c_1)^2 + (x_2 - c_2)^2}, x_3 - c_3) 
  @f}
*/
SphereCoord RectToSphere(RectCoord x , RectCoord c = RectCoord(), SurfaceCoord axis = SurfaceCoord(), double radius = 1.0)
{
    SphereCoord temp;

    RectCoord ephi(-sin(axis.phi), cos(axis.phi), 0);

    RectCoord xr = rotate(x - c, ephi, -axis.theta);
    
    temp.rho = sqrt(dot(xr, xr)) / radius;
    
    temp.s.phi = atan2(xr.y, xr.x);
    
    double planenorm = sqrt((xr.x) * (xr.x) + (xr.y) * (xr.y));

    temp.s.theta = atan2(planenorm, xr.z);
    
    temp.center = c;
    temp.northPole = axis;
    temp.radius = radius;

    return temp;

}
/**
 * @brief Converts a SphereCoord to a RectCoord.
 *
 * We use the surface conversion formulas:
 * 
 *  @f{eqnarray}{
  s_1 &=&  \rho \sin(\theta) \cos(\phi)\\
  s_2 &=& \rho \sin(\theta) \sin(\phi) \\
  s_3 &=& \rho \cos(\theta) 
  @f}
*
*  Then apply the scaling and translation to account for the properties of the 
*  given coordinate system to yield @f$ \mbox{radius} \cdot s + \mbox{center}@f$
*/
RectCoord SphereToRect(SphereCoord r)
{
    RectCoord temp;
    temp.x = r.rho * sin(r.s.theta) * cos(r.s.phi);
    temp.y = r.rho * sin(r.s.theta) * sin(r.s.phi);
    temp.z = r.rho * cos(r.s.theta);

    RectCoord ephi(-sin(r.northPole.phi), cos(r.northPole.phi), 0);
    temp = rotate(temp, ephi, r.northPole.theta);

    temp = r.radius * temp + r.center;

    return temp;
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

/**
 * @brief Converts a SphereCoord to a SphereCoord corresponding to the same point, in a different spherical coordinate system.
 * @param x The point to convert.
 * @param c The new center.
 * @param r the new radius.
 * @return The converted point.
*/
SphereCoord recenter(SphereCoord x, RectCoord c = RectCoord(),SurfaceCoord axis = SurfaceCoord(), double r = 1.0)
{
    return RectToSphere(SphereToRect(x), c,axis ,r);
}

/**
 * @brief Adds two points in spherical coordinate systems. The sum is equivalent to summing their RectCoord equivalents.
*/
RectCoord operator +(const SphereCoord & r, const SphereCoord& q)
{
    return SphereToRect(r) + SphereToRect(q) ;
}

/**
 * @brief Subtracts two points in spherical coordinate systems. The difference is equivalent to subtracting their RectCoord equivalents.
*/
RectCoord operator -(const SphereCoord & r, const SphereCoord & q)
{
    return SphereToRect(r) - SphereToRect(q);
}

/**
 * @brief Subtracts two pointS. The difference is equivalent to subtracting their RectCoord equivalents.
*/
RectCoord operator -(const SphereCoord& r, const RectCoord& q)
{
    return SphereToRect(r) - q;
}

/**
 * @brief Subtracts two points. The difference is equivalent to subtracting their RectCoord equivalents.
*/
RectCoord operator -(const RectCoord& r, const SphereCoord& q)
{
    return r - SphereToRect(q);
}

/**
 * @brief Adds a RectCoord to a SphereCoord. The sum is equivalent to summing their RectCoord equivalents.
*/
RectCoord operator +(const RectCoord & x,const SphereCoord & s)
{
    return x + SphereToRect(s);
}
/**
 * @brief Adds a RectCoord to a SphereCoord. The sum is equivalent to summing their RectCoord equivalents.
*/
RectCoord operator +(const SphereCoord & s, const RectCoord & x )
{
    return x + SphereToRect(s);
}

/**
 * @brief Scales a SphereCoord.
*/
SphereCoord operator *(const double & a, const SphereCoord & x)
{
    return RectToSphere(a * SphereToRect(x) , x.center , x.northPole,x.radius);
}

/**
 * @brief Scales a SphereCoord.
*/
SphereCoord operator *(const SphereCoord & x,const double & a)
{
    return RectToSphere(a * SphereToRect(x) , x.center , x.northPole,x.radius);
}

/**
 * @brief Computes the dot product of two SphereCoord objects. The sum is equivalent to summing their RectCoord equivalents.
*/
double dot(const SphereCoord & x,const SphereCoord & y)
{
    return dot(SphereToRect(x), SphereToRect(y));
}

/**
 * @brief Computes the cross product of two SphereCoord objects. The sum is equivalent to summing their RectCoord equivalents.
*/
RectCoord cross(const SphereCoord & x,const SphereCoord & y)
{
    return cross(SphereToRect(x), SphereToRect(y));
}

/**
 * @brief computes the norm of a SphereCoord in the eulerian description.
*/
double norm(const SphereCoord & x)
{
    RectCoord temp = SphereToRect(x);
    return sqrt(dot(temp , temp));
}

/**
 * @brief Computes the Spherical Basis vector @f$ e_r @f$ corresponding to a point.
 * 
 * We ue the formula @f$ e_r(x) = \frac{x}{ r} @f$
*/
RectCoord e_r(const SphereCoord & x)
{
    RectCoord r =x - x.center;
    return r / norm(r);
}


/**
 * @brief Computes the Spherical Basis vector @f$ e_r @f$ corresponding to a point.
 *
 * We ue the formula @f$ e_\phi(x) = [ -\sin\phi , \cos\phi , 0] @f$
*/
RectCoord e_phi(const SphereCoord & x)
{
    RectCoord temp = RectCoord(-sin(x.s.phi), cos(x.s.phi), 0);
    RectCoord ephiPole = RectCoord(-sin(x.northPole.phi), cos(x.northPole.phi), 0.0);
    return rotate(temp , ephiPole , x.northPole.theta);
}


/**
 * @brief Computes the Spherical Basis vector @f$ e_r @f$ corresponding to a point.
 *
 * We ue the formula @f$ e_\phi(x) = [ \cos\theta \cos\phi , \cos\theta \sin\phi , -\sin\theta] @f$
*/
RectCoord e_theta(const SphereCoord & x)
{
    RectCoord temp = RectCoord(cos(x.s.theta) * cos(x.s.phi), cos(x.s.theta) * sin(x.s.phi), -sin(x.s.theta));
    RectCoord ephiPole = RectCoord(-sin(x.northPole.phi), cos(x.northPole.phi), 0.0);
    return rotate(temp, ephiPole, x.northPole.theta);
}



/**
 * @brief The Functor template class. allows us to create objects which behave syntactically like functions, but can be passed as data.
 * @tparam DOM The domain type of the functor.
 * @tparam RAN The range type of the functor.
 * 
 * The mathematical functions we will work with are derived classes of Functors. We take advantage of polymorphism mostly for testing purposes.
*/
template<class DOM , class RAN>
class Functor
{
private:
    RAN (*f)(const DOM & s);
public:
    std::string name; /**< the name of the function.*/
    /**
     * @brief Default cuntructor. creates a constant functor which returns the default value of RAN.
    */
    Functor() : f(nullptr) {}

    /**
     * @brief Function pointer constructor. The functor will call the passed function and take on it's behavior.
     * @param f0 The function to emulate.
     * @param n The name.
    */
    Functor(RAN (*f0)(const DOM& s) , std::string n = "") : f(f0) , name(n) {}

    /**
     * @brief The evaluation method for the Functor.
     * @param s the domain element which the functor will evaluate at.
     * @return In base class, the pointer passed on constructor. 
     * 
     * This function will be redefined to have whatever behavior we choose through inheritence.
    */
    virtual  RAN operator()(const DOM & s)const  
    {
        if (f != nullptr)
            return (*f)(s);
        else
            return RAN();
    }

};




/**
 * @brief Type corresponding to functions which map a point on a spherical surface to a real value.
*/
typedef Functor<SurfaceCoord, double> SurfaceScalarFunction;

/**
 * @brief The numerical surface gradient of a function evaluated at a point.
 * @param u the function to compute the gradient of.
 * @param s the point at which we evaluate the gradient.
 * @param stepsize the numerical parameter for computing partial derivates using finite differences.
 * @note the numerical surface gradient is only used for testing. Exact relations for surface gradients are used in solving the mobility problem.
 * @attention The Surface gradient is singular at the poles @f$ \theta = 0 @f$ and @f$ \theta = \pi @f$.
 *
 *The surface gradient is defined by the formula: 
  @f{equation}{
  \nabla_\Gamma u = \frac{\partial u}{\partial \theta} e_\theta + \frac{1}{\sin\theta} \frac{\partial u}{\partial\phi}e_\phi
  @f}
  *
  * and the partial derivatives are evaulated using the finite difference, where @f$ h = @f$ stepsize:
  @f{equation}{
  \frac{\partial u}{\partial \theta} = \frac{u(\theta + h , \phi) - u(\theta - h , \phi)}{2h}
  @f}
  * and the same numerical scheme is used for @f$ \frac{\partial}{\partial \phi}. @f$
*/
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
    temp = RectToSphere(temp + dtheta * etheta);
    return temp;

    

}

/**
 * @brief Functor where the input is a Spherical Coordinate and the output is real.
*/
typedef Functor<SphereCoord, double> SphericalScalarFunction;

/**
 * @brief Functor where the input is a SphereCoord and the output is a RectCoord.
*/
typedef Functor<SphereCoord, RectCoord> SphericalVectorField;


/**
 * @brief A special type of functor which evaulates the sum of functions accumulated through an FunctorSum::append() function.
 * @tparam DOM The common domain of the functions to sum.
 * @tparam RAN The common range of the functions to sum. Ran must have a constructor RAN(double ).
*/
template<class DOM, class RAN>
class FunctorSum: public Functor<DOM , RAN>
{
private:


    FunctorSum(const FunctorSum<DOM , RAN> & f){}
    FunctorSum(const FunctorSum<DOM , RAN> && f) noexcept{}
    FunctorSum & operator =(const FunctorSum<DOM , RAN> & f){}
    FunctorSum & operator =(const FunctorSum<DOM , RAN> && f) noexcept {}

    int numterms;
    std::vector<Functor<DOM,RAN>*> fs;
public:

    /**
     * @brief Constructs an empty sum.
    */
    FunctorSum()
    {
        numterms = 0;
        fs.clear();
    }
    /**
     * @brief pushes the function f onto the sum.
    */
    void append(Functor<DOM , RAN>* f)
    {
        fs.push_back(f);
        numterms++;
    }

    void clear()
    {
        fs.clear();
        numterms = 0;
    }

    /**
     * @brief Evaluates the sum of all functions added via the append() method.
     * @param x 
     * @return 
    */
    virtual RAN operator()(const DOM & x) const
    {
        RAN temp = RAN(0.0);
        for (auto f : fs )
        {
            temp = temp + (*f)(x);
        }

        return temp;
    }
};

typedef FunctorSum<SphereCoord, double>       SphereScalFunctionSum;

typedef FunctorSum<SphereCoord, RectCoord>  VectorFieldSum;


/**
 * @brief Computes the @f$ L^2(\Gamma) @f$ difference of the functions @f$ f_1 @f$ and @f$ f_2 @f$ over an arbitrary spherical surface.
 * @param radius The radius of the sphere to integrate over.
 * @param center  The center of the sphere to integrate over.
 * @param debugOut Determines whether messages will print when a pointwise large error occurs.
 * @return The integral @f$\left(\int_\Gamma |f_1(x) - f_2(x)|^2 dS_x\right)^{1/2} @f$
*/
double L2Difference(SphericalScalarFunction* f1, SphericalScalarFunction* f2, double radius = 1.0 , RectCoord center = RectCoord() , bool debugOut = false)
{
   
    double total = 0;
for (int i = 0; i < NUMGLNODES; i++)
    for (int p = 0; p < NUMTRAPNODES; p++) 
        {
            SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
            SphereCoord x(radius, s , center);
            double f = (*f1)(x);
            double g = (*f2)(x);
            double diff =  f - g ;
            total += radius * radius * sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
            if ((abs(diff) > 1e-4)&&debugOut)
            {
                std::cout << "Large error between " << f1->name << " and " << f2->name << std::endl;
                std::cout << "absolute error of " << diff << " at " << x << std::endl;
                std::cout << "ratio of " << f / g << std::endl;
                std::cout << "value of f: " << f << std::endl;
                std::cout << "value of g: " << g << std::endl;

            }
            
        }

    return sqrt(total);
}

/**
 * @brief Computes the @f$ L^2(\Gamma) @f$ inner productof the functions @f$ f_1 @f$ and @f$ f_2 @f$ over an arbitrary spherical surface.
 * @param center  The center of the sphere to integrate over.
 * @param debugOut Determines whether messages will print when a pointwise large error occurs.
 * @return The integral @f$\int_\Gamma f_1(x) \cdot f_2(x) dS_x @f$
*/
double L2InnerProduct(SphericalVectorField* f1, SphericalVectorField* f2,  RectCoord c = RectCoord())
{  

    double total = 0;

    for (int p = 0; p < NUMTRAPNODES; p++)
        for (int i = 0; i < NUMGLNODES; i++)
        {
            SurfaceCoord s( PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
            SphereCoord x(s , c);
            total += sin(s.theta) * GLweights[i] * dot((*f1)(x), (*f2)(x)) *  (PI / NUMTRAPNODES) * PI;
        }

    return total;
}
/**
 * @brief Computes the @f$ L^2(\Gamma) @f$ inner productof the functions @f$ f_1 @f$ and @f$ f_2 @f$ over an arbitrary spherical surface.
 * @param center  The center of the sphere to integrate over.
 * @param debugOut Determines whether messages will print when a pointwise large error occurs.
 * @return The integral @f$\int_\Gamma f_1(x) \cdot f_2(x) dS_x @f$
*/
double L2InnerProduct(SphericalScalarFunction* f1, SphericalScalarFunction* f2, RectCoord c = RectCoord())
{

    double total = 0;

    for (int p = 0; p < NUMTRAPNODES; p++)
        for (int i = 0; i < NUMGLNODES; i++)
        {
            SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
            SphereCoord x(s, c);
            total += sin(s.theta) * GLweights[i] * (*f1)(x) * (*f2)(x) * (PI / NUMTRAPNODES) * PI;
        }

    return total;
}



/**
 * @brief Vector storing the discretization of a SphericalVectorField.
*/
typedef std::vector<std::vector<RectCoord>> SphereData;


SphereData operator +(const SphereData& a, const SphereData& b)
{
    SphereData temp;
    temp.resize(NUMGLNODES);
    for (int t = 0; t < NUMGLNODES; t++)
    {
        temp[t].resize(NUMTRAPNODES);
        for (int p = 0; p < NUMTRAPNODES; p++)
            temp[t][p] = a[t][p] + b[t][p];
    }

    return temp;
}

SphereData operator -(const SphereData& a, const SphereData& b)
{
    SphereData temp;
    temp.resize(NUMGLNODES);
    for (int t = 0; t < NUMGLNODES; t++)
    {
        temp[t].resize(NUMTRAPNODES);
        for (int p = 0; p < NUMTRAPNODES; p++)
            temp[t][p] = a[t][p] - b[t][p];
    }

    return temp;
}

/**
 * @brief Vector storing the discretization of a SphericalScalarFunction.
*/
typedef std::vector<std::vector<double>> ScalSphereData;

/**
 * @brief Computes the L2 inner product of a function with a discrete set of data located on the quadrature nodes.
 * @param c center of the sphere to integrate over.
 * @return The integral @f$\int_\Gamma data(x) \cdot f_2(x) dS_x @f$
*/
double L2InnerProductDiscrete(SphereData& data, SphericalVectorField* f2, RectCoord c = RectCoord())
{

    double total = 0;

    for (int p = 0; p < NUMTRAPNODES; p++)
        for (int i = 0; i < NUMGLNODES; i++)
        {
            SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
            SphereCoord x(s, c);
            total += sin(s.theta) * GLweights[i] * dot(data[i][p], (*f2)(x)) * (PI / NUMTRAPNODES) * PI;
        }

    return total;
}
/**
 * @brief Computes the L2 inner product of a function with a discrete set of data located on the quadrature nodes.
 * @param c center of the sphere to integrate over.
 * @return The integral @f$\int_\Gamma data1(x) \cdot data2(x) dS_x @f$
*/
double L2InnerProductDiscrete(SphereData& data1, SphereData& data2, RectCoord c = RectCoord())
{

    double total = 0;

    for (int p = 0; p < NUMTRAPNODES; p++)
        for (int i = 0; i < NUMGLNODES; i++)
        {
            SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
            SphereCoord x(s, c);
            total += sin(s.theta) * GLweights[i] * dot(data1[i][p], data2[i][p]) * (PI / NUMTRAPNODES) * PI;
        }

    return total;
}

/**
 * @brief Computes the @f$ L^2(\Gamma) @f$ difference of the functions @f$ f_1 @f$ and @f$ f_2 @f$ over an arbitrary spherical surface.
 * @param radius The radius of the sphere to integrate over.
 * @param center  The center of the sphere to integrate over.
 * @param debugOut Determines whether messages will print when a pointwise large error occurs.
 * @return The integral @f$\left(\int_\Gamma |f_1(x) - f_2(x)|^2 dS_x\right)^{1/2} @f$
*/
double L2Difference(SphericalVectorField* f1, SphericalVectorField* f2,  double radius = 1, RectCoord center = RectCoord(), bool debugOut = false)
{
    
    double total = 0;

    for (int p = 0; p < NUMTRAPNODES; p++)
        for (int i = 0; i < NUMGLNODES; i++)
        {
            SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
            SphereCoord x(radius , s , center);
            RectCoord diff = (*f1)(x) - (*f2)(x);
            total += radius * radius * sin(s.theta) * GLweights[i] * dot(diff, diff) * (PI / NUMTRAPNODES) * PI;
            if ((dot(diff, diff) > 1e-4 )|| isnan(dot(diff,diff))&& debugOut)
            {
                std::cout << " error of " << dot(diff, diff) << " at " << x << std::endl;
                std::cout << "x ratio of " << (*f1)(x).x / (*f2)(x).x << "\n";
                std::cout << "y ratio of " << (*f1)(x).y / (*f2)(x).y << "\n";
                std::cout << "z ratio of " << (*f1)(x).z / (*f2)(x).z << "\n";
                std::cout << "value of f: " << (*f1)(x) << std::endl;
                std::cout << "value of g: " << (*f2)(x) << std::endl;

            }
            
        }

    return sqrt(total);
}

/**
 * @brief Computes the integral of the function @f$ f @f$  over an arbitrary spherical surface.
 * @param radius The radius of the sphere to integrate over.
 * @param center  The center of the sphere to integrate over.
 * @param debugOut Determines whether messages will print when a pointwise large error occurs.
 * @return The integral @f$\int_\Gamma f(x) dS_x @f$
*/
RectCoord Integrate(SphericalVectorField* f, double radius = 1.0, RectCoord center = RectCoord())
{

    RectCoord total = 0.0;

    for (int p = 0; p < NUMTRAPNODES; p++)
        for (int i = 0; i < NUMGLNODES; i++)
        {
            SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
            SphereCoord x(radius, s, center);
            total =  total + radius * radius * sin(s.theta) * GLweights[i] * ((*f)(x)) * (PI / NUMTRAPNODES) * PI;
            

        }

    return total;
}

/**
 * @brief Computes the integral of the function @f$ f @f$  over an arbitrary spherical surface.
 * @param radius The radius of the sphere to integrate over.
 * @param center  The center of the sphere to integrate over.
 * @param debugOut Determines whether messages will print when a pointwise large error occurs.
 * @return The integral @f$\int_\Gamma f(x) dS_x @f$
*/
double Integrate(SphericalScalarFunction* f, double radius = 1.0, RectCoord center = RectCoord())
{

    double total = 0.0;

    for (int p = 0; p < NUMTRAPNODES; p++)
        for (int i = 0; i < NUMGLNODES; i++)
        {
            SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
            SphereCoord x(radius, s, center);
            total = total + radius * radius * sin(s.theta) * GLweights[i] * ((*f)(x)) * (PI / NUMTRAPNODES) * PI;


        }

    return total;
}

/**
 * @brief Computes the sum of the data over an arbitrary spherical surface.
 * @param radius The radius of the sphere to integrate over.
 * @param center  The center of the sphere to integrate over.
 * @param debugOut Determines whether messages will print when a pointwise large error occurs.
 * @return The integral @f$\int_\Gamma data(x) dS_x @f$
*/
RectCoord IntegrateDiscrete(const SphereData& data,  double radius = 1, RectCoord center = RectCoord())
{

    RectCoord total = 0.0;

    for (int p = 0; p < NUMTRAPNODES; p++)
        for (int i = 0; i < NUMGLNODES; i++)
        {
            SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
            SphereCoord x(radius, s, center);
            total = total + radius * radius * sin(s.theta) * GLweights[i] * data[i][p] * (PI / NUMTRAPNODES) * PI;


        }

    return total;
}

/**
 * @brief Computes the sum of the data over an arbitrary spherical surface.
 * @param radius The radius of the sphere to integrate over.
 * @param center  The center of the sphere to integrate over.
 * @param debugOut Determines whether messages will print when a pointwise large error occurs.
 * @return The integral @f$\int_\Gamma data(x) dS_x @f$
*/
double IntegrateDiscrete(const ScalSphereData& data,  double radius = 1)
{

    double total = 0.0;

    for (int p = 0; p < NUMTRAPNODES; p++)
        for (int i = 0; i < NUMGLNODES; i++)
        {
            SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
            total = total + radius * radius * sin(s.theta) * GLweights[i] * data[i][p] * (PI / NUMTRAPNODES) * PI;
        }

    return total;
}

/**
 * @brief computes the max value in the coordinates of x.
*/
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

/**
 * @brief Computes the @f$ L^\infty(\Gamma) @f$ difference of the functions @f$ f_1 @f$ and @f$ f_2 @f$ over an arbitrary spherical surface.
 * @param radius The radius of the sphere to integrate over.
 * @param center  The center of the sphere to integrate over.
 * @param debugOut Determines whether messages will print when a pointwise large error occurs.
 * @return  @f$ \sup_{x \in \Gamma} |f_1(x) - f_2(x)| @f$
*/
double LInfDifference(SphericalVectorField* f1, SphericalVectorField* f2, double r = 1.0, RectCoord center = RectCoord() , SurfaceCoord axis = SurfaceCoord())
{
   
    double sup = 0;
    double curr = 0;

    for (int i = 0; i < NUMGLNODES; i++)
        for (double j = 0; j < NUMTRAPNODES; j++)
        {
            SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)j / (double)NUMTRAPNODES);
            SphereCoord x(s , center ,axis, r);
            RectCoord diff = (*f1)(x) - (*f2)(x);
            curr = max(diff);
            if (curr > sup)
                sup = curr;
            
        }

    return sup;
}

/**
 * @brief Functor corresponding to the @f$ e_r @f$ component of a given function.
*/
class dot_eR : public SphericalScalarFunction
{

private:
    SphericalVectorField* rho;
    RectCoord center;
public:

    /**
     * @brief Contructs the functor.
     * @param u The function to extract the component of.
    */
    dot_eR(SphericalVectorField* u , RectCoord c = RectCoord())
    {
        rho = u;
        name = u->name + " dot eR";
        center = c;
    }

    /**
     * @brief Evaluates the function.
     * @return @f$ u(x) \cdot e_r(x) @f$
    */
    double const operator()(const SphereCoord & x)
    {
        SphereCoord temp = recenter(x, center);
        return dot((*rho)(x), e_r(temp));
    }
};

/**
 * @brief Functor corresponding to the @f$ e_\theta @f$ component of a given function.
*/
class dot_eTheta : public SphericalScalarFunction
{

private:
    SphericalVectorField* rho;
    RectCoord center;

public:

    /**
     * @brief Contructs the functor.
     * @param u The function to extract the component of.
    */
    dot_eTheta(SphericalVectorField* u , RectCoord c = RectCoord())
    {
        rho = u;
        name = u->name + " dot eTheta";
        center = c;
    }

    /**
     * @brief Evaluates the function.
     * @return @f$ u(x) \cdot e_\theta(x) @f$
    */
    double const operator()(const SphereCoord &  x)
    {
        SphereCoord temp = recenter(x, center);
        return dot((*rho)(x), e_theta(temp));
    }
};

/**
 * @brief Functor corresponding to the @f$ e_\phi @f$ component of a given function.
*/
class dot_ePhi : public SphericalScalarFunction
{

private:
    SphericalVectorField* rho;
    RectCoord center;

public:

    /**
     * @brief Contructs the functor.
     * @param u The function to extract the component of.
    */
    dot_ePhi(SphericalVectorField* u , RectCoord c = RectCoord())
    {
        rho = u;
        name = u->name + " dot ePhi";
        center = c;
    }

    /**
     * @brief Evaluates the function.
     * @return @f$ u(x) \cdot e_\phi(x) @f$
    */
    double const operator()(const SphereCoord & x)
    {
        SphereCoord temp = recenter(x, center);
        return dot((*rho)(x), e_phi(temp));
    }
};

/**
 * @brief Computes the numerical partial derivative of a function with respect to @f$ r @f$.
*/
class NdR : public SphericalScalarFunction
{
private:
    SphericalScalarFunction* rho;
    RectCoord center;

public:
    /**
     * @brief Constructs the functor.
     * @param u The function to differentiate.
    */
    NdR(SphericalScalarFunction* u , RectCoord c = RectCoord())
    {
        rho = u;
        center = c;
        name = u->name + " dR";
    }

    /**
     * @brief Evaluates the derivative using finite difference.
     * @return @f$ \frac{\partial u}{\partial r}\left(x\right) @f$
    */
    double const operator()(const SphereCoord & x)
    {

        SphereCoord temp = recenter(x, center);
        double step = 1e-3;
        SphereCoord xplus = SphereCoord(temp.rho + step, temp.s , center);
        SphereCoord xminus = SphereCoord(temp.rho - step, temp.s, center);

        return ((*rho)(xplus) - (*rho)(xminus)) / (2 * step);
    }

};

/**
 * @brief Computes the numerical partial derivative of a function with respect to @f$ \theta @f$.
*/
class NdTheta : public SphericalScalarFunction
{
private:
    SphericalScalarFunction* rho;
    RectCoord center;


public:
    /**
    * @brief Constructs the functor.
    * @param u The function to differentiate.
   */
    NdTheta(SphericalScalarFunction* u, RectCoord c = RectCoord())
    {
        rho = u;
        center = c;
        name = u->name + " dTheta";
    }

    /**
     * @brief Evaluates the derivative using finite difference.
     * @return @f$ \frac{\partial u}{\partial \theta}\left(x\right) @f$
    */
    double  operator()(const SphereCoord & x) const
    {
        SphereCoord temp = recenter(x, center);
        double step = 1e-3;
        SphereCoord xplus = SphereCoord(temp.rho, temp.s + SurfaceCoord(step , 0) , center);
        SphereCoord xminus = SphereCoord(temp.rho, temp.s + SurfaceCoord(-step , 0), center);

        return ((*rho)(xplus) - (*rho)(xminus)) / (2 * step);
    }

};

/**
 * @brief Computes the numerical partial derivative of a function with respect to @f$ \phi @f$.
*/
class NdPhi : public SphericalScalarFunction
{
private:
    SphericalScalarFunction* rho;
    RectCoord center;

public:
    /**
    * @brief Constructs the functor.
    * @param u The function to differentiate.
   */
    NdPhi(SphericalScalarFunction* u, RectCoord c = RectCoord())
    {
        rho = u;
        center = c;
        name = u->name + " dPhi";
    }

    /**
     * @brief Evaluates the derivative using finite difference.
     * @return @f$ \frac{\partial u}{\partial \phi}\left(x\right) @f$
    */
    double const operator()(const SphereCoord & x)
    {
        SphereCoord temp = recenter(x, center);
        double step = 1e-3;
        SphereCoord xplus = SphereCoord(temp.rho, temp.s + SurfaceCoord(0, step), center);
        SphereCoord xminus = SphereCoord(temp.rho, temp.s + SurfaceCoord(0,-step),center);

        return ((*rho)(xplus) - (*rho)(xminus)) / (2 * step);
    }

};

/**
 * @brief Creates a SphereData from a function.
 * @param f the function to discretize.
 * @param center canter of the sphere to evaluate over.
 * @return a pointer to SphereData whos entries are the values of f at the quadrature nodes.
*/
SphereData discretize(const SphericalVectorField* f , RectCoord center = RectCoord())
{
    SphereData temp;

    temp.resize((int)NUMGLNODES);

    for (int t = 0; t < (int)NUMGLNODES; t++)
    {
        temp[t].resize((int)NUMTRAPNODES);
        for (int p = 0; p < (int)NUMTRAPNODES; p++)
        {
            SurfaceCoord s(PI / 2.0 * (GLnodes[t] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
            SphereCoord x(s, center);
            temp[t][p] = (*f)(x);
        }
    }
    return temp;
}