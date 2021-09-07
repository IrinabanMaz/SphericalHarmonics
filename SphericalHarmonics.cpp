// SphericalHarmonics.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>
#include <vector>

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
    YReal(int a, int b )
    {
        m = a;
        n = b;
       
    }


    //Implements the function Y 
    double operator()(SurfaceCoord s)
    {
    
            
       
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;
        
        coef *= poly(cos(s.theta) , m, n);

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

        coef *= poly(cos(s.theta),m,n);

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
    VReal(int a, int b, double steps)
    {
        m = a;
        n = b;
        h = steps;
    }

    SphereCoord operator()(SphereCoord x)
    {
        SphereCoord temp;
        YReal y(m, n );

        temp = -((double)n+1.0) * y(x.s) * SphereCoord(1, x.s);
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
    VImag(int a, int b, double steps)
    {
        m = a;
        n = b;
        h = steps;
    }

    SphereCoord operator()(SphereCoord x)
    {
        SphereCoord temp;
        YImag y(m, n);

        temp = -((double)n + 1.0) * y(x.s) * SphereCoord(1, x.s);
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

    WReal(int a, int b, double steps)
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

    WImag(int a, int b, double steps)
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

    XReal(int a, int b, double steps)
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

    XImag(int a, int b, double steps)
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
    public:

    double coefr;
    double coefi;
    double denom;
    SphericalVectorField* fr;
    SphericalVectorField* fi;


public:
    SphericalHarmonic() { coefr = 1; fr = nullptr; coefi = 1; fi = nullptr;}

    SphericalHarmonic(double cr, double ci ,  SphericalVectorField& f0r, SphericalVectorField& f0i)
    {
        coefr = cr;
        coefi = ci;
        fr = &f0r;
        fi = &f0i;

        //denom = L2InnerProduct(fr, fr, 20) + L2InnerProduct(fi , fi , 20);
    }

    void setpointers(SphericalVectorField* r, SphericalVectorField* i)
    {
        fr = r;
        fi = i;
    }

    SphereCoord operator()(SphereCoord x)
    {
        SphereCoord frx = (*fr)(x);
        SphereCoord Frnormalized; 
        if (denom < 1e-16)
            Frnormalized = SphereCoord(0.0, SurfaceCoord(0.0, 0.0));
        else
            Frnormalized = (1.0 / denom)* frx;

        SphereCoord fix = (*fi)(x);
        SphereCoord Finormalized;
        if (denom < 1e-16)
            Finormalized = SphereCoord(0.0, SurfaceCoord(0.0, 0.0));
        else
            Finormalized = (1.0 / denom) * fix;

        return coefr * Frnormalized + coefi * Finormalized;
    }

};

class VSHSeries : public VectorFieldSum
{

private:

    

    Legendre P;

    VReal** vr;
    WReal** wr;
    XReal** xr;
    VImag** vi;
    WImag** wi;
    XImag** xi;

public:
    int N;
    int NUMGRIDS;
    double GRADSTEP;
    SphericalHarmonic** Vterm = nullptr;
    SphericalHarmonic** Wterm = nullptr;
    SphericalHarmonic** Xterm = nullptr;




    VSHSeries(int n, int numgrids, double gradstep)
    {
        N = n;
        NUMGRIDS = numgrids;
        GRADSTEP = gradstep;

        vr = new VReal * [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            vr[i] = new VReal[N + 1];

        wr = new WReal * [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            wr[i] = new WReal[N + 1];

        xr = new XReal * [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            xr[i] = new XReal[N + 1];

        vi = new VImag * [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            vi[i] = new VImag[N + 1];

        wi = new WImag * [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            wi[i] = new WImag[N + 1];

        xi = new XImag * [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            xi[i] = new XImag[N + 1];

        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                vr[m][n] = VReal(m, n, GRADSTEP);
                vr[m + N][n] = VReal(-m, n, GRADSTEP);
                wr[m][n] = WReal(m, n, GRADSTEP);
                wr[m + N][n] = WReal(-m, n, GRADSTEP);
                xr[m][n] = XReal(m, n, GRADSTEP);
                xr[m + N][n] = XReal(-m, n, GRADSTEP);

                vi[m][n] = VImag(m, n, GRADSTEP);
                vi[m + N][n] = VImag(-m, n, GRADSTEP);
                wi[m][n] = WImag(m, n, GRADSTEP);
                wi[m + N][n] = WImag(-m, n, GRADSTEP);
                xi[m][n] = XImag(m, n, GRADSTEP);
                xi[m + N][n] = XImag(-m, n, GRADSTEP);
            }

        Vterm = new SphericalHarmonic * [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            Vterm[i] = new SphericalHarmonic[N + 1];

        
      
        Wterm = new SphericalHarmonic * [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            Wterm[i] = new SphericalHarmonic[N + 1];
        
        
        
        Xterm = new SphericalHarmonic * [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            Xterm[i] = new SphericalHarmonic[N + 1];
        


    }

    void approximate(SphericalVectorField r)
    {

        //calculate the inner products;

    //for V
        double** rhohatVr = new double* [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            rhohatVr[i] = new double[N + 1];

        double** rhohatVi = new double* [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            rhohatVi[i] = new double[N + 1];

        std::cout << "Computing inner product coefficiencts for V... " << std::endl;
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

            }

        std::cout << "Done!" << std::endl << std::endl;
        static auto offset = SphericalVectorField([](SphereCoord s) { return -1.0 * VectorFieldSum()(s); });
        append(offset);

        for (int i = 0; i < 2 * N + 1; i++)
        {
            delete[] rhohatVi[i];
            delete[] rhohatVr[i];
            delete[] rhohatWr[i];
            delete[] rhohatWi[i];
            delete[] rhohatXr[i];
            delete[] rhohatXi[i];
        }

        delete[] rhohatVi;
        delete[] rhohatVr;
        delete[] rhohatWr;
        delete[] rhohatWi;
        delete[] rhohatXr;
        delete[] rhohatXi;

    }

    void operator ~()
    {
        
        for (int i = 0; i < 2 * N + 1; i++)
            delete[] Vterm[i];

        delete[] Vterm;

        
        for (int i = 0; i < 2 * N + 1; i++)
            delete[] Wterm[i];

        delete[] Wterm;

        
        for (int i = 0; i < 2 * N + 1; i++)
            delete[] Xterm[i];

        delete[] Xterm;

        vr = new VReal * [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            delete[] vr[i];

        wr = new WReal * [2 * N + 1];
        for (int i = 0; i < 2 * N + 1; i++)
            delete[] wr[i];

        
        for (int i = 0; i < 2 * N + 1; i++)
            delete[] xr[i];
        
        delete[] xr;


        
        for (int i = 0; i < 2 * N + 1; i++)
            delete[] vi[i];
        
        delete[] vi;

        
        for (int i = 0; i < 2 * N + 1; i++)
            delete[] wi[i];
        
        delete[] wi;
        
        for (int i = 0; i < 2 * N + 1; i++)
            delete[] xi[i];

        delete[] xi;
    }

    VSHSeries operator =(VSHSeries vshs)
    {
        *this = VSHSeries(vshs.N, vshs.NUMGRIDS, vshs.GRADSTEP);

        if (Vterm != nullptr)
        {
            for (int n = 0; n <= N; n++)
                for (int m = 0; m <= n; m++)
                {

                    Vterm[m][n] = vshs.Vterm[m][n];
                    Vterm[m][n].setpointers(&vr[m][n], &vi[m][n]);
                    if (m > 0)
                    {
                        Vterm[m + N][n] = vshs.Vterm[m + N][n];
                        Vterm[m + N][n].setpointers(&vr[m+N][n], &vi[m+N][n]);
                    }

                    Wterm[m][n] = vshs.Wterm[m][n];
                    Wterm[m][n].setpointers(&wr[m][n], &wi[m][n]);
                    if (m > 0)
                    {
                        Wterm[m + N][n] = vshs.Wterm[m + N][n];
                        Wterm[m + N][n].setpointers(&wr[m + N][n], &wi[m + N][n]);
                    }

                    Xterm[m][n] = vshs.Xterm[m][n];
                    Xterm[m][n].setpointers(&xr[m][n], &xi[m][n]);
                    if (m > 0)
                    {
                        Xterm[m + N][n] = vshs.Xterm[m + N][n];
                        Xterm[m + N][n].setpointers(&xr[m + N][n], &xi[m + N][n]);
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


                   

                }

            static auto offset = SphericalVectorField([](SphereCoord s) {return -1.0 * VectorFieldSum()(s); });
            append(offset);

        }

        return *this;
    }

    VSHSeries(const VSHSeries & vshs)
    {
        *this = vshs;
    }

    SphereCoord operator()(SphereCoord x)
    {
        P.populate(cos(x.s.theta), N, N);

        return VectorFieldSum::operator()(x);
    }
};

double fV(double r, int n, double rhohatV , double rhohatW)
{

    if (n == 0)
        return 0;

    double m = double(n);

    double temp1 = m / ((2.0 * m + 1.0) * (2.0 * m + 3.0) * pow(r, n + 2));
    temp1 *= rhohatV / ((2.0 * m + 1.0) * (m + 1.0));
    
    double temp2 = m / (4.0 * m + 2.0) * (1.0 / pow(r, n + 2) - 1.0 / pow(r, n));
    temp2 *= rhohatW / ((2.0 * m + 1.0) * m);

    return temp1 + temp2;
}

double fW(double r, int n, double rhohatW)
{
    if (n == 0)
        return 0;

    double m = double(n);

    double temp1 = (m + 1.0) / ((2.0 * m + 1.0) * (2.0 * m - 1.0) * pow(r, n));
    temp1 *= rhohatW / ((2.0 * m + 1.0) * m);

    

   

    return temp1;

}

double fX(double r, int n, double rhohatX)
{
    if (n == 0)
        return 0;

    double m = double(n);

    double temp = 1.0 / ((2.0 * m + 1.0) * pow(r, n + 1));
    temp *= rhohatX / (m * (m + 1.0));

    return temp;
}

double g(double r, int n, double rhohatW)
{
    double m = double(n);

    double temp = 1.0 / pow(r, n + 1);
    temp *= rhohatW / (2.0 * m + 1.0);

    return temp;
}

class StokesFlowTerm : public SphericalVectorField
{
private:

    int m_;
    int n_;

    double rhohatVr;
    double rhohatVi;
    double rhohatWr;
    double rhohatWi;
    double rhohatXr;
    double rhohatXi;

    SphericalVectorField* Vr;
    SphericalVectorField* Vi;
    SphericalVectorField* Wr;
    SphericalVectorField* Wi;
    SphericalVectorField* Xr;
    SphericalVectorField* Xi;

public:

    StokesFlowTerm(){}
    StokesFlowTerm(VSHSeries& series, int m, int n)
    {
    
        m_ = m;
        n_ = n;
        Vr = series.Vterm[m][n].fr;
        Vi = series.Vterm[m][n].fi;
        Wr = series.Wterm[m][n].fr;
        Wi = series.Wterm[m][n].fi;
        Xr = series.Xterm[m][n].fr;
        Xi = series.Xterm[m][n].fi;

        rhohatVr = series.Vterm[m][n].coefr;
        rhohatVi = series.Vterm[m][n].coefi;
        rhohatWr = series.Wterm[m][n].coefr;
        rhohatWi = series.Wterm[m][n].coefi;
        rhohatXr = series.Xterm[m][n].coefr;
        rhohatXi = series.Xterm[m][n].coefi;

    }

    SphereCoord operator()(SphereCoord s)
    {
        SphereCoord temp = fV(s.rho, n_, rhohatVr, rhohatWr) * (*Vr)(s) + fV(s.rho, n_, rhohatVi , rhohatWi) * (*Vi)(s);
        temp = temp + fW(s.rho, n_, rhohatWr) * (*Wr)(s) + fW(s.rho, n_, rhohatWi) * (*Wi)(s);
        temp = temp + fX(s.rho, n_, rhohatXr) * (*Xr)(s) + fX(s.rho, n_, rhohatXi) * (*Xi)(s);

        return temp;
    }

};

class StokesFlow : public VectorFieldSum
{
private:

    Legendre P;
    std::vector<std::vector<StokesFlowTerm>> terms;

public:
    StokesFlow(int n)
    {
        terms.resize(2 * n + 1);

        for (int i = 0; i < terms.size(); i++)
            terms[i].resize(n + 1);
    }

    StokesFlow(VSHSeries& series)
    {
        terms.resize(2 * series.N + 1);
        for (int i = 0; i < 2 *series.N + 1; i++)
            terms[i].resize(series.N + 1);
        for (int n = 0; n <= series.N; n++)
           for (int m = 0; m <= n; m++)
            {
                terms[m][n] = StokesFlowTerm(series, m, n);
                terms[m + series.N][n] = StokesFlowTerm(series, m + series.N, n);
                append(terms[m][n]);
                if (m > 0)
                    append(terms[m + series.N][n]);
            }


    }

    void solve(VSHSeries& series)
    {
        terms.resize(2 * series.N + 1);
        for (int i = 0; i < series.N; i++)
            terms[i].resize(series.N);

        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[m][n] = StokesFlowTerm(series, m, n);
                terms[m + series.N][n] = StokesFlowTerm(series, m + series.N, n);
                append(terms[m][n]);
                if (m > 0)
                    append(terms[m + series.N][n]);
            }
    }

    void solve(SphericalVectorField rho , int N, int numgrids  , double gradstep )
    {
        VSHSeries series(N, numgrids, gradstep);
        series.approximate(rho);
        
        solve(series);
    }

    SphereCoord operator()(SphereCoord s)
    {
        int N = terms[0].size();
        P.populate(cos(s.s.theta), N, N);
        return VectorFieldSum::operator()(s);

    }

};

class StokesPressureTerm : public SphericalScalarFunction
{
private:

    int m_;
    int n_;

    double rhohatWr;
    double rhohatWi;

    YReal Yr;
    YImag Yi;

public:

    StokesPressureTerm() {}
    StokesPressureTerm(VSHSeries& series, int m, int n)
    {

        m_ = m;
        n_ = n;
        Yr = YReal(m, n);
        Yi = YImag(m, n);

        rhohatWr = series.Wterm[m][n].coefr;
        rhohatWi = series.Wterm[m][n].coefi;

    }

    double operator()(SphereCoord s)
    {
        double temp =  g(s.rho, n_, rhohatWr) * Yr(s.s) + g(s.rho, n_, rhohatWi) * Yi(s.s);

        return temp;
    }

};

class StokesPressure : public SphereScalFunctionSum
{
private:

    Legendre P;
    std::vector<std::vector<StokesPressureTerm>> terms;

public:
    StokesPressure(int n)
    {
        terms.resize(2 * n + 1);

        for (int i = 0; i < terms.size(); i++)
            terms[i].resize(n);
    }

    StokesPressure(VSHSeries& series)
    {
        terms.resize(2 * series.N + 1);
        for (int i = 0; i < series.N; i++)
            terms[i].resize(series.N);

        for (int m = 0; m < 2 * series.N + 1; m++)
            for (int n = 0; n < series.N; n++)
            {
                terms[m][n] = StokesPressureTerm(series, m, n);
                append(terms[m][n]);
            }

    }

    void solve(VSHSeries& series)
    {
        terms.resize(2 * series.N + 1);
        for (int i = 0; i < series.N; i++)
            terms[i].resize(series.N);

        for (int m = 0; m < 2 * series.N + 1; m++)
            for (int n = 0; n < series.N; n++)
            {
                terms[m][n] = StokesPressureTerm(series, m, n);
                append(terms[m][n]);
            }
    }

    void solve(SphericalVectorField rho, int N, int numgrids, double gradstep)
    {
        VSHSeries series(N, numgrids, gradstep);
        series.approximate(rho);

        solve(series);
    }

    double operator()(SphereCoord s)
    {
        int N = terms[0].size();
        P.populate(cos(s.s.theta), N, N);
        SphereScalFunctionSum* ptr = this;
        ptr->operator()(s);

    }

};

int main()
{   
    //numerical parameters.
    
    //series truncation
    const int N = 5;
    
    //step size for gradient.
    const double GRADSTEP = 1e-5;

    //number of panel points per dimension in integrals.
    const int NUMGRIDS = 22;

    SphericalVectorField boundaryData([](SphereCoord s)
        {return cos(s.s.theta) * e_r(s) - sin(s.s.theta) * e_theta(s); });

    SphericalVectorField TrueSolution([](SphereCoord s)
        {
            SphereCoord temp;
            temp = RectToSphere(e_r(s)* cos(s.s.theta) * (3.0 / (2.0 * s.rho) - 1.0 / (2.0 * s.rho * s.rho * s.rho)));
            temp = temp - RectToSphere(sin(s.s.theta) * (3.0 / (4.0 * s.rho) + 1.0 / (4.0 * s.rho * s.rho * s.rho)) * e_theta(s));
            return temp;

        });

    VSHSeries bDApprox(N, NUMGRIDS, GRADSTEP);

    bDApprox.approximate(boundaryData);

    StokesFlow approxSolution(bDApprox);


    std::cout << L2Difference(&approxSolution, &TrueSolution, NUMGRIDS, 1.0);

   
    return 0;

    
}

