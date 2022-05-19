#pragma once
#include "Matrix.h"
#include "SphericalHarmonics.h"


double fV(double r, int n, double rhohatV, double rhohatW)
{

    if (n == 0)
        return 0;

    double m = double(n);

    double temp1 = m / ((2.0 * m + 1.0) * (2.0 * m + 3.0) * pow(r, n + 2));
    temp1 *= rhohatV;

    double temp2 = m / (4.0 * m + 2.0) * (1.0 / pow(r, n + 2) - 1.0 / pow(r, n));
    temp2 *= rhohatW;

    return temp1 + temp2;
}

double fW(double r, int n, double rhohatW)
{
    if (n == 0)
        return 0;

    double m = double(n);

    double temp1 = (m + 1.0) / ((2.0 * m + 1.0) * (2.0 * m - 1.0) * pow(r, n));
    temp1 *= rhohatW;





    return temp1;

}

double fX(double r, int n, double rhohatX)
{
    if (n == 0)
        return 0;

    double m = double(n);

    double temp = 1.0 / ((2.0 * m + 1.0) * pow(r, n + 1));
    temp *= rhohatX;

    return temp;
}

double g(double r, int n, double rhohatW)
{
    double m = double(n);

    double temp = m / pow(r, n + 1);
    temp *= rhohatW;

    return temp;
}

double fVprime(double r, int n, double rhohatV, double rhohatW)
{

    if (n == 0)
        return 0;

    double m = double(n);

    double temp1 = -m * (m + 2) / ((2.0 * m + 1.0) * (2.0 * m + 3.0) * pow(r, n + 3));
    temp1 *= rhohatV;

    double temp2 = m / (4.0 * m + 2.0) * (-(m + 2) / pow(r, n + 3) + m / pow(r, n + 1));
    temp2 *= rhohatW;

    return temp1 + temp2;
}

double fWprime(double r, int n, double rhohatW)
{
    if (n == 0)
        return 0;

    double m = double(n);

    double temp1 = -(m + 1.0) * m / ((2.0 * m + 1.0) * (2.0 * m - 1.0) * pow(r, n + 1));
    temp1 *= rhohatW;





    return temp1;

}

double fXprime(double r, int n, double rhohatX)
{
    if (n == 0)
        return 0;

    double m = double(n);

    double temp = -(m + 1) / ((2.0 * m + 1.0) * pow(r, n + 2));
    temp *= rhohatX;

    return temp;
}



class SingleLayerHarmonicTerm : public SphericalVectorField
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

    VReal vr;
    VImag vi;
    WReal wr;
    WImag wi;
    XReal xr;
    XImag xi;

public:

    SingleLayerHarmonicTerm() {}
    SingleLayerHarmonicTerm(VSHSeries& series, int m, int n)
    {

        m_ = m;
        n_ = n;
        vr.reset(m, n);
        vi.reset(m, n);
        wr.reset(m, n);
        wi.reset(m, n);
        xr.reset(m, n);
        xi.reset(m, n);

        rhohatVr = series.terms[n][m].rhohatVr;
        rhohatVi = series.terms[n][m].rhohatVi;
        rhohatWr = series.terms[n][m].rhohatWr;
        rhohatWi = series.terms[n][m].rhohatWi;
        rhohatXr = series.terms[n][m].rhohatXr;
        rhohatXi = series.terms[n][m].rhohatXi;

    }

    SphereCoord operator()(SphereCoord s)
    {
        SphereCoord temp = fV(s.rho, n_, rhohatVr, rhohatWr) * vr(s) + fV(s.rho, n_, rhohatVi, rhohatWi) * vi(s);
        temp = temp + fW(s.rho, n_, rhohatWr) * wr(s) + fW(s.rho, n_, rhohatWi) * wi(s);
        temp = temp + fX(s.rho, n_, rhohatXr) * xr(s) + fX(s.rho, n_, rhohatXi) * xi(s);

        return temp;
    }

    double eR(SphereCoord s)
    {
        double temp = fV(s.rho, n_, rhohatVr, rhohatWr) * vr.eR(s) + fV(s.rho, n_, rhohatVi, rhohatWi) * vi.eR(s);
        temp = temp + fW(s.rho, n_, rhohatWr) * wr.eR(s) + fW(s.rho, n_, rhohatWi) * wi.eR(s);
        temp = temp + fX(s.rho, n_, rhohatXr) * xr.eR(s) + fX(s.rho, n_, rhohatXi) * xi.eR(s);

        return temp;
    }

    double eTheta(SphereCoord s)
    {
        double temp = fV(s.rho, n_, rhohatVr, rhohatWr) * vr.eTheta(s) + fV(s.rho, n_, rhohatVi, rhohatWi) * vi.eTheta(s);
        temp = temp + fW(s.rho, n_, rhohatWr) * wr.eTheta(s) + fW(s.rho, n_, rhohatWi) * wi.eTheta(s);
        temp = temp + fX(s.rho, n_, rhohatXr) * xr.eTheta(s) + fX(s.rho, n_, rhohatXi) * xi.eTheta(s);

        return temp;
    }

    double ePhi(SphereCoord s)
    {
        double temp = fV(s.rho, n_, rhohatVr, rhohatWr) * vr.ePhi(s) + fV(s.rho, n_, rhohatVi, rhohatWi) * vi.ePhi(s);
        temp = temp + fW(s.rho, n_, rhohatWr) * wr.ePhi(s) + fW(s.rho, n_, rhohatWi) * wi.ePhi(s);
        temp = temp + fX(s.rho, n_, rhohatXr) * xr.ePhi(s) + fX(s.rho, n_, rhohatXi) * xi.ePhi(s);

        return temp;
    }

    double eR_dR(SphereCoord s)
    {
        double temp = fVprime(s.rho, n_, rhohatVr, rhohatWr) * vr.eR(s) + fVprime(s.rho, n_, rhohatVi, rhohatWi) * vi.eR(s);
        temp = temp + fWprime(s.rho, n_, rhohatWr) * wr.eR(s) + fWprime(s.rho, n_, rhohatWi) * wi.eR(s);
        temp = temp + fXprime(s.rho, n_, rhohatXr) * xr.eR(s) + fXprime(s.rho, n_, rhohatXi) * xi.eR(s);

        return temp;
    }

    double eTheta_dR(SphereCoord s)
    {
        double temp = fVprime(s.rho, n_, rhohatVr, rhohatWr) * vr.eTheta(s) + fVprime(s.rho, n_, rhohatVi, rhohatWi) * vi.eTheta(s);
        temp = temp + fWprime(s.rho, n_, rhohatWr) * wr.eTheta(s) + fWprime(s.rho, n_, rhohatWi) * wi.eTheta(s);
        temp = temp + fXprime(s.rho, n_, rhohatXr) * xr.eTheta(s) + fXprime(s.rho, n_, rhohatXi) * xi.eTheta(s);

        return temp;
    }

    double ePhi_dR(SphereCoord s)
    {
        double temp = fVprime(s.rho, n_, rhohatVr, rhohatWr) * vr.ePhi(s) + fVprime(s.rho, n_, rhohatVi, rhohatWi) * vi.ePhi(s);
        temp = temp + fWprime(s.rho, n_, rhohatWr) * wr.ePhi(s) + fWprime(s.rho, n_, rhohatWi) * wi.ePhi(s);
        temp = temp + fXprime(s.rho, n_, rhohatXr) * xr.ePhi(s) + fXprime(s.rho, n_, rhohatXi) * xi.ePhi(s);

        return temp;
    }

    double eR_dTheta(SphereCoord s)
    {
        double temp = fV(s.rho, n_, rhohatVr, rhohatWr) * vr.eR_dTheta(s) + fV(s.rho, n_, rhohatVi, rhohatWi) * vi.eR_dTheta(s);
        temp = temp + fW(s.rho, n_, rhohatWr) * wr.eR_dTheta(s) + fW(s.rho, n_, rhohatWi) * wi.eR_dTheta(s);
        temp = temp + fX(s.rho, n_, rhohatXr) * xr.eR_dTheta(s) + fX(s.rho, n_, rhohatXi) * xi.eR_dTheta(s);

        return temp;
    }

    double eTheta_dTheta(SphereCoord s)
    {
        double temp = fV(s.rho, n_, rhohatVr, rhohatWr) * vr.eTheta_dTheta(s) + fV(s.rho, n_, rhohatVi, rhohatWi) * vi.eTheta_dTheta(s);
        temp = temp + fW(s.rho, n_, rhohatWr) * wr.eTheta_dTheta(s) + fW(s.rho, n_, rhohatWi) * wi.eTheta_dTheta(s);
        temp = temp + fX(s.rho, n_, rhohatXr) * xr.eTheta_dTheta(s) + fX(s.rho, n_, rhohatXi) * xi.eTheta_dTheta(s);

        return temp;
    }

    double ePhi_dTheta(SphereCoord s)
    {
        double temp = fV(s.rho, n_, rhohatVr, rhohatWr) * vr.ePhi_dTheta(s) + fV(s.rho, n_, rhohatVi, rhohatWi) * vi.ePhi_dTheta(s);
        temp = temp + fW(s.rho, n_, rhohatWr) * wr.ePhi_dTheta(s) + fW(s.rho, n_, rhohatWi) * wi.ePhi_dTheta(s);
        temp = temp + fX(s.rho, n_, rhohatXr) * xr.ePhi_dTheta(s) + fX(s.rho, n_, rhohatXi) * xi.ePhi_dTheta(s);

        return temp;
    }

    double eR_dPhi(SphereCoord s)
    {
        double temp = fV(s.rho, n_, rhohatVr, rhohatWr) * vr.eR_dPhi(s) + fV(s.rho, n_, rhohatVi, rhohatWi) * vi.eR_dPhi(s);
        temp = temp + fW(s.rho, n_, rhohatWr) * wr.eR_dPhi(s) + fW(s.rho, n_, rhohatWi) * wi.eR_dPhi(s);
        temp = temp + fX(s.rho, n_, rhohatXr) * xr.eR_dPhi(s) + fX(s.rho, n_, rhohatXi) * xi.eR_dPhi(s);

        return temp;
    }

    double eTheta_dPhi(SphereCoord s)
    {
        double temp = fV(s.rho, n_, rhohatVr, rhohatWr) * vr.eTheta_dPhi(s) + fV(s.rho, n_, rhohatVi, rhohatWi) * vi.eTheta_dPhi(s);
        temp = temp + fW(s.rho, n_, rhohatWr) * wr.eTheta_dPhi(s) + fW(s.rho, n_, rhohatWi) * wi.eTheta_dPhi(s);
        temp = temp + fX(s.rho, n_, rhohatXr) * xr.eTheta_dPhi(s) + fX(s.rho, n_, rhohatXi) * xi.eTheta_dPhi(s);

        return temp;
    }

    double ePhi_dPhi(SphereCoord s)
    {
        double temp = fV(s.rho, n_, rhohatVr, rhohatWr) * vr.ePhi_dPhi(s) + fV(s.rho, n_, rhohatVi, rhohatWi) * vi.ePhi_dPhi(s);
        temp = temp + fW(s.rho, n_, rhohatWr) * wr.ePhi_dPhi(s) + fW(s.rho, n_, rhohatWi) * wi.ePhi_dPhi(s);
        temp = temp + fX(s.rho, n_, rhohatXr) * xr.ePhi_dPhi(s) + fX(s.rho, n_, rhohatXi) * xi.ePhi_dPhi(s);

        return temp;
    }

};

class SingleLayerHarmonic : public VectorFieldSum
{
private:

    Legendre P;
    std::vector<std::vector<SingleLayerHarmonicTerm>> terms;

public:
    SingleLayerHarmonic(int n)
    {
        terms.resize(n + 1);

        for (int i = 0; i < terms.size(); i++)
            terms[i].resize(2 * i + 1);
    }

    SingleLayerHarmonic(VSHSeries& series)
    {
        terms.resize(series.N + 1);
        for (int i = 0; i <= series.N; i++)
            terms[i].resize(2 * i + 1);
        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = SingleLayerHarmonicTerm(series, m, n);
                append(terms[n][m]);
            }


    }

    void solve(VSHSeries& series)
    {
        terms.resize(series.N + 1);
        for (int i = 0; i <= series.N; i++)
            terms[i].resize(2 * i + 1);
        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = SingleLayerHarmonicTerm(series, m, n);
                append(terms[m][n]);
            }
    }

    void solve(SphericalVectorField rho, int N, int numgrids, double gradstep)
    {
        VSHSeries series(N, numgrids);
        series.approximate(&rho);

        solve(series);
    }

    SphereCoord operator()(SphereCoord s)
    {
        int N = terms[0].size();
        P.populate(cos(s.s.theta), N, N);
        return VectorFieldSum::operator()(s);

    }



    double eR(SphereCoord s)
    {
        double temp = 0.0;

        for(auto i : terms)
            for (auto j : i)
            {
                temp += j.eR(s);
            }

        return temp;
    }

    double eTheta(SphereCoord s)
    {
        double temp = 0.0;

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.eTheta(s);
            }

        return temp;
    }

    double ePhi(SphereCoord s)
    {
        double temp = 0.0;

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.ePhi(s);
            }

        return temp;
    }

    double eR_dR(SphereCoord s)
    {
        double temp = 0.0;

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.eR_dR(s);
            }

        return temp;
    }

    double eTheta_dR(SphereCoord s)
    {
        double temp = 0.0;

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.eTheta_dR(s);
            }

        return temp;
    }

    double ePhi_dR(SphereCoord s)
    {
        double temp = 0.0;

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.ePhi_dR(s);
            }

        return temp;
    }

    double eR_dTheta(SphereCoord s)
    {
        double temp = 0.0;

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.eR_dTheta(s);
            }

        return temp;
    }

    double eTheta_dTheta(SphereCoord s)
    {
        double temp = 0.0;

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.eTheta_dTheta(s);
            }

        return temp;
    }

    double ePhi_dTheta(SphereCoord s)
    {
        double temp = 0.0;

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.ePhi_dTheta(s);
            }

        return temp;
    }

    double eR_dPhi(SphereCoord s)
    {
        double temp = 0.0;

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.eR_dPhi(s);
            }

        return temp;
    }

    double eTheta_dPhi(SphereCoord s)
    {
        double temp = 0.0;

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.eTheta_dPhi(s);
            }

        return temp;
    }

    double ePhi_dPhi(SphereCoord s)
    {
        double temp = 0.0;

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.ePhi_dPhi(s);
            }

        return temp;
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

        rhohatWr = series.terms[n][m].rhohatWr;
        rhohatWi = series.terms[n][m].rhohatWi;

    }

    double operator()(SphereCoord s)
    {
        double temp = g(s.rho, n_, rhohatWr) * Yr(s.s) + g(s.rho, n_, rhohatWi) * Yi(s.s);

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
        terms.resize(n);

        for (int i = 0; i < terms.size(); i++)
            terms[i].resize(2 * i + 1);
    }

    StokesPressure(VSHSeries& series)
    {
        terms.resize(series.N + 1);
        for (int i = 0; i <= series.N; i++)
            terms[i].resize(2 * i + 1);

        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = StokesPressureTerm(series, m, n);
                append(terms[n][m]);
            }

    }

    void solve(VSHSeries& series)
    {
        terms.resize(series.N + 1);
        for (int i = 0; i <= series.N; i++)
            terms[i].resize(2 * i + 1);

        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = StokesPressureTerm(series, m, n);
                append(terms[n][m]);
            }
    }

    void solve(SphericalVectorField rho, int N, int numgrids)
    {
        VSHSeries series(N, numgrids);
        series.approximate(&rho);

        solve(series);
    }

    double operator()(SphereCoord s)
    {
        int N = terms[0].size();
        P.populate(cos(s.s.theta), N, N);
        return SphereScalFunctionSum::operator()(s);

    }

};

vec<3> RecttoVec(RectCoord x)
{
    vec<3> temp;
    temp[0] = x.x;
    temp[1] = x.y;
    temp[2] = x.z;
    return temp;
}
vec<3> SpheretoVec(SphereCoord s)
{
    RectCoord x(s);
    return RecttoVec(x);

}

SphereCoord VectoSphere(vec<3> v)
{
    RectCoord temp;

    temp.x = v[0];
    temp.y = v[1];
    temp.z = v[2];

    return RectToSphere(temp);
}
class SingleLayerDirect : public SphericalVectorField
{
private:

    int n;
    SphericalVectorField* u;

   Matrix<3, 3> Kernel(SphereCoord x, SphereCoord y)
   {
       Matrix<3, 3> temp;

       //convert to the argument of the stokeslet.
       vec<3> r = SpheretoVec(x - y);
       double R = sqrt(dot(r,r));

       //iterate over entries over the tensor.
       for(int i = 0; i < 3; i++)
           for (int j = 0; j < 3; j++)
           {
               //formula for the i,jth entry of the stokeslet
               temp[i][j] = (double)(i ==j)/ R + r[i] * r[j] / (R * R * R);
               temp[i][j] /= 8.0 * PI;
           }

       return temp;
   } 
public:

    SingleLayerDirect(SphericalVectorField *rho, int N) 
    {
        n = N;
        u = rho;
    }


    SphereCoord operator()(SphereCoord x)
    {
        /*const double GLnodes[16] = {-0.0950125098376374  , 0.0950125098376374  , -0.2816035507792589 , 0.2816035507792589 ,  -0.4580167776572274  , 0.4580167776572274,
                         -0.6178762444026438  , 0.6178762444026438 , -0.7554044083550030 , 0.7554044083550030 , -0.8656312023878318 , 0.8656312023878318  ,
                         -0.9445750230732326 , 0.9445750230732326  , -0.9894009349916499  , 0.9894009349916499 };

   */
        //Gauss Legendre nodes
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

        //Gauss Legendre weights
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

        vec<3> total;
        //loop over trapezoidal(outside) and legendre(inside) nodes.
        for (int p = 0; p < n; p++)
            for (int i = 0; i < 32; i++)
            {
                //generate point on unit sphere, note theta is the first argument, phi the second.
                SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1) , 2.0 * PI * (double)p / (double)n);
                //convert to spherical coordinate.
                SphereCoord y(s);
                //add contribution at y to integral.
                total = total +  dot(Kernel(x , y), SpheretoVec((*u)(y))) *sin(s.theta) * GLweights[i] * (PI / n) * PI;
            }

        return VectoSphere(total);
    }
};

class StokesTractionDirect : public SphericalVectorField
{
private:

    int N;
    SphericalVectorField* u;

    Matrix<3, 3> Kernel(SphereCoord x, SphereCoord y)
    {
        Matrix<3, 3> temp;

        //convert to the argument of the stokeslet.
        vec<3> r = SpheretoVec(x - y);
        double R = sqrt(dot(r, r));
        double rn = dot(r,SpheretoVec(e_r(x)));
        //iterate over entries over the tensor.
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                //formula for the i,jth entry of the stokeslet
                temp[i][j] = r[i] * r[j] * rn / pow(R , 5);
                temp[i][j] *= -0.75 / PI;
            }

        return temp;
    }
public:

    StokesTractionDirect(SphericalVectorField* rho, int Nn)
    {
        N = Nn;
        u = rho;
    }


    SphereCoord operator()(SphereCoord x)
    {
        /*const double GLnodes[16] = {-0.0950125098376374  , 0.0950125098376374  , -0.2816035507792589 , 0.2816035507792589 ,  -0.4580167776572274  , 0.4580167776572274,
                         -0.6178762444026438  , 0.6178762444026438 , -0.7554044083550030 , 0.7554044083550030 , -0.8656312023878318 , 0.8656312023878318  ,
                         -0.9445750230732326 , 0.9445750230732326  , -0.9894009349916499  , 0.9894009349916499 };

   */
   //Gauss Legendre nodes
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

        //Gauss Legendre weights
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

        vec<3> total;
        //loop over trapezoidal(outside) and legendre(inside) nodes.
        for (int p = 0; p < N; p++)
            for (int i = 0; i < 32; i++)
            {
                //generate point on unit sphere, note theta is the first argument, phi the second.
                SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)N);
                //convert to spherical coordinate.
                SphereCoord y(s);
                //add contribution at y to integral.
                total = total + dot(Kernel(x, y), SpheretoVec((*u)(y))) * sin(s.theta) * GLweights[i] * (PI / N) * PI;
            }

        return VectoSphere(total);
    }
};


class StokesTractionHarmonic : public SphericalVectorField
{
private:
    SingleLayerHarmonic* slh;
    StokesPressure* sp;
    Matrix<3, 3> kernel;


public:

    StokesTractionHarmonic(SingleLayerHarmonic* u, StokesPressure* p)
    {
        slh = u;
        sp = p;
    }

    SphereCoord operator()(SphereCoord s)
    {
        kernel[0][0] = slh->eR_dR(s) - sp->operator()(s);
        kernel[1][0] = kernel[0][1] = (s.rho * slh->eTheta_dR(s) - slh->eTheta(s)) / (2.0 * s.rho) + slh->eR_dTheta(s) / (2.0 * s.rho);
        kernel[2][0] = kernel[0][2] = slh->eR_dPhi(s) / (2.0 * s.rho * sin(s.s.theta)) + (s.rho * slh->ePhi_dR(s) - slh->ePhi(s)) / (2.0 * s.rho);
        kernel[1][1] = slh->eTheta_dTheta(s) / s.rho + slh->eR(s) / s.rho - sp->operator()(s);
        kernel[2][1] = kernel[1][2] = (sin(s.s.theta) * slh->ePhi_dTheta(s) - cos(s.s.theta) * slh->ePhi(s)) / (2.0 * s.rho * sin(s.s.theta)) + slh->eTheta_dPhi(s) / (2.0 * s.rho * sin(s.s.theta));
        kernel[2][2] = slh->ePhi_dPhi(s) / (s.rho * sin(s.s.theta)) + slh->eR(s) / s.rho + slh->eTheta(s) * cos(s.s.theta) / (s.rho * sin(s.s.theta)) - sp->operator()(s);

        vec<3> normal = SpheretoVec(e_r(s));

        return VectoSphere(dot(kernel, normal));

    }

};
