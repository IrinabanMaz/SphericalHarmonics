#pragma once
#include "Matrix.h"
#include "SphericalHarmonics.h"
#include <cmath>


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

    if (isnan(rhohatW))
        return 0.0;

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

    RectCoord operator()(SphereCoord s)
    {
        RectCoord temp = fV(s.rho, n_, rhohatVr, rhohatWr) * vr(s) + fV(s.rho, n_, rhohatVi, rhohatWi) * vi(s);
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

        NdR epdr(&dep);
        NdR etdr(&det);
        NdR erdr(&der);



        double ephierr = 0.0;
        double ethetaerr = 0.0;
        double ererr = 0.0;
        double ephidphierr = 0.0;
        double ephidthetaerr = 0.0;
        double ephidrerr = 0.0;
        double erdphierr = 0.0;
        double ethetadthetaerr = 0.0;
        double ethetadphierr = 0.0;
        double ethetadrerr = 0.0;
        double erdthetaerr = 0.0;
        double erdrerr = 0.0;


        for (int p = 0; p < NUMTRAPNODES; p++)
            for (int i = 0; i < NUMGLNODES; i++)
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

                diff = ePhi_dR(x) - epdr(x);
                ephidrerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;


                diff = eTheta_dTheta(x) - etdt(x);
                ethetadthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dPhi(x) - etdp(x.s);
                ethetadphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eTheta_dR(x) - etdr(x);
                ethetadrerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;


                diff = eR_dTheta(x) - erdt(x);
                erdthetaerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dPhi(x) - erdp(x);
                erdphierr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;

                diff = eR_dR(x) - erdr(x);
                erdrerr += sin(s.theta) * GLweights[i] * diff * diff * (PI / NUMTRAPNODES) * PI;
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
        std::cout << "Error in computing ephidr: " << ephidrerr << "\n";

        std::cout << "Error in computing ethetadtheta: " << ethetadthetaerr << "\n";
        std::cout << "Error in computing ethetadphi: " << ethetadphierr << "\n";
        std::cout << "Error in computing ethetadr: " << ethetadrerr << "\n";

        std::cout << "Error in computing erdtheta: " << erdthetaerr << "\n";
        std::cout << "Error in computing erdphi: " << erdphierr << "\n";
        std::cout << "Error in computing erdr: " << erdrerr << "\n";

    }


};

class SingleLayerHarmonic : public VectorFieldSum
{
private:

    Legendre P;
    
    std::vector<std::vector<SingleLayerHarmonicTerm>> terms;

public:
    
    RectCoord center;
    SingleLayerHarmonic(int n , RectCoord c = RectCoord())
    {

        center = c;
        terms.resize(n + 1);

        for (int i = 0; i < terms.size(); i++)
            terms[i].resize(i + 1);
    }

    SingleLayerHarmonic(VSHSeries& series)
    {
        terms.resize(series.N + 1);
        for (int i = 0; i <= series.N; i++)
            terms[i].resize( i + 1);
        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = SingleLayerHarmonicTerm(series, m, n);
                append(&terms[n][m]);
            }

        center = series.center;

    }

    void solve(VSHSeries& series)
    {
        terms.clear();
        terms.resize(series.N + 1);
        for (int i = 0; i <= series.N; i++)
            terms[i].resize(i + 1);
        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = SingleLayerHarmonicTerm(series, m, n);
                append(&terms[n][m]);
            }

        center = series.center;
    }

    void solve(SphericalVectorField* rho, int N, int numgrids, double gradstep , RectCoord c = RectCoord())
    {
        VSHSeries series(N, numgrids , c);
        series.approximate(rho);

        solve(series);
    }

    RectCoord operator()(SphereCoord s)
    {

        SphereCoord temp = recenter(s, center);
        int N = terms.size();
        P.populate(cos(temp.s.theta), N, N);
        return VectorFieldSum::operator()(temp);

    }



    double eR(SphereCoord s)
    {
        double temp = 0.0;

        SphereCoord Temp = recenter(s, center);

        for(auto i : terms)
            for (auto j : i)
            {
                temp += j.eR(Temp);
            }

        return temp;
    }

    double eTheta(SphereCoord s)
    {
        double temp = 0.0;

        SphereCoord Temp = recenter(s, center);

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.eTheta(Temp);
            }

        return temp;
    }

    double ePhi(SphereCoord s)
    {
        double temp = 0.0;

        SphereCoord Temp = recenter(s, center);

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.ePhi(Temp);
            }

        return temp;
    }

    double eR_dR(SphereCoord s)
    {
        double temp = 0.0;

        SphereCoord Temp = recenter(s, center);

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.eR_dR(Temp);
            }

        return temp;
    }

    double eTheta_dR(SphereCoord s)
    {
        double temp = 0.0;

        SphereCoord Temp = recenter(s, center);

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.eTheta_dR(Temp);
            }

        return temp;
    }

    double ePhi_dR(SphereCoord s)
    {
        double temp = 0.0;

        SphereCoord Temp = recenter(s, center);

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.ePhi_dR(Temp);
            }

        return temp;
    }

    double eR_dTheta(SphereCoord s)
    {
        double temp = 0.0;

        SphereCoord Temp = recenter(s, center);

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.eR_dTheta(Temp);
            }

        return temp;
    }

    double eTheta_dTheta(SphereCoord s)
    {
        double temp = 0.0;

        SphereCoord Temp = recenter(s, center);

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.eTheta_dTheta(Temp);
            }

        return temp;
    }

    double ePhi_dTheta(SphereCoord s)
    {
        double temp = 0.0;

        SphereCoord Temp = recenter(s, center);

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.ePhi_dTheta(Temp);
            }

        return temp;
    }

    double eR_dPhi(SphereCoord s)
    {
        double temp = 0.0;

        SphereCoord Temp = recenter(s, center);

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.eR_dPhi(Temp);
            }

        return temp;
    }

    double eTheta_dPhi(SphereCoord s)
    {
        double temp = 0.0;

        SphereCoord Temp = recenter(s, center);

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.eTheta_dPhi(Temp);
            }

        return temp;
    }

    double ePhi_dPhi(SphereCoord s)
    {
        double temp = 0.0;

        SphereCoord Temp = recenter(s, center);

        for (auto i : terms)
            for (auto j : i)
            {
                temp += j.ePhi_dPhi(Temp);
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

        //std::cout << temp << "\n";

        return temp;
    }

    

};

class StokesPressure : public SphereScalFunctionSum
{
private:

    Legendre P;
    std::vector<std::vector<StokesPressureTerm>> terms;

public:

    RectCoord center;

    StokesPressure(int n , RectCoord c = RectCoord())
    {

        center = c;
        terms.resize(n+1);

        for (int i = 0; i < terms.size(); i++)
            terms[i].resize(i + 1);
    }

    StokesPressure(VSHSeries& series , std::string n = "")
    {
        name = n;
        center = series.center;
        terms.resize(series.N + 1);
        for (int i = 0; i <= series.N; i++)
            terms[i].resize(i + 1);

        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = StokesPressureTerm(series, m, n);
                append(&terms[n][m]);
            }

    }

    void solve(VSHSeries& series)
    {
        terms.clear();
        terms.resize(series.N + 1);
        for (int i = 0; i <= series.N; i++)
            terms[i].resize(i + 1);

        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = StokesPressureTerm(series, m, n);
                append(&terms[n][m]);
            }

        center = series.center;
    }

    void solve(SphericalVectorField* rho, int N, int numgrids , RectCoord center = RectCoord())
    {
        VSHSeries series(N, numgrids , center);
        series.approximate(rho);

        solve(series);
    }

    double operator()(SphereCoord s)
    {

        SphereCoord temp = recenter(s, center);
        int N = terms.size();
        P.populate(cos(temp.s.theta), N, N);
        return SphereScalFunctionSum::operator()(temp);

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

RectCoord VectoRect(vec<3> temp)
{
    RectCoord x;
    x.x = temp[0];
    x.y = temp[1];
    x.z = temp[2];
    return x;
}
class SingleLayerDirect : public SphericalVectorField
{
private:

    int n;
    RectCoord center;
    SphericalVectorField* u;

   Matrix<3, 3> Kernel(SphereCoord x, SphereCoord y)
   {
       Matrix<3, 3> temp;

       //convert to the argument of the stokeslet.
       vec<3> r = RecttoVec(x - y);
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

    SingleLayerDirect(SphericalVectorField *rho, int N, RectCoord c = RectCoord(), std::string na = "Direct Single Layer")
    {
        n = N;
        u = rho;
        name = na;
        center = c;
    }


    RectCoord operator()(SphereCoord x)
    {

        

        vec<3> total = 0.0;
        //loop over trapezoidal(outside) and legendre(inside) nodes.
        for (int p = 0; p < n; p++)
            for (int i = 0; i < NUMGLNODES; i++)
            {
                //generate point on unit sphere, note theta is the first argument, phi the second.
                SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1) , 2.0 * PI * (double)p / (double)n);
                //convert to spherical coordinate.
                SphereCoord y(s , center);
                //add contribution at y to integral.
                total = total +  dot(Kernel(x , y), RecttoVec((*u)(y))) *sin(s.theta) * GLweights[i] * (PI / n) * PI;
            }

        return VectoRect(total);
    }
};



class StokesTractionDirect : public SphericalVectorField
{
private:

    int N;
    SphericalVectorField* u;
    RectCoord center;

    Matrix<3, 3> Kernel(SphereCoord x, SphereCoord y)
    {
        Matrix<3, 3> temp;


        //convert to the argument of the stokeslet.
        vec<3> r = RecttoVec(x - y);
        double R = sqrt(dot(r, r));
        double rn = dot(r,RecttoVec(e_r(x)));
        //iterate over entries over the tensor.
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                //formula for the i,jth entry of the stokeslet
                temp[i][j] = r[i] * r[j] * rn / pow(R , 5);
                temp[i][j] *= - 0.75 / PI;
            }

        return temp;
    }
public:

    StokesTractionDirect(SphericalVectorField* rho, int Nn, RectCoord c = RectCoord(), std::string na = "Traction Direct")
    {
        N = Nn;
        u = rho;
        name = na;
        center = c;
    }


    RectCoord operator()(SphereCoord x)
    {



        vec<3> total = 0.0;
        //loop over trapezoidal(outside) and legendre(inside) nodes.
        for (int p = 0; p < N; p++)
            for (int i = 0; i < NUMGLNODES; i++)
            {
                //generate point on unit sphere, note theta is the first argument, phi the second.
                SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)N);
                //convert to spherical coordinate.
                SphereCoord y(s , center);
                //add contribution at y to integral.
                total = total + dot(Kernel(x, y), RecttoVec((*u)(y))) * sin(s.theta) * GLweights[i] * (PI / N) * PI;
                if (isnan(total[0]) || isnan(total[1]) || isnan(total[2]))
                {
                    std::cout << "nan detected from function " << u->name << " at " << y <<"\n";
                    std::cout << u->name << "(y) = " << ((*u)(y)) << "\n";
                }
            }

        return VectoRect(total);
    }
};


class StokesTractionHarmonic : public SphericalVectorField
{
private:
    SingleLayerHarmonic* slh;
    StokesPressure* sp;


public:

    StokesTractionHarmonic(SingleLayerHarmonic* u, StokesPressure* p ,std::string na = "Traction Harmonic")
    {
        slh = u;
        sp = p;
        name = na;
    }

    RectCoord operator()(SphereCoord s)
    {
        
        SphereCoord temp = recenter(s , slh->center);

       
        double tr = 2.0 * slh->eR_dR(s) -(*sp)(s);
        double ttheta = slh->eTheta_dR(s) + slh->eR_dTheta(s) / temp.rho - slh->eTheta(s) / temp.rho;
        double tphi = slh->ePhi_dR(s) + slh->eR_dPhi(s) / (temp.rho * sin(temp.s.theta)) - slh->ePhi(s) / temp.rho;

        return tr * e_r(s) + ttheta * e_theta(s) + tphi * e_phi(s);

    }

};

class StokesPressurefromTraction : public SphericalScalarFunction
{
private:
    SingleLayerHarmonic* slh;
    StokesTractionDirect* t;

public:
    StokesPressurefromTraction(SingleLayerHarmonic* slh0 , StokesTractionDirect* t0): slh(slh0) , t(t0){}

    double operator()(SphereCoord x)
    {
        return 2.0 * slh->eR_dR(x) - dot((*t)(x) , e_r(x));
    }
};