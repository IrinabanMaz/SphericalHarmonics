// SphericalHarmonics.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "StokesOperators.h"




typedef std::vector<std::vector<double>> bundle;

class TestSeriesApprox : public SphericalVectorField
{
private:
    VReal vr;
    VImag vi;
    WReal wr;
    WImag wi;
    XReal xr;
    XImag xi;

    
        
    int N;
    void generateCoefficients(bundle & coef , int N)
    {
        std::default_random_engine gen;
        std::uniform_real_distribution<double> r(0.0, 10.0);

        coef.resize(N + 1);
        for (int i = 0; i < N + 1; i++)
        {
            coef[i].resize(2 * i + 1); \
            for (int j = 0; j < 2 * i + 1; j++)
                coef[i][j] = r(gen) * exp(-i);
        }
    }

public:
    bundle VRealcoefs, VImagcoefs,
        WRealcoefs, WImagcoefs,
        XRealcoefs, XImagcoefs;

    TestSeriesApprox(int n)
    {
        N = n;
        generateCoefficients(VRealcoefs, N);
        generateCoefficients(VImagcoefs, N);
        generateCoefficients(WRealcoefs, N);
        generateCoefficients(WImagcoefs, N);
        generateCoefficients(XRealcoefs, N);
        generateCoefficients(XImagcoefs, N);
        for (int n = 0; n <= N; n++)
        {
            VImagcoefs[n][0] = 0.0;
            WImagcoefs[n][0] = 0.0;
            XRealcoefs[n][0] = 0.0;
            XImagcoefs[n][0] = 0.0;
        }
    }

    SphereCoord operator()(SphereCoord x)
    {
        SphereCoord temp;
        for (int i = 0; i <= N; i++)
        {
            
            for (int j = 0; j <= i; j++)
            {

                vr.reset(j, i);
                vi.reset(j, i);
                wr.reset(j, i);
                wi.reset(j, i);
                xr.reset(j, i);
                xi.reset(j, i);

                if (i == 0)
                {
                    temp = temp + VRealcoefs[i][j] * vr(x) + VImagcoefs[i][j] * vi(x);
                }
                else 
                {
                    temp = temp + VRealcoefs[i][j] * vr(x)  + VImagcoefs[i][j] * vi(x) 
                        + WRealcoefs[i][j]* wr(x)  + WImagcoefs[i][j] * wi(x) 
                        + XRealcoefs[i][j]* xr(x)  + XImagcoefs[i][j] * xi(x) ;
                }
                
            }
        }
       
        return temp;
    }
};

int main()
{   


    
    //numerical parameters.
    
    //series truncation
    const int N = 5;
    

    //number of panel points per dimension in integrals.
    const int NUMGRIDS = 20;

    std::string Zcoords;
    std::string Ycoords;
    std::string UZcoords;
    std::string VYcoords;
    std::string Rcoords;

    std::string fileoutput;

    std::fstream file;
    file.open("plotdata.py");

    TestSeriesApprox test(N);
    /*
    SphericalVectorField test([](SphereCoord x)
        {
            return VReal(2, 3)(x) + 2.4 * VImag(1, 4)(x)
                + 12.2 * WReal(3, 3)(x) + 5.5 * WImag(3, 3)(x)
                + 18.0 * XReal(2, 5)(x) + 3.3 * XImag(4, 4)(x);
        });
     */
    VSHSeries testapprox(N, NUMGRIDS);

    testapprox.approximate(&test);

    std::cout << "Max Error in approximating by VSH: " << LInfDifference(&test, &testapprox, NUMGRIDS) << std::endl;

    std::cout << "Errors in approximating fourier coefficients for V: " << std::endl;
    
   

    for (int n = 0; n <= N; n++)
        for (int m = 0; m <= n; m++)
            std::cout << " m = " << m << ", n = " << n << std::endl
            << "Real part: " << testapprox.terms[n][m].rhohatVr - test.VRealcoefs[n][m] << std::endl
            << "Imaginary part: " << testapprox.terms[n][m].rhohatVr - test.VRealcoefs[n][m] << std::endl;

    std::cout << "Errors in approximating fourier coefficients for W: " << std::endl;

    
    for (int n = 0; n <= N; n++)
        for (int m = 0; m <= n; m++)
            std::cout << " m = " << m << ", n = " << n << std::endl
            << "Real part: " << testapprox.terms[n][m].rhohatWr - test.WRealcoefs[n][m]  << std::endl
            << "Imaginary part: " << testapprox.terms[n][m].rhohatWi -test.WImagcoefs[n][m]<< std::endl;

    std::cout << "Errors in approximating fourier coefficients for X: " << std::endl;

    
    for (int n = 0; n <= N; n++)
        for (int m = 0; m <= n; m++)
            std::cout << " m = " << m << ", n = " << n << std::endl
            << "Real part: " << testapprox.terms[n][m].rhohatXr - test.XRealcoefs[n][m]<< std::endl
            << "Imaginary part: " << testapprox.terms[n][m].rhohatXi - test.XRealcoefs[n][m]<< std::endl;



    SingleLayerDirect directSLapprox(&test, NUMGRIDS);
    SingleLayerHarmonic harmonicSL(testapprox);
    StokesPressure harmonicP(testapprox);

    std::cout << "Error in applying single layer operator using both methods: "
        << L2Difference(&directSLapprox, &harmonicSL, N, 2.0) << std::endl;

     /*
    SphericalVectorField boundaryData([](SphereCoord s)
        {return cos(s.s.theta) * e_r(s) - sin(s.s.theta) * e_theta(s); });

    SphericalVectorField TrueSolution([](SphereCoord s)
        {
            SphereCoord temp;
            temp = RectToSphere(e_r(s)* cos(s.s.theta) * (3.0 / (2.0 * s.rho) - 1.0 / (2.0 * s.rho * s.rho * s.rho)));
            temp = temp - RectToSphere(sin(s.s.theta) * (3.0 / (4.0 * s.rho) + 1.0 / (4.0 * s.rho * s.rho * s.rho)) * e_theta(s));
            return 2.0 / 3.0 *temp;

        });
        
    VSHSeries bDApprox(N, NUMGRIDS);

    bDApprox.approximate(&boundaryData);

    SingleLayerHarmonic approxSolution(bDApprox);

    std::cout << L2Difference(&approxSolution, &TrueSolution, NUMGRIDS, 1.0);
   */

    StokesTractionDirect tractionDirect(&test , NUMGRIDS);
    StokesTractionHarmonic tractionHarmonic(&harmonicSL, &harmonicP);

    SphericalVectorField zero([](SphereCoord x) { return SphereCoord(0.0); });

    //std::cout << "Norm in computing traction directly: " << L2Difference(&tractionDirect, &zero, NUMGRIDS, 2.0) << std::endl;
    std::cout << "Norm in computing traction with harmonics: " << L2Difference(&tractionHarmonic, &zero, NUMGRIDS, 2.0) << std::endl;


    /*
    Zcoords.append("Z = np.linspace(-4.0 , 4.0 , 400) ");
    Ycoords.append("Y = np.linspace(-4.0 , 4.0 , 400) ");
    UZcoords.append("UZ = [[");
    VYcoords.append("UY = [[");
    Rcoords.append("R = [[");

    for(int i = 0; i < 400; i++)
        for (int j = 0; j < 400; j++)
        {
            double zcoord = -4.0 + .02 * j;
            double ycoord = -4.0 + .02 * i;
            RectCoord evaluation;
            SphereCoord evaluationr;
            if (zcoord * zcoord + ycoord * ycoord > 1)
            {
                RectCoord point(0, ycoord, zcoord);
                SphereCoord spoint = RectToSphere(point);
                evaluationr = approxSolution(spoint);
                evaluation = RectCoord(evaluationr);
            }
            else
            {
                
                evaluation = RectCoord(0.0, 0.0, 1.0);
                evaluationr = RectToSphere(evaluation);
            }
           if (j < 399)
           {
                    //Zcoords.append(std::to_string(zcoord) + ", ");
                    //Ycoords.append(std::to_string(ycoord) + ", ");
                    UZcoords.append(std::to_string(evaluation.z - 1.0) + ", ");
                    VYcoords.append(std::to_string(evaluation.y) + ", ");
                    Rcoords.append(std::to_string(evaluationr.rho) + ", ");
           }
           else if (i < 399)
           {
                    //Zcoords.append(std::to_string(zcoord) + ", ");
                    //Ycoords.append(std::to_string(ycoord) + ", ");
                    UZcoords.append(std::to_string(evaluation.z) + "], [");
                    VYcoords.append(std::to_string(evaluation.y) + "], [");
                    Rcoords.append(std::to_string(evaluationr.rho) + "], [");
           }
           else
           {
                    //Zcoords.append(std::to_string(zcoord) + "]");
                    //Ycoords.append(std::to_string(ycoord) + "]");
                    UZcoords.append(std::to_string(evaluation.z) + "]]");
                    VYcoords.append(std::to_string(evaluation.y) + "]]");
                    Rcoords.append(std::to_string(evaluationr.rho) + "]]");
           }
            
        }

    fileoutput = Zcoords + "\n" + Ycoords + "\n" + UZcoords + "\n" + VYcoords + "\n" + Rcoords;

    file << fileoutput;
   */
    return 0;

    
}

