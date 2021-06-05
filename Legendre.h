#include <vector>
#include<cmath>

#pragma once
int factorial(int n)
{
    int ret = 1;
    if (n < 0)
    {
        std::cout << "Error: negative integer passed to factorial";
        return -1;
    }
    else for (int i = 2; i <= n; i++)
        ret *= i;

    return ret;
}

//Class which computes the value of the Associated Legendre Functions using recurrence relations. The class stores all intermediate values (over m,n)for a fixed input
// for potential later use.
class Legendre
{


private:
    //the value of the input variable to the function.
    double input = 100;
    //the maximum value of the parameter m computed for this value of the input.
    int maxorder = 0;
    //the maximum value of the degree n computed for this value of the input.
    int maxdegree = 0;
    std::vector<std::vector<double>> values;
    ;

public:

    Legendre()
    {
        input = 100;
        maxorder = -1;
        maxdegree = -1;
    }

    ///Populate the values of the legendre polynomials at the gicen point values.
    ///-double val: The point we are evaluating at.
    ///-int M: the maximum order of the legendre polynomials
    ///int L: the maximum degree of the polynomials.
    void populate(double val, int M, int L)
    {


        input = val;
        maxorder = abs(M);
        maxdegree = L;
        double s = sqrt(1 - input * input) + 10e-12;

        


        values.resize(abs(M) + 1);

        for (int m = 0; m < abs(M) + 1; m++)
            values[m].resize(std::max(L + 1 , 2));
        
        //set values in case of boundary condition.

        //set initial values.
        values[0][0] = 1.0;
        values[0][1] = input;

        //set values for m = 0 using recurrence in terms of two previous i values.
        for (int i = 1; i < L; i++)
            values[0][i + 1] = ((2.0 * (double)i + 1) * input * values[0][i] - (double)i * values[0][i - 1]) / ((double)i + 1);


        //computes values for larger values of m in terms of previous value of m. for i = 0, uses the values from i = 1, i = 2, instead.
        for (int m = 0; m < abs(M); m++)
        {

            for (int i = 1; i < L + 1; i++)
                values[m + 1][i] = ((((double)i - (double)m)) * input * values[m][i] - ((double)i + (double)m) * values[m][i - 1]) / s;

            values[m + 1][0] = ((-(double)m + 1.0) * values[m][1] + -((double)m + 1.0) * values[m][0]) / s;

        }


    }

    //Accesses the value of the legendre functions for specified values of m,i.
    double getValue(int m, int i)
    {
        if (m >= 0)
            return values[m][i];
        else 
            return ((m % 2 == 0) ? 1.0 : -1.0) * (double)factorial(i + m) / (double)factorial(i - m) * values[-m][i];

        
    }

};
double testpoly(double x)
{
    double xsq = 1 - x * x;
    xsq = pow(xsq, 2.5);
    return -10395.0 *x * xsq;
}

double testerr()
{
    int m = 5;
    int l = 6;

    Legendre P;
    
    double GLnodes[16] = { -0.0950125098376374  , 0.0950125098376374  , -0.2816035507792589 , 0.2816035507792589 ,  -0.4580167776572274  , 0.4580167776572274
                          - 0.6178762444026438  , 0.6178762444026438 , -0.7554044083550030 , 0.7554044083550030 , -0.8656312023878318 , 0.8656312023878318  ,
                          -0.9445750230732326 , 0.9445750230732326  , -0.9894009349916499  , 0.9894009349916499 };

    double GLweights[16] = { 0.1894506104550685 , 0.1894506104550685 , 0.1826034150449236  , 0.1826034150449236  , 0.1691565193950025  , 0.1691565193950025 ,
                            0.1495959888165767  , 0.1495959888165767 , 0.1246289712555339  , 0.1246289712555339  , 0.0951585116824928  , 0.0951585116824928 ,
                            0.0622535239386479  , 0.0622535239386479 , 0.0271524594117541  , 0.0271524594117541 };
    double total = 0.0;
    for (int i = 0; i < 16; i++)
    {
        P.populate(GLnodes[i], m, l);
        double err = P.getValue(m, l) - testpoly(GLnodes[i]);
        total += GLweights[i] * err * err;
    }

    return total;
}