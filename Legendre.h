#pragma once
#include <vector>
#include<cmath>
#include<iostream>

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
    inline static double input;
    //the maximum value of the parameter m computed for this value of the input.
    inline static int maxorder;
    //the maximum value of the degree n computed for this value of the input.
    inline static int maxdegree;
    inline static std::vector<std::vector<double>> values;
    

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

        if (abs(val - input) < 1e-16)
            if (maxdegree > L)
                if (maxorder > L)
                    return;

        input = val;
        maxorder = std::max(L, 2);
        maxdegree = std::max(L , 2);
        double s = sqrt(1 - input * input) + 1e-16;

        


        values.resize(maxdegree + 1);

        for (int l = 0; l <=maxdegree; l++)
            values[l].resize(2*l + 1);
        
        //set values in case of boundary condition.

        //set initial values.
        values[0][0] = 1.0;
        values[1][0] = input;

        //set values for m = 0 using recurrence in terms of two previous i values.
        for (int i = 1; i < L; i++)
            values[i + 1][0] = ((2.0 * (double)i + 1) * input * values[i][0] - (double)i * values[i - 1][0]) / ((double)i + 1);


        //computes values for larger values of m in terms of previous value of m. for i = 0, uses the values from i = 1, i = 2, instead.
        for (int i = 1; i <=L; i++)
           for (int m = 0; m <i; m++)
                values[i][m + 1] = ((((double)i - (double)m)) * input * values[i][m] - ((double)i + (double)m) * values[i - 1][m]) / s;

            

        
        
    }

    //Accesses the value of the legendre functions for specified values of m,i.
    double getValue(int m, int i)
    {
        if (m >= 0)
            return values[i][m];
        else 
            return ((m % 2 == 0) ? 1.0 : -1.0) * (double)factorial(i + m) / (double)factorial(i - m) * values[i][-m];

        
    }

    double operator()(double x, int m, int l)
    {
        populate(x, m, l);
        return getValue(m, l);
    }

    double dtheta(int m, int n, double theta, int order = 1)
    {
        populate(cos(theta), m + 1, n + 1);

        double temp = 0.0;
        if (order == 0)
            return getValue(m, n);
        else if(order == 1)
        {
            if (n == 0)
                temp = 0;
            else if (m >= 0)
                temp = -(double)((n + m) * (n - m + 1)) * getValue(m - 1, n)
                - (double)m * cos(theta) / sin(theta) * getValue(m, n);
            else
                temp = getValue(m + 1, n) + (double)m * cos(theta) / sin(theta) * getValue(m, n);


            return temp;
        }
        else if (order == 2)
        {
            if (n == 0)
                return 0;
            else if (m >= 0)
            {
                return -(double)(n + m) * (n - m + 1) * dtheta(m - 1, n, theta, 1) + (double)m * cos(theta) / (sin(theta) * sin(theta)) * getValue(m, n)
                    - (double)m * cos(theta) / sin(theta) * dtheta(m, n, theta, 1);
            }
            else
            {
                return -(double)(n + m) * (n - m + 1) * dtheta(m + 1, n, theta, 1) - (double)m * cos(theta) / (sin(theta) * sin(theta)) * getValue(m, n)
                    + (double)m * cos(theta) / sin(theta) * dtheta(m, n, theta, 1);
            }
        }

    }


};
double testpoly(double x)
{
    double s = sqrt(1 - x * x);
    return 1.0;
}

double testerr()
{
    int m = 0;
    int l = 0;

    Legendre P;
    
    double GLnodes[16] = { -0.0950125098376374  , 0.0950125098376374  , -0.2816035507792589 , 0.2816035507792589 ,  -0.4580167776572274  , 0.4580167776572274,
                          - 0.6178762444026438  , 0.6178762444026438 , -0.7554044083550030 , 0.7554044083550030 , -0.8656312023878318 , 0.8656312023878318  ,
                          -0.9445750230732326 , 0.9445750230732326  , -0.9894009349916499  , 0.9894009349916499 };

    double GLweights[16] = { 0.1894506104550685 , 0.1894506104550685 , 0.1826034150449236  , 0.1826034150449236  , 0.1691565193950025  , 0.1691565193950025 ,
                            0.1495959888165767  , 0.1495959888165767 , 0.1246289712555339  , 0.1246289712555339  , 0.0951585116824928  , 0.0951585116824928 ,
                            0.0622535239386479  , 0.0622535239386479 , 0.0271524594117541  , 0.0271524594117541 };
    double total = 0.0;
    for (int i = 0; i < 16; i++)
    {
        
        double err = P(GLnodes[i],m, l) - testpoly(GLnodes[i]);
        total += GLweights[i] * err * err;
    }

    return total;
}
