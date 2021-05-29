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
        double s = sqrt(1 - val * val) + 10e-12;

        


        values.resize(abs(M) + 1);

        for (int m = 0; m < abs(M) + 1; m++)
            values[m].resize(std::max(L + 1 , 2));
        
        //set values in case of boundary condition.
        if (1.0 - val < 10e-6)
            for (int i = 0; i < L; i++)
                for (int m = 0; m < M; m++)
                    values[m][i] = 1.0;

        //set initial values.
        values[0][0] = 1;
        values[0][1] = input;

        //set values for m = 0 using recurrence in terms of two previous i values.
        for (int i = 1; i < L; i++)
            values[0][i + 1] = ((2 * (double)i + 1) * input * values[0][i] - i * values[0][i - 1]) / ((double)i + 1);


        //computes values for larger values of m in terms of previous value of m. for i = 0, uses the values from i = 1, i = 2, instead.
        for (int m = 0; m < abs(M); m++)
        {

            for (int i = 1; i < L + 1; i++)
                values[m + 1][i] = ((((double)i - (double)m)) * input * values[m][i] - ((double)i + (double)m) * values[m][i - 1]) / s;

            values[m + 1][0] = ((-(double)m + 1) * values[m][1] + -((double)m + 1) * values[m][0]) / s;

        }


    }

    //Accesses the value of the legendre functions for specified values of m,i.
    double getValue(int m, int i)
    {
        if (m >= 0)
            return values[m][i];
        else if (m < 0)
            return ((m % 2 == 0) ? 1 : -1) * (double)factorial(i + m) / factorial(i - m) * values[-m][i];

        else return -100;
    }

};
double testpoly(double x)
{
    double xsq = 1 - x * x;
    return 105 * xsq * xsq;
}