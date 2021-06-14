#include <iostream>

/**
 * Program demonstrates the Runge phenomenon in Lagrange polynomial interpolation.
 * Interpolated function: f(x) = 1/(1+10x^6).
 * Function defined on the interval [-1, 1].
 * The construction of interpolating polynomials uses a Newton's basis with equidistant nodes and Chebyshev nodes.
 * */

const double PI = 3.1415926535897932384;
const int N = 10; // nr of nodes
const double X_MIN =  -1.0;
const double X_MAX =  1.0;

double interpolated_func(double x);
void interpolation_polynomials_construction();

int main()
{
    interpolation_polynomials_construction();

    return 0;
}

void interpolation_polynomials_construction()
{
    std::cout << PI;
}

double interpolated_func(double x)
{
    return 1.0 / (1.0 + 10.0 * x * x * x * x * x * x);
}
