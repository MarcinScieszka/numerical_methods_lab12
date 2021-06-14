#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

/**
 * Program designed to demonstrate the Runge phenomenon in Lagrange polynomial interpolation.
 * Interpolated function: f(x) = 1/(1+10x^6).
 * Function defined on the interval [-1, 1].
 * The construction of interpolating polynomials uses a Newton's basis with equidistant nodes and Chebyshev nodes.
 * */

const double PI = 3.141592653589793238;
const double X_MIN =  -1.0; // beginning of the specified interval
const double X_MAX =  1.0; // the end of the specified interval
const double plot_h = 0.01; // step used for plotting

double interpolated_func(double x);
void interpolation_polynomials_construction();
void print_array(double *array, int size, int print_width=12);
void copy_array(double *dest_array, const double *source_array, int size);

double get_value_of_interpolation_polynomial(double x, const double *xn, const double *cn, int d);

int main()
{
    interpolation_polynomials_construction();

    return 0;
}

void interpolation_polynomials_construction()
{
    int d; // degree of the polynomial
    double interpolated_func_val; // return value of interpolated function for a given x
    double plot_x; // current position on x axis
    double step_sum; // current sum of steps
    double ksi; // root of the Chebyshev polynomial
    double h; // equidistant node step
    int i, j;

    std::cout.precision(5);
    std::cout.setf(std::ios::fixed);

//    for (d=5; d<21; d+=5)
    for (d=15; d<16; d+=1)
    {
        std::string results_Lag_New_equi_filename = std::to_string(d) + "_d_results_Lag_New_equidistant.txt";
        std::string results_Lag_New_Che_filename = std::to_string(d) + "_d_results_Lag_New_Chebyshev.txt";

        std::ofstream results_Lag_New_equi(results_Lag_New_equi_filename);
        std::ofstream results_Lag_New_Che(results_Lag_New_Che_filename);
        if(!results_Lag_New_equi || !results_Lag_New_Che)
        {
            std::cerr << "Error: file could not be opened" << std::endl;
            exit(1);
        }

        results_Lag_New_equi << "x analytical_sol interpolation\n";
        results_Lag_New_Che << "x analytical_sol interpolation\n";

        auto *xe = new double[d]; // equidistant nodes
        auto *ye = new double[d]; // f(xe)
        auto *xc = new double[d]; // Chebyshev nodes
        auto *yc = new double[d]; // f(xc)
        auto *cc = new double[d]; // Chebyshev coefficients
        auto *ce= new double[d]; // equidistant coefficients

        // determining equidistant nodes
        step_sum = 0.0;
        h = (X_MAX - X_MIN) / (d - 1.0);

        for (i = 0; i < d; i++)
        {
            xe[i] = X_MIN + step_sum;
            ye[i] = interpolated_func(xe[i]);
            step_sum += h;
        }

        // determining Chebyshev nodes
        for (i = 0; i < d; i++)
        {
            ksi = cos(((2.0*i + 1.0)/(2.0*d + 2.0)) * PI);
            xc[i] = ((X_MAX + X_MIN)/2.0) + ((X_MAX - X_MIN)/2.0)*ksi;
            yc[i] = interpolated_func(xc[i]);
        }

        print_array(xc, d);

        // determining equidistant c coefficients using difference quotients
        copy_array(ce, ye, d);
        for (i = 1; i < d; i++)
        {
            j = d - 1;
            while (j >= i)
            {
                ce[j] = (ce[j] - ce[j-1]) / (xe[j] - xe[j-i]);
                j--;
            }
        }
        print_array(ce, d);

        // determining Chebyshev c coefficients using difference quotients -
        copy_array(cc, yc, d);
        for (i = 1; i < d; i++)
        {
            j = d - 1;
            while (j >= i)
            {
                cc[j] = (cc[j] - cc[j-1]) / (xc[j] - xc[j-i]);
                j--;
            }
        }
        print_array(cc, d);

        // saving results to file
        plot_x = X_MIN;
        while (plot_x < X_MAX)
        {
            interpolated_func_val = interpolated_func(plot_x);
            results_Lag_New_equi << plot_x << " " << interpolated_func_val << " " << get_value_of_interpolation_polynomial(plot_x, xe, ce, d) << "\n";
            results_Lag_New_Che << plot_x << " " << interpolated_func_val << " " << get_value_of_interpolation_polynomial(plot_x, xc, cc, d) << "\n";
            plot_x += plot_h;
        }

        delete[] xc;
        delete[] yc;
        delete[] xe;
        delete[] ye;
        delete[] cc;
        delete[] ce;
        results_Lag_New_Che.close();
        results_Lag_New_equi.close();
    }
}

double interpolated_func(double x)
{
    return 1.0 / (1.0 + 10.0 * x * x * x * x * x * x);
}

double get_value_of_interpolation_polynomial(double x, const double *xn, const double *cn, int d)
{
    double ip_value; // value of the interpolation polynomial
    int i;

    ip_value = cn[d - 1];
    for (i = d - 1; i > 0; i--)
    {
        ip_value = ip_value * (x - xn[i - 1]) + cn[i - 1];
    }
    return ip_value;
}

void copy_array(double *dest_array, const double *source_array, int size)
{
    int i;
    for (i = 0; i < size; i++)
    {
        dest_array[i] = source_array[i];
    }
}

void print_array(double *array, int size, int print_width)
{
    int i;
    std::cout << "[ ";
    for(i=0; i<size; i++)
    {
        std::cout << array[i]  << std::setw(print_width);
    }
    std::cout << "] \n";
}