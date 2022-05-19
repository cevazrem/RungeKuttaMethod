#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;
const double eps = 0.01; // epsilon
const int N = 3; // number of equations

// ----- Function of calculation right part -----
// Y     - input :  array of variables
// M    - input : parameter
// D     - input : parameter
// FRight - output : array of right part
double* F(double* Y, double M, double D) {
    double* FRight = new double[N];
    FRight[0] = (Y[0] - D/2*Y[1]+Y[1]*(Y[2]+Y[0]*Y[0]));
    FRight[1] = (D/2*Y[0]+Y[1]+Y[0]*(3*Y[2]-Y[0]*Y[0]));
    FRight[2] = ((-2)*Y[2]*(M+Y[0]*Y[1]));

    return FRight;
}

// ----- Function Method Runge-Kytta -----
// Y     - input :  array of variables
// M    - input : parameter
// D     - input : parameter
// h     - input : step
// NextY - output : array of variables on next iteration
double* MethodRungeKytta(double M, double D, double* Y, double h) {
    double* ArrHalf = new double[N]; // additional array for intermediate calculations 
    double* NextY = new double[N];
    double* k1, * k2, * k3, * k4; // arrays of coefficients

    k1 = F(Y, M, D);

    ArrHalf[0] = Y[0] + (h / 2) * k1[0];
    ArrHalf[1] = Y[1] + (h / 2) * k1[1];
    ArrHalf[2] = Y[2] + (h / 2) * k1[2];
    k2 = F(ArrHalf, M, D);

    ArrHalf[0] = Y[0] + (h / 2) * k2[0];
    ArrHalf[1] = Y[1] + (h / 2) * k2[1];
    ArrHalf[2] = Y[2] + (h / 2) * k2[2];
    k3 = F(ArrHalf, M, D);

    ArrHalf[0] = Y[0] + h * k3[0];
    ArrHalf[1] = Y[1] + h * k3[1];
    ArrHalf[2] = Y[2] + h * k3[2];
    k4 = F(ArrHalf, M, D);

    NextY[0] = Y[0] + (h / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    NextY[1] = Y[1] + (h / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    NextY[2] = Y[2] + (h / 6) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);

    return NextY;
}


// ----- Function Calculate some same things -----
// M    - input : parameter
// A     - input : parameter
double calcRoot(double M, double D) {
	double result = 0;
	result = (-D*M + sqrt(D*D*M*M+3*M*M-4*M*M*M))/(3-4*M);
	return result;
}

int main()
{
    double* Y = new double[N]; // array of variables

    double D; // parameter
    double M; // parameter
    double hStep = 0.01; // step
    double tMax;

    cout << "D: ";
    cin >> D;
    cout << "M: ";
    cin >> M;
    cout << "tMax: ";
    cin >> tMax;

    ofstream result("Result.txt"); // open file for recording
    // check file opening and print error
    if (!result.is_open())
    {
        cout << "Oops! Couldn't open the file, bro... :с" << endl;
    }

    // setting the initial conditions
    Y[0] = sqrt(calcRoot(M, D))+ eps; 
    cout << Y[0] << endl;
    Y[1] = (-M)/sqrt(calcRoot(M, D)) + eps;
    cout << Y[1] << endl;
    Y[2] = D/2 + calcRoot(M, D)*(1/M - 1)+ eps;
    cout << Y[2];
    result << "# D = " << D << " M = " << M << " Eps = " << eps << endl;
    result << "# Time" << ' ' << "Y_1" << ' ' << "Y_2" << ' ' << "Y_3" << ' ' << "Y_2-Y_1" << ' ' << "Y_3-Y_1" << ' ' << "Y_3-Y_2" << ' ' << "Y_3-Y_2-Y_1" << endl;

    double t = 0;
    while (t < tMax) {
        for (int j = 0; j < 1; j++) {
            Y = MethodRungeKytta(M, D, Y, hStep);
            t += hStep;
        }
        
        result << t << ' ' << Y[0] << ' ' << Y[1] << ' ' << Y[2] << ' ' << Y[1]-Y[0] << ' ' << Y[2]-Y[0] << ' ' << Y[2]-Y[1] << ' ' << Y[2]-Y[1]-Y[0] << endl;
    }

    return 0;
}
