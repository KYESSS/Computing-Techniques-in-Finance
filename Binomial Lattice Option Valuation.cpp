// Tian (1993) based on Espen Haug
 
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <iostream>
#include <vector>
#include<iomanip>
#include <string>
 
#include <time.h>
 
//using namespace System;
using namespace std;
 
// Function for binomial tree
double Binomial(int n, double S, double K, double r, double q, double v, double T, char PutCall, char OpStyle) {
	int i, j;
 
	double dt, u, d, p;
	// new for Tian
	double M, H;
	int z;
 
	// Quantities for the tree
	dt = T / n;
	
	M = exp((r-q)*dt);
	H = exp(v*v*dt);
	u = M*H / 2 * ((H + 1) + sqrt(H*H + 2 * H - 3));
	d = M*H / 2 * ((H + 1) - sqrt(H*H + 2 * H - 3));
	p = (M - d) / (u - d);
 
	if (PutCall == 'C')
	{
		z = 1;
	}
	else if (PutCall == 'P')
	{
		z = -1;
	}
 
	vector<double> OptionValue;
 
	//resize the column
	OptionValue.resize(n + 1);
 
	for (i = 0; i <= n; i++) {
 
	   OptionValue[i] = max(z*(S*pow(u, i)*pow(d, n - i) - K), 0.0);
	}
 
	// Backward recursion through the tree
	for (j = n - 1; j >= 0; j--)
		for (i = 0; i <= j; i++) {
			if (OpStyle == 'E')
				OptionValue[i] = exp(-r*dt)*(p*(OptionValue[i + 1]) + (1.0 - p)*(OptionValue[i]));
			else {
 
				OptionValue[i] = max(z*(S*pow(u, i)*pow(d, j - i) - K), exp(-r*dt)*(p*(OptionValue[i + 1]) + (1.0 - p)*(OptionValue[i])));
				//OptionValue[i] = max(z*(S*pow(u, (2 * i - j)) - K), exp(-r*dt)*(p*(OptionValue[i + 1]) + (1.0 - p)*(OptionValue[i])));
			}
		}
	// Return the option price
	return OptionValue[0];
}
 
int main() {
	//double S, K, T, v, r;
	//char PutCall, OpStyle;
	//int n;
	clock_t start_time, end_time;
	start_time = clock();
 
 
	int n = 500;							// Number of steps
	double S = 100.0;						// Spot Price
	double K = 100.0;						// Strike Price
	double T = 3;							// Years to maturity
	double r = 0.03;						// Risk Free Rate
	double q = 0.07;						//double q = 			
	double v = 0.20;
	char PutCall = 'C';
	char OpStyle = 'A';
 
 
    cout << setprecision(10);
	cout << "The binomial price is " << Binomial(n, S, K, r, q, v, T, PutCall, OpStyle) << endl;
 
	end_time = clock();
	cout << " " << start_time << " " << end_time << " " << (end_time - start_time) << endl;
	cout << " " << start_time << " " << end_time << " " << (end_time - start_time) / (double)CLOCKS_PER_SEC << " seconds" << endl;
 
	//system("PAUSE");
}
