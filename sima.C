#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

double sim_anneal(double t0, int num_temperatures, int num_iterations,
		  double start_value, vector<double>* chain);
double cauchy_llh(double alpha);
double std_normal();
  
int main()
{
  srandom(time(0));  // Seed the random number generator

  cout << "Enter t0, the starting temperature: ";
  double t0;
  cin >> t0;

  cout << "Enter num_temperatures, the number of temperatures: ";
  int num_temperatures;
  cin >> num_temperatures;

  cout << "Enter num_iterations, the number of iterations per temperature: ";
  int num_iterations;
  cin >> num_iterations;

  cout << "Enter start_value, the starting value of the time 0 chain: ";
  double start_value;
  cin >> start_value;

  vector<double> chain = vector<double>(num_iterations);
  double min = sim_anneal(
      t0, num_temperatures, num_iterations, start_value, &chain);

  /*
    Uncomment if printing out the last iteration
    for (int i = 0; i < num_iterations; ++i)
    cout << chain[i] << endl;
  */

  cout << "The estimate of the minimum is: " << min << endl;

  return 0;
}

double sim_anneal(double t0, int num_temperatures, int num_iterations,
		  double start_value, vector<double>* chain) {
  double oldval = start_value;
  for (int n = 0; n < num_temperatures; ++n) {
    double currtemp = t0 / (n + 1);
    for (int iter = 0; iter < num_iterations; ++iter) {
      (*chain)[iter] = std_normal() * sqrt(currtemp) + oldval;
      double delta = exp(
          (cauchy_llh(oldval) - cauchy_llh((*chain)[iter])) / currtemp);
      if (delta < 1) {
        double unif = 1.0 * random() / RAND_MAX;
        if (unif > delta)
          (*chain)[iter] = oldval;
      }
      oldval = (*chain)[iter];
    }
    oldval = (*chain)[num_iterations - 1];
  }
  return (*chain)[num_iterations - 1];
}

double cauchy_llh(double alpha) {
  const double beta = .1;
  const double points[] = {-1.4, -1.1, -.8, 0, 1, 1.4, 2, 2.5};
  
  double val = 0;
  for (int i = 0; i < 8; ++i) {
    val += log(beta * beta + (points[i] - alpha) * (points[i] - alpha));
  }
  return val;
}

double std_normal() {
  static bool toggler = false;
  static double val2;

  if ((toggler = !toggler)) {
    double u = 1.0 * random() / RAND_MAX;
    double v = 1.0 * random() / RAND_MAX;
    double val1 = sqrt(-2 * log(u)) * sin(2 * M_PI * v);
    val2 = sqrt(-2 * log(u)) * cos(2 * M_PI * v);
    return val1;
  } else {
    return val2;
  }
}
