#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>
#include <random>
#include <string>
#include <algorithm>
#include <cassert>

using namespace std;

// Constants for the simulation
const double MASS = 2.0;                      // [pg] Mass of the particle
const double TOTAL_TIME = 1e10;               // [µs] Total simulation time
const long int NUM_STEPS = 1e12;             // Number of simulation steps
const double DAMPING_PARAM = 14.0;           // [pN·µs/nm] Damping parameter
const double KBT = 4.1;                      // [pN·nm] Boltzmann constant times temperature
const double TIME_STEP_INIT = TOTAL_TIME / NUM_STEPS; // Initial time step
const double DIM_C = 0.5 * DAMPING_PARAM * TIME_STEP_INIT / MASS;
const double DIM_B = 1.0 / (1.0 + DIM_C);
const double DIM_A = DIM_B * (1.0 - DIM_C);
const double PRECISION = 0.0001;             // Precision for Newton's method

// Polynomial coefficients for potential energy
const double POLY_COEFFS[5] = {64.0, 0.0, -32.0, 0.0, 4.0};
const double X_LEFT = -4.0;                  // Left boundary of the potential
const double X_RIGHT = 4.0;                  // Right boundary of the potential
const int NUM_DIVISIONS = 1000;              // Number of divisions for numerical integration
const double DX = (X_RIGHT - X_LEFT) / NUM_DIVISIONS; // Step size for numerical integration

// Random number generators setup
mt19937 rng(random_device{}());
auto normal_dist = bind(normal_distribution<double>(0, 1), rng);
auto uniform_dist = bind(uniform_real_distribution<double>(0, 1), rng);
auto velocity_dist = bind(normal_distribution<double>(0, sqrt(KBT / MASS)), rng);

/**
 * Computes the external pulling force at a given time.
 * @param time Current simulation time.
 * @param pulling_rate The rate of pulling force [pN/µs].
 * @return The external force value.
 */
double external_force(double time, double pulling_rate) {
    return pulling_rate * time;
}

/**
 * Computes the kinetic energy of a particle.
 * @param velocity Current velocity of the particle.
 * @return Kinetic energy value.
 */
double kinetic_energy(double velocity) {
    return 0.5 * MASS * velocity * velocity;
}

/**
 * Computes the potential energy of the particle at a given position and time.
 * @param position Current position of the particle.
 * @param time Current simulation time.
 * @param pulling_rate The rate of pulling force [pN/µs].
 * @return Potential energy value.
 */
double potential_energy(double position, double time, double pulling_rate) {
    return POLY_COEFFS[0] + position * (POLY_COEFFS[1] - external_force(time, pulling_rate) + position * (POLY_COEFFS[2] + position * (POLY_COEFFS[3] + position * POLY_COEFFS[4])));
}

/**
 * Computes the first derivative of the potential energy with respect to position.
 * @param position Current position of the particle.
 * @param time Current simulation time.
 * @param pulling_rate The rate of pulling force [pN/µs].
 * @return Derivative of the potential energy.
 */
double potential_derivative(double position, double time, double pulling_rate) {
    return POLY_COEFFS[1] - external_force(time, pulling_rate) + position * (2 * POLY_COEFFS[2] + position * (3 * POLY_COEFFS[3] + position * 4 * POLY_COEFFS[4]));
}

/**
 * Computes the second derivative of the potential energy with respect to position.
 * @param position Current position of the particle.
 * @return Second derivative of the potential energy.
 */
double potential_second_derivative(double position) {
    return 2 * POLY_COEFFS[2] + position * (6 * POLY_COEFFS[3] + position * 12 * POLY_COEFFS[4]);
}

/**
 * Integration function for computing the partition function.
 * @param position Current position.
 * @param pulling_rate The rate of pulling force [pN/µs].
 * @return Value of the integrand.
 */
double integration_function(double position, double pulling_rate) {
    return exp(-potential_energy(position, 0, pulling_rate) / KBT);
}

/**
 * Computes the partition function using the trapezoidal rule.
 * @param pulling_rate The rate of pulling force [pN/µs].
 * @return Partition function value.
 */
double trapezoidal_integration(double pulling_rate) {
    const size_t num_points = 500;
    double step_size = (X_RIGHT - X_LEFT) / (num_points - 1);
    double value = 0.5 * (integration_function(X_LEFT, pulling_rate) + integration_function(X_RIGHT, pulling_rate));

    for (size_t i = 1; i < num_points - 1; ++i) {
        value += integration_function(X_LEFT + i * step_size, pulling_rate);
    }

    return step_size * value;
}

/**
 * Computes the number density of the particle.
 * @param position Current position.
 * @param pulling_rate The rate of pulling force [pN/µs].
 * @return Number density value.
 */
double number_density(double position, double pulling_rate) {
    double normalization_factor = 1.0 / trapezoidal_integration(pulling_rate);
    return normalization_factor * integration_function(position, pulling_rate);
}

/**
 * Solves for the roots of the potential derivative using Newton's method.
 * @param position Reference to the position to adjust.
 * @param time Current simulation time.
 * @param pulling_rate The rate of pulling force [pN/µs].
 */
void newton_method(double &position, double time, double pulling_rate) {
    double delta;
    do {
        delta = potential_derivative(position, time, pulling_rate) / potential_second_derivative(position);
        position -= delta;
    } while (fabs(delta) >= PRECISION);
}

/**
 * Updates the position and velocity using the modified Verlet algorithm.
 * @param position Reference to the current position.
 * @param velocity Reference to the current velocity.
 * @param time Current simulation time.
 * @param time_step Reference to the current time step.
 * @param pulling_rate The rate of pulling force [pN/µs].
 */
void update_modified_verlet(double &position, double &velocity, double time, double &time_step, double pulling_rate) {
    double random_force = sqrt(2 * KBT * DAMPING_PARAM * time_step) * normal_dist();
    double force1 = -potential_derivative(position, time, pulling_rate);
    position += (velocity + 0.5 * (force1 * time_step + random_force) / MASS) * DIM_B * time_step;
    double force2 = -potential_derivative(position, time, pulling_rate);
    velocity = DIM_A * velocity + 0.5 * time_step / MASS * (DIM_A * force1 + force2) + DIM_B / MASS * random_force;
}

/**
 * Finds the extrema of the potential energy function.
 * @param xl Reference to the left minimum.
 * @param xb Reference to the barrier position.
 * @param xr Reference to the right minimum.
 * @param time Current simulation time.
 * @param pulling_rate The rate of pulling force [pN/µs].
 */
void find_extrema(double &xl, double &xb, double &xr, double time, double pulling_rate) {
    int j = 0;
    assert(potential_derivative(X_LEFT, time, pulling_rate) < 0);

    for (; j < NUM_DIVISIONS; ++j) {
        double position = X_LEFT + j * DX;
        if (potential_derivative(position, time, pulling_rate) > 0) {
            xl = position;
            break;
        }
    }

    for (; j < NUM_DIVISIONS; ++j) {
        double position = X_LEFT + j * DX;
        if (potential_derivative(position, time, pulling_rate) < 0) {
            xb = position;
            break;
        }
    }

    for (; j < NUM_DIVISIONS; ++j) {
        double position = X_LEFT + j * DX;
        if (potential_derivative(position, time, pulling_rate) > 0) {
            xr = position;
            break;
        }
    }

    newton_method(xl, time, pulling_rate);
    newton_method(xb, time, pulling_rate);
    newton_method(xr, time, pulling_rate);
}

/**
 * Generates a random initial position based on the number density.
 * @param pulling_rate The rate of pulling force [pN/µs].
 * @return Initial position value.
 */
double random_initial_position(double pulling_rate) {
    double random_value = uniform_dist();
    double cumulative = 0.0;
    const double delta_x = 0.05;

    for (size_t i = 0; i <= 500; ++i) {
        double position = X_LEFT + i * delta_x;
        cumulative += delta_x * number_density(position, pulling_rate);
        if (cumulative > random_value) {
            return position;
        }
    }

    assert(false);
    return -99999.0;
}

/**
 * Main simulation program.
 * Runs simulations for multiple pulling rates and outputs results to a file.
 */
int main() {
    vector<double> pulling_rates = {0.0001, 0.001, 0.01, 0.1};
    ofstream output_file("results.dat", ios::app);
    output_file << "Run_Num Pulling_Rate External_Force" << endl;

    for (double pulling_rate : pulling_rates) {
        for (int run = 1; run <= 5; ++run) {
            double xl, xb, xr;
            find_extrema(xl, xb, xr, 0, pulling_rate);

            double initial_position = random_initial_position(pulling_rate);
            while (initial_position >= 0.5 * (xl + xb)) {
                initial_position = random_initial_position(pulling_rate);
            }

            double initial_velocity = velocity_dist();
            double time = 0;
            double position = initial_position;
            double velocity = initial_velocity;
            double time_step = TIME_STEP_INIT;

            while (time < TOTAL_TIME) {
                time += time_step;
                update_modified_verlet(position, velocity, time, time_step, pulling_rate);

                if (position > 0.001 * (xl + xb)) {
                    time_step = min(time_step, DX / fabs(velocity));
                }

                newton_method(xl, time, pulling_rate);
                newton_method(xb, time, pulling_rate);
                newton_method(xr, time, pulling_rate);

                double potential_barrier = potential_energy(xb, time, pulling_rate) - max(potential_energy(xl, time, pulling_rate), potential_energy(xr, time, pulling_rate));
                if (potential_barrier < 1E-20) {
                    break;
                }

                if (position > xb) {
                    output_file << run << " " << pulling_rate << " " << external_force(time, pulling_rate) << endl;
                    break;
                }
            }
        }
    }

    output_file.close();
    return 0;
}

