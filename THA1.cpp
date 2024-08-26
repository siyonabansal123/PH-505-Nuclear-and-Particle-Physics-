#include <iostream>
#include <cmath>
#include <random>
#include <vector>

// constants
const double V0 = 40.0;  // MeV
const double r0 = 2.1;   // fm

// Potential function
double potential(double r) {
    return -V0 * exp(-pow(r / r0, 2));
}

// Trial wavefunction (Gaussian)
double trial_wavefunction(double r, double alpha) {
    return exp(-alpha * r * r);
}

// Local energy for given r and alpha
double local_energy(double r, double alpha) {

    double kinetic = (alpha * (2 * alpha * r * r - 3)) / (2.0);
    double potential_energy = potential(r);
    return kinetic + potential_energy;
}

// Monte Carlo integration for energy calculation
double monte_carlo_integration(double alpha, int num_samples) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 5);  // Random samples from 0 to 5 fm
    
    double energy_sum = 0.0;
    
    for (int i = 0; i < num_samples; ++i) {
        double r = dis(gen);
        double wf = trial_wavefunction(r, alpha);
        double energy = local_energy(r, alpha);
        energy_sum += energy * wf * wf;  // Importance sampling
    }
    
    return energy_sum / num_samples;
}

int main() {

    const int num_samples = 100000;  // Number of Monte Carlo samples
    double alpha = 0.5;
    double min_energy = 1e10;
    double best_alpha = alpha;
    
    for (double a = 0.1; a <= 2.0; a += 0.1) {
        double energy = monte_carlo_integration(a, num_samples);
        std::cout << "Alpha: " << a << ", Energy: " << energy << std::endl;
        if (energy < min_energy) {
            min_energy = energy;
            best_alpha = a;
        }
    }
    
    std::cout << "Minimum Energy of " << min_energy << " found at alpha = " << best_alpha << std::endl;
    
    return 0;
}
