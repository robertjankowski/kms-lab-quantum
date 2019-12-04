#include "simulation.hpp"

int main() {
    Simulation simulation("../input/parameters.txt");

    const std::string metricsFile = "../output/metrics_n=1_kappa=8_omega=32.csv";
    const std::string densityFile = "../output/density__n=1_kappa=8_omega=32.csv";
    simulation.run(metricsFile, densityFile);
    return 0;
}