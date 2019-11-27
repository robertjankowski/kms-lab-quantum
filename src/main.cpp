#include "simulation.hpp"

int main() {
    const std::string metricsFile = "../output/metrics.csv";
    const std::string densityFile = "../output/density.csv";

    Simulation simulation("../input/parameters.txt");
    simulation.run(metricsFile, densityFile);
    return 0;
}