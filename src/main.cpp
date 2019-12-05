#include "simulation.hpp"
#include <sstream>

std::string getCustomFilename(std::string_view basePath, std::string_view filename, int omega, int n) {
    std::ostringstream oss;
    oss << basePath << "/" << filename << "_omega=" << omega << "_n=" << n << ".csv";
    return oss.str();
}

int main() {
    Simulation simulation("../input/parameters.txt");

    std::string_view basePath = "../output";
    for (int i = 91; i < 96; ++i) {
        const std::string metricsFile = getCustomFilename(basePath, "resonance_metrics_kappa=3_", i, 9);
        const std::string densityFile = getCustomFilename(basePath, "resonance_density_kappa=3_", i, 9);
        simulation.setOmega(i);
        simulation.run(metricsFile, densityFile);
    }
    return 0;
}