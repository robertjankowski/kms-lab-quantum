#include "simulation.hpp"
#include <sstream>

std::string getCustomFilename(std::string_view basePath, std::string_view filename, int kappa, int omega, int n) {
    std::ostringstream oss;
    oss << basePath << "/" << filename << "_kappa=" << kappa << "_omega=" << omega << "_n=" << n << ".csv";
    return oss.str();
}

int main() {
    std::string_view basePath = "../output";
    const auto kappa = 9;
    for (int i = 5; i < 50; ++i) {
        Simulation simulation("../input/parameters.txt");
        const std::string metricsFile = getCustomFilename(basePath, "new_resonance_metrics", kappa, i, 9);
        const std::string densityFile = getCustomFilename(basePath, "new_resonance_density", kappa, i, 9);
        simulation.setOmega(i);
        simulation.run(metricsFile, densityFile);
    }
    return 0;
}