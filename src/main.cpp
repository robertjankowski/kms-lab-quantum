#include "simulation.hpp"
#include <sstream>

std::string getCustomFilename(std::string_view basePath, std::string_view filename, int kappa, double omega, int n) {
    std::ostringstream oss;
    oss << basePath << "/" << filename << "_kappa=" << kappa << "_omega=" << omega << "_n=" << n << ".csv";
    return oss.str();
}

int main() {
    std::string_view basePath = "../output";
    const auto kappa = 1;
    auto step{14.8};
    for (int i = 0; i < 1; ++i) {
        Simulation simulation("../input/parameters.txt");
        const std::string metricsFile = getCustomFilename(basePath, "metrics_resonance", kappa, step, 1);
        const std::string densityFile = getCustomFilename(basePath, "density_resonance", kappa, step, 1);
        simulation.setOmega(step);
        simulation.run(metricsFile, densityFile);
        step += 0.1;
    }
    return 0;
}