#include <iostream>
#include <algorithm>
#include "simulation.hpp"

int main() {
    constexpr auto N = 100;
    constexpr auto tau = 10;
    constexpr auto deltaTau = 1e-4;
    constexpr auto kappa = 0;
    constexpr auto omega = 0;
    Simulation simulation(N, tau, deltaTau, kappa, omega);
    simulation.run();

    const auto rhos = simulation.getRho();
    std::for_each(rhos.begin(), rhos.end(), [](const auto &rho) { std::cout << rho << '\n'; });

    const auto m = simulation.getMetrics();
    std::cout << m.N << ' ' << m.xMean << ' ' << m.epsilon << '\n';
    return 0;
}