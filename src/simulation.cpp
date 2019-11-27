#include "simulation.hpp"
#include <cmath>

Simulation::Simulation(uint N, double tau, double kappa, double omega) : m_N(N), m_tau(tau), m_kappa(kappa),
                                                                         m_omega(omega) {
    initialize(N);
}

void Simulation::initialize(uint N) {
    const auto n = N + 1;
    m_xK.reserve(n);
    m_phiImg.reserve(n);
    m_phiReal.reserve(n);
    m_Hreal.reserve(n);
    m_Himg.reserve(n);

    const auto sq2 = std::sqrt(2);
    for (uint i = 0; i <= N; ++i) {
        const auto xk = static_cast<double>(i) / N;
        m_xK.push_back(xk);  // (27)
        m_phiReal.push_back(sq2 * std::sin(i * M_PI * xk))  // (29)
        m_phiImg.push_back(0.0);  // (29)
    }

    m_Hreal.at(0) = 0;
    m_Himg.at(0) = 0;

    for (uint i = 1; i < N; ++i) {
        // update hamiltonian
    }

}