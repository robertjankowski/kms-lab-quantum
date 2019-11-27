#include "simulation.hpp"
#include <cmath>
#include <iostream>

Simulation::Simulation(uint N, double tau, double deltaTau, double kappa, double omega) : m_N(N), m_tau(tau),
                                                                                          m_deltaTau(deltaTau),
                                                                                          m_kappa(kappa),
                                                                                          m_omega(omega) {
    initialize(N);
}

void Simulation::initialize(uint N) {
    const auto n = N + 1;
    m_Hreal.resize(n);
    m_Himg.resize(n);

    const auto sq2 = std::sqrt(2);
    for (uint i = 0; i < n; ++i) {
        const auto xk = static_cast<double>(i) / N;
        m_xK.push_back(xk);  // (27)
        m_phiReal.push_back(sq2 * std::sin(m_n * M_PI * xk)); // (29)
        m_phiImg.push_back(0.0);  // (29)
    }

    std::fill(m_Hreal.begin(), m_Hreal.end(), 0.0);
    std::fill(m_Himg.begin(), m_Himg.end(), 0.0);
    setHamiltonian(m_Hreal, m_phiReal);
    setHamiltonian(m_Himg, m_phiImg);
}

void Simulation::setHamiltonian(std::vector<double> &hamiltonian, const std::vector<double> &phi) {
    for (uint i = 1; i < m_N; ++i) {
        hamiltonian.at(i) =
                -0.5 * (phi.at(i + 1) + phi.at(i - 1) - 2 * phi.at(i)) * m_N * m_N +
                m_kappa * (m_xK.at(i) - 0.5) * phi.at(i) * std::sin(m_omega * m_tau);
    }
}

double Simulation::getRho(uint k) {
    const auto phiReal = m_phiReal.at(k);
    const auto phiImg = m_phiImg.at(k);
    return phiReal * phiReal + phiImg * phiImg;
}

std::vector<double> Simulation::getRho() {
    std::vector<double> rhos;
    rhos.reserve(m_N);
    for (uint i = 0; i < m_N + 1; ++i)
        rhos.push_back(getRho(i));
    return rhos;
}

void Simulation::setMetrics() {
    double N{0.0};
    double xMean{0.0};
    double epsilon{0.0};
    for (uint i = 0; i < m_N + 1; ++i) {
        const auto phiReal = m_phiReal.at(i);
        const auto phiImg = m_phiImg.at(i);
        const auto Hreal = m_Hreal.at(i);
        const auto Himg = m_Himg.at(i);

        N += phiReal * phiReal + phiImg * phiImg;
        xMean += m_xK.at(i) * (phiReal * phiReal + phiImg * phiImg);
        epsilon += phiReal * Hreal + phiImg * Himg;
    }
    const auto deltaX = 1.0 / m_N;
    m_metrics = Metrics{N * deltaX, xMean * deltaX, epsilon * deltaX};
}


void Simulation::run() {
    const auto nStep = static_cast<uint>(m_tau / m_deltaTau);
    double step{0.0};
    for (uint i = 0; i < nStep; ++i) {
        for (uint k = 0; k < m_N; ++k)
            m_phiReal.at(k) += m_Himg.at(k) * m_deltaTau / 2; // (32)

        setHamiltonian(m_Hreal, m_phiReal);
        for (uint k = 0; k < m_N; ++k)
            m_phiImg.at(k) -= m_Hreal.at(k) * m_deltaTau; // (33)

        setHamiltonian(m_Himg, m_phiImg);
        setHamiltonian(m_Hreal, m_phiReal);

        for (uint k = 0; k < m_N; ++k)
            m_phiReal.at(k) += m_Himg.at(k) * m_deltaTau / 2; // (34)

        step += m_deltaTau;
        const auto m = getMetrics();
        std::cout << "Step: " << step << "\tX_mean: " << m.xMean << "\tEpsilon: " << m.epsilon << '\n';
    }
}