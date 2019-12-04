#include "simulation.hpp"
#include <cmath>
#include <iostream>

Simulation::Simulation(uint N, double tau, double deltaTau, double kappa, double omega, uint nStepMetrics,
                       uint nStepDensity) : m_N(N), m_tau(tau),
                                            m_deltaTau(deltaTau),
                                            m_kappa(kappa),
                                            m_omega(omega),
                                            m_nStepMetrics(nStepMetrics),
                                            m_nStepDensity(nStepDensity) {
    initialize(N);
}

Simulation::Simulation(const std::string &filename) {
    std::ifstream file(filename);
    checkIfFileExists(file);

    std::string line;
    m_N = getNextParameter<uint>(file, line);
    m_tau = getNextParameter<double>(file, line);
    m_deltaTau = getNextParameter<double>(file, line);
    m_kappa = getNextParameter<double>(file, line);
    m_n = getNextParameter<double>(file, line);
    m_nStepMetrics = getNextParameter<uint>(file, line);
    m_nStepDensity = getNextParameter<uint>(file, line);
    // Simulation with electric field (page 6)
    // omega = {3pi^2/2, 4pi*2/2, 8pi^2/2, ...}
    m_omega = 32 * M_PI * M_PI / 2;
    initialize(m_N);
    file.close();
}


void Simulation::setHamiltonian(std::vector<double> &hamiltonian, const std::vector<double> &phi, double tau) noexcept {
    for (uint i = 1; i < m_N; ++i) {
        hamiltonian.at(i) =
                -0.5 * (phi.at(i + 1) + phi.at(i - 1) - 2 * phi.at(i)) * m_N * m_N +
                m_kappa * (m_xK.at(i) - 0.5) * phi.at(i) * std::sin(m_omega * tau);
    }
}

void Simulation::initialize(uint N) {
    m_metrics = Metrics{0, 0, 0};

    const auto n = N + 1;
    m_xK.reserve(n);
    m_phiReal.reserve(n);
    m_phiImg.reserve(n);
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
    setHamiltonian(m_Hreal, m_phiReal, 0);
    setHamiltonian(m_Himg, m_phiImg, 0);
}

double Simulation::getRho(uint k) const {
    const auto phiReal = m_phiReal.at(k);
    const auto phiImg = m_phiImg.at(k);
    return phiReal * phiReal + phiImg * phiImg;
}

std::vector<double> Simulation::getRho() const noexcept {
    std::vector<double> rhos;
    rhos.reserve(m_N);
    for (uint i = 0; i < m_N + 1; ++i)
        rhos.emplace_back(getRho(i));
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


void Simulation::run(const std::string &metricsFilename, const std::string &densityFilename) {
    std::ofstream metricsFile(metricsFilename);
    std::ofstream densityFile(densityFilename);
    checkIfFileExists(metricsFile);
    checkIfFileExists(densityFile);

    auto meanEnergy{0.0};
    const auto nStep = static_cast<uint>(m_tau / m_deltaTau);
    auto step{0.0};
    for (uint i = 0; i < nStep; ++i) {
        for (uint k = 0; k < m_N; ++k)
            m_phiReal.at(k) += m_Himg.at(k) * m_deltaTau / 2; // (32)

        setHamiltonian(m_Hreal, m_phiReal, step);
        for (uint k = 0; k < m_N; ++k)
            m_phiImg.at(k) -= m_Hreal.at(k) * m_deltaTau; // (33)

        setHamiltonian(m_Himg, m_phiImg, step);
        setHamiltonian(m_Hreal, m_phiReal, step);

        for (uint k = 0; k < m_N; ++k)
            m_phiReal.at(k) += m_Himg.at(k) * m_deltaTau / 2; // (34)

        if (i % m_nStepMetrics == 0)
            saveMetrics(step, metricsFile);
        if (i % m_nStepDensity == 0)
            saveDensity(densityFile);
        step += m_deltaTau;
        const auto[N, xMean, epsilon] = getMetrics();
        meanEnergy += epsilon;
    }
    meanEnergy /= nStep;
    std::cout << m_omega << "," << meanEnergy << std::endl;
    metricsFile.close();
    densityFile.close();
}

void Simulation::saveMetrics(double step, std::ofstream &file) {
    const auto[N, xMean, epsilon] = getMetrics();
    file << step << ',' << N << ',' << xMean << ',' << epsilon << '\n';
}

void Simulation::saveDensity(std::ofstream &file) {
    for (const auto &rho : getRho())
        file << rho << ',';
    file << '\n';
}
