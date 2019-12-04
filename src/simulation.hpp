#ifndef KMS_LAB_QUANTUM_SIMULATION_HPP
#define KMS_LAB_QUANTUM_SIMULATION_HPP

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

using uint = unsigned int;

template<typename T>
T split(const std::string &s, char delim);

struct Metrics {
    double N, xMean, epsilon;
};

class Simulation {
    std::vector<double> m_xK;
    std::vector<double> m_phiReal;
    std::vector<double> m_phiImg;
    std::vector<double> m_Hreal;
    std::vector<double> m_Himg;
    uint m_N;
    double m_tau;
    double m_deltaTau;
    double m_kappa;
    double m_omega;
    double m_n{1};
    uint m_nStepMetrics;
    uint m_nStepDensity;
    Metrics m_metrics{};
public:
    Simulation(uint N, double tau, double deltaTau, double kappa, double omega, uint nStepMetrics, uint nStepDensity);

    explicit Simulation(const std::string &filename);

    [[nodiscard]] std::vector<double> getRho() const noexcept;

    [[nodiscard]] Metrics getMetrics() {
        setMetrics();
        return m_metrics;
    };

    void run(const std::string &metricsFilename, const std::string &densityFilename);

    void setOmega(double omega) { m_omega = omega; }

private:
    void initialize(uint N);

    void setHamiltonian(std::vector<double> &hamiltonian, const std::vector<double> &phi, double tau) noexcept;

    [[nodiscard]] double getRho(uint k) const;

    void setMetrics();

    template<typename T>
    T getNextParameter(std::ifstream &file, std::string &line, char delim = '#') {
        std::getline(file, line);
        return split<T>(line, delim);
    }

    template<typename T>
    void checkIfFileExists(T &file) {
        static_assert(std::is_base_of<std::ios, T>::value, "T must inherit from ios");
        if (file.fail()) {
            std::cerr << "Unable to open file\n";
            exit(0);
        }
    }

    void saveMetrics(double step, std::ofstream &file);

    void saveDensity(std::ofstream &file);
};

template<typename T>
T split(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> elems;
    while (std::getline(ss, item, delim)) {
        elems.push_back(std::move(item));
    }
    std::stringstream val(elems.at(0));
    T value;
    val >> value;
    return value;
}

#endif //KMS_LAB_QUANTUM_SIMULATION_HPP
