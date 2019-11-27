#ifndef KMS_LAB_QUANTUM_SIMULATION_HPP
#define KMS_LAB_QUANTUM_SIMULATION_HPP

#include <vector>

using uint = unsigned int;

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
    Metrics m_metrics;
public:
    Simulation(uint N, double tau, double deltaTau, double kappa, double omega);

    std::vector<double> getRho();

    [[nodiscard]] Metrics getMetrics() {
        setMetrics();
        return m_metrics;
    };

    void run();

private:
    void initialize(uint N);

    void setHamiltonian(std::vector<double> &hamiltonian, const std::vector<double> &phi);

    double getRho(uint k);

    void setMetrics();
};


#endif //KMS_LAB_QUANTUM_SIMULATION_HPP
