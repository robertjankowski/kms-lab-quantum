#ifndef KMS_LAB_QUANTUM_SIMULATION_HPP
#define KMS_LAB_QUANTUM_SIMULATION_HPP

#include <vector>

using uint = unsigned int;

class Simulation {
    std::vector<double> m_xK;
    std::vector<double> m_phiReal;
    std::vector<double> m_phiImg;
    std::vector<double> m_Hreal;
    std::vector<double> m_Himg;
    uint m_N;
    double m_tau{0.001};
    double m_kappa;
    double m_omega;
public:
    Simulation(uint N, double tau, double kappa, double omega);

private:
    void initialize(uint N);

};


#endif //KMS_LAB_QUANTUM_SIMULATION_HPP
