#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#define M_PI 3.14159265358979323846

int main(){
    int tmin = 0, R = 100, i, j;
    double t = 0.0001, L = 0.1, C = 0.001;
    double I0 = 0, Q0 = 0, lambda = -1;
    double omega_0 = 1/sqrt(L*C);
    double T_0 = 2*M_PI/omega_0;
    double tmax = 4*T_0;
    double V, V_1, V_2;
    double omega_v[4] = {0.5 * omega_0, 0.8 * omega_0, 1.0 * omega_0, 1.2 * omega_0};


    std::vector<std::vector<double>> I{{I0}, {I0}, {I0}, {I0}};
    std::vector<std::vector<double>> Q{{Q0}, {Q0}, {Q0}, {Q0}};
    std::ofstream plik("wyniki4_1.dat");

    //obliczenia

    for(j = 0; j<4; j++){
        int krok = 1;
        do{
            V = 10.0 * sin(omega_v[j] * t * krok);
            V_1 = 10.0 * sin(omega_v[j] * t * (krok+1/2.0));
            V_2 = 10.0 * sin(omega_v[j] * t * (krok+1));
            double k1_Q = I[j][krok-1];
            double k1_I = V/L - 1/(L*C)*Q[j][krok-1] - R/L * I[j][krok-1];
            double k2_Q = I[j][krok-1] + t/2 * k1_I;
            double k2_I = V_1/L - 1/(L*C)*(Q[j][krok-1] + t/2 * k1_Q) - R/L * (I[j][krok-1] + t/2 * k1_I);
            double k3_Q = I[j][krok-1] + t/2 *k2_I;
            double k3_I = V_1/L - 1/(L*C)*(Q[j][krok-1] + t/2 * k2_Q) - R/L * (I[j][krok-1] + t/2 * k2_I);
            double k4_Q = I[j][krok-1] + t * k3_I;
            double k4_I = V_2/L - 1/(L*C)*(Q[j][krok-1] + t * k3_Q) - R/L * (I[j][krok-1] + t * k3_I);
            Q[j].push_back(Q[j][krok-1] + t/6.0*(k1_Q + 2.0*k2_Q + 2.0*k3_Q + k4_Q));
            I[j].push_back(I[j][krok-1] + t/6.0*(k1_I + 2.0*k2_I + 2.0*k3_I + k4_I));
            krok++;
        }while(krok*t<=tmax);
    }

    //wypisanie I i Q

    for(j = 0; j<4; j++){
        for(i=0; i<I[j].size(); i++){
            plik << std::fixed << std::setprecision(10)
                << std::setw(12) << i*t << ' '
                << std::setw(12) <<  I[j][i] << ' '
                << std::setw(12) << Q[j][i] << std::endl;
        }
        plik << std::endl << std::endl;
    }

    plik.close();

    return 0;
}