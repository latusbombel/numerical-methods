#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

double funkcja(double alpha, double beta, double y){
    return alpha*y - beta*y*y;
}

void met_newtona(double U[], double TOL, int max_iter, int i, double t, std::vector<double> y, double alpha, double beta){
    int iter = 0;
    double F1, F2;
    double U1 = U[0];
    double U2 = U[1];
    double delta_U1, delta_U2;

    

    do{
        double U_old1 = U1;
        double U_old2 = U2;
        double Alpha[2][2] = {
            {1/4.0, 1/4.0 - sqrt(3)/6},
            {1/4.0 + sqrt(3)/6, 1/4.0}
                };

        double m[2][2] = {
            {1 - t*Alpha[0][0] * (alpha - 2 * beta * U1), -t * Alpha[0][1]*(alpha - 2*beta*U2)},
            {-t * Alpha[1][0]*(alpha - 2*beta*U1), 1 - t*Alpha[1][1] * (alpha - 2 * beta * U2)}
                };
        F1 = U1 - y[i] - t*(Alpha[0][0] * (alpha * U1 - beta * U1*U1) + Alpha[0][1] * (alpha * U2 - beta * U2 *U2));
        F2 = U2 - y[i] - t*(Alpha[1][0] * (alpha * U1 - beta * U1*U1) + Alpha[1][1] * (alpha * U2 - beta * U2 *U2));
        delta_U1 = (F2 * m[0][1] - F1 * m[1][1]) / (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
        delta_U2 = (F1 * m[1][0] - F2 * m[0][0]) / (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
       
        U1 = U1 + delta_U1;
        U2 = U2 + delta_U2;
        iter++;
        U[0] = U1;
        U[1] = U2;
        if(fabs(U_old1 - U1) < TOL && fabs(U_old2 - U2) < TOL) break;
    }while(iter < max_iter);
}


int main(){
    int tmax = 100, N = 500, max_iter = 20;
    double y0 = 1, beta = 0.001, gamma = 0.1, t = 0.1, TOL = pow(10, -6);
    double alpha = beta * N - gamma;

    //metoda RK2

    double U[2];

    double delta_U1, delta_U2;

    //inicjalizacja wektora wyników

    std::vector<double> y{y0};
    std::vector<double> z{N-y0};

    //obliczenia

    for(int i = 0; ; i++){
        y.push_back(y[i]);
        U[0] = y[i];
        U[1] = y[i];
        met_newtona(U, TOL, max_iter, i, t, y, alpha, beta);
        y[i+1] = y[i] + t*(1/2.0 * funkcja(alpha, beta, U[0]) + 1/2.0 * funkcja(alpha, beta, U[1]));
        z.push_back(N-y[i+1]);
        if((i+1)*t>=tmax) break;
    }

    //zapis wyników

    std::ofstream plik ("wyniki3.dat");
    for(int i = 0; i < y.size(); i++){
        plik << std::fixed << std::setprecision(4) 
            << std::setw(8) << i*t << ' ' 
            << std::setw(8) << y[i] << ' '
            << std::setw(8) << z[i] << std::endl;
    }

    plik.close();
    return 0;
}