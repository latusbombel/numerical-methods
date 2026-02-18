#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

double funkcja(double alpha, double beta, double y){
    return alpha*y - beta*y*y;
}

double met_picarda(std::vector<double> y, int max_iter, double t, int i, double y1, double y2, double TOL){
    int iter = 0;
    do{
        double y_old = y[i+1];
        y[i+1] = y[i] + t/2*(y1 + y2);
        if(fabs(y_old - y[i+1]) < TOL) break;
        iter++;
    }while(iter <= max_iter);
    return y[i+1];
}

int main(){
    int tmax = 100, N = 500, max_iter = 20;
    double y0 = 1, beta = 0.001, gamma = 0.1, t = 0.1, TOL = pow(10, -6);
    double alpha = beta * N - gamma;

    //inicjalizacja wektora wyników

    std::vector<double> y{y0};
    std::vector<double> z{N-y0};

    //obliczenia

    for(int i = 0; ; i++){
        y.push_back(y[i]);
        y[i+1] = met_picarda(y, max_iter, t, i, funkcja(alpha, beta, y[i]), funkcja(alpha, beta, y[i+1]), TOL);
        z.push_back(N-y[i+1]);
        if((i+1)*t>=tmax) break;
    }

    //zapis wyników

    std::ofstream plik ("wyniki1.dat");
    for(int i = 0; i < y.size(); i++){
        plik << std::fixed << std::setprecision(4) 
            << std::setw(8) << i*t << ' ' 
            << std::setw(8) << y[i] << ' '
            << std::setw(8) << z[i] << std::endl;
    }

    plik.close();
    return 0;
}