#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

double przyspieszenie(int i, double x, double xF, double t, double tMax,double delta, double alfa, double beta, double dt, std::vector<double> uNew, std::vector<double> uOld){
    double deltaKron;
    if(fabs(x - xF) < 1e-9) deltaKron = 1.;
    else deltaKron = 0.;
    double aF = cos(50./tMax * t) * deltaKron;
    return (uNew[i+1] - 2*uNew[i] + uNew[i-1])/(delta*delta) - beta*(uNew[i] - uOld[i])/dt + alfa * aF;
}

int main(){
    // inicjalizacja danych
    int nx = 150, nt = 1000;
    double delta = 0.1, dt = 0.05, t = 0., sigma = 0.5, xA = 7.5, alfa = 1., beta = 1., xF = 2.5;
    double E, x;

    std::vector<double> uOld;
    uOld.resize(nx+1);
    std::vector<double> vOld;
    vOld.resize(nx+1);
    std::vector<double> uNew;
    uNew.resize(nx+1);
    std::vector<double> a;
    a.resize(nx+1);
    std::vector<double> vNew;
    vNew.resize(nx+1);

    std::ofstream Eplik("energia3.dat");
    std::ofstream uPlik("wychylenie3.dat");

    // warunki brzegowe
    uNew[0] = 0.0; uNew[nx] = 0.0;
    uOld[0] = 0.0; uOld[nx] = 0.0;
    vNew[0] = 0.0; vNew[nx] = 0.0;
    vOld[0] = 0.0; vOld[nx] = 0.0;

    // warunki początkowe
    for(int i = 1; i<nx; i++){
        double x = i*delta;
        uNew[i] = 0.0; //exp(-pow((x-xA)/sigma,2)*0.5);
        uOld[i] = uNew[i];
        vOld[i] = 0.0;
        vNew[i] = 0.0;
    }

    // obliczenia
    int tMax = dt * nt;
    for(int i=1; i<nx; i++){
        double x = i*delta;
        a[i] = przyspieszenie(i, x, xF, t, tMax, delta, alfa, beta, dt, uNew, uOld);
    }

    for(int n = 1; n< nt; n++){
        t = dt*n;
        for(int i = 1; i<nx; i++){
            vOld[i] = vNew[i] + dt/2.*a[i];
            uOld[i] = uNew[i];
            uNew[i] = uNew[i] + dt*vOld[i];
        }
        for(int i = 1; i<nx; i++){
            x = i*delta;
            a[i] = przyspieszenie(i, x, xF, t, tMax, delta, alfa, beta, dt, uNew, uOld);
            vNew[i] = vOld[i] + dt/2.*a[i];
        }

        E = delta/4.*(pow((uNew[1] - uNew[0])/delta, 2) + pow((uNew[nx] - uNew[nx-1])/delta, 2));
        for(int i = 1; i<nx; i++){
            E += delta/2.*(vNew[i]*vNew[i] + pow((uNew[i+1] - uNew[i-1])/(2*delta),2));
        }
        Eplik << std::fixed << std::setprecision(2) << std::setw(3) << t << " " << 
        std::fixed << std::setprecision(5) << std::setw(6) << E << std::endl;

        for(int i = 0; i<=nx; i++){
            x = i*delta;
            uPlik << std::fixed << std::setprecision(2) << std::setw(3) << t << " " << x << " " << 
            std::fixed << std::setprecision(5) << std::setw(6) << 
            uNew[i] << std::endl;
        }
        // uPlik << std::endl;
    }

    return 0;
}