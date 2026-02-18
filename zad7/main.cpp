#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

int main(){
    int nx = 200, ny = 90, i1 = 50, j1 = 55, j2 = 57, itMax = 20000;
    double Qwyjsciowe, Omega, delta = 0.01, rho = 1, mi = 1, Qwejsciowe = -1000;

    double** psi = new double* [ny+1];
    double** zeta = new double* [ny+1];
    for(int i = 0; i<ny+1; i++){
        psi[i] = new double[nx+1];
        zeta[i] = new double[nx+1];
        for(int j = 0 ; j<nx+1; j++){
            psi[i][j] = 0.0;
            zeta[i][j] = 0.0;
        }
    }

    std::ofstream plik ("psi.dat");
    std::ofstream plik2 ("zeta.dat");


    double yj1 = delta * j1;
    double yny = delta * ny;
    Qwyjsciowe = Qwejsciowe * (pow(yny, 3) - pow(delta * j1, 3) - 3*yj1*pow(yny, 2) + 3*pow(delta * j1, 2)* yny)/pow(yny, 3);

    // wawrunki brzegowe
    // dla psi

    for(int i = j1; i<ny+1; i++){
        psi[i][0] = Qwejsciowe/(2. * mi)*(pow(delta * i, 3)/3. - pow(delta * i, 2)/2.*(yj1 + yny) + delta*i * yj1 * yny);
    }

    for(int i = 0; i<=ny; i++){
        psi[i][nx] = Qwyjsciowe/(2. * mi)*(pow(delta * i, 3)/3. - pow(delta * i, 2)/2. * yny) + Qwejsciowe * pow(yj1, 2)*((-yj1)+ 3* yny)/(12. * mi);
    }

    for(int j = 1; j<nx; j++){
        psi[ny][j] = psi[ny][0];
    }

    for(int j = i1; j<nx; j++){
        psi[0][j] = psi[j1][0];
    }


    for(int i = 1; i<= j1; i++){
        psi[i][i1] = psi[j1][0];
    }

    for(int j = 1; j <= i1; j++){
        psi[j1][j] = psi[j1][0];
    }

    //dla zeta

    for(int i = j1; i<=ny; i++){
        zeta[i][0] = Qwejsciowe/(2. * mi)* (2*delta*i - yj1 - yny);
    }

    for(int i = 0; i<=ny; i++){
        zeta[i][nx] = Qwyjsciowe/(2. * mi) * (2*delta*i - yny);
    }

    for(int j = 1; j < nx; j++){
        zeta[ny][j] = 2./(delta*delta) * (psi[ny-1][j] - psi[ny][j]);
    }

    for(int j = i1 + 1; j<nx; j++){
        zeta[0][j] = 2./(delta*delta) * (psi[1][j] - psi[0][j]);
    }

    for(int i = 1; i<j1; i++){
        zeta[i][i1] = 2./(delta*delta) * (psi[i][i1+1] - psi[i][i1]);
    }
    
    for(int j = 1; j<=i1; j++){
        zeta[j1][j] = 2./(delta*delta) * (psi[j1+1][j] - psi[j1][j]);
    }

    zeta[j1][i1] = 0.5*(zeta[j1-1][i1] + zeta[j1][i1-1]);

    //wyliczanie
    int it = 0;

    do{
        if(it < 200){ Omega = 0.0;}
        else {Omega = 1.0;}

        for(int i = 1; i<ny; i++){
            for(int j = 1; j<nx; j++){
                if(i <= j1 && j <= i1) continue;
                else{
                    psi[i][j] = 1/4. * (psi[i][j+1] + psi[i][j-1] + psi[i+1][j] + psi[i-1][j] - delta*delta*zeta[i][j]);
                    zeta[i][j] = 1/4. * (zeta[i][j+1] + zeta[i][j-1] + zeta[i+1][j] + zeta[i-1][j]) - Omega*rho/(16. * mi)*((psi[i+1][j] - psi[i-1][j])*(zeta[i][j+1] - zeta[i][j-1]) - (psi[i][j+1] - psi[i][j-1])*(zeta[i+1][j]-zeta[i-1][j]));
                }
            }
        }

        for(int i = j1; i<=ny; i++){
            zeta[i][0] = Qwejsciowe/(2. * mi)* (2*delta*i - yj1 - yny);
        }

        for(int i = 0; i<=ny; i++){
            zeta[i][nx] = Qwyjsciowe/(2. * mi) * (2*delta*i - yny);
        }

        for(int j = 1; j < nx; j++){
            zeta[ny][j] = 2./(delta*delta) * (psi[ny-1][j] - psi[ny][j]);
        }

        for(int j = i1 + 1; j<nx; j++){
            zeta[0][j] = 2./(delta*delta) * (psi[1][j] - psi[0][j]);
        }

        for(int i = 1; i<j1; i++){
            zeta[i][i1] = 2./(delta*delta) * (psi[i][i1+1] - psi[i][i1]);
        }
        
        for(int j = 1; j<=i1; j++){
            zeta[j1][j] = 2./(delta*delta) * (psi[j1+1][j] - psi[j1][j]);
        }

        zeta[j1][i1] = 0.5*(zeta[j1-1][i1] + zeta[j1][i1-1]);

        it++;
    }while(it <= itMax);

    //wypisanie danych

    
    // std::ofstream plik3 ("potencjal.dat");

    for(int i = 0; i<ny+1; i++){
        for(int j = 0 ; j<nx+1; j++){
            plik << std::fixed << std::setprecision(5) 
            << std::setw(6) << psi[i][j] << ' ';
            plik2 << std::fixed << std::setprecision(5) 
            << std::setw(6) << zeta[i][j] << ' ';
        }
        plik << std::endl;
        plik2 << std::endl;
    }

    //zwalnianie pamieci

    for(int i = 0; i<ny; i++){
        delete[] psi[i];
        delete[] zeta[i];
    }
    delete[] psi;
    delete[] zeta;

    return 0;
}