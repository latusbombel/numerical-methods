#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

double obliczCalke(int nx, int ny, double delta, std::vector<std::vector <double>> u){
    double rezultat = 0.;
    for(int i = 0; i<=nx; i++){
        for(int j = 0; j<=ny; j++){
            rezultat += u[i][j] *delta*delta;
        }
    }
    return rezultat;
}

double obliczXsr(int nx, int ny, double delta, std::vector<std::vector <double>> u){
    double rezultat = 0.;
    double x;
    for(int i = 0; i<=nx; i++){
        for(int j = 0; j<=ny; j++){
            x = delta*i;
            rezultat += u[i][j] *delta*delta * x;
        }
    }
    return rezultat;
}

int main(){
    // inicjalizacja danych

    int nx = 400, ny = 90, i1 = 200, i2 = 210, j1 = 50;
    double delta = 0.01, sigma = 10*delta, xA = 0.45, yA = 0.45, D = 0.1;

    std::vector<std::vector <double>> uOld(nx+1, std::vector<double>(ny+1));
    std::vector<std::vector <double>> uNew(nx+1, std::vector<double>(ny+1));
    std::vector<std::vector <double>> vx(nx+1, std::vector<double>(ny+1));
    std::vector<std::vector <double>> vy(nx+1, std::vector<double>(ny+1));
    std::vector<double> calka;
    std::vector<double> xSr;
    
    // wczytanie pliku Psi
    std::ifstream psiPlik("psi.dat");

    // wczytanie plików wynikowych

    std::ofstream calkiPlik("calki1.dat");

    std::vector<std::vector <double>> PsiWektor;
    PsiWektor.reserve(nx + 100);
    std::vector<double> rezerwacja;
    rezerwacja.reserve(ny + 100);
    int temp_i, temp_j;
    double temp_psi;

    while(psiPlik >> temp_i >> temp_j >> temp_psi){
        if (temp_j == 0)PsiWektor.push_back(rezerwacja);
        PsiWektor[temp_i].push_back(temp_psi);
    }

    // pole prędkości

    for(int i = 1; i< nx; i++){
        for(int j = 1; j < ny; j++){
            vx[i][j] = (PsiWektor[i][j+1] - PsiWektor[i][j-1])/(2*delta);
            vy[i][j] = -(PsiWektor[i+1][j] - PsiWektor[i-1][j])/(2*delta);
        }
        vx[i][0] = 0.;
        vy[i][ny] = 0.;
    }

    for(int i = i1; i<=i2; i++){
        for(int j = 0; j<=j1; j++){
            vx[i][j] = 0.;
        }
    }

    for(int j = 0; j<=ny; j++){
        vx[0][j] = vx[1][j];
        vx[nx][j] = vx[nx-1][j];
    }

    // vmax i krok czasowy dt

    double vMax, dt;
    vMax = sqrt(vx[0][0]*vx[0][0] + vy[0][0]*vy[0][0]);

    for(int i = 1; i<= nx; i++){
        for(int j = 1; j<= ny; j++){
            double vTemp = sqrt(vx[i][j]*vx[i][j] + vy[i][j]*vy[i][j]);
            if(vTemp > vMax) vMax = vTemp;
        }
    }

    dt = delta/(4*vMax);

    // inicjalizacja gestości u

    for(int i = 0; i<=nx; i++){
        for(int j = 0; j<=ny; j++){
            uOld[i][j] = 1/(2. * M_PI * sigma*sigma) * exp(-1* ((delta*i - xA)*(delta*i - xA) + (delta*j - yA)*(delta*j - yA))/(2*sigma*sigma));
        }
    }

    // obliczenia
    int it = 0, itMax = 10000;
    calka.reserve(itMax);
    xSr.reserve(itMax);
    // std::ofstream uPlik("gestosc.dat");
    // int ktoryPlik = 0;

    do{
        //iteracja Picarda
        for(int i = 0; i<=nx; i++){
                for(int j = 0; j<=ny; j++){
                    uNew[i][j] = uOld[i][j];
                }
            }
        int k = 0;
        do{  
            for(int i = 0; i<=nx; i++){
                for(int j = 0; j<=ny; j++){
                    if((i > i1 && i < i2) && (j < j1)) continue;

                    else if(i == 0){
                        uNew[i][j] = (1/(1+2*D*dt/(delta*delta)))*(uOld[i][j] - dt/2. * vx[i][j] * ((uOld[i+1][j] - uOld[nx][j])/(2*delta) + (uNew[i+1][j] - uNew[nx][j])/(2*delta)) - dt/2. * vy[i][j] * ((uOld[i][j+1] - uOld[i][j-1])/(2*delta) + (uNew[i][j+1] - uNew[i][j-1])/(2*delta)) + dt/2. * D * ((uOld[i+1][j] + uOld[nx][j] + uOld[i][j+1] + uOld[i][j-1] - 4* uOld[i][j])/(delta*delta) +  (uNew[i+1][j] + uNew[nx][j] + uNew[i][j+1] + uNew[i][j-1])/(delta*delta)));
                    }
                    else if(i == nx){
                        uNew[i][j] = (1/(1+2*D*dt/(delta*delta)))*(uOld[i][j] - dt/2. * vx[i][j] * ((uOld[ 0][j] - uOld[i-1][j])/(2*delta) + (uNew[ 0][j] - uNew[i-1][j])/(2*delta)) - dt/2. * vy[i][j] * ((uOld[i][j+1] - uOld[i][j-1])/(2*delta) + (uNew[i][j+1] - uNew[i][j-1])/(2*delta)) + dt/2. * D * ((uOld[0][j] + uOld[i-1][j] + uOld[i][j+1] + uOld[i][j-1] - 4* uOld[i][j])/(delta*delta) +  (uNew[0][j] + uNew[i-1][j] + uNew[i][j+1] + uNew[i][j-1])/(delta*delta)));
                    }
                    else{
                        uNew[i][j] = (1/(1+2*D*dt/(delta*delta)))*(uOld[i][j] - dt/2. * vx[i][j] * ((uOld[i+1][j] - uOld[i-1][j])/(2*delta) + (uNew[i+1][j] - uNew[i-1][j])/(2*delta)) - dt/2. * vy[i][j] * ((uOld[i][j+1] - uOld[i][j-1])/(2*delta) + (uNew[i][j+1] - uNew[i][j-1])/(2*delta)) + dt/2. * D * ((uOld[i+1][j] + uOld[i-1][j] + uOld[i][j+1] + uOld[i][j-1] - 4* uOld[i][j])/(delta*delta) +  (uNew[i+1][j] + uNew[i-1][j] + uNew[i][j+1] + uNew[i][j-1])/(delta*delta)));
                    }
                }
            }
            
            k++;
        }while(k<20);
        for(int i = 0; i<=nx; i++){
            for(int j = 0; j<=ny; j++){
                uOld[i][j] = uNew[i][j];
            }
        }
        calka.push_back(obliczCalke(nx, ny, delta, uNew));
        xSr.push_back(obliczXsr(nx, ny, delta, uNew));

        double t = dt*it;
        calkiPlik << std::fixed << std::setprecision(5) << std::setw(6) << 
        t << " " << xSr[it] << " " << calka[it] << std::endl;


        // if(it % (itMax/50) == 0){
        //     std::string nazwa = "gestosc" + std::to_string(ktoryPlik) + ".dat";

        //     std::ofstream uPlik(nazwa.c_str());
        //     for(int i = 0; i<=nx; i++){
        //         for(int j = 0; j<=ny; j++){
        //             uPlik << i*delta << " " << j*delta << " " << 
        //             std::fixed << std::setprecision(5) << std::setw(6) <<
        //             uNew[i][j] << std::endl;
        //         }
        //         uPlik << std::endl;
        //     }
        //     ktoryPlik ++;
        //     uPlik.close();
        // }
        it++;
    }while(it < itMax);

    // std::ofstream vPlik("predkosc1.dat");

    // for(int i = 0; i<=nx; i++){
    //     for(int j = 0; j<=ny; j++){
    //         vPlik << i*delta << " " << j*delta << " " << 
    //         std::fixed << std::setprecision(5) << std::setw(6) << 
    //         vx[i][j] << " " << vy[i][j] << std::endl;
    //     }
    // }

    // int m = 0;
    // do{
    //     std::cout << m <<  " martynka" << std::endl;
    //     m++;
    // }while(m<100);

    return 0;
}