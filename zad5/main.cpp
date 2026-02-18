#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

double S(const double V[][129], int k){
    int nx = 128, ny = 128;
    double S = 0, delta = 0.2;

    for(int i = 0; i<nx-k; i=i+k){
        for(int j = 0; j<ny-k; j=j+k){
            S += pow(k*delta, 2)/2.0*(pow((V[i+k][j] - V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i][j+k])/(2*k*delta),2)+1/2.0*pow((V[i][j+k] - V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i][j+k])/(2*k*delta),2));
        }
    }
    return S;
}

int main(){
    const int nx = 128, ny = 128;
    double  x, y;
    double delta = 0.2;
    double x_max = delta*nx;
    double y_max = delta*ny;
    int k_tab[5] = {16,8,4,2,1};
    double TOL = 10e-8;
    double S_new, S_old;

    double V[nx+1][ny+1];

    std::ofstream plik ("akcja.dat");
    std::ofstream plik2 ("potencjal.dat");

    for(int i = 0; i<nx+1; i++){
        x = delta * i;
        V[i][nx] = (-1.0)*sin(2*M_PI*x/x_max);
        V[i][0] = 1.0*sin(2*M_PI*x/x_max);
    }

    for(int j=0; j<ny+1;j++){
        y = delta*j;
        V[0][j] = 1.0*sin(M_PI*y/y_max);
        V[ny][j] = 1.0*sin(M_PI*y/y_max);
    }

    for(int i = 1; i<nx; i++){
        for(int j=1; j<ny;j++){
            V[i][j] = 0.0;
        }
    }

    for(int k : k_tab){
        int iter = 0;
        do{
            S_old = S(V,k);
            //obliczenie w punktach
            for(int i = k; i<=nx-k;i = i+k){
                for(int j = k; j<=ny-k; j = j+k){
                    V[i][j] = 1/4.0*(V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);
                }
            }
            //oblizcenie sredniej w punktach nieoblicznaych
            if(k != 1){
                for(int i = k; i<nx-k; ){
                    for(int j = k; j<ny-k; ){
                        V[i+k/2][j+k/2] = 1/2.0*(V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]);
                        V[i+k][j+k/2] = 1/2.0*(V[i+k][j]+ V[i+k][j+k]);
                        V[i+k/2][j+k] = 1/2.0*(V[i][j+k] + V[i+k][j+k]);
                        V[i+k/2][j] = 1/2.0*(V[i][j] + V[i+k][j]);
                        V[i][j+k/2] = 1/2.0*(V[i][j] + V[i][j+k]);
                        j = j+k;
                    }
                    i = i+k;
                }
            }
            
            S_new = S(V,k);
            plik << std::fixed << std::setprecision(4) 
            << std::setw(2) << iter << ' ' << std::setprecision(8) 
            << std::setw(8) << S_new
            << std::endl;
            iter++;
        }while(fabs((S_new-S_old)/S_old) > TOL);
        plik << std::endl << std::endl;

        //wypisanie danych
        x = 0;
        for(int i = 0; i<=nx-k; i=i+k){
            x = i * delta;
            y = 0;
            for(int j=0; j<=ny-k;j=j+k){
                y = j * delta;
                plik2 << std::fixed << std::setprecision(10) 
                << std::setw(4) << x << ' ' << ' ' << ' ' << std::setw(4) << y << ' ' << ' ' << ' '
                << std::setw(8) << V[i][j] << std:: endl;
            }
            plik2 << std::endl;
        }
        plik2 << std::endl << std::endl;
    }
        
    return 0;
}