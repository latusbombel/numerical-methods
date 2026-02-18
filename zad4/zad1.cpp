#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

double rho_1(double x, double y){
    double x_max = 150*0.1;
    double y_max = 100*0.1;
    double sigma_x = 0.1*x_max;
    double sigma_y = 0.1*y_max;
    return std::exp(-std::pow((x-0.35*x_max)/sigma_x,2)-std::pow((y-0.5*y_max)/sigma_y,2));
}

double rho_2(double x, double y){
    double x_max = 150*0.1;
    double y_max = 100*0.1;
    double sigma_y = 0.1*y_max;
    double sigma_x = 0.1*x_max;
    return (-1)*std::exp(-std::pow((x-0.65*x_max)/sigma_x,2)-std::pow((y-0.5*y_max)/sigma_y,2));
}

double rho_obl(double rho_1, double rho_2){
    return rho_1+rho_2;
}

double S(const double V[][101], const double rho[][101]){
    int nx = 150, ny = 100;
    double S = 0, delta = 0.1;

    for(int i = 0; i<nx; i++){
        for(int j = 0; j<ny; j++){
            S += pow(delta, 2)*(1/2.0*pow((V[i+1][j] - V[i][j])/delta,2)+1/2.0*pow((V[i][j+1] - V[i][j])/delta,2) - rho[i][j] * V[i][j]);
        }
    }
    return S;
}

int main(){
    const int eps = 1, nx = 150, ny = 100;
    double  V1 = 10;
    double V2 = 0;
    double delta = 0.1;
    double x_max = delta*nx;
    double y_max = delta*ny;
    double omega_tab[2] = {0.6, 1.0};
    double TOL = 10e-8;

    double Vn[nx+1][ny+1];
    double Vs[nx+1][ny+1];
    

    double S_old;
    double S_new;
    double rho[nx+1][ny+1];


    std::ofstream plik ("akcja.dat");
    std::ofstream plik2 ("potencjal.dat");
    std::ofstream plik3 ("blad.dat");

    double x = 0;
    double y = 0;
    for(int i = 0; i<nx+1; i++){
        x = delta*i;
        // plik2 << x;
        y = 0;
        for(int j=0; j<ny+1;j++){
            
            y = delta*j;
            double rho1 = rho_1(x, y);
            double rho2 = rho_2(x, y);
            rho[i][j] = rho_obl(rho1, rho2);
            // plik2 << y << ' ';
        }
        // plik2 << std::endl;
    }

    double blad[nx+1][ny+1];

    for(double omega : omega_tab){
        double it = 0;
        double war;
        for(int i = 1; i<nx; i++){
            for(int j = 1; j<ny; j++){
                Vs[i][j] = 0.0;
            }   
        }

        for(int i = 0; i<nx+1; i++){
            Vs[i][0] = V1;
            Vs[i][ny] = V2;
        }

        do{
            S_old = S(Vs, rho);
            for(int i = 1; i<nx; i++){              
                for(int j = 1; j<ny; j++){                                      
                    Vn[i][j] = 1/4.0*(Vs[i+1][j] + Vs[i-1][j] + Vs[i][j-1] + Vs[i][j+1] + pow(delta,2)/eps*rho[i][j]);
                    Vn[0][j] = Vn[1][j];
                    Vn[nx][j] = Vn[nx-1][j];
                }   
            }

            for(int i = 0; i<nx+1; i++){
                for(int j = 1; j<ny; j++){
                    Vs[i][j] = (1-omega)*Vs[i][j] + omega*Vn[i][j];
                    
                }   
            }

            S_new = S(Vs, rho); 
            // wypisanie całki
            plik << std::fixed << std::setprecision(20) 
            << std::setw(8) << it << ' ' 
            << std::setw(8) << S_new
            << std::endl;
            war = fabs((S_new - S_old)/S_old);
            it++;
        }while(war > TOL);


        for(int i = 1; i<nx; i++){
            for(int j=1; j<ny;j++){
                blad[i][j] = (Vs[i][j+1] - 2*Vs[i][j] + Vs[i][j-1])/pow(delta,2) + (Vs[i+1][j] - 2*Vs[i][j] + Vs[i-1][j])/pow(delta,2) + rho[i][j]/eps;
            }
        }

        // wypisanie potencjału
        for(int i = 0; i<ny+1; i++){
            for(int j=0; j<nx+1;j++){
                plik2 << std::fixed << std::setprecision(8) 
                << std::setw(8) << Vs[j][i] << ' ';
            }
            plik2 << std::endl;
        }
        plik2 << std::endl << std::endl;
        plik << std::endl << std::endl;
           
    
    for(int i = 1; i<ny; i++){
            for(int j=1; j<nx;j++){
                plik3 << std::fixed << std::setprecision(8) 
                << std::setw(8) << blad[j][i] << ' ';
            }
            plik3 << std::endl;
        }
        plik3 << std::endl << std::endl;
    }     
    
    return 0;
}
