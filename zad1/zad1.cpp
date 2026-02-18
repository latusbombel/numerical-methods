#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

int main(){
    int tmax = 5, tmin = 0, i, j;
    int n[3];
    double t[3] = {0.01, 0.1, 1.0};
    double y0 = 1, lambda = -1;

    std::vector<std::vector<double>> y{{y0}, {y0}, {y0}};
    /*std::vector<double> y1{y0};
    std::vector<double> y2{y0};*/

    std::ofstream plik("wyniki1.txt");
    
    
    for(j = 0; j<3; j++){
        int krok = 1;
        do{
            y[j].push_back(y[j][krok-1] + lambda*t[j]*y[j][krok-1]);
            krok++;
        }while(krok*t[j]<=tmax);
    }
    
    for(j = 0; j<3; j++){
        for(i=0; i<y[j].size(); i++){
            double y_anal = exp(lambda * i * t[j]);
            plik << std::fixed << std::setprecision(4)
                << std::setw(6) << i*t[j] << ' '
                << std::setw(6) <<  y[j][i] << ' ' 
                << std::setw(6) << y_anal << ' ' << std::endl;
        }
        plik << std::endl << std::endl;
    }

    
    std::ofstream plik1 ("wyniki2.dat");


    for(j = 0; j<3; j++){
        for(i=0; i<y[j].size(); i++){
            double y_anal = exp(lambda * i * t[j]);
            plik1 << std::fixed << std::setprecision(4)
                << std::setw(6) << i*t[j] << ' '
                << std::setw(6) <<  y[j][i]-y_anal << std::endl;
        }
        plik1 << std::endl << std::endl;
    }
    
    plik1.close();
    plik.close();
    return 0;
}