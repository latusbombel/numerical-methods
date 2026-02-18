#include <iostream>
#include <cmath>
#include "mgmres.h"
#include <fstream>
#include <iomanip>

int eps1 = 1;
int eps2 = 10;

int eps(int l){
    int nx = 100;
    int j = floor(l/(nx+1));
    int i = l - j * (nx+1);
    if(i <= nx/2){
        return eps1;
    }
    else return eps2;
}

double rho1(int i, int j, double delta, double x_max, double y_max){
    double x = i * delta;
    double y = j * delta;
    double sigma = x_max/10.0;
    return exp(-(x-0.25*x_max)*(x-0.25*x_max)/(sigma*sigma) - (y-0.5*y_max)*(y-0.5*y_max)/(sigma*sigma));
}

double rho2(int i, int j, double delta, double x_max, double y_max){
    double x = i * delta;
    double y = j * delta;
    double sigma = x_max/10.0;
    return  (-1)*exp(-(x-0.75*x_max)*(x-0.75*x_max)/(sigma*sigma) - (y-0.5*y_max)*(y-0.5*y_max)/(sigma*sigma));
}


int main(){
    double delta = 0.1, V1 = 0., V2 = 0., V3 = 0., V4 = 0.;
    int nx = 100, ny = 100;
    double x_max = nx*delta;
    double y_max = ny*delta;
    // double rho1 = 0., rho2 = 0.;
    int N = (nx+1)*(ny+1);
    double* a = new double[5*N];
    int* ja = new int[5*N];
    int* ia = new int[N+1];

    for(int i = 0; i<N+1; i++){
        ia[i] = -1;
    }

    double b[N];
    double V[N];
    int k = -1;

    //zapis do plików

    std::ofstream plik ("macierz_A.dat");
    std::ofstream plik2 ("wektor_B.dat");
    std::ofstream plik3 ("potencjal.dat");

    //UZALEŻNIĆ EPS OD L

    for(int l=0; l<N; l++){
        int brzeg =0; // wskaźnik położenia: 0-środek obszaru; 1-brzeg
        double vb=0; // potencjal na brzegu
        int j = floor(l/(nx+1));
        int i = l - j * (nx+1);
        if(i==0){ //lewy brzeg
            brzeg=1;
            vb=V1;
        }
        if(j==ny){ //górny brzeg
            brzeg=1;
            vb=V2;
        }
        if(i==nx){ //prawy brzeg
            brzeg=1;
            vb=V3;
        }
        if(j==0) { //dolny brzeg
            brzeg=1;
            vb=V4;
        }
        // wypełniamy od razu wektor wyrazów wolnych
        b[l]=-(rho1(i, j, delta, x_max, y_max) + rho2(i, j, delta, x_max, y_max)); //jeśli w środku jest gęstość
        if(brzeg ==1) {
        b[l]=vb; // wymuszamy wartość potencjału na brzegu
        }
        // wypełniamy elementy macierzy A
        
        ia[l]=-1; // wskaźnik dla pierwszego el. w wierszu
        //lewa skrajna przekatna
        if(l-nx -1>=0 && brzeg ==0 ){
            k++;
            if(ia[l]<0) ia[l]=k;
            a[k]= eps(l)/(delta*delta);
            ja[k]=l - nx - 1;
        }
        // poddiagonala
        if(l-1>=0 && brzeg ==0 ){
            k++;
            if(ia[l]<0)ia[l]=k;
            a[k]=eps(l)/(delta*delta);
            ja[k]=l - 1;
        }
        // diagonala
        k++;
        if(ia[l]<0)ia[l]=k;
        if(brzeg ==0) {
            a[k]= - (2*eps(l)+eps(l+1)+eps(l+nx+1))/(delta*delta);}
        else{
            a[k]=1;
        }
        ja[k]=l;
        // naddiagonala
        if(l<N && brzeg ==0 ){
            k++;
            a[k]=eps(l+1)/(delta*delta);
            ja[k]= l + 1;
        }
        // prawa skrajna przekątna
        if(l<N -nx -1 && brzeg ==0 ){
            k++;
            a[k]=eps(l+nx+1)/(delta*delta);
            ja[k]=l + nx + 1;
        }

        plik << std::fixed << std::setprecision(4) 
            << std::setw(2) << l << ' ' 
            << std::setw(2) << i << ' '
            << std::setw(2) << j << ' '
            << std::setprecision(8) 
            << std::setw(8) << a[l]
            << std::endl;

        plik2 << std::fixed << std::setprecision(4) 
            << std::setw(2) << l << ' ' 
            << std::setw(2) << i << ' '
            << std::setw(2) << j << ' '
            << std::setprecision(8) 
            << std::setw(8) << b[l]
            << std::endl;

    }
    int nz_num = k+1;
    ia[N] = nz_num;
    int itr_max = 500;
    int mr = 500;
    double total_abs = 10e-8, tot_rel = 10e-8;

    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, total_abs, tot_rel);

    for (int l = 0; l<N;l++){
        int j = floor(l/(nx+1));
        int i = l - j * (nx+1);
        plik3 << std::fixed << std::setprecision(4) 
            << std::setw(2) << l << ' ' 
            << std::setw(2) << i << ' '
            << std::setw(2) << j << ' '
            << std::setprecision(8) 
            << std::setw(8) << V[l]
            << std::endl;
    }


    delete[] a;
    delete[] ja;
    delete[] ia;
    
    return 0;
}
