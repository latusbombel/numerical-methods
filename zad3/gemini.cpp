#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <algorithm>

#define M_PI 3.14159265358979323846
#define G (4*M_PI*M_PI)

// --- Funkcje Sił i Pochodnych ---
double ax(double x, double y){
    double denom = pow(x*x + y*y, 1.5);
    return -G * x / denom;
}

double ay(double x, double y){
    double denom = pow(x*x + y*y, 1.5);
    return -G * y / denom;
}

double axx(double x, double y){
    double r2 = x*x + y*y;
    return -G * (y*y - 2*x*x) / pow(r2, 2.5);
}

double ayy(double x, double y){
    double r2 = x*x + y*y;
    return -G * (x*x - 2*y*y) / pow(r2, 2.5);
}

double axy(double x, double y){
    double r2 = x*x + y*y;
    return 3 * G * x * y / pow(r2, 2.5);
}

// --- Implicit RK2 Solver (Newton Iteration) ---
void oblicz_U(double Ux[], double Uy[], double Uvx[], double Uvy[], double A[2][2], double dt, double x, double y, double vx, double vy){
    double TOL = 1e-6; // Zwiększona precyzja dla Newtona
    int iter = 0;
    int MaxIter = 50; // 200 to zazwyczaj za dużo, Newton jest szybki

    // Inicjalizacja (Guess) - używamy wartości z poprzedniego kroku
    Ux[0] = x; Ux[1] = x;
    Uy[0] = y; Uy[1] = y;
    Uvx[0] = vx; Uvx[1] = vx;
    Uvy[0] = vy; Uvy[1] = vy;

    while(iter < MaxIter){
        // Obliczanie F (Residuum)
        // F = U - w_n - dt * sum(A * k)
        // Zauważ: W definicji problemu k to pochodne, czyli v dla pozycji i a dla prędkości.
        
        double F[8];
        // Równania dla X1, Y1
        F[0] = Ux[0] - x - dt*(A[0][0]*Uvx[0] + A[0][1]*Uvx[1]);
        F[2] = Uy[0] - y - dt*(A[0][0]*Uvy[0] + A[0][1]*Uvy[1]);
        
        // Równania dla Vx1, Vy1
        F[1] = Uvx[0] - vx - dt*(A[0][0]*ax(Ux[0], Uy[0]) + A[0][1]*ax(Ux[1], Uy[1]));
        F[3] = Uvy[0] - vy - dt*(A[0][0]*ay(Ux[0], Uy[0]) + A[0][1]*ay(Ux[1], Uy[1]));

        // Równania dla X2, Y2
        F[4] = Ux[1] - x - dt*(A[1][0]*Uvx[0] + A[1][1]*Uvx[1]);
        F[6] = Uy[1] - y - dt*(A[1][0]*Uvy[0] + A[1][1]*Uvy[1]);
        
        // Równania dla Vx2, Vy2
        F[5] = Uvx[1] - vx - dt*(A[1][0]*ax(Ux[0], Uy[0]) + A[1][1]*ax(Ux[1], Uy[1]));
        F[7] = Uvy[1] - vy - dt*(A[1][0]*ay(Ux[0], Uy[0]) + A[1][1]*ay(Ux[1], Uy[1]));

        // Jakobian J[i][j] = dFi / dUj
        // Kolejność zmiennych w U: Ux1, Uvx1, Uy1, Uvy1, Ux2, Uvx2, Uy2, Uvy2
        // Ale uwaga: Twój kod używa kolejności indeksów [0..7] mapowanej na {Ux[0], Uvx[0], Uy[0], Uvy[0]...}
        // Utrzymuję Twoją konwencję mapowania wektora U i F:
        // 0:x1, 1:vx1, 2:y1, 3:vy1, 4:x2, 5:vx2, 6:y2, 7:vy2

        double J[8][8] = {0}; // Zerowanie

        // Blok 1 (dla etapu 1)
        J[0][0]=1; J[0][1]=-dt*A[0][0]; J[0][5]=-dt*A[0][1]; // dFx1

        J[1][0] = -dt*A[0][0]*axx(Ux[0],Uy[0]); 
        J[1][1] = 1; 
        J[1][2] = -dt*A[0][0]*axy(Ux[0],Uy[0]);
        J[1][4] = -dt*A[0][1]*axx(Ux[1],Uy[1]);
        J[1][6] = -dt*A[0][1]*axy(Ux[1],Uy[1]);

        J[2][2]=1; J[2][3]=-dt*A[0][0]; J[2][7]=-dt*A[0][1]; // dFy1

        J[3][0] = -dt*A[0][0]*axy(Ux[0],Uy[0]);
        J[3][2] = -dt*A[0][0]*ayy(Ux[0],Uy[0]);
        J[3][3] = 1;
        J[3][4] = -dt*A[0][1]*axy(Ux[1],Uy[1]);
        J[3][6] = -dt*A[0][1]*ayy(Ux[1],Uy[1]);

        // Blok 2 (dla etapu 2)
        J[4][1]=-dt*A[1][0]; J[4][4]=1; J[4][5]=-dt*A[1][1]; // dFx2

        J[5][0] = -dt*A[1][0]*axx(Ux[0],Uy[0]);
        J[5][2] = -dt*A[1][0]*axy(Ux[0],Uy[0]);
        J[5][5] = 1;
        J[5][4] = -dt*A[1][1]*axx(Ux[1],Uy[1]);
        J[5][6] = -dt*A[1][1]*axy(Ux[1],Uy[1]);

        J[6][3]=-dt*A[1][0]; J[6][6]=1; J[6][7]=-dt*A[1][1]; // dFy2

        J[7][0] = -dt*A[1][0]*axy(Ux[0],Uy[0]);
        J[7][2] = -dt*A[1][0]*ayy(Ux[0],Uy[0]);
        J[7][4] = -dt*A[1][1]*axy(Ux[1],Uy[1]);
        J[7][6] = -dt*A[1][1]*ayy(Ux[1],Uy[1]);
        J[7][7] = 1;

        // Rozwiązanie układu J * dU = -F
        Eigen::Matrix<double, 8, 8> J_mat;
        Eigen::VectorXd F_vec(8);
        for(int i=0; i<8; i++) {
            F_vec(i) = -F[i]; // Przeniesienie F na prawą stronę
            for(int j=0; j<8; j++) J_mat(i,j) = J[i][j];
        }

        Eigen::VectorXd dU = J_mat.fullPivLu().solve(F_vec);

        double max_dU = 0;
        for(int i=0; i<8; i++){
            // Update zmiennych
            double val = dU(i);
            if(i==0) Ux[0]+=val; else if(i==1) Uvx[0]+=val;
            else if(i==2) Uy[0]+=val; else if(i==3) Uvy[0]+=val;
            else if(i==4) Ux[1]+=val; else if(i==5) Uvx[1]+=val;
            else if(i==6) Uy[1]+=val; else if(i==7) Uvy[1]+=val;
            
            if(std::abs(val) > max_dU) max_dU = std::abs(val);
        }

        iter++;
        if(max_dU < TOL) break;
    }
}

// Funkcja wykonująca jeden krok RK2 (zwraca wynik przez referencję)
void perform_step_IRK2(double dt, double x, double y, double vx, double vy, 
                      double &x_out, double &y_out, double &vx_out, double &vy_out,
                      double A[2][2], double b1, double b2) {
    
    double Ux[2], Uy[2], Uvx[2], Uvy[2];
    oblicz_U(Ux, Uy, Uvx, Uvy, A, dt, x, y, vx, vy);

    // Korektor (Wzór 33-36 w instrukcji lub standardowy wzór RK)
    // x_n+1 = x_n + dt * sum(b_i * k_i)
    // gdzie k_i to prędkości (dla x) lub przyspieszenia (dla v)
    x_out = x + dt * (b1 * Uvx[0] + b2 * Uvx[1]);
    y_out = y + dt * (b1 * Uvy[0] + b2 * Uvy[1]);
    vx_out = vx + dt * (b1 * ax(Ux[0], Uy[0]) + b2 * ax(Ux[1], Uy[1]));
    vy_out = vy + dt * (b1 * ay(Ux[0], Uy[0]) + b2 * ay(Ux[1], Uy[1]));
}

int main(){
    // Parametry metody IRK2 (Gauss-Legendre 4th order)
    double A[2][2] = {
            {1.0/4.0, 1.0/4.0 - sqrt(3)/6.0},
            {1.0/4.0 + sqrt(3)/6.0, 1.0/4.0}
    };
    double b1 = 0.5, b2 = 0.5;

    // Warunki początkowe
    double x = 1.0, y = 0.0;
    double vx = 0.0, vy = sqrt(G);
    double t = 0.0;
    double t_max = 2.0;

    // Parametry kontroli kroku
    double dt = 0.01;
    double S = 0.75;
    double tol = 1e-6;
    double p = 4.0; // Rząd metody

    // Plik wyjściowy
    std::ofstream plik("wyniki_irk2.dat");
    
    // Zapis stanu początkowego
    double E_tot = (vx*vx + vy*vy)/2.0 - G/sqrt(x*x + y*y);
    double L = std::abs(x*vy - y*vx);
    plik << std::fixed << std::setprecision(6) 
         << t << "\t" << x << "\t" << y << "\t" << vx << "\t" << vy << "\t" 
         << E_tot << "\t" << L << "\t" << dt << "\n";

    std::cout << "Start symulacji..." << std::endl;

    while(t < t_max){
        // 1. Wykonaj dwa kroki z dt
        double x1, y1, vx1, vy1; // Po pierwszym małym kroku
        double x2, y2, vx2, vy2; // Po drugim małym kroku
        
        perform_step_IRK2(dt, x, y, vx, vy, x1, y1, vx1, vy1, A, b1, b2);
        perform_step_IRK2(dt, x1, y1, vx1, vy1, x2, y2, vx2, vy2, A, b1, b2);

        // 2. Wykonaj jeden duży krok z 2*dt
        double x_big, y_big, vx_big, vy_big;
        perform_step_IRK2(2*dt, x, y, vx, vy, x_big, y_big, vx_big, vy_big, A, b1, b2);

        // 3. Oblicz błąd (Wzór 52)
        // Używamy abs() żeby uniknąć błędów z liczbami ujemnymi!
        double Ex = std::abs(x2 - x_big) / (pow(2, p) - 1);
        double Ey = std::abs(y2 - y_big) / (pow(2, p) - 1);
        double Evx = std::abs(vx2 - vx_big) / (pow(2, p) - 1);
        double Evy = std::abs(vy2 - vy_big) / (pow(2, p) - 1);

        double max_error = std::max({Ex, Ey, Evx, Evy});

        // Zabezpieczenie przed dzieleniem przez zero
        if(max_error < 1e-16) max_error = 1e-16;

        // 4. Decyzja
        if(max_error < tol){
            // AKCEPTACJA
            t += 2*dt;
            x = x2;
            y = y2;
            vx = vx2;
            vy = vy2;

            // Zapis wyników
            E_tot = (vx*vx + vy*vy)/2.0 - G/sqrt(x*x + y*y);
            L = std::abs(x*vy - y*vx);
            plik << t << "\t" << x << "\t" << y << "\t" << vx << "\t" << vy << "\t" 
                 << E_tot << "\t" << L << "\t" << dt << "\n";
                 
            // std::cout << "t: " << t << " dt: " << dt << std::endl; 
        } else {
            // ODSTAWIENIE (REJECTION)
            // Nie zmieniamy t, x, y... po prostu pętla przejdzie jeszcze raz z mniejszym dt
            // std::cout << "Krok odrzucony przy t=" << t << ", error=" << max_error << std::endl;
        }

        // 5. Nowy krok czasowy (Wzór 51)
        double new_dt = dt * pow( (S * tol / max_error), 1.0/(p+1) );
        
        // Opcjonalne ograniczniki, żeby dt nie zmieniało się zbyt gwałtownie
        // if(new_dt > 5*dt) new_dt = 5*dt;
        // if(new_dt < 0.1*dt) new_dt = 0.1*dt;
        
        dt = new_dt;
        
        // Zabezpieczenie przed wyjściem poza t_max w ostatnim kroku
        if(t + 2*dt > t_max && t < t_max) {
             // Można tu dostosować dt, ale zazwyczaj algorytm sam dojedzie
        }
    }

    plik.close();
    std::cout << "Koniec. Wyniki w wyniki_irk2.dat" << std::endl;
    return 0;
}