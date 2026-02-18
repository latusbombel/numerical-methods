#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <algorithm>
#define M_PI 3.14159265358979323846
#define G 4*M_PI*M_PI

double ax(double x, double y){
    double x0 = 0;
    double y0 = 0;
    return -1*G/pow(((x-x0)*(x-x0)+ (y-y0)*(y-y0)),3/2.0)*(x-x0);
}

double ay(double x, double y){
    double x0 = 0;
    double y0 = 0;
    return -1*G/pow(((x-x0)*(x-x0)+ (y-y0)*(y-y0)),3/2.0)*(y-y0);
}

double axx(double x, double y){
    return -G*(y*y - 2*x*x)/pow((x*x + y*y),5/2.0);
}

double ayy(double x, double y){
    return -G*(x*x - 2*y*y)/pow((x*x + y*y),5/2.0);
}

double axy(double x, double y){
    return 3*G*x*y/pow((x*x + y*y),5/2.0);
}

void oblicz_U(double Ux[], double Uy[], double Uvx[], double Uvy[], double A[2][2], double dt, double x, double y, double vx, double vy){
    double Fx[2], Fy[2], Fvx[2], Fvy[2];
    double TOL = pow(10, -6);
    int iter = 0;
    int MaxIter = 200;
    do{
        Fx[0] = Ux[0] - x - dt*(A[0][0]*Uvx[0] + A[0][1]*Uvx[1]);
        Fy[0] = Uy[0] - y - dt*(A[0][0]*Uvy[0] + A[0][1]*Uvy[1]);

        Fvx[0] = Uvx[0] - vx - dt*(A[0][0]*ax(Ux[0], Uy[0]) + A[0][1]*ax(Ux[1],Uy[1]));
        Fvy[0] = Uvy[0] - vy - dt*(A[0][0]*ay(Ux[0], Uy[0]) + A[0][1]*ay(Ux[1],Uy[1]));

        Fx[1] = Ux[1] - x - dt*(A[1][0]*Uvx[0] + A[1][1]*Uvx[1]);
        Fy[1] = Uy[1] - y - dt*(A[1][0]*Uvy[0] + A[1][1]*Uvy[1]);    

        Fvx[1] = Uvx[1] - vx - dt*(A[1][0]*ax(Ux[0], Uy[0]) + A[1][1]*ax(Ux[1],Uy[1]));
        Fvy[1] = Uvy[1] - vy - dt*(A[1][0]*ay(Ux[0], Uy[0]) + A[1][1]*ay(Ux[1],Uy[1]));
        double U[8] = {Ux[0], Uvx[0], Uy[0], Uvy[0], Ux[1], Uvx[1], Uy[1], Uvy[1]};
        double F[8] = {-1*Fx[0], -1*Fvx[0], -1*Fy[0], -1*Fvy[0], -1*Fx[1], -1*Fvx[1], -1*Fy[1], -1* Fvy[1]};
        double dU[8];
        double J[8][8] = {
            // Jakobian J[i][j] = dF_i/dU_j
                    {1, -dt*A[0][0], 0, 0, 0, -dt*A[0][1], 0, 0},
                    {-dt*A[0][0]*axx(Ux[0],Uy[0]), 1, -dt*A[0][0]*axy(Ux[0],Uy[0]), 0, -dt*A[0][1]*axx(Ux[1], Uy[1]), 0, -dt*A[0][1]*axy(Ux[1], Uy[1]), 0},
                    {0, 0, 1, -dt*A[0][0], 0, 0, 0, -dt*A[0][1]},
                    {-dt*A[0][0]*axy(Ux[0],Uy[0]), 0, -dt*A[0][0]*ayy(Ux[0],Uy[0]), 1, -dt*A[0][1]*axy(Ux[1],Uy[1]), 0, -dt*A[0][1]*ayy(Ux[1],Uy[1]), 0},
                    {0, -dt*A[1][0], 0, 0, 1, -dt*A[1][1], 0, 0},
                    {-dt*A[1][0]*axx(Ux[0],Uy[0]), 0, -dt*A[1][0]*axy(Ux[0],Uy[0]), 0, -dt*A[1][1]*axx(Ux[1],Uy[1]), 1, -dt*A[1][1]*axy(Ux[1],Uy[1]), 0},
                    {0, 0, 0, -dt*A[1][0], 0, 0, 1, -dt*A[1][1]},
                    {-dt*A[1][0]*axy(Ux[0],Uy[0]), 0, -dt*A[1][0]*ayy(Ux[0],Uy[0]), 0, -dt*A[1][1]*axy(Ux[1],Uy[1]), 0, -dt*A[1][1]*ayy(Ux[1],Uy[1]), 1},
                };
        
                // Używamy Eigen do rozwiązania układu równań
        Eigen::Matrix<double, 8, 8> J_mat;
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                J_mat(i, j) = J[i][j];
            }
        }
        Eigen::VectorXd F_vec(8);
        for (int i = 0; i < 8; i++) {
            F_vec(i) = F[i];
        }

        Eigen::VectorXd dU_vec = J_mat.fullPivLu().solve(F_vec);

        for (int i = 0; i < 8; i++) {
            dU[i] = dU_vec(i);
        }

        //Gauss_Elimination(J, F, dU);

        for(int i = 0; i<8; i++){
            U[i] += dU[i];
        }

        double maxVal = dU[0];

        for(int i = 1; i<8; i++){
            if(fabs(dU[i]) > fabs(maxVal)){
                maxVal = dU[i];
            }
        }

        Ux[0] = U[0]; 
        Uvx[0] = U[1];
        Uy[0] = U[2];
        Uvy[0] = U[3];
        Ux[1] =  U[4];
        Uvx[1] = U[5];
        Uy[1] = U[6];
        Uvy[1] = U[7];
        iter++;
        if(fabs(maxVal) < TOL) break;
    }while(iter < MaxIter);   
}

void NRK2(double A[2][2], double b1, double b2, double dt, double x, double y, double vx, double vy, double& x_new, double& y_new, double& vx_new, double& vy_new){
    double Ux[2], Uy[2], Uvx[2], Uvy[2];

    Ux[0] = x;
    Ux[1] = x;
    Uy[0] = y;
    Uy[1] = y;

    Uvx[0] = vx;
    Uvx[1] = vx;
    Uvy[0] = vy;
    Uvy[1] = vy;

    oblicz_U(Ux, Uy, Uvx, Uvy, A, dt, x, y, vx ,vy);

    x_new = x + dt*(b1*Uvx[0] + b2*Uvx[1]);
    y_new = y + dt*(b1*Uvy[0] + b2*Uvy[1]);

    vx_new = vx + dt*(b1*ax(Ux[0], Uy[0]) + b2*ax(Ux[1], Uy[1]));
    vy_new = vy + dt*(b1*ay(Ux[0], Uy[0]) + b2*ay(Ux[1], Uy[1]));
}

int main(){
    double A[2][2] = {
            {1/4.0, 1/4.0 - sqrt(3)/6},
            {1/4.0 + sqrt(3)/6, 1/4.0}
                };
    
    double x = 1, y = 0, vx = 0, vy = sqrt(G);
    double dt = 0.01;
    double t_max = 2.0;
    double S = 0.75;
    double tol = pow(10,-6);
    double b1 = 0.5, b2 = 0.5;
    double t = 0;

    long double E = (vx*vx + vy*vy)/2 - G/sqrt(x*x+y*y);

    
    long double Ex, Ey, Evx, Evy;
    std::ofstream plik ("wyniki2.dat"); 
    do{
        double x1_temp, y1_temp, vx1_temp, vy1_temp;
        double x2_temp, y2_temp, vx2_temp, vy2_temp;
        double x2_new, y2_new, vx2_new, vy2_new;
        
        NRK2(A, b1, b2, 2*dt, x, y, vx, vy, x1_temp, y1_temp, vx1_temp, vy1_temp);

        NRK2(A, b1, b2, dt, x, y, vx, vy, x2_temp, y2_temp, vx2_temp, vy2_temp);
        NRK2(A, b1, b2, dt, x2_temp, y2_temp, vx2_temp, vy2_temp, x2_new, y2_new, vx2_new, vy2_new);
        
        Ex = (x2_new - x1_temp)/(pow(2, 4) - 1);
        Ey = (y2_new - y1_temp)/(pow(2, 4) - 1);
        Evx = (vx2_new - vx1_temp)/(pow(2, 4) - 1);
        Evy = (vy2_new - vy1_temp)/(pow(2, 4) - 1);
        long double E[4] = {Ex, Ey, Evx, Evy};
        long double maxVal = E[0];

        for(int j = 0; j<4; j++){
            if(fabs(E[j]) > fabs(maxVal)){
                maxVal = fabs(E[j]);
            }
        }

        if(maxVal < tol){
            x = x2_new;
            y = y2_new;
            vx = vx2_new;
            vy = vy2_new;
            t += 2*dt;
            long double E2 = (vx*vx + vy*vy)/2 - G/sqrt(x*x+y*y);
            double L = std::abs(x*vy - y*vx);
            plik << std::fixed << std::setprecision(20) 
            << std::setw(8) << t << ' ' 
            << std::setw(8) << x << ' '
            << std::setw(8) << y << ' '
            << std::setw(8) << vx << ' '
            << std::setw(8) << vy << ' '
            << std::setw(20) << E2 << ' '
            << std::setw(20) << L << ' '
            << std::setw(8) << dt << std::endl;
        }

        dt = pow((S*tol/maxVal),1/5.0)*dt;
    }while(t <= t_max);
    


    plik.close();
    return 0;
}