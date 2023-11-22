#include <iostream>
#include <math.h>
#include <omp.h>
#include <chrono>
#include<fstream>

using namespace std;

int n=101;
double length = 1;
double ds=(length)/(n-1);
double KinVisc = 1e-2;
double KinVisc_1 = 1.0/KinVisc;
double U_inf=1;
double Re= (length)*(U_inf)/(KinVisc);
double zeta = 1.0/(2.0 * ds), alpha = 2.0 * U_inf / ds;

double** psi;
double** u;
double** v;
double** W;
double** W_old;

void allocatememory()
{
    psi = new double*[n];
    W = new double*[n];
    u = new double*[n];
    v = new double*[n];
    W_old = new double*[n];
    for (int k = 0; k < n; k++)
    {
        psi[k] = new double[n];
        W[k] = new double[n];
        u[k] = new double[n];
        v[k] = new double[n];
        W_old[k] = new double[n];
    }
}

int main()
{
    allocatememory();
    // cout<<"trace1"<<endl;

    // initialzing arrays (W and psi)
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            W[i][j] = 0.0;
            psi[i][j] = 0.0;
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            W_old[i][j] = 0.0;
        }
    }

    // //upper wall boundary condition
    // for(int i=1; i<n; i++)
    // {
    //     u[i][n-1]=U_inf;
    //     v[i][n-1]=0.0;
    // }
    //iteration starts
    int itr;
    double start, end;
    start = omp_get_wtime();
    for(itr=0;itr<1000000; itr++)
    {

        //computing psi
        #pragma omp parallel for
        for(int i=1; i<n-1; i++)
        {
            for(int j=1; j<n-1; j++)
            {
            psi[i][j]=0.25*(psi[i+1][j]+psi[i-1][j]+psi[i][j-1]+psi[i][j+1]+ ds * ds * W[i][j]);
            }
        }

        //computing velocity
        #pragma omp parallel for collapse(2)
        for(int i=1; i<n-1; i++)
        {
            for(int j=1; j<n-1; j++)
            {
                u[i][j]=zeta*(psi[i][j+1]-psi[i][j-1]);
                v[i][j]=-zeta*(psi[i+1][j]-psi[i-1][j]);
                // cout<<u[i][j]<<" "<<v[i][j]<<endl;
            }
        }
        // cout<<"trace 2"<<endl;

        //upper wall boundary condition
        #pragma omp simd
        for(int i=0; i<n; i++)
        {
            u[i][n-1] = U_inf;
            // psi[i][n-1]=U_inf*ds + psi[i-1][n-2];
            W[i][n-1]=(8 * zeta * zeta)*(psi[i][n-1]-psi[i][n-2]) - alpha;
            // psi[i-1][n-2]=-U_inf*ds + psi[i][n-1];
        }

        //left wall BC
        #pragma omp simd
        for(int j=0; j<n; j++)
        {
            // psi[1][j-1]=psi[0][j];
            W[0][j]=(8 * zeta * zeta)*(psi[0][j]-psi[1][j]);
        }

        //right wall BC
        #pragma omp simd
        for(int j=0; j<n; j++)
        {
            // psi[n-2][j]=psi[n-1][j-1];
            W[n-1][j]=(8 * zeta * zeta)*(psi[n-1][j]-psi[n-2][j]);
        }

        //bottom wall BC
        #pragma omp simd
        for(int i=0; i<n; i++)
        {
            // psi[i-1][1]=psi[i][0];
            W[i][0]=(8 * zeta * zeta)*(psi[i][0]-psi[i][1]);
        }

        //computing vorticity transport equation
        //discretizing using upwind scheme & central differencing depedning on local peclet number
        #pragma omp parallel for collapse(2) 
        for(int i=1; i<n-1; i++)
        {
            for(int j=1; j<n-1; j++)
            {
                double A_1x = 0.0, A_2x = 0.0, B_1x = 0.0 , B_2x = 0.0;
                double Re_x = (u[i][j]*ds) * (KinVisc_1);
                double Re_y = (v[i][j]*ds) * (KinVisc_1);
                
                if(Re_x <= 2.0 && Re_x >= 0)
                {
                   A_1x = 1.0/2.0 , A_2x = 1.0/2.0;
                }
                else if(Re_x > 2.0 )
                {
                    A_1x = 0.0 , A_2x = -1.0;
                }
                else if(Re_x<0 && Re_x>= -2.0)
                {
                    A_1x = 1.0/2.0 , A_2x = 1.0/2.0;
                }
                else if(Re_x < -2.0)
                {
                    A_1x = 1.0 , A_2x = 0.0;
                } //if for x direction concludes

                if(Re_y <= 2.0 && Re_y >= 0)
                {
                   B_1x= 1.0/2.0 , B_2x= 1.0/2.0 ;
                }
                else if(Re_y > 2.0 )
                {
                    B_1x= 0.0 , B_2x= -1.0 ;
                }
                else if(Re_y<0 && Re_y>= -2.0)
                {
                    B_1x= 1.0/2.0 , B_2x= 1.0/2.0 ;
                }
                else if(Re_y < -2.0)
                {
                    B_1x= 1.0 , B_2x= 0.0 ;
                } //if for y direction concludes

                double aE = (0.25) - (0.25)*(u[i][j]*A_1x*Re*ds);
                double aW = (0.25) + (0.25)*(u[i][j]*A_2x*Re*ds);
                double aN = (0.25) - (0.25)*(v[i][j]*B_1x*Re*ds);
                double aS = (0.25) + (0.25)*(v[i][j]*B_2x*Re*ds);

                W[i][j] = aE*W[i+1][j] + aW*W[i-1][j] + aN*W[i][j+1] + aS*W[i][j-1]; //vorticity transport equation

            }
        }
        
        double W_error_num=0.0;
        double W_error_den=0.0;
        #pragma omp parallel for collapse(2) reduction(+:W_error_num)
        for(int i=0; i<n; i++)
        {
            for(int j=0; j<n; j++)
            {
                W_error_num += (W_old[i][j] - W[i][j]) * (W_old[i][j] - W[i][j]);
            }
        }
        double W_error = sqrt((W_error_num));

        if(W_error < 1e-7)
        {
            break;
        }

        cout<<"iteration ="<<itr<<"   |   "<<"error ="<<W_error<<endl;


        for(int i=0; i<n; i++)
        {
            for(int j=0; j<n; j++)
            {
                W_old[i][j]=W[i][j];
            }
        }

    }
    end = omp_get_wtime();

    cout<<"\n"<<"time taken = "<<" "<< (end-start)<<endl;


    
ofstream f;
f.open("output_parallel.dat", ios::out);
      for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            f<<i<<" "<<j<<" "<<psi[i][j]<<" "<<W[i][j]<<" "<<u[i][j]<<" "<<v[i][j]<<" "<<sqrt(pow(u[i][j],2)+pow(v[i][j],2))<<endl;
        }
    } 
f.close();
//output file
f.open("U-Re_x=0.5_parallel.dat", ios::out);
     
    double y[n] ;
    y[0]=0;
    for(int j=0; j<n; j++)
    {
        y[j] = y[j-1]+ds;
    }
    for (int j = 0; j < n; j++)
    {
        f<<y[j]<<" "<<u[((n-1)/(2))][j]<<endl;
    }

f.close();
//u file for x=1.0/2.0
f.open("V-Re_x=0.5_parallel.dat", ios::out);
     
     double x[n] ;
    x[0]=0;
    for(int j=0; j<n; j++)
    {
        x[j] = x[j-1]+ds;
    }
    for (int j = 0; j < n; j++)
    {
        f<<x[j]<<" "<<v[((n-1)/(2))][j]<<endl;
    }
f.close();
//v for x=1.0/2.0
}
