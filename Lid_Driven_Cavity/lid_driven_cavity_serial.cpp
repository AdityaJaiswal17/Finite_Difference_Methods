#include <iostream>
#include <math.h>
#include <chrono>
#include<fstream>

using namespace std;

int n=101;
double length = 1;
double ds=(length)/(n-1);
double KinVisc = 0.01;
double U_inf=1.0;
double Re= (length)*(U_inf)/(KinVisc);

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
    clock_t start, end;
    start = clock();
    for(itr=0;itr<1000000; itr++)
    {

        //computing psi
        for(int i=1; i<n-1; i++)
        {
            for(int j=1; j<n-1; j++)
            {
            psi[i][j]=(1.0/4.0)*(psi[i+1][j]+psi[i-1][j]+psi[i][j-1]+psi[i][j+1]+ pow(ds,2)*W[i][j]);
            }
        }

        //computing velocity
        for(int i=1; i<n-1; i++)
        {
            for(int j=1; j<n-1; j++)
            {
                u[i][j]=(1.0/(2.0*ds))*(psi[i][j+1]-psi[i][j-1]);
                v[i][j]=(-1.0/(2.0*ds))*(psi[i+1][j]-psi[i-1][j]);
                // cout<<u[i][j]<<" "<<v[i][j]<<endl;
            }
        }
        // cout<<"trace 2"<<endl;

        //upper wall boundary condition
        for(int i=0; i<n; i++)
        {
            u[i][n-1] = U_inf;
            // psi[i][n-1]=U_inf*ds + psi[i-1][n-2];
            W[i][n-1]=(2.0/pow(ds,2))*(psi[i][n-1]-psi[i][n-2]) - (2.0*U_inf/ds);
            // psi[i-1][n-2]=-U_inf*ds + psi[i][n-1];
        }

        //left wall BC
        for(int j=0; j<n; j++)
        {
            // psi[1][j-1]=psi[0][j];
            W[0][j]=(2.0/pow(ds,2))*(psi[0][j]-psi[1][j]);
        }

        //right wall BC
        for(int j=0; j<n; j++)
        {
            // psi[n-2][j]=psi[n-1][j-1];
            W[n-1][j]=(2.0/pow(ds,2))*(psi[n-1][j]-psi[n-2][j]);
        }

        //bottom wall BC
        for(int i=0; i<n; i++)
        {
            // psi[i-1][1]=psi[i][0];
            W[i][0]=(2.0/pow(ds,2))*(psi[i][0]-psi[i][1]);
        }

        double K = 2.0/(Re * ds);

        //computing vorticity transport equation
        //discretizing using upwind scheme & central differencing depedning on local peclet number
        for(int i=1; i<n-1; i++)
        {
            for(int j=1; j<n-1; j++)
            {
                double A_1x = 0.0, A_2x = 0.0, B_1x = 0.0 , B_2x = 0.0;
                double Re_x = (u[i][j]*ds)/(KinVisc);
                double Re_y = (v[i][j]*ds)/(KinVisc);
                
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

                double aE = (1.0/4.0) - (1.0/4.0)*(u[i][j]*A_1x*Re*ds);
                double aW = (1.0/4.0) + (1.0/4.0)*(u[i][j]*A_2x*Re*ds);
                double aN = (1.0/4.0) - (1.0/4.0)*(v[i][j]*B_1x*Re*ds);
                double aS = (1.0/4.0) + (1.0/4.0)*(v[i][j]*B_2x*Re*ds);

                W[i][j] = aE*W[i+1][j] + aW*W[i-1][j] + aN*W[i][j+1] + aS*W[i][j-1]; //vorticity transport equation

            }
        }
        
        double W_error_num=0.0;
        double W_error_den=0.0;
        for(int i=0; i<n; i++)
        {
            for(int j=0; j<n; j++)
            {
                W_error_num += pow(fabs(W_old[i][j] - W[i][j]),2);
            }
        }
        double W_error = sqrt((W_error_num));

        cout<<"iteration ="<<itr<<"   |   "<<"error ="<<W_error<<endl;

        if(W_error < 1e-7)
        {
            break;
        }

        for(int i=0; i<n; i++)
        {
            for(int j=0; j<n; j++)
            {
                W_old[i][j]=W[i][j];
            }
        }

    }
    end = clock();

    cout<<"\n"<<"time taken = "<<" "<< ((end-start)/double(CLOCKS_PER_SEC))<<endl;


    
ofstream f;
f.open("output_serial.dat", ios::out);
      for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            f<<i<<" "<<j<<" "<<psi[i][j]<<" "<<W[i][j]<<" "<<u[i][j]<<" "<<v[i][j]<<" "<<sqrt(pow(u[i][j],2)+pow(v[i][j],2))<<endl;
        }
    } 
f.close();
//output file
f.open("U-Re_x=0.5_serial.dat", ios::out);
     
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
//u file for x=0.5
f.open("V-Re_x=0.5_serial.dat", ios::out);
     
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
//v for x=0.5
}
