#include<iostream>
#include<math.h>
#include<chrono>
#include<omp.h>
#include<fstream>

using namespace std;

int LX=101, LY=101;
int n_iter=1000000;
double length= 1.0;
double delX=(length)/(LX);
double delY=(length)/(LY);
double zeta= pow((delX/delY),2);
double time_parallel;
double time_serial;

double** temp;
double** temp_old;
double** temp_series;
double** temp_parallel;
double** gradTx;
double** gradTy;

void memoryalloc()
{
    temp = new double*[LX];
    temp_old = new double*[LX];
    temp_series = new double*[LX];
    temp_parallel = new double*[LX];
    gradTx = new double*[LX];
    gradTy = new double*[LY];
    for(int i=0; i<LX; i++)
    {
        temp[i]=new double[LY];
        temp_old[i]=new double[LY];
        temp_series[i]=new double[LY];
        temp_parallel[i]=new double[LY];
        gradTx[i]= new double[LY];
        gradTy[i]= new double[LY];
    }
}

void parallel(double** temp)
{
     for(int i=0; i<LX; i++)
    {
        for(int j=0; j<LY; j++)
        {
        temp[i][j]=0;
        temp_old[i][j]=0;
        }
    }

    //initialization
    for(int j=0; j<LY; j++)
    {
        temp[0][j]=100.0;
        temp_old[0][j]=100.0;
    }
    for(int i=0; i<LX; i++)
    {
        temp[i][LY-1]=100.0;
        temp_old[i][LY-1]=100.0;
    }
    for(int i=0; i<LX; i++)
    {
        temp[i][1]=temp[i][0]; 
        temp_old[i][1]=temp_old[i][0]; 
    }
    for(int j=0; j<LY; j++)
    {
        temp[LX-2][j]=temp[LX-1][j];
        temp_old[LX-2][j]=temp_old[LX-1][j];
    }
    

    double start = omp_get_wtime();
    // #pragma omp parallel for collapse(3)
    for(int n=0; n<n_iter; n++)
    {
        #pragma omp parallel for
        for(int i=0; i<LX; i++)
        {
            temp[i][1]=temp[i][0]; 
        }
        #pragma omp parallel for
        for(int j=0; j<LY; j++)
        {
            temp[LX-2][j]=temp[LX-1][j];
        }
        #pragma omp parallel for collapse(2)
        for(int i=1; i<LX-1; i++)
        {
            for(int j=1; j<LY-1; j++)
            {
                temp[i][j] = ((temp[i+1][j]+temp[i-1][j]) + zeta*(temp[i][j+1]+temp[i][j-1]))/(2*(1+zeta));
            }
        }

    float num_sum=0.0;
    float deno_sum=0.0;

        #pragma omp parallel for collapse(2) reduction(+:num_sum, deno_sum)
        for(int i=1; i<LX-1; i++)
		{ 
            for(int j=1; j<LY-1; j++)
            {
                num_sum += pow((temp[i][j]-temp_old[i][j]),2);
                deno_sum += pow((temp_old[i][j]),2);
            }
		}
		float error = sqrt(num_sum/deno_sum);

        cout<<"Parallel iteration ="<<" "<<n<<"   |   "<<error<<endl;
		if (error<1e-12)
		break;
		
        #pragma omp parallel for collapse(2)
		for(int i=1; i<LX-1; i++)
		{
            for(int j=1; j<LY-1; j++)
            {
			    temp_old[i][j]=temp[i][j];
            }
		}
    }
    double end = omp_get_wtime();

    #pragma omp parallel for
    for(int j=1; j<LY-1; j++)
    {
        for(int i=1; i<LX-1; i++)
        {
        gradTx[i][j] = (temp[i+1][j]-temp[i-1][j])/(2*delX); 
        }
    }

    #pragma omp parallel for
     for(int i=1; i<LX-1; i++)
    {
        for(int j=1; j<LY-1; j++)
        {
        gradTy[i][j] = (temp[i][j+1]-temp[i][j-1])/(2*delY); 
        }
    }



    time_parallel = (end-start);
}

void serial(double** temp)
{
    for(int i=0; i<LX; i++)
    {
        for(int j=0; j<LY; j++)
        {
        temp[i][j]=0;
        }
    }

    //initialization
    for(int j=0; j<LY; j++)
    {
        temp[0][j]=100.0;
    }
    for(int i=0; i<LX; i++)
    {
        temp[i][LY-1]=100.0;
    }
    for(int i=0; i<LX; i++)
    {
        temp[i][1]=temp[i][0]; 
    }
    for(int j=0; j<LY; j++)
    {
        temp[LX-2][j]=temp[LX-1][j];
    }

    clock_t start = clock();
    for(int n=0; n<n_iter; n++)
    {
        for(int i=1; i<LX-1; i++)
        {
            for(int j=1; j<LY-1; j++)
            {
                temp[i][j] = ((temp[i+1][j]+temp[i-1][j]) + zeta*(temp[i][j+1]+temp[i][j-1]))/(2*(1+zeta));
            }
        }
    
    float num_sum=0.0;
    float deno_sum=0.0;

        for(int i=1; i<LX-1; i++)
		{ 
            for(int j=1; j<LY-1; j++)
            {
                num_sum += pow((temp[i][j]-temp_old[i][j]),2);
                deno_sum += pow((temp_old[i][j]),2);
            }
		}
		float error = sqrt(num_sum/deno_sum);
        cout<<"Serial iteration ="<<" "<<n<<"   |   "<<error<<endl;
		if (error<1e-12)
        break;

		for(int i=1; i<LX-1; i++)
		{
            for(int j=1; j<LY-1; j++)
            {
			    temp_old[i][j]=temp[i][j];
            }
		}
    }
   
    clock_t end= clock();

    time_serial = (end-start)/double(CLOCKS_PER_SEC);
}

int main()
{
    memoryalloc();

    //CHECKING ARRAY
    // for(int i=0; i<LX; i++)
    // {
    //     for(int j=0; j<LY; j++)
    //     {
    //     temp[i][j]=i+j;
    //     }
    // }

    // for(int i=0; i<LX; i++)
    // {
    //     for(int j=0; j<LY; j++)
    //     {
    //     cout<<i<<" "<<j<<" "<<temp[i][j]<<endl;
    //     }
    // }

    parallel(temp_parallel);
    serial(temp_series);    

    ofstream f;
    f.open("output_parallel.dat", ios::out);
    for(int i=0; i<LX; i++)
    {
        for(int j=0; j<LY; j++)
        {
            // cout<<i<<" "<<j<<" "<<temp_parallel[i][j]<<" "<<temp_series[i][j]<<"   |   "<<"Temp difference= "<<temp_parallel[i][j]-temp_series[i][j]<<endl;
            f<<i<<" "<<j<<" "<<temp_parallel[i][j]<<" "<<gradTx[i][j]<<" "<<gradTy[i][j]<<endl;
        }
        
    }

    fstream g;
    g.open("output_serial.dat", ios::out);
    for(int i=0; i<LX; i++)
    {
        for(int j=0; j<LY; j++)
        {
            // cout<<i<<" "<<j<<" "<<temp_parallel[i][j]<<" "<<temp_series[i][j]<<"   |   "<<"Temp difference= "<<temp_parallel[i][j]-temp_series[i][j]<<endl;
            g<<i<<" "<<j<<" "<<temp_series[i][j]<<" "<<gradTx[i][j]<<" "<<gradTy[i][j]<<endl;
        }
        
    }
    cout<<"\n"<<endl;
    cout<<"Parallel runtime="<<" "<<time_parallel<<" "<<   "|"   <<"Serial runtime="<<" "<<time_serial<<endl;
//     cout<<"\n"<<endl;
//     cout<<"Serial runtime="<<" "<<time_serial<<endl;
}
