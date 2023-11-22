#include<iostream>
#include<chrono>
#include<time.h>
#include<math.h>
#include<omp.h>
#include<fstream>

using namespace std;

int LX=501;
double length = 1;
double density=1, thermal_cond=100, specific_heat=10;
double delt=1;
double delX=length/LX;
int n_serial = 0.0, n_parallel = 0.0;
double time_serial = 0.0, time_parallel = 0.0;
int n_iter=100000; //number of iterations

// void parameters(void)
// {
// 	double *temp = new double[LX];
// 	double *temp_old= new double[LX];
//  	double alpha=(density*specific_heat)/(thermal_cond);
// 	double error=0;
// 	double num_sum=0.0 ,deno_sum=0.0;
// }

void parallel(double* temp)
{
	double *temp_old= new double[LX];
 	double alpha=(density*specific_heat)/(thermal_cond);
	double error=0;
	double num_sum=0.0 ,deno_sum=0.0;

 	//test for array population 

 	// for(int i=0; i<LX; i++)
 	// {
 	// 	temp[i]=i;
 	// }
 	// for(int i=0; i<LX; i++)
 	// {
 	// 	std::cout<<temp[i]<<"\n"<<std::endl;
 	// }

	//initializaiton

 	for(int i=0; i<LX; i++)
 	{
 		temp[i]=0;
		temp_old[i]=0;
 	}

 	temp[0]=100.0;
 	temp[LX-1]=1000.0;
	temp_old[0]=100.0;
 	temp_old[LX-1]=1000.0;

 	// for(int n_iter=0; n_iter<10000; n_iter++)
	// {
	// 	for(int i=1; i<LX-1; i++)
	// 	{
	// 		temp[i]=0.5*(temp[i+1]+temp[i-1]);
	// 	}

	// }
	int n_iter=0;
 	double start = omp_get_wtime();
	// clock_t start = clock();
	
	for(n_iter=0; n_iter<1000000; n_iter++)
	{
	    num_sum = 0.0;
		deno_sum = 0.0;
		#pragma omp parallel for
		for(int i=1; i<LX-1; i++)
		{
			temp[i]=0.5*(temp[i+1]+temp[i-1]);
		}

		#pragma omp parallel for reduction(+:num_sum, deno_sum)
		for(int j=1; j<LX-1; j++)
		{ 
			num_sum += pow((temp[j]-temp_old[j]),2);
			deno_sum += pow((temp_old[j]),2);
		}
		error = sqrt(num_sum/deno_sum);
		if (error<1e-12)
		break;
		
		#pragma omp parallel for
		for(int k=1; k<LX-1; k++)
		{
			temp_old[k]=temp[k];
		}
	}
	double end = omp_get_wtime();
	// clock_t end = clock();

	// for(int i=0; i<LX; i++)
	// {
	// 	cout<<i<<" "<<temp[i]<<endl;
	// }
	// cout<<"\n"<<"parallel iteration no. ="<<" "<<n_iter<<endl;
	// cout<<endl<<"Parallel Run time = "<<(end-start)/double(CLOCKS_PER_SEC)<<endl;
	double time = (end-start);
	n_parallel = n_iter;
	time_parallel = time;
}


void serial(double* temp)
{
	double *temp_old= new double[LX];
 	double alpha=(density*specific_heat)/(thermal_cond);
	double error=0;
	double num_sum=0.0 ,deno_sum=0.0;
 	//test for array population 

 	// for(int i=0; i<LX; i++)
 	// {
 	// 	temp[i]=i;
 	// }
 	// for(int i=0; i<LX; i++)
 	// {
 	// 	std::cout<<temp[i]<<"\n"<<std::endl;
 	// }

	//initialization

 	for(int i=0; i<LX; i++)
 	{
 		temp[i]=0;
		temp_old[i]=0;
 	}

 	temp[0]=100.0;
 	temp[LX-1]=1000.0;
	temp_old[0]=100.0;
 	temp_old[LX-1]=1000.0;
	clock_t start, end;


 	// for(int n_iter=0; n_iter<10000; n_iter++)
	// {
	// 	for(int i=1; i<LX-1; i++)
	// 	{
	// 		temp[i]=0.5*(temp[i+1]+temp[i-1]);
	// 	}

	// }
	int n_iter=0;
 	start = clock();
	for(n_iter=0; n_iter<1000000; n_iter++)
	{
		num_sum = 0.0;
		deno_sum = 0.0;
		for(int i=1; i<LX-1; i++)
		{
			temp[i]=0.5*(temp[i+1]+temp[i-1]);
		}

		for(int j=1; j<LX-1; j++)
		{ 
			num_sum += pow((temp[j]-temp_old[j]),2);
			deno_sum += pow((temp_old[j]),2);
		}
		error = sqrt(num_sum/deno_sum);
		if (error<1e-12)
		break;
		
		for(int k=1; k<LX-1; k++)
		{
			temp_old[k]=temp[k];
		}
	}
	end = clock();

	// cout<<"\n"<<"serial iteration no. ="<<" "<<n_iter<<endl;
	// cout<<endl<<"Serial Run time = "<<(end-start)/double(CLOCKS_PER_SEC)<<endl;
	double time = (end-start)/double(CLOCKS_PER_SEC);
	n_serial = n_iter;
	time_serial = time;
}


int main()
{
	double* temp_serial = new double[LX];
	double* temp_parallel = new double[LX];
	double diff[LX];

	serial(temp_serial);
	parallel(temp_parallel);

	for(int i=0; i<LX; i++)
	{
		diff[i]= (fabs(temp_parallel[i] - temp_serial[i]) / temp_serial[i]) * 100;
	}

	ofstream f;
	f.open("output.dat", ios::out);
		for(int i=0; i<LX; i++)
		{
			f<<i<<" "<<(temp_parallel[i])<<" "<<(temp_serial[i])<<" "<<diff[i]<<endl;
		}
	f.close();	

	// for(int i=0; i<LX; i++)
	// {
	// // cout<<(temp_parallel[i])<<" | "<<(temp_serial[i])<<" | "<<"difference in temp="<<" "<<diff[i]<<endl; //comment for gnuplot style data
	// cout<<i<<" "<<(temp_parallel[i])<<" "<<(temp_serial[i])<<" "<<"difference in temp="<<" "<<diff[i]<<endl;
	// }

	//comment all below for gnuplot style data
	// cout<<"\n"<<endl;
	// cout<<"Parallel Iteration no. ="<<" "<<n_parallel<<"\n"<<"Parallel runtime="<<time_parallel<<endl;
	// cout<<"\n"<<endl;
	// cout<<"Serial Iteration no. ="<<" "<<n_serial<<"\n"<<"Serial runtime="<<time_serial<<endl;

}

