/*
Raupova Mokhinur
*/

#include <stdio.h>
#include <mpi.h> 
#include <complex.h> 
#include <math.h>
#include "timer.h" 

#define PI 3.14159265
#define bigN 16384
#define howmanytimesavg 3

int main()
{
	int my_rank,comm_sz;
	MPI_Init(NULL,NULL); //запустить MPI
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);   
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);  
	double start,finish;
	double avgtime = 0;
	FILE *outfile;
	int h;
	if(my_rank == 0) 
	{
		outfile = fopen("ParallelVersionOutput.txt", "w"); 
	}
	for(h = 0; h < howmanytimesavg; h++) 
	{
		if(my_rank == 0)
		{	
			start = MPI_Wtime();
		}
		int i,k,n,j; 

		double complex evenpart[(bigN / comm_sz / 2)]; 
		double complex oddpart[(bigN / comm_sz / 2)]; 
		double complex evenpartmaster[ (bigN / comm_sz / 2) * comm_sz]; 
		double complex oddpartmaster[ (bigN / comm_sz / 2) * comm_sz]; 
		double storeKsumreal[bigN]; 
		double storeKsumimag[bigN]; 
		
		double subtable[(bigN / comm_sz)][3]; 
		
		double table[bigN][3] = 
							{
							 0,3.6,2.6,
							 1,2.9,6.3,
							 2,5.6,4.0,
							 3,4.8,9.1,
							 4,3.3,0.4,
							 5,5.9,4.8,
							 6,5.0,2.6,
							 7,4.3,4.1,
							 };
			if(bigN > 8)  
			{
				for(i = 8; i < bigN; i++)
				{
					table[i][0] = i;
					for(j = 1; j < 3;j++)
					{
						table[i][j] = 0.0; 
					}
				}
			}
		int sendandrecvct = (bigN / comm_sz) * 3; 
		MPI_Scatter(table,sendandrecvct,MPI_DOUBLE,subtable,sendandrecvct,MPI_DOUBLE,0,MPI_COMM_WORLD); 
		for (k = 0; k < bigN / 2; k++) 
		{
					
			double sumrealeven = 0.0; 
			double sumimageven = 0.0;
			double sumrealodd = 0.0; 
			double sumimagodd = 0.0; 
			
			for(i = 0; i < (bigN/comm_sz)/2; i++) 
			{
				double factoreven , factorodd = 0.0;
				int shiftevenonnonzeroP = my_rank * subtable[2*i][0];
				int shiftoddonnonzeroP = my_rank * subtable[2*i + 1][0]; 
								
				double realeven = subtable[2*i][1]; 
				double complex imaginaryeven = subtable[2*i][2];
				double complex componeeven = (realeven + imaginaryeven * I); 
				if(my_rank == 0) 
				{
					factoreven = ((2*PI)*((2*i)*k))/bigN; 				
				}
				else
				{
					factoreven = ((2*PI)*((shiftevenonnonzeroP)*k))/bigN; 
				}
				double complex comptwoeven = (cos(factoreven) - (sin(factoreven)*I)); 
				
				evenpart[i] = (componeeven * comptwoeven); 
				
				double realodd = subtable[2*i + 1][1]; 
				double complex imaginaryodd = subtable[2*i + 1][2]; 
				double complex componeodd = (realodd + imaginaryodd * I); 
				if (my_rank == 0)
				{
					factorodd = ((2*PI)*((2*i+1)*k))/bigN;
				}
				else
				{
					factorodd = ((2*PI)*((shiftoddonnonzeroP)*k))/bigN;
				}
							
				double complex comptwoodd = (cos(factorodd) - (sin(factorodd)*I));

				oddpart[i] = (componeodd * comptwoodd);
				
			}
			MPI_Gather(evenpart,(bigN / comm_sz / 2),MPI_DOUBLE_COMPLEX,evenpartmaster,(bigN / comm_sz / 2), MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
			MPI_Gather(oddpart,(bigN / comm_sz / 2),MPI_DOUBLE_COMPLEX,oddpartmaster,(bigN / comm_sz / 2), MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);

			if(my_rank == 0)
			{
				for(i = 0; i < (bigN / comm_sz / 2) * comm_sz; i++)
				{
					sumrealeven += creal(evenpartmaster[i]);
					sumimageven += cimag(evenpartmaster[i]); 
					sumrealodd += creal(oddpartmaster[i]); 
					sumimagodd += cimag(oddpartmaster[i]); 
				}
				storeKsumreal[k] = sumrealeven + sumrealodd; 
				storeKsumimag[k]  = sumimageven + sumimagodd;
				storeKsumreal[k + bigN/2] = sumrealeven - sumrealodd; 
				storeKsumimag[k + bigN/2] = sumimageven - sumimagodd; 
				if(k <= 10) 
				{
					if(k == 0)
					{
						fprintf(outfile," \n\n TOTAL PROCESSED SAMPLES : %d\n",bigN);
					}
					fprintf(outfile,"================================\n");
					fprintf(outfile,"XR[%d]: %.4f XI[%d]: %.4f \n",k,storeKsumreal[k],k,storeKsumimag[k]);
					fprintf(outfile,"================================\n");
				}
			}
		}
		if(my_rank == 0)
		{
			GET_TIME(finish);
			double timeElapsed = finish-start; 
			avgtime = avgtime + timeElapsed;
			fprintf(outfile,"Time Elaspsed on Iteration %d: %f Seconds\n", (h+1),timeElapsed);
		}
	}
	if(my_rank == 0)
	{
		avgtime = avgtime / howmanytimesavg; 
		fprintf(outfile,"\nAverage Time Elaspsed: %f Seconds", avgtime);
		fclose(outfile); 
	}
	MPI_Barrier(MPI_COMM_WORLD); 
	MPI_Finalize(); 
	return 0;
}
