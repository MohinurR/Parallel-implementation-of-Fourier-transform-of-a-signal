/*
Raupova Mokhinur
Moscow State University in Tashkent
FFT PARALLEL IMPLEMETATION
*/
/* 
	Эта программа вычисляет БПФ параллельно с использованием mpi.
	Запускает программу X раз и вычисляет время для каждой итерации. Затем усредняет затраченное время
	Каждый процесс получает подтаблицу из таблицы. Затем вычислите четное и нечетное для каждого коэффициента X, тогда P0 соберет эти значения и распечатает их в файл
*/
#include <stdio.h>
#include <mpi.h> //Чтобы использовать MPI
#include <complex.h> //Чтобы использовать комплексные числа
#include <math.h>	//для cos() и sin()
#include "timer.h" //для использованию таймера

#define PI 3.14159265
#define bigN 16384 //Размер задачи
#define howmanytimesavg 3 //Сколько раз я хочу бегать в AVG?

int main()
{
	int my_rank,comm_sz;
	MPI_Init(NULL,NULL); //запустить MPI
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);   ///сколько процессов мы используем?
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);  
	double start,finish;
	double avgtime = 0;
	FILE *outfile;
	int h;
	if(my_rank == 0) //если процесс 0 открыть файл
	{
		outfile = fopen("ParallelVersionOutput.txt", "w"); //открыть из текущего каталога
	}
	for(h = 0; h < howmanytimesavg; h++) //цикл для выполнения несколько раз за время AVG.
	{
		if(my_rank == 0) //Если 0, запускается таймер
		{	
			start = MPI_Wtime();
		}
		int i,k,n,j; //Основные переменные цикла

		double complex evenpart[(bigN / comm_sz / 2)]; //массив для сохранения данных для EVENHALF
		double complex oddpart[(bigN / comm_sz / 2)]; //массив для сохранения данных для ODDHALF
		double complex evenpartmaster[ (bigN / comm_sz / 2) * comm_sz]; //массив для сохранения данных для EVENHALF
		double complex oddpartmaster[ (bigN / comm_sz / 2) * comm_sz]; //массив для сохранения данных для ODDHALF
		double storeKsumreal[bigN]; //store the K real variable so we can abuse symmerty
		double storeKsumimag[bigN]; //store the K imaginary variable so we can abuse symmerty
		
		double subtable[(bigN / comm_sz)][3]; //Each process owns a subtable from the table below 
		
		double table[bigN][3] = //TABLE of numbers to use
							{
							 0,3.6,2.6, //n, Real,Imaginary CREATES TABLE
							 1,2.9,6.3,
							 2,5.6,4.0,
							 3,4.8,9.1,
							 4,3.3,0.4,
							 5,5.9,4.8,
							 6,5.0,2.6,
							 7,4.3,4.1,
							 };
			if(bigN > 8)  //Everything after row 8 is all 0's
			{
				for(i = 8; i < bigN; i++)
				{
					table[i][0] = i;
					for(j = 1; j < 3;j++)
					{
						table[i][j] = 0.0; //set to 0.0
					}
				}
			}
		int sendandrecvct = (bigN / comm_sz) * 3; //how much to send and recieve??
		MPI_Scatter(table,sendandrecvct,MPI_DOUBLE,subtable,sendandrecvct,MPI_DOUBLE,0,MPI_COMM_WORLD); //scatter the table to subtables
		for (k = 0; k < bigN / 2; k++) //K coeffiencet Loop 
		{
					/* Variables used for the computation */
			double sumrealeven = 0.0; //sum of real numbers for even
			double sumimageven = 0.0; //sum of imaginary numbers for even
			double sumrealodd = 0.0; //sum of real numbers for odd
			double sumimagodd = 0.0; //sum of imaginary numbers for odd
			
			for(i = 0; i < (bigN/comm_sz)/2; i++) //Sigma loop EVEN and ODD
			{
				double factoreven , factorodd = 0.0;
				int shiftevenonnonzeroP = my_rank * subtable[2*i][0]; //used to shift index numbers for correct results for EVEN.
				int shiftoddonnonzeroP = my_rank * subtable[2*i + 1][0]; //used to shift index numbers for correct results for ODD.
				
				/* -------- EVEN PART -------- */
				double realeven = subtable[2*i][1]; //Access table for real number at spot 2i
				double complex imaginaryeven = subtable[2*i][2]; //Access table for imaginary number at spot 2i
				double complex componeeven = (realeven + imaginaryeven * I); //Create the first component from table
				if(my_rank == 0) //if proc 0, dont use shiftevenonnonzeroP
				{
					factoreven = ((2*PI)*((2*i)*k))/bigN; //Calculates the even factor for Cos() and Sin()										
								//   *********Reduces computational time*********
				}
				else //use shiftevenonnonzeroP
				{
					factoreven = ((2*PI)*((shiftevenonnonzeroP)*k))/bigN; //Calculates the even factor for Cos() and Sin()										
								//   *********Reduces computational time*********
				}
				double complex comptwoeven = (cos(factoreven) - (sin(factoreven)*I)); //Create the second component
				
				evenpart[i] = (componeeven * comptwoeven); //store in the evenpart array
				
				/* -------- ODD PART -------- */
				double realodd = subtable[2*i + 1][1]; //Access table for real number at spot 2i+1
				double complex imaginaryodd = subtable[2*i + 1][2]; //Access table for imaginary number at spot 2i+1
				double complex componeodd = (realodd + imaginaryodd * I); //Create the first component from table
				if (my_rank == 0)//if proc 0, dont use shiftoddonnonzeroP
				{
					factorodd = ((2*PI)*((2*i+1)*k))/bigN;//Calculates the odd factor for Cos() and Sin()										
								// *********Reduces computational time*********	
				}
				else //use shiftoddonnonzeroP
				{
					factorodd = ((2*PI)*((shiftoddonnonzeroP)*k))/bigN;//Calculates the odd factor for Cos() and Sin()										
								// *********Reduces computational time*********
				}
							
				double complex comptwoodd = (cos(factorodd) - (sin(factorodd)*I));//Create the second component

				oddpart[i] = (componeodd * comptwoodd); //store in the oddpart array
				
			}
			/*Process ZERO gathers the even and odd part arrays and creates a evenpartmaster and oddpartmaster array*/
			MPI_Gather(evenpart,(bigN / comm_sz / 2),MPI_DOUBLE_COMPLEX,evenpartmaster,(bigN / comm_sz / 2), MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
			MPI_Gather(oddpart,(bigN / comm_sz / 2),MPI_DOUBLE_COMPLEX,oddpartmaster,(bigN / comm_sz / 2), MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);

			if(my_rank == 0)
			{
				for(i = 0; i < (bigN / comm_sz / 2) * comm_sz; i++) //loop to sum the EVEN and ODD parts
				{
					sumrealeven += creal(evenpartmaster[i]); //sums the realpart of the even half
					sumimageven += cimag(evenpartmaster[i]); //sums the imaginarypart of the even half
					sumrealodd += creal(oddpartmaster[i]); //sums the realpart of the odd half
					sumimagodd += cimag(oddpartmaster[i]); //sums the imaginary part of the odd half
				}
				storeKsumreal[k] = sumrealeven + sumrealodd; //add the calculated reals from even and odd
				storeKsumimag[k]  = sumimageven + sumimagodd; //add the calculated imaginary from even and odd
				storeKsumreal[k + bigN/2] = sumrealeven - sumrealodd; //ABUSE symmetry Xkreal + N/2 = Evenk - OddK
				storeKsumimag[k + bigN/2] = sumimageven - sumimagodd; //ABUSE symmetry Xkimag + N/2 = Evenk - OddK
				if(k <= 10) //Do the first 10 K's
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
			GET_TIME(finish); //stop timer
			double timeElapsed = finish-start; //Time for that iteration
			avgtime = avgtime + timeElapsed; //AVG the time 
			fprintf(outfile,"Time Elaspsed on Iteration %d: %f Seconds\n", (h+1),timeElapsed);
		}
	}
	if(my_rank == 0)
	{
		avgtime = avgtime / howmanytimesavg; //get avg time
		fprintf(outfile,"\nAverage Time Elaspsed: %f Seconds", avgtime);
		fclose(outfile); //CLOSE file ONLY proc 0 can.
	}
	MPI_Barrier(MPI_COMM_WORLD); //wait to all proccesses to catch up before finalize
	MPI_Finalize(); //End MPI
	return 0;
}
