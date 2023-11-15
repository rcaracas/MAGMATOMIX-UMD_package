#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


//These scripts are meant to calculate the msd in slightly different ways. Their functionment is very similar.

void free_memory(double *array){
	free(array);
}

double* compute_msd(const double *pos, const int hh, const int vv, const int ballistic, const int nitmax)//Calculates the msd from the positions pos.
{

	double* msd;
	msd = calloc((nitmax/vv),sizeof(double));

        for (int i=ballistic+hh ; i<nitmax ; i=i+hh){
		//We extract the position of the atom at time i
            const double Xorigin = pos[3*i];
            const double Yorigin = pos[3*i+1];
            const double Zorigin = pos[3*i+2];

            for (int j=0 ; j<nitmax ; j=j+vv){
			//And at time j
                const double Xfuture = pos[3*(i+j)];
                const double Yfuture = pos[3*(i+j)+1];
                const double Zfuture = pos[3*(i+j)+2];

                msd[(j/vv)] = msd[(j/vv)] + pow((Xfuture-Xorigin),2) + pow((Yfuture-Yorigin),2) + pow((Zfuture-Zorigin),2);//MSD calculation

            }
        }
	return msd;
}

double* compute_msd_tilted(const double *pos, const int hh, const int vv, const int ballistic, const int nitmax, const double *Axes)//Computes the msd along particular axes
{
	double* msd;
	msd = calloc(3*(nitmax/vv),sizeof(double));

	double Xorigin, Yorigin, Zorigin, Xfuture, Yfuture, Zfuture, deltX, deltY, deltZ, deltA1, deltA2, deltA3;

        for (int i=ballistic+hh ; i<nitmax ; i=i+hh){
			
        	Xorigin = pos[3*i];
        	Yorigin = pos[3*i+1];
        	Zorigin = pos[3*i+2];

            	for (int j=0 ; j<nitmax ; j=j+vv){

                	Xfuture = pos[3*(i+j)];
                	Yfuture = pos[3*(i+j)+1];
                	Zfuture = pos[3*(i+j)+2];

			deltX = Xfuture-Xorigin;
			deltY = Yfuture-Yorigin;
			deltZ = Zfuture-Zorigin;
	
			deltA1 = deltX*Axes[0]+deltY*Axes[1]+deltZ*Axes[2];//Projection of the position difference vector on the choosen axes
			deltA2 = deltX*Axes[3]+deltY*Axes[4]+deltZ*Axes[5];
			deltA3 = deltX*Axes[6]+deltY*Axes[7]+deltZ*Axes[8];

                	msd[(j/vv)] = msd[(j/vv)] + pow(deltA1,2);
                	msd[(j/vv)+(nitmax/vv)] = msd[(j/vv)+(nitmax/vv)] + pow(deltA2,2);
                	msd[(j/vv)+2*(nitmax/vv)] = msd[(j/vv)+2*(nitmax/vv)] + pow(deltA3,2);
            }
        }
return msd;

}

double* compute_msd_tilted_multi(const double *pos, const int hh, const int vv, const int ballistic, const int nitmax, const double *Axes, const int nAxes)//Computes the msd along an arbitraty number of particular axes
{
	double* msd;
	double deltAxe;

	msd = calloc(nAxes*(nitmax/vv),sizeof(double));

	double Xorigin, Yorigin, Zorigin, Xfuture, Yfuture, Zfuture, deltX, deltY, deltZ, deltA1, deltA2, deltA3;

        for (int i=ballistic+hh ; i<nitmax ; i=i+hh){
			
        	Xorigin = pos[3*i];
        	Yorigin = pos[3*i+1];
        	Zorigin = pos[3*i+2];

            	for (int j=0 ; j<nitmax ; j=j+vv){

                	Xfuture = pos[3*(i+j)];
                	Yfuture = pos[3*(i+j)+1];
                	Zfuture = pos[3*(i+j)+2];

			deltX = Xfuture-Xorigin;
			deltY = Yfuture-Yorigin;
			deltZ = Zfuture-Zorigin;

			for(int k=0 ; k< nAxes ; k++){
				deltAxe=deltX*Axes[3*k]+deltY*Axes[3*k+1]+deltZ*Axes[3*k+2];//Projection of the position difference vector on the choosen axes
				msd[(j/vv) +k*nitmax/vv] = msd[(j/vv) +k*nitmax/vv] + pow(deltAxe,2);
			}	
            }
        }
return msd;
}