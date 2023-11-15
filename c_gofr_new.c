#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

void free_memory(double *tab){
	free(tab);
}


void compute_gofr(int *gofr, double *MySnapshot, int ntypes, int *Types, double discrete, double *acell, double *rvec, int ndivx, double maxlength, int diag){
	double dist, deltX, deltY, deltZ, Xne, Yne, Zne, Xat, Yat, Zat, rda, rdb, rdc;
	int tmin, tmax, adj_tmin, adj_tmax, col;

	for(int el=0 ; el<ntypes ; el++){//We go element by element 
		
		tmin = Types[el];//And determine the first and last index in the atom list corresponting to this element
		tmax = Types[el+1];

		for(int at = tmin ; at<tmax ; at++){//parsing the atoms list for central atoms			
	
			Xat = MySnapshot[3*at];
			Yat = MySnapshot[3*at+1];
			Zat = MySnapshot[3*at+2];
			for(int adj_el = el ; adj_el<ntypes ; adj_el++){//We go element by element for their neighbors				

				//Definition of the first and last indexes in the atom list corresponding to this element

				if(adj_el==el){//If the 2 elements are the same, we adjust the first neighbor atom's index
					if(at==tmax-1){
						adj_el++; 
						adj_tmin=Types[adj_el];
						adj_tmax=Types[adj_el+1];
					}
					else{
						adj_tmin = at+1;
						adj_tmax = Types[adj_el+1];
					}
				}
				else{
					adj_tmin = Types[adj_el];
					adj_tmax = Types[adj_el+1];
				}								

				for(int ne=adj_tmin ; ne<adj_tmax ; ne++){//And for their neighbors
		
					Xne = MySnapshot[3*ne];
					Yne = MySnapshot[3*ne+1];
					Zne = MySnapshot[3*ne+2];

					deltX = Xat-Xne;
					deltY = Yat-Yne;
					deltZ = Zat-Zne;
					//Definition of the distance along x, y and z for the 2 atoms
					if(diag){//Orthogonal cell case
						if(fabs(deltX/acell[0])>0.5){
							deltX = acell[0]-fabs(deltX);
						}
						if(fabs(deltY/acell[1])>0.5){
							deltY = acell[1]-fabs(deltY);
						}
						if(fabs(deltZ/acell[2])>0.5){
							deltZ = acell[2]-fabs(deltZ);
						}
					}
					else{//Not orthogonal case
					    	rda = deltX * rvec[0] + deltY * rvec[1] + deltZ * rvec[2];//We project the position difference in the 3 directions of the cell thanks to reciprocal vectors (remark : 1 > (rda, rdb, rdc) > 0)
			    			rdb = deltX * rvec[3] + deltY * rvec[4] + deltZ * rvec[5];
			    			rdc = deltX * rvec[6] + deltY * rvec[7] + deltZ * rvec[8];
				
					    	//We check whether finding the smallest distance between at and ne necessitates to fold the cell, for each 3 directions
					    	if(fabs(rda)>0.5){
							rda = -1*(rda/fabs(rda)*(1-fabs(rda)));//Folding
					    	}
					    	if(fabs(rda)>0.5){
							rdb = -1*(rdb/fabs(rdb)*(1-fabs(rdb)));
					    	}
					    	if(fabs(rdc)>0.5){
							rdc = -1*(rdc/fabs(rdc)*(1-fabs(rdc)));
					    	}
					
					    	deltX = rda * acell[0] + rdb * acell[3] + rdc * acell[6];//We rebuild the position difference in the cartesian space
					    	deltY = rda * acell[1] + rdb * acell[4] + rdc * acell[7];
					    	deltZ = rda * acell[2] + rdb * acell[5] + rdc * acell[8];
						//printf("%f %f %f",deltX,deltY,deltZ);
					}
					dist = sqrt(pow(deltX,2)+pow(deltY,2)+pow(deltZ,2));
					if(dist<=maxlength/2){//If the distance is in the volume we take in account, we incremente the gofr in the right bin
						col = (int)(dist/(discrete));
						if(el == adj_el){
							gofr[(ntypes*ntypes)*col+(el*ntypes + adj_el)] = gofr[(ntypes*ntypes)*col+(el*ntypes + adj_el)]+2;
						}
						else{
							gofr[(ntypes*ntypes)*col+(el*ntypes + adj_el)]++;
						}
					}
				}	
			}
		}	
	}
}	