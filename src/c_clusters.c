#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>
#include <string.h>


void free_memory_int(int *tab){
	free(tab);
}



int min(int a,int b){//returns the minimum of 2 numbers, except if one of them is negative (it will return the other).

	if(a>b||a<0){return b;}
	else{return a;}

}	


void reinit(int *Tab,int *Indexes,int n){//(re)initializes CentAts or OutAts with all atoms presents
	int imax=0;
	int imin=0;
	for(int el = 0 ; el<n ; el++){
		imin = Indexes[2*el];
		imax = Indexes[2*el+1];
		for(int i = imin ; i<imax+1 ; i++){
			Tab[i] = i;
		}	
	}
}

bool is_in(int atom,int *Indexes, int n){//"Is atom in the ensemble designed by the tab Indexes, who has a length of 2n ?"
	bool result = false;
	for(int i=0 ; i<n ; i++){
		if((atom>=Indexes[2*i] && atom<=Indexes[2*i+1])){
			result = true;
		}
	}
	return result;
}

int num_atoms(int* Indexes, int n){//"How many atoms are designed by the tab Indexes ?"
	int num = 0;
	for(int i=0 ; i<n ; i++){
		num = num + Indexes[2*i+1]-Indexes[2*i] + 1;
	}
	return num;
}

void clusteringall(int *neighbors,const int *SnapshotBonds, const int *indBonds,int *AllAtoms, const int atom){		//Clustering as an infinite polymerization when r=0

	int at;

	for(int i = indBonds[atom]; i<indBonds[atom+1]; i++){	//Indexes of all the atoms bound to atom
		at = SnapshotBonds[i];				//We extract these atoms one by one
		if(AllAtoms[at]!=-1){				//If it's in the list
			AllAtoms[at]=-1;			//We retire the atom from the AllAtoms list
			neighbors[neighbors[0]]=at;		//We add the atom in the cluster
			neighbors[0]++;
			clusteringall(neighbors,SnapshotBonds,indBonds,AllAtoms,at);	//We loop again until all the connected atoms are exhausted
		}
	}
}

void clusteringrec(int *neighbors,const int *SnapshotBonds, const int *indBonds, int *OutAts, int *OutIndexes, const int nOut, const int atom, const int r){	//

	int at=0;
	//printf("in clusteringrec, with rflag = %d that has indBonds[atom] %d\n",r,indBonds[atom]);
	if(r>0){	//The recursion stops when r becomes 0
		for(int i = indBonds[atom]; i<indBonds[atom+1]; i++){	//We loop over the coordinated atoms
			at = SnapshotBonds[i];
			//printf(" next bonded atom is: %d\n",at);
			if(is_in(at,OutIndexes,nOut)){	//If, of course, the coordinated atom is in the list of the outer atoms
				if(OutAts[at]!=-1){	//If at has not been count yet
					OutAts[at]=-1;	//We retire at from the OutAts
					neighbors[neighbors[0]]=at;
					neighbors[0]++;
				}
				clusteringrec(neighbors,SnapshotBonds,indBonds,OutAts,OutIndexes,nOut,at,r-1);
			}
		}
	}
}

int* fullclustering(const int *SnapshotBonds, const int *indBonds, const int natom, const int nAts, int *CentIndexes, int *OutIndexes, int nCent, int nOut, const int M, const int r){	//Uses the Bonds information to create the species
	
	int *AllAts = NULL;	//All the concerned atoms 	(used only if n=0)
	int *neighbors = NULL;	//Will contain the clusters
	int *CentAts = NULL;	//All the central atoms		(used only if n>0)
	int *OutAts = NULL;	//All the outer atoms		(used only if n>0)

	int maxClustSize = 0;
	int index = 0;
	int atom;
	int len;
	int prev_index=0;

	int numOut = num_atoms(OutIndexes,nOut);

    //printf("Starting the c function fullclustering \n");

	if(r>0){
		maxClustSize = min(pow(M,r)+1,numOut+2);	//The maximum size of a cluster, used as an estimation unit to increase the size of neighbors
		//printf("maxClustSize is %d\n",maxClustSize);
		len = indBonds[nAts]+2*nAts+1 ;
		CentAts = calloc(natom,sizeof(int));
		OutAts = calloc(natom,sizeof(int));
	

		if (CentAts == NULL||OutAts == NULL)
		{
		    printf("Memory Allocation failure for CentAts (length %d) or OutAts (length %d)\n,",natom,natom);
		    //return EXIT_FAILURE;
		}

		for(int i=0 ; i<natom ; i++){
			CentAts[i]=-1;
			OutAts[i]=-1;
		}
		reinit(CentAts,CentIndexes,nCent);
		reinit(OutAts,OutIndexes,nOut);
	}	
	else{
		len = 2*nAts+2;					
		
		AllAts = calloc(natom,sizeof(int));
		if (AllAts == NULL)
		{
		    printf("Memory Allocation failure for AllAts : len is %d\n",len);
		    //return EXIT_FAILURE;
		}
		for(int i=0 ; i<natom ; i++){
			AllAts[i]=-1;
		}
		reinit(AllAts,CentIndexes,nCent);
		reinit(AllAts,OutIndexes,nOut);
	}

    //printf("finished allocating memory for various strings\n");
	neighbors = calloc(len,sizeof(int));
	if (neighbors == NULL)
	{
	    printf("Memory Allocation failure for neighbors : len is %d\n",len);
	    //return EXIT_FAILURE;
	}

    //printf("Check in fullclustering: neighboring size is %ld and they are:\n", sizeof(neighbors));
    //for (int i = 0; i < sizeof(neighbors); i++) {
    //    printf("%d", neighbors[i]); // Print each integer
    //    if (i < sizeof(neighbors) - 1) {
    //        printf(", "); // Add a comma and space between numbers
    //    }
    //}


//Initialisation

	for(int i=0 ; i<len ; i++){
		neighbors[i]=-1;
	}
	neighbors[0]=1;
    //printf("new initialized neighbor list\n");
    //for (int i = 0; i < sizeof(neighbors); i++) {
    //    printf("%d", neighbors[i]); // Print each integer
    //    if (i < sizeof(neighbors) - 1) {
    //        printf(", "); // Add a comma and space between numbers
    //    }
    //}



//neighbors computation
    //printf("\nWe are starting to compute neighbors \n");
    //printf("There are natom %d\n",natom);
	if(r>0){
	    //printf("Will compute coordination spheres of order %d\n",r);
		while(index<natom){
    	    //printf("index, i.e. current atom, is %d\n",index);    //prints out the current atom, going through the entire list of atoms
    	    //printf("CentAts for current atom is %d\n",CentAts[index]);
		    //printf("Centats is %d\n",CentAts[index]);
			if(CentAts[index] != -1){	//Looping central atom by central atom
				atom=CentAts[index];
				//printf("Current central atom is %d\n",atom);        //prints only the central atom that counts
				CentAts[index]=-1;
				neighbors[neighbors[0]]=atom;
				neighbors[0]++;
				prev_index = neighbors[0];
				clusteringrec(neighbors,SnapshotBonds,indBonds,OutAts,OutIndexes,nOut,atom,r);
				//printf("First neighbor in the list is %d\n",neighbors[0]);
				if(neighbors[0]==prev_index){	//We separate the clusters by the number -1
					neighbors[0]--;
					neighbors[neighbors[0]]=-1;
				}
				else{
					neighbors[0]++;
				}
			}

            //printf("In the middle of the while loop. OutAts is \n");
            //for (int i = 0; i < sizeof(OutAts); i++) {
            //    printf("%d",OutAts[i]); // Print each integer
            //    if (i < sizeof(neighbors) - 1) {
            //        printf(", "); // Add a comma and space between numbers
            //    }
            //}
            //printf("\n and OutIndexes is \n");
            //for (int i = 0; i < sizeof(OutIndexes); i++) {
            //    printf("%d",OutIndexes[i]); // Print each integer
            //    if (i < sizeof(neighbors) - 1) {
            //        printf(", "); // Add a comma and space between numbers
            //    }
            //}
            //printf("\n");


			index++;
			reinit(OutAts,OutIndexes,nOut);
			if(neighbors[0]>len-2*maxClustSize){	//We reallocate if the size of neighbors is too small
				
				int *neighborstemp = (int *)realloc(neighbors,(len+2*maxClustSize+1)*sizeof(int));

				if(neighborstemp == NULL){printf("Memory allocation failure for neighbors \n");}
				// if(neighborstemp == NULL){printf("Memory allocation failure for neighbors \n");return EXIT_FAILURE;}

				neighbors = neighborstemp;
				len = len + 2*maxClustSize +1;
				for(int k=neighbors[0] ; k<len ; k++){
					neighbors[k] = -1;
				}													
			}
		}
		free(OutAts);
		free(CentAts);
		while(neighbors[neighbors[0]]==-1){
			neighbors[0]--;   //updates the size of actual neighbors
		}
		return neighbors;
	}
	else{	//If r=0
		while(index<natom){	//We parse the AllAts list
			if(AllAts[index] != -1){
				atom=AllAts[index];
				AllAts[index]=-1;
				neighbors[neighbors[0]]=atom;
				neighbors[0]++;
				clusteringall(neighbors,SnapshotBonds,indBonds,AllAts,atom);
				neighbors[0]++;
			}					
			index++;
		}
		while(neighbors[neighbors[0]]==-1){
			neighbors[0]--;
		}
		
		free(AllAts);
		return neighbors;
	}
}
	