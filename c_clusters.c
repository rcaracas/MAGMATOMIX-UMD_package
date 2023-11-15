#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>
#include <string.h>


void free_memory_int(int *tab){
	free(tab);
}


void free_memory_double(double *tab){
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

double min3(double a, double b, double c){//The minimum between 3 numbers.

	double minimum = a;

	if(minimum>b){minimum=b;}
	if(minimum>c){minimum=c;}

	return minimum;
}


double* angles(int *polyhedras, int nPol, int M, double *MySnapshot, int *CentIndexes, int nCent, int *OutIndexes, int nOut, double acell0, double acell1, double acell2){//Calculates the angle between the atoms of the type Out-Cent-Out in a snapshot
	double* AnglesList = NULL;//The list of numerical values of the angles
	int* OutList = NULL;//The list of outer type atoms
	double* CoordList = NULL;//The list which will (transitorily) contains the coordinates of the atoms linked to a given central atom
	int at,ne1,ne2;	
	double xCent,yCent,zCent,Deltx1,Delty1,Deltz1,Deltx2,Delty2,Deltz2,valx1,valy1,valz1,valx2,valy2,valz2,deltzz,deltxx,deltyy,valxx,valyy,valzz,dist1,dist2,distij,angle;
	int centat=1;
	int nbounds=0 ;
	int ncent=0;
	int indexPol=0;
	int flagcoord=0;

	CoordList = calloc(3*M,sizeof(double));
	OutList = calloc(M,sizeof(int));

	int numCent = num_atoms(CentIndexes,nCent);
	int numOut = num_atoms(OutIndexes,nOut);

	AnglesList = calloc(((M+1)*M/2*numCent+numCent),sizeof(double));//Maximum possible size

	if(CoordList == NULL || OutList == NULL || AnglesList == NULL){
		printf("Memory Allocation failure for CoordList (length %d), OutList (length %d) or AnglesList (length %d) in function 'angles'",3*M,M,(M+1)*M/2*numCent+numCent);
		return EXIT_FAILURE;
	}
	AnglesList[0]=1;//Index of the next angle to write in the tab
	while(ncent<=nPol){//We parse the polyhedras central atom by central atom

		indexPol++;
		at = polyhedras[indexPol];
		if(at==-1){//If we're at the end of the cluster, we mark it by writing a -1 and we update the control variables
			AnglesList[(int)AnglesList[0]]=-1;
			AnglesList[0]++;
			centat=1;
			ncent++;
			nbounds = 0;
		}		
		
		else if(centat==1){	//If the atom we're looking at is of the central type
			centat=0;
			//We extract the coordinates of said central atom
			xCent = MySnapshot[3*at];
			yCent = MySnapshot[3*at+1];
			zCent = MySnapshot[3*at+2];

		}
		else{	//If not, at is a coordinating atom
			flagcoord=1;
			while(flagcoord){	//we parse polyhedras to extract all the coordinating atoms of the cluster, count them and put them in OutList
				at = polyhedras[indexPol];				
				OutList[nbounds]=at;
				CoordList[3*nbounds]=MySnapshot[3*at];
				CoordList[3*nbounds+1]=MySnapshot[3*at+1];
				CoordList[3*nbounds+2]=MySnapshot[3*at+2];
				nbounds++;
				indexPol++;				
				if(polyhedras[indexPol]==-1){
					flagcoord=0;
					indexPol--;
				}

			}
			
			for(int i=0 ; i<nbounds ; i++){//Calculation of the angles ; double looping on the OutList
				
				Deltx1 = CoordList[3*i]-xCent;
				Delty1 = CoordList[3*i+1]-yCent;
				Deltz1 = CoordList[3*i+2]-zCent;

				for(int j=i+1 ; j<nbounds ; j++){
					

					Deltx2 = CoordList[3*j]-xCent;
					Delty2 = CoordList[3*j+1]-yCent;
					Deltz2 = CoordList[3*j+2]-zCent;

					deltxx = CoordList[3*i]-CoordList[3*j]; 					
					deltyy = CoordList[3*i+1]-CoordList[3*j+1]; 					
					deltzz = CoordList[3*i+2]-CoordList[3*j+2];					

					valx1=min3(pow(Deltx1,2),pow(Deltx1-acell0,2),pow(Deltx1+acell0,2));													
					valy1=min3(pow(Delty1,2),pow(Delty1-acell1,2),pow(Delty1+acell1,2));
					valz1=min3(pow(Deltz1,2),pow(Deltz1-acell2,2),pow(Deltz1+acell2,2));

					dist1=valx1+valy1+valz1;

					valx2=min3(pow(Deltx2,2),pow(Deltx2-acell0,2),pow(Deltx2+acell0,2));
					valy2=min3(pow(Delty2,2),pow(Delty2-acell1,2),pow(Delty2+acell1,2));
					valz2=min3(pow(Deltz2,2),pow(Deltz2-acell2,2),pow(Deltz2+acell2,2));

					dist2=valx2+valy2+valz2;
				
					valxx = min3(pow(deltxx,2),pow(deltxx-acell0,2),pow(deltxx+acell0,2));
					valyy = min3(pow(deltyy,2),pow(deltyy-acell1,2),pow(deltyy+acell1,2));
					valzz = min3(pow(deltzz,2),pow(deltzz-acell2,2),pow(deltzz+acell2,2));

					distij = valxx+valyy+valzz;
					AnglesList[(int)AnglesList[0]] = acos((dist1+dist2-distij)/(2*sqrt(dist1*dist2)))/3.141592653589793*180;	//We put the value in AnglesList
					AnglesList[0]++;
				}
			}
		}
	}

	free(OutList);
	free(CoordList);
	return AnglesList;
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
	if(r>0){	//The recursion stops when r becomes 0
		for(int i = indBonds[atom]; i<indBonds[atom+1]; i++){	//We loop over the coordinated atoms
			at = SnapshotBonds[i];
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
	int imax;
	int prev_index=0;

	int numCent = num_atoms(CentIndexes,nCent);
	int numOut = num_atoms(OutIndexes,nCent);

	int meanClustSize = indBonds[nAts]/(nAts);

	if(r>0){
		maxClustSize = min(pow(M,r)+1,numOut+2);	//The maximum size of a cluster, used as an estimation unit to increase the size of neighbors
		len = indBonds[nAts]+2*nAts+1 ;
		CentAts = calloc(natom,sizeof(int));
		OutAts = calloc(natom,sizeof(int));
	

		if (CentAts == NULL||OutAts == NULL)
		{
		    printf("Memory Allocation failure for CentAts (length %d) or Outats (length %d)\n,",natom,natom);
		    return EXIT_FAILURE;
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
		    return EXIT_FAILURE;
		}
		for(int i=0 ; i<natom ; i++){
			AllAts[i]=-1;
		}
		reinit(AllAts,CentIndexes,nCent);
		reinit(AllAts,OutIndexes,nOut);
	}

	neighbors = calloc(len,sizeof(int));
	if (neighbors == NULL)
	{
	    printf("Memory Allocation failure for neighbors : len is %d\n",len);
	    return EXIT_FAILURE;
	}


//Initialisation

	for(int i=0 ; i<len ; i++){
		neighbors[i]=-1;
	}
	neighbors[0]=1;



//neighbors computation

	if(r>0){
		while(index<natom){
			if(CentAts[index] != -1){	//Looping central atom by central atom
				atom=CentAts[index];
				CentAts[index]=-1;
				neighbors[neighbors[0]]=atom;
				neighbors[0]++;
				prev_index = neighbors[0];
				clusteringrec(neighbors,SnapshotBonds,indBonds,OutAts,OutIndexes,nOut,atom,r);
				if(neighbors[0]==prev_index){	//We separate the clusters by the number -1
					neighbors[0]--;
					neighbors[neighbors[0]]=-1;
				}
				else{
					neighbors[0]++;
				}

			}		


			index++;
			reinit(OutAts,OutIndexes,nOut);
			if(neighbors[0]>len-2*maxClustSize){	//We reallocate if the size of neighbors is too small
				
				int *neighborstemp = (int *)realloc(neighbors,(len+2*maxClustSize+1)*sizeof(int));

				if(neighborstemp == NULL){printf("Memory allocation failure for neighbors \n");return EXIT_FAILURE;}

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
			neighbors[0]--;
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
	