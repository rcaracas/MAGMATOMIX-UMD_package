#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

void free_memory_int(int *tab){
	free(tab);
}

void free_memory(int **tab){
	free(tab);
}


double min(double a, double b, double c){
	double m=a;
	if(b<m){
		m=b;
	}
	if(c<m){
		m=c;
	}
return m;
}

int maxval(const int tab[], int L){//Maximal value of a tab

    int val;
    val=tab[0];
    for(int i=0 ; i<L ; i++){
        if(val<tab[i]){
            val=tab[i];
        }
    }
return val;
}

bool specificity(double x, double y, double z, const double *specifics){//Check if a particular point meets the spatial selectivity criteria

	if(x>=specifics[0] && x<=specifics[0]+specifics[3] && y>=specifics[1] && y<=specifics[1]+specifics[4] && z>=specifics[2] && z<=specifics[2]+specifics[5]){
		return true;
	}
	else{
		return false;
	}
}

int** compute_Bonds_full(const double *MySnapshot, const double *BondTable, const int *CrystalTypes, const int nAtoms, const int ntypat, const double *cell, const double *rvec, const int nCells, const double *specifics, const int specification, const int diag){
    int* DicoAtoms;
    int* nMesh;
    int** Mesh;
    int* MeshIndex;
    int** Bonding;
    int* nBonding;

    double x,y,z;
    int nX,nY,nZ,cellnum;
    double Lx, Ly, Lz;
    double rda,rdb,rdc;

    if(diag==1){            //If the simulation cell is orthogonal
	Lx = cell[0]/nCells;//we define the dimension of a sub-cell along x
    	Ly = cell[1]/nCells;			      //along y
   	Lz = cell[2]/nCells;			      //along z
    }
    
    int ntotAts=0;
    int numCoord=0;
    int indexCoord=0;
    int index = 0;
    int cx,cy,cz;    
    int at,ne;
    double dx,dy,dz,valx,valy,valz,distij;


    DicoAtoms = calloc((3*nAtoms+1),sizeof(int));	   //Will contain the indexes of the sub-cell in which each atom is
    nMesh = calloc((nCells*nCells*nCells),sizeof(int));    //Will contain the number of atoms in each sub-cell
    MeshIndex = calloc((nCells*nCells*nCells),sizeof(int));//Same thing than nMesh but will be useful during it's filling
    nBonding = calloc(nAtoms , sizeof(int));		   //Will contain the number of bonds for each atom

    for(int i=0 ; i<3*nAtoms+1 ; i++){
	DicoAtoms[i]=-1;
    }
    for (int i=0 ; i<nAtoms ; i=i+1){//This loop is used to fill nMesh and DicoAtoms with the appropriate values

        x = MySnapshot[3*i];
        y = MySnapshot[3*i+1];
        z = MySnapshot[3*i+2];
	if(specification==0||specificity(x,y,z,specifics)){
		if(diag==1){//If the simulation cell is orthogonal
		        nX = ((int)floor(x/Lx)+nCells)%nCells;//We place each atom in the corresponding cell
        		nY = ((int)floor(y/Ly)+nCells)%nCells;
        		nZ = ((int)floor(z/Lz)+nCells)%nCells;
			
        		DicoAtoms[3*i] = nX;
        		DicoAtoms[3*i+1] = nY;
        		DicoAtoms[3*i+2] = nZ;
		if(nX<0||nY<0||nZ<0){printf("A cell index is < 0. (atom %d) : nX %d ; nY %d ; nZ%d /x %f; y %f; z %f. Please check your UMD file ; XCart should always be >0.\n",i,nX,nY,nZ,x,y,z);}
		}
		else{//If it is not, the sub-cells aren't either
			rda = x * rvec[0] + y * rvec[1] + z * rvec[2];//We project the position in the 3 directions of the cell thanks to reciprocal vectors
			rdb = x * rvec[3] + y * rvec[4] + z * rvec[5];
			rdc = x * rvec[6] + y * rvec[7] + z * rvec[8];

        		nX = (int)floor(rda * nCells);//We place each atom in the corresponding cell
        		nY = (int)floor(rdb * nCells);
        		nZ = (int)floor(rdc * nCells);

        		DicoAtoms[3*i] = nX;
        		DicoAtoms[3*i+1] = nY;
        		DicoAtoms[3*i+2] = nZ;
		}
		nMesh[nX*nCells*nCells+nY*nCells+nZ]++;//We then increment the atom count of said cell
	}

    }


    int M = maxval(nMesh,nCells*nCells*nCells);          //maximum number of atoms that a cell contains in this snapshot
    Mesh = (int **)malloc(nCells*nCells*nCells*sizeof(int*)); //Mesh will explicitely contain the atoms of each sub-cell
    Bonding = (int **)malloc(nAtoms*sizeof(int*));        //Bonding will contain the atoms that are bound to another. Each sub-tab corresponds to a single (central) atom. All the values will be initialized to -1, except the first of each sub-tab which will be its size.
    int Ats[M];                                          //table transitorily containing the atoms in a central cell
    int Coord[M*27];                                     //table transitorily containing the coordinating atoms in the 27 surrounding cells

    for(int i=0 ; i<nCells*nCells*nCells ; i++){
	Mesh[i] = (int *)malloc(nMesh[i]*sizeof(int));
    }

    for(int i=0; i<nAtoms ; i++){
	Bonding[i]=(int *)malloc((2)*sizeof(int));

	Bonding[i][1]=-1;
	Bonding[i][0]=1;
    }

    for(int iat = 0 ; iat<nAtoms ; iat++){//Filling the Mesh tab using the information previousely stored in DicoAtoms

	if(DicoAtoms[3*iat] != -1){
	        nX = DicoAtoms[3*iat];
        	nY = DicoAtoms[3*iat+1];
        	nZ = DicoAtoms[3*iat+2];
        	cellnum = nX*nCells*nCells+nY*nCells+nZ;
        	index = MeshIndex[cellnum];
   		Mesh[cellnum][index] = iat;
		MeshIndex[cellnum]++;  				//Counts the number of atoms presents in each cell
	}
    }

    for(int iat=0 ; iat<nAtoms ; iat++) {//For each atom, we will determine the atoms to which it is bound
//	printf("atom %d\n",iat);

        numCoord = 0;

        if(DicoAtoms[3*iat]!=-1){//If we did not already directly examine the links between this atom and all the other atoms in it's sub-cell (or the atom is not in the specified zone)

            //Extraction of the relevant information

	    nX = DicoAtoms[3*iat];
            nY = DicoAtoms[3*iat+1];
            nZ = DicoAtoms[3*iat+2];
            cellnum = nX*nCells*nCells+nY*nCells+nZ;
            index = MeshIndex[cellnum];
            
	    for(int iatom=0 ; iatom<index ; iatom++){//We note that all the atoms in the sub-cell will be examined in their potential bonding relationships to each other

              
	        Ats[iatom]=Mesh[cellnum][iatom];
		DicoAtoms[3*Ats[iatom]]=-1;
            }

            for(int i=-1 ; i<2 ; i++){
                    cx = (nX+i+nCells)%nCells;
                    for(int j=-1 ; j<2 ; j++){
                        cy = (nY+j+nCells)%nCells;
                        for(int k=-1 ; k<2 ; k++){
                            cz = (nZ+k+nCells)%nCells;
                            cellnum = cx*nCells*nCells+cy*nCells+cz;
                            indexCoord = MeshIndex[cellnum];
                         
			    for(int jatom = 0 ; jatom<indexCoord ; jatom++){
	
                                Coord[numCoord]=Mesh[cellnum][jatom];//We fill the Coord tab with all the atoms in the surrounding cells + the central cell
                                numCoord++;			      //Counting the atoms in the tab
			
			       }
                            }
                        }
                    }


                 for(int iatom = 0 ; iatom<index ; iatom++){//We parse the Ats tab

                    at = Ats[iatom];
		    
                    for(int jatom = 0 ; jatom<numCoord ; jatom++){//And we examine their relationships to the atoms in the Coord tab

                        ne = Coord[jatom];

                        if(at<ne){
	    
	                    dx = MySnapshot[at*3]-MySnapshot[ne*3];
        	            dy = MySnapshot[at*3+1]-MySnapshot[ne*3+1];
                            dz = MySnapshot[at*3+2]-MySnapshot[ne*3+2];

			    if(diag==1){
					
				    if(fabs(dx/cell[0])>0.5){
				    	dx=cell[0]-fabs(dx);    
				    }
				    if(fabs(dy/cell[1])>0.5){
				    	dy=cell[1]-fabs(dy);    
				    }
				    if(fabs(dz/cell[2])>0.5){
				    	dz=cell[2]-fabs(dz);    
				    }

        	                    distij = pow(dx,2)+pow(dy,2)+pow(dz,2);
			    }

			    else{
			    	rda = dx * rvec[0] + dy * rvec[1] + dz * rvec[2];//We project the position difference in the 3 directions of the cell thanks to reciprocal vectors (remark : 1 > (rda, rdb, rdc) > 0)
			    	rdb = dx * rvec[3] + dy * rvec[4] + dz * rvec[5];
			    	rdc = dx * rvec[6] + dy * rvec[7] + dz * rvec[8];
	
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
			
			    	dx = rda * cell[0] + rdb * cell[3] + rdc * cell[6];//We rebuild the position difference in the cartesian space
			    	dy = rda * cell[1] + rdb * cell[4] + rdc * cell[7];
			    	dz = rda * cell[2] + rdb * cell[5] + rdc * cell[8];

                            	valx = pow(dx,2);
                            	valy = pow(dy,2);
                            	valz = pow(dz,2);
	
                            	distij = valx+valy+valz;
			    }

                            if(distij<BondTable[CrystalTypes[at]*ntypat+CrystalTypes[ne]]){//If at and ne are effectively linked to each other, we update Bonding and nBonding
				if(Bonding[at][0]<=nBonding[at]+1){
					Bonding[at] = (int *)realloc(Bonding[at],(Bonding[at][0]+2)*sizeof(int));
					Bonding[at][0] = Bonding[at][0]+2;
				}
				if(Bonding[ne][0]<=nBonding[ne]+1){
					Bonding[ne] = (int *)realloc(Bonding[ne],(Bonding[ne][0]+2)*sizeof(int));
					Bonding[ne][0] = Bonding[ne][0]+2;
				}

                                nBonding[ne]++;
                                nBonding[at]++;
				Bonding[at][nBonding[at]]=ne;
                                Bonding[ne][nBonding[ne]]=at;

				ntotAts = ntotAts+2;
                            }
                        }
                    }
                }
            }	
        }

	for(int i=0 ; i<nAtoms ; i++){//Now the first number of each tab is the number of coordinating atoms
		Bonding[i][0]=nBonding[i];
	}

	for(int i=0 ; i<nCells*nCells*nCells ; i++){free(Mesh[i]);}

	free(Mesh);
	free(nMesh);
	free(MeshIndex);
	free(DicoAtoms);
	free(nBonding);

	return Bonding;
}