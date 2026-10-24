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

bool specificity(double x, double y, double z, const double *specifics){//Check if a particular point meets the spatial selectivity criteria

	if(x>=specifics[0] && x<=specifics[0]+specifics[3] && y>=specifics[1] && y<=specifics[1]+specifics[4] && z>=specifics[2] && z<=specifics[2]+specifics[5]){
		return true;
	}
	else{
		return false;
	}
}

static void add_bond(int **Bonding, int *nBonding, int at, int ne){//Appends ne to the list of at, Bonding[at][0] being its allocated length during the search
	if(Bonding[at][0]<=nBonding[at]+1){
		Bonding[at] = (int *)realloc(Bonding[at],(Bonding[at][0]+4)*sizeof(int));
		Bonding[at][0] += 4;
	}
	nBonding[at]++;
	Bonding[at][nBonding[at]] = ne;
}

/* Bonds of one snapshot, for any cell (orthogonal or not).
   MySnapshot   : cartesian coordinates, nAtoms x 3
   BondTable    : squared cutoff of each pair of types, ntypat x ntypat
   cell         : lattice vectors, rows a, b, c (cartesian = s0 a + s1 b + s2 c)
   rvec         : reciprocal vectors, rows a*, b*, c* (s_k = x.rvec_k), the transpose of the inverse of cell
   nCells       : number of sub-cells along a, b, c ; the width of a sub-cell (distance between its opposite faces) must not be
                  smaller than the largest cutoff
   Each atom is put in a sub-cell from its reduced coordinates, and compared with the atoms of the 27 surrounding sub-cells,
   each counted once even when there are fewer than 3 sub-cells along a direction (then -1 and +1 are the same sub-cell).
   The nearest image is taken in reduced coordinates (ds - rint(ds)), exact for the distances shorter than half the smallest
   width of the cell. Returns for each atom a list whose first number is the number of bound atoms, followed by them. */
int** compute_Bonds_full(const double *MySnapshot, const double *BondTable, const int *CrystalTypes, const int nAtoms, const int ntypat, const double *cell, const double *rvec, const int *nCells, const double *specifics, const int specification){
	int nTot = nCells[0]*nCells[1]*nCells[2];
	double *red = malloc(3*nAtoms*sizeof(double));//reduced coordinates, in [0,1)
	int *cellOf = malloc(nAtoms*sizeof(int));      //sub-cell of each atom, -1 outside the specified zone
	int *head = malloc(nTot*sizeof(int));          //first atom of each sub-cell, then next[] chains the others
	int *next = malloc(nAtoms*sizeof(int));
	int *nBonding = calloc(nAtoms,sizeof(int));
	int **Bonding = (int **)malloc(nAtoms*sizeof(int*));
	int offsets[3][3], nOffsets[3];

	for(int k=0 ; k<3 ; k++){//Distinct neighbouring sub-cells along each direction
		if(nCells[k]>=3){
			nOffsets[k] = 3;
			offsets[k][0] = -1; offsets[k][1] = 0; offsets[k][2] = 1;
		}
		else{
			nOffsets[k] = nCells[k];
			for(int o=0 ; o<nCells[k] ; o++){offsets[k][o] = o;}
		}
	}
	for(int c=0 ; c<nTot ; c++){head[c] = -1;}

	for(int i=nAtoms-1 ; i>=0 ; i--){//Backwards, so that each sub-cell lists its atoms in increasing order
		const double *x = MySnapshot+3*i;
		Bonding[i] = (int *)malloc(4*sizeof(int));
		Bonding[i][0] = 4;
		cellOf[i] = -1;
		if(specification==0 || specificity(x[0],x[1],x[2],specifics)){
			int n[3];
			for(int k=0 ; k<3 ; k++){
				double s = x[0]*rvec[3*k]+x[1]*rvec[3*k+1]+x[2]*rvec[3*k+2];
				s -= floor(s);
				red[3*i+k] = s;
				n[k] = (int)(s*nCells[k]);
				if(n[k]>=nCells[k]){n[k] = nCells[k]-1;}
			}
			cellOf[i] = (n[0]*nCells[1]+n[1])*nCells[2]+n[2];
			next[i] = head[cellOf[i]];
			head[cellOf[i]] = i;
		}
	}

	for(int at=0 ; at<nAtoms ; at++){
		if(cellOf[at]<0){continue;}
		int c = cellOf[at];
		int n0 = c/(nCells[1]*nCells[2]), n1 = (c/nCells[2])%nCells[1], n2 = c%nCells[2];
		const double *cut = BondTable+CrystalTypes[at]*ntypat;
		for(int i=0 ; i<nOffsets[0] ; i++){
			int cx = (n0+offsets[0][i]+nCells[0])%nCells[0];
			for(int j=0 ; j<nOffsets[1] ; j++){
				int cy = (n1+offsets[1][j]+nCells[1])%nCells[1];
				for(int k=0 ; k<nOffsets[2] ; k++){
					int cz = (n2+offsets[2][k]+nCells[2])%nCells[2];
					for(int ne=head[(cx*nCells[1]+cy)*nCells[2]+cz] ; ne>=0 ; ne=next[ne]){
						if(ne<=at){continue;}//each pair once
						double da = red[3*ne]-red[3*at], db = red[3*ne+1]-red[3*at+1], dc = red[3*ne+2]-red[3*at+2];
						da -= rint(da);
						db -= rint(db);
						dc -= rint(dc);
						double dx = da*cell[0]+db*cell[3]+dc*cell[6];
						double dy = da*cell[1]+db*cell[4]+dc*cell[7];
						double dz = da*cell[2]+db*cell[5]+dc*cell[8];
						if(dx*dx+dy*dy+dz*dz < cut[CrystalTypes[ne]]){
							add_bond(Bonding,nBonding,at,ne);
							add_bond(Bonding,nBonding,ne,at);
						}
					}
				}
			}
		}
	}

	for(int i=0 ; i<nAtoms ; i++){//Now the first number of each list is the number of bound atoms
		Bonding[i][0] = nBonding[i];
	}

	free(red);
	free(cellOf);
	free(head);
	free(next);
	free(nBonding);

	return Bonding;
}
