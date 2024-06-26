#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>
#include <string.h>

//This script is intended to work with the python script umd_processes_fast.py
//It's main purpose is to read a snapshot and to extract the wanted information and only this

int min(int a,int b){

	if(a>b){return b;}
	else{return a;}

}	

double min3(double a, double b, double c){

	double minimum = a;

	if(minimum>b){minimum=b;}
	if(minimum>c){minimum=c;}

	return minimum;

}


void free_memory(double *Tab){
	free(Tab);
}


double* read_umd_values(char *MySnapshot, int nAtoms, int X){//reads some specifics values of a snapshot and returns it in a tab

	int len = 0;
	int depth = 0;
	int atom = 0;
	int flagcoord = 0;//If we're on the right place on the line to read the values
	int flagexp = 0;//If we're reading an exponent
	int sign = 1;//Takes in account the fact that the value can be negative
	int signexp = 1;//And the exponent too
	double number = 0;
	int exp = 0;
	double* Values;//This tab will contain the values and will be returned at the end of the function
	int lineIndex=0;//index on the line of the number we are looking at
	
	if(X!=-1){
		Values = calloc(3*nAtoms,sizeof(double));//3 values per atom
		 }
	else{
		Values = calloc(12*nAtoms,sizeof(double));//Here, all the values (xred, xcart, absxcart and vels) will be read and returned
	    }
	len = strlen(MySnapshot);
	for(int i=0 ; i<len ; i++){

		if(MySnapshot[i]==' '){	//If the next symbol is a space, it means we're at the end of the number
			if(flagcoord){
				if(atom<nAtoms){
					if(X!=-1){//Different ways of ordering in function of the value of X
						Values[3*atom+lineIndex-X]=number*sign*pow(10,exp*signexp);
					}
					else if(lineIndex<12){
						Values[12*atom+lineIndex]=number*sign*pow(10,exp*signexp);
					}
				}
				flagcoord=0;
			}
			//We reinitialize all the values of the variables for the next number

			exp=0;
			flagexp=0;
			number=0;
			sign = 1;
			signexp = 1;

			lineIndex++;//We increment the indice by 1
			depth=0;//The next digits will be in the integer part
		}
		else if(MySnapshot[i]=='-'){//If there's a minus sign, we update the information
			if(flagexp){//If we're on the exponent part
				signexp=-1;
			}
			else{
				sign = -1;
			}
		}
		else if((lineIndex==X||lineIndex==X+1||lineIndex==X+2)||(X==-1 && lineIndex<12)){//If we are on the right coordinates
			flagcoord = 1;
			if(isdigit(MySnapshot[i])){
				if(flagexp){//We built the exponent
					exp = exp*10;
					exp = exp + (MySnapshot[i]-'0');
				}
				else if(depth == 0){//We build the int part of the number
					number=number*10;
					number=number+(MySnapshot[i]-'0');
				}
			
				else{//We're on the decimal part, so we divide the number by a power of 10
					number=number+(double)(MySnapshot[i]-'0')/depth;
					depth=depth*10;

				}
			}
			else if(MySnapshot[i]=='.'){
				depth=10;//We're in the decimal part now
			}
			else if(MySnapshot[i]=='e'){
				exp=0;
				flagexp = 1;//From now the digits we encounter are the exponent's
			}
			

		}
		else if(MySnapshot[i]=='\n'){
			//We increment the atom index and reinitialize all the other variables for the next line
			atom++;
			lineIndex=0;
			exp=0;
			sign=1;
			signexp=1;
			depth=0;
			number=0;
			flagexp=0;
			flagcoord = 0;
		}
	}
return Values;
}


bool is_in(int atom,int *Indexes, int n){//Is an atom included in the ensemble described by a list of indexes ?
	bool result = false;
	for(int i=0 ; i<n ; i++){
		if((atom>=Indexes[2*i] && atom<=Indexes[2*i+1])){
			result = true;
		}
	}
	return result;
}

int num_atoms(int* Indexes, int n){//Return the number of atoms in the aforementioned ensemble
	int num = 0;
	for(int i=0 ; i<n ; i++){
		num = num + Indexes[2*i+1]-Indexes[2*i] + 1;
	}
	return num;
}

int* defineBonds(const char *lines,int len,int *CentIndexes,int *OutIndexes, int nCent, int nOut,int natom){//Read the bonds of a bond file and return the info in a tab
	
	int N=0;
	int number=0;
	int nAts=0;
	int nMax=0;
	bool ligand = false ;//true when the examined atom is on the list of the elements we're looking
	bool newline = true ;//true when we're looking at the first number of the line
	int *BondsList;	     //Contains the atoms bound to the other atoms. BondsList[0] will contain the index in BondsList of the last atom written.
	int *BondIndexes;    //BondsIndexes[X] contains the first index in the BondsList pertaining to the atom X for 0<X<natom. BondsIndexes[-1] contains the maximum number of coordinating atoms for any atom. BondsIndexes[0] contains the indexe in BondsIndexes of the last number written.
	int numCent = 0;     //Will be initialized to the number of central atoms
	int numOut = 0;	     //Will be initialized to the number of outer atoms
	int iatom = 0;
	
	bool OnAtoms = true;	

	numCent = num_atoms(CentIndexes,nCent);
	numOut = num_atoms(OutIndexes,nOut);

	//Max possible sizes
	BondsList = calloc((long long int)(numCent+numOut)*(numCent+numOut+2),sizeof(int));
	BondIndexes = calloc(natom+3,sizeof(int));
	
	BondIndexes[0]=1;
	BondIndexes[1]=0;
	BondsList[0]=1;
	for(int i=0 ; i<len ; i++){


		if(isdigit(lines[i])){number = number*10; number = number + (lines[i]-'0');}//If the character is a digit, we built the number corresponding to the atom

		else if(lines[i] == '\t' && OnAtoms){				//If it's a tab, we treat the lastly built number as the atom it embodies
			if(newline||ligand){//If the atom is the first of the line (central) OR if it's among the coordinating atoms when the first of the line is a ligand
				if(is_in(number,CentIndexes,nCent)||is_in(number,OutIndexes,nOut)){	//If the atom is a ligand
					if(newline){				//If it's the first of the line, all the others of this one will be examined
						newline = false;
						ligand = true;
					}		
					else{					//If it's not, it's bound to the first one ; we update the BondsList.
						BondsList[BondsList[0]]=number; 
						BondsList[0]++;
						nAts++;				//Counts the number of atoms that are bound to the first one.
					}
				}	
				else{newline = false ; number=0;}		//We won't examine the non-ligands atoms.
			}
			number=0;//reinitialization of the number

		}

		else if(lines[i] == '\n' && OnAtoms){//If we're going toward the next line and we're still on the "atoms" part
			
			iatom  ++;
			if(iatom == natom){OnAtoms = false;}	//If we've examined all the atoms, the next lines won't contain any relevant information
			if(nAts>nMax){nMax=nAts;}		//Updating the "maximum number of bonds for any atoms" variable
			BondIndexes[BondIndexes[0]+1]=BondIndexes[BondIndexes[0]]+nAts;//And the indexes
			BondIndexes[0]++;

			//Reinitializing variables for the next line
			nAts=0;
			ligand = false;
			newline = true;
			number=0;
		}
	}

	
	BondIndexes[BondIndexes[0]+1]=nMax;
	BondIndexes[0]++;

	int *BList;//Will contain all the information in a single list.

	BList = malloc((BondsList[0]+BondIndexes[0]+2)*sizeof(int));

	for(int i=1 ; i<BondsList[0]+1 ; i++){
		BList[i+1]=BondsList[i];
	}

	for(int i=1 ; i<BondIndexes[0]+1 ; i++){
		BList[BondsList[0]+1+i]=BondIndexes[i];
	}
	
	BList[0]=BondsList[0];
	BList[1]=BondIndexes[0];

	free(BondsList);
	free(BondIndexes);

	return BList ;
}
