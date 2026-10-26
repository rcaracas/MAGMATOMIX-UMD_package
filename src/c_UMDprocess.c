#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>

//Readers of UMD files and bonding files, used by umd_processes_fast.py (ctypes)


void free_memory(void *Tab){
	free(Tab);
}


static long read_line(FILE *f, char **buf, size_t *size){//Reads one line of any length into *buf (grown as needed) ; returns its length, -1 at the end of the file
	size_t len;
	if(fgets(*buf,(int)*size,f)==NULL){return -1;}
	len = strlen(*buf);
	while(len>0 && (*buf)[len-1]!='\n'){
		*size *= 2;
		*buf = realloc(*buf,*size);
		if(fgets(*buf+len,(int)(*size-len),f)==NULL){break;}
		len += strlen(*buf+len);
	}
	return (long)len;
}


static const char* first_word(const char *line){//The line without its leading blanks
	while(*line==' ' || *line=='\t'){line++;}
	return line;
}


static bool starts_with(const char *line, const char *word){//The first word of the line is word
	size_t n = strlen(word);
	return strncmp(line,word,n)==0 && (line[n]==' ' || line[n]=='\t' || line[n]=='\n' || line[n]=='\r' || line[n]=='\0');
}


static int read_numbers(const char *line, int n, double *values){//The first n numbers after the first word of the line ; returns how many were read
	char *end;
	int k;
	while(*line && !isspace((unsigned char)*line)){line++;}
	for(k=0 ; k<n ; k++){
		values[k] = strtod(line,&end);
		if(end==line){break;}
		line = end;
	}
	return k;
}


long umd_count(const char *filename, long natom, long *total){//Number of complete snapshots of a UMD file (an "atoms:" line followed by natom lines) ;
	FILE *f = fopen(filename,"r");                            //total = number of "atoms:" lines (the last snapshot may be incomplete) ; -1 if the file cannot be opened
	size_t size = 65536;
	char *buf;
	long complete = 0;
	long len, iatom;
	*total = 0;
	if(f==NULL){return -1;}
	buf = malloc(size);
	while((len=read_line(f,&buf,&size))>=0){
		if(starts_with(first_word(buf),"atoms:")){
			(*total)++;
			for(iatom=0 ; iatom<natom ; iatom++){
				if(read_line(f,&buf,&size)<0){break;}
			}
			if(iatom==natom){complete++;}
		}
	}
	free(buf);
	fclose(f);
	return complete;
}


long umd_read(const char *filename, long natom, long first, long last, long sampling, int col0, int ncol, double *out, double *cells){
	//Snapshots first, first+sampling, ... before last (counted from 0, complete snapshots only) : columns col0 to col0+ncol-1 of the atom lines
	//(reduced*3, cartesian*3, unwrapped*3, velocity*3, ...) into out (snapshot, atom, column), and if cells is not NULL the lattice vectors of each
	//of these snapshots (rows a, b, c) from its rprimd lines, else from its acell (orthogonal cell), else NaN.
	//Returns the number of snapshots stored ; -1 if the file cannot be opened ; -2-k if an atom line of the k-th stored snapshot has too few numbers
	FILE *f = fopen(filename,"r");
	size_t size = 65536;
	char *buf;
	const char *line;
	double rprimd[9], acell[3];
	bool hasrprimd[3] = {false,false,false};
	bool hasacell = false;
	long isnap = 0, nstored = 0, len, iatom;
	int k, j, row;
	if(f==NULL){return -1;}
	buf = malloc(size);
	while(isnap<last && (len=read_line(f,&buf,&size))>=0){
		line = first_word(buf);
		if(starts_with(line,"timestep")){
			hasrprimd[0] = hasrprimd[1] = hasrprimd[2] = false;
			hasacell = false;
		}
		else if(cells!=NULL && (starts_with(line,"rprimd_a") || starts_with(line,"rprimd_b") || starts_with(line,"rprimd_c"))){
			row = line[7]-'a';
			hasrprimd[row] = read_numbers(line,3,rprimd+3*row)==3;
		}
		else if(cells!=NULL && starts_with(line,"acell")){
			hasacell = read_numbers(line,3,acell)==3;
		}
		else if(starts_with(line,"atoms:")){
			bool wanted = isnap>=first && (isnap-first)%sampling==0;
			for(iatom=0 ; iatom<natom ; iatom++){
				if(read_line(f,&buf,&size)<0){break;}
				if(!wanted || ncol==0){continue;}
				char *p = buf, *end;
				double *values = out + (nstored*natom+iatom)*ncol;
				for(k=0 ; k<col0 ; k++){//skip the columns before col0
					while(isspace((unsigned char)*p)){p++;}
					while(*p && !isspace((unsigned char)*p)){p++;}
				}
				for(k=0 ; k<ncol ; k++){
					values[k] = strtod(p,&end);
					if(end==p){
						free(buf);
						fclose(f);
						return -2-nstored;
					}
					p = end;
				}
			}
			if(iatom<natom){break;}//incomplete last snapshot
			if(wanted){
				if(cells!=NULL){
					for(j=0 ; j<9 ; j++){
						if(hasrprimd[0] && hasrprimd[1] && hasrprimd[2]){cells[9*nstored+j] = rprimd[j];}
						else if(hasacell){cells[9*nstored+j] = (j%4==0) ? acell[j/4] : 0.0;}
						else{cells[9*nstored+j] = NAN;}
					}
				}
				nstored++;
			}
			isnap++;
		}
	}
	free(buf);
	fclose(f);
	return nstored;
}


int* defineBonds(const char *lines,int len,int *CentIndexes,int *OutIndexes, int nCent, int nOut,int natom){//Read the bonds of a snapshot of a bonding file
	//(one line per atom : the atom, then its bound atoms, each followed by a tab) between the central and outer atoms (ranges of indexes, inclusive).
	//Returns BList (to free with free_memory) : BList[0] = 1 + the number of bonds stored, BList[1] = natom + 2 ; the bound atoms of the central and
	//outer atoms, atom after atom ; the index in this list where the bound atoms of each atom begin (natom+1 values) and the largest number of bound atoms
	int number=0;
	int nAts=0;
	int nMax=0;
	bool ligand = false ;//true when the examined atom is on the list of the elements we're looking
	bool newline = true ;//true when we're looking at the first number of the line
	int *BondsList;	     //Contains the atoms bound to the other atoms. BondsList[0] will contain the index in BondsList of the last atom written.
	int *BondIndexes;    //BondsIndexes[X] contains the first index in the BondsList pertaining to the atom X for 0<X<natom. BondsIndexes[-1] contains the maximum number of coordinating atoms for any atom. BondsIndexes[0] contains the indexe in BondsIndexes of the last number written.
	int iatom = 0;
	bool OnAtoms = true;

	//Every bound atom takes at least two characters (a digit and a tab)
	BondsList = calloc((size_t)len/2+2,sizeof(int));
	BondIndexes = calloc(natom+3,sizeof(int));

	BondIndexes[0]=1;
	BondIndexes[1]=0;
	BondsList[0]=1;
	for(int i=0 ; i<len ; i++){

		if(isdigit((unsigned char)lines[i])){number = number*10; number = number + (lines[i]-'0');}//If the character is a digit, we built the number corresponding to the atom

		else if(lines[i] == '\t' && OnAtoms){				//If it's a tab, we treat the lastly built number as the atom it embodies
			if(newline||ligand){//If the atom is the first of the line (central) OR if it's among the coordinating atoms when the first of the line is a ligand
				bool in = false;
				for(int k=0 ; k<nCent ; k++){if(number>=CentIndexes[2*k] && number<=CentIndexes[2*k+1]){in = true;}}
				for(int k=0 ; k<nOut ; k++){if(number>=OutIndexes[2*k] && number<=OutIndexes[2*k+1]){in = true;}}
				if(in){	//If the atom is a ligand
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
