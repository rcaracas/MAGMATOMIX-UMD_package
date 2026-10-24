#include <math.h>

/* Histograms of the interatomic distances, for any cell (orthogonal or not), over a block of snapshots.
   red     : reduced coordinates, nsnap x natom x 3
   cells   : lattice vectors of each snapshot, nsnap x 3 x 3, rows a, b, c (cartesian = s0 a + s1 b + s2 c)
   weights : weight of each snapshot in wcounts (its volume, for the normalization of g(r) when the cell changes)
   counts, wcounts : ntypes x ntypes x nbins, filled by the caller with zeros ; every pair is counted in both orders (i,j) and (j,i)
   The nearest image is taken in reduced coordinates (ds - rint(ds)), which is exact for all the distances shorter than half the
   smallest perpendicular width of the cell : nbins*dr must not exceed it. */
void gofr_histogram(const double *red, int nsnap, int natom, const int *typat, int ntypes, const double *cells,
                    const double *weights, double dr, int nbins, double *counts, double *wcounts){
	double rmax2 = (nbins*dr)*(nbins*dr);
	for(int s=0 ; s<nsnap ; s++){
		const double *r = red + (long)s*natom*3;
		const double *m = cells + 9*s;
		double w = weights[s];
		for(int i=0 ; i<natom-1 ; i++){
			double xi = r[3*i], yi = r[3*i+1], zi = r[3*i+2];
			int ti = typat[i];
			for(int j=i+1 ; j<natom ; j++){
				double da = r[3*j]-xi, db = r[3*j+1]-yi, dc = r[3*j+2]-zi;
				da -= rint(da);
				db -= rint(db);
				dc -= rint(dc);
				double dx = da*m[0]+db*m[3]+dc*m[6];
				double dy = da*m[1]+db*m[4]+dc*m[7];
				double dz = da*m[2]+db*m[5]+dc*m[8];
				double d2 = dx*dx+dy*dy+dz*dz;
				if(d2 < rmax2){
					int bin = (int)(sqrt(d2)/dr);
					if(bin >= nbins) continue;
					int tj = typat[j];
					long k1 = ((long)ti*ntypes+tj)*nbins+bin;
					long k2 = ((long)tj*ntypes+ti)*nbins+bin;
					counts[k1] += 1.0;
					counts[k2] += 1.0;
					wcounts[k1] += w;
					wcounts[k2] += w;
				}
			}
		}
	}
}
