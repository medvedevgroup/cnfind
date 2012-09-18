/*Copyright 2012, Paul Medvedev

This file is part of Cnfind

Cnfind is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cnfind is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with cnfind (see file gpl.txt). If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <fcntl.h>
#include <sys/mman.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "include/dbtypes.h"
#include "include/defs.h"

int num_bins;
string sequence;
string out_file;
int window_size;
double * scov;
char * masks;

int whats_my_bin(int gc_count,int at_count) {
	if (gc_count==0 && at_count==0) {
		return -1;
	}
	int bin = ((gc_count * num_bins) / (at_count+gc_count + 1));
	if (bin >= num_bins || bin < 0) {
		fprintf(stderr,"Window error!");
	}
	assert(bin >= 0 && bin < num_bins);
	return bin;
}


void usage(int argc, char** argv) {
	cout << "Usage: " << argv[0] << " fasta_filename masks_filename scov_filename out_gcbins_file window_size num_bins\n";
	exit(1);
}

void build_gc_bins(double bins_lambdas[]) {

	int i;
	double bins_hits[num_bins];
	double bins_length[num_bins];

	for (i=0; i<num_bins; i++) {
		bins_hits[i]=0.0;
		bins_length[i]=0.0;
		bins_lambdas[i]=0.0;
	}

	for (i = 0; i < sequence.length(); ) {
		int left = sequence.length() - i;
		if (left > window_size) left = window_size;

		int gc_content=0;
		int at_content=0;
		int repeats_inside=0;

		int j;
		for (j=0; j<left; j++) {
			char c=sequence[i+j];
			if (c=='c' || c=='C' || c=='g' || c=='G') {
				gc_content+=1;
			}
			if (c=='a' || c=='A' || c=='t' || c=='T') {
				at_content+=1;
			}
			if (masks[i+j] == '1' || c=='N' || c=='n')  repeats_inside++; 
		}	

		int bin = whats_my_bin(gc_content,at_content);

		/*for (j=0; j<left; j++) {
			if (masks[i+j] == '0') {
				bins_hits[bin] += scov[i+j];
				bins_length[bin]++; 
			}
		}
		*/

		if (repeats_inside==0) {
		  for (j=0; j<left; j++) bins_hits[bin] += scov[i+j];
		  bins_length[bin] += window_size;
		}
	
		i += window_size;
	}

	//compute bins_lambdas:  the arrival rate for each bin, and write out
	ofstream outf;
	open_file(outf, out_file);
	outf << "Bin\tRangeStart\tLambda\tHits\tLen\n";
	for (i=0; i<num_bins; i++) {
		if (bins_length[i] == 0) {
			bins_lambdas[i] = -2.0;
		} else {
			bins_lambdas[i] = bins_hits[i] / bins_length[i];
		}
		double bin_size = 1.0 / num_bins;
		outf << i << "\t" << bin_size * i << "\t" << bins_lambdas[i] << "\t" << bins_hits[i] << "\t"  << bins_length[i] << endl;
	}
	outf.close();
}

int main (int argc, char** argv) {
	if (argc != 7) usage(argc, argv);

	string fasta_filename    = argv[1];
	string masks_filename    = argv[2];
	string scov_filename     = argv[3];
	out_file                 = argv[4];
	window_size              = atoi(argv[5]);
	num_bins                 = atoi(argv[6]);

	//fprintf(stderr,"Using window_size %d, num_bins %d\n",window_size,num_bins);

	double bins_lambdas[num_bins];
	sequence = read_genome(fasta_filename); //read in fasta file
	//printf("Read in %d chars from %s\n",sequence_length,fasta_filename);

	ifstream inf;
	open_file_binary(inf, scov_filename);
	scov = new double[sequence.length()];
	inf.read((char *) scov, sizeof(double) * sequence.length());
	inf.close();

	//read in the masks
	open_file_binary(inf, masks_filename);
	masks = new char[sequence.length()];
	if (masks ==NULL) {
		fprintf(stderr,"Problem with allocating memory!\n");
		exit(1);
	}
	inf.read(masks, sequence.length());
	inf.close();


	build_gc_bins(bins_lambdas);
	delete[] scov;
	delete[] masks;

	return 0;

}
