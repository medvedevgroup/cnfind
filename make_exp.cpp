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
string gcbins_filename;
string exp_filename;
int window_size;
double * scov;
char * masks;

void usage(int argc, char** argv) {
	cout << "Usage: " << argv[0] << " fasta_filename masks_filename scov_filename gc_filename window_size num_bins exp_filename\n";
	exit(1);
}


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


void read_gc_bins( double bins_lambdas[], string filename) {
	ifstream inf;
	string sbuf;
	open_file(inf, filename);
	getline(inf, sbuf); //skip_header
	int num_read = 0;
	vector<string> row;
	while (get_row(inf, row)) {
		num_read++;
		if (num_read > num_bins) {
			cerr << "Attempting to read in more bins than allocated...error!\n";
			exit(1);
		}
		bins_lambdas[num_read - 1] = atof(row[2].c_str());
	}
	assert(num_read == num_bins);
	//cerr << "bins are: \n"; for (int i = 0; i < num_bins; i++) { cerr << bins_lambdas[i] << "\n"; }


}


int main (int argc, char** argv) {
	if (argc!=8) usage(argc, argv);

	string fasta_filename    = argv[1];
	string masks_filename    = argv[2];
	string scov_filename     = argv[3];
	gcbins_filename          = argv[4];
	window_size              = atoi(argv[5]);
	num_bins                 = atoi(argv[6]);
	exp_filename             = argv[7];


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


	int i;

	read_gc_bins(bins_lambdas, gcbins_filename);

	//calculate pos2bins, i.e. the gc content at every position.
	int * pos2bin=(int*)malloc(sizeof(int)*sequence.length());
	for (i=0; i<sequence.length(); i++) pos2bin[i]=-1;

	for (i = 0; i < sequence.length(); ) {
		int left = sequence.length() - i;
		if (left > window_size) left = window_size;

		int gc_content=0;
		int at_content=0;

		int j;
		for (j=0; j<left; j++) {
			char c = sequence[i+j];
			if (c=='c' || c=='C' || c=='g' || c=='G')  gc_content+=1;
			if (c=='a' || c=='A' || c=='t' || c=='T') at_content+=1;
		}	

		int bin = whats_my_bin(gc_content,at_content);
		for (j=0; j<left; j++) {
			char c=sequence[i+j];
			if (masks[i+j]  == '1' || c=='N' || c=='n') {
				pos2bin[i+j] = -1;
			} else {
				pos2bin[i+j] = bin;
			}
		}

		i += window_size;
	}

	//calculate pos2lambda : expected arrivals per position
	delete[] masks;
	delete[] scov;
	double* pos2lambda = (double*)malloc(sizeof(double)*(sequence.length()));
	if (pos2lambda == NULL) {
		fprintf(stderr,"Trouble allocating pos2lambda\n");
		exit(1);
	}
	unsigned int repeats = 0;
	unsigned int not_repeats=0;
	unsigned int shifted_positions=0;
	for(i = 0; i < sequence.length(); i++) {
		int bin = pos2bin[i];
		//cout << "Looking at position " << i << endl;
		if (bin >= 0) {
			not_repeats++;
			if (bins_lambdas[bin] == 0.0 || bins_lambdas[bin] == -2.0) {
				shifted_positions++;
			}
			//printf("Thinking about bin %d, %f\n",bin,bins_lambdas[bin]);
			if (bin > (num_bins/2)) {
				while(bins_lambdas[bin] == 0.0 || bins_lambdas[bin] == -2.0)  bin--;
			} else {
				while(bins_lambdas[bin] == 0.0 || bins_lambdas[bin] == -2.0)  bin++;
			}

			if ( bin >= num_bins  || bin < 0 ) { 
				pos2lambda[i] = 0.0;
			} else {
				pos2lambda[i] = bins_lambdas[bin];
				assert(pos2lambda[i] > 0.0);
			}
		} else {
			assert(bin == -1);
			repeats++;
			pos2lambda[i] = -1.0;
		}
	}

	//write out gc file
	ofstream outf;
	open_file_binary(outf, exp_filename);
	outf.write((char *) pos2lambda, sizeof(double) * sequence.length());
	outf.close();

	//printf("Repeats: %d, Un-Repeats: %d\n",repeats, not_repeats);
	//printf("Shifted positions: %d\n",shifted_positions);

	free(pos2lambda);
	return 0;

}
