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
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>
#include "include/interval.h"

string sbuf;
FILE *f;
FILE *g;
string outputFormat = "full";

bool process_interval(Interval i) {
	long range = i.end - i.start + 1;
	if (range <= 0) return false;


	double * fs= (double*)malloc(sizeof(double) * range);
	if (fs==NULL) {
		printf("Cannot malloc mem\n");
		exit(1);
	}
	double * gs= (double*)malloc(sizeof(double) * range);
	if (gs==NULL) {
		printf("Cannot malloc mem\n");
		exit(1);
	}


	long go_to = (i.start - 1) * sizeof(double);
	int res = fseek(f, go_to, SEEK_SET);
	if (res != 0) {
		cerr << "doc_walker: error at interval " << i << endl;
		exit(1);
	}
	res = fseek(g, go_to, SEEK_SET);
	if (res != 0) {
		cerr << "doc_walker: error at interval " << i << endl;
		exit(1);
	}
	res = fread(fs, sizeof(double) * range, 1, f);
	if (res != 1) {
		cerr << "doc_walker: error at interval " << i << endl;
		exit(1);
	}
	res = fread(gs, sizeof(double) * range, 1, g);
	if (res != 1) {
		cerr << "doc_walker: error at interval " << i << endl;
		exit(1);
	}
	double sum_unmasked = 0.0;
	double sum_masked = 0.0;
	double sum_expected = 0.0;
	int masked = 0;
	int index = 0;
	while (index < range) {
		//printf("%d\t%.4f\t%.4f\n",index+start,fs[index],gs[index]);
		sum_unmasked += fs[index];
		if (gs[index] > 0.0) {
			sum_masked += fs[index];
			sum_expected += gs[index];
		} else {
			masked++;
		}
		index++;
	}
	double doc_ratio = (double) sum_masked / (double) sum_expected;
	cout << i << "\t";
	if (outputFormat == "concise") {
		printf ("%.4f\n", doc_ratio);
	} else if (outputFormat == "full") {
		printf ("DOC_ratio:\t%.4f\tObserved_coverage:\t%.4f\tExpected_coverage:\t%.4f\tNum_masked_positions:\t%d\tObserved_coverage_including_masked_regions:\t%.4f\n", doc_ratio, sum_masked, sum_expected, masked, sum_unmasked);
	} else if (outputFormat == "internal") {
		printf("%.4f\t%.4f\t%d\n", sum_masked, sum_expected, range - masked);
	} else {
		cerr << "Unknown output format: " << outputFormat << endl;
		exit(1);
	}


	delete fs;
	delete gs;
	return true;

}

int main(int argc, char** argv) {
	unsigned int start = -1;
	unsigned int end = -1;
	if (argc == 5) {
		start = atoi(argv[3]);
		end   = atoi(argv[4]);
	} else if (argc == 6) {
		start = atoi(argv[3]);
		end   = atoi(argv[4]);
		outputFormat = argv[5];
	} else if (argc == 4) {
		outputFormat = argv[3];
	} else if (argc != 3) {
		printf("%s scov_file gc_file [start end] [concise|internal] \n", argv[0]);
		printf("The .scov and .gc files can be found in the work_dir.\n");
		printf("If start/end is not specified, an interval file is taken as an input.\n");
		exit(1);
	}

	f = fopen(argv[1],"rb");
	if (f==NULL){
		printf("Error opening file\n");
		exit(1);
	}
	g = fopen(argv[2],"rb");
	if (g==NULL){
		printf("Error opening file\n");
		exit(1);
	}

	if (start == -1) {
		while (getline(cin, sbuf)) {
			Interval i = read_interval(sbuf);
			if (!process_interval(i)) cout << sbuf << endl;
		}
	} else {
		process_interval(Interval("-1", start, end));
	}
	fclose(f);
	fclose(g);
	return 0;
}
