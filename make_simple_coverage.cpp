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
//#include <fcntl.h>
//#include <sys/mman.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
//#include <sys/types.h>
//#include <sys/stat.h>

#include "include/dbtypes.h"
#include "include/defs.h"

void usage(int argc, char** argv) {
	cout << "Usage: " << argv[0] << " rmap_file output_scov_file ref_name ref_length\n";
	exit(1);
}


int main (int argc, char** argv) {
	if (argc != 5) usage(argc, argv);

	string rmap_filename = argv[1];
	string scov_filename = argv[2];
	string ref_name = argv[3];
	uint64_t ref_len = atol(argv[4]);
	
	double * sum_normodds = new double[ref_len];
	for (int i = 0; i < ref_len; i++) sum_normodds[i] = 0;

	//Read in rmap file
	ifstream inf;
	open_file_binary(inf, rmap_filename);
	Map_t mapt;
	while (inf) {
		inf.read((char*) & mapt, sizeof(mapt));
		if (mapt.ref_name == ref_name) {
			if (mapt.ref_pos > ref_len) {
				cerr << argv[0] << ": Error, out of range position " << mapt.ref_pos << " on ref " << ref_name << " of size " << ref_len << ".\n";
				exit(1);
			}
			sum_normodds[mapt.ref_pos - 1] += mapt.normodds;
		}
	}
	inf.close();

	//Write out scov file
	ofstream outf;
	open_file_binary(outf, scov_filename);
	outf.write((char*) sum_normodds, ref_len * sizeof(double));
	if (outf.bad()) {
		cerr << argv[0] << ": Error writing scov file " << scov_filename << ".\n";
		exit(1);
	}
	outf.close();
	delete[] sum_normodds;
	return 0;
}
