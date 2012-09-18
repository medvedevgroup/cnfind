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

#include<cassert>
#include<cmath>
#include<cstring>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<list>
#include<queue>
#include<deque>
#include<cstdarg>
#include<algorithm>
#include<stdint.h>

using namespace std;
string sbuf;

#include "include/defs.h"
#include "include/interval.h"


/* Takes an interval query and outputs the gcContent of that range (%gc out of non-ambig bases, and %non-ambig bases)
   If no inteval query is given, takes list of intervals from stdin
   */

string genome;

void process_interval(const Interval & t) {
	long nonambig = 0;
	long gc = 0;
	for (int i = t.start; i <= t.end; i++) {
		char c = toupper(genome[i]);
		if (c == 'A' || c == 'T') {
			nonambig++;
		} else if ( c == 'C' || c == 'G') {
			nonambig++;
			gc++;
		}
	}
	cout << t << "\t" << gc / (double) nonambig << "\t" << nonambig / (double) (t.end - t.start + 1) << endl;
}

void usage(char* argv[]) {
	cerr << "Usage: " << argv[0] << " fasta_file [interval] " << endl;
	cerr << "Takes an interval query and outputs the gcContent of that range (%gc out of non-ambig bases, and %non-ambig bases). " << endl;
    cerr << "If no inteval query is given, takes list of intervals from stdin." << endl;
	cerr << "Intervals are closed zero-based." << endl;
}

int main(int argc, char * argv[]) {
	if ((argc != 2) && (argc != 4)) {
		usage(argv);
		exit(1);
	}
	string filename = argv[1];
	uint64_t start = -1;
	uint64_t end = -1;
	if (argc == 4) {
		start = atol(argv[2]);
		end   = atol(argv[3]);
	}
		
	genome = read_genome(filename);

	if (start == -1) {
		while (getline(cin, sbuf)) {
			Interval i = read_interval(sbuf);
			process_interval(i);
		}
	} else {
		process_interval(Interval("-1", start, end));
	}


	return 0;
}


