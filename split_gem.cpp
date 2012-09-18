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
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<list>
#include<queue>
#include<cstdarg>
#include<algorithm>

using namespace std;
ostringstream oMsg;
string sbuf;
#include "include/defs.h"

/***********************************/
/* this program converts a gem mappability file into one with stricter criteria (unique I think, but you should check the code */
/***********************************/

void usage(int argc, char * argv[]) {
	cerr << "Usage: " << argv[0] << " <db> -g genome -k kmersize -o outfile" << endl;
	cerr << "Program descrption.\n";
	exit(1);
}


void replace(string & sbuf) {
	for (int i = 0; i < sbuf.length(); i++) {
		if (sbuf[i] != '!') sbuf[i] = ' ';
	}
}

int main(int argc, char * argv[]) {

	ifstream inf;
	open_file(inf, argv[1]);
	string target = argv[2];
	bool skipping = false;
	bool inHeader = true;
	while (getline(inf, sbuf)) {
		if (sbuf.substr(0,4) == "~chr") { 
			inHeader=false;
			if (sbuf.substr(1,6) == target) {
				skipping = false;
			} else {
				skipping = true;
			}
		} 
		if (!skipping) {
			if (!inHeader) replace(sbuf);
			cout << sbuf << endl;
		}
	}

	
	/*
	char ch;
	while ((ch = getopt(argc, argv, "g:k:o:")) != -1) {
		switch (ch) {
			case 'g':
				genome = read_genome(optarg);
				break;
			case 'k':
				kmersize = atoi(optarg);
				break;
			case 'o':
				base = optarg;
				break;
		}
	}
	*/
	return 0;
}


