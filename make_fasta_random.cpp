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



void usage(int argc, char * argv[]) {
	cerr << "Usage: " << argv[0] << " " << endl;
	cerr << "Makes the nucleotides random and writes to stdout\n";
	exit(1);
}


int main(int argc, char * argv[]) {

	if (argc != 1) usage(argc, argv);
	
	//output header
	getline(cin, sbuf);
	cout << sbuf << endl;

	char c;
	while (cin.get(c)) {
		if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
			c = toupper(num2nt(randNumber(0,3)));
		} else if (c == 'a' || c == 'c' || c == 'g' || c == 't') {
			c = tolower(num2nt(randNumber(0,3)));
		} 
		cout << c;
	}

	return 0;
}


	
