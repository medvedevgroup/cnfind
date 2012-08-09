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


	
