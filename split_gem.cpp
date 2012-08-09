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


