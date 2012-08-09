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


string sbuf;
const int maxChrLen = 400000000;

void convert(char * masks, int & curPos, string line) {

	for (int i = 0; i < line.size(); i++) {
		if (line[i] == '!') {
			masks[curPos] = '0';
		} else {
			masks[curPos] = '1';
		}
		curPos++;
	}

}

void usage(int argc, char** argv) {
	cout << "Usage: " << argv[0] << " gem_filename output_prefix \n";
	exit(1);
}

int main (int argc, char** argv) {
	if (argc!=3) usage(argc, argv);

	string gem_filename      = argv[1];
	string output_dir        = argv[2];


	ifstream inf;
	open_file(inf, gem_filename);
	char * masks = NULL;
	int curPos;
	ofstream outf;
	bool inHeader = true;

	while (getline(inf, sbuf)) {
		if (sbuf.substr(0,4) == "~chr") { //do startup or cleanup 

			//write masks to file
			if (masks != NULL) {
				outf.write((char *) masks, curPos);
				delete masks;
				outf.close();
			}

			//start new chromosome
			inHeader=false;
			string chr_name = sbuf.substr(1,6);
			string masks_filename = output_dir + chr_name + ".mask";
			open_file_binary(outf, masks_filename);
			masks = new char[maxChrLen];
			curPos = 0;
		} else if (!inHeader) {
			convert(masks, curPos, sbuf);
		}
	}

	//cleanup after last
	if (masks != NULL) {
		outf.write((char *) masks, curPos);
		delete masks;
		outf.close();
	}

}
