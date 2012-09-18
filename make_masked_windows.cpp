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


void usage(int argc, char** argv) {
	cout << "Usage: " << argv[0] << " masks_filename unmasked_size max_total chr_name\n";
	exit(1);
}

int main (int argc, char** argv) {
	if (argc != 6) usage(argc, argv);

	string masks_filename           = argv[1];
	int unmaskedSize                = atoi(argv[2]);
	int maxTotal                    = atoi(argv[3]);
	string chrName                  = argv[4];

	struct stat filestatus;
	stat( masks_filename.c_str(), &filestatus );
	off_t chrLen = filestatus.st_size;

	ifstream inf;
	open_file_binary(inf, masks_filename);
	char * masks = new char[chrLen];
	inf.read(masks, chrLen);
	inf.close();

	int tot_unmasked = 0;
	int startPos = 0;

	for (int pos = 0; pos < chrLen; pos++) {
		if (masks[pos] == '0') tot_unmasked++;
		if (tot_unmasked == unmaskedSize) { //close off window
			if ((pos - startPos) < maxTotal) {
				cout << chrName << "\t" << startPos + 1 << "\t" << pos + 1 << "\n";
			}
			tot_unmasked = 0;
			startPos = pos + 1;
		}
	}
	
	//close last window
	if ((chrLen - startPos) < maxTotal) {
		cout << chrName << "\t" << startPos + 1 << "\t" << chrLen  << "\n";
	}

	delete[] masks;
	return 0;
}
