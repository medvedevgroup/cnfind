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

#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


#include "include/dbtypes.h"
#include "include/defs.h"
#include "include/InputReader.h"
//#include "bamtools/BamReader.h" #include "bamtools/BamAux.h"

//using namespace BamTools;
string sbuf;
vector<string> ref_names;
typedef pair<string,int> ref_pair;
typedef map<string,int> ref_map_t;
ref_map_t ref_names_map;
vector<ofstream *> outf;
ofstream readNamesIdx;
string short_chr;

// append every mapping with the norm odds, modified chromosome name, and dump it
void process_block(deque<Map_t> & block) {
	if (block.size() == 0 || 0 == strcmp(block[0].ref_name, "-1")) return;

	//note normodds calculation uses all of the mappings, not just to the current chromosomes.
	double normodds = 1.0 / block.size(); 
	for (int i = 0; i < block.size(); i++) {
		block[i].normodds = normodds;
		string ref_name = block[i].ref_name;
		if (short_chr == "short_chr") {
			ref_name = "chr" + ref_name;
		}
		ref_map_t::iterator it = ref_names_map.find(ref_name); //Identify reference index
		if (it != ref_names_map.end()) {
			int file_index = it->second; 
			outf[file_index]->write((char *) & block[i], sizeof(block[i]));
		}
	}
	return;
}

void usage( int argc, char ** argv) {
	printf("Usage: %s map_list ref_names input_format \n", argv[0]);
	cout << "\tmap_list     : file with names of mapping files.\n";
	cout << "\tref_names    : file with names of reference sequences.\n";
	cout << "\tinput_format : format of mapping files, either \"bam\" or \"txt\".\n";
	cout << "\tshort_chr    : either \"short_chr\" or something else.\n";
	cout << endl;
	exit(1);
}

int main ( int argc, char ** argv) {
	if (argc != 5) usage(argc, argv);
	vector<string> map_list;
	read_file_into_vector(argv[1], map_list);
	read_file_into_vector(argv[2], ref_names);
	string input_format = argv[3];
	short_chr    = argv[4];
	long read_id = 0;

	//open outputfiles and initialize ref_names_map
	outf.resize(ref_names.size());
	for (int i = 0; i < ref_names.size(); i++)  {
		ref_names_map.insert(make_pair(ref_names[i], i));
		outf[i] = new ofstream;
		open_file_binary(*outf[i], ref_names[i]  + ".rmap");
	}

	open_file_binary(readNamesIdx, "readNames.idx");

	//read in the read mappings, grouping them into blocks with the same read_id.  
	//Each block is then handled in process_block
	uint64_t lastReadId = -1;
	for (int map_file_index = 0; map_file_index < map_list.size(); map_file_index++) {
		InputReader reader(input_format, map_list[map_file_index], lastReadId + 1);
		Map_t mapt;
		deque<Map_t> block;
		while (reader.getNext(mapt)) {
			if (lastReadId != mapt.read_id) {
				string bam_readname = reader.getLastReadName();
				assert(bam_readname.size() < MAX_BAM_READNAME);
				char name[MAX_BAM_READNAME];
				strcpy(name, bam_readname.c_str());
				readNamesIdx.write(bam_readname.c_str(), MAX_BAM_READNAME);
				process_block(block);  
				block.clear();
			} 
			lastReadId = mapt.read_id;
			block.push_back(mapt);
		}
		process_block(block);
		block.clear();
	}

	//close files
	for (int i = 0; i < outf.size(); i++) {
		outf[i]->close();
		delete outf[i];
	}

	readNamesIdx.close();

	return 0;	
}
