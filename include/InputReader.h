#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


#include "dbtypes.h"
#include "defs.h"
#include "../bamtools/BamReader.h"
#include "../bamtools/BamAux.h"

using namespace BamTools;

class InputReader {
public:
	InputReader(string _input_format, string _filename, uint64_t _first_read_id = 0) : input_format(_input_format), filename(_filename), lastReadName("2paclives"), lastReadId(_first_read_id - 1) {
		if (input_format == "bam"){
			bamreader = new BamReader();
			bamreader->Open(filename);
			refVector = bamreader->GetReferenceData();
		} else if (input_format == "txt" || input_format == "bowtie_concise") {
			inf = new ifstream();
			open_file(*inf, filename);
		} else {
			cerr << "InputReader: unknown input_format " << input_format << endl;
			exit(1);
		}
	}

	bool getNext(Map_t & mapt) {
		string sbuf;
		if (input_format == "bam") {
			BamAlignment bammap;
			if (!bamreader->GetNextAlignment(bammap)) return false;
			mapt = Bam2mapt(bammap);
			mapt.ref_pos++; //cnver uses 1-based indexing, while bam uses 0-based
		} else if (input_format == "bowtie_concise") {
			if (!getline(*inf, sbuf)) return false;
			int ref_num;
			int read_in = sscanf(sbuf.c_str(),"%u%c:<%d,%u,%*u>",&(mapt.read_id),&(mapt.orientation),&ref_num,&(mapt.ref_pos));
			assert(read_in==4);
			if (ref_num< 22) { //Bowtie chr index correction
				ref_num++;
			} else if (ref_num == 22) {
				ref_num = 0;
			}
			strcpy(mapt.ref_name, ("chr" + make_string(ref_num)).c_str());
			mapt.ref_pos++; //cnver uses 1-based indexing, while bowtie is 0-based.
		} else if (input_format == "txt"){
			if (!getline(*inf, sbuf)) return false;
			istringstream line(sbuf);
			line >> mapt.read_id >> mapt.orientation >> mapt.ref_name >> mapt.ref_pos;
		}
		return true;
	}

	// checks if no refNames start with "chr"
	bool isRefShort() {
		if (input_format != "bam") {
			cerr << "InputReader: isRefShort function is only supported in bam mode\n";
			exit(1);
		}
		for (int i = 0; i < refVector.size(); i++) {
			if (refVector[i].RefName.substr(0,3) == "chr") return false;
		}
		return true;
	}

	~InputReader() {
		if (input_format == "bam"){
			delete bamreader;
		} else if (input_format == "txt" || input_format == "bowtie_concise") {
			inf->close();
			delete inf;
		}
	}
	string getLastReadName() { return lastReadName; } //this is only supported with bam input

private:
	string input_format;
	string filename;
	BamReader * bamreader;
	ifstream * inf;
	string lastReadName;
	uint64_t lastReadId;
	BamTools::RefVector refVector;

	Map_t Bam2mapt(const BamAlignment & bammap ) {
		Map_t mapt;
		if (lastReadName != bammap.Name) {
			lastReadName = bammap.Name;
			lastReadId++;
		}
		mapt.read_id = lastReadId;
		if (bammap.IsReverseStrand()) {
			mapt.orientation = '-';
		} else {
			mapt.orientation = '+';
		}
		if (bammap.RefID == -1) { //unaligned
			strcpy(mapt.ref_name, "-1");
		} else {
			strcpy(mapt.ref_name, refVector.at(bammap.RefID).RefName.c_str());
		}
		mapt.ref_pos = bammap.Position;
		mapt.normodds = 0;
		return mapt;
	}
};

