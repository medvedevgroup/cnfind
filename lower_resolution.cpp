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
string sbuf;
#include "include/defs.h"
#include "include/interval.h"


vector<Interval> ints;

void usage(int argc, char * argv[]) {
	cerr << "Usage: " << argv[0] << " windows maxgap" << endl;
	cerr << "Takes intervals in stdin, and outputs...\n";
	exit(1);
}

int calc_gap(int first, int last) {
	assert(first <= last);
	if (first == last) return  0;

	int totGap = 0;
	for (int i = first; i < last; i++) {
		totGap += ints[i + 1].start - ints[i].end - 1;
	}
	return totGap;
}



vector<double> add_rows(vector<double> & row1, vector<double> & row2) {
	assert(row1.size() == row2.size());
	vector<double> retval(row1.size());

	for (int i = 0; i < row1.size(); i++) {
		retval[i] = row1[i] + row2[i];
	}
	return retval;
}


vector<double> parse_label(string & label) {
	istringstream lineStream(label);
	string s;
	vector<double> retval;
	while (getline(lineStream, s, '\t')) {
		retval.push_back(atof(s.c_str()));
	}
	return retval;
}



string merge_labels(int first, int last) {
	assert(first <= last);
	if (first == last) return ints[first].label;
	vector<double> row = parse_label(ints[first].label);
	for (int i = first + 1; i <= last; i++) {
		vector<double> row2 = parse_label(ints[i].label);
		row = add_rows(row, row2); 
	}
	ostringstream oStr;
	for (int i = 0; i < row.size() - 1; i++) {
		oStr << make_string(row[i]) << "\t";
	}
	oStr << make_string(row[row.size() - 1]);
	return oStr.str();
}


int main(int argc, char * argv[]) {

	if (argc != 3) usage(argc, argv);
	int window = atoi(argv[1]);
	int maxGap = atoi(argv[2]);

	read_intervals(cin, ints);

	//we assume its all on the same chromosome

	assert (ints.size() > 0);

	for (int i = 0; i < ints.size() - window + 1; i++) {
		int gap = calc_gap(i, i + window - 1);
		if (gap < maxGap) {
			string joint_label = merge_labels(i, i + window - 1);
			cout << ints[i].chr << "\t" << ints[i].start << "\t" << ints[i+window-1].end << "\t" << joint_label << endl;
		}
	}

	//forget about the last window (on purpose)

}

