
#ifndef INTERVAL_H 
#define INTERVAL_H


#include<cassert>
#include<cmath>
#include<cstring>
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


class Link {
public:
	string chr;
	int from;
	int to;
	int type;
	string label;

	Link() {}
	Link(string _chr, int _from, int _to, int _type, string _label = "") : chr(_chr), from(_from), to(_to), type(_type), label(_label)  {}

};

ostream & operator << (ostream & out, const Link & i) {
	out << i.chr << "\t" << i.from << "\t" << i.to << "\t" << i.type;
	if (i.label != "") out << "\t" << i.label;
	return out;
}

Link read_link(string sbuf) {
	istringstream line(sbuf);
	Link i;
	line >> i.chr >> i.from >> i.to >> i.type;
	getline(line, i.label); //the rest of the line is a label
	size_t pos = i.label.find_first_not_of(" \t"); //trim leading whitespace
	if (pos != string::npos) i.label = i.label.substr(pos);
	return i;
}


template<class T>
int read_links(istream & in, T & links) {
	string sbuf;
	int num_read = 0;
	while (getline(in, sbuf)) {
		links.push_back(read_link(sbuf));
		num_read++;
	}
	return num_read;
}

class Interval {  //the intervals are assumed to contain the endpoints.
	public:
		string chr;
		int start;
		int end;
		string label;

		Interval() {}
		Interval(string _chr, int _start, int _end, string _label = "") : chr(_chr), start(_start), end(_end), label(_label)  {}

		bool operator< (const Interval & i) const {
			if (chr == i.chr) return start < i.start;
			return chr < i.chr;
		}
		/*bool operator<= (const Interval & i) const {
		  if (chr == i.chr) return start <= i.start;
		  return chr < i.chr;
		  }*/

		bool contains(const Interval & c1) const {
			if (c1.chr != chr) return false;
			if ((start <= c1.start) && (c1.end <= end)) return true;
			return false;
		}

		bool contains(const string _chr, const int point) {
			if (_chr != chr) return false;
			if ((start <= point) && (point <= end)) return true;
			return false;
		}

		bool overlaps(const Interval &c) { // note containment counts as overlaping
			if (c.chr != chr) return false;
			if ((c.start <= start) && (c.end >= start)) return true;
			if ((c.start <= end)   && (c.end >= end)) return true;
			if ((start <= c.start) && (c.end <= end)) return true;
			return false;
		}
		int amountThatOverlaps(const Interval &c) const {
			if (c.chr != chr) return 0;
			if ((c.start <= start) && (start <= c.end)) return min(end - start + 1, c.end - start + 1);
			if ((c.start <= end)   && (end <= c.end)) return min(end - start + 1, end - c.start + 1);
			if ((start <= c.start) && (c.end <= end)) return c.end - c.start + 1;
			return 0;
		}

		bool fullyLeftOf(const Interval &c) {
			if (chr != c.chr) return (chr < c.chr);
			if (end < c.start) return true;
			return false;
		}
		bool fullyRightOf(const Interval &c) {
			if (chr != c.chr) return (chr > c.chr);
			if (c.end < start) return true;
			return false;
		}

};

bool comp_le (const Interval &i1, const Interval &i2) {
	if (i1.chr == i2.chr) return i1.start <= i2.start;
	return i1.chr < i2.chr;
}
/*bool comp_le_by_right (const Interval &i1, const Interval &i2) {
  if (i1.chr == i2.chr) return i1.end <= i2.end;
  return i1.chr < i2.chr;
  }
  bool fullyLeftOf(const Interval &i1, const Interval &i2) {
  if (i1.chr != i2.chr) return i1.chr < i2.chr;
  return i1.end < i2.start;
  }
 */

ostream & operator << (ostream & out, const Interval & i) {
	//out << i.chr << " " << i.start << " " << i.end;
	//if (i.label != "") out << " " << i.label;
	out << i.chr << "\t" << i.start << "\t" << i.end;
	if (i.label != "") out <<  "\t" << i.label;
	return out;
}

Interval read_interval(string sbuf) {
	istringstream line(sbuf);
	Interval i;
	line >> i.chr >> i.start >> i.end;
	getline(line, i.label); //the rest of the line is a label
	//trim leading whitespace in label
	size_t pos = i.label.find_first_not_of(" \t");
	if (pos != string::npos) i.label = i.label.substr(pos);
	return i;
}


template<class T>
int read_intervals(istream & in, T & intervals) {
	string sbuf;
	int num_read = 0;
	while (getline(in, sbuf)) {
		intervals.push_back(read_interval(sbuf));
		num_read++;
	}
	return num_read;
}


//range is a closed interval
//points should be sorted
//this ignores the chromosome for now
void searchRange(vector<int> & points, Interval range, int & start, int & end) {
	vector<int>::iterator it = lower_bound(points.begin(), points.end(), range.start);
	if (it == points.end() || *it > range.end) {
		start = 0;
		end = 0;
		return;
	}
	start = it - points.begin();
	end = start + 1;
	it++;
	while (it != points.end() && *it <= range.end) {
		it++;
		end++;
	}
}

#endif
