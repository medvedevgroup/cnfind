#ifndef __DBTYPES_H__
#define __DBTYPES_H__

#include <stdint.h>
#include <fstream>
#include <iostream>

using namespace std;

/* Read mapping struct */
//#define READNAME_LEN    48
//#define CONTIGNAME_LEN  32
//#define EDITSTRING_LEN  32
#define MAX_BAM_READNAME 60
//#define CONTIGNAME_LEN_SMALL 6
#define CONTIGNAME_LEN_SMALL 10

typedef struct {
	uint64_t 	read_id;
	char 		orientation;
	char        ref_name[CONTIGNAME_LEN_SMALL];
	uint64_t    ref_pos;
	double      normodds;
} Map_t;

ostream & operator<< (ostream & out, const Map_t & mapt) {
	out << mapt.read_id << "\t" << mapt.orientation << "\t" << mapt.ref_name << "\t" << mapt.ref_pos << "\t" << mapt.normodds;
	return out;
}

istream & operator>> (istream & in, Map_t & mapt) {
	in >> mapt.read_id >> mapt.orientation >> mapt.ref_name >> mapt.ref_pos >> mapt.normodds;
	return in;
}


typedef struct {
	uint64_t	dist;
	char        ref_name[CONTIGNAME_LEN_SMALL];
	uint64_t    left_pos;
	uint64_t    right_pos;
	uint64_t    read_id;
	uint64_t	type;
	char        l_orientation;
	char        r_orientation;
} Mmap_t;

ostream & operator<< (ostream & out, const Mmap_t & mmap) {
	out << mmap.dist << "\t";
	out << mmap.ref_name << "\t";
	out << mmap.left_pos << "\t";
	out << mmap.right_pos << "\t";
	out << mmap.read_id << "\t";
	out << mmap.type << "\t";
	out << mmap.l_orientation << "\t";
	out << mmap.r_orientation;
	return out;
}

istream & operator>> (istream & in, Mmap_t & mmap) {
	in >> mmap.dist >> mmap.ref_name >> mmap.left_pos >>  mmap.right_pos >> mmap.read_id >> mmap.type >> mmap.l_orientation >>  mmap.r_orientation;
	return in;
}



typedef struct {
	double      normodds;
} coverage_t_small;


#endif /* !defined(__DBTYPES_H__) */
