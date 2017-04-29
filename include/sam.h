#ifndef __SAM_H__
#define __SAM_H__

using namespace std;

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

#include <miblab/stringtokenizer.h>
#include <miblab/strconv.h>

#include "./transcript.h"

struct sam_record {
	string str_read_name;
	int i_flag;
	string str_ref_name;
	unsigned long ul_map_pos;
	int i_map_qual;
	string str_cigar;
	string str_mate_ref_name;
	unsigned long ul_mate_map_pos;
	int i_tlen;	// total template length of pair
	int i_end_tlen;	// template length of one end, if paired
	int i_strand;

	vector<string> vec_sam;	// all records
};

/*
struct reference_coords {
	string str_chromosome;
	unsigned long ul_start_coord;
	unsigned long ul_end_coord;
	string str_transcript_id;
};
*/

bool split_cigar_string(string str_cigar, vector<int>& vec_cigar_numbers, vector<char>& vec_cigar_codes);
bool translate_cigar_string(string str_cigar, Transcript& trans, unsigned long ul_trans_coord, int i_exon_number, string& str_new_cigar, int& i_end_tlen);
void parse_sam_record(char* ch_buf, struct sam_record& sr);
void print_sam_record(struct sam_record& sr);
bool translate_sam_record(struct sam_record& sr, map<string,Transcript>& map_transcript_structure);
void fix_template_length(struct sam_record& sr1, struct sam_record& sr2);

#endif
