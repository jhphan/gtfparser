#include "../include/sam.h"

bool split_cigar_string(string str_cigar, vector<int>& vec_cigar_numbers, vector<char>& vec_cigar_codes) {

	int i_length = str_cigar.length();


	int i_cur_pos = 0;
	int i_start_pos = 0;
	string str_substr;
	while (i_cur_pos < i_length) {
		char ch_cur = str_cigar[i_cur_pos];
		switch (ch_cur) {
			case 'M':
			case 'I':
			case 'D':
			case 'N':
			case 'S':
			case 'H':
			case 'P':
			case '=':
			case 'X': 
				//
				if (i_cur_pos == i_start_pos) {
					// no number prior to code
					// this shouldn't happen
					cerr << "!error: split_cigar_string(): missing cigar number: " << str_cigar << endl;
					return false;
				}
				str_substr = str_cigar.substr(i_start_pos, i_cur_pos-i_start_pos);
				int i_cigar_number = atoi(str_substr.c_str());
				vec_cigar_numbers.push_back(i_cigar_number);
				vec_cigar_codes.push_back(ch_cur);
				i_start_pos = i_cur_pos+1;
				break;
		}
		i_cur_pos++;
	}

	return true;
}

bool translate_cigar_string(string str_cigar, Transcript& trans, unsigned long ul_trans_coord, int i_exon_number, string& str_new_cigar, int& i_end_tlen) {
	if (i_exon_number >= trans.get_num_exons()) {
		cerr << "!error: translate_cigar_string(): exon number out of range: " << i_exon_number << endl;
		return false;
	}

	// split cigar into pieces
	vector<int> vec_cigar_numbers;
	vector<char> vec_cigar_codes;
	if (!split_cigar_string(str_cigar, vec_cigar_numbers, vec_cigar_codes)) {
		cerr << "!error: translate_cigar_string(): cannot split cigar string: " << str_cigar << endl;
		return false;
	}
	
	int i_num_cigar = vec_cigar_numbers.size();
	str_new_cigar = "";
	unsigned long ul_cur_pos = ul_trans_coord;
	i_end_tlen = 0;

	bool b_reached_end = false;

	for (int i=0; i<i_num_cigar; i++) {

		if (b_reached_end) {
			
		}

		switch (vec_cigar_codes[i]) {
			case 'S':
			case 'H': {
				// ignore hard or soft clipped regions (don't change position
				str_new_cigar += conv2str(vec_cigar_numbers[i]);
				str_new_cigar += vec_cigar_codes[i];
				break;
			}
			case 'M':
			case 'D': {
				// introns can be inserted into the middle of M's (matches/mismatches) or D's (deletions)
				if (ul_cur_pos > trans.get_exon_rel_end(i_exon_number)) {
					// reached the end of the exon, try to go to the next
					if (i_exon_number+1 > trans.get_num_exons()) {
						b_reached_end = true;
						//cerr << "!error: translate_cigar_string(): read maps to position outside of known exon regions (M,D)" << endl;
						//cerr << trans << endl;
						//return false;
					}
					// add intron
					i_exon_number++;
					int i_intron_size = trans.get_exon_start(i_exon_number)-trans.get_exon_end(i_exon_number-1)-1;
					str_new_cigar += conv2str(i_intron_size);
					str_new_cigar += 'N';
					i_end_tlen += i_intron_size;
				}

				int i_cigar_number_remaining = vec_cigar_numbers[i];
				while (i_cigar_number_remaining > 0) {
					if (ul_cur_pos+i_cigar_number_remaining-1>trans.get_exon_rel_end(i_exon_number)) {
						// 'M' or 'D' cross exon junction
						//cerr << "exon number: " << i_exon_number << endl;
						//cerr << "max exon number: " << trans.get_num_exons() << endl;
						//cerr << "exon rel end: " << trans.get_exon_rel_end(i_exon_number) << endl;
						//cerr << "ul_cur_pos: " << ul_cur_pos << endl;
						int i_num_M = trans.get_exon_rel_end(i_exon_number)-ul_cur_pos+1;
						//cerr << "M: " << i_num_M << endl;
						i_cigar_number_remaining -= i_num_M;
						//cerr << "cigar num remain: " << i_cigar_number_remaining << endl;
						ul_cur_pos += i_num_M;
						if (i_exon_number+1 >= trans.get_num_exons()) {
							b_reached_end = true;
							//cerr << "!error: translate_cigar_string(): read maps to position outside of known exon regions (M,D)" << endl;
							//cerr << trans << endl;
							//return false;
							
							// force a softclip
							//cerr << "!warning: translate_cigar_string(): read maps to position outside of known exon regions (M,D), forcing a soft clip" << endl;
							//cerr << trans << endl;

							//cerr << "M: " << i_num_M << endl;
							//cerr << "cigar num remain: " << i_cigar_number_remaining << endl;

							str_new_cigar += conv2str(i_num_M);
							//cerr << str_new_cigar << endl;
							str_new_cigar += vec_cigar_codes[i];
							//cerr << str_new_cigar << endl;
							i_end_tlen += i_num_M;
							str_new_cigar += conv2str(i_cigar_number_remaining);
							//cerr << str_new_cigar << endl;
							str_new_cigar += string("S");
							//cerr << str_new_cigar << endl;
							i_cigar_number_remaining = 0;
						} else {
							str_new_cigar += conv2str(i_num_M);
							str_new_cigar += vec_cigar_codes[i];
							i_end_tlen += i_num_M;
							i_exon_number++;
							int i_intron_size = trans.get_exon_start(i_exon_number)-trans.get_exon_end(i_exon_number-1)-1;
							str_new_cigar += conv2str(i_intron_size);
							str_new_cigar += 'N';
							i_end_tlen += i_intron_size;
						}
					} else {
						// 'M' or 'D' doesn't cross exon junction
						str_new_cigar += conv2str(i_cigar_number_remaining);
						str_new_cigar += vec_cigar_codes[i];
						i_end_tlen += i_cigar_number_remaining;
						ul_cur_pos += i_cigar_number_remaining;
						i_cigar_number_remaining = 0;
					}
				}
				break;
			}
			case 'I': {
				if (ul_cur_pos > trans.get_exon_rel_end(i_exon_number)) {
					// if 'I' starts at exon boundary, it comes before the intron
					str_new_cigar += conv2str(vec_cigar_numbers[i]);
					str_new_cigar += vec_cigar_codes[i];
					if (i_exon_number+1 > trans.get_num_exons()) {
						cerr << "!error: translate_cigar_string(): read maps to position outside of known exon regions (I)" << endl;
						cerr << trans << endl;
						return false;
					}
					i_exon_number++;
					int i_intron_size = trans.get_exon_start(i_exon_number)-trans.get_exon_end(i_exon_number-1)-1;
					str_new_cigar += conv2str(i_intron_size);
					str_new_cigar += 'N';
					i_end_tlen += i_intron_size;
				} else {
					str_new_cigar += conv2str(vec_cigar_numbers[i]);
					str_new_cigar += vec_cigar_codes[i];
					i_end_tlen -= vec_cigar_numbers[i];
				}
				break;
			}
		}
	}

	return true;
}

void parse_sam_record(char* ch_buf, struct sam_record& sr) {
	sr.vec_sam = StringTokenizer(ch_buf, '\t').split();

	sr.str_read_name = sr.vec_sam[0];
	sr.i_flag = atoi(sr.vec_sam[1].c_str());
	sr.str_ref_name = sr.vec_sam[2];
	sr.ul_map_pos = atol(sr.vec_sam[3].c_str());
	sr.i_map_qual = atoi(sr.vec_sam[4].c_str());
	sr.str_cigar = sr.vec_sam[5];
	sr.str_mate_ref_name = sr.vec_sam[6];
	sr.ul_mate_map_pos = atol(sr.vec_sam[7].c_str());
	sr.i_tlen = atoi(sr.vec_sam[8].c_str());
	sr.i_end_tlen = 0;
	sr.i_strand = 1;
}

void print_sam_record(struct sam_record& sr) {
	cout << sr.str_read_name;
	cout << "\t" << sr.i_flag;
	cout << "\t" << sr.str_ref_name;
	cout << "\t" << sr.ul_map_pos;
	cout << "\t" << sr.i_map_qual;
	cout << "\t" << sr.str_cigar;
	cout << "\t" << sr.str_mate_ref_name;
	cout << "\t" << sr.ul_mate_map_pos;
	cout << "\t" << sr.i_tlen;
	for (int i=9; i<sr.vec_sam.size(); i++)
		cout << "\t" << sr.vec_sam[i];
	cout << "\t" << "XS:A:";
	if (sr.i_strand > 0) cout << "+"; else cout << "-";
	cout << endl;
}

bool translate_sam_record(struct sam_record& sr, map<string,Transcript>& map_transcript_structure) {
	// check if transcript exists
	map<string,Transcript>::iterator it_trans;
	it_trans = map_transcript_structure.find(sr.str_ref_name);
	if (it_trans == map_transcript_structure.end()) {
		cerr << "!error: translate_sam_record(): transcript not found: " << sr.str_ref_name << endl;
		return false;
	}

	// translate map position
	unsigned long ul_genomic_map_pos = 0;
	int i_exon_number = 0;
	if (!it_trans->second.translate_to_genomic_coord(sr.ul_map_pos, i_exon_number, ul_genomic_map_pos)) {
		cerr << "!error: translate_sam_record(): cannot translate to genomic coord: " << sr.ul_map_pos << endl;
		return false;
	}

	// translate cigar string
	string str_new_cigar;
	int i_end_tlen = 0;
	if (!translate_cigar_string(sr.str_cigar, it_trans->second, sr.ul_map_pos, i_exon_number, str_new_cigar, i_end_tlen)) {
		cerr << "!error: translate_sam_record(): cannot translate cigar string: " << sr.str_cigar << ", " << sr.str_read_name << endl;
 		return false;
	}

	// check mate information
	unsigned long ul_mate_genomic_map_pos;
	string str_mate_ref_name = sr.str_mate_ref_name;
	if (sr.str_mate_ref_name == "*") {
		// mate not mapped
		ul_mate_genomic_map_pos = 0;
	} else if (sr.str_mate_ref_name == "=") {
		// mate mapped to same ref
		if (!it_trans->second.translate_to_genomic_coord(sr.ul_mate_map_pos, i_exon_number, ul_mate_genomic_map_pos)) {
			cerr << "!error: translate_sam_record(): cannot translate to genomic coord: " << sr.ul_mate_map_pos << endl;
			return false;
		}
	} else {
		// mate mapped to a different ref
		map<string,Transcript>::iterator it_mate_trans;
		it_mate_trans = map_transcript_structure.find(sr.str_mate_ref_name);
		if (it_mate_trans == map_transcript_structure.end()) {
			cerr << "!error: translate_sam_record(): mate transcript not found: " << sr.str_mate_ref_name << endl;
			return false;
		}

		str_mate_ref_name = it_mate_trans->second.get_chromosome();
		if (!it_mate_trans->second.translate_to_genomic_coord(sr.ul_mate_map_pos, i_exon_number, ul_mate_genomic_map_pos)) {
			cerr << "!error: translate_sam_record(): cannot translate mate to genomic coord: " << sr.ul_mate_map_pos << endl;
			return false;
		}
	}

	// update record
	sr.str_ref_name = it_trans->second.get_chromosome();
	sr.ul_map_pos = ul_genomic_map_pos;
	sr.str_cigar = str_new_cigar;
	sr.str_mate_ref_name = str_mate_ref_name;
	sr.ul_mate_map_pos = ul_mate_genomic_map_pos;
	sr.i_end_tlen = i_end_tlen;
	sr.i_strand = it_trans->second.get_strand();
	
	return true;
}

void fix_template_length(struct sam_record& sr1, struct sam_record& sr2) {
	if (sr1.ul_map_pos < sr2.ul_map_pos) {
		sr1.i_tlen = sr2.ul_map_pos + sr2.i_end_tlen - sr1.ul_map_pos;
		sr2.i_tlen = -sr1.i_tlen;
	} else {
		sr2.i_tlen = sr1.ul_map_pos + sr1.i_end_tlen - sr2.ul_map_pos;
		sr1.i_tlen = -sr2.i_tlen;
	}
}

