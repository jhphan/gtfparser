#include <stdio.h>

#include <iostream>
#include <vector>
#include <map>
#include <bitset>

#include <miblab/commandline.h>
#include <miblab/stringtokenizer.h>
#include <miblab/strconv.h>

#define MAX_BUF_LEN 1048576

using namespace std;

int read_line(FILE* f, char* ch_line) {
	int i_size = 0;
	char ch;

	ch = fgetc(f);
	while ((ch != '\n') && (ch != EOF) && (i_size<MAX_BUF_LEN-1)) {
		ch_line[i_size] = ch;
		ch = fgetc(f);
		i_size++;
	}
	ch_line[i_size] = 0;
	if (ch == EOF) return 0;
	return i_size+1;
}

void usage() {
	cout << "filter_duplicates" << endl;
	cout << "\t-sam [sam/bam file]" << endl;
}

int get_pair_flag(int i_flag, map<int,int>& map_pairs) {
	if (i_flag > 256) {
		map<int,int>::iterator it = map_pairs.find(i_flag-256);
		if (it == map_pairs.end())
			return 0;	//  unpaired secondary
		return it->second+256;
	}
	map<int,int>::iterator it = map_pairs.find(i_flag);
	if (it == map_pairs.end())
		return 0;	// unpaired primary
	return it->second;
}

bool filter_sam_pairs(vector<string>& vec_sam_lines, vector<vector<string> >& vec_sam_lines_filtered) {

	map<int,int> map_pairs;
	map_pairs[83] = 163;	map_pairs[163] = 83;
	map_pairs[99] = 147;	map_pairs[147] = 99;
	map_pairs[67] = 131;	map_pairs[131] = 67;
	map_pairs[115] = 179;	map_pairs[179] = 115;
	map_pairs[81] = 161;	map_pairs[161] = 81;
	map_pairs[97] = 145;	map_pairs[145] = 97;
	map_pairs[65] = 129;	map_pairs[129] = 65;
	map_pairs[113] = 177;	map_pairs[177] = 113;

	int i_num_lines = vec_sam_lines.size();

	map<string,string> map_reverse_key;
	map<string,map<int,vector<string> > > map_sam;
	
	for (int i=0; i<i_num_lines; i++) {
		vector<string> vec_sam = StringTokenizer(vec_sam_lines[i],'\t').split();
		string str_end1_seq = vec_sam[2];
		string str_end2_seq;
		if (vec_sam[6] == "=") 
			str_end2_seq = vec_sam[2];
		else
			str_end2_seq = vec_sam[6];
		string str_end1_coord = vec_sam[3];
		string str_end2_coord = vec_sam[7];

		string str_coord_key = str_end1_seq+string(":")+str_end1_coord+string("-")+str_end2_seq+string(":")+str_end2_coord;
		string str_match_key = str_end2_seq+string(":")+str_end2_coord+string("-")+str_end1_seq+string(":")+str_end1_coord;
		int i_flag = atoi(vec_sam[1].c_str());

		bool b_reverse = false;
		if (str_coord_key == str_match_key) {	// duplicate reads in the pair
			if (i_flag < 256) {
				if (i_flag > 128) b_reverse = true;
			} else {
				if (i_flag > 384) b_reverse = true;
			}
			if (b_reverse) {
				str_coord_key = str_coord_key+string("r");
			} else {
				str_match_key = str_match_key+string("r");
			}
		}

		map_reverse_key[str_coord_key] = str_match_key;
		map_sam[str_coord_key][i_flag] = vec_sam;

		//cout << vec_sam_lines[i] << endl;
	}

	// make sure all mappings are properly paired
	/*
	map<string,map<int,string> >::iterator it;
	for (it=map_sam.begin(); it!=map_sam.end(); it++) {
		map<int,string>::iterator it2;
		for (it2=it->second.begin(); it2!=it->second.end(); it2++) {
			map<int,string>::iterator it2_search;
			map<int,string> map_reverse = map_sam[map_reverse_key[it->first]];
			it2_search = map_reverse.find(get_pair_flag(it2->first));
			if (it2_search == map_reverse.end()) {
				cerr << "!error: filter_duplicates: filter_sam_pairs(): incorrent sam pairing" << endl;

				for (int i=0; i<i_num_lines; i++) {
					cerr << vec_sam_lines[i] << endl;
				}

				return false;
			}
		}
	}
	*/

	// remove duplicates
	map<string,map<int,vector<string> > >::iterator it;
	for (it=map_sam.begin(); it!=map_sam.end(); it++) {
		vector<int> vec_flags;
		int i_smallest = 500;
		map<int,vector<string> >::iterator it2;
		for (it2=it->second.begin(); it2!=it->second.end(); it2++) {
			vec_flags.push_back(it2->first);
			if (i_smallest > it2->first) i_smallest = it2->first;
		}
		for (int i=0; i<vec_flags.size(); i++) {
			if (vec_flags[i] > i_smallest) { // not the primary alignment
				// remove its pair, if it exists
				map_sam[map_reverse_key[it->first]].erase(get_pair_flag(vec_flags[i], map_pairs));
				// remove it
				it->second.erase(vec_flags[i]);
			}
		}
	}
	
	// copy remaining lines to new list
	vec_sam_lines_filtered = vector<vector<string> >();
	for (it=map_sam.begin(); it!=map_sam.end(); it++) {
		map<int,vector<string> >::iterator it2;
		for (it2=it->second.begin(); it2!=it->second.end(); it2++) {
			if (it2->second.size() > 0) {
				// find pair if it exists
				int i_pair_flag = get_pair_flag(it2->first, map_pairs);
				if (i_pair_flag > 0) {
					// read should be paired
					string str_reverse_key = map_reverse_key[it->first];
					map<int,vector<string> >::iterator it2_search;
					it2_search = map_sam[str_reverse_key].find(i_pair_flag);
					if (it2_search == map_sam[str_reverse_key].end()) {
						// other end of pair not found
						string str_seq1 = it2->second[2];
						string str_seq2 = it2->second[6];
						int i_new_flag = it2->first;
						bitset<8*sizeof(int)> bit(i_new_flag);
						if ( ((str_seq1 != str_seq2) && (str_seq2 != "=")) || (!bit[1]) ) {
							// change flag to singleton
							//cerr << "flag: " << i_new_flag << endl;
							//cerr << "bits: " << bit << endl;
							if (bit[5]) {
								i_new_flag -= 32;
								//cerr << "new flag: " << i_new_flag << endl;
							}
							i_new_flag += 8;
							//cerr << "final flag: " << i_new_flag << endl;
							it2->second[1] = conv2str(i_new_flag); 
						}
						
						vec_sam_lines_filtered.push_back(it2->second);
						it2->second.clear();
					} else {
						// print pair
						vec_sam_lines_filtered.push_back(it2->second);
						it2->second.clear();
						vec_sam_lines_filtered.push_back(it2_search->second);
						it2_search->second.clear();
					}
				} else {
					vec_sam_lines_filtered.push_back(it2->second);
					it2->second.clear();
				}
			}
		}
	}
	
	return true;
}

bool filter_duplicates(string str_sam_file) {

	FILE* file_sam;
	string str_cmd = string("java -jar /data/ngs/bin/picard/ViewSam.jar I=")+str_sam_file+string(" VALIDATION_STRINGENCY=SILENT");
	file_sam = popen(str_cmd.c_str(), "r");
	if (file_sam == NULL) {
		cerr << "!error: filter_duplicates: filter_duplicates(): cannot open: " << str_cmd << endl;
		return false;
	}


	// current list of all reads
	vector<string> vec_sam_lines = vector<string>();
	string str_cur_read;

	char ch_buf[MAX_BUF_LEN];
	while (read_line(file_sam, ch_buf)) {
		if (ch_buf[0] == '@') {
			// header line
			cout << ch_buf << endl;
		} else {
			string str_sam = string(ch_buf);
			string str_read_name;
			StringTokenizer(str_sam, '\t').getNextToken(str_read_name);
			if (str_read_name != str_cur_read) {
				if (vec_sam_lines.size() > 0) {
					// process this list for duplicates
					vector<vector<string> > vec_sam_lines_filtered;
					if (!filter_sam_pairs(vec_sam_lines, vec_sam_lines_filtered)) {
						cerr << "!error: filter_duplicates: filter_duplicates(): cannot filter sam pairs" << endl;
						return false;
					}
					// print list
					for (int i=0; i<vec_sam_lines_filtered.size(); i++) {
						for (int j=0; j<vec_sam_lines_filtered[i].size(); j++) {
							if (j > 0) cout << "\t";
							cout << vec_sam_lines_filtered[i][j];
						}
						cout << endl;
					}
				}
				// new read
				str_cur_read = str_read_name;
				// reset list
				vec_sam_lines = vector<string>(1);
				vec_sam_lines[0] = str_sam;
			} else {
				// current line = current read
				vec_sam_lines.push_back(str_sam);
			}
		}
	}

	// reached end of file
	if (vec_sam_lines.size() > 0) {
		// process this list for duplicates
		vector<vector<string> > vec_sam_lines_filtered;
		if (!filter_sam_pairs(vec_sam_lines, vec_sam_lines_filtered)) {
			cerr << "!error: filter_duplicates: filter_duplicates(): cannot filter sam pairs" << endl;
			return false;
		}
		// print list
		for (int i=0; i<vec_sam_lines_filtered.size(); i++) {
			for (int j=0; j<vec_sam_lines_filtered[i].size(); j++) {
				if (j > 0) cout << "\t";
				cout << vec_sam_lines_filtered[i][j];
			}
			cout << endl;
		}
	}

	fclose(file_sam);

	return true;
}


int main(int argc, char* argv[]) {

	string str_sam_file;
	CommandLine cl(argc, argv);
	if (!cl.getArg(str_sam_file, "sam")) {
		usage();
		return 1;
	}

	// read sam file
	if (!filter_duplicates(str_sam_file)) {
		cerr << "!error: filter_duplicates: cannot filter file" << endl;
		return 1;
	}

	return 0;
}
