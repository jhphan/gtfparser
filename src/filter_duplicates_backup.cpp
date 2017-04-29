#include <stdio.h>

#include <iostream>
#include <vector>
#include <map>

#include <miblab/commandline.h>
#include <miblab/stringtokenizer.h>

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

int get_pair_flag(int i_flag) {
	if (i_flag == 83) return 163;
	if (i_flag == 163) return 83;
	if (i_flag == 99) return 147;
	if (i_flag == 147) return 99;
	if (i_flag == 355) return 403;
	if (i_flag == 403) return 355;
	if (i_flag == 339) return 419;
	if (i_flag == 419) return 339;
	return 0;
}

bool filter_sam_pairs(vector<string>& vec_sam_lines, vector<string>& vec_sam_lines_filtered) {

	int i_num_lines = vec_sam_lines.size();

	map<string,string> map_reverse_key;
	map<string,map<int,string> > map_sam;
	
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

		map_reverse_key[str_coord_key] = str_match_key;
		map_sam[str_coord_key][i_flag] = vec_sam_lines[i];
	}

	// make sure all mappings are properly paired
	map<string,map<int,string> >::iterator it;
	for (it=map_sam.begin(); it!=map_sam.end(); it++) {
		map<int,string>::iterator it2;
		for (it2=it->second.begin(); it2!=it->second.end(); it2++) {
			map<int,string>::iterator it2_search;
			map_reverse = map_sam[map_reverse_key[it->first]];
			it2_search = map_reverse.find(get_pair_flag(it2->first));
			if (it2_search == map_reverse.end()) {
				cerr << "!error: filter_duplicates: filter_sam_pairs(): incorrent sam pairing" << endl;
				return false;
			}
		}
	}

	// remove duplicates
	map<string,map<int,string> >::iterator it;
	for (it=map_sam.begin(); it!=map_sam.end(); it++) {
		vector<int> vec_flags;
		int i_smallest = 500;
		map<int,string>::iterator it2;
		for (it2=it->second.begin(); it2!=it->second.end(); it2++) {
			vec_flags.push_back(it2->first);
			if (i_smallest > it2->first) i_smallest = it2->first;
		}
		for (int i=0; i<vec_flags.size(); i++) {
			if (vec_flags[i] > i_smallest) { // not the primary alignment
				// remove its pair
				map_sam[map_reverse_key[it->first]].erase(get_pair_flag(vec_flags[i]));
				// remove it
				it->second.erase(vec_flags[i]);
			}
		}
	}
	
	// copy remaining lines to new list
	vec_sam_lines_filtered = vector<string>();
	map<string,map<int,string> >::iterator it;
	for (it=map_sam.begin(); it!=map_sam.end(); it++) {
		map<int,string>::iterator it2;
		for (it2=it->second.begin(); it2!=it->second.end(); it2++) {
			vec_sam_lines_filtered.push_back(it2->second);
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
					vector<string> vec_sam_lines_filtered;
					if (!filter_sam_pairs(vec_sam_lines, vec_sam_lines_filtered)) {
						cerr << "!error: filter_duplicates: filter_duplicates(): cannot filter sam pairs" << endl;
						return false;
					}
					// print list
					for (int i=0; i<vec_sam_lines_filtered.size(); i++)
						cout << vec_sam_lines_filtered[i] << endl;
				} else {
					// empty list
				}
				// new read
				str_cur_read = str_read_name;
				// reset list
				vec_sam_lines = vector<string>();
			} else {
				// current line = current read
				vec_sam_lines.push_back(str_sam);
			}
		}
	}

	// reached end of file
	if (vec_sam_lines.size() > 0) {
		// process this list for duplicates
		vector<string> vec_sam_lines_filtered;
		if (!filter_sam_pairs(vec_sam_lines, vec_sam_lines_filtered)) {
			cerr << "!error: filter_duplicates: filter_duplicates(): cannot filter sam pairs" << endl;
			return false;
		}
		// print list
		for (int i=0; i<vec_sam_lines_filtered.size(); i++)
			cout << vec_sam_lines_filtered[i] << endl;
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
