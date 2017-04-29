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
	cout << "read_sam_flags" << endl;
	cout << "\t[ -sam <sam file> || -bam <bam file> ]" << endl;
}

bool read_file(map<vector<int>,unsigned int>& map_sam_flags, string str_sam_file, string file_type) {

	FILE* file_sam;
	string str_cmd;
	if (!file_type.compare("sam")) {
		str_cmd = string("/data/ngs/bin/samtools/samtools view -S -h ")+str_sam_file;
	}
	else {
		str_cmd = string("/data/ngs/bin/samtools/samtools view -h ")+str_sam_file;
	}

	file_sam = popen(str_cmd.c_str(), "r");
	if (file_sam == NULL) {
		cerr << "!error: read_sam_flags: read_file: cannot open: " << str_cmd << endl;
		return false;
	}

	char ch_buf[MAX_BUF_LEN];
	while (read_line(file_sam, ch_buf)) {

		if (ch_buf[0] == '@') {
			// header line
			string sch_buf = ch_buf;
			cout << sch_buf.erase(sch_buf.find("query")) + string("unsorted") << endl;
			cout << ch_buf[0-2] << endl;
			//cout << ch_buf << endl;
getchar();
		} else {
			vector<string> vec_sam_1 = StringTokenizer(ch_buf, '\t').split();
			// read next line
			if (!read_line(file_sam, ch_buf)) {
				cerr << "!error: read_sam_flags: read_file(): odd number of reads?" << endl;
				return false;
			}
			vector<string> vec_sam_2 = StringTokenizer(ch_buf, '\t').split();
			if (vec_sam_1[0] != vec_sam_2[0]) {
				cerr << "!error: read_sam_flags: read_file(): adjacent reads don't match" << endl;
				return false;
			}
			int i_flag_1 = atoi(vec_sam_1[1].c_str());
			int i_flag_2 = atoi(vec_sam_2[1].c_str());

			vector<int> vec_flag_pair = vector<int>(2,0);
			vec_flag_pair[0] = i_flag_1;
			vec_flag_pair[1] = i_flag_2;

			map<vector<int>,unsigned int>::iterator it;
			it = map_sam_flags.find(vec_flag_pair);
			if (it == map_sam_flags.end()) {
				// new pair
				map_sam_flags[vec_flag_pair] = 1;
			} else {
				it->second++;
			}			
		}
	}

	fclose(file_sam);

	return true;
}


int main(int argc, char* argv[]) {

	string str_sam_file;
	string str_bam_file;
	string file_type = "sam";

	CommandLine cl(argc, argv);

	if (!cl.getArg(str_sam_file, "sam")) {
		file_type = "bam";
		if (!cl.getArg(str_bam_file, "bam")) {
			
			usage();
			return 1;
		}
	}

	// read sam file
	map<vector<int>,unsigned int> map_sam_flags;

	if (!str_sam_file.empty()){
		if (!read_file(map_sam_flags, str_sam_file, file_type)) {
			cerr << "!error: wham_fixing: cannot read sam file" << endl;
			return 1;
		}
	}
	else
	{
		if (!read_file(map_sam_flags, str_bam_file, file_type)) {
			cerr << "!error: wham_fixing: cannot read bam file" << endl;
			return 1;
		}
	}

	unsigned int ui_total_pairs = 0;

	map<vector<int>,unsigned int>::iterator it;
	for (it = map_sam_flags.begin(); it != map_sam_flags.end(); it++) {
		cout << "(" << it->first[0] << "," << it->first[1] << ") => " << it->second << endl;
		ui_total_pairs += it->second;
	}

	cout << "Total Pairs: " << ui_total_pairs << endl;
	
	return 0;
}
