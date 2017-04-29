#include <stdio.h>

#include <iostream>
#include <vector>

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
	cout << "filter_sam" << endl;
	cout << "\t-sam [sam/bam file]" << endl;
	cout << "\t-filter [mapped/unmapped]" << endl;
}

bool filter_sam(string str_sam_file, bool b_mapped) {

	FILE* file_sam;
	string str_cmd = string("java -jar /data/ngs/bin/picard/ViewSam.jar I=")+str_sam_file+string(" VALIDATION_STRINGENCY=SILENT");
	file_sam = popen(str_cmd.c_str(), "r");
	if (file_sam == NULL) {
		cerr << "!error: filter_sam: filter_sam(): cannot open: " << str_cmd << endl;
		return false;
	}

	char ch_buf1[MAX_BUF_LEN];
	char ch_buf2[MAX_BUF_LEN];
	while (read_line(file_sam, ch_buf1)) {
		if (ch_buf1[0] == '@') {
			// header line
			cout << ch_buf1 << endl;
		} else {
			vector<string> vec_sam_1 = StringTokenizer(ch_buf1, '\t').split();
			// read next line
			if (!read_line(file_sam, ch_buf2)) {
				cerr << "!error: filter_sam: filter_sam(): odd number of reads?" << endl;
				return false;
			}
			vector<string> vec_sam_2 = StringTokenizer(ch_buf2, '\t').split();
			if (vec_sam_1[0] != vec_sam_2[0]) {
				cerr << "!error: filter_sam: filter_sam(): adjacent reads don't match, must be sorted by queryname" << endl;
				return false;
			}
			int i_flag_1 = atoi(vec_sam_1[1].c_str());
			int i_flag_2 = atoi(vec_sam_2[1].c_str());
			if (
				(i_flag_1 == 83 && i_flag_2 == 163) ||
				(i_flag_1 == 163 && i_flag_2 == 83) ||
				(i_flag_1 == 99 && i_flag_2 == 147) ||
				(i_flag_1 == 147 && i_flag_2 == 99)
			) {
				// correctly paired and mapped
				if (b_mapped) {
					cout << ch_buf1 << endl;
					cout << ch_buf2 << endl;
				}
			} else {
				// unmapped
				if (!b_mapped) {
					cout << ch_buf1 << endl;
					cout << ch_buf2 << endl;
				}
			}
		}
	}

	fclose(file_sam);

	return true;
}


int main(int argc, char* argv[]) {

	string str_sam_file;
	string str_filter;
	bool b_mapped;

	CommandLine cl(argc, argv);

	if (!cl.getArg(str_sam_file, "sam")) {
		usage();
		return 1;
	}
	if (!cl.getArg(str_filter, "filter")) {
		usage();
		return 1;
	}

	if (str_filter == "mapped") {
		b_mapped = true;
	} else if (str_filter == "unmapped") {
		b_mapped = false;
	} else {
		usage();
		return 1;
	}

	// read sam file
	if (!filter_sam(str_sam_file, b_mapped)) {
		cerr << "!error: filter_sam: cannot read sam file" << endl;
		return 1;
	}

	return 0;
}
