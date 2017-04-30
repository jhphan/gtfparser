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

bool pair_consistent(vector<string>& vec_sam_1, vector<string>& vec_sam_2) {
	if (vec_sam_1[2] != vec_sam_2[2]) {
		// pair doesn't map to the same reference sequence
		return false;
	}
	if (vec_sam_1[3] != vec_sam_2[7]) {
		// first end and second end mate pos do not match
		return false;
	}
	if (vec_sam_2[3] != vec_sam_1[7]) {
		// second end and first end mate pos do not match
		return false;
	}
	return true;
}

bool filter_sam(string str_sam_file, bool b_mapped) {

	FILE* file_sam;
	string str_cmd = string("java -jar /data/ngs/bin/picard/ViewSam.jar I=")+str_sam_file+string(" VALIDATION_STRINGENCY=SILENT");
	file_sam = popen(str_cmd.c_str(), "r");
	if (file_sam == NULL) {
		cerr << "!error: filter_sam: filter_sam(): cannot open: " << str_cmd << endl;
		return false;
	}

	string str_sam_1;
	string str_sam_2;
	vector<string> vec_sam_1 = vector<string>();
	vector<string> vec_sam_2 = vector<string>();
	int i_flag_1 = 0;
	int i_flag_2 = 0;

	char ch_buf[MAX_BUF_LEN];
	while (read_line(file_sam, ch_buf)) {
		if (ch_buf[0] == '@') {
			// header line
			cout << ch_buf << endl;
		} else {
			if (i_flag_1 == 0) {
				// first of pair has not been found
				str_sam_1 = string(ch_buf);
				vec_sam_1 = StringTokenizer(str_sam_1, '\t').split();
				i_flag_1 = atoi(vec_sam_1[1].c_str());
				
				bool b_valid = false;
				if (
					i_flag_1 == 83 || i_flag_1 == 163
					|| i_flag_1 == 99 || i_flag_1 == 147
					|| i_flag_1 == 355 || i_flag_1 == 403
					|| i_flag_1 == 339 || i_flag_1 == 419
				) b_valid = true;


				if (!b_valid && i_flag_1 > 255) {
					i_flag_1 = 0;
				}
//				if ( 
//					(!b_valid && b_mapped)
//					|| (!b_mapped && i_flag_1 > 255)
//				) {
//					i_flag_1 = 0;
//				}
			} else {
				if (i_flag_2 == 0) {
					// first of pair has been found, but not second
					str_sam_2 = string(ch_buf);
					vec_sam_2 = StringTokenizer(str_sam_2, '\t').split();
					i_flag_2 = atoi(vec_sam_2[1].c_str());

					bool b_valid = false;
					if (
						i_flag_2 == 83 || i_flag_2 == 163
						|| i_flag_2 == 99 || i_flag_2 == 147
						|| i_flag_2 == 355 || i_flag_2 == 403
						|| i_flag_2 == 339 || i_flag_2 == 419
					) b_valid = true;

					if (!b_valid && i_flag_2 > 255) {
						i_flag_2 = 0;
					}
//					if (
//						(!b_valid && b_mapped)
//						|| (!b_mapped && i_flag_2 > 255)
//					) {
//						i_flag_2 = 0;
//					}
				}
			}

			if (i_flag_1 > 0 && i_flag_2 > 0) {
				// found a pair
				if (vec_sam_1[0] != vec_sam_2[0]) {
					cerr << "!error: filter_sam: filter_sam(): adjacent reads don't match" << endl;
					cerr << str_sam_1 << endl;
					cerr << str_sam_2 << endl;
					return false;
				}

				// check if it's a correct pair
				if (
					(i_flag_1 == 83 && i_flag_2 == 163) ||
					(i_flag_1 == 163 && i_flag_2 == 83) ||
					(i_flag_1 == 99 && i_flag_2 == 147) ||
					(i_flag_1 == 147 && i_flag_2 == 99) ||
					(i_flag_1 == 355 && i_flag_2 == 403) ||
					(i_flag_1 == 403 && i_flag_2 == 355) ||
					(i_flag_1 == 339 && i_flag_2 == 419) ||
					(i_flag_1 == 419 && i_flag_2 == 339)
				) {
					if (pair_consistent(vec_sam_1, vec_sam_2)) {
						if (b_mapped) {
							cout << str_sam_1 << endl;
							cout << str_sam_2 << endl;
						}
					} else {
						if (!b_mapped) {
							if (i_flag_1 < 256) { // ignore secondary hits
								cout << str_sam_1 << endl;
								cout << str_sam_2 << endl;
							}
						}
					}
				} else {
					if (!b_mapped) {
						cout << str_sam_1 << endl;
						cout << str_sam_2 << endl;
					}
				}
				i_flag_1 = 0;
				i_flag_2 = 0;
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
