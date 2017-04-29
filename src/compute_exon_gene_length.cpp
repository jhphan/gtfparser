#include <stdio.h>

#include <iostream>
#include <vector>
#include <map>

#include <miblab/commandline.h>
#include <miblab/stringtokenizer.h>

#include "../include/transcript.h"
#include "../include/sam.h"

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
	cout << "compute_exon_gene_length" << endl;
	cout << "\t-f [sorted gene to exon file]" << endl;
}

bool load_transcript_structure(map<string,Transcript>& map_gene_structure, string str_file) {

	// use the transcript structure to store all exons for a gene

	FILE* file;
	file = fopen(str_file.c_str(), "r");
	if (file == NULL) {
		cerr << "!error: compute_exon_gene_length: load_transcript_structure: cannot open file" << endl;
		return false;
	}

	char ch_buf[MAX_BUF_LEN];
	while (read_line(file, ch_buf)) {
		vector<string> vec_line = StringTokenizer(ch_buf, '\t').split();

		string str_gene_id = vec_line[0];
		vector<string> vec_chr = StringTokenizer(vec_line[1], ':').split();
		vector<string> vec_coords = StringTokenizer(vec_chr[1], '-').split();
		
		string str_chromosome = vec_chr[0];
		unsigned long ul_start_coord = atol(vec_coords[0].c_str());
		unsigned long ul_end_coord = atol(vec_coords[1].c_str());

		map<string,Transcript>::iterator it;
		it = map_gene_structure.find(str_gene_id);
		if (it == map_gene_structure.end()) {
			// new gene
			Transcript trans(str_gene_id, str_gene_id, str_chromosome, 1);
			trans.add_exon(ul_start_coord, ul_end_coord);
			map_gene_structure[str_gene_id] = trans;
		} else {
			// add exon to existing gene
			if (!it->second.add_exon(ul_start_coord, ul_end_coord)) {
				cout << str_gene_id << endl;
				return false;
			}
		}
	}

	fclose(file);
	return true;

}

int main(int argc, char* argv[]) {

	string str_file;

	CommandLine cl(argc, argv);

	if (!cl.getArg(str_file, "f")) {
		usage();
		return 1;
	}

	// load transcriptome structure
	map<string,Transcript> map_gene_structure;
	if (!load_transcript_structure(map_gene_structure, str_file)) {
		cerr << "!error: compute_exon_gene_length: cannot load transcript structure" << endl;
		return 1;
	}


	map<string,Transcript>::iterator it;
	for (it=map_gene_structure.begin(); it != map_gene_structure.end(); it++) {
		cout << it->first << "\t" << it->second.get_length() << endl;
	}

	return 0;
}
