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
	cout << "merge_exons" << endl;
	cout << "\t-gtf [gtf file]" << endl;
}

bool load_transcript_structure(map<string,Transcript>& map_transcript_structure, string str_gtf_file) {
	FILE* file_gtf;
	file_gtf = fopen(str_gtf_file.c_str(), "r");
	if (file_gtf == NULL) {
		cerr << "!error: transcriptome_to_genome: load_transcript_structure: cannot open gtf file" << endl;
		return false;
	}

	char ch_buf[MAX_BUF_LEN];
	while (read_line(file_gtf, ch_buf)) {
		vector<string> vec_gtf = StringTokenizer(ch_buf, '\t').split();
		if (vec_gtf.size() != 9) {
			cerr << "!error: transcriptome_to_genome: load_transcriptome_structure: malformed gtf: " << ch_buf << endl;
			return false;
		}

		if (vec_gtf[2] == "exon") {
			string str_chromosome = vec_gtf[0];
			unsigned long ul_start_coord = atol(vec_gtf[3].c_str());
			unsigned long ul_end_coord = atol(vec_gtf[4].c_str());
			int i_strand = 0;
			if (vec_gtf[6] == "+") {
				i_strand = 1;
			} else if (vec_gtf[6] == "-") {
				i_strand = -1;
			} else {
				cerr << "!error: transcriptome_to_genome: load_transcriptome_structure: invalid strand: " << ch_buf << endl;
				return false;
			}
			vector<string> vec_attrib = StringTokenizer(vec_gtf[8], ' ').split();
			string str_transcript_id = vec_attrib[1].substr(1, vec_attrib[1].size()-3);
			string str_gene_id = vec_attrib[3].substr(1, vec_attrib[3].size()-3);
			
			map<string,Transcript>::iterator it;
			it = map_transcript_structure.find(str_transcript_id);
			if (it == map_transcript_structure.end()) {
				// new transcript
				Transcript trans(str_transcript_id, str_gene_id, str_chromosome, i_strand);
				trans.add_exon(ul_start_coord, ul_end_coord);
				map_transcript_structure[str_transcript_id] = trans;
			} else {
				// add exon to existing transcript
				it->second.add_exon(ul_start_coord, ul_end_coord);
			}
		}
	}

	fclose(file_gtf);
	return true;

}

int main(int argc, char* argv[]) {

	string str_gtf_file;

	CommandLine cl(argc, argv);

	if (!cl.getArg(str_gtf_file, "gtf")) {
		usage();
		return 1;
	}

	// load transcriptome structure
	map<string,Transcript> map_transcript_structure;
	if (!load_transcript_structure(map_transcript_structure, str_gtf_file)) {
		cerr << "!error: transcriptome_to_genome: cannot load transcript structure" << endl;
		return 1;
	}


	map<string,Transcript>::iterator it;
	for (it=map_transcript_structure.begin(); it != map_transcript_structure.end(); it++) {
		int i_num_exons = it->second.get_num_exons();
		for (int i=0; i<i_num_exons; i++) {
			cout << it->second.get_chromosome();
			cout << "\tMerged";
			cout << "\texon";
			cout << "\t" << it->second.get_exon_start(i);
			cout << "\t" << it->second.get_exon_end(i);
			cout << "\t.";
			if (it->second.get_strand() > 0) {
				cout << "\t+";
			} else {
				cout << "\t-";
			}
			cout << "\t.";
			cout << "\ttranscript_id \"" << it->second.get_name() << "\"; ";
			cout << "gene_id \"" << it->second.get_gene_name() << "\";";
			cout << endl;
		}
	}

	return 0;
}
