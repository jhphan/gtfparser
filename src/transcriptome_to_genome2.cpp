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
	cout << "transcriptome_to_genome" << endl;
	cout << "\t-sam [sam/bam file]" << endl;
	cout << "\t-gtf [gtf file]" << endl;
}
/*
bool load_reference_coords(map<string,reference_coords>& map_reference_coords, string str_ref_file) {
	FILE* file_ref;
	file_ref = fopen(str_ref_file.c_str(), "r");
	if (file_ref == NULL) {
		cerr << "!error: transcriptome_to_genome: load_reference_coords(): cannot open reference coord file" << endl;
		return false;
	}

	char ch_buf[MAX_BUF_LEN];
	while (read_line(file_ref, ch_buf)) {
		vector<string> vec_header_pieces = StringTokenizer(ch_buf, ' ').split();
		if (vec_header_pieces.size() != 3) {
			cerr << "!error: transcriptome_to_genome: load_reference_coords(): malformed header: " << ch_buf << endl;
			return false;
		}
		
		struct reference_coords rc;
		string str_ref_id = vec_header_pieces[0].substr(1,vec_header_pieces[0].size()-1);
		StringTokenizer st(vec_header_pieces[1], ':');
		if (!st.getNextToken(rc.str_chromosome)) {
			cerr << "!error: transcriptome_to_genome: load_reference_coords(): cannot get chromosome name: " << ch_buf << endl;
			return false;
		}
		string str_chromosome_coords;
		if (!st.getNextToken(str_chromosome_coords)) {
			cerr << "!error: transcriptome_to_genome: load_reference_coords(): cannot get chromosome coords: " << ch_buf << endl;
			return false;
		}
		st = StringTokenizer(str_chromosome_coords, '-');
		int i_coord;
		if (!st.getNextToken(i_coord)) {
			cerr << "!error: transcriptome_to_genome: load_reference_coords(): cannot get chromosome start: " << ch_buf << endl;
			return false;
		}
		rc.ul_start_coord = (unsigned long)i_coord;
		if (!st.getNextToken(i_coord)) {
			cerr << "!error: transcriptome_to_genome: load_reference_coords(): cannot get chromosome end: " << ch_buf << endl;
			return false;
		}
		rc.ul_end_coord = (unsigned long)i_coord;
		rc.str_transcript_id = vec_header_pieces[2];

		map<string,reference_coords>::iterator it;
		it = map_reference_coords.find(str_ref_id);
		if (it == map_reference_coords.end()) {
			// add reference element
			map_reference_coords[str_ref_id] = rc;
		} else {
			cerr << "!error: transcriptome_to_genome: load_reference_coords(): duplicate reference id: " << str_ref_id << endl;
			return false;
		}
	}

	fclose(file_ref);

	return true;
}
*/
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

bool translate_sam_file(map<string,Transcript>& map_transcript_structure, string str_sam_file) {

	FILE* file_sam;
	string str_cmd = string("java -jar /data/ngs/bin/picard/ViewSam.jar I=")+str_sam_file+string(" VALIDATION_STRINGENCY=SILENT");
	file_sam = popen(str_cmd.c_str(), "r");
	if (file_sam == NULL) {
		cerr << "!error: transcriptome_to_genome: translate_sam_file(): cannot open: " << str_cmd << endl;
		return false;
	}

	char ch_buf[MAX_BUF_LEN];
	struct sam_record sr1;
	struct sam_record sr2;
	bool done = false;
	bool read_first_line = true;
	bool b_error_first = false;
	bool b_error_second = false;
	while (!done) {
		if (read_first_line) {
			// find first of pair
			if (!read_line(file_sam, ch_buf)) break;
			while (ch_buf[0] == '@') {
				// header
				cout << ch_buf << endl;
				if (!read_line(file_sam, ch_buf)) {
					done = true;
					break;
				}
			}
			if (done) break;
			parse_sam_record(ch_buf, sr1);
			if (!translate_sam_record(sr1, map_transcript_structure)) {
				cerr << "!error: transcriptome_to_genome: translate_sam_file(): cannot translate sam record" << endl;
				cerr << ch_buf << endl;
				b_error_first = true;
				//return false;
			}
		}

		// find second of pair (if paired)
		if (!read_line(file_sam, ch_buf)) {
			if (!b_error_first) print_sam_record(sr1);	// not paired
			break;
		}
		while (ch_buf[0] == '@') {
			cout << ch_buf << endl;
			if (!read_line(file_sam, ch_buf)) {
				if (!b_error_first) print_sam_record(sr1);	// not paired
				done = true;
				break;
			}
		}
		if (done) break;
		parse_sam_record(ch_buf, sr2);
		if (!translate_sam_record(sr2, map_transcript_structure)) {
			cerr << "!error: transcriptome_to_genome: translate_sam_file(): cannot translate sam record" << endl;
			cerr << ch_buf << endl;
			b_error_second = true;
			//return false;
		}

		if (sr1.str_read_name == sr2.str_read_name) {
			if (!b_error_first && !b_error_second) {
				// if read names match
				// this is a paired set of reads
				fix_template_length(sr1, sr2);
				print_sam_record(sr1);
				print_sam_record(sr2);
			}
			b_error_first = false;
			b_error_second = false;
			read_first_line = true;
		} else {
			// not paired
			if (!b_error_first) print_sam_record(sr1);
			sr1 = sr2;
			b_error_first = b_error_second;
			b_error_second = false;
			read_first_line = false;
		}
	}

	fclose(file_sam);

	return true;
}


int main(int argc, char* argv[]) {

	string str_sam_file;
	string str_gtf_file;

	CommandLine cl(argc, argv);

	if (!cl.getArg(str_sam_file, "sam")) {
		usage();
		return 1;
	}
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

	// read sam file
	if (!translate_sam_file(map_transcript_structure, str_sam_file)) {
		cerr << "!error: transcriptome_to_genome: cannot translate sam file" << endl;
		return 1;
	}
	
	return 0;
}
