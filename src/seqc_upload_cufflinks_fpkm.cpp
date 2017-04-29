#include <iostream>
#include <stdio.h>
#include <miblab/matrix.h>
#include <string>
#include <map>
#include <miblab/datafileheader.h>
#include <miblab/commandline.h>
#include <miblab/stringtokenizer.h>
#include <mysql/mysql.h>
#include "../include/seqc_db.h"

#define MAX_BUF_LEN 1048576

using namespace std;

void usage() {
	cout << "seqc_upload_cufflinks_fpkm" << endl;
	cout << "\t-dbhost [db hostname]" << endl;
	cout << "\t-dbuser [db user name]" << endl;
	cout << "\t-dbpass [db password]" << endl;
	cout << "\t-dbname [db name]" << endl;
	cout << endl;
	cout << "\t-id [isoform or other expression id]" << endl;
	cout << "\t-template [template name]" << endl;
	cout << "\t-f [input cufflinks file]" << endl;
}


int read_line(FILE* f, char* line) {
	int size;
	char ch;
	size = 0;

	ch = fgetc(f);
	while ((ch != '\n') && (ch != EOF) && (size<MAX_BUF_LEN-1)) {
		line[size] = ch;
		ch = fgetc(f);
		size++;
	}
	line[size] = '\0';
	if (ch == EOF) return 0;
	return size+1;
}


int main(int argc, char* argv[]) {

	// process command line parameters
	// database parameters
	string str_db_host;
	string str_db_user;
	string str_db_pass;
	string str_db_name;
	string str_id;
	string str_template_name;
	string str_file;
	
	CommandLine cl(argc, argv);

	if (!cl.getArg(str_db_host, "dbhost")) {
		usage();
		return 1;
	}
	if (!cl.getArg(str_db_user, "dbuser")) {
		usage();
		return 1;
	}
	if (!cl.getArg(str_db_pass, "dbpass")) {
		usage();
		return 1;
	}
	if (!cl.getArg(str_db_name, "dbname")) {
		usage();
		return 1;
	}
	if (!cl.getArg(str_id, "id")) {
		usage();
		return 1;
	}
	if (!cl.getArg(str_template_name, "template")) {
		usage();
		return 1;
	}
	if (!cl.getArg(str_file, "f")) {
		usage();
		return 1;
	}
	

	// open database connection
	mysql_library_init(0, NULL, NULL);
	MYSQL* mysql = mysql_init(NULL);
	if (mysql == NULL) {
		cerr << "!error: seqc_upload_cufflinks_fpkm: cannot init mysql" << endl;
		mysql_library_end();
		return 1;
	}
	if (!mysql_real_connect(mysql, str_db_host.c_str(), str_db_user.c_str(), str_db_pass.c_str(), str_db_name.c_str(), 0, NULL, 0)) {
		cerr << "!error: seqc_upload_cufflinks_fpkm: cannot connect to mysql database" << endl;
		mysql_library_end();
		return 1;
	}

	//cout << "get ordered template sequence annotation" << endl;
	
	vector<string> vec_template_annotations;
	if (!get_ordered_template_sequence_annotation1(
		mysql,
		str_template_name,
		vec_template_annotations
	)) {
		cerr << "!error: seqc_upload_cufflinks_fpkm: error getting template annotations" << endl;
		mysql_close(mysql);
		mysql_library_end();
		return 1;
	}

	int i_num_annotations = vec_template_annotations.size();
	//cout << "num annotations: " << i_num_annotations << endl;
	map<string,float> map_fpkm;
	for (int i=0; i<i_num_annotations; i++) {
	//	cout << vec_template_annotations[i] << endl;
		map_fpkm[vec_template_annotations[i]] = 0;
	}

	//cout << "reading file..." << endl;

	// read cufflinks file
	FILE* file;
	file = fopen(str_file.c_str(), "r");
	if (file == NULL) {
		cerr << "!error: seqc_upload_cufflinks_fpkm: cannot open file" << endl;
		mysql_close(mysql);
		mysql_library_end();
		return 1;
	}
	char buf[MAX_BUF_LEN];
	// read header line
	read_line(file, buf);
	//cout << "header: " << buf << endl;
	while (read_line(file, buf)) {
		string str_buf = string(buf);
		//cout << str_buf << endl;
		StringTokenizer st(str_buf,0,true);
		vector<string> vec_str_buf = st.split();
		//cout << "size: " << vec_str_buf.size() << endl;
		//cout << vec_str_buf[6] << "\t" << vec_str_buf[9] << endl;
		string str_map_loc = vec_str_buf[6];
		StringTokenizer st_map_loc(str_map_loc,':',true);
		string str_ercc;
		st_map_loc.getNextToken(str_ercc);
		map<string,float>:: iterator it;
		it = map_fpkm.find(str_ercc);
		if (it == map_fpkm.end()) {
			cerr << "!error: seqc_upload_cufflinks_fpkm: mapped location not found in db" << endl;
			fclose(file);
			mysql_close(mysql);
			mysql_library_end();
			return 1;
		}
		float f_val = atof(vec_str_buf[9].c_str());
		if (f_val > map_fpkm[str_ercc]) {
			map_fpkm[str_ercc] = f_val;
		}
	}
	
	fclose(file);

	Matrix<float> mat_x(i_num_annotations, 1, 0);
	for (int i=0; i<i_num_annotations; i++) {
		//cout << map_fpkm[vec_template_annotations[i]] << endl;
		mat_x[i][0] = map_fpkm[vec_template_annotations[i]];
	}

	string str_table = string("other_expression");
	if (!put_sample(
		mysql,
		str_id,
		str_table,
		mat_x,
		0
	)) {
		cerr << "!error: seqc_upload_cufflinks_fpkm: cannot upload sample to db" << endl;
		mysql_close(mysql);
		mysql_library_end();
		return 1;
	}

	// close database connection
	mysql_close(mysql);
	mysql_library_end();

	return 0;
}
