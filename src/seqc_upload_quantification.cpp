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
	cout << "seqc_upload_quantification" << endl;
	cout << "\t-dbhost [db hostname]" << endl;
	cout << "\t-dbuser [db user name]" << endl;
	cout << "\t-dbpass [db password]" << endl;
	cout << "\t-dbname [db name]" << endl;
	cout << endl;
	cout << "\t-id [unique id of record with data field]" << endl;
	cout << "\t-type [gene, isoform, exon, or ercc]" << endl;
	cout << "\t-f [input file]" << endl;
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
	string str_type;
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
	if (!cl.getArg(str_type, "type")) {
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
		cerr << "!error: seqc_upload_quantification: cannot init mysql" << endl;
		mysql_library_end();
		return 1;
	}
	if (!mysql_real_connect(mysql, str_db_host.c_str(), str_db_user.c_str(), str_db_pass.c_str(), str_db_name.c_str(), 0, NULL, 0)) {
		cerr << "!error: seqc_upload_quantification: cannot connect to mysql database" << endl;
		mysql_library_end();
		return 1;
	}

	//cout << "get ordered template sequence annotation" << endl;
	
	vector<string> vec_feature_list;
	if (!get_ordered_feature_list(
		mysql,
		str_type,
		vec_feature_list
	)) {
		cerr << "!error: seqc_upload_quantification: error getting ordered feature list" << endl;
		mysql_close(mysql);
		mysql_library_end();
		return 1;
	}

	int i_num_features = vec_feature_list.size();
	map<string,float> map_vals;
	for (int i=0; i<i_num_features; i++) {
		map_vals[vec_feature_list[i]] = 0;
	}

	// read cufflinks file
	FILE* file;
	file = fopen(str_file.c_str(), "r");
	if (file == NULL) {
		cerr << "!error: seqc_upload_quantification: cannot open file" << endl;
		mysql_close(mysql);
		mysql_library_end();
		return 1;
	}
	char buf[MAX_BUF_LEN];
	// read header line
	while (read_line(file, buf)) {
		string str_buf = string(buf);
		StringTokenizer st(str_buf,0,true);
		vector<string> vec_str_buf = st.split();
		string str_feature = vec_str_buf[0];

		map<string,float>:: iterator it;
		it = map_vals.find(str_feature);
		if (it == map_vals.end()) {
			cerr << "!warning: seqc_upload_quantification: extra feature: " << str_feature << endl;
		} else {
			float f_val = atof(vec_str_buf[1].c_str());
			map_vals[str_feature] = f_val;
		}
	}
	
	fclose(file);

	Matrix<float> mat_x(i_num_features, 1, 0);
	for (int i=0; i<i_num_features; i++) {
		mat_x[i][0] = map_vals[vec_feature_list[i]];
	}

	string str_table = string("quantification");
	if (!put_sample(
		mysql,
		str_id,
		str_table,
		mat_x,
		0
	)) {
		cerr << "!error: seqc_upload_quantification: cannot upload sample to db" << endl;
		mysql_close(mysql);
		mysql_library_end();
		return 1;
	}

	// close database connection
	mysql_close(mysql);
	mysql_library_end();

	return 0;
}
