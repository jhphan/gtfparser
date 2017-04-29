#include <stdlib.h>
#include <iostream>
#include <string>
#include <mysql/mysql.h>
#include <miblab/matrix.h>
#include <miblab/strconv.h>
#include <miblab/stringtokenizer.h>

struct sequence {
	string str_name;
	int i_total_reads;
	string str_platform_id;
	string str_site;
	string str_sample;
	string str_replicate;
	string str_lane;
	string str_barcode;
	string str_flowcell_id;
	string str_path1;
	string str_path2;
	string str_path3;
	string str_path4;
};

bool get_ordered_feature_list(
	MYSQL* mysql,
	string& str_type,
	vector<string>& vec_str_list
);

bool put_sample(
	MYSQL* mysql,
	string& str_id,
	string& str_table,
	Matrix<float>& mat_x,
	int i_col
);

bool get_sample(
	MYSQL* mysql,
	string& str_id,
	Matrix<float>& mat_x,
	int i_col
);

