#include "../include/seqc_db.h"

bool get_ordered_feature_list(
	MYSQL* mysql,
	string& str_type,
	vector<string>& vec_str_list
) {

	// get annotation names
	string str_query
		= string("select id from ")+str_type+string(" order by order_index asc");

	if (mysql_real_query(
		mysql,
		str_query.c_str(),
		str_query.size()
	)) {
		cerr << "!error: get_ordered_feature_list(): query failed" << endl;
		return false;
	}

	MYSQL_RES* mysql_result = mysql_store_result(mysql);

	// make sure database is consistent
	int i_num_rows = mysql_num_rows(mysql_result);
	vec_str_list = vector<string>(i_num_rows);

	MYSQL_ROW mysql_row;

	for (int i=0; i<i_num_rows; i++) {
		mysql_row = mysql_fetch_row(mysql_result);
		vec_str_list[i] = string(mysql_row[0]);
	}

	mysql_free_result(mysql_result);

	return true;
}

bool put_sample(
	MYSQL* mysql,
	string& str_id,
	string& str_table,
	Matrix<float>& mat_x,
	int i_col=0
) {

	int i_dimensions = mat_x.getHeight();

	int i_datasize = sizeof(float)*i_dimensions;

	string str_query_pre;
	str_query_pre = string("update ")+str_table+string(" set data = '");
	int i_presize = str_query_pre.size();

	string str_query_post;
	str_query_post = string("' where id = '")+str_id+string("'");
	int i_postsize = str_query_post.size();

	char* ch_all_data = new char[i_datasize];
	char* ch_all_data_escape = new char[i_datasize*2];

	for (int i=0; i<i_dimensions; i++)
		memcpy(
			ch_all_data+i*sizeof(float),
			&(mat_x[i][i_col]),
			sizeof(float)
		);

	int i_datasize_escape = mysql_real_escape_string(
		mysql,
		ch_all_data_escape,
		ch_all_data,
		i_datasize
	);
	delete[] ch_all_data;

	char* ch_query
		= new char[i_presize+i_postsize+(i_datasize_escape)+1];
	memcpy(
		ch_query,
		str_query_pre.c_str(),
		i_presize
	);
	
	memcpy(
		ch_query+i_presize,
		ch_all_data_escape,
		i_datasize_escape
	);
	delete[] ch_all_data_escape;

	memcpy(
		ch_query+i_presize+i_datasize_escape,
		str_query_post.c_str(),
		i_postsize
	);
	ch_query[i_presize+i_datasize_escape+i_postsize] = 0;

	int i_result = mysql_real_query(
		mysql,
		ch_query,
		i_presize+i_datasize_escape+i_postsize
	);
	
	delete[] ch_query;

	if (i_result != 0) {
		cerr << mysql_error(mysql) << endl;
		cerr << "!error: put_sample(): mysql_real_query() failed" << endl;
		return false;
	}

	return true;
}


bool get_sample(
	MYSQL* mysql,
	string& str_id,
	Matrix<float>& mat_x,
	int i_col
) {

	string str_query 
		= string("select data from quantification where id = '")
			+str_id
			+string("'");

	if (mysql_real_query(
		mysql,
		str_query.c_str(),
		str_query.size()
	)) {
		cerr << "!error: get_sample(): query failed" << endl;
		return false;
	}

	MYSQL_RES* mysql_result = mysql_store_result(mysql);
	MYSQL_ROW mysql_row = mysql_fetch_row(mysql_result);
	if (!mysql_row) {
		cerr << "!error: get_sample(): sample "
			<< str_id << endl;
		mysql_free_result(mysql_result);
		return false;
	}

	if (!mysql_row[0]) {
		cerr << "!error: get_sample(): data null, "
			<< str_id << endl;
		mysql_free_result(mysql_result);
		return false;
	}

	unsigned long* l_lengths = mysql_fetch_lengths(mysql_result);

	if (l_lengths[0]%4 > 0) {
		cerr << "!error: get_sample(): corrupt data" << endl;
		mysql_free_result(mysql_result);
		return false;
	}
	if (l_lengths[0]/4 != mat_x.getHeight()) {
		cerr << "!error: get_sample(): matrix and data sizes do not agree" << endl;
		mysql_free_result(mysql_result);
		return false;
	}

	// copy data to matrix
	for (int i=0; i<l_lengths[0]/4; i++) {
		memcpy(&(mat_x[i][i_col]), mysql_row[0]+i*4, 4);
	}
	
	mysql_free_result(mysql_result);
	return true;
}
