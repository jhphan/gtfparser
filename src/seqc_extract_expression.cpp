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
	cout << "seqc_extract_expression" << endl;
	cout << "\t-dbhost [db hostname]" << endl;
	cout << "\t-dbuser [db user name]" << endl;
	cout << "\t-dbpass [db password]" << endl;
	cout << "\t-dbname [db name]" << endl;
	cout << endl;
	cout << "\t-f [output file]" << endl;
	cout << "\t-mapping_pipeline [name]" << endl;
	cout << "\t-mapping_pipeline_stage [stage]" << endl;
	cout << "\t-ref_type [gene,isoform,exon,ercc]" << endl;
	cout << "\t-algorithm [ht-seq,cufflinks,rsem]" << endl;
	cout << "\t-normalization [count,fpkm,fpm]" << endl;
	cout << "\t-platforms [platform ids: ILM,LIF,ROC]" << endl;
	cout << "\t-sites [comma separated list of sites]" << endl;
	cout << "\t-samples [comma separated list of samples]" << endl;
	cout << "\t-replicates [comma separated list of replicates]" << endl;
	cout << "\t-flowlanes [comma separated list of flowcell/lane combinations, e.g.: A:L01-L02_B:L01]" << endl;
}

bool get_sequence_data(
	MYSQL* mysql,
	vector<string>& vec_platforms,
	vector<string>& vec_sites,
	vector<string>& vec_samples,
	vector<string>& vec_replicates,
	vector<string>& vec_flowlanes,
	map<string,map<string,vector<struct sequence> > >& map_sequence_data,
	map<string,int>& map_flowlane_count
) {

	// construct sql query
	string str_query = "select * from sequence where";

	string str_tmp;
	// platform
	for (int i=0; i<vec_platforms.size(); i++) {
		if (i > 0) str_tmp += string(" or ");
		str_tmp += string("platform_id = '")+vec_platforms[i]+string("'");
	}
	str_query += string(" (")+str_tmp+string(")");

	str_tmp = "";
	// site
	for (int i=0; i<vec_sites.size(); i++) {
		if (i > 0) str_tmp += string(" or ");
		str_tmp += string("site = '")+vec_sites[i]+string("'");
	}
	str_query += string(" and (")+str_tmp+string(")");

	str_tmp = "";
	// samples
	for (int i=0; i<vec_samples.size(); i++) {
		if (i > 0) str_tmp += string(" or ");
		str_tmp += string("sample = '")+vec_samples[i]+string("'");
	}
	str_query += string(" and (")+str_tmp+string(")");

	str_tmp = "";
	// replicates
	for (int i=0; i<vec_replicates.size(); i++) {
		if (i > 0) str_tmp += string(" or ");
		str_tmp += string("replicate = '")+vec_replicates[i]+string("'");
	}
	str_query += string(" and (")+str_tmp+string(")");

	for (int i=0; i<vec_flowlanes.size(); i++) {
		int i_flowlane_count = 0;

		// split at _
		vector<string> vec_flowcell_split = StringTokenizer(vec_flowlanes[i],'_').split();
		
		string str_query_post = "";
		
		for (int j=0; j<vec_flowcell_split.size(); j++) {
			// split at :
			vector<string> vec_flowlane_split = StringTokenizer(vec_flowcell_split[j],':').split();
			if ( (vec_flowlane_split[0] != "A") & (vec_flowlane_split[0] != "B") ) {
				cerr << "!error: seqc_extract_expression: get_sequence_data(): invalid flowcell" << endl;
				return false;
			}
			if (vec_flowlane_split.size() != 2) {
				cerr << "!error: seqc_extract_expression: get_sequence_data(): invalid flowlane format" << endl;
				return false;
			}
			// split at '-'
			vector<string> vec_lane_split = StringTokenizer(vec_flowlane_split[1],'-').split();
			for (int k=0; k<vec_lane_split.size(); k++) {
				if (str_query_post.length() > 0) str_query_post += string(" or ");
				str_query_post 
					+= string("(flowcell_id like '")
					+vec_flowlane_split[0]
					+string("%' and lane = '")
					+vec_lane_split[k]
					+string("')");
				i_flowlane_count++;
			}
		}
		map_flowlane_count[vec_flowlanes[i]] = i_flowlane_count;

		string str_final = str_query + string(" and (")+str_query_post+string(")");
		//cout << str_final << endl;
	
		if (mysql_real_query(
			mysql,
			str_final.c_str(),
			str_final.size()
		)) {
			cerr << "!error: seqc_extract_expression: get_sequence_data(): cannot query mysql" << endl;
			return false;
		}
		MYSQL_RES* mysql_result = mysql_store_result(mysql);
		int i_num_rows = mysql_num_rows(mysql_result);
	
		MYSQL_ROW mysql_row;
		for (int j=0; j<i_num_rows; j++) {
			mysql_row = mysql_fetch_row(mysql_result);
			struct sequence seq;
			seq.str_name = string(mysql_row[0]);
			seq.i_total_reads = atoi(mysql_row[1]);
			seq.str_platform_id = string(mysql_row[2]);
			seq.str_site = string(mysql_row[3]);
			seq.str_sample = string(mysql_row[4]);
			seq.str_replicate = string(mysql_row[5]);
			seq.str_lane = string(mysql_row[6]);
			seq.str_barcode = string(mysql_row[7]);
			seq.str_flowcell_id = string(mysql_row[8]);
			seq.str_path1 = string(mysql_row[9]);
			seq.str_path2 = string(mysql_row[10]);
			seq.str_path3 = string(mysql_row[11]);
			seq.str_path4 = string(mysql_row[12]);
			
			string str_key = seq.str_platform_id
				+string("_")+seq.str_site
				+string("_")+seq.str_sample
				+string("_")+seq.str_replicate;
			//cout << str_key << endl;

			map_sequence_data[str_key][vec_flowlanes[i]].push_back(seq);
		}
		mysql_free_result(mysql_result);
	}

	return true;
}

bool get_mapping_pipeline_stage_id(
	MYSQL* mysql,
	string str_mapping_pipeline,
	string str_mapping_pipeline_stage,
	string& str_mapping_pipeline_stage_id
) {
	string str_query = 
		string("select mapping_pipeline_stage.id ")
		+string("from mapping_pipeline_stage,mapping_pipeline ")
		+string("where mapping_pipeline_stage.mapping_pipeline_id = mapping_pipeline.id ")
			+string("and mapping_pipeline.name = '")+str_mapping_pipeline+string("' ")
			+string("and mapping_pipeline_stage.number = ")+str_mapping_pipeline_stage;

	
	if (mysql_real_query(
		mysql,
		str_query.c_str(),
		str_query.size()
	)) {
		cerr << "!error: seqc_extract_expression: get_mapping_pipeline_stage_id(): cannot query mysql" << endl;
		return false;
	}
	MYSQL_RES* mysql_result = mysql_store_result(mysql);
	int i_num_rows = mysql_num_rows(mysql_result);
	if (i_num_rows == 0) {
		cerr << "!error: seqc_extract_expression: get_mapping_pipeline_stage_id(): mapping pipeline not found" << endl;
		return false;
	}
	if (i_num_rows != 1) {
		cerr << "!error: seqc_extract_expression: get_mapping_pipeline_stage_id(): multiple mapping pipelines found" << endl;
		return false;
	}

	MYSQL_ROW mysql_row = mysql_fetch_row(mysql_result);
	str_mapping_pipeline_stage_id = string(mysql_row[0]);
	mysql_free_result(mysql_result);

	return true;
}

bool get_mapping_id(
	MYSQL* mysql,
	string str_sequence_name,
	string str_mapping_pipeline_stage_id,
	string& str_mapping_id
) {
	string str_query = 
		string("select id from mapping ")
		+string("where sequence_name = '")+str_sequence_name+string("' ")
		+string("and mapping_pipeline_stage_id = '")+str_mapping_pipeline_stage_id+string("'");
	//cout << str_query << endl;
	if (mysql_real_query(
		mysql,
		str_query.c_str(),
		str_query.size()
	)) {
		cerr << "!error: seqc_extract_expression: get_mapping_id(): cannot query mysql" << endl;
		return false;
	}
	MYSQL_RES* mysql_result = mysql_store_result(mysql);
	int i_num_rows = mysql_num_rows(mysql_result);
	if (i_num_rows == 0) {
		cerr << "!error: seqc_extract_expression: get_mapping_id(): mapping not found" << endl;
		return false;
	}
	if (i_num_rows != 1) {
		cerr << "!error: seqc_extract_expression: get_mapping_id(): multiple mappings found" << endl;
		return false;
	}

	MYSQL_ROW mysql_row = mysql_fetch_row(mysql_result);
	str_mapping_id = string(mysql_row[0]);
	mysql_free_result(mysql_result);
	return true;
}

bool get_quantification_id(
	MYSQL* mysql,
	string str_ref_type,
	string str_algorithm,
	string str_normalization,
	vector<string>& vec_mapping_ids,
	string& str_quantification_id
) {
	string str_query =
		string("select distinct quantification.id from quantification,quantification_mapping_junction ")
		+string("where quantification.type = '")+str_ref_type+string("' ")
		+string("and quantification.algorithm = '")+str_algorithm+string("' ")
		+string("and quantification.normalization = '")+str_normalization+string("' ")
		+string("and quantification.id = quantification_mapping_junction.quantification_id");

	string str_query_post = "";
	for (int i=0; i<vec_mapping_ids.size(); i++) {
		if (i > 0) str_query_post += string(" or ");
		str_query_post += string("quantification_mapping_junction.mapping_id = '")+vec_mapping_ids[i]+string("'");
	}
	string str_final = str_query+string(" and (")+str_query_post+string(")");
	//cout << str_final << endl;
	if (mysql_real_query(
		mysql,
		str_final.c_str(),
		str_final.size()
	)) {
		cerr << "!error: seqc_extract_expression: get_quantification_id(): cannot query mysql" << endl;
		return false;
	}
	MYSQL_RES* mysql_result = mysql_store_result(mysql);
	int i_num_rows = mysql_num_rows(mysql_result);
	if (i_num_rows == 0) {
		cerr << "!error: seqc_extract_expression: get_quantification_id(): quantification id not found" << endl;
		return false;
	}
	vector<string> vec_quantification_ids(i_num_rows);
	for (int i=0; i<i_num_rows; i++) {
		MYSQL_ROW mysql_row = mysql_fetch_row(mysql_result);
		vec_quantification_ids[i] = string(mysql_row[0]);
	}
	mysql_free_result(mysql_result);

	// check each quantification id to make sure it's linked to the right number of mappings
	for (int i=0; i<i_num_rows; i++) {
		str_query = string("select count(*) from quantification_mapping_junction ")
			+string("where quantification_id = '")+vec_quantification_ids[i]+string("'");
		if (mysql_real_query(
			mysql,
			str_query.c_str(),
			str_query.size()
		)) {
			cerr << "!error: seqc_extract_expression: get_quantification_id(): cannot query for count" << endl;
			return false;
		}
		mysql_result = mysql_store_result(mysql);
		MYSQL_ROW mysql_row = mysql_fetch_row(mysql_result);
		if (atoi(mysql_row[0]) == vec_mapping_ids.size()) {
			mysql_free_result(mysql_result);
			str_quantification_id = vec_quantification_ids[i];
			return true;
		}
		mysql_free_result(mysql_result);
	}

	str_quantification_id = "";
	return false;
}

int main(int argc, char* argv[]) {

	// process command line parameters
	// database parameters
	string str_db_host;
	string str_db_user;
	string str_db_pass;
	string str_db_name;

	string str_file;

	string str_mapping_pipeline;
	string str_mapping_pipeline_stage;
	string str_ref_type;
	string str_algorithm;
	string str_normalization;

	string str_platforms;
	vector<string> vec_platforms;
	string str_samples;
	vector<string> vec_samples;
	string str_sites;
	vector<string> vec_sites;
	string str_replicates;
	vector<string> vec_replicates;
	string str_flowlanes;
	vector<string> vec_flowlanes;
	
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
	if (!cl.getArg(str_file, "f")) {
		usage();
		return 1;
	}
	if (!cl.getArg(str_mapping_pipeline, "mapping_pipeline")) {
		usage();
		return 1;
	}
	if (!cl.getArg(str_mapping_pipeline_stage, "mapping_pipeline_stage")) {
		usage();
		return 1;
	}
	if (!cl.getArg(str_ref_type, "ref_type")) {
		usage();
		return 1;
	}
	if (!cl.getArg(str_algorithm, "algorithm")) {
		usage();
		return 1;
	}
	if (!cl.getArg(str_normalization, "normalization")) {
		usage();
		return 1;
	}
	if (!cl.getArg(str_platforms, "platforms")) {
		usage();
		return 1;
	}
	vec_platforms = StringTokenizer(str_platforms,',').split();
	if (!cl.getArg(str_sites, "sites")) {
		usage();
		return 1;
	}
	vec_sites = StringTokenizer(str_sites,',').split();
	if (!cl.getArg(str_samples, "samples")) {
		usage();
		return 1;
	}
	vec_samples = StringTokenizer(str_samples,',').split();
	if (!cl.getArg(str_replicates, "replicates")) {
		usage();
		return 1;
	}
	vec_replicates = StringTokenizer(str_replicates,',').split();
	if (!cl.getArg(str_flowlanes, "flowlanes")) {
		usage();
		return 1;
	}
	vec_flowlanes = StringTokenizer(str_flowlanes,',').split();

	// open database connection
	mysql_library_init(0, NULL, NULL);
	MYSQL* mysql = mysql_init(NULL);
	if (mysql == NULL) {
		cerr << "!error: seqc_extract_expression: cannot init mysql" << endl;
		mysql_library_end();
		return 1;
	}
	if (!mysql_real_connect(mysql, str_db_host.c_str(), str_db_user.c_str(), str_db_pass.c_str(), str_db_name.c_str(), 0, NULL, 0)) {
		cerr << "!error: seqc_extract_expression: cannot connect to mysql database" << endl;
		mysql_library_end();
		return 1;
	}

	// get sequence information
	map<string,map<string,vector<struct sequence> > > map_sequence_data;
	map<string,int> map_flowlane_count;
	if (!get_sequence_data(
		mysql,
		vec_platforms,
		vec_sites,
		vec_samples,
		vec_replicates,
		vec_flowlanes,
		map_sequence_data,
		map_flowlane_count
	)) {
		cerr << "!error: seqc_extract_expression: cannot get sequence data" << endl;
		mysql_close(mysql);
		mysql_library_end();
		return 1;
	}

	// get mapping pipeline stage id
	string str_mapping_pipeline_stage_id;
	if (!get_mapping_pipeline_stage_id(
		mysql,
		str_mapping_pipeline,
		str_mapping_pipeline_stage,
		str_mapping_pipeline_stage_id
	)) {
		cerr << "!error: seqc_extract_expression: cannot get mapping pipeline stage id" << endl;
		mysql_close(mysql);
		mysql_library_end();
		return 1;
	}


	// verify that the correct number of sequences were found per flowlane combination
	// 	and get mapping, and quantification ids
	map<string,map<string,vector<string> > > map_mapping_data;
	map<string,map<string,string> > map_quantification_data;

	int i_num_columns = 0;
	
	map<string,map<string,vector<struct sequence> > >::iterator it;
	for (it=map_sequence_data.begin(); it!=map_sequence_data.end(); it++) {
		//cout << it->first << endl;
		map<string,vector<struct sequence> >::iterator it2;
		for (it2=it->second.begin(); it2!=it->second.end(); it2++) {
			//cout << "\t" << it2->first << endl;
			// check count of sequences
			if (it2->second.size() != map_flowlane_count[it2->first]) {
				cerr << "!error: seqc_extract_expression: invalid number of sequences found for flowlane combination" << endl;
				mysql_close(mysql);
				mysql_library_end();
				return 1;
			}
			vector<string> vec_mapping_ids(it2->second.size());

			// get mapping ids
			for (int i=0; i<it2->second.size(); i++) {
				string str_mapping_id;
				if (!get_mapping_id(
					mysql,
					it2->second[i].str_name,
					str_mapping_pipeline_stage_id,
					str_mapping_id
				)) {
					cerr << "!error: seqc_extract_expression: cannot get mapping id" << endl;
					mysql_close(mysql);
					mysql_library_end();
					return 1;
				}
				vec_mapping_ids[i] = str_mapping_id;
			}
			map_mapping_data[it->first][it2->first] = vec_mapping_ids;

			// get quantification id
			string str_quantification_id;
			if (!get_quantification_id(
				mysql,
				str_ref_type,
				str_algorithm,
				str_normalization,
				vec_mapping_ids,
				str_quantification_id
			)) {
				cerr << "!error: seqc_extract_expression: cannot get quantification id" << endl;
				mysql_close(mysql);
				mysql_library_end();
				return 1;
			}
			map_quantification_data[it->first][it2->first] = str_quantification_id;
			i_num_columns++;
			//cout << "quantifcation id: " << str_quantification_id << endl;
		}
	}

	vector<string> vec_feature_list;
	if (!get_ordered_feature_list(
		mysql,
		str_ref_type,
		vec_feature_list
	)) {
		cerr << "!error: seqc_extract_expression: cannot get feature list" << endl;
		mysql_close(mysql);
		mysql_library_end();
		return 1;
	}

	//vector<string> vec_column_headers(i_num_columns);
	Matrix<string> mat_column_headers(6, i_num_columns+1);
	mat_column_headers[1][0] = "Platform";
	mat_column_headers[2][0] = "Site";
	mat_column_headers[3][0] = "Sample";
	mat_column_headers[4][0] = "Replicate";
	mat_column_headers[5][0] = "Flowcells/Lanes";

	// download data
	int i_cur_column = 0;
	Matrix<float> mat_data(vec_feature_list.size(), i_num_columns, 0);
	map<string,map<string,string> >::iterator qit;
	for (qit=map_quantification_data.begin(); qit!=map_quantification_data.end(); qit++) {
		map<string,string>::iterator qit2;
		for (qit2=qit->second.begin(); qit2!=qit->second.end(); qit2++) {
			if (!get_sample(
				mysql,
				qit2->second,
				mat_data,
				i_cur_column
			)) {
				cerr << "!error: seqc_extract_expression: cannot get sample" << endl;
				mysql_close(mysql);
				mysql_library_end();
				return 1;
			}
			//vec_column_headers[i_cur_column] = qit->first+string("|")+qit2->first;
			mat_column_headers[0][i_cur_column+1] = qit->first+string("|")+qit2->first;
			struct sequence seq = map_sequence_data[qit->first][qit2->first][0];
			mat_column_headers[1][i_cur_column+1] = seq.str_platform_id;
			mat_column_headers[2][i_cur_column+1] = seq.str_site;
			mat_column_headers[3][i_cur_column+1] = seq.str_sample;
			mat_column_headers[4][i_cur_column+1] = seq.str_replicate;
			mat_column_headers[5][i_cur_column+1] = qit2->first;
			i_cur_column++;
		}
	}

	DataFileHeader<string> dfh_headers(str_file, '\t');
	dfh_headers.writeDataFile(mat_column_headers, false, false);

	DataFileHeader<float> dfh(str_file, '\t');
	dfh.setRowHeaders(vec_feature_list);
	dfh.writeDataFile(mat_data, false, true, ofstream::out|ofstream::app);

	// close database connection
	mysql_close(mysql);
	mysql_library_end();

	return 0;
}
