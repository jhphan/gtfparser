#ifndef __TRANSCRIPT_H__
#define __TRANSCRIPT_H__

#include <iostream>
#include <vector>
#include <string>

using namespace std;

class Transcript;

ostream& operator<< (ostream& ostr, const Transcript& trans);

class Transcript {

private:
	string _str_name;
	string _str_gene_name;
	string _str_chromosome;	
	int _i_strand; 		// +1 or -1
	unsigned long _ul_start;	// transcript start coordinate
	unsigned long _ul_end;	// transcript end coordinate
	int _i_num_exons;
	
	vector<unsigned long> _vec_ul_exon_start;	// genomic start coord
	vector<unsigned long> _vec_ul_exon_end;	// genomic end coord
	vector<unsigned long> _vec_ul_exon_rel_start;	// start coord relative to transcript
	vector<unsigned long> _vec_ul_exon_rel_end;	// end coord relative to transcript
	
	void copy(const Transcript& trans);

public:

	Transcript(string str_name = "transcript", string str_gene_name = "gene", string str_chromosome = "chr1", int i_strand = 1);
	Transcript(const Transcript& trans);

	Transcript& operator=(const Transcript& trans);

	string get_name() const { return _str_name; };
	string get_gene_name() const { return _str_gene_name; };
	string get_chromosome() const { return _str_chromosome; };
	int get_strand() const { return _i_strand; };
	unsigned long get_start() const { return _ul_start; };
	unsigned long get_end() const { return _ul_end; };
	int get_num_exons() const { return _i_num_exons; };

	void set_name(string str_name);
	void set_gene_name(string str_gene_name);
	void set_chromosome(string str_chromosome);
	void set_strand(int i_strand);
	void set_start(unsigned long ul_start);
	void set_end(unsigned long ul_end);

	void clear_exons();
	bool add_exon(unsigned long ul_exon_start, unsigned long ul_exon_end);
	bool set_exon_start(int i_exon_num, unsigned long ul_exon_start);
	bool set_exon_end(int i_exon_num, unsigned long ul_exon_end);
	bool set_exon_rel_start(int i_exon_num, unsigned long ul_exon_rel_start);
	bool set_exon_rel_end(int i_exon_num, unsigned long ul_exon_rel_end);
	bool set_exon_coords(int i_exon_num, unsigned long ul_exon_start, unsigned long ul_exon_end, unsigned long ul_exon_rel_start, unsigned long ul_exon_rel_end);

	unsigned long get_exon_start(int i_exon_num) const;
	unsigned long get_exon_end(int i_exon_num) const;
	unsigned long get_exon_rel_start(int i_exon_num) const;
	unsigned long get_exon_rel_end(int i_exon_num) const;

	bool translate_to_genomic_coord(unsigned long ul_trans_coord, int& i_exon_number, unsigned long& ul_genomic_coord);

	unsigned long get_length();

	~Transcript();

	friend ostream& operator<< (ostream& ostr, const Transcript& trans);
};

#endif
