#include "../include/transcript.h"

Transcript::Transcript(string str_name, string str_gene_name, string str_chromosome, int i_strand) {
	set_name(str_name);
	set_gene_name(str_gene_name);
	set_chromosome(str_chromosome);
	set_strand(i_strand);
	set_start(0);
	set_end(0);
	clear_exons();
}

Transcript::Transcript(const Transcript& trans) {
	copy(trans);
}

void Transcript::copy(const Transcript& trans) {
	set_name(trans.get_name());
	set_gene_name(trans.get_gene_name());
	set_chromosome(trans.get_chromosome());
	set_strand(trans.get_strand());
	set_start(trans.get_start());
	set_end(trans.get_end());

	_i_num_exons = trans.get_num_exons();
	_vec_ul_exon_start = vector<unsigned long>(_i_num_exons, 0);
	_vec_ul_exon_end = vector<unsigned long>(_i_num_exons, 0);
	_vec_ul_exon_rel_start = vector<unsigned long>(_i_num_exons, 0);
	_vec_ul_exon_rel_end = vector<unsigned long>(_i_num_exons, 0);

	for (int i=0; i<_i_num_exons; i++) {
		_vec_ul_exon_start[i] = trans.get_exon_start(i);
		_vec_ul_exon_end[i] = trans.get_exon_end(i);
		_vec_ul_exon_rel_start[i] = trans.get_exon_rel_start(i);
		_vec_ul_exon_rel_end[i] = trans.get_exon_rel_end(i);
	}
}

Transcript& Transcript::operator=(const Transcript& trans) {
	copy(trans);
	return *this;
}

void Transcript::set_name(string str_name) {
	_str_name = str_name;
}

void Transcript::set_gene_name(string str_gene_name) {
	_str_gene_name = str_gene_name;
}

void Transcript::set_chromosome(string str_chromosome) {
	_str_chromosome = str_chromosome;
}

void Transcript::set_strand(int i_strand) {
	_i_strand = i_strand;
	if (_i_strand == 0) _i_strand = 1;
}

void Transcript::set_start(unsigned long ul_start) {
	_ul_start = ul_start;
}

void Transcript::set_end(unsigned long ul_end) {
	_ul_end = ul_end;
}

void Transcript::clear_exons() {
	_i_num_exons = 0;

	_vec_ul_exon_start = vector<unsigned long>();
	_vec_ul_exon_end = vector<unsigned long>();
	_vec_ul_exon_rel_start = vector<unsigned long>();
	_vec_ul_exon_rel_end = vector<unsigned long>();	
}


bool Transcript::add_exon(
	unsigned long ul_exon_start,
	unsigned long ul_exon_end
) {
	if (ul_exon_start > ul_exon_end) {
		cerr << "!error: Transcript::add_exon(): start and end coordinates out of order" << endl;
		return false;
	}

	if (_i_num_exons == 0) {
		// first exon
		_i_num_exons = 1;
		_vec_ul_exon_start.push_back(ul_exon_start);
		_vec_ul_exon_end.push_back(ul_exon_end);
		_vec_ul_exon_rel_start.push_back(1);
		_vec_ul_exon_rel_end.push_back(ul_exon_end-ul_exon_start+1);
		return true;
	}
	int i_cur_exon_num = _i_num_exons-1;
	if (ul_exon_start < _vec_ul_exon_start[i_cur_exon_num]) {
		cerr << "!error: Transcript::add_exon(): exons must be added in order by start coord" << endl;
		return false;
	}
	if (ul_exon_start <= _vec_ul_exon_end[i_cur_exon_num]) {
		// overlapping exon
		if (ul_exon_end <= _vec_ul_exon_end[i_cur_exon_num]) {
			// exon falls completely within previous exon
			// don't need to do anything
			return true;
		} else {
			// merge the exons by extending previous exon end
			_vec_ul_exon_rel_end[i_cur_exon_num] += (ul_exon_end-_vec_ul_exon_end[i_cur_exon_num]);
			_vec_ul_exon_end[i_cur_exon_num] = ul_exon_end;
			return true;
		}
	} else {
		// new exon
		_i_num_exons++;
		_vec_ul_exon_start.push_back(ul_exon_start);
		_vec_ul_exon_end.push_back(ul_exon_end);
		int i_rel_start = _vec_ul_exon_rel_end[i_cur_exon_num]+1;
		int i_rel_end = i_rel_start+ul_exon_end-ul_exon_start;
		_vec_ul_exon_rel_start.push_back(i_rel_start);
		_vec_ul_exon_rel_end.push_back(i_rel_end);
	}
	return true;
}

bool Transcript::set_exon_start(int i_exon_num, unsigned long ul_exon_start) {
	if (i_exon_num >= _i_num_exons) {
		cerr << "!error: Transcript::set_exon_start(): exon number out of range" << endl;
		return false;
	}
	_vec_ul_exon_start[i_exon_num] = ul_exon_start;
	return true;
}

bool Transcript::set_exon_end(int i_exon_num, unsigned long ul_exon_end) {
	if (i_exon_num >= _i_num_exons) {
		cerr << "!error: Transcript::set_exon_end(): exon number out of range" << endl;
		return false;
	}
	_vec_ul_exon_end[i_exon_num] = ul_exon_end;
	return true;
}

bool Transcript::set_exon_rel_start(int i_exon_num, unsigned long ul_exon_rel_start) {
	if (i_exon_num >= _i_num_exons) {
		cerr << "!error: Transcript::set_exon_rel_start(): exon number out of range" << endl;
		return false;
	}
	_vec_ul_exon_rel_start[i_exon_num] = ul_exon_rel_start;
	return true;
}

bool Transcript::set_exon_rel_end(int i_exon_num, unsigned long ul_exon_rel_end) {
	if (i_exon_num >= _i_num_exons) {
		cerr << "!error: Transcript::set_exon_rel_end(): exon number out of range" << endl;
		return false;
	}
	_vec_ul_exon_rel_end[i_exon_num] = ul_exon_rel_end;
	return true;
}
	
bool Transcript::set_exon_coords(
	int i_exon_num,
	unsigned long ul_exon_start,
	unsigned long ul_exon_end,
	unsigned long ul_exon_rel_start,
	unsigned long ul_exon_rel_end
) {
	if (i_exon_num >= _i_num_exons) {
		cerr << "!error: Transcript::set_exon_coords(): exon number out of range" << endl;
		return false;
	}
	_vec_ul_exon_start[i_exon_num] = ul_exon_start;
	_vec_ul_exon_end[i_exon_num] = ul_exon_end;
	_vec_ul_exon_rel_start[i_exon_num] = ul_exon_rel_start;
	_vec_ul_exon_rel_end[i_exon_num] = ul_exon_rel_end;
	return true;
}

unsigned long Transcript::get_exon_start(int i_exon_num) const {
	if (i_exon_num >= _i_num_exons) return 0;
	return _vec_ul_exon_start[i_exon_num];
}

unsigned long Transcript::get_exon_end(int i_exon_num) const {
	if (i_exon_num >= _i_num_exons) return 0;
	return _vec_ul_exon_end[i_exon_num];
}

unsigned long Transcript::get_exon_rel_start(int i_exon_num) const {
	if (i_exon_num >= _i_num_exons) return 0;
	return _vec_ul_exon_rel_start[i_exon_num];
}

unsigned long Transcript::get_exon_rel_end(int i_exon_num) const {
	if (i_exon_num >= _i_num_exons) return 0;
	return _vec_ul_exon_rel_end[i_exon_num];
}


// given a transcriptome map position, get the corresponding exon number and 
// genomic coord
bool Transcript::translate_to_genomic_coord(unsigned long ul_trans_coord, int& i_exon_number, unsigned long& ul_genomic_coord) {
	int i_cur_exon_number = 0;

	// loop through exons and check if trans_coord falls within relative start/end pos
	while (i_cur_exon_number < _i_num_exons) {
		if ( 
			(_vec_ul_exon_rel_start[i_cur_exon_number] <= ul_trans_coord)
			&& (_vec_ul_exon_rel_end[i_cur_exon_number] >= ul_trans_coord)
		) {
			break;
		}
		i_cur_exon_number++;
	}

	if (i_cur_exon_number == _i_num_exons) {
		// trans_coord is outside of transcript range or in intronic region
		return false;
	}

	i_exon_number = i_cur_exon_number;

	ul_genomic_coord = _vec_ul_exon_start[i_exon_number]+(ul_trans_coord-_vec_ul_exon_rel_start[i_exon_number]);

	return true;
}

unsigned long Transcript::get_length() {
	unsigned long ul_length = 0;
	for (int i=0; i<_i_num_exons; i++)
		ul_length += _vec_ul_exon_rel_end[i]-_vec_ul_exon_rel_start[i]+1;
	return ul_length;
}

Transcript::~Transcript() {
}

ostream& operator<< (ostream& ostr, const Transcript& trans) {
	ostr << trans.get_name() << ":" << trans.get_gene_name() << ":" << trans.get_chromosome() << ":";
	if (trans.get_strand() > 0) {
		ostr << "+:";
	} else {
		ostr << "-:";
	}
	ostr << trans.get_start() << "-" << trans.get_end() << ":";
	ostr << trans.get_num_exons() << endl;
	int i_num_exons = trans.get_num_exons();
	for (int i=0; i<i_num_exons; i++) {
		ostr << "\t" 
			<< trans.get_exon_start(i) 
			<< " (" << trans.get_exon_rel_start(i) << ") - "
			<< trans.get_exon_end(i)
			<< " (" << trans.get_exon_rel_end(i) << ")" << endl;
	}
	return ostr;
}
