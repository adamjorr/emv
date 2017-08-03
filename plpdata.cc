#include "plpdata.h"
#include <iostream>
#include <vector>

Pileupdata::Pileupdata(std::string filename, std::string refname, std::string region) : plp(filename, refname, region), data() {
	populate_data();
}

Pileupdata::Pileupdata(std::string filename, std::string refname) : plp(filename, refname), data() {
	populate_data();
}

Pileupdata::Pileupdata(Pileup p) : plp(p), data() {
	populate_data();
}

Pileupdata::Pileupdata(std::vector<char> x, char ref, std::vector<char> quals) : plp(), data(){
	populate_data(x,ref,quals);
}

Pileupdata::Pileupdata(std::vector<char> x) : Pileupdata(x, x[0], x) {
}

//I'll need to think of something better; this will break if the pileup isn't completely contiguous (ie multiple ranges)
void Pileupdata::populate_data(){
	std::map<std::string,int> chrs = plp.get_name_map();
	int val;
	while((val = plp.next()) != 0){
		if (val == 1){
			int tid = plp.get_tid();
			char ref_char = plp.ref_char;		
			data.resize(tid + 1);
			data[tid].push_back(std::make_tuple(plp.alleles,plp.counts,plp.qual,ref_char,plp.readgroups));
			++ref_counts[ref_char];
		}
	}
}

void Pileupdata::populate_data(std::vector<char> x, char ref, std::vector<char> quals){
	std::map<char,int> counts = {{'A',0},{'T',0},{'G',0},{'C',0}};
	std::vector<std::string> rgs(x.size(),"RG0");
	for (char i : x){
		++counts[i];
	}
	data.resize(1);
	std::vector<pileuptuple_t> v(1,std::make_tuple(x,counts,quals,ref,rgs));
	data.push_back(v);
	++ref_counts[ref];
}


std::vector<char> Pileupdata::bases_at(int tid, int pos){
	return std::get<0>(data[tid][pos]);
}

int Pileupdata::depth_at(int tid, int pos){
	return Pileupdata::bases_at(tid,pos).size();
}

int Pileupdata::num_base(int tid, int pos, char base){
	return std::get<1>(data[tid][pos])[base];
}

std::map<std::string,int> Pileupdata::get_name_map(){
	return plp.get_name_map();
}

pileupdata_t Pileupdata::get_data(){
	return data;
}

