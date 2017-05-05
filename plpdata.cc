#include "plpdata.h"
#include "iostream"

Pileupdata::Pileupdata(std::string filename, std::string refname, std::string region) : plp(filename, refname, region), data(1, std::vector<pileuptuple_t>(1)) {
	populate_data();
}

Pileupdata::Pileupdata(std::string filename, std::string refname) : plp(filename, refname), data(1, std::vector<pileuptuple_t>(1)) {
	populate_data();
}

Pileupdata::Pileupdata(Pileup p) : plp(p), data() {
	populate_data();
}

void Pileupdata::populate_data(){
	std::map<std::string,int> chrs = plp.get_name_map();
	data.reserve(chrs.size());
	while(plp.next() != 0){
		int tid = plp.get_tid();
		int pos = plp.get_pos();
		data[tid].reserve(pos + 1);
		std::cout << "tid is: " << tid << " and pos is: " << pos << std::endl;
		data[tid][pos] = std::make_tuple(plp.alleles,plp.counts,plp.qual,plp.readgroups);
		std::cout << "made it out alive" << std::endl;
	}
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

