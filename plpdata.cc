#include "plpdata.h"

Pileupdata::Pileupdata(std::string filename, std::string refname, std::string region) : plp(filename, refname, region), data() {
	populate_data();
}

Pileupdata(std::string filename, std::string refname) : plp(filename, refname), data() {
	populate_data();
}

Pileupdata(Pileup p) plp(p), data() {
	populate_data();
}

void Pileupdata::populate_data(){
	while(plp.next() != 0){
		data[plp.get_tid()][plp.get_pos()] = std::make_tuple(plp.alleles,plp.qual,plp.readgroups);
	}
}

vector<char> Pileupdata::bases_at(int tid, int pos){
	return std::get<0>(data[tid][pos]);
}

int Pileupdata::depth_at(int tid, int pos){
	return Pileupdata::bases_at(tid,pos).size();
}

int Pileupdata::num_base(int tid, int pos, char base){
	std::map<char, int> counts{{"A",0},{"T",0},{"G",0},{"C",0}};
	vector<char> alleles = bases_at(tid,pos);
	for(auto it = alleles.begin(); it != alleles.end(); ++it){
		counts[*it] += counts[*it];
	}
}

std::map<std::string,int> Pileupdata::get_name_map(){
	return plp.get_name_map();
}
