#ifndef __MEEP_PLPDATA_INCLUDED__
#define __MEEP_PLPDATA_INCLUDED__

typedef pileuptuple_t std::tuple<vector<char>,std::map<char,int>,vector<char>,vector<std::string>> //(bases, counts, qualities, readgroups)
typedef vector<vector<pileuptuple_t>> pileupdata_t; //data[tid][pos] = (bases, counts, qualities, readgroups)


//class for slurping in pileup data
class Pileupdata{
protected:
	Pileup plp;
	pileupdata_t data; //data[tid][pos] = (bases, counts, qualities, readgroups)
	void populate_data();
public:
	vector<char> bases_at(int tid, int pos);
	int depth_at(int tid, int pos);
	int num_base(int tid, int pos, char base);
	pileupdata_t get_data();
	std::map<std::string,int> get_name_map();
	Pileupdata(std::string filename, std::string refname, std::string region);
	Pileupdata(std::string filename, std::string refname);
	Pileupdata(Pileup p);
};




#endif
