#ifndef __MEEP_PLPDATA_INCLUDED__
#define __MEEP_PLPDATA_INCLUDED__

#include <tuple>
#include <vector>
#include <map>
#include "pileup.h"

typedef std::tuple<std::vector<char>,std::map<char,int>,std::vector<char>,char,std::vector<std::string> > pileuptuple_t; //(bases, counts, qualities, readgroups)
typedef std::vector<std::vector<pileuptuple_t> > pileupdata_t; //data[tid][pos] = (bases, counts, qualities, ref, readgroups)


//class for slurping in pileup data
class Pileupdata{
protected:
	Pileup plp;
	pileupdata_t data; //data[tid][pos] = (bases, counts, qualities, ref, readgroups)
	std::map<char,int> ref_counts;
	void populate_data();
	void populate_data(std::vector<char> x, char ref, std::vector<char> quals);
public:
	std::vector<char> bases_at(int tid, int pos);
	int depth_at(int tid, int pos);
	int num_base(int tid, int pos, char base);
	pileupdata_t get_data();
	std::map<std::string,int> get_name_map();
	std::map<char,int> get_ref_counts();
	Pileupdata(std::string filename, std::string refname, std::string region);
	Pileupdata(std::string filename, std::string refname);
	Pileupdata(Pileup p);
	Pileupdata(std::vector<char> x, char ref, std::vector<char> quals);
	Pileupdata(std::vector<char> x);
};




#endif
