#ifndef __MEEP_PLPDATA_INCLUDED__
#define __MEEP_PLPDATA_INCLUDED__

//class for slurping in pileup data
class Pileupdata{
protected:
	Pileup plp;
	vector<vector<std::tuple<vector<char>,vector<char>,vector<std::string>>>> data; //data[tid][pos] = (bases, qualities, readgroups)
	void populate_data();
public:
	vector<char> bases_at(int tid, int pos);
	int depth_at(int tid, int pos);
	int num_base(int tid, int pos, char base);
	vector<vector<std::tuple<vector<char>,vector<char>,vector<std::string>>>> get_data();
	std::map<std::string,int> get_name_map();
	Pileupdata(std::string filename, std::string refname, std::string region);
	Pileupdata(std::string filename, std::string refname);
	Pileupdata(Pileup p);
};




#endif
