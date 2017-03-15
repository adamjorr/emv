#ifndef __ERRLUMINA_VCALLER_INCLUDED__
#define __ERRLUMINA_VCALLER_INCLUDED__

class VCaller{
protected:
	SamReader reader;
	Reftype ref;
	VCFwriter writer;
	std::vector<int> cutoffs;
	long double seq_error;
	long double cutoff_error;
public:
	VCaller(std::string samfile, std::string reffile, std::string outfile);
	std::vector<char> call_variants();
	void pileup_and_call();
	int get_cutoff(int coverage);
};









#endif