#ifndef __MEEP_PILEUP_INCLUDED__
#define __MEEP_PILEUP_INCLUDED__

#include "samio.h"
#include "reftype.h"
#include <htslib/sam.h>
#include <string>
#include <vector>

class Pileup{
protected:
	SamReader reader;
	Reftype ref;
	int tid;
	int pos;
	int cov;
	const bam_pileup1_t *pileup;
	bam_plp_t iter;
public:
	Pileup(std::string samfile, std::string reffile);
	Pileup(std::string samfile, std::string reffile, std::string region);
	Pileup();
	~Pileup();
	std::vector<char> alleles;
	std::vector<char> qual;
	std::vector<std::string> names;
	std::vector<std::string> readgroups;
	std::map<char,int> counts;
	char ref_char;
	static int plp_get_read(void *data, bam1_t *b);
	int next();
	int get_tid();
	int get_pos();
	int get_ref_tid(std::string name);
	std::map<std::string,int> get_name_map();
	std::string get_chr_name(int tid);
};


#endif
