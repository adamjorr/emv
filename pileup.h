#ifndef __MEEP_PILEUP_INCLUDED__
#define __MEEP_PILEUP_INCLUDED__

#include "samio.h"
#include "reftype.h"
#include <htslib/sam.h>
#include <string>

class Pileup{
protected:
	SamReader reader;
	Reftype ref;
	int tid;
	int pos;
	int cov;
	bam_pileup1_t *pileup;
	bam_plp_t iter;
public:
	Pileup(std::string samfile, std::string reffile);
	Pileup(std::string samfile, std::string reffile, std::string region);
	~Pileup();
	std::vector<char> alleles;
	std::vector<char> qual;
	std::vector<std::string> names;
	std::vector<std::string> readgroups;
	std::vector<std::string> samples;
	int plp_get_read(void *data, bam1_t *b);
	int next();
};



#endif
