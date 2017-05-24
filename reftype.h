#ifndef __ERRLUMINA_REFTYPE_INCLUDED__
#define __ERRLUMINA_REFTYPE_INCLUDED__

#include <htslib/faidx.h>
#include <string>

class Reftype{
protected:
	faidx_t* faidx_p;
	char* ref_p;
	std::string ref;
	int ref_len;
	std::string region;
public:
	Reftype(std::string reference_name);
	Reftype(faidx_t* faidx_p);
	Reftype();
	~Reftype();
	std::string get_ref(std::string region); //update ref if necessary, otherwise do nothing. then return ref. throws.
	int get_ref_len();
};



#endif