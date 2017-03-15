#ifndef __ERRLUMINA_REFTYPE_INCLUDED__
#define __ERRLUMINA_REFTYPE_INCLUDED__

class Reftype{
protected:
	faidx_t* faidx_p;
	std::string ref;
	int ref_len;
	std::string region;
public:
	Reftype(std::string reference_name);
	Reftype(faidx_t* faidx_p);
	~Reftype();
	std::string get_ref(std::string region); //update ref if necessary, otherwise do nothing. then return ref. throws.
	int get_ref_len();
};



#endif