#ifndef __ERRLUMINA_VCF_IO_INCLUDED__
#define __ERRLUMINA_VCF_IO_INCLUDED__

#include <htslib/vcf.h>
#include <string>

class VCFwriter{
protected:
	htsFile* outfh;
	bcf_hdr_t* header;
public:
	VCFWriter();
	VCFWriter(std::string);
	VCFWriter(std::string, bcf_hdr_t*);
	VCFWriter(bcf_hdr_t*);
	~VCFWriter();
	void check_open_success(); //void, throws error if not open
	void append_to_header(std::string s); // void, call before writing header. throws on error.
	void append_to_headef(); //void, throws on error
	void write_header(); //throws error if header is null or fails to write
	void write_variant(bcf1_t);
};



#endif