#ifndef __ERRLUMINA_HTS_API_INCLUDED__
#define __ERRLUMINA_HTS_API_INCLUDED__

#include <memory>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <string>

class SamReader{
protected:
	htsFile* in;
	bam_hdr_t* header;
	hts_idx_t* idx;
	hts_itr_t* iter;
	std::string region;
	bool region_exists;
	void open(std::string);
	void open(std::string, std::string);
public:
	bool has_region();
	bam_hdr_t* get_header();
	int next(bam1_t *b);
	std::string get_ref_name(bam1_t* b);
	std::string get_ref_name(int tid);
	SamReader();
	SamReader(const std::string filename);
	SamReader(const std::string filename, const std::string region);
	~SamReader();
};

class SamWriter{
protected:
	htsFile* outfh;
	bam_hdr_t* header;
public:
	SamWriter(std::string, bam_hdr_t*);
	SamWriter(std::string);
	SamWriter(bam_hdr_t*);
	SamWriter();
	~SamWriter();
	void check_open_success(); //void, throws error if not open
	void write_header(); //throws error if header is null or fails to write
	void write_read(const bam1_t*);
};

#endif