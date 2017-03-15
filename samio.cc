#include "samio.h"
#include <stdexcept>

SamReader::SamReader() : in(), header(), idx(), iter(), region(), region_exists(){
	open("-");
	this->region_exists = false;
}

SamReader::SamReader(std::string filename) : in(), header(), idx(), iter(), region(), region_exists(){
	open(filename);
	this->region_exists = false;
}

SamReader::SamReader(std::string filename, std::string region) : in(), header(), idx(), iter(), region(), region_exists(){
	open(filename, region);
	this->region_exists = true;
}

SamReader::~SamReader(){
	if (in != nullptr){
		sam_close(in);
	}
	if (idx != nullptr){
		hts_idx_destroy(idx);
	}
	if (iter != nullptr){
		hts_itr_destroy(iter);
	}
	if (header != nullptr){
		bam_hdr_destroy(header);
	}
}

//set in and header
void SamReader::open(std::string filename_in){
	htsFile* htsin = NULL;
	bam_hdr_t* htsheader = NULL;

	htsin = hts_open(filename_in.c_str(), "r");
	if (htsin == NULL) {
    	//error
        throw std::runtime_error("error, fail to open");
    }
    this->in = htsin;

    htsheader = sam_hdr_read(htsin);
    if (htsheader == NULL || htsheader->n_targets == 0){
    	//error
    	throw std::runtime_error("error, no header");
    }
    this->header = htsheader;
}

//set in and header and idx and iter
void SamReader::open(std::string filename_in, std::string region){
	open(filename_in);
    hts_idx_t* idx = sam_index_load(this->in, filename_in.c_str()); //&hts_idx_destroy
    if (idx != nullptr){
    	//error
    	throw std::runtime_error("error");
    }
    this->idx = idx;

    hts_itr_t* iter = sam_itr_querys(idx, header, region.c_str()); //&hts_itr_destroy
    if (iter != nullptr){
    	throw std::runtime_error("error");
    }
    this->iter = iter;
}

//updates b with next read; returns >=0 on success
int SamReader::next(bam1_t* b){
	if (this->has_region()){
		return sam_itr_next(this->in, this->iter, b);
	}
	else{
		return sam_read1(this->in, this->header, b);
	}
}

bool SamReader::has_region(){
	return this->region_exists;
}

bam_hdr_t* SamReader::get_header(){
	return this->header;
}

std::string SamReader::get_ref_name(bam1_t* b){
	return std::string(get_ref_name(b->core.tid));
}

std::string SamReader::get_ref_name(int tid){
	return std::string(header->target_name[tid]);
}

//SamWriter class
void SamWriter::check_open_success(){
	if (this->outfh == 0) {
		throw std::runtime_error("error opening file");
	}
}

SamWriter::SamWriter(){
	this->outfh = hts_open("-", "w");
	this->header = nullptr;
	check_open_success();
}

SamWriter::SamWriter(std::string filename){
	this->outfh = hts_open(filename.c_str(), "w");
	this->header = nullptr;
	check_open_success();
}

SamWriter::SamWriter(bam_hdr_t* h){
	this->outfh = hts_open("-", "w");
	this->header = h;
	check_open_success();
}

SamWriter::SamWriter(std::string filename, bam_hdr_t* h){
	this->outfh = hts_open(filename.c_str(), "w");
	this->header = h;
	check_open_success();
}

SamWriter::~SamWriter(){
	if (outfh != nullptr){
		sam_close(outfh);
	}
}

void SamWriter::write_header(){
	if (this->header == nullptr){
		throw std::runtime_error("error writing header: null header");
	}
	else{
		if (sam_hdr_write(this->outfh, this->header) != 0){
			throw std::runtime_error("error writing header");
		}
	}
}

void SamWriter::write_read(const bam1_t* b){
	if (sam_write1(this->outfh, this->header, b) < 0){
		throw std::runtime_error("error writing read");
	}
}

