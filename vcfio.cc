#include "vcfio.h"
#include <stdexcept>

// VCFWriter::VCFWriter() : outfh(), header() {
// 	this->outfh = hts_open("-", "w");
// 	check_open_success();
// }

// VCFWriter::VCFWriter(std::string filename) : outfh(), header() {
// 	this->outfh = hts_open(filename.c_str(), "w");
// 	check_open_success();
// }

// VCFWriter::VCFWriter(std::string filename, bcf_hdr_t* h) : outfh(), header(h) {
// 	this->outfh = hts_open(filename.c_str(), "w");
// 	check_open_success();
// }

// VCFWriter::VCFWriter(bcf_hdr_t* h) : outfh(), header(h) {
// 	this->outfh = hts_open("-", "w");
// 	check_open_success();
// }

VCFWriter::VCFWriter(std::string filename, bcf_hdr_t* h) : outfh(), header(h) {
	this->outfh = hts_open(filename.c_str(), "w");
	check_open_success();
}

VCFWriter::VCFWriter(std::string filename) : VCFWriter(filename, nullptr) {
	this->header = bcf_hdr_init("w");
}

VCFWriter::VCFWriter(bcf_hdr_t* h) : VCFWriter("-", h) {}
VCFWriter::VCFWriter() : VCFWriter("-") {}


VCFWriter::~VCFWriter(){
	if (outfh != nullptr){
		hts_close(outfh);
	}
	if (header != nullptr){
		bcf_hdr_destroy(header);
	}
}

VCFWriter::check_open_success(){
	if (outfh == nullptr){
		throw std::runtime_error("error opening vcf for writing");
	}
}

VCFWriter::append_to_header(std::string s){
	if (bcf_hdr_append(header, s.c_str()) != 0){
		throw std::runtime_error("error appending to header");
	}
}

VCFWriter::write_header(){
	if (this->header == nullptr){
		throw std::runtime_error("error writing header: null header");
	}
	else{
		if (bcf_hdr_write(this->outfh, this->header) != 0){
			throw std::runtime_error("error writing header");
		}
	}
}

VCFWriter::write_variant(bcf1_t *v){
	if (bcf_write(this->outfh, this->header, v) < 0){
		throw std::runtime_error("error writing read");
	}
}






