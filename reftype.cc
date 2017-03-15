#include "reftype.h"

Reftype::Reftype(std::string reference_name) : {
	faidx_t* faidx = fai_load(reference_name.c_str());
	if (faidx == nullptr){
		throw std::runtime_error("error loading reference");
	}
	faidx_p = faidx;
}

Reftype::Reftype(faix_t* faidx_p) : faidx_p(faidx_p) {
}

Reftype::~Reftype(){
	if (faidx_p != nullptr){
		fai_destroy(faidx_p);
	}
	if (ref != nullptr){
		free(ref);
	}
}

std::string Reftype::get_ref(std::string region){
	if (region != this->region){
		this->region = region;
		ref = std::string(fai_fetch_seq(faidx_p,region.c_str(),&ref_len));
		if (ref == nullptr){
			throw std::runtime_error("error getting ref");
		}
	}
	return ref;
}

int Reftype::get_ref_len(){
	return ref_len;
}



