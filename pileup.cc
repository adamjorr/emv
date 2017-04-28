#include "pileup.h"
#include <htslib/sam.h>
#include "samio.h"

Pileup::Pileup(std::string samfile, std::string reffile): reader(samfile), ref(reffile), tid(), pos(), cov(), pileup(nullptr), iter(), alleles(), qual(), names(), readgroups(), samples(), counts({{'A',0},{'T',0},{'G',0},{'C',0}})  {
	iter = bam_plp_init(&Pileup::plp_get_read, &reader);
}

Pileup::Pileup(std::string samfile, std::string reffile, std::string region): reader(samfile, region), ref(reffile), tid(), pos(), cov(), pileup(nullptr), iter() {
	iter = bam_plp_init(&Pileup::plp_get_read, &reader);
}

Pileup::~Pileup(){
	bam_plp_destroy(iter);
}

//typedef int (*bam_plp_auto_f)(void *data, bam1_t *b);
int Pileup::plp_get_read(void *data, bam1_t *b){
	SamReader *reader = (SamReader*)data;
	reader->next(b);
}

//possible optimization: store sequence strings in a hash w/ alignment, throw out of hash once no longer in pileup
int Pileup::next(){
	if((pileup = bam_plp_auto(iter, &tid, &pos, &cov)) != nullptr){ //successfully pile up new position
		alleles.clear(); qual.clear(); names.clear(); readgroups.clear(); samples.clear(); counts.clear();
		alleles.reserve(cov); qual.reserve(cov); names.reserve(cov); readgroups.reserve(cov); samples.reserve(cov);
		for (int i = 0; i < cov; i++){
			bam1_t* alignment = pileup[i].b;
			uint8_t *seq = bam_get_seq(alignment);
			int qpos = pileup[i].qpos;
			int baseint = bam_seqi(seq,qpos);
			char allele = seq_nt16_str[baseint];
			alleles[i] = allele;
			++counts[allele];
			qual[i] = bam_get_qual(alignment)[qpos];
			names[i] = std::string(bam_get_qname(alignment));
			readgroups[i] = std::to_string(*bam_aux_get(alignment, "RG"));
			samples[i] = std::to_string(*bam_aux_get(alignment, "SM"));
		}
		return 1;
	} else {
		return 0;
	}
}

int Pileup::get_tid(){
	return tid;
}

int Pileup::get_pos(){
	return pos;
}

std::string Pileup::chr_name(){
	return reader.get_ref_name(tid);
}

int Pileup::get_ref_tid(std::string name){
	return reader.get_ref_tid(name);
}

std::map<std::string,int> Pileup::get_name_map(){
	return reader.get_name_map();
}

