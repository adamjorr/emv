#include "pileup.h"
#include <htslib/sam.h>
#include "samio.h"

Pileup::Pileup(std::string samfile, std::string reffile): reader(samfile), ref(reffile), tid(), pos(), cov(), pileup(nullptr), iter() {
	iter = bam_plp_init(&plp_get_read, &reader);
}

Pileup::Pileup(std::string samfile, std::string reffile, std::string region): reader(samfile, region), ref(reffile), tid(), pos(), cov(), pileup(nullptr), iter() {
	iter = bam_plp_init(&plp_get_read, &reader);
}

//typedef int (*bam_plp_auto_f)(void *data, bam1_t *b);
int plp_get_read(void *data, bam1_t *b){
	SamReader *reader = data;
	reader->next(b);
}

//possible optimization: store sequence strings in a hash w/ alignment, throw out of hash once no longer in pileup
int next(){
	if((pileup = bam_plp_auto(iter, &tid, &pos, &cov)) != nullptr){ //successfully pile up new position
		alleles.clear(); qual.clear(); names.clear(); readgroups.clear(); samples.clear();
		alleles.reserve(cov); qual.reserve(cov); names.reserve(cov); readgroups.reserve(cov); samples.reserve(cov);
		for (i = 0; i < cov; i++){
			bam1_t* alignment = pileup[i].b;
			uint8_t *seq = bam_get_seq(alignment);
			int qpos = pileup[i].qpos;
			int baseint = bam_seqi(seq,qpos)
			alleles[i] = seq_nt16_str[baseint];
			qual[i] = bam_get_qual(alignment)[qpos];
			names[i] = std::string(bam_get_qname(alignment));
			readgroups[i] = bam_aux_get(alignment, "RG");
			samples[i] = bam_aux_get(alignment, "SM");
		}
		return 1;
	} else {
		return 0;
	}
}

void clear_vectors(){
	alleles.clear();
	qual.clear();
	names.clear();
}

void reserve_vectors(){

}


