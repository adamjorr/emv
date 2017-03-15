#include "mismatchfinder.h"
#include <stdexcept>

MismatchFinder::MismatchFinder(SamReader* r, std::string ref_name) : tid(-2), ref_len(-1), ref_t(ref_name), reader(r) {
}

// int len;
    // ref_file = argc > optind + 1 ? argv[optind+1] : ga.reference;
    // fai = fai_load(ref_file);
//     ref = fai_fetch(fai, header->target_name[b->core.tid], &len);

//to get ref and len:
//https://github.com/samtools/htslib/blob/develop/htslib/faidx.h
void MismatchFinder::dump_mismatches(std::string fileout){
	SamWriter w(this->reader->get_header());
	w.write_header();
	bam1_t *b = NULL;
	b = bam_init1();
	int r;

	while((r = this->reader->next(b)) >= 0){
		if (has_mismatch(b, ref_t.get_ref(reader->get_ref_name(b)), ref_t.get_ref_len())){
			w.write_read(b);
		}
	}
}

void MismatchFinder::dump_mismatches(){
	dump_mismatches("-");
}

void MismatchFinder::dump_locations(std::ostream *os){
	bam1_t *b = NULL;
	b = bam_init1();
	int r;
	int location;

	while((r = this->reader->next(b)) >= 0){
		location = mismatch_location(b,ref_t.get_ref(get_region(b)).c_str(),ref_t.get_ref_len());
		if (location >= 0){
			*os << location << "\t" << get_cigar_str(b) << std::endl;
		}
	}
}

void MismatchFinder::dump_locations(std::string fileout){
	if (fileout == "-"){
		dump_locations();
	}
	else{
		std::ofstream ofs(fileout, std::ofstream::out);
		dump_locations(&ofs);
		ofs.close();
	}
}

void MismatchFinder::dump_locations(){
	MismatchFinder::dump_locations(&std::cout);
}


//see https://github.com/samtools/samtools/blob/1bae2510c2f58e0332b84780b3c6bd438c58ed3c/bam_md.c#L58-L90
//this is simplified from samtools code, so we'll leave ref as a c_str
static bool has_mismatch(bam1_t *b, char *ref, int ref_len){
	uint8_t *seq = bam_get_seq(b);
	uint32_t *cigar = bam_get_cigar(b);
	bam1_core_t *c = &b->core;
	int i, x, y = 0;
	for (i = y = 0, x = c->pos; i < c->n_cigar; ++i){
		int j, l = cigar[i]>>4,  op = cigar[i]&0xf;
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF){
			for (j = 0; j < l; ++j){
				int c1, c2, z = y + j;
				if (x + j >= ref_len || ref[x + j] == '\0') break; //out of bounds
				c1 = bam_seqi(seq, z);
				c2 = seq_nt16_table[(int)ref[x+j]];
				if (!((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0)){ //doesn't match
					return true;
				}
			}
			if (j < l) break; //break again if we have to break because we're out of bounds
			x += l;
			y += l;
		}
		else{
			return true;
		}
	}
	return false;
}

//see https://github.com/samtools/samtools/blob/1bae2510c2f58e0332b84780b3c6bd438c58ed3c/bam_md.c#L58-L90
//this is simplified from samtools code, so we'll leave ref as a c_str
static int mismatch_location(bam1_t *b, char *ref, int ref_len){
	uint8_t *seq = bam_get_seq(b);
	uint32_t *cigar = bam_get_cigar(b);
	bam1_core_t *c = &b->core;
	int i, x, y = 0;
	for (i = y = 0, x = c->pos; i < c->n_cigar; ++i){
		int j, l = cigar[i]>>4,  op = cigar[i]&0xf;
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF){
			for (j = 0; j < l; ++j){
				int c1, c2, z = y + j;
				if (x + j >= ref_len || ref[x + j] == '\0') break; //out of bounds
				c1 = bam_seqi(seq, z);
				c2 = seq_nt16_table[(int)ref[x+j]];
				if (!((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0)){ //doesn't match
					return z;
				}
			}
			if (j < l) break; //break again if we have to break because we're out of bounds
			x += l;
			y += l;
		}
		else{
			return y;
		}
	}
	return -1;
}

std::string get_sequence(bam1_t *b){
	uint8_t *seq = bam_get_seq(b);
	int32_t seqlen = b->core.l_qseq;
	int baseint;
	char basech;
	std::string s;

	for (int i = 0; i < seqlen; i++){
		baseint = bam_seqi(seq, i);
		basech = seq_nt16_str[baseint];
		s.push_back(basech);
	}
	return s;
}

std::string get_cigar_str(bam1_t *b){
	uint32_t *cigar = bam_get_cigar(b);
	uint32_t cigarlen = b->core.n_cigar;
	std::string s;

	for (int i = 0; i < cigarlen; i++){
		uint32_t c = cigar[i];
		int len = bam_cigar_oplen(c);
		char op = bam_cigar_opchr(c);
		s += std::to_string(len);
		s.push_back(op);
	}
	return s;
}



