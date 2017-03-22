#ifndef __ERRLUMINA_VCALLER_INCLUDED__
#define __ERRLUMINA_VCALLER_INCLUDED__

class VCaller{
protected:
	SamReader reader;
	Reftype ref;
	VCFwriter writer;
	std::vector<int> cutoffs;
	long double seq_error;
	long double cutoff_error;
public:
	VCaller(std::string samfile, std::string reffile, std::string outfile);
	std::vector<char> call_variants(std::vector<char>);
	void pileup_and_call();
	static long double binomial_cdf(int successes, int trials, long double p);
	static long double binomial_pdf(int successes, int trials, long double p);
	static int binomial_coeff(int n, int k);
	static long double p_incorrect_call(int count, int coverage, long double error);
	static int p_threshold(int coverage, long double error, long double cutoff_error);
	int get_cutoff(int coverage);
	std::vector<char> get_piled_alleles(bam_pileup1_t* pileup, int cov);
	bcf1_t* create_vcf_line(bcf_hdr_t *hdr, char ref, std::vector<char> variants);
	int plp_get_read(void *data, bam1_t *b);
};









#endif