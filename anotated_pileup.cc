//https://github.com/samtools/samtools/blob/9df632114a537f57928e4b2c17fbb1928679f877/bam_plcmd.c#L317-758

typedef struct {
    int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag, all;
    int rflag_require, rflag_filter;
    int openQ, extQ, tandemQ, min_support; // for indels
    double min_frac; // for indels
    char *reg, *pl_list, *fai_fname, *output_fname;
    faidx_t *fai;
    void *bed, *rghash;
    int argc;
    char **argv;
    sam_global_args ga;
} mplp_conf_t;

mplp_conf_t mplp;
memset(&mplp, 0, sizeof(mplp_conf_t));
mplp.min_baseQ = 13;
mplp.capQ_thres = 0;
mplp.max_depth = 250; mplp.max_indel_depth = 250;
mplp.openQ = 40; mplp.extQ = 20; mplp.tandemQ = 100;
mplp.min_frac = 0.002; mplp.min_support = 1;
mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_SMART_OVERLAPS;
mplp.argc = argc; mplp.argv = argv;
mplp.rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
mplp.output_fname = NULL;
mplp.all = 0;
sam_global_args_init(&mplp.ga);
mpileup(&mplp, argc - optind, argv + optind);

typedef struct {
    samFile *fp;
    hts_itr_t *iter;
    bam_hdr_t *h;
    mplp_ref_t *ref;
    const mplp_conf_t *conf;
} mplp_aux_t;

static int mplp_func(void *data, bam1_t *b)
{
    char *ref;
    mplp_aux_t *ma = (mplp_aux_t*)data;
    int ret, skip = 0, ref_len;
    do {
        int has_ref;
        ret = ma->iter? sam_itr_next(ma->fp, ma->iter, b) : sam_read1(ma->fp, ma->h, b);
        if (ret < 0) break;
        // The 'B' cigar operation is not part of the specification, considering as obsolete.
        //  bam_remove_B(b);
        if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) { // exclude unmapped reads
            skip = 1;
            continue;
        }
        if (ma->conf->rflag_require && !(ma->conf->rflag_require&b->core.flag)) { skip = 1; continue; }
        if (ma->conf->rflag_filter && ma->conf->rflag_filter&b->core.flag) { skip = 1; continue; }
        if (ma->conf->bed && ma->conf->all == 0) { // test overlap
            skip = !bed_overlap(ma->conf->bed, ma->h->target_name[b->core.tid], b->core.pos, bam_endpos(b));
            if (skip) continue;
        }
        if (ma->conf->rghash) { // exclude read groups
            uint8_t *rg = bam_aux_get(b, "RG");
            skip = (rg && khash_str2int_get(ma->conf->rghash, (const char*)(rg+1), NULL)==0);
            if (skip) continue;
        }
        if (ma->conf->flag & MPLP_ILLUMINA13) {
            int i;
            uint8_t *qual = bam_get_qual(b);
            for (i = 0; i < b->core.l_qseq; ++i)
                qual[i] = qual[i] > 31? qual[i] - 31 : 0;
        }

        if (ma->conf->fai && b->core.tid >= 0) {
            has_ref = mplp_get_ref(ma, b->core.tid, &ref, &ref_len);
            if (has_ref && ref_len <= b->core.pos) { // exclude reads outside of the reference sequence
                fprintf(stderr,"[%s] Skipping because %d is outside of %d [ref:%d]\n",
                        __func__, b->core.pos, ref_len, b->core.tid);
                skip = 1;
                continue;
            }
        } else {
            has_ref = 0;
        }

        skip = 0;
        if (has_ref && (ma->conf->flag&MPLP_REALN)) sam_prob_realn(b, ref, ref_len, (ma->conf->flag & MPLP_REDO_BAQ)? 7 : 3);
        if (has_ref && ma->conf->capQ_thres > 10) {
            int q = sam_cap_mapq(b, ref, ref_len, ma->conf->capQ_thres);
            if (q < 0) skip = 1;
            else if (b->core.qual > q) b->core.qual = q;
        }
        if (b->core.qual < ma->conf->min_mq) skip = 1;
        else if ((ma->conf->flag&MPLP_NO_ORPHAN) && (b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) skip = 1;
    } while (skip);
    return ret;
}

static int mpileup(mplp_conf_t *conf, int n, char **fn)
{
    extern void *bcf_call_add_rg(void *rghash, const char *hdtext, const char *list);
    extern void bcf_call_del_rghash(void *rghash);
    mplp_aux_t **data;
    int i, tid, pos, *n_plp, beg0 = 0, end0 = INT_MAX, tid0 = 0, ref_len, max_depth, max_indel_depth;
    const bam_pileup1_t **plp;
    mplp_ref_t mp_ref = MPLP_REF_INIT;
    bam_mplp_t iter;
    bam_hdr_t *h = NULL; /* header of first file in input list */
    char *ref;
    void *rghash = NULL;
    FILE *pileup_fp = NULL;

    bcf_callaux_t *bca = NULL;
    bcf_callret1_t *bcr = NULL;
    bcf_call_t bc;
    htsFile *bcf_fp = NULL;
    bcf_hdr_t *bcf_hdr = NULL;

    bam_sample_t *sm = NULL;
    kstring_t buf;
    mplp_pileup_t gplp;

    memset(&gplp, 0, sizeof(mplp_pileup_t));
    memset(&buf, 0, sizeof(kstring_t));
    memset(&bc, 0, sizeof(bcf_call_t));
    data = calloc(n, sizeof(mplp_aux_t*));
    plp = calloc(n, sizeof(bam_pileup1_t*));
    n_plp = calloc(n, sizeof(int));
    sm = bam_smpl_init();

    if (n == 0) {
        fprintf(stderr,"[%s] no input file/data given\n", __func__);
        exit(EXIT_FAILURE);
    }

    // read the header of each file in the list and initialize data
    for (i = 0; i < n; ++i) {
        bam_hdr_t *h_tmp;
        data[i] = calloc(1, sizeof(mplp_aux_t));
        data[i]->fp = sam_open_format(fn[i], "rb", &conf->ga.in);
        if ( !data[i]->fp )
        {
            fprintf(stderr, "[%s] failed to open %s: %s\n", __func__, fn[i], strerror(errno));
            exit(EXIT_FAILURE);
        }
        if (hts_set_opt(data[i]->fp, CRAM_OPT_DECODE_MD, 0)) {
            fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
            exit(EXIT_FAILURE);
        }
        if (conf->fai_fname && hts_set_fai_filename(data[i]->fp, conf->fai_fname) != 0) {
            fprintf(stderr, "[%s] failed to process %s: %s\n",
                    __func__, conf->fai_fname, strerror(errno));
            exit(EXIT_FAILURE);
        }
        data[i]->conf = conf;
        data[i]->ref = &mp_ref;
        h_tmp = sam_hdr_read(data[i]->fp);
        if ( !h_tmp ) {
            fprintf(stderr,"[%s] fail to read the header of %s\n", __func__, fn[i]);
            exit(EXIT_FAILURE);
        }
        bam_smpl_add(sm, fn[i], (conf->flag&MPLP_IGNORE_RG)? 0 : h_tmp->text);
        // Collect read group IDs with PL (platform) listed in pl_list (note: fragile, strstr search)
        rghash = bcf_call_add_rg(rghash, h_tmp->text, conf->pl_list);
        if (conf->reg) {
            hts_idx_t *idx = sam_index_load(data[i]->fp, fn[i]);
            if (idx == NULL) {
                fprintf(stderr, "[%s] fail to load index for %s\n", __func__, fn[i]);
                exit(EXIT_FAILURE);
            }
            if ( (data[i]->iter=sam_itr_querys(idx, h_tmp, conf->reg)) == 0) {
                fprintf(stderr, "[E::%s] fail to parse region '%s' with %s\n", __func__, conf->reg, fn[i]);
                exit(EXIT_FAILURE);
            }
            if (i == 0) beg0 = data[i]->iter->beg, end0 = data[i]->iter->end, tid0 = data[i]->iter->tid;
            hts_idx_destroy(idx);
        }
        else
            data[i]->iter = NULL;

        if (i == 0) h = data[i]->h = h_tmp; // save the header of the first file
        else {
            // FIXME: check consistency between h and h_tmp
            bam_hdr_destroy(h_tmp);

            // we store only the first file's header; it's (alleged to be)
            // compatible with the i-th file's target_name lookup needs
            data[i]->h = h;
        }
    }
    // allocate data storage proportionate to number of samples being studied sm->n
    gplp.n = sm->n;
    gplp.n_plp = calloc(sm->n, sizeof(int));
    gplp.m_plp = calloc(sm->n, sizeof(int));
    gplp.plp = calloc(sm->n, sizeof(bam_pileup1_t*));

    fprintf(stderr, "[%s] %d samples in %d input files\n", __func__, sm->n, n);
    //if (BCF output) { ... }
    else {
        pileup_fp = conf->output_fname? fopen(conf->output_fname, "w") : stdout;

        if (pileup_fp == NULL) {
            fprintf(stderr, "[%s] failed to write to %s: %s\n", __func__, conf->output_fname, strerror(errno));
            exit(EXIT_FAILURE);
        }
    }

    // init pileup
    iter = bam_mplp_init(n, mplp_func, (void**)data); // looks like mplp_func is called on each read (alignment) to create an iterable over the reads
    if ( conf->flag & MPLP_SMART_OVERLAPS ) bam_mplp_init_overlaps(iter);
    max_depth = conf->max_depth;
    if (max_depth * sm->n > 1<<20)
        fprintf(stderr, "(%s) Max depth is above 1M. Potential memory hog!\n", __func__);
    if (max_depth * sm->n < 8000) {
        max_depth = 8000 / sm->n;
        fprintf(stderr, "<%s> Set max per-file depth to %d\n", __func__, max_depth);
    }
    max_indel_depth = conf->max_indel_depth * sm->n;
    bam_mplp_set_maxcnt(iter, max_depth);
    bcf1_t *bcf_rec = bcf_init1();
    int ret;
    int last_tid = -1, last_pos = -1;

    // begin pileup
    while ( (ret=bam_mplp_auto(iter, &tid, &pos, n_plp, plp)) > 0) {
        if (conf->reg && (pos < beg0 || pos >= end0)) continue; // out of the region requested
        mplp_get_ref(data[0], tid, &ref, &ref_len);
        //printf("tid=%d len=%d ref=%p/%s\n", tid, ref_len, ref, ref);
        if (conf->flag & MPLP_BCF) {
            int total_depth, _ref0, ref16;
            if (conf->bed && tid >= 0 && !bed_overlap(conf->bed, h->target_name[tid], pos, pos+1)) continue;
            for (i = total_depth = 0; i < n; ++i) total_depth += n_plp[i];
            group_smpl(&gplp, sm, &buf, n, fn, n_plp, plp, conf->flag & MPLP_IGNORE_RG);
            _ref0 = (ref && pos < ref_len)? ref[pos] : 'N';
            ref16 = seq_nt16_table[_ref0];
            bcf_callaux_clean(bca, &bc);
            for (i = 0; i < gplp.n; ++i)
                bcf_call_glfgen(gplp.n_plp[i], gplp.plp[i], ref16, bca, bcr + i);
            bc.tid = tid; bc.pos = pos;
            bcf_call_combine(gplp.n, bcr, bca, ref16, &bc);
            bcf_clear1(bcf_rec);
            bcf_call2bcf(&bc, bcf_rec, bcr, conf->fmt_flag, 0, 0);
            bcf_write1(bcf_fp, bcf_hdr, bcf_rec);
            // call indels; todo: subsampling with total_depth>max_indel_depth instead of ignoring?
            if (!(conf->flag&MPLP_NO_INDEL) && total_depth < max_indel_depth && bcf_call_gap_prep(gplp.n, gplp.n_plp, gplp.plp, pos, bca, ref, rghash) >= 0)
            {
                bcf_callaux_clean(bca, &bc);
                for (i = 0; i < gplp.n; ++i)
                    bcf_call_glfgen(gplp.n_plp[i], gplp.plp[i], -1, bca, bcr + i);
                if (bcf_call_combine(gplp.n, bcr, bca, -1, &bc) >= 0) {
                    bcf_clear1(bcf_rec);
                    bcf_call2bcf(&bc, bcf_rec, bcr, conf->fmt_flag, bca, ref);
                    bcf_write1(bcf_fp, bcf_hdr, bcf_rec);
                }
            }
        } else {
            fprintf(pileup_fp, "%s\t%d\t%c", h->target_name[tid], pos + 1, (ref && pos < ref_len)? ref[pos] : 'N'); //first, second, third column (chr, pos, ref nucleotide)
            for (i = 0; i < n; ++i) {
                int j, cnt;
                for (j = cnt = 0; j < n_plp[i]; ++j) {
                    const bam_pileup1_t *p = plp[i] + j;
                    int c = p->qpos < p->b->core.l_qseq
                             ? bam_get_qual(p->b)[p->qpos]
                             : 0;
                    if (c >= conf->min_baseQ) ++cnt;
                }
                fprintf(pileup_fp, "\t%d\t", cnt); //depth, 4th column

                //handle bases
                if (n_plp[i] == 0) {
                    fputs("*\t*", pileup_fp);
                    if (conf->flag & MPLP_PRINT_MAPQ) fputs("\t*", pileup_fp);
                    if (conf->flag & MPLP_PRINT_POS) fputs("\t*", pileup_fp);
                } else {
                    int n = 0;
                    for (j = 0; j < n_plp[i]; ++j) {
                        const bam_pileup1_t *p = plp[i] + j;
                        int c = p->qpos < p->b->core.l_qseq
                            ? bam_get_qual(p->b)[p->qpos]
                            : 0;
                        if (c >= conf->min_baseQ)
                            n++, pileup_seq(pileup_fp, plp[i] + j, pos, ref_len, ref);
                    }
                    if (!n) putc('*', pileup_fp);

                    n = 0;
                    putc('\t', pileup_fp);

                    //handle base quality
                    for (j = 0; j < n_plp[i]; ++j) {
                        const bam_pileup1_t *p = plp[i] + j;
                        int c = p->qpos < p->b->core.l_qseq
                            ? bam_get_qual(p->b)[p->qpos]
                            : 0;
                        if (c >= conf->min_baseQ) {
                            c = c + 33 < 126? c + 33 : 126;
                            putc(c, pileup_fp);
                            n++;
                        }
                    }
                    if (!n) putc('*', pileup_fp);

                    //optionally handle alignment quality
                    if (conf->flag & MPLP_PRINT_MAPQ) {
                        n = 0;
                        putc('\t', pileup_fp);
                        for (j = 0; j < n_plp[i]; ++j) {
                            const bam_pileup1_t *p = plp[i] + j;
                            int c = bam_get_qual(p->b)[p->qpos];
                            if ( c < conf->min_baseQ ) continue;
                            c = plp[i][j].b->core.qual + 33;
                            if (c > 126) c = 126;
                            putc(c, pileup_fp);
                            n++;
                        }
                        if (!n) putc('*', pileup_fp);
                    }

                    //optionally print base positions on reads
                    if (conf->flag & MPLP_PRINT_POS) {
                        n = 0;
                        putc('\t', pileup_fp);
                        for (j = 0; j < n_plp[i]; ++j) {
                            const bam_pileup1_t *p = plp[i] + j;
                            int c = bam_get_qual(p->b)[p->qpos];
                            if ( c < conf->min_baseQ ) continue;

                            if (n > 0) putc(',', pileup_fp);
                            fprintf(pileup_fp, "%d", plp[i][j].qpos + 1); // FIXME: printf() is very slow...
                            n++;
                        }
                        if (!n) putc('*', pileup_fp);
                    }
                }
            } //end of iteration through positions
            putc('\n', pileup_fp); // eol
        } //end of pileup format output
    } //end of while loop

    if (conf->all && !(conf->flag & MPLP_BCF)) {
        // Handle terminating region
        if (last_tid < 0 && conf->reg && conf->all > 1) {
            last_tid = tid0;
            last_pos = beg0-1;
            mplp_get_ref(data[0], tid0, &ref, &ref_len);
        }
       while (last_tid >= 0 && last_tid < h->n_targets) {
            while (++last_pos < h->target_len[last_tid]) {
                if (last_pos >= end0) break;
                if (conf->bed && bed_overlap(conf->bed, h->target_name[last_tid], last_pos, last_pos + 1) == 0)
                    continue;
                print_empty_pileup(pileup_fp, conf, h->target_name[last_tid], last_pos, n, ref, ref_len);
            }
            last_tid++;
            last_pos = -1;
            if (conf->all < 2 || conf->reg)
                break;
        }
    }

    // clean up
    free(bc.tmp.s);
    bcf_destroy1(bcf_rec);
    if (bcf_fp)
    {
        hts_close(bcf_fp);
        bcf_hdr_destroy(bcf_hdr);
        bcf_call_destroy(bca);
        free(bc.PL);
        free(bc.DP4);
        free(bc.ADR);
        free(bc.ADF);
        free(bc.fmt_arr);
        free(bcr);
    }
    if (pileup_fp && conf->output_fname) fclose(pileup_fp);
    bam_smpl_destroy(sm); free(buf.s);
    for (i = 0; i < gplp.n; ++i) free(gplp.plp[i]);
    free(gplp.plp); free(gplp.n_plp); free(gplp.m_plp);
    bcf_call_del_rghash(rghash);
    bam_mplp_destroy(iter);
    bam_hdr_destroy(h);
    for (i = 0; i < n; ++i) {
        sam_close(data[i]->fp);
        if (data[i]->iter) hts_itr_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data); free(plp); free(n_plp);
    free(mp_ref.ref[0]);
    free(mp_ref.ref[1]);
    return ret;
}