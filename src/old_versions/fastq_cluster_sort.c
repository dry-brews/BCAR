#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <getopt.h>
#include <time.h>
#include <sys/stat.h>
#include <errno.h>
#include <limits.h>

#include "sort_module.h"

/* ------------------------------------------------------------------ */
/*  Program parameters                                                */
/* ------------------------------------------------------------------ */

typedef struct {
    char **in_files;
    int    n_in;
    char  *out_file;
    char  *fail_file;
    char  *tmp_dir;
    long   chunk_mem;
    int    max_open_files;
    int    max_mismatches;
    int    max_bc_indels;
    int    max_n;
    extraction_config_t extract;
} params_t;

/* ------------------------------------------------------------------ */
/*  Globals                                                           */
/* ------------------------------------------------------------------ */

/* cleanup state */
static char g_tmp_dir[PATH_MAX] = {0};
static char **g_chunk_paths = NULL;
static int    g_n_chunks = 0;

/* ------------------------------------------------------------------ */
/*  Timing utilities                                                  */
/* ------------------------------------------------------------------ */

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void format_elapsed(double seconds, char *buf, int buflen) {
    int h = (int)(seconds / 3600);
    int m = (int)((seconds - h * 3600) / 60);
    int s = (int)(seconds - h * 3600 - m * 60);
    if (h > 0)
        snprintf(buf, buflen, "%dh %dm %ds", h, m, s);
    else if (m > 0)
        snprintf(buf, buflen, "%dm %ds", m, s);
    else
        snprintf(buf, buflen, "%ds", s);
}

/* ------------------------------------------------------------------ */
/*  Context pattern parsing                                           */
/* ------------------------------------------------------------------ */

/* Parse --context string like "^ATGCGTNNNNNN$" into context_pattern_t.
   Returns 0 on success, -1 on error. */
static int parse_context(const char *pattern, context_pattern_t *out) {
    memset(out, 0, sizeof(*out));
    const char *p = pattern;
    int slen = (int)strlen(pattern);

    /* Check for anchors */
    if (p[0] == '^') {
        out->anchored_5 = 1;
        p++;
        slen--;
    }
    /* Work on a mutable copy for trailing $ */
    char *work = xstrdup(p, "parse_context");
    int wlen = (int)strlen(work);
    if (wlen > 0 && work[wlen - 1] == '$') {
        out->anchored_3 = 1;
        work[wlen - 1] = '\0';
        wlen--;
    }

    /* Find first and last N */
    int first_n = -1, last_n = -1;
    for (int i = 0; i < wlen; i++) {
        if (work[i] == 'N' || work[i] == 'n') {
            if (first_n < 0) first_n = i;
            last_n = i;
        }
    }
    if (first_n < 0) {
        fprintf(stderr, "Error: --context pattern must contain at least one N\n");
        free(work);
        return -1;
    }

    /* Validate N's are contiguous */
    for (int i = first_n; i <= last_n; i++) {
        if (work[i] != 'N' && work[i] != 'n') {
            fprintf(stderr, "Error: N's in --context pattern must be contiguous\n");
            free(work);
            return -1;
        }
    }

    out->bc_len = last_n - first_n + 1;
    if (out->bc_len > MAX_BC_LEN) {
        fprintf(stderr, "Error: barcode length %d exceeds max %d\n", out->bc_len, MAX_BC_LEN);
        free(work);
        return -1;
    }

    /* Extract prefix and suffix */
    out->prefix_len = first_n;
    if (out->prefix_len > 0) memcpy(out->prefix, work, out->prefix_len);
    out->prefix[out->prefix_len] = '\0';

    out->suffix_len = wlen - last_n - 1;
    if (out->suffix_len > 0) memcpy(out->suffix, work + last_n + 1, out->suffix_len);
    out->suffix[out->suffix_len] = '\0';

    free(work);
    return 0;
}

/* ------------------------------------------------------------------ */
/*  CLI parsing                                                       */
/* ------------------------------------------------------------------ */

static void print_usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options]\n\n"
        "Required:\n"
        "  --in <file[,file,...]>           FASTQ input(s), comma-separated, gzip OK\n\n"
        "Barcode extraction (choose mode):\n"
        "  --bc-start <int>                 Fixed position start (0-based, default: 0)\n"
        "  --bc-len <int>                   Fixed barcode length (default: 18)\n"
        "  --context <pattern>              Context pattern with N's for barcode positions\n"
        "                                     e.g. ATGCGTNNNNNNTTGACT or ^NNNNNNNNNNNNNNNNNN$\n"
        "  --max-context-mismatches <int>   Allowed mismatches in flanking seqs (default: 1)\n\n"
        "Clustering:\n"
        "  --max-mismatches <int>           Hamming distance radius (default: 1)\n"
        "  --max-bc-indels <int>            Max indels in clustering/extraction (0-2, default: 0)\n\n"
        "N handling:\n"
        "  --max-n <int>                    Max N bases in barcode (default: 2)\n\n"
        "Output:\n"
        "  --out <file>                     Sorted output (default: sorted_out.fastq)\n"
        "  --fail <file>                    Reads where barcode not found\n\n"
        "Resources:\n"
        "  --chunk-mem <bytes>              Sort buffer size (default: 4294967296 = 4GB)\n"
        "  --max-open-files <int>           Max files during merge (default: 240)\n"
        "  --tmp-dir <path>                 Temp directory (default: .)\n\n",
        prog);
    exit(EXIT_FAILURE);
}

static params_t parse_args(int argc, char *argv[]) {
    params_t p;
    memset(&p, 0, sizeof(p));
    p.out_file = "sorted_out.fastq";
    p.fail_file = NULL;
    p.tmp_dir = ".";
    p.chunk_mem = 4LL * 1024 * 1024 * 1024;
    p.max_open_files = 240;
    p.max_mismatches = 1;
    p.max_bc_indels = 0;
    p.max_n = 2;
    p.extract.mode = 0;
    p.extract.bc_start = 0;
    p.extract.bc_len = 18;
    p.extract.max_context_mm = 1;
    p.extract.max_bc_indels = 0;

    static struct option long_options[] = {
        {"in",                     required_argument, 0, 'a'},
        {"bc-start",               required_argument, 0, 's'},
        {"bc-len",                 required_argument, 0, 'l'},
        {"context",                required_argument, 0, 'c'},
        {"max-context-mismatches", required_argument, 0, 'x'},
        {"max-mismatches",         required_argument, 0, 'm'},
        {"max-bc-indels",          required_argument, 0, 'i'},
        {"max-n",                  required_argument, 0, 'n'},
        {"out",                    required_argument, 0, 'o'},
        {"fail",                   required_argument, 0, 'f'},
        {"chunk-mem",              required_argument, 0, 'k'},
        {"max-open-files",         required_argument, 0, 'u'},
        {"tmp-dir",                required_argument, 0, 't'},
        {0, 0, 0, 0}
    };

    int c, option_index = 0;
    char *in_arg = NULL;
    char *context_arg = NULL;
    while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
        switch (c) {
            case 'a': in_arg = optarg; break;
            case 's': p.extract.bc_start = atoi(optarg); break;
            case 'l': p.extract.bc_len = atoi(optarg); break;
            case 'c': context_arg = optarg; break;
            case 'x': p.extract.max_context_mm = atoi(optarg); break;
            case 'm': p.max_mismatches = atoi(optarg); break;
            case 'i': p.max_bc_indels = atoi(optarg); break;
            case 'n': p.max_n = atoi(optarg); break;
            case 'o': p.out_file = optarg; break;
            case 'f': p.fail_file = optarg; break;
            case 'k': p.chunk_mem = atol(optarg); break;
            case 'u': p.max_open_files = atoi(optarg); break;
            case 't': p.tmp_dir = optarg; break;
            default: print_usage(argv[0]);
        }
    }

    if (!in_arg) {
        fprintf(stderr, "Error: --in is required\n");
        print_usage(argv[0]);
    }

    /* Split comma-separated input files */
    char *in_copy = xstrdup(in_arg, "parse_args in");
    int n_files = 0;
    char *tok = strtok(in_copy, ",");
    while (tok) { n_files++; tok = strtok(NULL, ","); }
    free(in_copy);

    p.in_files = xmalloc(n_files * sizeof(char *), "parse_args in_files");
    p.n_in = n_files;
    in_copy = xstrdup(in_arg, "parse_args in2");
    tok = strtok(in_copy, ",");
    for (int i = 0; i < n_files; i++) {
        p.in_files[i] = xstrdup(tok, "parse_args file");
        tok = strtok(NULL, ",");
    }
    free(in_copy);

    /* Context mode */
    if (context_arg) {
        p.extract.mode = 1;
        if (parse_context(context_arg, &p.extract.ctx) < 0)
            exit(EXIT_FAILURE);
        p.extract.bc_len = p.extract.ctx.bc_len;
    }

    /* Propagate indel setting */
    p.extract.max_bc_indels = p.max_bc_indels;

    if (p.max_bc_indels > 2) {
        fprintf(stderr, "Error: --max-bc-indels must be 0, 1, or 2\n");
        exit(EXIT_FAILURE);
    }

    return p;
}

/* ------------------------------------------------------------------ */
/*  Main                                                              */
/* ------------------------------------------------------------------ */

int main(int argc, char *argv[]) {
    params_t params = parse_args(argc, argv);
    double t_start = now_sec();

    /* Print run info */
    fprintf(stderr, "\n=== BCAR Cluster Sort ===\n");
    fprintf(stderr, "Input files: %d\n", params.n_in);
    for (int i = 0; i < params.n_in; i++)
        fprintf(stderr, "  %s\n", params.in_files[i]);
    fprintf(stderr, "Output: %s\n", params.out_file);
    if (params.fail_file) fprintf(stderr, "Fail file: %s\n", params.fail_file);
    if (params.extract.mode == 0) {
        fprintf(stderr, "Extraction: fixed position (start=%d, len=%d)\n",
                params.extract.bc_start, params.extract.bc_len);
    } else {
        fprintf(stderr, "Extraction: context pattern (bc_len=%d, prefix='%s', suffix='%s'",
                params.extract.ctx.bc_len, params.extract.ctx.prefix, params.extract.ctx.suffix);
        if (params.extract.ctx.anchored_5) fprintf(stderr, ", anchored 5'");
        if (params.extract.ctx.anchored_3) fprintf(stderr, ", anchored 3'");
        fprintf(stderr, ")\n");
    }
    fprintf(stderr, "Clustering: max_mismatches=%d, max_bc_indels=%d, max_n=%d\n",
            params.max_mismatches, params.max_bc_indels, params.max_n);
    fprintf(stderr, "Chunk memory: %ld bytes\n", params.chunk_mem);
    fprintf(stderr, "\n");

    /* Create temp directory */
    char tmp_template[PATH_MAX];
    snprintf(tmp_template, PATH_MAX, "%s/bcar_sort_XXXXXX", params.tmp_dir);
    char *tmp_result = mkdtemp(tmp_template);
    if (!tmp_result) {
        fprintf(stderr, "Fatal: cannot create temp directory in %s: %s\n",
                params.tmp_dir, strerror(errno));
        exit(EXIT_FAILURE);
    }
    strncpy(g_tmp_dir, tmp_result, PATH_MAX - 1);
    fprintf(stderr, "Temp directory: %s\n\n", g_tmp_dir);

    /* ============================================================== */
    /*  Phase A: Count unique barcodes                                */
    /* ============================================================== */
    fprintf(stderr, "Phase A: Counting unique barcodes...\n");
    double t_phase = now_sec();

    barcode_ht_t ht;
    ht_init(&ht, HT_INITIAL_CAP);

    long total_reads = 0;
    long n_with_n = 0;
    long n_failed = 0;
    long reads_per_file[params.n_in];
    memset(reads_per_file, 0, sizeof(reads_per_file));

    char hdr[MAX_LINE], seq[MAX_LINE], plus[MAX_LINE], qual[MAX_LINE];

    for (int fi = 0; fi < params.n_in; fi++) {
        gz_reader_t *reader = gz_reader_open(params.in_files[fi]);
        while (read_fastq(reader, hdr, seq, plus, qual)) {
            total_reads++;
            reads_per_file[fi]++;

            extraction_result_t ext;
            int slen = (int)strlen(seq);
            if (extract_barcode(seq, slen, &params.extract, &ext) < 0) {
                n_failed++;
                continue;
            }

            barcode_t bc;
            int nn = encode_barcode(ext.barcode, ext.bc_len, &bc);
            if (nn < 0) {
                n_failed++;
                continue;
            }
            if (nn > params.max_n) {
                n_failed++;
                continue;
            }
            if (nn > 0) {
                n_with_n++;
                continue; /* N-containing: counted but not inserted */
            }

            bc_key_t key;
            barcode_to_key(&bc, &key);
            ht_insert(&ht, &key);

            if (total_reads % 10000000 == 0)
                fprintf(stderr, "  %ld reads scanned...\n", total_reads);
        }
        gz_reader_close(reader);
    }

    char elapsed[64];
    format_elapsed(now_sec() - t_phase, elapsed, sizeof(elapsed));
    fprintf(stderr, "\nPhase A complete (%s)\n", elapsed);
    fprintf(stderr, "  Total reads: %ld\n", total_reads);
    for (int i = 0; i < params.n_in; i++)
        fprintf(stderr, "  %s: %ld reads\n", params.in_files[i], reads_per_file[i]);
    fprintf(stderr, "  Unique barcodes: %lu\n", (unsigned long)ht.count);
    fprintf(stderr, "  Reads with N's in barcode: %ld\n", n_with_n);
    fprintf(stderr, "  Failed extraction: %ld\n\n", n_failed);

    /* ============================================================== */
    /*  Phase B: Cluster barcodes                                     */
    /* ============================================================== */
    fprintf(stderr, "Phase B: Clustering barcodes...\n");
    t_phase = now_sec();

    /* Assign uf_idx to each occupied slot */
    union_find_t uf;
    uf_init(&uf, (uint32_t)ht.count);
    {
        uint32_t idx = 0;
        for (uint64_t i = 0; i < ht.capacity; i++) {
            if (ht.slots[i].occupied) {
                ht.slots[i].uf_idx = idx++;
            }
        }
    }

    cluster_barcodes(&ht, &uf, params.max_mismatches, params.max_bc_indels);

    uint64_t *ubid_map = NULL;
    uint64_t n_clusters = assign_ubids(&ht, &uf, &ubid_map);

    format_elapsed(now_sec() - t_phase, elapsed, sizeof(elapsed));
    fprintf(stderr, "\nPhase B complete (%s)\n", elapsed);
    fprintf(stderr, "  Clusters: %lu\n\n", (unsigned long)n_clusters);

    uf_free(&uf);

    /* ============================================================== */
    /*  Phase C: Re-read, assign UBIDs, sort chunks                   */
    /* ============================================================== */
    fprintf(stderr, "Phase C: Assigning UBIDs and creating sorted chunks...\n");
    t_phase = now_sec();

    FILE *fail_fp = NULL;
    if (params.fail_file) {
        fail_fp = fopen(params.fail_file, "w");
        if (!fail_fp) {
            fprintf(stderr, "Fatal: cannot open fail file: %s\n", params.fail_file);
            exit(EXIT_FAILURE);
        }
    }

    /* Chunk file tracking */
    int chunk_cap = 256;
    char **chunk_paths = xmalloc(chunk_cap * sizeof(char *), "chunk paths");
    int n_chunks = 0;

    chunk_buffer_t cbuf;
    chunk_init(&cbuf, 100000);

    long assigned_reads = 0;
    long fail_reads = 0;

    for (int fi = 0; fi < params.n_in; fi++) {
        gz_reader_t *reader = gz_reader_open(params.in_files[fi]);
        while (read_fastq(reader, hdr, seq, plus, qual)) {
            int slen = (int)strlen(seq);
            extraction_result_t ext;
            uint64_t ubid = 0;

            if (extract_barcode(seq, slen, &params.extract, &ext) < 0) {
                /* Failed extraction */
                if (fail_fp)
                    fprintf(fail_fp, "%s\n%s\n%s\n%s\n", hdr, seq, plus, qual);
                fail_reads++;
                continue;
            }

            barcode_t bc;
            int nn = encode_barcode(ext.barcode, ext.bc_len, &bc);
            if (nn < 0 || nn > params.max_n) {
                if (fail_fp)
                    fprintf(fail_fp, "%s\n%s\n%s\n%s\n", hdr, seq, plus, qual);
                fail_reads++;
                continue;
            }

            if (nn == 0) {
                /* Concrete barcode: direct lookup */
                bc_key_t key;
                barcode_to_key(&bc, &key);
                ht_entry_t *e = ht_lookup(&ht, &key);
                if (e) {
                    uint64_t slot_idx = (uint64_t)(e - ht.slots);
                    ubid = ubid_map[slot_idx];
                }
            } else {
                /* N-containing: enumerate possibilities */
                ubid = resolve_n_barcode(&bc, &ht, ubid_map);
            }

            if (ubid == 0) {
                /* No cluster match */
                if (fail_fp)
                    fprintf(fail_fp, "%s\n%s\n%s\n%s\n", hdr, seq, plus, qual);
                fail_reads++;
                continue;
            }

            /* Append UBID to header */
            char new_hdr[MAX_LINE];
            snprintf(new_hdr, MAX_LINE, "%s;UBID=%lu", hdr, (unsigned long)ubid);

            chunk_push(&cbuf, new_hdr, seq, plus, qual, ubid);
            assigned_reads++;

            if (cbuf.mem_used >= params.chunk_mem) {
                if (n_chunks >= chunk_cap) {
                    chunk_cap *= 2;
                    chunk_paths = xrealloc(chunk_paths, chunk_cap * sizeof(char *), "chunk paths grow");
                }
                chunk_paths[n_chunks] = write_chunk(&cbuf, n_chunks, g_tmp_dir);
                n_chunks++;
                chunk_clear(&cbuf);
            }

            if ((assigned_reads + fail_reads) % 10000000 == 0)
                fprintf(stderr, "  %ld reads processed...\n", assigned_reads + fail_reads);
        }
        gz_reader_close(reader);
    }

    /* Flush remaining chunk */
    if (cbuf.n > 0) {
        if (n_chunks >= chunk_cap) {
            chunk_cap *= 2;
            chunk_paths = xrealloc(chunk_paths, chunk_cap * sizeof(char *), "chunk paths grow");
        }
        chunk_paths[n_chunks] = write_chunk(&cbuf, n_chunks, g_tmp_dir);
        n_chunks++;
    }
    chunk_free(&cbuf);

    if (fail_fp) fclose(fail_fp);

    /* Store for cleanup */
    g_chunk_paths = chunk_paths;
    g_n_chunks = n_chunks;

    format_elapsed(now_sec() - t_phase, elapsed, sizeof(elapsed));
    fprintf(stderr, "\nPhase C complete (%s)\n", elapsed);
    fprintf(stderr, "  Assigned reads: %ld\n", assigned_reads);
    fprintf(stderr, "  Failed/unmatched reads: %ld\n", fail_reads);
    fprintf(stderr, "  Chunks created: %d\n\n", n_chunks);

    /* Free hash table and UBID map (no longer needed) */
    ht_free(&ht);
    free(ubid_map);

    /* ============================================================== */
    /*  Phase D: K-way merge                                          */
    /* ============================================================== */
    fprintf(stderr, "Phase D: Merging sorted chunks...\n");
    t_phase = now_sec();

    merge_sorted_chunks(chunk_paths, n_chunks, params.out_file,
                        params.max_open_files, g_tmp_dir);

    format_elapsed(now_sec() - t_phase, elapsed, sizeof(elapsed));
    fprintf(stderr, "\nPhase D complete (%s)\n", elapsed);

    /* Cleanup temp files */
    for (int i = 0; i < n_chunks; i++) {
        remove(chunk_paths[i]);
        free(chunk_paths[i]);
    }
    free(chunk_paths);
    g_chunk_paths = NULL;
    g_n_chunks = 0;
    rmdir(g_tmp_dir);

    /* Free input file list */
    for (int i = 0; i < params.n_in; i++) free(params.in_files[i]);
    free(params.in_files);

    /* Summary */
    format_elapsed(now_sec() - t_start, elapsed, sizeof(elapsed));
    fprintf(stderr, "\n=== Sort complete ===\n");
    fprintf(stderr, "Total time: %s\n", elapsed);
    fprintf(stderr, "Total reads: %ld\n", total_reads);
    fprintf(stderr, "Assigned: %ld (%.1f%%)\n", assigned_reads,
            total_reads > 0 ? 100.0 * assigned_reads / total_reads : 0.0);
    fprintf(stderr, "Clusters: %lu\n", (unsigned long)n_clusters);
    fprintf(stderr, "Output: %s\n\n", params.out_file);

    return 0;
}
