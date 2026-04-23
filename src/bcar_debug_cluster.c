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
#include <unistd.h>

#include "sort_module.h"

/* ------------------------------------------------------------------ */
/*  Timing                                                            */
/* ------------------------------------------------------------------ */

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ------------------------------------------------------------------ */
/*  Program parameters                                                */
/* ------------------------------------------------------------------ */

typedef struct {
    char **in_files;
    int    n_in;
    char  *log_dir;
    int    max_mismatches;
    int    max_bc_indels;
    int    max_n;
    extraction_config_t extract;
} params_t;

/* ------------------------------------------------------------------ */
/*  Context pattern parsing (identical to bcar_pipeline)              */
/* ------------------------------------------------------------------ */

static int parse_context(const char *pattern, context_pattern_t *out) {
    memset(out, 0, sizeof(*out));
    const char *p = pattern;
    int slen = (int)strlen(pattern);

    if (p[0] == '^') {
        out->anchored_5 = 1;
        p++;
        slen--;
    }
    char *work = xstrdup(p, "parse_context");
    int wlen = (int)strlen(work);
    if (wlen > 0 && work[wlen - 1] == '$') {
        out->anchored_3 = 1;
        work[wlen - 1] = '\0';
        wlen--;
    }

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
        "  Headers must contain 'barcode=<seq>' giving the ground truth barcode.\n\n"
        "Barcode extraction (choose mode):\n"
        "  --bc-start <int>                 Fixed position start (0-based, default: 0)\n"
        "  --bc-len <int>                   Fixed barcode length (default: 15)\n"
        "  --context <pattern>              Context pattern with N's for barcode positions\n"
        "  --max-context-mismatches <int>   Allowed mismatches in flanking seqs (default: 1)\n\n"
        "Clustering:\n"
        "  --max-mismatches <int>           Max substitutions between clustered barcodes (default: 2)\n"
        "  --max-bc-indels <int>            Max indels between clustered barcodes (0-2, default: 1)\n"
        "  --max-n <int>                    Max N bases in barcode (default: 2)\n\n"
        "Output:\n"
        "  --log-dir <path>                 Directory for log file (default: .)\n\n",
        prog);
    exit(EXIT_FAILURE);
}

static params_t parse_args(int argc, char *argv[]) {
    params_t p;
    memset(&p, 0, sizeof(p));
    p.log_dir           = ".";
    p.max_mismatches    = 2;
    p.max_bc_indels     = 1;
    p.max_n             = 2;
    p.extract.mode      = 0;
    p.extract.bc_start  = 0;
    p.extract.bc_len    = 15;
    p.extract.max_context_mm  = 1;
    p.extract.max_bc_indels   = 0;

    static struct option long_options[] = {
        {"in",                     required_argument, 0, 'a'},
        {"bc-start",               required_argument, 0, 's'},
        {"bc-len",                 required_argument, 0, 'l'},
        {"context",                required_argument, 0, 'c'},
        {"max-context-mismatches", required_argument, 0, 'x'},
        {"max-mismatches",         required_argument, 0, 'm'},
        {"max-bc-indels",          required_argument, 0, 'i'},
        {"max-n",                  required_argument, 0, 'n'},
        {"log-dir",                required_argument, 0, 'd'},
        {0, 0, 0, 0}
    };

    int c, option_index = 0;
    char *in_arg = NULL;
    char *context_arg = NULL;

    while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
        switch (c) {
            case 'a': in_arg        = optarg;          break;
            case 's': p.extract.bc_start         = atoi(optarg); break;
            case 'l': p.extract.bc_len           = atoi(optarg); break;
            case 'c': context_arg   = optarg;          break;
            case 'x': p.extract.max_context_mm   = atoi(optarg); break;
            case 'm': p.max_mismatches           = atoi(optarg); break;
            case 'i': p.max_bc_indels            = atoi(optarg); break;
            case 'n': p.max_n                    = atoi(optarg); break;
            case 'd': p.log_dir     = optarg;          break;
            default:  print_usage(argv[0]);
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

    if (context_arg) {
        p.extract.mode = 1;
        if (parse_context(context_arg, &p.extract.ctx) < 0)
            exit(EXIT_FAILURE);
        p.extract.bc_len = p.extract.ctx.bc_len;
    }

    p.extract.max_bc_indels = p.max_bc_indels;

    if (p.max_bc_indels > 2) {
        fprintf(stderr, "Error: --max-bc-indels must be 0, 1, or 2\n");
        exit(EXIT_FAILURE);
    }

    return p;
}

/* ------------------------------------------------------------------ */
/*  Read pair: (ground truth barcode, assigned cluster UBID)          */
/* ------------------------------------------------------------------ */

typedef struct {
    char     true_bc[MAX_BC_LEN + 1];
    uint64_t ubid;
} read_pair_t;

/* Composite sorts for the two analyses */
static int cmp_by_true_bc_ubid(const void *a, const void *b) {
    const read_pair_t *A = a, *B = b;
    int c = strcmp(A->true_bc, B->true_bc);
    if (c) return c;
    if (A->ubid < B->ubid) return -1;
    if (A->ubid > B->ubid) return  1;
    return 0;
}

static int cmp_by_ubid_true_bc(const void *a, const void *b) {
    const read_pair_t *A = a, *B = b;
    if (A->ubid < B->ubid) return -1;
    if (A->ubid > B->ubid) return  1;
    return strcmp(A->true_bc, B->true_bc);
}

/* Parse "barcode=<seq>" from a FASTQ header line.
   Writes the sequence into out (must be MAX_BC_LEN+1 bytes).
   Returns 0 on success, -1 if the tag is absent or empty. */
static int parse_true_bc(const char *hdr, char *out) {
    const char *p = strstr(hdr, "barcode=");
    if (!p) return -1;
    p += 8;
    int i = 0;
    while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r'
           && i < MAX_BC_LEN) {
        out[i++] = *p++;
    }
    out[i] = '\0';
    return i > 0 ? 0 : -1;
}

/* ------------------------------------------------------------------ */
/*  Main                                                              */
/* ------------------------------------------------------------------ */

int main(int argc, char *argv[]) {
    params_t params = parse_args(argc, argv);
    double t_start = now_sec();

    /* Build uniquely named log file */
    char log_path[PATH_MAX];
    {
        time_t now = time(NULL);
        struct tm *tm_info = localtime(&now);
        snprintf(log_path, PATH_MAX,
                 "%s/bcar_debug_%04d%02d%02d_%02d%02d%02d_%d.log",
                 params.log_dir,
                 tm_info->tm_year + 1900, tm_info->tm_mon + 1, tm_info->tm_mday,
                 tm_info->tm_hour, tm_info->tm_min, tm_info->tm_sec,
                 (int)getpid());
    }

    FILE *log = fopen(log_path, "w");
    if (!log) {
        fprintf(stderr, "Fatal: cannot open log file %s: %s\n",
                log_path, strerror(errno));
        exit(EXIT_FAILURE);
    }

    /* Log parameters */
    fprintf(log, "=== BCAR Debug Cluster ===\n");
    fprintf(log, "Input files: %d\n", params.n_in);
    for (int i = 0; i < params.n_in; i++)
        fprintf(log, "  %s\n", params.in_files[i]);
    if (params.extract.mode == 0) {
        fprintf(log, "Extraction:               fixed position (start=%d, len=%d)\n",
                params.extract.bc_start, params.extract.bc_len);
    } else {
        fprintf(log, "Extraction:               context pattern (bc_len=%d, prefix='%s', suffix='%s'",
                params.extract.ctx.bc_len,
                params.extract.ctx.prefix,
                params.extract.ctx.suffix);
        if (params.extract.ctx.anchored_5) fprintf(log, ", anchored 5'");
        if (params.extract.ctx.anchored_3) fprintf(log, ", anchored 3'");
        fprintf(log, ")\n");
        fprintf(log, "Max context mismatches:   %d\n", params.extract.max_context_mm);
    }
    fprintf(log, "max_mismatches:           %d\n", params.max_mismatches);
    fprintf(log, "max_bc_indels:            %d\n", params.max_bc_indels);
    fprintf(log, "max_n:                    %d\n\n", params.max_n);

    /* ============================================================== */
    /*  Phase A: Count unique barcodes                                */
    /* ============================================================== */
    fprintf(stderr, "Phase A: Counting unique barcodes...\n");
    double t_phase = now_sec();

    barcode_ht_t ht;
    ht_init(&ht, HT_INITIAL_CAP);

    long total_reads = 0;
    long n_with_n    = 0;
    long n_failed    = 0;

    char hdr[MAX_LINE], seq[MAX_LINE], plus[MAX_LINE], qual[MAX_LINE];

    for (int fi = 0; fi < params.n_in; fi++) {
        gz_reader_t *reader = gz_reader_open(params.in_files[fi]);
        while (read_fastq(reader, hdr, seq, plus, qual)) {
            total_reads++;

            extraction_result_t ext;
            int slen = (int)strlen(seq);
            if (extract_barcode(seq, slen, &params.extract, &ext) < 0) {
                n_failed++;
                continue;
            }

            barcode_t bc;
            int nn = encode_barcode(ext.barcode, ext.bc_len, &bc);
            if (nn < 0)              { n_failed++; continue; }
            if (nn > params.max_n)   { n_failed++; continue; }
            if (nn > 0)              { n_with_n++;  continue; }

            bc_key_t key;
            barcode_to_key(&bc, &key);
            ht_insert(&ht, &key);

            if (total_reads % 10000000 == 0)
                fprintf(stderr, "  %ld reads scanned...\n", total_reads);
        }
        gz_reader_close(reader);
    }

    double t_a = now_sec() - t_phase;
    fprintf(log, "=== Phase A: Count unique barcodes (%.1fs) ===\n", t_a);
    fprintf(log, "Total reads:              %ld\n",              total_reads);
    fprintf(log, "Unique barcodes:          %lu\n",  (unsigned long)ht.count);
    fprintf(log, "Reads with N in barcode:  %ld\n",              n_with_n);
    fprintf(log, "Failed extraction:        %ld\n\n",            n_failed);
    fprintf(stderr, "Phase A complete (%.1fs): %lu unique barcodes\n\n",
            t_a, (unsigned long)ht.count);

    /* ============================================================== */
    /*  Phase B: Cluster                                              */
    /* ============================================================== */
    fprintf(stderr, "Phase B: Clustering barcodes...\n");
    t_phase = now_sec();

    union_find_t uf;
    uf_init(&uf, (uint32_t)ht.count);
    {
        uint32_t idx = 0;
        for (uint64_t i = 0; i < ht.capacity; i++) {
            if (ht.slots[i].occupied)
                ht.slots[i].uf_idx = idx++;
        }
    }

    cluster_barcodes(&ht, &uf, params.max_mismatches, &params.extract);

    uint64_t *ubid_map  = NULL;
    uint64_t  n_clusters = assign_ubids(&ht, &uf, &ubid_map);

    uf_free(&uf);

    double t_b = now_sec() - t_phase;
    fprintf(log, "=== Phase B: Clustering (%.1fs) ===\n", t_b);
    fprintf(log, "Clusters:                 %lu\n\n", (unsigned long)n_clusters);
    fprintf(stderr, "Phase B complete (%.1fs): %lu clusters\n\n",
            t_b, (unsigned long)n_clusters);

    /* ============================================================== */
    /*  Phase C: Cluster size distribution                            */
    /* ============================================================== */
    fprintf(stderr, "Phase C: Computing cluster size distribution...\n");

    /* Per-cluster counts — indexed 1..n_clusters (UBID 0 = unassigned) */
    uint64_t *uniq_bc_count  = xcalloc(n_clusters + 1, sizeof(uint64_t), "uniq_bc_count");
    uint64_t *read_count     = xcalloc(n_clusters + 1, sizeof(uint64_t), "read_count");

    for (uint64_t i = 0; i < ht.capacity; i++) {
        if (!ht.slots[i].occupied) continue;
        uint64_t ubid = ubid_map[i];
        if (ubid == 0 || ubid > n_clusters) continue;
        uniq_bc_count[ubid]++;
        read_count[ubid] += ht.slots[i].count;
    }

    /* Histogram by unique barcodes per cluster: sizes 1..12, then >12 */
    #define N_SIZE_BINS 13  /* bins 0..11 = sizes 1..12; bin 12 = >12 */
    long size_hist[N_SIZE_BINS] = {0};

    uint64_t min_reads = UINT64_MAX, max_reads = 0;
    double   sum_reads = 0.0;
    long     n_nonempty = 0;

    for (uint64_t u = 1; u <= n_clusters; u++) {
        uint64_t sz = uniq_bc_count[u];
        if (sz == 0) continue;
        int bin = (sz <= 12) ? (int)sz - 1 : 12;
        size_hist[bin]++;

        uint64_t rc = read_count[u];
        if (rc < min_reads) min_reads = rc;
        if (rc > max_reads) max_reads = rc;
        sum_reads += (double)rc;
        n_nonempty++;
    }

    free(uniq_bc_count);
    free(read_count);

    fprintf(log, "=== Phase C: Cluster size distribution ===\n");
    fprintf(log, "Unique BCs/cluster    Clusters\n");
    fprintf(log, "-------------------   --------\n");
    for (int b = 0; b < 12; b++)
        fprintf(log, "%-21d %ld\n", b + 1, size_hist[b]);
    fprintf(log, "%-21s %ld\n", ">12", size_hist[12]);
    fprintf(log, "\n");
    if (n_nonempty > 0) {
        if (min_reads == UINT64_MAX) min_reads = 0;
        fprintf(log, "Read count per cluster: min=%lu, max=%lu, mean=%.1f\n\n",
                (unsigned long)min_reads, (unsigned long)max_reads,
                sum_reads / n_nonempty);
    }

    fprintf(stderr, "Phase C complete\n\n");

    /* ============================================================== */
    /*  Phase D: Ground truth evaluation                              */
    /* ============================================================== */
    fprintf(stderr, "Phase D: Ground truth evaluation...\n");
    t_phase = now_sec();

    long pair_cap  = 1024 * 1024;
    long pair_n    = 0;
    read_pair_t *pairs = xmalloc(pair_cap * sizeof(read_pair_t), "pairs");

    long n_no_truth = 0;
    long n_no_ubid  = 0;

    for (int fi = 0; fi < params.n_in; fi++) {
        gz_reader_t *reader = gz_reader_open(params.in_files[fi]);
        while (read_fastq(reader, hdr, seq, plus, qual)) {
            char true_bc[MAX_BC_LEN + 1];
            if (parse_true_bc(hdr, true_bc) < 0) {
                n_no_truth++;
                continue;
            }

            uint64_t ubid = 0;
            extraction_result_t ext;
            int slen = (int)strlen(seq);
            if (extract_barcode(seq, slen, &params.extract, &ext) >= 0) {
                barcode_t bc;
                int nn = encode_barcode(ext.barcode, ext.bc_len, &bc);
                if (nn >= 0 && nn <= params.max_n) {
                    if (nn == 0) {
                        bc_key_t key;
                        barcode_to_key(&bc, &key);
                        ht_entry_t *e = ht_lookup(&ht, &key);
                        if (e) {
                            uint64_t slot = (uint64_t)(e - ht.slots);
                            ubid = ubid_map[slot];
                        }
                    } else {
                        ubid = resolve_n_barcode(&bc, &ht, ubid_map);
                    }
                }
            }

            if (ubid == 0) {
                n_no_ubid++;
                continue;
            }

            if (pair_n >= pair_cap) {
                pair_cap *= 2;
                pairs = xrealloc(pairs, pair_cap * sizeof(read_pair_t), "pairs grow");
            }
            strncpy(pairs[pair_n].true_bc, true_bc, MAX_BC_LEN);
            pairs[pair_n].true_bc[MAX_BC_LEN] = '\0';
            pairs[pair_n].ubid = ubid;
            pair_n++;
        }
        gz_reader_close(reader);
    }

    /* --- Splitting: one true barcode → multiple clusters --- */
    /* Sort by (true_bc, ubid) so distinct UBIDs within a true_bc run are adjacent */
    qsort(pairs, pair_n, sizeof(read_pair_t), cmp_by_true_bc_ubid);

    long n_true_bcs = 0;
    long split_hist[6] = {0};  /* index k = (k+1) distinct UBIDs; index 5 = 6+ */

    for (long i = 0; i < pair_n; ) {
        long j = i;
        while (j < pair_n && strcmp(pairs[j].true_bc, pairs[i].true_bc) == 0)
            j++;
        /* Count distinct UBIDs in [i, j): entries are sorted by ubid */
        int n_distinct = 1;
        for (long k = i + 1; k < j; k++) {
            if (pairs[k].ubid != pairs[k - 1].ubid)
                n_distinct++;
        }
        n_true_bcs++;
        int bin = (n_distinct <= 5) ? n_distinct - 1 : 5;
        split_hist[bin]++;
        i = j;
    }

    /* --- Merging: one cluster ← multiple true barcodes --- */
    /* Sort by (ubid, true_bc) so distinct true_bcs within a ubid run are adjacent */
    qsort(pairs, pair_n, sizeof(read_pair_t), cmp_by_ubid_true_bc);

    long n_pop_clusters = 0;
    long merge_hist[6] = {0};  /* index k = (k+1) distinct true barcodes; index 5 = 6+ */

    for (long i = 0; i < pair_n; ) {
        long j = i;
        while (j < pair_n && pairs[j].ubid == pairs[i].ubid)
            j++;
        /* Count distinct true_bcs in [i, j): entries are sorted by true_bc */
        int n_distinct = 1;
        for (long k = i + 1; k < j; k++) {
            if (strcmp(pairs[k].true_bc, pairs[k - 1].true_bc) != 0)
                n_distinct++;
        }
        n_pop_clusters++;
        int bin = (n_distinct <= 5) ? n_distinct - 1 : 5;
        merge_hist[bin]++;
        i = j;
    }

    double t_d = now_sec() - t_phase;
    fprintf(log, "=== Phase D: Ground truth evaluation (%.1fs) ===\n", t_d);
    fprintf(log, "Reads missing 'barcode=' tag:  %ld\n", n_no_truth);
    fprintf(log, "Reads without cluster assignment: %ld\n\n", n_no_ubid);

    /* Splitting report */
    long split_count = 0;
    for (int b = 1; b < 6; b++) split_count += split_hist[b];

    fprintf(log, "--- Splitting (true barcode spread across >1 cluster) ---\n");
    fprintf(log, "Total true barcodes seen: %ld\n", n_true_bcs);
    for (int b = 0; b < 5; b++)
        fprintf(log, "  %d cluster(s):  %ld true barcodes\n", b + 1, split_hist[b]);
    fprintf(log, "  6+ clusters:   %ld true barcodes\n", split_hist[5]);
    if (n_true_bcs > 0)
        fprintf(log, "Split rate: %.4f%% (%ld / %ld)\n\n",
                100.0 * split_count / n_true_bcs, split_count, n_true_bcs);
    else
        fprintf(log, "Split rate: N/A\n\n");

    /* Merging report */
    long merge_count = 0;
    for (int b = 1; b < 6; b++) merge_count += merge_hist[b];

    fprintf(log, "--- Merging (cluster contains >1 true barcode) ---\n");
    fprintf(log, "Total populated clusters: %ld\n", n_pop_clusters);
    for (int b = 0; b < 5; b++)
        fprintf(log, "  %d true barcode(s): %ld clusters\n", b + 1, merge_hist[b]);
    fprintf(log, "  6+ true barcodes:   %ld clusters\n", merge_hist[5]);
    if (n_pop_clusters > 0)
        fprintf(log, "Merge rate: %.4f%% (%ld / %ld)\n\n",
                100.0 * merge_count / n_pop_clusters, merge_count, n_pop_clusters);
    else
        fprintf(log, "Merge rate: N/A\n\n");

    /* Summary */
    double t_total = now_sec() - t_start;
    fprintf(log, "=== Summary ===\n");
    fprintf(log, "Total time:               %.1fs\n", t_total);
    fprintf(log, "Log:                      %s\n", log_path);

    fprintf(stderr, "Phase D complete (%.1fs)\n", t_d);
    fprintf(stderr, "\nDone. Log: %s\n", log_path);

    /* Cleanup */
    free(pairs);
    ht_free(&ht);
    free(ubid_map);
    fclose(log);
    for (int i = 0; i < params.n_in; i++) free(params.in_files[i]);
    free(params.in_files);

    return 0;
}
