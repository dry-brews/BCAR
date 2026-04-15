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
#include <math.h>
#include <pthread.h>

#include "sort_module.h"
#include "seq_module.h"

/* Tell the linker we use fatal_alloc from sort_module.c */
#define FATAL_ALLOC_PROVIDED

/* ------------------------------------------------------------------ */
/*  Program parameters                                                */
/* ------------------------------------------------------------------ */

typedef struct {
    /* Sort parameters */
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
    /* Consensus parameters */
    int    num_threads;
    bool   no_alignment;
} params_t;

/* ------------------------------------------------------------------ */
/*  Globals                                                           */
/* ------------------------------------------------------------------ */

/* cleanup state */
static char g_tmp_dir[PATH_MAX] = {0};
static char **g_chunk_paths = NULL;
static int    g_n_chunks = 0;

/* consensus globals (used by seq_module) */
double gap_pen = DEFAULT_GAP_SCORE;
int max_line_len = DEFAULT_MAX_LINE_LEN;
double call_total = 0.0;
double call_missense = 0.0;
double call_indel = 0.0;
long total_barcodes = 0;

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
        "Barcode extraction (choose mode):\n"
        "  --bc-start <int>                 Fixed position start (0-based, default: 0)\n"
        "  --bc-len <int>                   Fixed barcode length (default: 18)\n"
        "  --context <pattern>              Context pattern with N's for barcode positions\n"
        "  --max-context-mismatches <int>   Allowed mismatches in flanking seqs (default: 1)\n\n"
        "Clustering:\n"
        "  --max-mismatches <int>           Hamming distance radius (default: 1)\n"
        "  --max-bc-indels <int>            Max indels in clustering/extraction (0-2, default: 0)\n"
        "  --max-n <int>                    Max N bases in barcode (default: 2)\n\n"
        "Consensus:\n"
        "  --threads <int>                  Number of worker threads (default: 1)\n"
        "  --gap <float>                    Gap penalty for alignment (default: -1.0)\n"
        "  --max-len <int>                  Maximum read length (default: 100,000)\n"
        "  --no-alignment                   Skip alignment\n\n"
        "Output:\n"
        "  --out <file>                     Output file (default: bcar_out.fastq)\n"
        "  --fail <file>                    Reads where barcode not found (if context used)\n\n"
        "Resources:\n"
        "  --chunk-mem <bytes>              Sort buffer size (default: 4GB)\n"
        "  --max-open-files <int>           Max files during merge (default: 240)\n"
        "  --tmp-dir <path>                 Temp directory (default: .)\n\n",
        prog);
    exit(EXIT_FAILURE);
}

static params_t parse_args(int argc, char *argv[]) {
    params_t p;
    memset(&p, 0, sizeof(p));
    p.out_file = "bcar_out.fastq";
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
    p.num_threads = 1;
    p.no_alignment = false;

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
        {"threads",                required_argument, 0, 'T'},
        {"gap",                    required_argument, 0, 'g'},
        {"max-len",                required_argument, 0, 'M'},
        {"no-alignment",           no_argument,       0, 'N'},
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
            case 'T': p.num_threads = atoi(optarg); break;
            case 'g': gap_pen = atof(optarg); break;
            case 'M': max_line_len = atoi(optarg); break;
            case 'N': p.no_alignment = true; break;
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

    if (p.num_threads < 1) {
        fprintf(stderr, "Error: --threads must be >= 1\n");
        exit(EXIT_FAILURE);
    }

    return p;
}

/* ------------------------------------------------------------------ */
/*  Threading structures (from bc_merger)                             */
/* ------------------------------------------------------------------ */

typedef struct WorkQueueNode {
    BarcodeGroup *group;
    struct WorkQueueNode *next;
} WorkQueueNode;

typedef struct {
    WorkQueueNode *head;
    WorkQueueNode *tail;
    int size;
    pthread_mutex_t mutex;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
    bool producer_done;
} WorkQueue;

typedef struct {
    FILE *fwd_out;
    FILE *fail_out;
    pthread_mutex_t mutex;
} OutputFiles;

typedef struct {
    WorkQueue *queue;
    OutputFiles *output;
    int thread_id;
    bool no_alignment;
    const extraction_config_t *extract;
} WorkerArgs;

/* Context passed to merge callback */
typedef struct {
    WorkQueue *queue;
} merge_context_t;

/* ------------------------------------------------------------------ */
/*  Work queue operations                                             */
/* ------------------------------------------------------------------ */

static void init_work_queue(WorkQueue *queue) {
    queue->head = NULL;
    queue->tail = NULL;
    queue->size = 0;
    queue->producer_done = false;
    pthread_mutex_init(&queue->mutex, NULL);
    pthread_cond_init(&queue->not_empty, NULL);
    pthread_cond_init(&queue->not_full, NULL);
}

static void queue_push(WorkQueue *queue, BarcodeGroup *group) {
    WorkQueueNode *node = xmalloc(sizeof(WorkQueueNode), "queue_push node");
    node->group = group;
    node->next = NULL;

    pthread_mutex_lock(&queue->mutex);
    while (queue->size >= 10000) {
        pthread_cond_wait(&queue->not_full, &queue->mutex);
    }
    if (queue->tail) {
        queue->tail->next = node;
    } else {
        queue->head = node;
    }
    queue->tail = node;
    queue->size++;
    pthread_cond_signal(&queue->not_empty);
    pthread_mutex_unlock(&queue->mutex);
}

static BarcodeGroup *queue_pop(WorkQueue *queue) {
    pthread_mutex_lock(&queue->mutex);
    while (queue->size == 0) {
        if (queue->producer_done) {
            pthread_mutex_unlock(&queue->mutex);
            return NULL;
        }
        pthread_cond_wait(&queue->not_empty, &queue->mutex);
    }
    WorkQueueNode *node = queue->head;
    BarcodeGroup *group = node->group;
    queue->head = node->next;
    if (!queue->head) queue->tail = NULL;
    queue->size--;
    pthread_cond_signal(&queue->not_full);
    pthread_mutex_unlock(&queue->mutex);
    free(node);
    return group;
}

/* ------------------------------------------------------------------ */
/*  Consensus processing (from bc_merger)                             */
/* ------------------------------------------------------------------ */

static int cmp_idxsim_desc(const void *a, const void *b) {
    const IdxSim *A = a;
    const IdxSim *B = b;
    if (A->sim < B->sim) return 1;
    if (A->sim > B->sim) return -1;
    return 0;
}

static void free_barcode_group(BarcodeGroup *group) {
    for (int i = 0; i < group->count; i++) {
        free(group->fwd_reads[i]);
        free(group->fwd_quals[i]);
    }
    free(group->fwd_reads);
    free(group->fwd_quals);
    free(group);
}

/* Compute minimum and mean Phred quality score (integers) across a quality string */
static void compute_q_stats(const char *qual, int *min_q, int *mean_q) {
    int len = (int)strlen(qual);
    if (len <= 0) { *min_q = 0; *mean_q = 0; return; }
    int mn = 255;
    long sum = 0;
    for (int i = 0; i < len; i++) {
        int q = (int)(unsigned char)qual[i] - 33;
        if (q < 0) q = 0;
        if (q < mn) mn = q;
        sum += q;
    }
    *min_q = mn;
    *mean_q = (int)(sum / len);
}

/* Build a formatted FASTQ entry string for the consensus output */
static char *build_fwd_entry(const char *bc_detected, const char *bcid,
                             int count, int min_q, int mean_q,
                             double minor_mean, double minor_max,
                             const char *seq, const char *qual) {
    const char *fmt = "@bc=%s;bcid=%s;count=%d;min_q=%d;mean_q=%d;minor_frac_mean=%.3f;minor_frac_max=%.3f\n%s\n+\n%s\n";
    int len = snprintf(NULL, 0, fmt, bc_detected, bcid, count, min_q, mean_q,
                       minor_mean, minor_max, seq, qual) + 1;
    char *entry = malloc(len);
    snprintf(entry, len, fmt, bc_detected, bcid, count, min_q, mean_q,
             minor_mean, minor_max, seq, qual);
    return entry;
}

/* Build a fail entry in the same 4-line FASTQ format used during Phase C */
static char *build_fail_entry(const char *bcid, int count,
                              const char *seq, const char *qual) {
    const char *fmt = "@bcid=%s;count=%d\n%s\n+\n%s\n";
    int len = snprintf(NULL, 0, fmt, bcid, count, seq, qual) + 1;
    char *entry = malloc(len);
    snprintf(entry, len, fmt, bcid, count, seq, qual);
    return entry;
}

static void process_barcode(BarcodeGroup *group, char **fwd_entry, char **fail_entry,
                            bool no_alignment, const extraction_config_t *extract) {
    *fwd_entry = NULL;
    *fail_entry = NULL;

    if (group->count == 0) return;
    if (group->count == 1) {
        const char *seq = group->fwd_reads[0];
        const char *qual = group->fwd_quals[0];

        extraction_result_t ext;
        if (extract_barcode(seq, (int)strlen(seq), extract, &ext) < 0) {
            *fail_entry = build_fail_entry(group->bc, group->count, seq, qual);
            return;
        }

        int min_q = 0, mean_q = 0;
        compute_q_stats(qual, &min_q, &mean_q);

        *fwd_entry = build_fwd_entry(ext.barcode, group->bc, group->count,
                                     min_q, mean_q, 0.0, 0.0, seq, qual);
        total_barcodes++;
        return;
    }

    /* Get max length of reads */
    int L_fwd = 0;
    for (int i = 0; i < group->count; i++) {
        int curr_fwd_len = strlen(group->fwd_reads[i]);
        if (curr_fwd_len > L_fwd) L_fwd = curr_fwd_len;
    }

    /* Convert reads to arrays */
    SeqArray *fwd_arrays = malloc(group->count * sizeof(SeqArray));
    if (!fwd_arrays) fatal_alloc("Failed to allocate fwd_arrays");
    for (int i = 0; i < group->count; i++) {
        fwd_arrays[i] = seq_to_array(group->fwd_reads[i], group->fwd_quals[i], L_fwd);
    }

    /* Build unaligned consensus */
    SeqArray fwd_consensus = build_unaligned_consensus(fwd_arrays, group->count);

    if (no_alignment) {
        char *fwd_seq = NULL, *fwd_qual = NULL;
        double minor_max = 0.0, minor_mean = 0.0;
        array_to_seq(&fwd_consensus, &fwd_seq, &fwd_qual, &minor_max, &minor_mean);

        extraction_result_t ext;
        if (extract_barcode(fwd_seq, (int)strlen(fwd_seq), extract, &ext) < 0) {
            *fail_entry = build_fail_entry(group->bc, group->count, fwd_seq, fwd_qual);
        } else {
            int min_q = 0, mean_q = 0;
            compute_q_stats(fwd_qual, &min_q, &mean_q);
            *fwd_entry = build_fwd_entry(ext.barcode, group->bc, group->count,
                                         min_q, mean_q, minor_mean, minor_max,
                                         fwd_seq, fwd_qual);
            total_barcodes++;
        }

        free(fwd_seq);
        free(fwd_qual);
        free_seq_array(&fwd_consensus);
        for (int i = 0; i < group->count; i++)
            free_seq_array(&fwd_arrays[i]);
        free(fwd_arrays);
        return;
    }

    /* Compare reads and sort by similarity to consensus */
    double *similarities = malloc(group->count * sizeof(double));
    if (!similarities) fatal_alloc("Failed to allocate similarities");
    IdxSim *pairs = malloc(group->count * sizeof(IdxSim));
    if (!pairs) fatal_alloc("Failed to allocate pairs");

    for (int i = 0; i < group->count; i++) {
        similarities[i] = compare_seqs(&fwd_arrays[i], &fwd_consensus);
        pairs[i].idx = i;
        pairs[i].sim = similarities[i];
    }

    qsort(pairs, group->count, sizeof(IdxSim), cmp_idxsim_desc);

    int *indices = malloc(group->count * sizeof(int));
    if (!indices) fatal_alloc("Failed to allocate indices");
    for (int i = 0; i < group->count; ++i) indices[i] = pairs[i].idx;
    free(pairs);

    /* Merge sequences into consensus */
    SeqArray merged_fwd = create_seq_array(fwd_arrays[indices[0]].length);
    memcpy(merged_fwd.positions, fwd_arrays[indices[0]].positions,
           fwd_arrays[indices[0]].length * sizeof(Position));

    for (int i = 1; i < group->count; i++) {
        SeqArray temp_fwd = merge_seqs(&merged_fwd, &fwd_arrays[indices[i]]);
        free_seq_array(&merged_fwd);
        merged_fwd = temp_fwd;
    }

    /* Convert to FASTQ */
    char *fwd_seq, *fwd_qual;
    double minor_max = 0.0, minor_mean = 0.0;
    array_to_seq(&merged_fwd, &fwd_seq, &fwd_qual, &minor_max, &minor_mean);

    extraction_result_t ext;
    if (extract_barcode(fwd_seq, (int)strlen(fwd_seq), extract, &ext) < 0) {
        *fail_entry = build_fail_entry(group->bc, group->count, fwd_seq, fwd_qual);
    } else {
        int min_q = 0, mean_q = 0;
        compute_q_stats(fwd_qual, &min_q, &mean_q);
        *fwd_entry = build_fwd_entry(ext.barcode, group->bc, group->count,
                                     min_q, mean_q, minor_mean, minor_max,
                                     fwd_seq, fwd_qual);
        total_barcodes++;
    }

    free(fwd_seq);
    free(fwd_qual);
    free(similarities);
    free(indices);
    free_seq_array(&fwd_consensus);
    free_seq_array(&merged_fwd);

    for (int i = 0; i < group->count; i++)
        free_seq_array(&fwd_arrays[i]);
    free(fwd_arrays);
}

/* ------------------------------------------------------------------ */
/*  Worker thread                                                     */
/* ------------------------------------------------------------------ */

static void *worker_thread(void *arg) {
    WorkerArgs *args = (WorkerArgs *)arg;
    WorkQueue *queue = args->queue;
    OutputFiles *output = args->output;
    bool no_alignment = args->no_alignment;
    const extraction_config_t *extract = args->extract;

    while (true) {
        BarcodeGroup *group = queue_pop(queue);
        if (!group) break;

        char *fwd_entry = NULL;
        char *fail_entry = NULL;
        process_barcode(group, &fwd_entry, &fail_entry, no_alignment, extract);
        if (fwd_entry || (fail_entry && output->fail_out)) {
            pthread_mutex_lock(&output->mutex);
            if (fwd_entry) fprintf(output->fwd_out, "%s", fwd_entry);
            if (fail_entry && output->fail_out) fprintf(output->fail_out, "%s", fail_entry);
            pthread_mutex_unlock(&output->mutex);
        }
        free(fwd_entry);
        free(fail_entry);
        free_barcode_group(group);
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  Merge callback: converts merge output to BarcodeGroup             */
/* ------------------------------------------------------------------ */

static void merge_group_handler(const char *ubid_str,
                                char **seqs, char **quals,
                                int count, void *user_data) {
    merge_context_t *ctx = (merge_context_t *)user_data;

    BarcodeGroup *group = xmalloc(sizeof(BarcodeGroup), "merge_group_handler group");
    strncpy(group->bc, ubid_str, MAX_BARCODE_LEN - 1);
    group->bc[MAX_BARCODE_LEN - 1] = '\0';
    group->count = count;
    group->capacity = count;
    /* Transfer ownership of the string arrays */
    group->fwd_reads = seqs;
    group->fwd_quals = quals;

    queue_push(ctx->queue, group);
}

/* ------------------------------------------------------------------ */
/*  Main                                                              */
/* ------------------------------------------------------------------ */

int main(int argc, char *argv[]) {
    params_t params = parse_args(argc, argv);
    double t_start = now_sec();

    /* Print run info */
    fprintf(stderr, "\n=== BCAR Pipeline ===\n");
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
    fprintf(stderr, "Consensus: threads=%d, no_alignment=%s\n",
            params.num_threads, params.no_alignment ? "yes" : "no");
    fprintf(stderr, "Chunk memory: %ld bytes\n", params.chunk_mem);
    fprintf(stderr, "\n");

    /* Create temp directory */
    char tmp_template[PATH_MAX];
    snprintf(tmp_template, PATH_MAX, "%s/bcar_XXXXXX", params.tmp_dir);
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
            if (nn < 0) { n_failed++; continue; }
            if (nn > params.max_n) { n_failed++; continue; }
            if (nn > 0) { n_with_n++; continue; }

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

    union_find_t uf;
    uf_init(&uf, (uint32_t)ht.count);
    {
        uint32_t idx = 0;
        for (uint64_t i = 0; i < ht.capacity; i++) {
            if (ht.slots[i].occupied)
                ht.slots[i].uf_idx = idx++;
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
                bc_key_t key;
                barcode_to_key(&bc, &key);
                ht_entry_t *e = ht_lookup(&ht, &key);
                if (e) {
                    uint64_t slot_idx = (uint64_t)(e - ht.slots);
                    ubid = ubid_map[slot_idx];
                }
            } else {
                ubid = resolve_n_barcode(&bc, &ht, ubid_map);
            }

            if (ubid == 0) {
                if (fail_fp)
                    fprintf(fail_fp, "%s\n%s\n%s\n%s\n", hdr, seq, plus, qual);
                fail_reads++;
                continue;
            }

            char new_hdr[MAX_LINE + 32];
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
    /*  Phase D: Merge + Consensus (combined)                         */
    /* ============================================================== */
    fprintf(stderr, "Phase D: Merging sorted chunks and building consensus...\n");
    t_phase = now_sec();

    /* Set up output */
    OutputFiles output = {NULL, NULL, PTHREAD_MUTEX_INITIALIZER};
    output.fwd_out = fopen(params.out_file, "w");
    if (!output.fwd_out) {
        fprintf(stderr, "Fatal: cannot open output file: %s\n", params.out_file);
        exit(EXIT_FAILURE);
    }
    if (params.fail_file) {
        output.fail_out = fopen(params.fail_file, "a");
        if (!output.fail_out) {
            fprintf(stderr, "Fatal: cannot open fail file for append: %s\n", params.fail_file);
            exit(EXIT_FAILURE);
        }
    }

    /* Set up work queue */
    WorkQueue work_queue;
    init_work_queue(&work_queue);

    /* Create worker threads */
    pthread_t *threads = xmalloc(params.num_threads * sizeof(pthread_t), "threads");
    WorkerArgs *thread_args = xmalloc(params.num_threads * sizeof(WorkerArgs), "thread_args");

    for (int i = 0; i < params.num_threads; i++) {
        thread_args[i].queue = &work_queue;
        thread_args[i].output = &output;
        thread_args[i].thread_id = i;
        thread_args[i].no_alignment = params.no_alignment;
        thread_args[i].extract = &params.extract;
        if (pthread_create(&threads[i], NULL, worker_thread, &thread_args[i]) != 0) {
            fprintf(stderr, "Fatal: failed to create worker thread %d\n", i);
            exit(EXIT_FAILURE);
        }
    }

    /* Run merge with callback that feeds groups into the work queue */
    merge_context_t merge_ctx = { .queue = &work_queue };
    merge_sorted_chunks_with_callback(chunk_paths, n_chunks,
                                      params.max_open_files, g_tmp_dir,
                                      merge_group_handler, &merge_ctx);

    /* Signal workers that producer is done */
    pthread_mutex_lock(&work_queue.mutex);
    work_queue.producer_done = true;
    pthread_cond_broadcast(&work_queue.not_empty);
    pthread_mutex_unlock(&work_queue.mutex);

    /* Wait for workers to finish */
    for (int i = 0; i < params.num_threads; i++)
        pthread_join(threads[i], NULL);

    /* Cleanup threading */
    pthread_mutex_destroy(&work_queue.mutex);
    pthread_cond_destroy(&work_queue.not_empty);
    pthread_cond_destroy(&work_queue.not_full);
    pthread_mutex_destroy(&output.mutex);
    fclose(output.fwd_out);
    if (output.fail_out) fclose(output.fail_out);
    free(threads);
    free(thread_args);

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
    fprintf(stderr, "\n=== Pipeline complete ===\n");
    fprintf(stderr, "Total time: %s\n", elapsed);
    fprintf(stderr, "Total reads: %ld\n", total_reads);
    fprintf(stderr, "Assigned: %ld (%.1f%%)\n", assigned_reads,
            total_reads > 0 ? 100.0 * assigned_reads / total_reads : 0.0);
    fprintf(stderr, "Clusters: %lu\n", (unsigned long)n_clusters);
    fprintf(stderr, "Consensus barcodes: %ld\n", total_barcodes);
    if (call_total > 0) {
        fprintf(stderr, "%% Indels: %.4f%%\n", (call_indel / call_total) * 100.0);
        fprintf(stderr, "%% Missense: %.4f%%\n", (call_missense / call_total) * 100.0);
    }
    fprintf(stderr, "Output: %s\n\n", params.out_file);

    return 0;
}
