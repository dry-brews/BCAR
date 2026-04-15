#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <pthread.h>
#include <getopt.h>
#include "seq_module.h"

// Structures for threading (not shared with seq_module)
typedef struct {
    double** data;
    int rows;
    int cols;
} Matrix;

typedef struct WorkQueueNode {
    BarcodeGroup* group;
    struct WorkQueueNode* next;
} WorkQueueNode;

typedef struct {
    WorkQueueNode* head;
    WorkQueueNode* tail;
    int size;
    pthread_mutex_t mutex;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
    bool producer_done;
} WorkQueue;

typedef struct {
    FILE* fwd_out;
    pthread_mutex_t mutex;
} OutputFiles;

typedef struct {
    WorkQueue* queue;
    OutputFiles* output;
    int thread_id;
    bool paired_end_mode;
    bool no_alignment;
} WorkerArgs;

// Global variables
double gap_pen = DEFAULT_GAP_SCORE;
int max_line_len = DEFAULT_MAX_LINE_LEN;
double call_total = 0.0;
double call_missense = 0.0;
double call_indel = 0.0;
double mix_frac_thresh = 0.15; // min fraction of evidence for 2nd base to flag a position
int mix_pos_thresh = 3;        // min flagged positions to call a read as mixed
long mixed_count = 0;
long total_barcodes = 0;
FILE* g_fwd_file = NULL;
OutputFiles* g_output_ptr = NULL;
pthread_t* g_threads_ptr = NULL;
WorkerArgs* g_thread_args_ptr = NULL;

// Function prototypes
void init_work_queue(WorkQueue* queue);
void queue_push(WorkQueue* queue, BarcodeGroup* group);
BarcodeGroup* queue_pop(WorkQueue* queue);
void* worker_thread(void* arg);
void process_barcode(BarcodeGroup* group, char** fwd_entry, bool no_alignment);
void free_barcode_group(BarcodeGroup* group);
void cleanup_on_error(FILE* fwd_file, OutputFiles* output, pthread_t* threads, WorkerArgs* thread_args);
void print_usage(const char* program_name);
static int cmp_idxsim_desc(const void* a, const void* b);

// Main
int main(int argc, char *argv[]) {
    // Default values
    const char *in1_file = NULL;
    const char *out1_file = "seq_merge_out1.fastq";
    int num_threads = 1;
    bool no_alignment = false;
    pthread_t* threads = NULL;
    WorkerArgs* thread_args = NULL;
    OutputFiles output = {NULL, PTHREAD_MUTEX_INITIALIZER};

    WorkQueue work_queue;
    init_work_queue(&work_queue);

    // Command line options
    static struct option long_options[] = {
        {"in", required_argument, 0, 'a'},
        {"max-len", required_argument, 0, 'm'},
        {"gap", required_argument, 0, 'g'},
        {"out", required_argument, 0, '1'},
        {"threads", required_argument, 0, 't'},
        {"no-alignment", no_argument, 0, 'n'},
        {"mix-frac", required_argument, 0, 'f'},
        {"mix-pos", required_argument, 0, 'p'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;
    while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
        switch (c) {
            case 'a': in1_file = optarg; break;
            case 'm': max_line_len = atoi(optarg); break;
            case 'g': gap_pen = atof(optarg); break;
            case '1': out1_file = optarg; break;
            case 't': num_threads = atoi(optarg); break;
            case 'n': no_alignment = true; break;
            case 'f': mix_frac_thresh = atof(optarg); break;
            case 'p': mix_pos_thresh = atoi(optarg); break;
            default: print_usage(argv[0]);
        }
    }

    // Validate arguments
    if (!in1_file) {
        fprintf(stderr, "Error: Read 1 file is required (--in)\n");
        print_usage(argv[0]);
    }

    if (num_threads < 1) {
        fprintf(stderr, "Error: Invalid thread count\n");
        exit(EXIT_FAILURE);
    }

    // Print run information
    fprintf(stderr, "\nInitializing sequence processing...\n");
    fprintf(stderr, "Threads: %d\n", num_threads);
    fprintf(stderr, "Input file(s):\n  Read 1: %s\n", in1_file);
    fprintf(stderr, "Output file(s):\n  Read 1: %s\n", out1_file);
    fprintf(stderr, "\n");

    // Initialize input files
    FILE* fwd_file = fopen(in1_file, "r");

    // Set global pointers for fatal cleanup helper
    g_fwd_file = fwd_file;
    g_output_ptr = &output;

    // Initialize output files
    output.fwd_out = fopen(out1_file, "w");

    // Create worker threads
    threads = malloc(num_threads * sizeof(pthread_t));
    if (!threads) {
        fatal_alloc("Failed to allocate thread handles");
    }
    // set global pointer now that threads pointer exists for cleanup
    g_threads_ptr = threads;

    thread_args = malloc(num_threads * sizeof(WorkerArgs));
    if (!thread_args) {
        fatal_alloc("Failed to allocate thread args");
    }
    g_thread_args_ptr = thread_args;

    for (int i = 0; i < num_threads; i++) {
        thread_args[i].queue = &work_queue;
        thread_args[i].output = &output;
        thread_args[i].thread_id = i;
        thread_args[i].no_alignment = no_alignment;
        if (pthread_create(&threads[i], NULL, worker_thread, &thread_args[i]) != 0) {
            fatal_alloc("Failed to create worker thread");
        }
    }

    // Producer: Read input files and create work
    BarcodeGroup* current_group = NULL;
    char line_fwd[max_line_len];

    while (1) {
        // Read header line
        if (!fgets(line_fwd, max_line_len, fwd_file)) {
            if (current_group) {
                queue_push(&work_queue, current_group);
            }
            break;
        }

        // Parse UBID from header (expects ";UBID=<number>" tag)
        char bc[MAX_BARCODE_LEN];
        bc[0] = '\0';
        const char* ubid_tag = strstr(line_fwd, ";UBID=");
        if (ubid_tag) {
            const char* val = ubid_tag + 6; // skip ";UBID="
            int i = 0;
            while (val[i] && val[i] != ';' && val[i] != '\n' && val[i] != '\r'
                   && i < MAX_BARCODE_LEN - 1) {
                bc[i] = val[i];
                i++;
            }
            bc[i] = '\0';
        }

        // Read sequences
        (void)fgets(line_fwd, max_line_len, fwd_file);
        char* fwd_seq = strndup(line_fwd, strcspn(line_fwd, "\n"));

        // Skip '+' lines
        (void)fgets(line_fwd, max_line_len, fwd_file);

        // Read qualities
        (void)fgets(line_fwd, max_line_len, fwd_file);
        char* fwd_qual = strndup(line_fwd, strcspn(line_fwd, "\n"));

        // Skip reads with no UBID
        if (bc[0] == '\0') {
            free(fwd_seq);
            free(fwd_qual);
            continue;
        }

        // Handle group management
        if (!current_group) {
            current_group = malloc(sizeof(BarcodeGroup));
            if (!current_group) fatal_alloc("Failed to allocate BarcodeGroup");
            strcpy(current_group->bc, bc);
            current_group->count = 0;
            current_group->capacity = 16; // initial capacity
            current_group->fwd_reads = malloc(sizeof(char*) * current_group->capacity);
            if (!current_group->fwd_reads) fatal_alloc("Failed to allocate fwd_reads");
            current_group->fwd_quals = malloc(sizeof(char*) * current_group->capacity);
            if (!current_group->fwd_quals) fatal_alloc("Failed to allocate fwd_quals");
        }

        if (strcmp(current_group->bc, bc) != 0) {
            queue_push(&work_queue, current_group);

            current_group = malloc(sizeof(BarcodeGroup));
            if (!current_group) fatal_alloc("Failed to allocate BarcodeGroup");
            strcpy(current_group->bc, bc);
            current_group->count = 0;
            current_group->capacity = 16; // initial capacity
            current_group->fwd_reads = malloc(sizeof(char*) * current_group->capacity);
            if (!current_group->fwd_reads) fatal_alloc("Failed to allocate fwd_reads");
            current_group->fwd_quals = malloc(sizeof(char*) * current_group->capacity);
            if (!current_group->fwd_quals) fatal_alloc("Failed to allocate fwd_quals");
        }

        // Add to group
        if (current_group->count >= current_group->capacity) {
            current_group->capacity *= 2;
            char** tmp = realloc(current_group->fwd_reads, sizeof(char*) * current_group->capacity);
            if (!tmp) fatal_alloc("Failed to realloc fwd_reads");
            current_group->fwd_reads = tmp;

            tmp = realloc(current_group->fwd_quals, sizeof(char*) * current_group->capacity);
            if (!tmp) fatal_alloc("Failed to realloc fwd_quals");
            current_group->fwd_quals = tmp;
        }

        current_group->fwd_reads[current_group->count] = fwd_seq;
        current_group->fwd_quals[current_group->count] = fwd_qual;
        current_group->count++;
    }

    // Signal that producer is done
    pthread_mutex_lock(&work_queue.mutex);
    work_queue.producer_done = true;
    pthread_cond_broadcast(&work_queue.not_empty);
    pthread_mutex_unlock(&work_queue.mutex);

    // Wait for all worker threads to finish
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    // Cleanup
    pthread_mutex_destroy(&work_queue.mutex);
    pthread_cond_destroy(&work_queue.not_empty);
    pthread_cond_destroy(&work_queue.not_full);
    pthread_mutex_destroy(&output.mutex);

    free(threads);
    free(thread_args);

    fclose(fwd_file);
    fclose(output.fwd_out);

    // Print error rate statistics
    fprintf(stderr, "\nError rate statistics:\n");
    fprintf(stderr, "%% Indels: %.4f%%\n", (call_indel / call_total) * 100.0);
    fprintf(stderr, "%% Missense: %.4f%%\n", (call_missense / call_total) * 100.0);
    fprintf(stderr, "%% Mixed barcodes: %.2f%%  (%ld / %ld)\n",
            total_barcodes > 0 ? (mixed_count * 100.0 / total_barcodes) : 0.0,
            mixed_count, total_barcodes);

    return 0;
}

// Get reads associated with a barcode, merge, and print
void process_barcode(BarcodeGroup* group, char** fwd_entry, bool no_alignment) {
    *fwd_entry = NULL;

    // Step 0: Check group size and don't do full processing if not necessary
    if(group->count == 0) return;
    if(group->count == 1) {
        // Format FASTQ entry
        int fwd_len = snprintf(NULL, 0, "@fwd_read;bc=%s;count=%d\n%s\n+\n%s\n",
                              group->bc, group->count, group->fwd_reads[0], group->fwd_quals[0]) + 1;

        *fwd_entry = malloc(fwd_len);
        snprintf(*fwd_entry, fwd_len, "@fwd_read;bc=%s;count=%d\n%s\n+\n%s\n",
                group->bc, group->count, group->fwd_reads[0], group->fwd_quals[0]);
        total_barcodes++;
        return;
    }

    // Step 1: Get max length of reads
    int L_fwd = 0;
    for(int i = 0; i < group->count; i++) {
        int curr_fwd_len = strlen(group->fwd_reads[i]);
        if(curr_fwd_len > L_fwd) L_fwd = curr_fwd_len;
    }

    // Step 2: Convert reads to arrays
    SeqArray* fwd_arrays = malloc(group->count * sizeof(SeqArray));
    if (!fwd_arrays) fatal_alloc("Failed to allocate fwd_arrays");
    for(int i=0; i<group->count; i++) {
        fwd_arrays[i] = seq_to_array(group->fwd_reads[i], group->fwd_quals[i], L_fwd);
    }

    // Step 3: Build unaligned consensus
    SeqArray fwd_consensus = build_unaligned_consensus(fwd_arrays, group->count);

    // If --no-alignment option is set, skip compare/sort and merge, convert consensus directly
    if (no_alignment) {
        char *fwd_seq = NULL, *fwd_qual = NULL;
        bool is_mixed = false;
        double minor_frac = 0.0;
        array_to_seq(&fwd_consensus, &fwd_seq, &fwd_qual, &is_mixed, &minor_frac);

        // Format FASTQ entry
        int fwd_len;
        if (is_mixed)
            fwd_len = snprintf(NULL, 0, "@fwd_read;bc=%s;count=%d;mixed=1;minor_frac=%.3f\n%s\n+\n%s\n",
                               group->bc, group->count, minor_frac, fwd_seq, fwd_qual) + 1;
        else
            fwd_len = snprintf(NULL, 0, "@fwd_read;bc=%s;count=%d\n%s\n+\n%s\n",
                               group->bc, group->count, fwd_seq, fwd_qual) + 1;

        *fwd_entry = malloc(fwd_len);
        if (is_mixed)
            snprintf(*fwd_entry, fwd_len, "@fwd_read;bc=%s;count=%d;mixed=1;minor_frac=%.3f\n%s\n+\n%s\n",
                     group->bc, group->count, minor_frac, fwd_seq, fwd_qual);
        else
            snprintf(*fwd_entry, fwd_len, "@fwd_read;bc=%s;count=%d\n%s\n+\n%s\n",
                     group->bc, group->count, fwd_seq, fwd_qual);

        total_barcodes++;
        if (is_mixed) mixed_count++;

        // Cleanup
        free(fwd_seq);
        free(fwd_qual);
        free_seq_array(&fwd_consensus);
        for(int i=0; i<group->count; i++) {
            free_seq_array(&fwd_arrays[i]);
        }
        free(fwd_arrays);
        return;
    }

    // Step 4: Compare reads and sort
    double* similarities = malloc(group->count * sizeof(double));
    if (!similarities) fatal_alloc("Failed to allocate similarities");
    IdxSim* pairs = malloc(group->count * sizeof(IdxSim));
    if (!pairs) fatal_alloc("Failed to allocate pairs");

    for(int i=0; i<group->count; i++) {
        similarities[i] = compare_seqs(&fwd_arrays[i], &fwd_consensus);
        pairs[i].idx = i;
        pairs[i].sim = similarities[i];
    }

    qsort(pairs, group->count, sizeof(IdxSim), cmp_idxsim_desc);

    int* indices = malloc(group->count * sizeof(int));
    if (!indices) fatal_alloc("Failed to allocate indices");
    for (int i = 0; i < group->count; ++i) indices[i] = pairs[i].idx;

    free(pairs);

    // Step 5: Merge sequences into consensus
    SeqArray merged_fwd = create_seq_array(fwd_arrays[indices[0]].length);
    memcpy(merged_fwd.positions, fwd_arrays[indices[0]].positions,
           fwd_arrays[indices[0]].length * sizeof(Position));

    for(int i=1; i<group->count; i++) {
        SeqArray temp_fwd = merge_seqs(&merged_fwd, &fwd_arrays[indices[i]]);
        free_seq_array(&merged_fwd);
        merged_fwd = temp_fwd;
    }

    // Step 6: Convert to FASTQ and print
    char *fwd_seq, *fwd_qual;
    bool is_mixed = false;
    double minor_frac = 0.0;
    array_to_seq(&merged_fwd, &fwd_seq, &fwd_qual, &is_mixed, &minor_frac);

    // Format FASTQ entry
    int fwd_len;
    if (is_mixed)
        fwd_len = snprintf(NULL, 0, "@fwd_read;bc=%s;count=%d;mixed=1;minor_frac=%.3f\n%s\n+\n%s\n",
                           group->bc, group->count, minor_frac, fwd_seq, fwd_qual) + 1;
    else
        fwd_len = snprintf(NULL, 0, "@fwd_read;bc=%s;count=%d\n%s\n+\n%s\n",
                           group->bc, group->count, fwd_seq, fwd_qual) + 1;

    *fwd_entry = malloc(fwd_len);
    if (is_mixed)
        snprintf(*fwd_entry, fwd_len, "@fwd_read;bc=%s;count=%d;mixed=1;minor_frac=%.3f\n%s\n+\n%s\n",
                 group->bc, group->count, minor_frac, fwd_seq, fwd_qual);
    else
        snprintf(*fwd_entry, fwd_len, "@fwd_read;bc=%s;count=%d\n%s\n+\n%s\n",
                 group->bc, group->count, fwd_seq, fwd_qual);

    total_barcodes++;
    if (is_mixed) mixed_count++;

    // Cleanup
    free(fwd_seq);
    free(fwd_qual);
    free(similarities);
    free(indices);
    free_seq_array(&fwd_consensus);
    free_seq_array(&merged_fwd);

    for(int i=0; i<group->count; i++) {
        free_seq_array(&fwd_arrays[i]);
    }
    free(fwd_arrays);
}

// Organize workers
// Initialize work queue
void init_work_queue(WorkQueue* queue) {
    queue->head = NULL;
    queue->tail = NULL;
    queue->size = 0;
    queue->producer_done = false;
    pthread_mutex_init(&queue->mutex, NULL);
    pthread_cond_init(&queue->not_empty, NULL);
    pthread_cond_init(&queue->not_full, NULL);
}

// Add work to queue
void queue_push(WorkQueue* queue, BarcodeGroup* group) {
    WorkQueueNode* node = malloc(sizeof(WorkQueueNode));
    if (!node) {
        fatal_alloc("Failed to allocate WorkQueueNode");
    }
    node->group = group;
    node->next = NULL;

    pthread_mutex_lock(&queue->mutex);

    // Wait if queue is too full (prevent memory explosion)
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

// Get work from queue
BarcodeGroup* queue_pop(WorkQueue* queue) {
    pthread_mutex_lock(&queue->mutex);

    while (queue->size == 0) {
        if (queue->producer_done) {
            pthread_mutex_unlock(&queue->mutex);
            return NULL;
        }
        pthread_cond_wait(&queue->not_empty, &queue->mutex);
    }

    WorkQueueNode* node = queue->head;
    BarcodeGroup* group = node->group;
    queue->head = node->next;
    if (!queue->head) {
        queue->tail = NULL;
    }
    queue->size--;

    pthread_cond_signal(&queue->not_full);
    pthread_mutex_unlock(&queue->mutex);

    free(node);
    return group;
}

// Worker thread function
void* worker_thread(void* arg) {
    WorkerArgs* args = (WorkerArgs*)arg;
    WorkQueue* queue = args->queue;
    OutputFiles* output = args->output;
    bool no_alignment = args->no_alignment;

    while (true) {
        BarcodeGroup* group = queue_pop(queue);
        if (!group) {
            break;  // Queue is empty and producer is done
        }
        char *fwd_entry;
        process_barcode(group, &fwd_entry, no_alignment);
        if (fwd_entry) {
            pthread_mutex_lock(&output->mutex);
            fprintf(output->fwd_out, "%s", fwd_entry);
            pthread_mutex_unlock(&output->mutex);
        }
        free(fwd_entry);
        free_barcode_group(group);
    }
    return NULL;
}

// Helper functions
static int cmp_idxsim_desc(const void* a, const void* b) {
    const IdxSim* A = a;
    const IdxSim* B = b;
    if (A->sim < B->sim) return 1;
    if (A->sim > B->sim) return -1;
    return 0;
}

/* fatal_alloc is provided by sort_module.c when linked together,
   or defined here as a fallback for standalone builds. */
#ifndef FATAL_ALLOC_PROVIDED
void fatal_alloc(const char* where) {
    fprintf(stderr, "Fatal error: %s\n", where);
    cleanup_on_error(g_fwd_file, g_output_ptr, g_threads_ptr, g_thread_args_ptr);
}
#endif

void free_barcode_group(BarcodeGroup* group) {
    for (int i = 0; i < group->count; i++) {
        free(group->fwd_reads[i]);
        free(group->fwd_quals[i]);
    }
    free(group->fwd_reads);
    free(group->fwd_quals);
    free(group);
}

void cleanup_on_error(FILE* fwd_file, OutputFiles* output,
    pthread_t* threads, WorkerArgs* thread_args) {
    if (fwd_file) fclose(fwd_file);
    if (output->fwd_out) fclose(output->fwd_out);
    pthread_mutex_destroy(&output->mutex);
    free(threads);
    free(thread_args);
    exit(EXIT_FAILURE);
}

// Print usage
void print_usage(const char* program_name) {
    fprintf(stderr, "Usage: %s [options]\n", program_name);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  --in file          Input FASTQ file (required, sorted by UBID)\n");
    fprintf(stderr, "  --max-len int      Maximum read length (default: 131072)\n");
    fprintf(stderr, "  --gap float        Gap penalty for alignment (default: -1.0)\n");
    fprintf(stderr, "  --out file         Output file for consensus read 1\n");
    fprintf(stderr, "  --threads int      Number of threads (default: 1)\n");
    fprintf(stderr, "  --no-alignment     Skip alignment/merging; use unaligned consensus directly\n");
    fprintf(stderr, "  --mix-frac float   Min secondary-base fraction to flag a position as mixed (default: 0.15)\n");
    fprintf(stderr, "  --mix-pos int      Min flagged positions to mark a read as mixed (default: 3)\n");
    fprintf(stderr, "\nInput must have ;UBID=<id> tags in headers (from bcar_sort).\n");
    exit(EXIT_FAILURE);
}
