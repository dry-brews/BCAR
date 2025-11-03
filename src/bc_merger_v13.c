#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <pthread.h>
#include <getopt.h>

#define MAX_BARCODE_LEN 128
#define DEFAULT_MAX_LINE_LEN 131072
#define BASES "ACGT-"
#define DEFAULT_GAP_SCORE -1.0
#define VECTOR_LENGTH 5
#define q_correct -0.004575749 // = 1/10 * log10(0.9)
#define q_error -0.1 // = 1/10 * log10(0.1)

// Structures
// Structure for barcode group (all reads associated with one barcode)
typedef struct {
    char bc[MAX_BARCODE_LEN];
    char** fwd_reads;
    char** fwd_quals;
    int count;
    int capacity; // added capacity field
} BarcodeGroup;

// Structure for position (base scores for A, C, G, T, and gap)
typedef struct {
    int scores[5]; // A, C, G, T, -
} Position;

// Structure for sequence array (array of positions)
typedef struct {
    Position* positions;
    int length;
} SeqArray;

// Structure for alignment matrices
typedef struct {
    double** data;
    int rows;
    int cols;
} Matrix;

// Queue structure for work items
typedef struct WorkQueueNode {
    BarcodeGroup* group;
    struct WorkQueueNode* next;
} WorkQueueNode;

// Structure for work queue
typedef struct {
    WorkQueueNode* head;
    WorkQueueNode* tail;
    int size;
    pthread_mutex_t mutex;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
    bool producer_done;
} WorkQueue;

// Structure to hold output file handles
typedef struct {
    FILE* fwd_out;
    pthread_mutex_t mutex;
} OutputFiles;

// Structure to pass to worker threads
typedef struct {
    WorkQueue* queue;
    OutputFiles* output;
    int thread_id;
    bool paired_end_mode;
    bool no_alignment;
} WorkerArgs;

// Structure for sorting reads by similarity
typedef struct { int idx; double sim; } IdxSim;

// Global variables
const char* bases = BASES;
const int basemap[256] = { ['A']=0, ['C']=1, ['G']=2, ['T']=3, ['-']=4, ['N']=4 };
const int q_adj[41] = { // adjusted quality score mapping
    0, 0, 0, 0, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40, 41, 42
};
double gap_pen = DEFAULT_GAP_SCORE;
int max_line_len = DEFAULT_MAX_LINE_LEN;
double call_total = 0.0;
double call_missense = 0.0;
double call_indel = 0.0;
FILE* g_fwd_file = NULL;
OutputFiles* g_output_ptr = NULL;
pthread_t* g_threads_ptr = NULL;
WorkerArgs* g_thread_args_ptr = NULL;

// Function prototypes
// organize workers into a queue
void init_work_queue(WorkQueue* queue);
void queue_push(WorkQueue* queue, BarcodeGroup* group);
BarcodeGroup* queue_pop(WorkQueue* queue);
void* worker_thread(void* arg);
// process barcode group
void process_barcode(BarcodeGroup* group, char** fwd_entry, bool no_alignment);
SeqArray seq_to_array(const char* seq, const char* qual, int length);
double compare_positions(const Position* a, const Position* b);
double compare_seqs(const SeqArray* a, const SeqArray* b);
SeqArray build_unaligned_consensus(SeqArray* sequences, int count);
int* align_arrays(const SeqArray* ref, const SeqArray* query, int *out_len);
int* align_arrays_band(const SeqArray* ref, const SeqArray* query, int *out_len, int max_phase_diff);
SeqArray merge_seqs(const SeqArray* seq1, const SeqArray* seq2);
void array_to_seq(const SeqArray* array, char** seq_out, char** qual_out);
// create and free objects
SeqArray create_seq_array(int length);
void free_seq_array(SeqArray* arr);
void free_barcode_group(BarcodeGroup* group);
void cleanup_on_error(FILE* fwd_file, OutputFiles* output, pthread_t* threads, WorkerArgs* thread_args);
void fatal_alloc(const char* where);
// other
double log10addexp_vec(const double *vals, int n);
void print_usage(const char* program_name);
static int cmp_idxsim_desc(const void* a, const void* b);

// Main
int main(int argc, char *argv[]) {
    // Default values
    const char *in1_file = NULL;
    const char *in2_file = NULL;
    const char *out1_file = "seq_merge_out1.fastq";
    const char *out2_file = "seq_merge_out2.fastq";
    int bc_start = 0;
    int bc_len = 18;
    int num_threads = 1;
    bool no_alignment = false;
    pthread_t* threads = NULL;
    WorkerArgs* thread_args = NULL;
    OutputFiles output = {NULL, NULL, PTHREAD_MUTEX_INITIALIZER};

    WorkQueue work_queue;
    init_work_queue(&work_queue);

    // Command line options
    static struct option long_options[] = {
        {"in", required_argument, 0, 'a'},
        {"bc-start", required_argument, 0, 's'},
        {"bc-len", required_argument, 0, 'l'},
        {"max-len", required_argument, 0, 'm'},
        {"gap", required_argument, 0, 'g'},
        {"out", required_argument, 0, '1'},
        {"threads", required_argument, 0, 't'},
        {"no-alignment", no_argument, 0, 'n'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;
    while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
        switch (c) {
            case 'a': in1_file = optarg; break;
            case 's': bc_start = atoi(optarg); break;
            case 'l': bc_len = atoi(optarg); break;
            case 'm': max_line_len = atoi(optarg); break;
            case 'g': gap_pen = atof(optarg); break;
            case '1': out1_file = optarg; break;
            case '2': out2_file = optarg; break;
            case 't': num_threads = atoi(optarg); break;
            case 'n': no_alignment = true; break;
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
        // Check for end of file and read headers
        if (!fgets(line_fwd, max_line_len, fwd_file)) {
            if (current_group) {
                queue_push(&work_queue, current_group);
            }
            break;
        }

        // Read sequences
        (void)fgets(line_fwd, max_line_len, fwd_file);
        char* fwd_seq = strndup(line_fwd, strcspn(line_fwd, "\n"));

        // Skip '+' lines
        (void)fgets(line_fwd, max_line_len, fwd_file);

        // Read qualities
        (void)fgets(line_fwd, max_line_len, fwd_file);
        char* fwd_qual = strndup(line_fwd, strcspn(line_fwd, "\n"));

        // Extract barcode
        char bc[MAX_BARCODE_LEN];
        strncpy(bc, fwd_seq + bc_start, bc_len);
        bc[bc_len] = '\0';

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
        array_to_seq(&fwd_consensus, &fwd_seq, &fwd_qual);

        // Format FASTQ entry
        int fwd_len = snprintf(NULL, 0, "@fwd_read;bc=%s;count=%d\n%s\n+\n%s\n", 
                              group->bc, group->count, fwd_seq, fwd_qual) + 1;
        
        *fwd_entry = malloc(fwd_len);
        snprintf(*fwd_entry, fwd_len, "@fwd_read;bc=%s;count=%d\n%s\n+\n%s\n", 
                group->bc, group->count, fwd_seq, fwd_qual);

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
    array_to_seq(&merged_fwd, &fwd_seq, &fwd_qual);

    // Format FASTQ entry
    int fwd_len = snprintf(NULL, 0, "@fwd_read;bc=%s;count=%d\n%s\n+\n%s\n", 
                          group->bc, group->count, fwd_seq, fwd_qual) + 1;
    
    *fwd_entry = malloc(fwd_len);
    snprintf(*fwd_entry, fwd_len, "@fwd_read;bc=%s;count=%d\n%s\n+\n%s\n", 
            group->bc, group->count, fwd_seq, fwd_qual);

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

// Merge all reads associated with one barcode
SeqArray merge_seqs(const SeqArray* seq1, const SeqArray* seq2) {
    // seq1 = reference = consensus
    // seq2 = query = individual read

    int trace_len = 0;
    int *trace = NULL;

    int len1 = seq1->length;
    int len2 = seq2->length;
    int longer = (len1 > len2) ? len1 : len2;

    // Use full alignment for very short sequences
    if (longer <= 50) {
        trace = align_arrays(seq1, seq2, &trace_len);
    } else {
        // Try banded alignments, doubling max_phase_diff until success or limit reached
        int max_phase_diff = (int)floor(sqrt((double)longer * 0.04));
        if (max_phase_diff < 1) max_phase_diff = 1;
        while (true) {
            trace = align_arrays_band(seq1, seq2, &trace_len, max_phase_diff);
            if (trace != NULL && trace_len > 0) {
                // banded alignment succeeded and was guaranteed optimal
                break;
            }
            // If band returned NULL but trace_len == 0 (no band cells) or trace == NULL,
            // we'll try increasing band unless we've exceeded allowed limit.
            if (max_phase_diff > (longer / 2)) {
                // Give up on banded approach
                break;
            }
            max_phase_diff *= 2;
            if (max_phase_diff > (longer / 2)) {
                max_phase_diff = (longer / 2) + 1; // ensure loop termination
            }
        }

        // If banded alignment did not produce a guaranteed-optimal traceback, fall back
        if (trace == NULL || trace_len <= 0) {
            trace = align_arrays(seq1, seq2, &trace_len);
        }
    }

    SeqArray merged = create_seq_array(trace_len);
    int pos1 = 0, pos2 = 0;

    for(int i=0; i<trace_len; i++) {
        Position p1 = {0}, p2 = {0};

        if (pos1 < seq1->length) {
            memcpy(p1.scores, seq1->positions[pos1].scores, sizeof(p1.scores));
        }
        if (pos2 < seq2->length) {
            memcpy(p2.scores, seq2->positions[pos2].scores, sizeof(p2.scores));
        }

        switch(trace[i]) {
            case 1:  // match
                for(int j=0; j<5; j++) {
                    merged.positions[i].scores[j] = p1.scores[j] + p2.scores[j];
                }
                pos1++;
                pos2++;
                break;
            case 2:  // gap in query
                if (pos1 < seq1->length) {
                    memcpy(merged.positions[i].scores, p1.scores, 5*sizeof(int));
                    for(int j=0; j<5; j++) {
                        merged.positions[i].scores[4] += p2.scores[j];
                    }
                }
                pos1++;
                break;
            case 3:  // gap in reference
                if (pos2 < seq2->length) {
                    memcpy(merged.positions[i].scores, p2.scores, 5*sizeof(int));
                    for(int j=0; j<5; j++) {
                        merged.positions[i].scores[4] += p1.scores[j];
                    }
                }
                pos2++;
                break;
            default:
                // any unexpected code: treat as mismatch (advance both)
                for(int j=0; j<5; j++) {
                    merged.positions[i].scores[j] = p1.scores[j] + p2.scores[j];
                }
                pos1++; pos2++;
                break;
        }
    }

    free(trace);
    return merged;
}

// Global alignment within a band
int* align_arrays_band(const SeqArray* ref, const SeqArray* query, int *out_len, int max_phase_diff) {

    const float NEG_INF = -1e30f;
    int m = ref->length;    // x dimension (reference)
    int n = query->length;  // y dimension (query)
    int ylen = n + 1;
    int w = max_phase_diff;

    if (out_len) *out_len = 0;

    // Allocate row band descriptors
    int *row_xmin = malloc(ylen * sizeof(int));
    int *row_width = malloc(ylen * sizeof(int));
    if (!row_xmin || !row_width) {
        free(row_xmin); free(row_width);
        fatal_alloc("Failed to allocate band index arrays");
    }

    // First pass: compute per-row xmin/width and total cells
    size_t total_cells = 0;
    for (int y = 0; y < ylen; y++) {
        int xmin = y - w;
        if (xmin < 0) xmin = 0;
        int xmax = y + w;
        if (xmax > m) xmax = m;
        int width = (xmax >= xmin) ? (xmax - xmin + 1) : 0;
        row_xmin[y] = xmin;
        row_width[y] = width;
        total_cells += (size_t)width;
    }

    if (total_cells == 0) {
        free(row_xmin); free(row_width);
        if (out_len) *out_len = 0;
        return NULL;
    }

    // Allocate contiguous blocks for score and trace
    float *score_block = malloc(total_cells * sizeof(float));
    int *trace_block = malloc(total_cells * sizeof(int));
    if (!score_block || !trace_block) {
        free(score_block); free(trace_block);
        free(row_xmin); free(row_width);
        fatal_alloc("Failed to allocate banded matrix blocks");
    }

    // Row pointer arrays into contiguous blocks
    float **score = malloc(ylen * sizeof(float*));
    int **trace = malloc(ylen * sizeof(int*));
    if (!score || !trace) {
        free(score_block); free(trace_block);
        free(score); free(trace);
        free(row_xmin); free(row_width);
        fatal_alloc("Failed to allocate band row pointers");
    }

    // Set row pointers and initialize cells to NEG_INF/0
    size_t offset = 0;
    for (int y = 0; y < ylen; y++) {
        int wrow = row_width[y];
        if (wrow > 0) {
            score[y] = score_block + offset;
            trace[y] = trace_block + offset;
            for (int i = 0; i < wrow; i++) {
                score[y][i] = NEG_INF;
                trace[y][i] = 0;
            }
            offset += (size_t)wrow;
        } else {
            score[y] = NULL;
            trace[y] = NULL;
        }
    }

    #define IN_BAND(yy, xx) ((xx) >= row_xmin[(yy)] && (xx) < row_xmin[(yy)] + row_width[(yy)])
    #define IDX(yy, xx) ((xx) - row_xmin[(yy)])

    // Linear gap penalty (use gap_open as linear gap)
    float g = (float)gap_pen;

    // Initialize (0,0) if in band
    if (IN_BAND(0,0)) {
        int i0 = IDX(0,0);
        score[0][i0] = 0.0f;
        trace[0][i0] = 0;
    }

    // Initialize first row (y=0): propagate gaps to the right within band
    if (row_width[0] > 0) {
        for (int x = row_xmin[0]; x < row_xmin[0] + row_width[0]; x++) {
            if (x == 0) continue;
            if (!IN_BAND(0, x-1)) break; // can't propagate if left cell not in band
            int i = IDX(0, x);
            int left_i = IDX(0, x-1);
            score[0][i] = score[0][left_i] + g;
            trace[0][i] = 2; // gap in query (deletion)
        }
    }

    // Initialize first column (x=0): propagate gaps downward within band
    for (int y = 1; y < ylen; y++) {
        if (!IN_BAND(y, 0) || !IN_BAND(y-1, 0)) continue;
        int i = IDX(y, 0);
        int up_i = IDX(y-1, 0);
        score[y][i] = score[y-1][up_i] + g;
        trace[y][i] = 3; // gap in reference (insertion)
    }

    // Fill DP within band computing match on-the-fly
    for (int y = 1; y < ylen; y++) {
        if (row_width[y] == 0) continue;
        for (int x = row_xmin[y]; x < row_xmin[y] + row_width[y]; x++) {
            int xi = IDX(y, x);
            float best = NEG_INF;
            int t = 0;

            // Diagonal (match/mismatch)
            if (x > 0 && y > 0 && IN_BAND(y-1, x-1)) {
                int prev_i = IDX(y-1, x-1);
                float prev_score = score[y-1][prev_i];
                if (prev_score > NEG_INF/2.0f) {
                    double mscore_d = compare_positions(&ref->positions[x-1], &query->positions[y-1]);
                    float mscore = (float)mscore_d;
                    float diag = prev_score + mscore;
                    best = diag;
                    t = 1;
                }
            }

            // Left (deletion)
            if (x > 0 && IN_BAND(y, x-1)) {
                int left_i = IDX(y, x-1);
                float left_score = score[y][left_i];
                if (left_score > NEG_INF/2.0f) {
                    // if consensus is gapped, relax the penalty against putting a gap in the query
                    float g_left = g;
                    if (ref->positions[x-1].scores[4] > 0) {
                        int total_ref = 0;
                        for (int bb = 0; bb < 5; bb++) total_ref += ref->positions[x-1].scores[bb];
                        // total_ref > 0 is guaranteed because gap count > 0 and scores are non-negative
                        float gap_frac_ref = (float)ref->positions[x-1].scores[4] / (float)total_ref;
                        g_left = g * (1.0f - gap_frac_ref);
                    }
                    float cand = left_score + g_left;
                    if (cand > best) { best = cand; t = 2; }
                }
            }

            // Up (insertion)
            if (y > 0 && IN_BAND(y-1, x)) {
                int up_i = IDX(y-1, x);
                float up_score = score[y-1][up_i];
                if (up_score > NEG_INF/2.0f) {
                    float cand = up_score + g;
                    if (cand > best) { best = cand; t = 3; }
                }
            }

            score[y][xi] = best;
            trace[y][xi] = t;
        }
    }

    // Ensure bottom-right (n,m) is in band and get sband from that cell
    if (!IN_BAND(n, m)) {
        // cleanup
        free(score_block); free(trace_block);
        free(score); free(trace);
        free(row_xmin); free(row_width);
        if (out_len) *out_len = 0;
        return NULL;
    }
    float sband = score[n][IDX(n, m)];
    if (sband <= NEG_INF/2.0f) {
        // unreachable end cell
        free(score_block); free(trace_block);
        free(score); free(trace);
        free(row_xmin); free(row_width);
        if (out_len) *out_len = 0;
        return NULL;
    }

    // Compute bestout bound: enforce l1 >= l2 and use m=1.0
    double l1 = (double)((m > n) ? m : n);
    double l2 = (double)((m > n) ? n : m);
    double mm = 1.0;
    double bestout = (2.0 * (w + 1.0) - (l1 - l2)) * (double)gap_pen + (l2 - (w + 1.0)) * mm;

    if (!((double)sband >= bestout)) {
        // band did not prove optimal
        free(score_block); free(trace_block);
        free(score); free(trace);
        free(row_xmin); free(row_width);
        if (out_len) *out_len = -1;
        return NULL;
    }

    // Traceback starting at (m,n)
    int tx = m;
    int ty = n;
    int max_trace_len = m + n + 1;
    int *traceback = malloc(max_trace_len * sizeof(int));
    if (!traceback) {
        free(score_block); free(trace_block);
        free(score); free(trace);
        free(row_xmin); free(row_width);
        return NULL;
    }
    int idx_tb = 0;

    while (tx > 0 || ty > 0) {
        if (!IN_BAND(ty, tx)) {
            // If we're outside band, consume gaps to move back into band
            if (tx > row_xmin[ty] + row_width[ty] - 1) {
                traceback[idx_tb++] = 2; tx--; continue;
            } else if (tx < row_xmin[ty]) {
                traceback[idx_tb++] = 3; ty--; continue;
            } else {
                break;
            }
        }
        int ci = IDX(ty, tx);
        int t = trace[ty][ci];
        if (t == 0) break;
        traceback[idx_tb++] = t;
        if (t == 1) { tx--; ty--; }
        else if (t == 2) { tx--; }
        else { ty--; }
    }

    // reverse traceback
    for (int i = 0; i < idx_tb/2; i++) {
        int tmp = traceback[i];
        traceback[i] = traceback[idx_tb - i - 1];
        traceback[idx_tb - i - 1] = tmp;
    }
    if (out_len) *out_len = idx_tb;

    // cleanup banded storage
    free(score_block);
    free(trace_block);
    free(score);
    free(trace);
    free(row_xmin);
    free(row_width);

    return traceback;
}

// True global alignment
int* align_arrays(const SeqArray* ref, const SeqArray* query, int *out_len) {

    const float NEG_INF = -1e30f;
    int ylen = query->length + 1; // rows (y)
    int xlen = ref->length + 1;   // cols (x)

    if (out_len) *out_len = 0;

    if (ylen <= 0 || xlen <= 0) return NULL;

    size_t cells = (size_t)ylen * (size_t)xlen;

    // Allocate contiguous blocks
    float *score_block = malloc(cells * sizeof(float));
    int *trace_block = malloc(cells * sizeof(int));
    if (!score_block || !trace_block) {
        free(score_block);
        free(trace_block);
        return NULL;
    }

    // convenience macro for indexing
    #define IDX(y,x) ((size_t)(y) * (size_t)xlen + (size_t)(x))

    // Initialize score and trace
    for (size_t i = 0; i < cells; i++) {
        score_block[i] = NEG_INF;
        trace_block[i] = 0;
    }

    float g = (float)gap_pen; // linear gap penalty

    // (0,0)
    score_block[IDX(0,0)] = 0.0f;
    trace_block[IDX(0,0)] = 0;

    // First row (gaps in query -> deletions)
    for (int x = 1; x < xlen; x++) {
        score_block[IDX(0,x)] = score_block[IDX(0,x-1)] + g;
        trace_block[IDX(0,x)] = 2; // deletion (gap in query)
    }

    // First column (gaps in reference -> insertions)
    for (int y = 1; y < ylen; y++) {
        score_block[IDX(y,0)] = score_block[IDX(y-1,0)] + g;
        trace_block[IDX(y,0)] = 3; // insertion (gap in reference)
    }

    // Fill DP; compute match on-the-fly
    for (int y = 1; y < ylen; y++) {
        for (int x = 1; x < xlen; x++) {
            float best = NEG_INF;
            int t = 0;

            // Diagonal (match/mismatch)
            float diag = NEG_INF;
            double mscore_d = compare_positions(&ref->positions[x-1], &query->positions[y-1]);
            diag = score_block[IDX(y-1,x-1)] + (float)mscore_d;
            best = diag;
            t = 1;

            // Left (deletion: gap in query)
            float g_left = g;
            // if consensus is gapped, relax the penalty against putting a gap in the query
            if (ref->positions[x-1].scores[4] > 0) {
                int total_ref = 0;
                for (int bb = 0; bb < 5; bb++) total_ref += ref->positions[x-1].scores[bb];
                // total_ref > 0 is guaranteed because gap count > 0 and scores are non-negative
                float gap_frac_ref = (float)ref->positions[x-1].scores[4] / (float)total_ref;
                g_left = g * (1.0f - gap_frac_ref);
            }
            float left = score_block[IDX(y,x-1)] + g_left;
            if (left > best) { best = left; t = 2; }

            // Up (insertion: gap in reference)
            float up = score_block[IDX(y-1,x)] + g;
            if (up > best) { best = up; t = 3; }

            score_block[IDX(y,x)] = best;
            trace_block[IDX(y,x)] = t;
        }
    }

    // Traceback (same semantics as previous implementation)
    int max_trace_len = query->length + ref->length + 1;
    int *traceback = malloc(max_trace_len * sizeof(int));
    if (!traceback) {
        free(score_block);
        free(trace_block);
        return NULL;
    }
    int idx = 0;
    int x = xlen - 1;
    int y = ylen - 1;

    while (x > 0 && y > 0) {
        int tv = trace_block[IDX(y,x)];
        traceback[idx++] = tv;
        if (tv == 1) { x--; y--; }
        else if (tv == 2) { x--; }
        else { y--; }
    }
    while (x > 0) { traceback[idx++] = 2; x--; }
    while (y > 0) { traceback[idx++] = 3; y--; }

    // reverse traceback
    for (int i = 0; i < idx/2; i++) {
        int tmp = traceback[i];
        traceback[i] = traceback[idx - i - 1];
        traceback[idx - i - 1] = tmp;
    }
    if (out_len) *out_len = idx;

    // cleanup
    free(score_block);
    free(trace_block);

    #undef IDX

    return traceback;
}

// Convert a sequence and quality string into a 5xL array
SeqArray seq_to_array(const char* seq, const char* qual, int length) {
    // Initialize array
    SeqArray array = create_seq_array(length);

    size_t seqlen = strlen(seq);
    size_t quallen = qual ? strlen(qual) : 0;
    size_t n = seqlen < (size_t)length ? seqlen : (size_t)length;
    
    // Fill the array
    for (size_t i = 0; i < n; i++) {
        // Check valid base
        char base = seq[i];
        if (base == 'A' || base == 'C' || base == 'G' || base == 'T') {
            // Convert to adjusted quality score with bounds (guard qual length)
            int adjusted_qual_score = 0;
            if (i >= quallen) {
                adjusted_qual_score = 0;
            } else if (qual[i] < 33) {
                adjusted_qual_score = 0;
            } else if (qual[i] > 73) {
                adjusted_qual_score = 42;
            } else {
                adjusted_qual_score = q_adj[qual[i] - 33];
            }
            // Assign quality score to the corresponding base
            array.positions[i].scores[basemap[(int)base]] = adjusted_qual_score;
        }
    }
    return array;
}

// Convert consensus array back to a sequence + quality
void array_to_seq(const SeqArray* array, char** seq_out, char** qual_out) {
    int out_pos = 0;
    *seq_out = malloc(array->length + 1);
    *qual_out = malloc(array->length + 1);
    
    // Calculate sum of all cells in the array and average evidence per position
    long total_evidence = 0;
    int total_positions = array->length;
    
    for(int i = 0; i < array->length; i++) {
        const Position* pos = &array->positions[i];
        for(int j = 0; j < 5; j++) {
            total_evidence += pos->scores[j];
        }
    }
    
    double avg_evidence_per_pos = (double)total_evidence / total_positions;
    
    // Find last position with >50% of average evidence
    int last_good_pos = array->length - 1;
    double threshold = 0.5 * avg_evidence_per_pos;
    
    for(; last_good_pos >= 0; last_good_pos--) {
        int pos_evidence = 0;
        for(int j = 0; j < 5; j++) {
            pos_evidence += array->positions[last_good_pos].scores[j];
        }
        
        if(pos_evidence > threshold) {
            break;
        }
    }
    
    for(int i = 0; i <= last_good_pos; i++) {
        const Position* pos = &array->positions[i];

        // Check for weird numbers
        for(int j = 0; j < 5; j++) {
            if (pos->scores[j] < -1000 || pos->scores[j] > 10000000) {
                fprintf(stderr, "Warning: Potentially overflowed score detected at position %d, base %d: %d\n", 
                        i, j, pos->scores[j]);
            }
        }
        
        // 1. Find global maximum (including gap)
        int max_idx = 0;
        int max_val = pos->scores[0];
        int total_non_gap = 0; // sum of A,C,G,T evidence
        for(int j = 0; j < 5; j++) {
            if (pos->scores[j] > max_val) {
                max_val = pos->scores[j];
                max_idx = j;
            }
            if (j < 4) total_non_gap += pos->scores[j];
        }

        // If gap is most supported, skip this position
        if(max_idx == 4) continue;

        // If there is no non-gap evidence, produce N and Q=0
        if (total_non_gap == 0) {
            (*seq_out)[out_pos] = 'N';
            (*qual_out)[out_pos] = '!'; // Q=0
            out_pos++;
            continue;
        }

        // Tie handling: if two or more non-gap bases share the same top evidence, declare N with Q=0.
        //   - We only consider A,C,G,T for tie detection (not gap, which was already handled above).
        int top_count = 0;
        for (int j = 0; j < 4; j++) {
            if (pos->scores[j] == max_val) top_count++;
        }
        if (top_count > 1) {
            // ambiguous top between two or more bases -> emit N, Q=0
            (*seq_out)[out_pos] = 'N';
            (*qual_out)[out_pos] = '!';
            out_pos++;
            continue;
        }

        // 2. Compute per-class log10 probabilities (one value for each of A,C,G,T,gap).
        //     q_correct = (1/10) * log10(0.9)
        //     q_error   = (1/10) * log10(0.1)
        double log10p[5];
        int total_all = 0;
        for (int j = 0; j < 5; j++) total_all += pos->scores[j]; // calculate total scores
        for (int j = 0; j < 5; j++) {
            // for each base B, log10(B) = Q_B * q_correct + Q_!B * q_error
            log10p[j] = pos->scores[j] * q_correct + (double)(total_all - pos->scores[j]) * q_error;
        }

        // 3. Compute log10 of probability mass for (a) all bases and (b) the nonconsensus bases
        //      denom = log10( sum_{B in bases} 10^{log10p[B]} )
        //      other = log10( sum_{B != max_idx in bases} 10^{log10p[B]} )
        //    Then the Phred-like quality is:
        //      Q = -10 * other + 10 * denom
        double denom = log10addexp_vec(log10p, 5);

        // Build a small array of size 4 for the "other" classes (all except max_idx)
        double other_vals[5];
        int oth_count = 0;
        for (int j = 0; j < 5; j++) {
            if (j == max_idx) continue;
            other_vals[oth_count++] = log10p[j];
        }
        double other = log10addexp_vec(other_vals, oth_count);

        // 4. Convert to a Phred-like quality score
        double qual_d = -10.0 * other + 10.0 * denom;

        // Clamp quality values between 0 and 40
        int qual_int = (int)round(qual_d);
        if (qual_int < 0) qual_int = 0;
        if (qual_int > 40) qual_int = 40;

        // For very low Q, declare N
        if (qual_int < 4) {
            (*seq_out)[out_pos] = 'N';
        } else {
            (*seq_out)[out_pos] = BASES[max_idx];
        }
        (*qual_out)[out_pos] = (char)(qual_int + 33);

        // Update global error counters
        for (int j = 0; j < 4; j++) {
            call_total += pos->scores[j];
            if (j != max_idx) {
                call_missense += pos->scores[j];
            }
        }
        if (pos->scores[4] > 0) {
            call_indel += pos->scores[4];
        }

        out_pos++;
    }
    
    // Terminate strings
    (*seq_out)[out_pos] = '\0';
    (*qual_out)[out_pos] = '\0';
    
    // Reallocate to trimmed size
    *seq_out = realloc(*seq_out, out_pos + 1);
    *qual_out = realloc(*qual_out, out_pos + 1);
}

// Compare two positions using scaled cosine similarity
double compare_positions(const Position* vec1, const Position* vec2) {
    double dot = 0.0, mag1 = 0.0, mag2 = 0.0;
    
    for (int i = 0; i < VECTOR_LENGTH; i++) {
        double v1 = vec1->scores[i];
        double v2 = vec2->scores[i];
        dot += v1 * v2;
        mag1 += v1 * v1;
        mag2 += v2 * v2;
    }
    
    // heuristics to speed up return
    if (mag1 == 0.0 || mag2 == 0.0) return 0.0; // empty vectors

    // calculate cosine similarity
    return (2.0 * (dot / (sqrt(mag1) * sqrt(mag2)))) - 1.0;
}

// Compare two sequences by summed position similarity
double compare_seqs(const SeqArray* a, const SeqArray* b) {
    if(a->length != b->length || a->length == 0) return 0.0;
    
    double total = 0.0;
    for(int i=0; i<a->length; i++) {
        total += compare_positions(&a->positions[i], &b->positions[i]);
    }
    return total;
}

// Build consensus without alignment by summing position scores across sequences
SeqArray build_unaligned_consensus(SeqArray* sequences, int count) {
    if(count == 0) return create_seq_array(0);
    
    int length = sequences[0].length;
    SeqArray consensus = create_seq_array(length);
    
    for(int i=0; i<count; i++) {
        for(int pos=0; pos<length; pos++) {
            for(int base=0; base<5; base++) {
                consensus.positions[pos].scores[base] += 
                    sequences[i].positions[pos].scores[base];
            }
        }
    }
    return consensus;
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
// Helper used locally to sort reads by similarity
static int cmp_idxsim_desc(const void* a, const void* b) {
    const IdxSim* A = a;
    const IdxSim* B = b;
    if (A->sim < B->sim) return 1;
    if (A->sim > B->sim) return -1;
    return 0;
}

// Helper function to create a new SeqArray
SeqArray create_seq_array(int length) {
    SeqArray arr;
    arr.length = length;
    arr.positions = calloc(length, sizeof(Position));
    if (!arr.positions) {
        fatal_alloc("Failed to allocate SeqArray");
        arr.positions = NULL;
        arr.length = 0;
    }
    return arr;
}

// Helper function to print message and perform global cleanup (will exit)
void fatal_alloc(const char* where) {
    fprintf(stderr, "Fatal error: %s\n", where);
    cleanup_on_error(g_fwd_file, g_output_ptr, g_threads_ptr, g_thread_args_ptr);
}

// Helper function to avoid underflow on q-scores
double log10addexp_vec(const double *vals, int n) {
    // Handle empty input
    if (n <= 0 || vals == NULL) {
        return -INFINITY; // log10(0)
    }

    // Find maximum value to factor out for numerical stability:
    // log10(sum_i 10^{v_i}) = maxV + log10( sum_i 10^{v_i - maxV} )
    double maxV = vals[0];
    for (int i = 1; i < n; ++i) {
        if (vals[i] > maxV) maxV = vals[i];
    }

    // If maxV is -inf, all terms are zero
    if (!isfinite(maxV)) return -INFINITY;

    // Accumulate scaled sums in normal space, but values are 10^{v - maxV}
    double acc = 0.0;
    for (int i = 0; i < n; ++i) {
        double diff = vals[i] - maxV;
        // If diff is very small, pow(10, diff) will underflow to 0; that's fine.
        acc += pow(10.0, diff);
    }

    // Return log10 of the sum
    return maxV + log10(acc);
}

// Cleanup definitions
void free_barcode_group(BarcodeGroup* group) {
    for (int i = 0; i < group->count; i++) {
        free(group->fwd_reads[i]);
        free(group->fwd_quals[i]);
    }
    free(group->fwd_reads);
    free(group->fwd_quals);
    free(group);
}

void free_seq_array(SeqArray* arr) {
    free(arr->positions);
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
    fprintf(stderr, "  --in file          Input FASTQ file 1 (required)\n");
    fprintf(stderr, "  --bc-start int     Barcode start position (Zero-indexed, default: 0)\n");
    fprintf(stderr, "  --bc-len int       Barcode length (default: 18)\n");
    fprintf(stderr, "  --max-len int      Maximum read length (default: 131072)\n");
    fprintf(stderr, "  --gap float        Gap penalty for alignment (default: -1.0)\n");
    fprintf(stderr, "  --out file         Output file for consensus read 1\n");
    fprintf(stderr, "  --threads int      Number of threads (default: 1)\n");
    fprintf(stderr, "  --no-alignment     Skip alignment/merging; use unaligned consensus directly\n");
    exit(EXIT_FAILURE);
}