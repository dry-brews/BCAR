#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <pthread.h>
#include <getopt.h>

#define MAX_BARCODE_LEN 128
#define MAX_LINE_LEN 1024
#define BASES "ACGT-"
#define MAX_CONS_LEN 1024 * 4
#define DEFAULT_GAP_SCORE -3.0
#define VECTOR_LENGTH 5
#define q_correct -0.004575749 // = 1/10 * log10(0.9)
#define q_error -0.1 // = 1/10 * log10(0.1)

// Structure for barcode group (all reads associated with one barcode)
typedef struct {
    char bc[MAX_BARCODE_LEN];
    char** fwd_reads;
    char** fwd_quals;
    char** rev_reads;
    char** rev_quals;
    int count;
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
    FILE* rev_out;
    pthread_mutex_t mutex;
} OutputFiles;

// Structure to pass to worker threads
typedef struct {
    WorkQueue* queue;
    OutputFiles* output;
    int thread_id;
} WorkerArgs;

// Global variables
const char* bases = BASES;
const int basemap[256] = { ['A']=0, ['C']=1, ['G']=2, ['T']=3, ['-']=4, ['N']=4 };
const int q_adj[41] = { // adjusted quality score mapping
    0, -6, -2, 0, 2, 4, 6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 18, 20, 21, 22, 23, 24, 25, 26, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 41, 43, 44, 45, 46
};
double gap_open = DEFAULT_GAP_SCORE;
double gap_extend = DEFAULT_GAP_SCORE;

// Function prototypes
// organize workers into a queue
void init_work_queue(WorkQueue* queue);
void queue_push(WorkQueue* queue, BarcodeGroup* group);
BarcodeGroup* queue_pop(WorkQueue* queue);
void* worker_thread(void* arg);
// process barcode group
void process_barcode(BarcodeGroup* group, char** fwd_entry, char** rev_entry);
SeqArray seq_to_array(const char* seq, const char* qual, int length);
double compare_positions(const Position* a, const Position* b);
double compare_seqs(const SeqArray* a, const SeqArray* b);
SeqArray build_unaligned_consensus(SeqArray* sequences, int count);
char* align_arrays(const SeqArray* query, const SeqArray* ref);
SeqArray merge_seqs(const SeqArray* seq1, const SeqArray* seq2);
void array_to_seq(const SeqArray* array, char** seq_out, char** qual_out);
// create and free objects
SeqArray create_seq_array(int length);
Matrix create_matrix(int rows, int cols);
void free_matrix(Matrix mat);
void free_seq_array(SeqArray* arr);
void free_barcode_group(BarcodeGroup* group);
// other helper functions
double log10addexp(double A, double B);
void print_usage(const char* program_name);

int main(int argc, char *argv[]) {
    // Default values
    const char *in1_file = NULL;
    const char *in2_file = NULL;
    const char *out1_file = "seq_merge_out1.fastq";
    const char *out2_file = "seq_merge_out2.fastq";
    int bc_start = 0;
    int bc_len = 18;
    int num_threads = 1;

    // Command line options
    static struct option long_options[] = {
        {"read1", required_argument, 0, 'a'},
        {"read2", required_argument, 0, 'b'},
        {"bc-start", required_argument, 0, 's'},
        {"bc-len", required_argument, 0, 'l'},
        {"gap-open", required_argument, 0, 'g'},
        {"gap-extend", required_argument, 0, 'e'},
        {"out1", required_argument, 0, '1'},
        {"out2", required_argument, 0, '2'},
        {"threads", required_argument, 0, 't'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;
    while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
        switch (c) {
            case 'a': in1_file = optarg; break;
            case 'b': in2_file = optarg; break;
            case 's': bc_start = atoi(optarg); break;
            case 'l': bc_len = atoi(optarg); break;
            case 'g': gap_open = atof(optarg); break;
            case 'e': gap_extend = atof(optarg); break;
            case '1': out1_file = optarg; break;
            case '2': out2_file = optarg; break;
            case 't': num_threads = atoi(optarg); break;
            default: print_usage(argv[0]);
        }
    }

    // Validate arguments
    if (!in1_file || !in2_file) {
        fprintf(stderr, "Error: Missing required arguments\n");
        print_usage(argv[0]);
    }

    if (num_threads < 1) {
        fprintf(stderr, "Error: Invalid thread count\n");
        exit(EXIT_FAILURE);
    }

    // Initialize input files
    FILE* fwd_file = fopen(in1_file, "r");
    FILE* rev_file = fopen(in2_file, "r");
    if (!fwd_file || !rev_file) {
        perror("Error opening input files");
        exit(EXIT_FAILURE);
    }

    // Initialize output files and mutex
    OutputFiles output;
    output.fwd_out = fopen(out1_file, "w");
    output.rev_out = fopen(out2_file, "w");
    if (!output.fwd_out || !output.rev_out) {
        perror("Error opening output files");
        exit(EXIT_FAILURE);
    }
    pthread_mutex_init(&output.mutex, NULL);

    // Initialize work queue
    WorkQueue work_queue;
    init_work_queue(&work_queue);

    // Create worker threads
    pthread_t* threads = malloc(num_threads * sizeof(pthread_t));
    WorkerArgs* thread_args = malloc(num_threads * sizeof(WorkerArgs));
    
    for (int i = 0; i < num_threads; i++) {
        thread_args[i].queue = &work_queue;
        thread_args[i].output = &output;
        thread_args[i].thread_id = i;
        pthread_create(&threads[i], NULL, worker_thread, &thread_args[i]);
    }

    // Producer: Read input files and create work
    BarcodeGroup* current_group = NULL;
    char line_fwd[MAX_LINE_LEN];
    char line_rev[MAX_LINE_LEN];

    while (1) {
        // Check for end of file and read headers
        if (!fgets(line_fwd, MAX_LINE_LEN, fwd_file) || 
            !fgets(line_rev, MAX_LINE_LEN, rev_file)) {
            if (current_group) {
                queue_push(&work_queue, current_group);
            }
            break;
        }

        // Read sequences
        (void)fgets(line_fwd, MAX_LINE_LEN, fwd_file);
        (void)fgets(line_rev, MAX_LINE_LEN, rev_file);
        char* fwd_seq = strndup(line_fwd, strcspn(line_fwd, "\n"));
        char* rev_seq = strndup(line_rev, strcspn(line_rev, "\n"));

        // Skip '+' lines
        (void)fgets(line_fwd, MAX_LINE_LEN, fwd_file);
        (void)fgets(line_rev, MAX_LINE_LEN, rev_file);

        // Read qualities
        fgets(line_fwd, MAX_LINE_LEN, fwd_file);
        fgets(line_rev, MAX_LINE_LEN, rev_file);
        char* fwd_qual = strndup(line_fwd, strcspn(line_fwd, "\n"));
        char* rev_qual = strndup(line_rev, strcspn(line_rev, "\n"));

        // Extract barcode
        char bc[MAX_BARCODE_LEN];
        strncpy(bc, fwd_seq + bc_start, bc_len);
        bc[bc_len] = '\0';

        // Handle group management
        if (!current_group) {
            current_group = malloc(sizeof(BarcodeGroup));
            strcpy(current_group->bc, bc);
            current_group->count = 0;
            current_group->fwd_reads = malloc(sizeof(char*));
            current_group->rev_reads = malloc(sizeof(char*));
            current_group->fwd_quals = malloc(sizeof(char*));
            current_group->rev_quals = malloc(sizeof(char*));
        }

        if (strcmp(current_group->bc, bc) != 0) {
            queue_push(&work_queue, current_group);
            
            current_group = malloc(sizeof(BarcodeGroup));
            strcpy(current_group->bc, bc);
            current_group->count = 0;
            current_group->fwd_reads = malloc(sizeof(char*));
            current_group->rev_reads = malloc(sizeof(char*));
            current_group->fwd_quals = malloc(sizeof(char*));
            current_group->rev_quals = malloc(sizeof(char*));
        }

        // Add to group
        current_group->fwd_reads = realloc(current_group->fwd_reads, 
            sizeof(char*) * (current_group->count + 1));
        current_group->rev_reads = realloc(current_group->rev_reads, 
            sizeof(char*) * (current_group->count + 1));
        current_group->fwd_quals = realloc(current_group->fwd_quals, 
            sizeof(char*) * (current_group->count + 1));
        current_group->rev_quals = realloc(current_group->rev_quals, 
            sizeof(char*) * (current_group->count + 1));

        current_group->fwd_reads[current_group->count] = fwd_seq;
        current_group->rev_reads[current_group->count] = rev_seq;
        current_group->fwd_quals[current_group->count] = fwd_qual;
        current_group->rev_quals[current_group->count] = rev_qual;
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
    fclose(rev_file);
    fclose(output.fwd_out);
    fclose(output.rev_out);

    return 0;
}

void process_barcode(BarcodeGroup* group, char** fwd_entry, char** rev_entry) {
    *fwd_entry = NULL;
    *rev_entry = NULL;
    // Step 0: Check group size and don't do full processing if not necessary
    if(group->count == 0) return;
    if(group->count == 1) {
        // Format FASTQ entries
        int fwd_len = snprintf(NULL, 0, "@fwd_read;bc=%s;count=%d\n%s\n+\n%s\n", 
                              group->bc, group->count, group->fwd_reads[0], group->fwd_quals[0]) + 1;
        int rev_len = snprintf(NULL, 0, "@rev_read;bc=%s;count=%d\n%s\n+\n%s\n", 
                              group->bc, group->count, group->rev_reads[0], group->rev_quals[0]) + 1;
        
        *fwd_entry = malloc(fwd_len);
        *rev_entry = malloc(rev_len);
        snprintf(*fwd_entry, fwd_len, "@fwd_read;bc=%s;count=%d\n%s\n+\n%s\n", 
                group->bc, group->count, group->fwd_reads[0], group->fwd_quals[0]);
        snprintf(*rev_entry, rev_len, "@rev_read;bc=%s;count=%d\n%s\n+\n%s\n", 
                group->bc, group->count, group->rev_reads[0], group->rev_quals[0]);
        return;
    }

    // Step 1: Get max length of reads
    int L_fwd = 0, L_rev = 0;
    for(int i = 0; i < group->count; i++) {
        int curr_fwd_len = strlen(group->fwd_reads[i]);
        int curr_rev_len = strlen(group->rev_reads[i]);
        if(curr_fwd_len > L_fwd) L_fwd = curr_fwd_len;
        if(curr_rev_len > L_rev) L_rev = curr_rev_len;
    }

    // Step 2: Convert reads to arrays
    SeqArray* fwd_arrays = malloc(group->count * sizeof(SeqArray));
    SeqArray* rev_arrays = malloc(group->count * sizeof(SeqArray));
    
    for(int i=0; i<group->count; i++) {
        fwd_arrays[i] = seq_to_array(group->fwd_reads[i], group->fwd_quals[i], L_fwd);
        rev_arrays[i] = seq_to_array(group->rev_reads[i], group->rev_quals[i], L_rev);
    }

    // Step 3: Build unaligned consensus
    SeqArray fwd_consensus = build_unaligned_consensus(fwd_arrays, group->count);
    SeqArray rev_consensus = build_unaligned_consensus(rev_arrays, group->count);

    // Step 4: Compare reads and sort
    double* similarities = malloc(group->count * sizeof(double));
    int* indices = malloc(group->count * sizeof(int));
    
    for(int i=0; i<group->count; i++) {
        double fwd_sim = compare_seqs(&fwd_arrays[i], &fwd_consensus);
        double rev_sim = compare_seqs(&rev_arrays[i], &rev_consensus);
        similarities[i] = (fwd_sim + rev_sim);
        indices[i] = i;
    }

    // Sort indices based on similarity (descending order)
    for(int i=0; i<group->count-1; i++) {
        for(int j=i+1; j<group->count; j++) {
            if(similarities[indices[i]] < similarities[indices[j]]) {
                int temp = indices[i];
                indices[i] = indices[j];
                indices[j] = temp;
            }
        }
    }

    // Step 5: Merge sequences into consensus
    // Create a copy of the first array instead of using it directly
    SeqArray merged_fwd = create_seq_array(fwd_arrays[indices[0]].length);
    SeqArray merged_rev = create_seq_array(rev_arrays[indices[0]].length);
    
    // Copy initial arrays
    memcpy(merged_fwd.positions, fwd_arrays[indices[0]].positions, 
           fwd_arrays[indices[0]].length * sizeof(Position));
    memcpy(merged_rev.positions, rev_arrays[indices[0]].positions, 
           rev_arrays[indices[0]].length * sizeof(Position));
    
    for(int i=1; i<group->count; i++) {
        SeqArray temp_fwd = merge_seqs(&merged_fwd, &fwd_arrays[indices[i]]);
        SeqArray temp_rev = merge_seqs(&merged_rev, &rev_arrays[indices[i]]);
        
        free_seq_array(&merged_fwd);
        free_seq_array(&merged_rev);
        
        merged_fwd = temp_fwd;
        merged_rev = temp_rev;
    }

    // Step 6: Convert to FASTQ and print
    char *fwd_seq, *fwd_qual, *rev_seq, *rev_qual;
    array_to_seq(&merged_fwd, &fwd_seq, &fwd_qual);
    array_to_seq(&merged_rev, &rev_seq, &rev_qual);

    // Format FASTQ entries
    int fwd_len = snprintf(NULL, 0, "@fwd_read;bc=%s;count=%d\n%s\n+\n%s\n", 
                          group->bc, group->count, fwd_seq, fwd_qual) + 1;
    int rev_len = snprintf(NULL, 0, "@rev_read;bc=%s;count=%d\n%s\n+\n%s\n", 
                          group->bc, group->count, rev_seq, rev_qual) + 1;
    
    *fwd_entry = malloc(fwd_len);
    *rev_entry = malloc(rev_len);
    snprintf(*fwd_entry, fwd_len, "@fwd_read;bc=%s;count=%d\n%s\n+\n%s\n", 
            group->bc, group->count, fwd_seq, fwd_qual);
    snprintf(*rev_entry, rev_len, "@rev_read;bc=%s;count=%d\n%s\n+\n%s\n", 
            group->bc, group->count, rev_seq, rev_qual);

    // Cleanup
    free(fwd_seq); free(fwd_qual); free(rev_seq); free(rev_qual);
    free(similarities); 
    free(indices);
    
    // Free the consensus arrays
    free_seq_array(&fwd_consensus);
    free_seq_array(&rev_consensus);
    
    // Free the merged results
    free_seq_array(&merged_fwd);
    free_seq_array(&merged_rev);
    
    // Free the individual arrays
    for(int i=0; i<group->count; i++) {
        free_seq_array(&fwd_arrays[i]);
        free_seq_array(&rev_arrays[i]);
    }
    free(fwd_arrays);
    free(rev_arrays);
}

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
    
    while (true) {
        BarcodeGroup* group = queue_pop(queue);
        if (!group) {
            break;  // Queue is empty and producer is done
        }
        
        char *fwd_entry, *rev_entry;
        process_barcode(group, &fwd_entry, &rev_entry);
        
        if (fwd_entry && rev_entry) {
            pthread_mutex_lock(&output->mutex);
            fprintf(output->fwd_out, "%s", fwd_entry);
            fprintf(output->rev_out, "%s", rev_entry);
            pthread_mutex_unlock(&output->mutex);
        }
        
        free(fwd_entry);
        free(rev_entry);
        free_barcode_group(group);
    }
    
    return NULL;
}

// Helper function to create a new SeqArray
SeqArray create_seq_array(int length) {
    SeqArray arr;
    arr.length = length;
    arr.positions = malloc(length * sizeof(Position));
    for(int i=0; i<length; i++) {
        for(int j=0; j<5; j++) {
            arr.positions[i].scores[j] = 0;
        }
    }
    return arr;
}

// Function to convert a sequence and quality string into a 5xL array
SeqArray seq_to_array(const char* seq, const char* qual, int length) {
    // Initialize array
    SeqArray array = create_seq_array(length);
    
    // Fill the array
    for (int i = 0; i < strlen(seq); i++) {
        // Check valid base
        char base = seq[i];
        if (base == 'A' || base == 'C' || base == 'G' || base == 'T') {
            // Convert to adjusted quality score with bounds
            int adjusted_qual_score = 0;
            if (qual[i] < 33) {
                adjusted_qual_score = 0;
            } else if (qual[i] > 73) {
                adjusted_qual_score = 46;
            } else {
                adjusted_qual_score = q_adj[qual[i] - 33];
            }
            // Assign quality score to the corresponding base
            array.positions[i].scores[basemap[(int)base]] = adjusted_qual_score;
        }
    }
    return array;
}

// Compare two positions using cosine similarity
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

// Compare two sequences by average position similarity
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

// Create a matrix with given dimensions
Matrix create_matrix(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        mat.data[i] = malloc(cols * sizeof(double));
    }
    return mat;
}

// Global alignment with traceback
char* align_arrays(const SeqArray* ref, const SeqArray* query) {
    int ylen = query->length + 1; // y = query = individual read
    int xlen = ref->length + 1;   // x = reference = consensus
    
    // Create matrices
    Matrix score = create_matrix(ylen, xlen);     // Main score matrix
    Matrix trace = create_matrix(ylen, xlen);     // Traceback matrix
    Matrix match = create_matrix(query->length, ref->length); // Match scores
    
    // New matrices for affine gap penalties
    Matrix gap_x = create_matrix(ylen, xlen);     // Gap in x direction (deletion)
    Matrix gap_y = create_matrix(ylen, xlen);     // Gap in y direction (insertion)
    
    // Precompute match scores
    for (int x = 0; x < ref->length; x++) {
        for (int y = 0; y < query->length; y++) {
            match.data[y][x] = compare_positions(&ref->positions[x], &query->positions[y]);
        }
    }
    
    // Initialize first row and column with affine gap penalties
    score.data[0][0] = 0;
    gap_x.data[0][0] = gap_y.data[0][0] = -INFINITY; // Can't start with a gap extension
    
    // Initialize first column (gaps in reference)
    for (int y = 1; y < ylen; y++) {
        score.data[y][0] = -INFINITY;  // Can't reach these cells via match
        gap_x.data[y][0] = -INFINITY;  // Can't have deletion in first column
        
        // Gap in y-direction (insertion)
        if (y == 1) {
            gap_y.data[y][0] = gap_open;  // First gap: open penalty
        } else {
            gap_y.data[y][0] = gap_y.data[y-1][0] + gap_extend;  // Extend penalty
        }
        
        // Overall score comes from gap_y for first column
        score.data[y][0] = gap_y.data[y][0];
        trace.data[y][0] = 3;  // Insertion (gap in reference)
    }
    
    // Initialize first row (gaps in query)
    for (int x = 1; x < xlen; x++) {
        score.data[0][x] = -INFINITY;  // Can't reach these cells via match
        gap_y.data[0][x] = -INFINITY;  // Can't have insertion in first row
        
        // Gap in x-direction (deletion)
        if (x == 1) {
            gap_x.data[0][x] = gap_open;  // First gap: open penalty
        } else {
            gap_x.data[0][x] = gap_x.data[0][x-1] + gap_extend;  // Extend penalty
        }
        
        // Overall score comes from gap_x for first row
        score.data[0][x] = gap_x.data[0][x];
        trace.data[0][x] = 2;  // Deletion (gap in query)
    }
    
    // Fill score and trace matrices
    for (int y = 1; y < ylen; y++) {
        for (int x = 1; x < xlen; x++) {
            // Calculate score for match/mismatch
            double diag = score.data[y-1][x-1] + match.data[y-1][x-1];
            
            // Calculate score for gap in x (deletion)
            double open_x = score.data[y][x-1] + gap_open;       // Open new gap
            double extend_x = gap_x.data[y][x-1] + gap_extend;   // Extend existing gap
            gap_x.data[y][x] = (open_x > extend_x) ? open_x : extend_x;
            
            // Calculate score for gap in y (insertion)
            double open_y = score.data[y-1][x] + gap_open;       // Open new gap
            double extend_y = gap_y.data[y-1][x] + gap_extend;   // Extend existing gap
            gap_y.data[y][x] = (open_y > extend_y) ? open_y : extend_y;
            
            // Find the best score among the three options
            if (diag >= gap_x.data[y][x] && diag >= gap_y.data[y][x]) {
                score.data[y][x] = diag;
                trace.data[y][x] = 1;  // Match/mismatch
            } else if (gap_x.data[y][x] >= gap_y.data[y][x]) {
                score.data[y][x] = gap_x.data[y][x];
                // Record whether this is an open or extend gap
                trace.data[y][x] = (open_x > extend_x) ? 2 : 4;  // 2=open deletion, 4=extend deletion
            } else {
                score.data[y][x] = gap_y.data[y][x];
                // Record whether this is an open or extend gap
                trace.data[y][x] = (open_y > extend_y) ? 3 : 5;  // 3=open insertion, 5=extend insertion
            }
        }
    }
    
    // Traceback
    int max_trace_len = query->length + ref->length + 1;
    char* traceback = malloc(max_trace_len * sizeof(char));
    if (!traceback) {
        free_matrix(score);
        free_matrix(trace);
        free_matrix(match);
        free_matrix(gap_x);
        free_matrix(gap_y);
        return NULL;
    }
    
    int idx = 0;
    int x = xlen - 1;
    int y = ylen - 1;
    
    while (x > 0 && y > 0) {
        // Bounds checking
        if (x >= xlen || y >= ylen) break;
        
        int trace_val = trace.data[y][x];
        traceback[idx++] = trace_val;
        
        if (trace_val == 1) {        // Match/mismatch
            x--;
            y--;
        } else if (trace_val == 2 || trace_val == 4) {  // Deletion (gap in query)
            x--;
        } else {                      // Insertion (gap in reference)
            y--;
        }
    }
    
    // Handle remaining gaps
    while (x > 0) {
        traceback[idx++] = 2;  // Use open gap code for simplicity in remaining sequence
        x--;
    }
    while (y > 0) {
        traceback[idx++] = 3;  // Use open gap code for simplicity in remaining sequence
        y--;
    }
    
    // Reverse the traceback
    for (int i = 0; i < idx/2; i++) {
        char temp = traceback[i];
        traceback[i] = traceback[idx-i-1];
        traceback[idx-i-1] = temp;
    }
    traceback[idx] = '\0';
    
    // Cleanup
    free_matrix(score);
    free_matrix(trace);
    free_matrix(match);
    free_matrix(gap_x);
    free_matrix(gap_y);
    
    return traceback;
}

SeqArray merge_seqs(const SeqArray* seq1, const SeqArray* seq2) {
    // seq1 = reference = consensus
    // seq2 = query = individual read
    char* trace = align_arrays(seq1, seq2);
    int trace_len = strlen(trace);
    
    SeqArray merged = create_seq_array(trace_len);
    int pos1 = 0, pos2 = 0;
    
    for(int i=0; i<trace_len; i++) {
        Position p1 = {0}, p2 = {0};
        
        // Handle sequence boundaries
        memcpy(p1.scores, seq1->positions[pos1].scores, sizeof(p1.scores));
        memcpy(p2.scores, seq2->positions[pos2].scores, sizeof(p2.scores));
        
        switch(trace[i]) {
            case 1:  // match
                for(int j=0; j<5; j++) {
                    merged.positions[i].scores[j] = p1.scores[j] + p2.scores[j];
                }
                pos1++;
                pos2++;
                break;
            case 2:  // gap in query / individual read (seq2)
                memcpy(merged.positions[i].scores, p1.scores, 5*sizeof(int)); // Copy ACGT- from reference
                for(int j=0; j<5; j++) {
                    merged.positions[i].scores[4] += p2.scores[j]; // Sum query's bases for gap support
                }
                pos1++;
                break;
            case 3:  // gap in reference / consensus (seq1)
                memcpy(merged.positions[i].scores, p2.scores, 5*sizeof(int)); // Copy ACGT- from query
                for(int j=0; j<5; j++) {
                    merged.positions[i].scores[4] += p1.scores[j]; // Sum ref's bases for gap support
                }
                pos2++;
                break;
        }
    }
    
    free(trace);
    return merged;
}

// Helper function to avoid underflow on q-scores
double log10addexp(double A, double B) {
    // Ensure A is the larger of the two values to avoid underflow
    if (A < B) {
        double temp = A;
        A = B;
        B = temp;
    }

    // Compute log10(10^A + 10^B) = A + log10(1 + 10^(B - A))
    double diff = B - A;
    return A + log10(1.0 + pow(10.0, diff));
}

// Convert consensus array to sequence and quality strings
void array_to_seq(const SeqArray* array, char** seq_out, char** qual_out) {
    int out_pos = 0;
    *seq_out = malloc(array->length + 1);
    *qual_out = malloc(array->length + 1);
    
    for(int i = 0; i < array->length; i++) {
        const Position* pos = &array->positions[i];
        
        // 1. Find global maximum (including gap)
        int max_idx = 0;
        int max_val = pos->scores[0];
        int sum = -pos->scores[4];
        for(int j = 0; j < 5; j++) {
            sum += pos->scores[j];
            if(pos->scores[j] > max_val) {
                max_val = pos->scores[j];
                max_idx = j;
            }
        }
        
        // Skip position if gap is most supported
        if(max_idx == 4) continue;
        
        // 2. Calculate quality scores
        (*seq_out)[out_pos] = BASES[max_idx];
        
        if(sum == 0) {
            (*qual_out)[out_pos] = '!';  // Q=0
            (*seq_out)[out_pos] = 'N';   // N for no coverage
        } else {
            double max_score = max_val;
            double other_sum = sum - max_score;
            
            // Bayesian calculation of output quality
            double log10_p_correct = max_score * q_correct + other_sum * q_error;
            double log10_p_error = max_score * q_error + other_sum * q_correct;
            double qual_d = -10.0 * log10_p_error + 10.0 * log10addexp(log10_p_correct, log10_p_error);
            
            // Clamp quality values
            int qual_int = (int)round(qual_d);
            qual_int = qual_int < 0 ? 0 : (qual_int > 40 ? 40 : qual_int);
            if(qual_int < 3) {
                (*seq_out)[out_pos] = 'N';  // N for low quality
            }
            (*qual_out)[out_pos] = (char)(qual_int + 33);
        }
        
        out_pos++;
    }
    
    // Terminate strings initially
    (*seq_out)[out_pos] = '\0';
    (*qual_out)[out_pos] = '\0';
    
    // Check if read needs to be trimmed from 3' end
    if(out_pos > 0 && ((*qual_out)[out_pos-1] - 33) <= 5) {
        // Trim iteratively from the end until a good base is found
        int trim_pos = out_pos - 1;
        while(trim_pos >= 0) {
            int qual = (*qual_out)[trim_pos] - 33;
            if(qual > 5) break;
            trim_pos--;
        }
        
        // Update length after trimming
        trim_pos++;  // Move forward one to include the last good base
        (*seq_out)[trim_pos] = '\0';
        (*qual_out)[trim_pos] = '\0';
        out_pos = trim_pos;
    }
    
    // Reallocate to trimmed size
    *seq_out = realloc(*seq_out, out_pos + 1);
    *qual_out = realloc(*qual_out, out_pos + 1);
}

void print_usage(const char* program_name) {
    fprintf(stderr, "Usage: %s [options]\n", program_name);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  --read1=file       Input FASTQ file 1\n");
    fprintf(stderr, "  --read2=file       Input FASTQ file 2\n");
    fprintf(stderr, "  --bc-start=int     Barcode start position (Zero-indexed, default: 0)\n");
    fprintf(stderr, "  --bc-len=int       Barcode length (default: 18)\n");
    fprintf(stderr, "  --gap-open=float   Gap open penalty (default: 3.0)\n");
    fprintf(stderr, "  --gap-extend=float Gap extend penalty (default: 3.0)\n");
    fprintf(stderr, "  --out1=file        Output file for consensus read 1\n");
    fprintf(stderr, "  --out2=file        Output file for consensus read 2\n");
    fprintf(stderr, "  --threads=int      Number of threads (default: 1)\n");
    exit(EXIT_FAILURE);
}

// Cleanup definitions
void free_barcode_group(BarcodeGroup* group) {
    for (int i = 0; i < group->count; i++) {
        free(group->fwd_reads[i]);
        free(group->rev_reads[i]);
        free(group->fwd_quals[i]);
        free(group->rev_quals[i]);
    }
    free(group->fwd_reads);
    free(group->rev_reads);
    free(group->fwd_quals);
    free(group->rev_quals);
    free(group);
}

void free_seq_array(SeqArray* arr) {
    free(arr->positions);
}

void free_matrix(Matrix mat) {
    for (int i = 0; i < mat.rows; i++) {
        free(mat.data[i]);
    }
    free(mat.data);
}

