#ifndef SEQ_MODULE_H
#define SEQ_MODULE_H

#include <stdbool.h>

// Constants
#define MAX_BARCODE_LEN 64
#define DEFAULT_MAX_LINE_LEN 100000
#define BASES "ACGT-"
#define DEFAULT_GAP_SCORE -1.0
#define VECTOR_LENGTH 5

// Structures
typedef struct {
    int scores[5]; // A, C, G, T, -
} Position;

typedef struct {
    Position* positions;
    int length;
} SeqArray;

typedef struct {
    char bc[MAX_BARCODE_LEN];
    char** fwd_reads;
    char** fwd_quals;
    int count;
    int capacity;
} BarcodeGroup;

typedef struct { int idx; double sim; } IdxSim;

// Globals defined in bc_merger_v14.c, used by seq_module functions
extern double gap_pen;
extern int max_line_len;
extern double call_total;
extern double call_missense;
extern double call_indel;
// Seq module functions
SeqArray seq_to_array(const char* seq, const char* qual, int length);
void array_to_seq(const SeqArray* array, char** seq_out, char** qual_out,
                  double* minor_frac_out);
double compare_positions(const Position* a, const Position* b);
double compare_seqs(const SeqArray* a, const SeqArray* b);
SeqArray build_unaligned_consensus(SeqArray* sequences, int count);
int* align_arrays(const SeqArray* ref, const SeqArray* query, int *out_len);
int* align_arrays_band(const SeqArray* ref, const SeqArray* query, int *out_len, int max_phase_diff);
SeqArray merge_seqs(const SeqArray* seq1, const SeqArray* seq2);
SeqArray create_seq_array(int length);
void free_seq_array(SeqArray* arr);
double log10addexp_vec(const double *vals, int n);

// fatal_alloc: provided by sort_module.c (or bc_merger_v14.c in standalone builds)
void fatal_alloc(const char* where);

#endif // SEQ_MODULE_H
