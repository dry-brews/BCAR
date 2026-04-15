#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "seq_module.h"

// Constants used only within seq_module
#define q_correct -0.004575749 // = 1/10 * log10(0.9)
#define q_error -0.1 // = 1/10 * log10(0.1)

static const char* bases = BASES;
static const int basemap[256] = { ['A']=0, ['C']=1, ['G']=2, ['T']=3, ['-']=4, ['N']=4 };
static const int q_adj[41] = {
    0, 0, 0, 0, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40, 41, 42
};

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

void free_seq_array(SeqArray* arr) {
    free(arr->positions);
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

// Convert consensus array back to a sequence + quality.
// Sets *minor_frac_out to the highest minor allele frequency observed at any position.
void array_to_seq(const SeqArray* array, char** seq_out, char** qual_out,
                  double* minor_frac_out) {
    int out_pos = 0;
    double max_f2 = 0.0;
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

        // Mixture detection: check if a secondary non-gap base has substantial evidence.
        // Only evaluated at positions where non-gap evidence exists and gap is not dominant.
        if (max_idx != 4 && total_non_gap > 0) {
            int second_val = 0;
            for (int j = 0; j < 4; j++) {
                if (j == max_idx) continue;
                if (pos->scores[j] > second_val) second_val = pos->scores[j];
            }
            if (second_val > 0) {
                double f2 = (double)second_val / total_non_gap;
                if (f2 > max_f2) max_f2 = f2;
            }
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

    // Report highest minor allele frequency
    if (minor_frac_out)
        *minor_frac_out = max_f2;
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

    #undef IN_BAND
    #undef IDX

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
