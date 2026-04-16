#ifndef SORT_MODULE_H
#define SORT_MODULE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <zlib.h>

/* ------------------------------------------------------------------ */
/*  Constants                                                         */
/* ------------------------------------------------------------------ */

#define MAX_LINE 65536
#define GZ_BUF_SIZE (1024 * 1024)
#define HT_INITIAL_CAP (1 << 20)
#define HT_LOAD_FACTOR 0.7
#define MAX_BC_LEN 64

/* ------------------------------------------------------------------ */
/*  Structures                                                        */
/* ------------------------------------------------------------------ */

/* Full barcode with N-mask */
typedef struct {
    uint64_t hi;        /* 2-bit packed bases, positions 32-63 */
    uint64_t lo;        /* 2-bit packed bases, positions 0-31  */
    uint64_t n_mask_hi; /* 1 = position is N (upper half)      */
    uint64_t n_mask_lo; /* 1 = position is N (lower half)      */
    uint8_t  len;
} barcode_t;

/* Hash table key: concrete barcodes only (no N's) */
typedef struct {
    uint64_t hi;
    uint64_t lo;
    uint8_t  len;
} bc_key_t;

/* Hash table entry (Robin Hood open addressing) */
typedef struct {
    bc_key_t bc;
    uint32_t count;
    uint32_t uf_idx;
    uint16_t psl;       /* probe sequence length */
    uint8_t  occupied;
} ht_entry_t;

/* Hash table */
typedef struct {
    ht_entry_t *slots;
    uint64_t capacity;
    uint64_t count;
    uint64_t mask;
} barcode_ht_t;

/* Union-Find */
typedef struct {
    uint32_t *parent;
    uint32_t *rank;
    uint32_t *best_idx; /* ht slot index of highest-count member */
    uint32_t  n;
} union_find_t;

/* Context pattern for barcode extraction */
typedef struct {
    char prefix[256];
    char suffix[256];
    int  prefix_len;
    int  suffix_len;
    int  bc_len;        /* number of N's in pattern */
    int  anchored_5;
    int  anchored_3;
} context_pattern_t;

/* Barcode extraction configuration */
typedef struct {
    int mode;               /* 0 = fixed position, 1 = context */
    int bc_start;
    int bc_len;
    context_pattern_t ctx;
    int max_context_mm;
    int max_bc_indels;
} extraction_config_t;

/* Extraction result */
typedef struct {
    char barcode[MAX_BC_LEN + 1];
    int  bc_len;
    int  found;
} extraction_result_t;

/* Record for chunk sorting */
typedef struct {
    char    *header;
    char    *seq;
    char    *plus;
    char    *qual;
    uint64_t ubid;
} sort_record_t;

/* Chunk buffer */
typedef struct {
    sort_record_t *recs;
    long  n;
    long  cap;
    long  mem_used;
} chunk_buffer_t;

/* Merge stream for k-way merge */
typedef struct {
    FILE    *fp;
    char     hdr[MAX_LINE];
    char     seq[MAX_LINE];
    char     plus_line[MAX_LINE];
    char     qual[MAX_LINE];
    uint64_t ubid;
    int      exhausted;
} merge_stream_t;

/* Min-heap for k-way merge */
typedef struct {
    merge_stream_t **arr;
    int n;
    int cap;
} min_heap_t;

/* Buffered gzip reader */
typedef struct {
    gzFile  gz;
    char    buf[GZ_BUF_SIZE];
    int     pos;
    int     len;
    int     is_gz;
    FILE   *fp;
} gz_reader_t;

/* ------------------------------------------------------------------ */
/*  Memory utilities                                                  */
/* ------------------------------------------------------------------ */

void fatal(const char *msg);
void fatal_alloc(const char *where);
void *xmalloc(size_t size, const char *where);
void *xcalloc(size_t n, size_t size, const char *where);
void *xrealloc(void *ptr, size_t size, const char *where);
char *xstrdup(const char *s, const char *where);

/* ------------------------------------------------------------------ */
/*  Gzip-aware FASTQ reader                                           */
/* ------------------------------------------------------------------ */

gz_reader_t *gz_reader_open(const char *filename);
void gz_reader_close(gz_reader_t *r);
int gz_getline(gz_reader_t *r, char *dest, int maxlen);
int read_fastq(gz_reader_t *r, char *hdr, char *seq, char *plus, char *qual);

/* ------------------------------------------------------------------ */
/*  Barcode encoding / decoding                                       */
/* ------------------------------------------------------------------ */

int encode_barcode(const char *seq, int len, barcode_t *out);
void barcode_to_key(const barcode_t *bc, bc_key_t *key);
void decode_key(const bc_key_t *key, char *out);

/* ------------------------------------------------------------------ */
/*  Hash table                                                        */
/* ------------------------------------------------------------------ */

void ht_init(barcode_ht_t *ht, uint64_t cap);
void ht_free(barcode_ht_t *ht);
void ht_insert(barcode_ht_t *ht, const bc_key_t *bc);
ht_entry_t *ht_lookup(barcode_ht_t *ht, const bc_key_t *bc);

/* ------------------------------------------------------------------ */
/*  Union-Find                                                        */
/* ------------------------------------------------------------------ */

void uf_init(union_find_t *uf, uint32_t n);
void uf_free(union_find_t *uf);
uint32_t uf_find(union_find_t *uf, uint32_t x);

/* ------------------------------------------------------------------ */
/*  Barcode extraction                                                */
/* ------------------------------------------------------------------ */

int extract_barcode(const char *seq, int slen,
                    const extraction_config_t *cfg,
                    extraction_result_t *out);

/* ------------------------------------------------------------------ */
/*  Clustering                                                        */
/* ------------------------------------------------------------------ */

void cluster_barcodes(barcode_ht_t *ht, union_find_t *uf,
                      int max_mismatches, int max_bc_indels);
uint64_t assign_ubids(barcode_ht_t *ht, union_find_t *uf,
                      uint64_t **ubid_map_out);
void check_cluster_diameters(const barcode_ht_t *ht, union_find_t *uf,
                              const uint64_t *ubid_map,
                              int max_mismatches, long total_reads);

/* ------------------------------------------------------------------ */
/*  N-containing barcode resolution                                   */
/* ------------------------------------------------------------------ */

uint64_t resolve_n_barcode(const barcode_t *bc, const barcode_ht_t *ht,
                           const uint64_t *ubid_map);

/* ------------------------------------------------------------------ */
/*  Chunk buffer                                                      */
/* ------------------------------------------------------------------ */

void chunk_init(chunk_buffer_t *cb, long initial_cap);
void chunk_push(chunk_buffer_t *cb, const char *hdr, const char *seq,
                const char *plus, const char *qual, uint64_t ubid);
void chunk_clear(chunk_buffer_t *cb);
void chunk_free(chunk_buffer_t *cb);
char *write_chunk(chunk_buffer_t *cb, int chunk_num, const char *tmp_dir);

/* ------------------------------------------------------------------ */
/*  K-way merge                                                       */
/* ------------------------------------------------------------------ */

void merge_sorted_chunks(char **paths, int n_chunks,
                         const char *out_path, int max_open_files,
                         const char *tmp_dir);

/* Callback-based merge: delivers completed UBID groups instead of writing
   to a file. Callback receives arrays of strdup'd strings that the callback
   owns (must free). */
typedef void (*merge_group_cb_t)(const char *ubid_str,
                                 char **seqs, char **quals,
                                 int count, void *user_data);

void merge_sorted_chunks_with_callback(
    char **paths, int n_chunks,
    int max_open_files, const char *tmp_dir,
    merge_group_cb_t cb, void *user_data);

#endif /* SORT_MODULE_H */
