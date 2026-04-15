#include "sort_module.h"
#include <errno.h>

/* ------------------------------------------------------------------ */
/*  File-private constants                                            */
/* ------------------------------------------------------------------ */

static const int basemap[256] = {
    ['A'] = 0, ['a'] = 0,
    ['C'] = 1, ['c'] = 1,
    ['G'] = 2, ['g'] = 2,
    ['T'] = 3, ['t'] = 3
};
static const char decode_base[4] = "ACGT";

/* ------------------------------------------------------------------ */
/*  Memory utilities                                                  */
/* ------------------------------------------------------------------ */

void fatal(const char *msg) {
    fprintf(stderr, "Fatal: %s\n", msg);
    exit(EXIT_FAILURE);
}

void fatal_alloc(const char *where) {
    fprintf(stderr, "Fatal: memory allocation failed in %s\n", where);
    exit(EXIT_FAILURE);
}

void *xmalloc(size_t size, const char *where) {
    void *p = malloc(size);
    if (!p) fatal_alloc(where);
    return p;
}

void *xcalloc(size_t n, size_t size, const char *where) {
    void *p = calloc(n, size);
    if (!p) fatal_alloc(where);
    return p;
}

void *xrealloc(void *ptr, size_t size, const char *where) {
    void *p = realloc(ptr, size);
    if (!p) fatal_alloc(where);
    return p;
}

char *xstrdup(const char *s, const char *where) {
    char *p = strdup(s);
    if (!p) fatal_alloc(where);
    return p;
}

/* ------------------------------------------------------------------ */
/*  Gzip-aware FASTQ reader                                           */
/* ------------------------------------------------------------------ */

static int is_gzipped(const char *filename) {
    FILE *f = fopen(filename, "rb");
    if (!f) return 0;
    unsigned char magic[2] = {0, 0};
    fread(magic, 1, 2, f);
    fclose(f);
    return (magic[0] == 0x1f && magic[1] == 0x8b);
}

gz_reader_t *gz_reader_open(const char *filename) {
    gz_reader_t *r = xcalloc(1, sizeof(gz_reader_t), "gz_reader_open");
    r->is_gz = is_gzipped(filename);
    if (r->is_gz) {
        r->gz = gzopen(filename, "r");
        if (!r->gz) {
            fprintf(stderr, "Fatal: cannot open gzipped file: %s\n", filename);
            exit(EXIT_FAILURE);
        }
        gzbuffer(r->gz, GZ_BUF_SIZE);
    } else {
        r->fp = fopen(filename, "r");
        if (!r->fp) {
            fprintf(stderr, "Fatal: cannot open file: %s\n", filename);
            exit(EXIT_FAILURE);
        }
    }
    r->pos = 0;
    r->len = 0;
    return r;
}

void gz_reader_close(gz_reader_t *r) {
    if (!r) return;
    if (r->is_gz) gzclose(r->gz);
    else if (r->fp) fclose(r->fp);
    free(r);
}

static int gz_fill(gz_reader_t *r) {
    r->len = gzread(r->gz, r->buf, GZ_BUF_SIZE);
    r->pos = 0;
    return r->len > 0;
}

/* Read one line into dest (no newline). Returns 1 on success, 0 on EOF. */
int gz_getline(gz_reader_t *r, char *dest, int maxlen) {
    int out = 0;
    if (r->is_gz) {
        while (1) {
            if (r->pos >= r->len) {
                if (!gz_fill(r)) {
                    dest[out] = '\0';
                    return out > 0;
                }
            }
            char *start = r->buf + r->pos;
            int avail = r->len - r->pos;
            char *nl = memchr(start, '\n', avail);
            if (nl) {
                int chunk = (int)(nl - start);
                if (out + chunk < maxlen) {
                    memcpy(dest + out, start, chunk);
                    out += chunk;
                }
                r->pos = (int)(nl - r->buf) + 1;
                /* strip \r */
                if (out > 0 && dest[out - 1] == '\r') out--;
                dest[out] = '\0';
                return 1;
            }
            int chunk = avail;
            if (out + chunk >= maxlen) chunk = maxlen - out - 1;
            memcpy(dest + out, start, chunk);
            out += chunk;
            r->pos = r->len;
        }
    } else {
        if (!fgets(dest, maxlen, r->fp)) {
            dest[0] = '\0';
            return 0;
        }
        int len = (int)strlen(dest);
        if (len > 0 && dest[len - 1] == '\n') dest[--len] = '\0';
        if (len > 0 && dest[len - 1] == '\r') dest[--len] = '\0';
        return 1;
    }
}

/* Read one FASTQ record. Returns 1 on success, 0 on EOF. */
int read_fastq(gz_reader_t *r, char *hdr, char *seq, char *plus, char *qual) {
    if (!gz_getline(r, hdr, MAX_LINE)) return 0;
    if (!gz_getline(r, seq, MAX_LINE)) return 0;
    if (!gz_getline(r, plus, MAX_LINE)) return 0;
    if (!gz_getline(r, qual, MAX_LINE)) return 0;
    return 1;
}

/* ------------------------------------------------------------------ */
/*  Barcode encoding / decoding                                       */
/* ------------------------------------------------------------------ */

/* Encode ASCII barcode to barcode_t. Returns number of N's found, or
   -1 if the barcode contains invalid characters. */
int encode_barcode(const char *seq, int len, barcode_t *out) {
    memset(out, 0, sizeof(*out));
    out->len = (uint8_t)len;
    int n_count = 0;
    for (int i = 0; i < len; i++) {
        char c = seq[i];
        int base;
        if (c == 'N' || c == 'n') {
            base = 0; /* placeholder */
            n_count++;
            if (i < 32)
                out->n_mask_lo |= (1ULL << i);
            else
                out->n_mask_hi |= (1ULL << (i - 32));
        } else {
            base = basemap[(unsigned char)c];
            if (c != 'A' && c != 'a' && c != 'C' && c != 'c' &&
                c != 'G' && c != 'g' && c != 'T' && c != 't')
                return -1;
        }
        if (i < 32)
            out->lo |= ((uint64_t)base << (2 * i));
        else
            out->hi |= ((uint64_t)base << (2 * (i - 32)));
    }
    return n_count;
}

/* Convert barcode_t with no N's to bc_key_t */
void barcode_to_key(const barcode_t *bc, bc_key_t *key) {
    key->hi = bc->hi;
    key->lo = bc->lo;
    key->len = bc->len;
}

/* Decode bc_key_t to ASCII string */
void decode_key(const bc_key_t *key, char *out) {
    for (int i = 0; i < key->len; i++) {
        int base;
        if (i < 32)
            base = (key->lo >> (2 * i)) & 3;
        else
            base = (key->hi >> (2 * (i - 32))) & 3;
        out[i] = decode_base[base];
    }
    out[key->len] = '\0';
}

/* ------------------------------------------------------------------ */
/*  Hash table (Robin Hood open addressing)                           */
/* ------------------------------------------------------------------ */

static uint64_t hash_key(const bc_key_t *k) {
    uint64_t h = k->lo ^ (k->hi * 0x9E3779B97F4A7C15ULL) ^ ((uint64_t)k->len << 57);
    h ^= h >> 30; h *= 0xbf58476d1ce4e5b9ULL;
    h ^= h >> 27; h *= 0x94d049bb133111ebULL;
    h ^= h >> 31;
    return h;
}

static int key_equal(const bc_key_t *a, const bc_key_t *b) {
    return a->len == b->len && a->lo == b->lo && a->hi == b->hi;
}

static void ht_grow(barcode_ht_t *ht);

/* Insert a new entry (used during ht_grow). Does not increment;
   caller sets count/uf_idx on the returned pointer. */
static ht_entry_t *ht_insert_new(barcode_ht_t *ht, const bc_key_t *bc,
                                  uint32_t count, uint32_t uf_idx) {
    uint64_t idx = hash_key(bc) & ht->mask;
    ht_entry_t incoming;
    incoming.bc = *bc;
    incoming.count = count;
    incoming.uf_idx = uf_idx;
    incoming.psl = 0;
    incoming.occupied = 1;

    while (1) {
        ht_entry_t *slot = &ht->slots[idx];
        if (!slot->occupied) {
            *slot = incoming;
            ht->count++;
            return slot;
        }
        /* Robin Hood: steal from the rich */
        if (slot->psl < incoming.psl) {
            ht_entry_t tmp = *slot;
            *slot = incoming;
            incoming = tmp;
        }
        incoming.psl++;
        idx = (idx + 1) & ht->mask;
    }
}

void ht_init(barcode_ht_t *ht, uint64_t cap) {
    /* Round up to power of 2 */
    uint64_t c = 1;
    while (c < cap) c <<= 1;
    ht->capacity = c;
    ht->mask = c - 1;
    ht->count = 0;
    ht->slots = xcalloc(c, sizeof(ht_entry_t), "ht_init");
}

void ht_free(barcode_ht_t *ht) {
    free(ht->slots);
    ht->slots = NULL;
}

/* Insert or increment count for bc. */
void ht_insert(barcode_ht_t *ht, const bc_key_t *bc) {
    if ((double)(ht->count + 1) / ht->capacity > HT_LOAD_FACTOR)
        ht_grow(ht);

    /* First check if key already exists */
    ht_entry_t *existing = ht_lookup(ht, bc);
    if (existing) {
        existing->count++;
        return;
    }

    /* Key not found: insert new entry */
    ht_insert_new(ht, bc, 1, 0);
}

static void ht_grow(barcode_ht_t *ht) {
    uint64_t old_cap = ht->capacity;
    ht_entry_t *old_slots = ht->slots;

    ht->capacity = old_cap * 2;
    ht->mask = ht->capacity - 1;
    ht->count = 0;
    ht->slots = xcalloc(ht->capacity, sizeof(ht_entry_t), "ht_grow");

    for (uint64_t i = 0; i < old_cap; i++) {
        if (old_slots[i].occupied) {
            ht_insert_new(ht, &old_slots[i].bc,
                          old_slots[i].count, old_slots[i].uf_idx);
        }
    }
    free(old_slots);
}

/* Lookup. Returns NULL if not found. */
ht_entry_t *ht_lookup(barcode_ht_t *ht, const bc_key_t *bc) {
    uint64_t idx = hash_key(bc) & ht->mask;
    uint16_t dist = 0;
    while (1) {
        ht_entry_t *slot = &ht->slots[idx];
        if (!slot->occupied) return NULL;
        if (slot->psl < dist) return NULL;
        if (key_equal(&slot->bc, bc)) return slot;
        dist++;
        idx = (idx + 1) & ht->mask;
    }
}

/* ------------------------------------------------------------------ */
/*  Union-Find                                                        */
/* ------------------------------------------------------------------ */

void uf_init(union_find_t *uf, uint32_t n) {
    uf->n = n;
    uf->parent = xmalloc(n * sizeof(uint32_t), "uf_init parent");
    uf->rank = xcalloc(n, sizeof(uint32_t), "uf_init rank");
    uf->best_idx = xmalloc(n * sizeof(uint32_t), "uf_init best_idx");
    for (uint32_t i = 0; i < n; i++) {
        uf->parent[i] = i;
        uf->best_idx[i] = i;
    }
}

void uf_free(union_find_t *uf) {
    free(uf->parent);
    free(uf->rank);
    free(uf->best_idx);
}

uint32_t uf_find(union_find_t *uf, uint32_t x) {
    while (uf->parent[x] != x) {
        uf->parent[x] = uf->parent[uf->parent[x]]; /* path halving */
        x = uf->parent[x];
    }
    return x;
}

/* Union two sets. Updates best_idx to track highest-count member.
   counts[uf_idx] holds the count for each UF element. */
static void uf_union(union_find_t *uf, uint32_t x, uint32_t y,
                     const uint32_t *counts) {
    uint32_t rx = uf_find(uf, x);
    uint32_t ry = uf_find(uf, y);
    if (rx == ry) return;

    /* Union by rank */
    if (uf->rank[rx] < uf->rank[ry]) { uint32_t t = rx; rx = ry; ry = t; }
    uf->parent[ry] = rx;
    if (uf->rank[rx] == uf->rank[ry]) uf->rank[rx]++;

    /* Track highest-count representative */
    uint32_t bi_x = uf->best_idx[rx];
    uint32_t bi_y = uf->best_idx[ry];
    if (counts[bi_y] > counts[bi_x])
        uf->best_idx[rx] = bi_y;
}

/* ------------------------------------------------------------------ */
/*  Barcode extraction                                                */
/* ------------------------------------------------------------------ */

/* Count substitutions between text at position pos and pattern.
   Returns error count, short-circuits at limit+1. */
static int count_subst(const char *text, int pos, const char *pattern,
                       int plen, int limit) {
    int errs = 0;
    for (int i = 0; i < plen; i++) {
        char a = text[pos + i];
        char b = pattern[i];
        /* Case-insensitive comparison */
        if ((a | 0x20) != (b | 0x20)) {
            if (++errs > limit) return errs;
        }
    }
    return errs;
}

/* Find best position of pattern in text[search_start..search_end-plen].
   Returns position or -1 if not found within max_errors. */
static int fuzzy_find(const char *text, int tlen, int search_start, int search_end,
                      const char *pattern, int plen, int max_errors) {
    if (plen == 0) return search_start; /* empty pattern always matches */
    int best_pos = -1;
    int best_err = max_errors + 1;
    int limit = search_end - plen;
    if (limit >= tlen - plen) limit = tlen - plen;
    if (search_start < 0) search_start = 0;

    for (int pos = search_start; pos <= limit; pos++) {
        int e = count_subst(text, pos, pattern, plen, best_err - 1);
        if (e < best_err) {
            best_err = e;
            best_pos = pos;
            if (e == 0) break;
        }
    }
    return best_pos;
}

/* Extract barcode from read. Returns 0 on success, -1 on failure. */
int extract_barcode(const char *seq, int slen,
                    const extraction_config_t *cfg,
                    extraction_result_t *out) {
    out->found = 0;

    if (cfg->mode == 0) {
        /* Fixed position mode */
        int end = cfg->bc_start + cfg->bc_len;
        if (end > slen) return -1;
        memcpy(out->barcode, seq + cfg->bc_start, cfg->bc_len);
        out->barcode[cfg->bc_len] = '\0';
        out->bc_len = cfg->bc_len;
        out->found = 1;
        return 0;
    }

    /* Context mode */
    const context_pattern_t *ctx = &cfg->ctx;
    int prefix_pos = -1;
    int prefix_end;

    if (ctx->prefix_len > 0) {
        int s_start = 0;
        int s_end = slen;
        if (ctx->anchored_5) {
            s_start = 0;
            s_end = ctx->prefix_len + cfg->max_context_mm + 1;
            if (s_end > slen) s_end = slen;
        }
        prefix_pos = fuzzy_find(seq, slen, s_start, s_end,
                                ctx->prefix, ctx->prefix_len, cfg->max_context_mm);
        if (prefix_pos < 0) return -1;
        if (ctx->anchored_5 && prefix_pos != 0) return -1;
        prefix_end = prefix_pos + ctx->prefix_len;
    } else {
        if (ctx->anchored_5)
            prefix_end = 0;
        else
            return -1; /* need at least a prefix or anchored start */
    }

    int barcode_start = prefix_end;
    int expected_bc_len = ctx->bc_len;

    if (ctx->suffix_len > 0) {
        /* Search for suffix in a window around expected position */
        int expected_suffix_start = barcode_start + expected_bc_len;
        int search_lo = expected_suffix_start - cfg->max_bc_indels;
        int search_hi = expected_suffix_start + cfg->max_bc_indels + ctx->suffix_len;
        if (search_lo < barcode_start) search_lo = barcode_start;
        if (search_hi > slen) search_hi = slen;

        int suffix_pos = fuzzy_find(seq, slen, search_lo, search_hi,
                                    ctx->suffix, ctx->suffix_len, cfg->max_context_mm);
        if (suffix_pos < 0) return -1;
        if (ctx->anchored_3 && suffix_pos + ctx->suffix_len != slen) return -1;

        int actual_bc_len = suffix_pos - barcode_start;
        if (actual_bc_len < 1) return -1;
        int indel_diff = actual_bc_len - expected_bc_len;
        if (indel_diff < 0) indel_diff = -indel_diff;
        if (indel_diff > cfg->max_bc_indels) return -1;

        out->bc_len = actual_bc_len;
    } else {
        /* No suffix */
        if (ctx->anchored_3) {
            int actual_bc_len = slen - barcode_start;
            int indel_diff = actual_bc_len - expected_bc_len;
            if (indel_diff < 0) indel_diff = -indel_diff;
            if (indel_diff > cfg->max_bc_indels) return -1;
            out->bc_len = actual_bc_len;
        } else {
            out->bc_len = expected_bc_len;
        }
    }

    if (barcode_start + out->bc_len > slen) return -1;
    if (out->bc_len > MAX_BC_LEN) return -1;
    memcpy(out->barcode, seq + barcode_start, out->bc_len);
    out->barcode[out->bc_len] = '\0';
    out->found = 1;
    return 0;
}

/* ------------------------------------------------------------------ */
/*  Neighborhood enumeration and clustering                           */
/* ------------------------------------------------------------------ */

/* Enumerate substitution neighbors at distance 1 and union with matches */
static void enum_subs_d1(const bc_key_t *bc, uint32_t bc_uf,
                         barcode_ht_t *ht, union_find_t *uf,
                         const uint32_t *counts) {
    for (int i = 0; i < bc->len; i++) {
        int orig;
        if (i < 32)
            orig = (bc->lo >> (2 * i)) & 3;
        else
            orig = (bc->hi >> (2 * (i - 32))) & 3;

        for (int alt = 0; alt < 4; alt++) {
            if (alt == orig) continue;
            bc_key_t nb = *bc;
            /* Set position i to alt */
            if (i < 32)
                nb.lo = (bc->lo & ~(3ULL << (2 * i))) | ((uint64_t)alt << (2 * i));
            else
                nb.hi = (bc->hi & ~(3ULL << (2 * (i - 32)))) | ((uint64_t)alt << (2 * (i - 32)));
            ht_entry_t *e = ht_lookup(ht, &nb);
            if (e) uf_union(uf, bc_uf, e->uf_idx, counts);
        }
    }
}

/* Enumerate substitution neighbors at distance 2 and union with matches */
static void enum_subs_d2(const bc_key_t *bc, uint32_t bc_uf,
                         barcode_ht_t *ht, union_find_t *uf,
                         const uint32_t *counts) {
    int L = bc->len;
    for (int i = 0; i < L; i++) {
        int orig_i;
        if (i < 32) orig_i = (bc->lo >> (2 * i)) & 3;
        else orig_i = (bc->hi >> (2 * (i - 32))) & 3;

        for (int j = i + 1; j < L; j++) {
            int orig_j;
            if (j < 32) orig_j = (bc->lo >> (2 * j)) & 3;
            else orig_j = (bc->hi >> (2 * (j - 32))) & 3;

            for (int ai = 0; ai < 4; ai++) {
                if (ai == orig_i) continue;
                for (int aj = 0; aj < 4; aj++) {
                    if (aj == orig_j) continue;
                    bc_key_t nb = *bc;
                    if (i < 32)
                        nb.lo = (nb.lo & ~(3ULL << (2 * i))) | ((uint64_t)ai << (2 * i));
                    else
                        nb.hi = (nb.hi & ~(3ULL << (2 * (i - 32)))) | ((uint64_t)ai << (2 * (i - 32)));
                    if (j < 32)
                        nb.lo = (nb.lo & ~(3ULL << (2 * j))) | ((uint64_t)aj << (2 * j));
                    else
                        nb.hi = (nb.hi & ~(3ULL << (2 * (j - 32)))) | ((uint64_t)aj << (2 * (j - 32)));

                    ht_entry_t *e = ht_lookup(ht, &nb);
                    if (e) uf_union(uf, bc_uf, e->uf_idx, counts);
                }
            }
        }
    }
}

/* Helper: build a bc_key_t by removing the base at position del from bc */
static void make_deletion(const bc_key_t *bc, int del, bc_key_t *out) {
    out->len = bc->len - 1;
    out->lo = 0;
    out->hi = 0;
    int dst = 0;
    for (int i = 0; i < bc->len; i++) {
        if (i == del) continue;
        int base;
        if (i < 32) base = (bc->lo >> (2 * i)) & 3;
        else base = (bc->hi >> (2 * (i - 32))) & 3;
        if (dst < 32)
            out->lo |= ((uint64_t)base << (2 * dst));
        else
            out->hi |= ((uint64_t)base << (2 * (dst - 32)));
        dst++;
    }
}

/* Helper: build a bc_key_t by inserting base at position ins into bc */
static void make_insertion(const bc_key_t *bc, int ins, int base, bc_key_t *out) {
    out->len = bc->len + 1;
    out->lo = 0;
    out->hi = 0;
    int dst = 0;
    for (int i = 0; i < bc->len + 1; i++) {
        int b;
        if (i == ins) {
            b = base;
        } else {
            int src = (i < ins) ? i : i - 1;
            if (src < 32) b = (bc->lo >> (2 * src)) & 3;
            else b = (bc->hi >> (2 * (src - 32))) & 3;
        }
        if (dst < 32)
            out->lo |= ((uint64_t)b << (2 * dst));
        else
            out->hi |= ((uint64_t)b << (2 * (dst - 32)));
        dst++;
    }
}

/* Enumerate deletion neighbors (1 and 2 deletions) */
static void enum_deletions(const bc_key_t *bc, uint32_t bc_uf,
                           barcode_ht_t *ht, union_find_t *uf,
                           const uint32_t *counts, int max_indels) {
    if (bc->len < 2) return;

    /* 1 deletion */
    for (int i = 0; i < bc->len; i++) {
        bc_key_t nb;
        make_deletion(bc, i, &nb);
        ht_entry_t *e = ht_lookup(ht, &nb);
        if (e) uf_union(uf, bc_uf, e->uf_idx, counts);
    }

    /* 2 deletions */
    if (max_indels >= 2 && bc->len >= 3) {
        for (int i = 0; i < bc->len; i++) {
            bc_key_t nb1;
            make_deletion(bc, i, &nb1);
            for (int j = 0; j < nb1.len; j++) {
                bc_key_t nb2;
                make_deletion(&nb1, j, &nb2);
                ht_entry_t *e = ht_lookup(ht, &nb2);
                if (e) uf_union(uf, bc_uf, e->uf_idx, counts);
            }
        }
    }
}

/* Enumerate insertion neighbors (1 and 2 insertions) */
static void enum_insertions(const bc_key_t *bc, uint32_t bc_uf,
                            barcode_ht_t *ht, union_find_t *uf,
                            const uint32_t *counts, int max_indels) {
    if (bc->len >= MAX_BC_LEN) return;

    /* 1 insertion */
    for (int i = 0; i <= bc->len; i++) {
        for (int base = 0; base < 4; base++) {
            bc_key_t nb;
            make_insertion(bc, i, base, &nb);
            ht_entry_t *e = ht_lookup(ht, &nb);
            if (e) uf_union(uf, bc_uf, e->uf_idx, counts);
        }
    }

    /* 2 insertions */
    if (max_indels >= 2 && bc->len + 2 <= MAX_BC_LEN) {
        for (int i = 0; i <= bc->len; i++) {
            for (int b1 = 0; b1 < 4; b1++) {
                bc_key_t nb1;
                make_insertion(bc, i, b1, &nb1);
                for (int j = 0; j <= nb1.len; j++) {
                    for (int b2 = 0; b2 < 4; b2++) {
                        bc_key_t nb2;
                        make_insertion(&nb1, j, b2, &nb2);
                        ht_entry_t *e = ht_lookup(ht, &nb2);
                        if (e) uf_union(uf, bc_uf, e->uf_idx, counts);
                    }
                }
            }
        }
    }
}

/* Run full clustering. Builds a uf_idx-indexed counts array from the HT
   so that uf_union can track the highest-count representative. */
void cluster_barcodes(barcode_ht_t *ht, union_find_t *uf,
                      int max_mismatches, int max_bc_indels) {
    /* Build counts array indexed by uf_idx */
    uint32_t *counts = xmalloc(uf->n * sizeof(uint32_t), "cluster counts");
    for (uint64_t i = 0; i < ht->capacity; i++) {
        if (!ht->slots[i].occupied) continue;
        counts[ht->slots[i].uf_idx] = ht->slots[i].count;
    }

    fprintf(stderr, "  Enumerating substitution neighbors (d=1)...\n");
    for (uint64_t i = 0; i < ht->capacity; i++) {
        if (!ht->slots[i].occupied) continue;
        enum_subs_d1(&ht->slots[i].bc, ht->slots[i].uf_idx, ht, uf, counts);
    }

    if (max_mismatches >= 2) {
        fprintf(stderr, "  Enumerating substitution neighbors (d=2)...\n");
        for (uint64_t i = 0; i < ht->capacity; i++) {
            if (!ht->slots[i].occupied) continue;
            enum_subs_d2(&ht->slots[i].bc, ht->slots[i].uf_idx, ht, uf, counts);
        }
    }

    if (max_bc_indels > 0) {
        fprintf(stderr, "  Enumerating deletion neighbors...\n");
        for (uint64_t i = 0; i < ht->capacity; i++) {
            if (!ht->slots[i].occupied) continue;
            enum_deletions(&ht->slots[i].bc, ht->slots[i].uf_idx, ht, uf, counts, max_bc_indels);
        }
        fprintf(stderr, "  Enumerating insertion neighbors...\n");
        for (uint64_t i = 0; i < ht->capacity; i++) {
            if (!ht->slots[i].occupied) continue;
            enum_insertions(&ht->slots[i].bc, ht->slots[i].uf_idx, ht, uf, counts, max_bc_indels);
        }
    }

    free(counts);
}

/* Assign sequential UBIDs. Returns number of clusters.
   ubid_map[ht slot index] = UBID */
uint64_t assign_ubids(barcode_ht_t *ht, union_find_t *uf,
                      uint64_t **ubid_map_out) {
    uint64_t *ubid_map = xcalloc(ht->capacity, sizeof(uint64_t), "assign_ubids");
    /* First pass: assign UBID to each root */
    uint64_t next_ubid = 1; /* 0 reserved for "no match" */
    uint64_t *root_ubid = xcalloc(uf->n, sizeof(uint64_t), "assign_ubids root_ubid");

    for (uint64_t i = 0; i < ht->capacity; i++) {
        if (!ht->slots[i].occupied) continue;
        uint32_t root = uf_find(uf, ht->slots[i].uf_idx);
        if (root_ubid[root] == 0)
            root_ubid[root] = next_ubid++;
        ubid_map[i] = root_ubid[root];
    }
    free(root_ubid);
    *ubid_map_out = ubid_map;
    return next_ubid - 1;
}

/* ------------------------------------------------------------------ */
/*  N-containing barcode resolution                                   */
/* ------------------------------------------------------------------ */

/* For a barcode with N's, enumerate all 4^k concrete substitutions,
   look up each in ht, and return the UBID of the highest-count match.
   Returns 0 if no match found. */
uint64_t resolve_n_barcode(const barcode_t *bc, const barcode_ht_t *ht,
                           const uint64_t *ubid_map) {
    /* Collect N positions */
    int n_positions[MAX_BC_LEN];
    int n_count = 0;
    for (int i = 0; i < bc->len; i++) {
        uint64_t mask = (i < 32) ? bc->n_mask_lo : bc->n_mask_hi;
        int shift = (i < 32) ? i : (i - 32);
        if (mask & (1ULL << shift))
            n_positions[n_count++] = i;
    }
    if (n_count == 0) return 0;

    /* Enumerate all 4^k combinations */
    int combos = 1;
    for (int i = 0; i < n_count; i++) combos *= 4;

    uint32_t best_count = 0;
    uint64_t best_ubid = 0;

    bc_key_t key;
    key.hi = bc->hi;
    key.lo = bc->lo;
    key.len = bc->len;

    for (int c = 0; c < combos; c++) {
        /* Set each N position to the base indicated by combo index */
        bc_key_t trial = key;
        int tmp = c;
        for (int j = 0; j < n_count; j++) {
            int base = tmp & 3;
            tmp >>= 2;
            int pos = n_positions[j];
            if (pos < 32)
                trial.lo = (trial.lo & ~(3ULL << (2 * pos))) | ((uint64_t)base << (2 * pos));
            else
                trial.hi = (trial.hi & ~(3ULL << (2 * (pos - 32)))) | ((uint64_t)base << (2 * (pos - 32)));
        }
        ht_entry_t *e = ht_lookup((barcode_ht_t *)ht, &trial);
        if (e) {
            /* Find UBID for this entry by scanning for its index */
            uint64_t slot_idx = (uint64_t)(e - ht->slots);
            uint64_t ubid = ubid_map[slot_idx];
            if (e->count > best_count) {
                best_count = e->count;
                best_ubid = ubid;
            }
        }
    }
    return best_ubid;
}

/* ------------------------------------------------------------------ */
/*  Chunk buffer                                                      */
/* ------------------------------------------------------------------ */

void chunk_init(chunk_buffer_t *cb, long initial_cap) {
    cb->recs = xmalloc(initial_cap * sizeof(sort_record_t), "chunk_init");
    cb->n = 0;
    cb->cap = initial_cap;
    cb->mem_used = 0;
}

void chunk_push(chunk_buffer_t *cb, const char *hdr, const char *seq,
                const char *plus, const char *qual, uint64_t ubid) {
    if (cb->n >= cb->cap) {
        cb->cap = cb->cap * 2;
        cb->recs = xrealloc(cb->recs, cb->cap * sizeof(sort_record_t), "chunk_push");
    }
    sort_record_t *r = &cb->recs[cb->n++];
    r->header = xstrdup(hdr, "chunk_push hdr");
    r->seq = xstrdup(seq, "chunk_push seq");
    r->plus = xstrdup(plus, "chunk_push plus");
    r->qual = xstrdup(qual, "chunk_push qual");
    r->ubid = ubid;
    cb->mem_used += strlen(hdr) + strlen(seq) + strlen(plus) + strlen(qual) + 4;
}

void chunk_clear(chunk_buffer_t *cb) {
    for (long i = 0; i < cb->n; i++) {
        free(cb->recs[i].header);
        free(cb->recs[i].seq);
        free(cb->recs[i].plus);
        free(cb->recs[i].qual);
    }
    cb->n = 0;
    cb->mem_used = 0;
}

void chunk_free(chunk_buffer_t *cb) {
    chunk_clear(cb);
    free(cb->recs);
}

static int cmp_sort_record(const void *a, const void *b) {
    const sort_record_t *ra = (const sort_record_t *)a;
    const sort_record_t *rb = (const sort_record_t *)b;
    if (ra->ubid < rb->ubid) return -1;
    if (ra->ubid > rb->ubid) return 1;
    return 0;
}

/* Write sorted chunk to temp file. Returns path (caller frees). */
char *write_chunk(chunk_buffer_t *cb, int chunk_num, const char *tmp_dir) {
    qsort(cb->recs, cb->n, sizeof(sort_record_t), cmp_sort_record);

    int pathlen = (int)strlen(tmp_dir) + 64;
    char *path = xmalloc(pathlen, "write_chunk path");
    snprintf(path, pathlen, "%s/chunk_%05d.fastq", tmp_dir, chunk_num);

    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "Fatal: cannot open chunk file: %s\n", path);
        exit(EXIT_FAILURE);
    }

    for (long i = 0; i < cb->n; i++) {
        sort_record_t *r = &cb->recs[i];
        fprintf(fp, "%s\n%s\n%s\n%s\n", r->header, r->seq, r->plus, r->qual);
    }

    fclose(fp);
    fprintf(stderr, "  Chunk %d written (%ld reads)\n", chunk_num + 1, cb->n);
    return path;
}

/* ------------------------------------------------------------------ */
/*  K-way merge (min-heap)                                            */
/* ------------------------------------------------------------------ */

static uint64_t parse_ubid_from_header(const char *hdr) {
    const char *tag = strstr(hdr, ";UBID=");
    if (!tag) return 0;
    return strtoull(tag + 6, NULL, 10);
}

static int stream_advance(merge_stream_t *ms) {
    if (!fgets(ms->hdr, MAX_LINE, ms->fp))  { ms->exhausted = 1; return 0; }
    if (!fgets(ms->seq, MAX_LINE, ms->fp))   { ms->exhausted = 1; return 0; }
    if (!fgets(ms->plus_line, MAX_LINE, ms->fp)) { ms->exhausted = 1; return 0; }
    if (!fgets(ms->qual, MAX_LINE, ms->fp))  { ms->exhausted = 1; return 0; }
    /* Strip newlines */
    int l;
    l = (int)strlen(ms->hdr); if (l > 0 && ms->hdr[l-1] == '\n') ms->hdr[--l] = '\0';
    l = (int)strlen(ms->seq); if (l > 0 && ms->seq[l-1] == '\n') ms->seq[--l] = '\0';
    l = (int)strlen(ms->plus_line); if (l > 0 && ms->plus_line[l-1] == '\n') ms->plus_line[--l] = '\0';
    l = (int)strlen(ms->qual); if (l > 0 && ms->qual[l-1] == '\n') ms->qual[--l] = '\0';
    ms->ubid = parse_ubid_from_header(ms->hdr);
    return 1;
}

/* Min-heap operations */
static void heap_init(min_heap_t *h, int cap) {
    h->arr = xmalloc(cap * sizeof(merge_stream_t *), "heap_init");
    h->n = 0;
    h->cap = cap;
}

static void heap_free(min_heap_t *h) { free(h->arr); }

static void heap_siftup(min_heap_t *h, int i) {
    while (i > 0) {
        int parent = (i - 1) / 2;
        if (h->arr[parent]->ubid <= h->arr[i]->ubid) break;
        merge_stream_t *t = h->arr[parent];
        h->arr[parent] = h->arr[i];
        h->arr[i] = t;
        i = parent;
    }
}

static void heap_siftdown(min_heap_t *h, int i) {
    while (1) {
        int smallest = i;
        int l = 2 * i + 1;
        int r = 2 * i + 2;
        if (l < h->n && h->arr[l]->ubid < h->arr[smallest]->ubid) smallest = l;
        if (r < h->n && h->arr[r]->ubid < h->arr[smallest]->ubid) smallest = r;
        if (smallest == i) break;
        merge_stream_t *t = h->arr[smallest];
        h->arr[smallest] = h->arr[i];
        h->arr[i] = t;
        i = smallest;
    }
}

static void heap_push(min_heap_t *h, merge_stream_t *ms) {
    h->arr[h->n] = ms;
    heap_siftup(h, h->n);
    h->n++;
}

static merge_stream_t *heap_pop(min_heap_t *h) {
    merge_stream_t *top = h->arr[0];
    h->n--;
    if (h->n > 0) {
        h->arr[0] = h->arr[h->n];
        heap_siftdown(h, 0);
    }
    return top;
}

/* Merge a batch of chunk files into one output file */
static void merge_batch(char **paths, int n_paths, const char *out_path) {
    FILE *out = fopen(out_path, "w");
    if (!out) {
        fprintf(stderr, "Fatal: cannot open output file: %s\n", out_path);
        exit(EXIT_FAILURE);
    }

    min_heap_t heap;
    heap_init(&heap, n_paths);

    merge_stream_t *streams = xcalloc(n_paths, sizeof(merge_stream_t), "merge_batch");
    for (int i = 0; i < n_paths; i++) {
        streams[i].fp = fopen(paths[i], "r");
        if (!streams[i].fp) {
            fprintf(stderr, "Fatal: cannot open chunk file: %s\n", paths[i]);
            exit(EXIT_FAILURE);
        }
        streams[i].exhausted = 0;
        if (stream_advance(&streams[i]))
            heap_push(&heap, &streams[i]);
    }

    long merged = 0;
    while (heap.n > 0) {
        merge_stream_t *ms = heap_pop(&heap);
        fprintf(out, "%s\n%s\n%s\n%s\n", ms->hdr, ms->seq, ms->plus_line, ms->qual);
        if (stream_advance(ms))
            heap_push(&heap, ms);
        if (++merged % 5000000 == 0)
            fprintf(stderr, "  Merged %ld reads...\n", merged);
    }

    for (int i = 0; i < n_paths; i++)
        fclose(streams[i].fp);
    free(streams);
    heap_free(&heap);
    fclose(out);
}

/* Multi-pass merge if too many chunks */
void merge_sorted_chunks(char **paths, int n_chunks,
                         const char *out_path, int max_open_files,
                         const char *tmp_dir) {
    if (n_chunks == 0) {
        /* No chunks: create empty output */
        FILE *f = fopen(out_path, "w");
        if (f) fclose(f);
        return;
    }

    /* Working copies of paths */
    char **working = xmalloc(n_chunks * sizeof(char *), "merge working");
    for (int i = 0; i < n_chunks; i++)
        working[i] = xstrdup(paths[i], "merge working path");
    int n_working = n_chunks;

    int pass = 0;
    int inter_id = 0;
    char **intermediate = NULL;
    int n_intermediate = 0;

    while (n_working > max_open_files) {
        fprintf(stderr, "Merge pass %d: reducing %d chunks...\n", pass + 1, n_working);
        int n_next = 0;
        /* Count how many output files this pass will produce */
        for (int i = 0; i < n_working; i += max_open_files) n_next++;
        char **next = xmalloc(n_next * sizeof(char *), "merge next");
        int ni = 0;

        for (int i = 0; i < n_working; i += max_open_files) {
            int batch_end = i + max_open_files;
            if (batch_end > n_working) batch_end = n_working;
            int batch_n = batch_end - i;

            int pathlen = (int)strlen(tmp_dir) + 64;
            char *out = xmalloc(pathlen, "merge inter path");
            snprintf(out, pathlen, "%s/pass_%d_chunk_%05d.fastq", tmp_dir, pass, inter_id++);

            merge_batch(working + i, batch_n, out);
            next[ni++] = out;

            /* Track intermediate for cleanup */
            intermediate = xrealloc(intermediate, (n_intermediate + 1) * sizeof(char *), "merge inter");
            intermediate[n_intermediate++] = xstrdup(out, "merge inter track");
        }

        for (int i = 0; i < n_working; i++) free(working[i]);
        free(working);
        working = next;
        n_working = ni;
        pass++;
    }

    fprintf(stderr, "Final merge of %d chunks into output...\n", n_working);
    merge_batch(working, n_working, out_path);

    /* Cleanup */
    for (int i = 0; i < n_working; i++) free(working[i]);
    free(working);
    for (int i = 0; i < n_intermediate; i++) {
        remove(intermediate[i]);
        free(intermediate[i]);
    }
    free(intermediate);
}

/* ------------------------------------------------------------------ */
/*  Callback-based merge (for combined pipeline)                      */
/* ------------------------------------------------------------------ */

/* Flush accumulated group to callback */
static void flush_group(merge_group_cb_t cb, void *user_data,
                        char *ubid_str, char **seqs, char **quals,
                        int count) {
    if (count > 0 && cb) {
        cb(ubid_str, seqs, quals, count, user_data);
    } else {
        /* No callback or empty group: free the arrays */
        for (int i = 0; i < count; i++) {
            free(seqs[i]);
            free(quals[i]);
        }
    }
}

/* Merge a batch of chunk files, delivering UBID groups via callback */
static void merge_batch_with_callback(char **paths, int n_paths,
                                      merge_group_cb_t cb, void *user_data) {
    min_heap_t heap;
    heap_init(&heap, n_paths);

    merge_stream_t *streams = xcalloc(n_paths, sizeof(merge_stream_t), "merge_batch_cb");
    for (int i = 0; i < n_paths; i++) {
        streams[i].fp = fopen(paths[i], "r");
        if (!streams[i].fp) {
            fprintf(stderr, "Fatal: cannot open chunk file: %s\n", paths[i]);
            exit(EXIT_FAILURE);
        }
        streams[i].exhausted = 0;
        if (stream_advance(&streams[i]))
            heap_push(&heap, &streams[i]);
    }

    /* Group accumulation state */
    int grp_cap = 256;
    char **grp_seqs = xmalloc(grp_cap * sizeof(char *), "merge_cb grp_seqs");
    char **grp_quals = xmalloc(grp_cap * sizeof(char *), "merge_cb grp_quals");
    int grp_count = 0;
    uint64_t cur_ubid = 0;
    char cur_ubid_str[32] = {0};

    long merged = 0;
    while (heap.n > 0) {
        merge_stream_t *ms = heap_pop(&heap);

        /* If UBID changed, flush the previous group */
        if (ms->ubid != cur_ubid && grp_count > 0) {
            flush_group(cb, user_data, cur_ubid_str, grp_seqs, grp_quals, grp_count);
            /* Allocate fresh arrays (old ones owned by callback) */
            grp_seqs = xmalloc(grp_cap * sizeof(char *), "merge_cb grp_seqs");
            grp_quals = xmalloc(grp_cap * sizeof(char *), "merge_cb grp_quals");
            grp_count = 0;
        }

        cur_ubid = ms->ubid;
        snprintf(cur_ubid_str, sizeof(cur_ubid_str), "%lu", (unsigned long)ms->ubid);

        /* Grow arrays if needed */
        if (grp_count >= grp_cap) {
            grp_cap *= 2;
            grp_seqs = xrealloc(grp_seqs, grp_cap * sizeof(char *), "merge_cb grp_seqs grow");
            grp_quals = xrealloc(grp_quals, grp_cap * sizeof(char *), "merge_cb grp_quals grow");
        }

        /* Strip trailing whitespace from seq/qual and strdup */
        grp_seqs[grp_count] = xstrdup(ms->seq, "merge_cb seq");
        grp_quals[grp_count] = xstrdup(ms->qual, "merge_cb qual");
        grp_count++;

        if (stream_advance(ms))
            heap_push(&heap, ms);

        if (++merged % 5000000 == 0)
            fprintf(stderr, "  Merged %ld reads...\n", merged);
    }

    /* Flush final group */
    if (grp_count > 0) {
        flush_group(cb, user_data, cur_ubid_str, grp_seqs, grp_quals, grp_count);
    } else {
        free(grp_seqs);
        free(grp_quals);
    }

    for (int i = 0; i < n_paths; i++)
        fclose(streams[i].fp);
    free(streams);
    heap_free(&heap);
}

/* Multi-pass merge with callback for final pass */
void merge_sorted_chunks_with_callback(
    char **paths, int n_chunks,
    int max_open_files, const char *tmp_dir,
    merge_group_cb_t cb, void *user_data) {

    if (n_chunks == 0) return;

    /* Working copies of paths */
    char **working = xmalloc(n_chunks * sizeof(char *), "merge_cb working");
    for (int i = 0; i < n_chunks; i++)
        working[i] = xstrdup(paths[i], "merge_cb working path");
    int n_working = n_chunks;

    int pass = 0;
    int inter_id = 0;
    char **intermediate = NULL;
    int n_intermediate = 0;

    /* Intermediate passes: use file-based merge_batch (same as before) */
    while (n_working > max_open_files) {
        fprintf(stderr, "Merge pass %d: reducing %d chunks...\n", pass + 1, n_working);
        int n_next = 0;
        for (int i = 0; i < n_working; i += max_open_files) n_next++;
        char **next = xmalloc(n_next * sizeof(char *), "merge_cb next");
        int ni = 0;

        for (int i = 0; i < n_working; i += max_open_files) {
            int batch_end = i + max_open_files;
            if (batch_end > n_working) batch_end = n_working;
            int batch_n = batch_end - i;

            int pathlen = (int)strlen(tmp_dir) + 64;
            char *out = xmalloc(pathlen, "merge_cb inter path");
            snprintf(out, pathlen, "%s/pass_%d_chunk_%05d.fastq", tmp_dir, pass, inter_id++);

            merge_batch(working + i, batch_n, out);
            next[ni++] = out;

            intermediate = xrealloc(intermediate, (n_intermediate + 1) * sizeof(char *), "merge_cb inter");
            intermediate[n_intermediate++] = xstrdup(out, "merge_cb inter track");
        }

        for (int i = 0; i < n_working; i++) free(working[i]);
        free(working);
        working = next;
        n_working = ni;
        pass++;
    }

    /* Final pass: callback-based merge */
    fprintf(stderr, "Final merge of %d chunks (with consensus)...\n", n_working);
    merge_batch_with_callback(working, n_working, cb, user_data);

    /* Cleanup */
    for (int i = 0; i < n_working; i++) free(working[i]);
    free(working);
    for (int i = 0; i < n_intermediate; i++) {
        remove(intermediate[i]);
        free(intermediate[i]);
    }
    free(intermediate);
}
