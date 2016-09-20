/*
 *
 * Copyright (c) 2011, Jue Ruan <ruanjue@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __DBGCNS_CNS_RJ_H
#define __DBGCNS_CNS_RJ_H

#include "dna.h"
#include "list.h"
#include "hashset.h"
#include "thread.h"
#include "ksw.h"
#include "string.h"

#define DBGCNS_KMER_MAX_SIZE	16
#define DBGCNS_KMER_MAX_NODE_COV	0x3FF
#define DBGCNS_KMER_MAX_EDGE_COV	0xFF

typedef struct {
	uint64_t mer:32, cov:10, extra:22;
	uint8_t  edges[2][4];
} dbgcns_kmer_t;
define_list(dbgcnskmerv, dbgcns_kmer_t);
#define UD(E) ((dbgcnskmerv*)set->userdata)->buffer[E].mer
#define kmer_hashcode(E) u32hashcode(UD(E))
#define kmer_hashequals(E1, E2) UD(E1) == UD(E2)
define_hashset(dbgcnskmerhash, uint32_t, kmer_hashcode, kmer_hashequals);
#undef UD

typedef struct {
	uint32_t ksize, kmask;
	dbgcnskmerv    *kmers;
	dbgcnskmerhash *khash;
} DBG;

#define DBGCNS_DP_SCORE_MIN	-0x7FFFFFFF
#define DBGCNS_PATH_M	0
#define DBGCNS_PATH_X	1
#define DBGCNS_PATH_I	2
#define DBGCNS_PATH_D	3
#define DBGCNS_CNS_NON	0
#define DBGCNS_CNS_TIP	1
#define DBGCNS_CNS_CUT	2
#define DBGCNS_CNS_EXT	3
#define DBGCNS_CNS_HIT	4

typedef struct {
	union {
		struct { uint32_t kidx:30, path:2; uint32_t qpos; };
		uint64_t identifier;
	};
	int score;
	uint32_t bt_idx;
} dbgcns_dp_t;
define_list(dbgcnsdpv, dbgcns_dp_t);
#define DD(E) ((dbgcnsdpv*)set->userdata)->buffer[E]
#define dp_hashcode(E) u64hashcode(DD(E).identifier)
#define dp_hashequals(E1, E2) DD(E1).identifier == DD(E2).identifier
define_hashset(dbgcnsdphash, uint32_t, dp_hashcode, dp_hashequals);
#undef DD

typedef struct {uint64_t off:40, len:24;} blk_t;
define_list(blkv, blk_t);

typedef struct {
	DBG      *g;
	int      C, M, X, I, D, E;
	int      Z, W;

	u1v      *qseqs;
	blkv     *qblks;

	uint8_t  *qry;
	uint32_t qlen;

	dbgcnsdpv      *dps;
	u4v      *heap;
	dbgcnsdphash   *hash;

	b4v      *qmaxs;
	uint32_t qtop;
	int      max_score;
	uint32_t best_idx;

	String   *seq;
	u8list   *cns;
	int      alns[4];
} CNS;

static inline DBG* init_dbg(uint32_t ksize){
	DBG *g;
	if(ksize > DBGCNS_KMER_MAX_SIZE){
		fprintf(stderr, " -- ksize MUST be no greater than %d, but is %d in %s -- %s:%d --\n", DBGCNS_KMER_MAX_SIZE, ksize, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		exit(1);
	}
	g = malloc(sizeof(DBG));
	g->ksize = ksize;
	g->kmask = 0xFFFFFFFFU >> ((16 - ksize) << 1);
	g->kmers = init_dbgcnskmerv(32);
	next_ref_dbgcnskmerv(g->kmers);
	memset(g->kmers->buffer, 0, sizeof(dbgcns_kmer_t));
	g->khash = init_dbgcnskmerhash(1023);
	set_userdata_dbgcnskmerhash(g->khash, g->kmers);
	return g;
}

static inline void free_dbg(DBG *g){
	free_dbgcnskmerv(g->kmers);
	free_dbgcnskmerhash(g->khash);
	free(g);
}

static inline void clear_dbg(DBG *g){
	clear_dbgcnskmerv(g->kmers);
	next_ref_dbgcnskmerv(g->kmers);
	memset(g->kmers->buffer, 0, sizeof(dbgcns_kmer_t));
	clear_dbgcnskmerhash(g->khash);
}

static inline int add_kmer_dbg(DBG *g, uint32_t kmer, uint8_t fbase, uint8_t rbase){
	dbgcns_kmer_t *k;
	uint32_t krev, *u;
	int dir, exists;
	if(0){
		krev = dna_rev_seq(kmer, g->ksize);
		if(kmer < krev){ dir = 0; }
		else { kmer = krev; dir = 1; }
	} else {
		dir = 0;
	}
	g->kmers->buffer[0].mer = kmer;
	u = prepare_dbgcnskmerhash(g->khash, 0, &exists);
	if(exists){
		k = ref_dbgcnskmerv(g->kmers, *u);
	} else {
		*u = g->kmers->size;
		k = next_ref_dbgcnskmerv(g->kmers);
		memset(k, 0, sizeof(dbgcns_kmer_t));
		k->mer = kmer;
	}
	if(k->cov < DBGCNS_KMER_MAX_NODE_COV) k->cov ++;
	if(dir){
		if(fbase < 4){
			fbase = (~fbase) & 0x03;
			if(k->edges[1][fbase] < DBGCNS_KMER_MAX_EDGE_COV) k->edges[1][fbase] ++;
		}
		if(rbase < 4){
			if(k->edges[0][rbase] < DBGCNS_KMER_MAX_EDGE_COV) k->edges[0][rbase] ++;
		}
	} else {
		if(fbase < 4){
			if(k->edges[0][fbase] < DBGCNS_KMER_MAX_EDGE_COV) k->edges[0][fbase] ++;
		}
		if(rbase < 4){
			rbase = (~rbase) & 0x03;
			if(k->edges[1][rbase] < DBGCNS_KMER_MAX_EDGE_COV) k->edges[1][rbase] ++;
		}
	}
	return exists;
}

static void add_seq_dbg(DBG *g, uint8_t *seq, uint32_t len){
	uint32_t kmer, i;
	uint8_t b, f, r;
	kmer = 0;
	for(i=0;i<len;){
		b = seq[i];
		kmer = ((kmer << 2) | b) & g->kmask;
		i ++;
		if(i < g->ksize) continue;
		f = (i < len)? seq[i] : 4;
		r = (i > g->ksize)? seq[i - g->ksize] : 4;
		add_kmer_dbg(g, kmer, f, r);
	}
}

static inline void print_kmers_dbg(DBG *g, FILE *out){
	dbgcns_kmer_t *k;
	uint64_t i;
	char seq[DBGCNS_KMER_MAX_SIZE + 1];
	for(i=1;i<g->kmers->size;i++){
		k = ref_dbgcnskmerv(g->kmers, i);
		kmer2seq(seq, k->mer, g->ksize);
		fprintf(out, "%s\t%d\t%d,%d,%d,%d\t%d,%d,%d,%d\n", seq, k->cov,
			k->edges[0][0], k->edges[0][1], k->edges[0][2], k->edges[0][3],
			k->edges[1][0], k->edges[1][1], k->edges[1][2], k->edges[1][3]);
	}
}

static inline CNS* init_cns(uint32_t ksize, int Z, int W, int M, int X, int I, int D, int E){
	CNS *cns;
	cns = malloc(sizeof(CNS));
	cns->g = init_dbg(ksize);
	cns->qseqs = init_u1v(32);
	cns->qblks = init_blkv(32);
	cns->Z = Z;
	cns->W = W;
	cns->C = 1;
	cns->M = M;
	cns->X = X;
	cns->I = I;
	cns->D = D;
	cns->E = E;
	cns->qlen = 0;
	cns->dps  = init_dbgcnsdpv(32);
	next_ref_dbgcnsdpv(cns->dps);
	memset(cns->dps->buffer, 0, sizeof(dbgcns_dp_t)); // no need, it is always zero-filled
	cns->heap = init_u4v(32);
	cns->hash = init_dbgcnsdphash(1023);
	set_userdata_dbgcnsdphash(cns->hash, cns->dps);
	cns->qmaxs = init_b4v(32);
	cns->qtop = 0;
	cns->max_score = DBGCNS_DP_SCORE_MIN;
	cns->best_idx  = 0;
	cns->seq = init_string(32);
	cns->cns = init_u8list(32);
	return cns;
}

static inline void free_cns(CNS *cns){
	free_dbg(cns->g);
	free_u1v(cns->qseqs);
	free_blkv(cns->qblks);
	free_dbgcnsdpv(cns->dps);
	free_u4v(cns->heap);
	free_dbgcnsdphash(cns->hash);
	free_b4v(cns->qmaxs);
	free_string(cns->seq);
	free_u8list(cns->cns);
	free(cns);
}

static inline void reset_cns(CNS *cns){
	clear_dbg(cns->g);
	clear_u1v(cns->qseqs);
	clear_blkv(cns->qblks);
	cns->qry = NULL;
	cns->qlen = 0;
}

static inline void add_seq_cns(CNS *cns, char *seq, int len){
	blk_t *b;
	int i;
	b = next_ref_blkv(cns->qblks);
	b->off = cns->qseqs->size;
	b->len = len;
	for(i=0;i<len;i++) push_u1v(cns->qseqs, base_bit_table[(int)seq[i]]);
}

static inline void ready_cns(CNS *cns){
	UNUSED(cns);
	//blk_t *b;
	//u4i i;
	//for(i=0;i<cns->qblks->size;i++){
		//b = ref_blkv(cns->qblks, i);
		//add_seq_dbg(cns->g, cns->qseqs->buffer + b->off, b->len);
	//}
}

static inline int dbg_cns_core2(CNS *cns){
	dbgcns_dp_t *dp, *dp2;
	dbgcns_kmer_t *k;
	uint32_t i, dp_idx, kmer, knew, kidx, *u;
	int sum, cov, cut;
	uint8_t b, q;
	int exists, score;
	if(cns->heap->size == 0) return DBGCNS_CNS_NON;
	dp_idx = cns->heap->buffer[0]; //array_heap_pop(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
	encap_dbgcnsdpv(cns->dps, 9);
	dp = ref_dbgcnsdpv(cns->dps, dp_idx);
	if(dp->qpos >= cns->qlen) return DBGCNS_CNS_HIT;
	if(dp->qpos + cns->W < cns->qtop){
		array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, 0, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
		return DBGCNS_CNS_CUT;
	} else if(dp->qpos > cns->qtop){
		cns->qtop = dp->qpos;
		for(i=cns->heap->size;i>0;i--){
			if(cns->dps->buffer[cns->heap->buffer[i-1]].qpos + cns->W < cns->qtop){
				array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, i - 1, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
			}
		}
	}
	if(dp->score < cns->qmaxs->buffer[dp->qpos]){
		array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, 0, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
		return DBGCNS_CNS_CUT;
	} else if(dp->score + cns->Z * cns->X > cns->qmaxs->buffer[dp->qpos]){
		cns->qmaxs->buffer[dp->qpos] = dp->score + cns->Z * cns->X;
	}
	u = prepare_dbgcnsdphash(cns->hash, dp_idx, &exists);
	if(exists){
		dp2 = ref_dbgcnsdpv(cns->dps, *u);
		if(dp->score > dp2->score){
			*u = dp_idx;
		} else {
			array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, 0, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
			return DBGCNS_CNS_CUT;
		}
	} else {
		*u = dp_idx;
	}
	array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, 0, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
	k = ref_dbgcnskmerv(cns->g->kmers, dp->kidx);
	kmer = k->mer;
	q = cns->qry[dp->qpos];
	sum = k->edges[0][0] + k->edges[0][1] + k->edges[0][2] + k->edges[0][3];
	if(sum == 0) return DBGCNS_CNS_TIP;
	cut = num_max(sum / 3, 3);
	for(b=0;b<4;b++){
		if((cov = k->edges[0][b]) == 0) continue;
		if(b == q){
			score = dp->score + cns->M;
		} else {
			score = dp->score + cns->X;
		}
		score += (cov > 1)? (cov >= cut? 1 : 0) : -1;
		knew = ((kmer << 2) | b) & cns->g->kmask;
		cns->g->kmers->buffer[0].mer = knew;
		u = get_dbgcnskmerhash(cns->g->khash, 0);
		kidx = *u;
		dp2 = next_ref_dbgcnsdpv(cns->dps);
		dp2->kidx = kidx;
		dp2->path = (b == q)? DBGCNS_PATH_M : DBGCNS_PATH_X;
		dp2->qpos = dp->qpos + 1;
		dp2->score = score;
		dp2->bt_idx = dp_idx;
		array_heap_push(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, cns->dps->size - 1, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
		// deletion
		if(dp->path != DBGCNS_PATH_I){
			dp2 = next_ref_dbgcnsdpv(cns->dps);
			score = dp->score + (cns->E + (dp->path == DBGCNS_PATH_D? 0 : cns->D));
			score += (cov > 1)? (cov >= cut? 1 : 0) : -1;
			dp2->kidx = kidx;
			dp2->path = DBGCNS_PATH_D;
			dp2->qpos = dp->qpos;
			dp2->score = score;
			dp2->bt_idx = dp_idx;
			array_heap_push(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, cns->dps->size - 1, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
		}
	}
	// insertion
	if(dp->path != DBGCNS_PATH_D){
		dp2 = next_ref_dbgcnsdpv(cns->dps);
		score = dp->score + (cns->E + (dp->path == DBGCNS_PATH_I? 0 : cns->I));
		dp2->kidx = dp->kidx;
		dp2->path = DBGCNS_PATH_I;
		dp2->qpos = dp->qpos + 1;
		dp2->score = score;
		dp2->bt_idx = dp_idx;
		array_heap_push(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, cns->dps->size - 1, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
	}
	return DBGCNS_CNS_EXT;
}

static inline int dbg_cns_core(CNS *cns){
	dbgcns_dp_t *dp, *dp2;
	dbgcns_kmer_t *k;
	uint32_t i, dp_idx, kmer, knew, kidx, *u;
	int sum, cov, cut;
	uint8_t b, q;
	int exists, score, nadd, mat_only;
	if(cns->heap->size == 0) return DBGCNS_CNS_NON;
	dp_idx = cns->heap->buffer[0]; //array_heap_pop(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
	encap_dbgcnsdpv(cns->dps, 9);
	dp = ref_dbgcnsdpv(cns->dps, dp_idx);
	mat_only = 0;
	if(dp->qpos >= cns->qlen) return DBGCNS_CNS_HIT;
	if(dp->qpos + cns->W < cns->qtop){
		mat_only = 1;
		//array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, 0, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
		//return DBGCNS_CNS_CUT;
	} else if(dp->qpos > cns->qtop){
		cns->qtop = dp->qpos;
		for(i=cns->heap->size;i>0;i--){
			if(cns->dps->buffer[cns->heap->buffer[i-1]].qpos + cns->W < cns->qtop){
				array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, i - 1, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
			}
		}
	}
	if(dp->score < cns->qmaxs->buffer[dp->qpos]){
		mat_only = 1;
		//array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, 0, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
		//return DBGCNS_CNS_CUT;
	} else if(dp->score + cns->Z * cns->X > cns->qmaxs->buffer[dp->qpos]){
		cns->qmaxs->buffer[dp->qpos] = dp->score + cns->Z * cns->X;
	}
	u = prepare_dbgcnsdphash(cns->hash, dp_idx, &exists);
	if(exists){
		dp2 = ref_dbgcnsdpv(cns->dps, *u);
		if(dp->score > dp2->score){
			*u = dp_idx;
		} else {
			array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, 0, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
			return DBGCNS_CNS_CUT;
		}
	} else {
		*u = dp_idx;
	}
	array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, 0, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
	k = ref_dbgcnskmerv(cns->g->kmers, dp->kidx);
	kmer = k->mer;
	q = cns->qry[dp->qpos];
	sum = k->edges[0][0] + k->edges[0][1] + k->edges[0][2] + k->edges[0][3];
	if(sum == 0) return DBGCNS_CNS_TIP;
	cut = num_max(sum / 3, 3);
	nadd = 0;
	for(b=0;b<4;b++){
		if((cov = k->edges[0][b]) == 0) continue;
		if(b == q){
			score = dp->score + cns->M;
		} else if(mat_only){
			continue;
		} else {
			score = dp->score + cns->X;
		}
		score += (cov > 1)? (cov >= cut? 1 : 0) : -1;
		knew = ((kmer << 2) | b) & cns->g->kmask;
		cns->g->kmers->buffer[0].mer = knew;
		u = get_dbgcnskmerhash(cns->g->khash, 0);
		kidx = *u;
		dp2 = next_ref_dbgcnsdpv(cns->dps);
		dp2->kidx = kidx;
		dp2->path = (b == q)? DBGCNS_PATH_M : DBGCNS_PATH_X;
		dp2->qpos = dp->qpos + 1;
		dp2->score = score;
		dp2->bt_idx = dp_idx;
		array_heap_push(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, cns->dps->size - 1, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
		nadd ++;
		// deletion
		if(dp->path != DBGCNS_PATH_I){
			dp2 = next_ref_dbgcnsdpv(cns->dps);
			score = dp->score + (cns->E + (dp->path == DBGCNS_PATH_D? 0 : cns->D));
			score += (cov > 1)? (cov >= cut? 1 : 0) : -1;
			dp2->kidx = kidx;
			dp2->path = DBGCNS_PATH_D;
			dp2->qpos = dp->qpos;
			dp2->score = score;
			dp2->bt_idx = dp_idx;
			array_heap_push(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, cns->dps->size - 1, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
			nadd ++;
		}
	}
	// insertion
	if(mat_only == 0 && dp->path != DBGCNS_PATH_D){
		dp2 = next_ref_dbgcnsdpv(cns->dps);
		score = dp->score + (cns->E + (dp->path == DBGCNS_PATH_I? 0 : cns->I));
		dp2->kidx = dp->kidx;
		dp2->path = DBGCNS_PATH_I;
		dp2->qpos = dp->qpos + 1;
		dp2->score = score;
		dp2->bt_idx = dp_idx;
		array_heap_push(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, cns->dps->size - 1, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
		nadd ++;
	}
	return nadd? DBGCNS_CNS_EXT : DBGCNS_CNS_CUT;
}

static inline void ready_core_cns(CNS *cns, int candidate_mode){
	blk_t *b;
	u4i i;
	if(candidate_mode == 3) i = 1;
	else i = 0;
	for(;i<cns->qblks->size;i++){
		b = ref_blkv(cns->qblks, i);
		add_seq_dbg(cns->g, cns->qseqs->buffer + b->off, b->len);
	}
}

static inline int run_core_cns(CNS *cns, uint8_t *qry, uint32_t qlen){
	dbgcns_dp_t *dp;
	dbgcns_kmer_t *k;
	uint32_t i, kmer, kidx, dp_idx, *u;
	int status;
	{
		cns->qry  = qry;
		cns->qlen = qlen;
		// reset auxiliaries
		clear_dbgcnsdpv(cns->dps);
		next_ref_dbgcnsdpv(cns->dps);
		memset(cns->dps->buffer, 0, sizeof(dbgcns_dp_t));
		clear_u4v(cns->heap);
		clear_dbgcnsdphash(cns->hash);
		clear_b4v(cns->qmaxs);
		for(i=0;i<qlen;i++) push_b4v(cns->qmaxs, DBGCNS_DP_SCORE_MIN);
		cns->qtop = 0;
		cns->max_score= DBGCNS_DP_SCORE_MIN;
		cns->best_idx = 0;
		clear_u8list(cns->cns);
		clear_string(cns->seq);
		cns->alns[0] = cns->alns[1] = cns->alns[2] = cns->alns[3] = 0;
	}
	// set first kmer
	kmer = 0;
	for(cns->qtop=0;cns->qtop<cns->g->ksize;cns->qtop++){
		kmer = (kmer << 2) | qry[cns->qtop];
	}
	cns->g->kmers->buffer[0].mer = kmer;
	u = get_dbgcnskmerhash(cns->g->khash, 0);
	if(u == NULL) return 0;
	kidx = *u;
	dp = next_ref_dbgcnsdpv(cns->dps);
	dp->kidx = kidx;
	dp->path = DBGCNS_PATH_M;
	dp->qpos = cns->qtop;
	dp->score = 0;
	dp->bt_idx = 0;
	array_heap_push(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, cns->dps->size - 1, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
	// dbg traversing
	while(cns->heap->size){
		status = dbg_cns_core(cns);
		if(status == DBGCNS_CNS_HIT){
			dp_idx = cns->heap->buffer[0];
			array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, 0, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
			dp = ref_dbgcnsdpv(cns->dps, dp_idx);
			if(dp->score > cns->max_score){
				cns->max_score = dp->score;
				cns->best_idx = dp_idx;
			}
		}
	}
	if(cns->best_idx == 0) return 0;
	// traceback to get cns seq
	dp_idx = cns->best_idx;
	while(dp_idx){
		dp = ref_dbgcnsdpv(cns->dps, dp_idx);
		cns->alns[dp->path] ++;
		if(dp->path != DBGCNS_PATH_I){
			k = ref_dbgcnskmerv(cns->g->kmers, dp->kidx);
			push_u8list(cns->cns, k->mer & 0x03);
			add_char_string(cns->seq, bit_base_table[k->mer & 0x03]);
		}
		dp_idx = dp->bt_idx;
	}
	// first ksize - 1 bases may be not corrected, truncated
	reverse_string(cns->seq);
	reverse_u8list(cns->cns);
	return cns->seq->size;
}

static inline int run_cns(CNS *cns, int candidate_mode){
	blk_t ref;
	u4i i;
	if(candidate_mode == 1){ // longest
		ref.off = ref.len = 0;
		for(i=0;i<cns->qblks->size;i++){
			if(cns->qblks->buffer[i].len > ref.len) ref = cns->qblks->buffer[i];
		}
	} else if(candidate_mode == 2){ // shortest
		ref.off = 0; ref.len = cns->qseqs->size;
		for(i=0;i<cns->qblks->size;i++){
			if(cns->qblks->buffer[i].len < ref.len) ref = cns->qblks->buffer[i];
		}
	} else if(candidate_mode == 3){ // first and exclusive
		ref = cns->qblks->buffer[0];
	} else { // median
		ref = quick_median_array(cns->qblks->buffer, cns->qblks->size, blk_t, num_cmpgt(a.len, b.len));
	}
	ready_core_cns(cns, candidate_mode);
	return run_core_cns(cns, cns->qseqs->buffer + ref.off, ref.len);
}

#endif
