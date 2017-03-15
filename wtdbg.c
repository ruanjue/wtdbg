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

#include "dmo.h"

#define WT_MAX_RD			0x03FFFFFF
#define WT_MAX_RDLEN		0x0001FFFF
#define WT_MAX_NODE			0x000000FFFFFFFFFFLLU
#define WT_MAX_EDGE			0x000000FFFFFFFFFFLLU
#define WT_MAX_NODE_EDGES	0xFFFF
#define WT_MAX_EDGE_COV		0x7FFF

typedef struct {
	uint64_t node;
	uint64_t rid:26, dir:1, beg:18, end:18, closed:1;
} rd_reg_t;
define_list(rdregv, rd_reg_t);

typedef struct {
	uint64_t rid:26, dir:1, beg:18, end:18, closed:1;
	uint64_t read_link;
} rd_rep_t;
define_list(rdrepv, rd_rep_t);

typedef struct {
	uint64_t node;
	uint64_t rid:26, dir:1, beg:18, end:18, closed:1;
	uint64_t read_link;
} reg_t;
define_list(regv, reg_t);

typedef struct { uint32_t x, y; } xy_t;
define_list(xyv, xy_t);

typedef struct { uint64_t idx:46, cnt:18; } vec_ref_t;
typedef struct { uint64_t idx:46, cnt:18; } ptr_ref_t;
typedef struct { uint64_t idx:46, cnt:18, fix:1, rank:63; } rnk_ref_t;
static const vec_ref_t VEC_REF_NULL = (vec_ref_t){0, 0};
static const ptr_ref_t PTR_REF_NULL = (ptr_ref_t){0, 0};
static const rnk_ref_t RNK_REF_NULL = (rnk_ref_t){0, 0, 0, 0};
define_list(vecrefv, vec_ref_t);
define_list(ptrrefv, ptr_ref_t);
define_list(rnkrefv, rnk_ref_t);
#define ptrref_hashcode(E) u64hashcode((E).idx)
#define ptrref_hashequals(E1, E2) ((E1).idx == (E2).idx)
define_hashset(ptrrefhash, ptr_ref_t, ptrref_hashcode, ptrref_hashequals);

#define WT_EDGE_CLOSED_NULL	0
#define WT_EDGE_CLOSED_MASK	1
#define WT_EDGE_CLOSED_LESS	2
#define WT_EDGE_CLOSED_HARD	3
typedef struct {
	uint64_t node1:45, dir1:1, dir2:1, status:1, closed:2, cov:14;
	uint64_t node2:45; int64_t off:19;
} edge_t;
define_list(edgev, edge_t);

static inline uint64_t _edge_hashcode(edge_t e){
	const uint64_t m = 0xc6a4a7935bd1e995LLU;
	const int r = 47;
	uint64_t h = 1023 ^ (16 * m);
	uint64_t k = (e.node1 << 1) | e.dir1;
	k *= m;
	k ^= k >> r;
	k *= m;
	h ^= k;
	h *= m;
	k = (e.node2 << 1) | e.dir2;
	k *= m;
	k ^= k >> r;
	k *= m;
	h ^= k;
	h *= m;
	h ^= h >> r;
	h *= m;
	h ^= h >> r;
	return h;
}
#define EDGEHASH(idx) ((edgev *)set->userdata)->buffer[idx]
#define edge_hashcode(E) _edge_hashcode(EDGEHASH(E))
#define edge_hashequals(E1, E2) (EDGEHASH(E1).node1 == EDGEHASH(E2).node1 && EDGEHASH(E1).node2 == EDGEHASH(E2).node2 \
			&& EDGEHASH(E1).dir1 == EDGEHASH(E2).dir1 && EDGEHASH(E1).dir2 == EDGEHASH(E2).dir2)
define_hashset(edgehash, uint64_t, edge_hashcode, edge_hashequals);

typedef struct { uint64_t idx:63, flg:1; uint64_t next; } edge_ref_t;
static const edge_ref_t EDGE_REF_NULL = (edge_ref_t){0x7FFFFFFFFFFFFFFLLU, 1, 0};
define_list(edgerefv, edge_ref_t);

#define MAX_REP_IDX	0xFFFFFFFFFFFFLLU

typedef struct {
	unsigned long long rep_idx:48, unvisit:16;
	uint64_t closed:1, single_in:1, bt_visit:45, rep_dir:1, bt_idx:16, init_end:1;
	vec_ref_t regs;
	ptr_ref_t edges[2];
} node_t;
define_list(nodev, node_t);

typedef struct {
	uint32_t count:31, flag:1;
	ptr_ref_t regs;
} read_t;
define_list(readv, read_t);

typedef struct {
	u8i node;
	//node_t *n;
	edge_ref_t edges[2];
	u4i dir:2, cov:30;
	int off;
} trace_t;
define_list(tracev, trace_t);
#define WT_TRACE_MSG_ZERO	0
#define WT_TRACE_MSG_ONE	1
#define WT_TRACE_MSG_MORE	2
#define WT_TRACE_MSG_VISITED	3
#define WT_TRACE_MSG_UNDEF	4

typedef struct {
	union {
		struct {
			u4i frg1:31, dir1:1;
			u4i frg2:31, dir2:1;
		};
		u8i key;
	};
	u4i cov:16, flag:1, tidx:12, weak:1, closed:2;
	b4i off;
} lnk_t;
define_list(lnkv, lnk_t);
#define lnk_hashcode(E) u64hashcode(E.key)
#define lnk_hashequals(E1, E2) (E1).key == (E2).key
define_hashset(lnkhash, lnk_t, lnk_hashcode, lnk_hashequals);

typedef struct {
	u8i toff:46, tcnt:18; // extended traces and core traces
	u4i tx, ty; // core traces
	ptr_ref_t lnks[2];
	u4i len, length;
	unsigned long long rep_idx:48, unvisit:16;
	uint64_t closed:1, single_in:1, bt_visit:46, rep_dir:1, bt_idx:16;
} frg_t;
define_list(frgv, frg_t);

typedef struct {
	edge_ref_t lnks[2];
	u4i frg;
	u4i dir:2, cov:30;
	int off;
	u4i tx, ty;
} path_t;
define_list(pathv, path_t);
#define WT_PATH_MSG_ZERO	0
#define WT_PATH_MSG_ONE	1
#define WT_PATH_MSG_MORE	2
#define WT_PATH_MSG_VISITED	3
#define WT_PATH_MSG_UNDEF	4

typedef struct {
	DMO      *dmo;
	regv     *regs;
	readv    *reads;
	nodev    *nodes;
	edgev    *edges;
	edgehash *ehash;
	edgerefv *erefs;

	frgv     *frgs;
	lnkv     *lnks;
	edgerefv *lrefs;
	tracev   *traces;

	int      node_order;
	uint32_t n_fix, only_fix; // first n sequences are accurate contigs; only_fix means whether to include other pacbio sequenes
	uint32_t reglen, regoff, regovl;
	float    node_merge_cutoff;
	uint32_t max_node_cov, min_node_cov, min_edge_cov;
	u4i      max_node_cov_sg, max_sg_end;
	int      cut_tip;
	int      store_low_cov_edge;
	int      rep_filter, rep_detach;
	uint32_t bub_step, tip_step, rep_step;
	int min_ctg_len, min_ctg_nds;

	vplist *utgs;
	u4i    major_nctg;
	vplist *ctgs;
} Graph;

static const char *colors[2][2] = {{"blue", "green"}, {"red", "gray"}};

Graph* init_graph(DMO *dmo){
	Graph *g;
	g = malloc(sizeof(Graph));
	g->dmo = dmo;
	g->regs = init_regv(32);
	g->reads = init_readv(dmo->n_rd);
	g->reads->size = dmo->n_rd;
	g->nodes = init_nodev(32);
	g->edges = init_edgev(32);
	g->ehash = init_edgehash(1023);
	set_userdata_edgehash(g->ehash, g->edges);
	g->erefs = init_edgerefv(32);
	g->frgs = init_frgv(32);
	g->lnks = init_lnkv(32);
	g->lrefs = init_edgerefv(32);
	g->traces = init_tracev(32);
	g->node_order = 0;
	g->n_fix = 0;
	g->only_fix = 0;
	g->reglen = 1000;
	g->regoff = 1000;
	g->regovl = 200;
	g->node_merge_cutoff = 0.8;
	g->min_node_cov = 3;
	g->max_node_cov = 60;
	g->min_edge_cov = 3;
	g->max_node_cov_sg = 2;
	g->max_sg_end = 5;
	g->store_low_cov_edge = 1;
	g->rep_filter = 1;
	g->rep_detach = 1;
	g->bub_step = 10;
	g->tip_step = 5;
	g->rep_step = 20;
	g->cut_tip = 1;
	g->min_ctg_len = 10000;
	g->min_ctg_nds = 5;
	g->utgs = init_vplist(32);
	g->ctgs = init_vplist(32);
	g->major_nctg = 0;
	return g;
}

void free_graph(Graph *g){
	uint64_t i;
	free_regv(g->regs);
	free_readv(g->reads);
	free_nodev(g->nodes);
	free_edgev(g->edges);
	free_edgehash(g->ehash);
	free_edgerefv(g->erefs);
	free_frgv(g->frgs);
	free_lnkv(g->lnks);
	free_edgerefv(g->lrefs);
	free_tracev(g->traces);
	for(i=0;i<g->utgs->size;i++) free_tracev(g->utgs->buffer[i]);
	free_vplist(g->utgs);
	for(i=0;i<g->ctgs->size;i++) free_tracev(g->ctgs->buffer[i]);
	free_vplist(g->ctgs);
	free(g);
}

thread_beg_def(mdbg);
Graph *g;
reg_t reg;
wtovlv *hits;
int align_mode;
thread_end_def(mdbg);

thread_beg_func(mdbg);
Graph *g;
DMO *wt;
DMOAux *aux;
hzmhash *hash;
hzmrv *seeds;
volatile reg_t *reg;
wtovlv *hits;
wt_ovl_t *hit;
uint32_t i;
g = mdbg->g;
wt = g->dmo;
reg = (reg_t*)&mdbg->reg;
hits = mdbg->hits;
aux = init_dmoaux();
hash = init_hzmhash(1023);
seeds = init_hzmrv(32);

thread_beg_loop(mdbg);
if(reg->closed) continue;
clear_wtovlv(hits);
query_dmo(wt, 0xFFFFFFFFU, wt->rdseqs, wt->reads->buffer[reg->rid].rdoff + reg->beg, reg->end - reg->beg, mdbg->align_mode, hits, NULL, aux);
for(i=0;i<hits->size;i++){
	hit = ref_wtovlv(hits, i);
	hit->qb = num_max(hit->qb - hit->tb, 0);
	hit->qe = num_min(hit->qe + reg->end - reg->beg - hit->te, (int)wt->reads->buffer[hit->pb2].rdlen);
}
sort_array(hits->buffer, hits->size, wt_ovl_t, num_cmpgt(a.pb2, b.pb2));
for(i=1;i<hits->size;i++){
	if(hits->buffer[i-1].pb2 == hits->buffer[i].pb2){
		hits->buffer[i-1].score = WT_OVL_NULL_SCORE;
		hits->buffer[i  ].score = WT_OVL_NULL_SCORE;
	}
}
thread_end_loop(mdbg);

free_dmoaux(aux);
free_hzmhash(hash);
free_hzmrv(seeds);
thread_end_func(mdbg);

int update_regs_core_graph(Graph *g, rdregv *regs, rnk_ref_t *nd, u8v *ins){
	read_t *rd;
	reg_t *reg, *r;
	rd_reg_t *hit, *lst;
	node_t *n;
	uint64_t idx, nidx, pre;
	uint32_t i, x, y, conflit, fix, cnt;
	int ret;
	cnt = 0;
	lst = hit = NULL;
	fix = nd->fix;
	encap_regv(g->regs, nd->cnt);
	ret = 0;
	clear_u8v(ins);
	lst = NULL;
	for(i=0;i<nd->cnt;i++){
		hit = ref_rdregv(regs, i + nd->idx);
		if(hzm_debug){
			fprintf(hzm_debug_out, "update_regs:HIT\t%s\t%c\t%d\t%d\t%d\n", g->dmo->reads->buffer[hit->rid].tag, "+-"[hit->dir], g->dmo->reads->buffer[hit->rid].rdlen, hit->beg, hit->end);
		}
		conflit = 0;
		if(lst && hit->rid == lst->rid){
			lst->closed = hit->closed = 1;
			conflit = 1;
		}
		pre = 0;
		//check whether it has conflit with previous regs
		if(conflit == 0){
			rd = ref_readv(g->reads, hit->rid);
			idx = rd->regs.idx;
			while(idx){
				reg = ref_regv(g->regs, idx);
				if(hit->beg >= reg->beg) pre = idx;
				idx = reg->read_link;
				if(hit->end <= reg->beg) break;
				if(hit->beg >= reg->end) continue;
				if(hit->beg <= reg->beg){
					x = reg->beg;
					if(hit->end >= reg->end){
						y = reg->end;
					} else y = hit->end;
				} else {
					x = hit->beg;
					if(hit->end <= reg->end){
						y = hit->end;
					} else y = reg->end;
				}
				if(x < y){
					if(x + g->regovl < y) conflit = 1;
				}
				if(conflit){
					hit->closed = 1; // if conflit with non-major node, still discard this hit
					if(hzm_debug){
						fprintf(hzm_debug_out, "conflit_regs:%d\t%d\n", reg->beg, reg->end);
						node_t *n2;
						reg_t  *r2;
						uint32_t i2;
						n2 = ref_nodev(g->nodes, reg->node);
						for(i2=0;i2<n2->regs.cnt;i2++){
							r2 = ref_regv(g->regs, n2->regs.idx + i2);
							fprintf(hzm_debug_out, "\t%s\t%c\t%d\t%d\n", g->dmo->reads->buffer[r2->rid].tag, "+-"[r2->dir], r2->beg, r2->end);
						}
					}
				}
			}
		}
		push_u8v(ins, pre);
	}
	for(i=cnt=0;i<nd->cnt;i++) if(regs->buffer[i+nd->idx].closed == 0) cnt ++;
	if(hzm_debug){ fprintf(hzm_debug_out, "update_regs:%d\t%d\n", cnt, (int)nd->cnt); }
	if(!fix){
		if(cnt < g->min_node_cov) return 1;
		if((int)cnt < nd->cnt / 2) return 1;
	}
	nidx = g->nodes->size;
	n = NULL;
	for(i=0;i<nd->cnt;i++){
		hit = ref_rdregv(regs, i + nd->idx);
		if(hit->closed) continue;
		rd = ref_readv(g->reads, hit->rid);
		if(n == NULL){
			n = next_ref_nodev(g->nodes);
			n->unvisit = 0;
			n->closed = 0;
			n->single_in = 0;
			n->bt_visit = 0;
			n->bt_idx = 0;
			n->init_end = 0;
			n->regs.idx = g->regs->size;
			n->regs.cnt = 0;
			n->edges[0] = PTR_REF_NULL;
			n->edges[1] = PTR_REF_NULL;
		}
		n->regs.cnt ++;
		reg = next_ref_regv(g->regs);
		reg->node = nidx;
		reg->rid  = hit->rid;
		reg->dir  = hit->dir;
		reg->beg  = hit->beg;
		reg->end  = hit->end;
		reg->closed = 0;
		reg->read_link = 0;
		if(1){
			if(ins->buffer[i]){
				r = ref_regv(g->regs, ins->buffer[i]);
				reg->read_link = r->read_link;
				r->read_link = g->regs->size - 1;
			} else {
				reg->read_link = rd->regs.idx;
				rd->regs.idx = g->regs->size - 1;
			}
		} else {
			if(rd->regs.idx){
				r = ref_regv(g->regs, rd->regs.idx);
				if(r->beg > reg->beg){
					if(ins->buffer[i] != 0){
						fprintf(stderr, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						abort();
					}
					reg->read_link = rd->regs.idx;
					rd->regs.idx = g->regs->size - 1;
				} else {
					while(1){
						if(r->read_link == 0) break;
						if(g->regs->buffer[r->read_link].beg > reg->beg) break;
						r = ref_regv(g->regs, r->read_link);
					}
					if(ins->buffer[i] + g->regs->buffer != r){
						fprintf(stderr, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						abort();
					}
					reg->read_link = r->read_link;
					r->read_link = g->regs->size - 1;
				}
			} else {
				if(ins->buffer[i] != 0){
					fprintf(stderr, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
				reg->read_link = rd->regs.idx;
				rd->regs.idx = g->regs->size - 1;
			}
		}
		rd->regs.cnt ++;
	}
	return ret;
}

void _add_rdreps(ptr_ref_t *rds, rdrepv *reps, wtovlv *hits){
	ptr_ref_t *rd;
	rd_rep_t *rep, *r;
	wt_ovl_t *hit;
	uint32_t i;
	for(i=0;i<hits->size;i++){
		hit = ref_wtovlv(hits, i);
		if(hit->score == WT_OVL_NULL_SCORE) continue;
		rd = rds + hit->pb2;
		rep = next_ref_rdrepv(reps);
		rep->rid = hit->pb2;
		rep->dir = hit->dir2;
		rep->beg = hit->qb;
		rep->end = hit->qe;
		rep->closed = 0;
		rep->read_link = 0;
		if(rd->idx){
			r = ref_rdrepv(reps, rd->idx);
			if(r->beg > rep->beg){
				rep->read_link = rd->idx;
				rd->idx = reps->size - 1;
			} else {
				while(1){
					if(r->read_link == 0) break;
					if(reps->buffer[r->read_link].beg > rep->beg) break;
					r = ref_rdrepv(reps, r->read_link);
				}
				rep->read_link = r->read_link;
				r->read_link = reps->size - 1;
			}
		} else {
			rd->idx = reps->size - 1;
		}
		rd->cnt ++;
	}
}

int _conflict_rdreps(u8i *_idx, rdrepv *reps, uint32_t beg, uint32_t end, uint32_t repovl){
	uint64_t idx;
	uint32_t x, y;
	rd_rep_t *rep;
	idx = *_idx;
	while(idx){
		rep = ref_rdrepv(reps, idx);
		idx = rep->read_link;
		if(end <= rep->beg) return 0;
		if(beg >= rep->end){
			*_idx = idx;
			continue;
		}
		if(beg <= rep->beg){
			x = rep->beg;
			if(end >= rep->end) return 1;
			y = end;
		} else {
			x = beg;
			if(end <= rep->end) return 1;
			y = rep->end;
		}
		if(x + repovl < y) return 1;
	}
	return 0;
}

void build_nodes_graph(Graph *g, int align_mode, u8i maxbp, int ncpu, FILE *alno){
	wt_ovl_t *hit;
	u8v *ins;
	u8i idx, rep_idx, rank, kcnts[256], nskip, upds[3], nbp, fix_node;
	uint32_t nqry, rid, i, x, regoff, dep, beg, end, ib, ie, mi, max_node_cov, round;
	pbread_t *pb;
	read_t *rd;
	rd = NULL;
	read_t *wushigang = rd;
	read_t *tmp = wushigang;
	wushigang = tmp;
	ptr_ref_t *rds;
	rdrepv *reps;
	rdregv *regs;
	rnkrefv *nds;
	rnk_ref_t *nd;
	thread_preprocess(mdbg);
	if(hzm_debug) ncpu = 1;
	thread_beg_init(mdbg, ncpu);
	mdbg->g = g;
	memset((void*)&mdbg->reg, 0, sizeof(reg_t));
	mdbg->reg.closed = 1;
	mdbg->hits = init_wtovlv(32);
	mdbg->align_mode = align_mode;
	thread_end_init(mdbg);
	if(g->rep_filter) rds = calloc(g->reads->size, sizeof(ptr_ref_t));
	else rds = NULL;
	reps = init_rdrepv(32);
	regs = init_rdregv(1024);
	nds = init_rnkrefv(1024);
	clear_regv(g->regs);
	next_ref_regv(g->regs);
	memset(kcnts, 0, sizeof(uint64_t) * 256);
	regoff = g->regoff;
	if(regoff < 1) regoff = 1;
	nqry = g->only_fix? g->n_fix : g->reads->size;
	nskip = 0;
	ib = g->only_fix? g->n_fix : 0;
	mi = g->reads->size;
	ie = 0;
	if(maxbp < g->dmo->nbp){
		max_node_cov = 1.0 * g->max_node_cov * maxbp / g->dmo->nbp;
		if(max_node_cov < 10) max_node_cov = 10;
	} else max_node_cov = g->max_node_cov;
	if(g->rep_filter > 1) max_node_cov = 0; // will be much much faster
	rank = 0;
	fix_node = 0;
	round = 0;
	while(ib < mi){
		nbp = 0;
		ie = ib;
		while(ie < mi && nbp < maxbp) nbp += g->dmo->reads->buffer[ie++].rdlen;
		clear_index_dmo(g->dmo);
		fprintf(hzm_debug_out, "[%s] indexing reads[%u,%u] (%llu bp), %d threads\n", date(), ib, ie, nbp, ncpu); fflush(hzm_debug_out);
		index_dmo(g->dmo, ib, ie, ncpu);
		fprintf(hzm_debug_out, "[%s] Done\n", date()); fflush(hzm_debug_out);
		ib = ie;
		rank = 0;
		round ++;
		if(round == 1) fix_node = 0;
		for(rid=0;rid<nqry;rid++){
			if(!hzm_debug && (rid % 1000) == 0){ fprintf(hzm_debug_out, "\r%u|%llu", rid, (u8i)regs->size); fflush(hzm_debug_out); }
			pb = ref_pbreadv(g->dmo->reads, rid);
			rd = ref_readv(g->reads, rid);
			rep_idx = rds[rid].idx;// newly added rd_rep_t may be missed in confliction checking, but no worry for that checking is just a speedup strategy.
			for(x=0;x+g->reglen<=pb->rdlen;x+=regoff,rank++){
				if(g->rep_filter && _conflict_rdreps(&rep_idx, reps, x, x + g->reglen, g->regovl)){
					nskip ++;
					continue;
				}
				thread_wait_one(mdbg);
				if(mdbg->reg.closed == 0){
					//if(mdbg->hits->size > 1){
						for(i=0;i<mdbg->hits->size;i++){
							hit = ref_wtovlv(mdbg->hits, i);
							if(hit->score == WT_OVL_NULL_SCORE) continue;
							if(hit->dir2){
								beg = g->dmo->reads->buffer[hit->pb2].rdlen - hit->qe;
								end = g->dmo->reads->buffer[hit->pb2].rdlen - hit->qb;
								hit->qb = beg; hit->qe = end;
							}
							push_rdregv(regs, (rd_reg_t){mdbg->reg.node, hit->pb2, hit->dir2, hit->qb, hit->qe, 0});
						}
					//}
					if(g->rep_filter && mdbg->hits->size > max_node_cov) _add_rdreps(rds, reps, mdbg->hits);
					if(alno){
						for(i=0;i<mdbg->hits->size;i++){
							hit = ref_wtovlv(mdbg->hits, i);
							fprintf(alno, "%s\t+\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%d\n", g->dmo->reads->buffer[mdbg->reg.rid].tag, mdbg->reg.beg, mdbg->reg.end,
								g->dmo->reads->buffer[hit->pb2].tag, "+-"[hit->dir2], g->dmo->reads->buffer[hit->pb2].rdlen, hit->qb, hit->qe, hit->aln, hit->mat);
						}
					}
					if(hzm_debug){
						fprintf(hzm_debug_out, "QUERY: %s\t+\t%d\t%d\n", g->dmo->reads->buffer[mdbg->reg.rid].tag, mdbg->reg.beg, mdbg->reg.end);
						for(i=0;i<mdbg->hits->size;i++){
							hit = ref_wtovlv(mdbg->hits, i);
							fprintf(hzm_debug_out, "\t%s\t%c\t%d\t%d\t%d\t%d\t%d\n", g->dmo->reads->buffer[hit->pb2].tag, "+-"[hit->dir2], g->dmo->reads->buffer[hit->pb2].rdlen, hit->qb, hit->qe, hit->aln, hit->mat);
						}
					}
				}
				mdbg->reg = (reg_t){rank, rid, 0, x, x + g->reglen, 0, 0};
				clear_wtovlv(mdbg->hits);
				thread_wake(mdbg);
			}
			if(rid < g->n_fix && round == 1) fix_node = rank;
		}
		thread_wait_all(mdbg);
		thread_beg_iter(mdbg);
		if(mdbg->reg.closed == 0){
					//if(mdbg->hits->size > 1){
						for(i=0;i<mdbg->hits->size;i++){
							hit = ref_wtovlv(mdbg->hits, i);
							if(hit->score == WT_OVL_NULL_SCORE) continue;
							if(hit->dir2){
								beg = g->dmo->reads->buffer[hit->pb2].rdlen - hit->qe;
								end = g->dmo->reads->buffer[hit->pb2].rdlen - hit->qb;
								hit->qb = beg; hit->qe = end;
							}
							push_rdregv(regs, (rd_reg_t){mdbg->reg.node, hit->pb2, hit->dir2, hit->qb, hit->qe, 0});
						}
					//}
					if(g->rep_filter && mdbg->hits->size > max_node_cov) _add_rdreps(rds, reps, mdbg->hits);
					if(alno){
						for(i=0;i<mdbg->hits->size;i++){
							hit = ref_wtovlv(mdbg->hits, i);
							fprintf(alno, "%s\t+\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%d\n", g->dmo->reads->buffer[mdbg->reg.rid].tag, mdbg->reg.beg, mdbg->reg.end,
								g->dmo->reads->buffer[hit->pb2].tag, "+-"[hit->dir2], g->dmo->reads->buffer[hit->pb2].rdlen, hit->qb, hit->qe, hit->aln, hit->mat);
						}
					}
					if(hzm_debug){
						fprintf(hzm_debug_out, "QUERY: %s\t+\t%d\t%d\n", g->dmo->reads->buffer[mdbg->reg.rid].tag, mdbg->reg.beg, mdbg->reg.end);
						for(i=0;i<mdbg->hits->size;i++){
							hit = ref_wtovlv(mdbg->hits, i);
							fprintf(hzm_debug_out, "\t%s\t%c\t%d\t%d\t%d\t%d\t%d\n", g->dmo->reads->buffer[hit->pb2].tag, "+-"[hit->dir2], g->dmo->reads->buffer[hit->pb2].rdlen, hit->qb, hit->qe, hit->aln, hit->mat);
						}
					}
					mdbg->reg.closed = 1;
		}
		thread_end_iter(mdbg);
		if(ib < mi && !hzm_debug) fprintf(hzm_debug_out, "\r");
	}
	thread_beg_close(mdbg);
	free_wtovlv(mdbg->hits);
	thread_end_close(mdbg);
	if(!hzm_debug) fprintf(hzm_debug_out, "\r%u reads|%llu+%llu intervals|total hits %llu\n", nqry, (u8i)rank, nskip, (u8i)regs->size);
	if(rds) free(rds);
	if(reps) free_rdrepv(reps);
	reset_index_dmo(g->dmo);
	fprintf(hzm_debug_out, "[%s] sorting alignments\n", date());
	psort_array(regs->buffer, regs->size, rd_reg_t, ncpu, num_cmpgtxx(a.node, b.node, a.rid, b.rid, a.beg, b.beg));
	rank = 0xFFFFFFFFFFFFFFFFLLU;
	nd = NULL;
	for(idx=0;idx<regs->size;idx++){
		if(regs->buffer[idx].node != rank){
			nd = next_ref_rnkrefv(nds);
			nd->idx = idx;
			nd->rank = regs->buffer[idx].node;
			nd->fix = regs->buffer[idx].node < fix_node;
			nd->cnt = 1;
			rank = regs->buffer[idx].node;
		} else {
			nd->cnt ++;
		}
	}
	psort_array(nds->buffer, nds->size, rnk_ref_t, ncpu, num_cmpgtx(b.cnt, a.cnt, a.rank, b.rank));
	if(0){
		for(idx=0;idx<nds->size;idx++){
			dep = num_min(nds->buffer[idx].cnt, 255);
			kcnts[dep] ++;
		}
		for(i=1;i<51;i++){
			fprintf(hzm_debug_out, "%10llu", (long long unsigned int)kcnts[i]);
			if(((i - 1) % 10) == 9) fprintf(hzm_debug_out, "\n");
		}
		if(((i - 1) % 10) != 0) fprintf(hzm_debug_out, "\n");
	}
	fprintf(hzm_debug_out, "[%s] selecing important intervals\n", date());
	upds[0] = upds[1] = upds[2] = 0;
	ins = init_u8v(32);
	for(idx=0;idx<nds->size;idx++){
		nd = ref_rnkrefv(nds, idx);
		//if(nd->cnt > g->max_node_cov) continue;
		if(nd->cnt < 2) break;
		upds[update_regs_core_graph(g, regs, nd, ins)] ++;
	}
	free_u8v(ins);
	free_rdregv(regs);
	free_rnkrefv(nds);
	encap_regv(g->regs, 1);
	memset(g->regs->buffer + g->regs->size, 0xFF, sizeof(reg_t));
	if(!hzm_debug) fprintf(hzm_debug_out, "[%s] Intervals: kept %llu, discarded %llu\n", date(), upds[0], upds[1]);
	print_proc_stat_info(0);
}

int is_conflict_reg_graph(Graph *g, uint32_t rid,  uint32_t beg, uint32_t end, uint32_t repovl){
	uint64_t idx;
	uint32_t x, y;
	read_t *rd;
	reg_t *reg;
	rd = ref_readv(g->reads, rid);
	idx = rd->regs.idx;
	while(idx){
		reg = ref_regv(g->regs, idx);
		idx = reg->read_link;
		if(end <= reg->beg) return 0;
		if(beg >= reg->end) continue;
		if(beg <= reg->beg){
			x = reg->beg;
			if(end >= reg->end) return 1;
			y = end;
		} else {
			x = beg;
			if(end <= reg->end) return 1;
			y = reg->end;
		}
		if(x + repovl < y) return 1;
	}
	return 0;
}

void fast_build_nodes_graph(Graph *g, int align_mode, int ncpu, FILE *alno){
	wt_ovl_t *hit;
	u8v *ins;
	u8i rank, kcnts[256], nskip, maxbp, nbp, fix_node;
	fix_node = 0;
	u8i wushigang = fix_node;
	u8i tmp = wushigang;
	wushigang = tmp;
	uint32_t nqry, rid, i, x, regoff, dep, beg, end, ib, ie, mi, round, pblen, rankinc;
	rdregv *regs;
	rnk_ref_t *nd, ND;
	thread_preprocess(mdbg);
	if(hzm_debug) ncpu = 1;
	thread_beg_init(mdbg, ncpu);
	mdbg->g = g;
	memset((void*)&mdbg->reg, 0, sizeof(reg_t));
	mdbg->reg.closed = 1;
	mdbg->hits = init_wtovlv(32);
	mdbg->align_mode = align_mode;
	thread_end_init(mdbg);
	regs = init_rdregv(1024);
	ins = init_u8v(32);
	clear_regv(g->regs);
	next_ref_regv(g->regs);
	memset(kcnts, 0, sizeof(uint64_t) * 256);
	regoff = g->regoff;
	if(regoff < 1) regoff = 1;
	nqry = g->only_fix? g->n_fix : g->reads->size;
	nskip = 0;
	ib = g->only_fix? g->n_fix : 0;
	mi = g->reads->size;
	ie = 0;
	maxbp = g->dmo->nbp;
	rank = 0;
	fix_node = 0;
	round = 0;
	while(ib < mi){
		nbp = 0;
		ie = ib;
		while(ie < mi && nbp < maxbp) nbp += g->dmo->reads->buffer[ie++].rdlen;
		clear_index_dmo(g->dmo);
		fprintf(hzm_debug_out, "[%s] indexing reads[%u,%u] (%llu bp), %d threads\n", date(), ib, ie, nbp, ncpu); fflush(hzm_debug_out);
		index_dmo(g->dmo, ib, ie, ncpu);
		fprintf(hzm_debug_out, "[%s] Done\n", date()); fflush(hzm_debug_out);
		//print_proc_stat_info(0);
		ib = ie;
		rank = 0;
		round ++;
		if(round == 1) fix_node = 0;
		for(rid=0;rid<=nqry;rid++){
			if(rid < nqry){
				if(!hzm_debug && (rid % 1000) == 0){ fprintf(hzm_debug_out, "\r%u|%llu", rid, (u8i)g->regs->size); fflush(hzm_debug_out); }
				pblen = ref_pbreadv(g->dmo->reads, rid)->rdlen;
				rankinc = 1;
			} else {
				pblen = g->reglen + g->regoff * ncpu;
				rankinc = 0;
			}
			for(x=0;x+g->reglen<=pblen;x+=regoff,rank+=rankinc){
				if(rid < nqry && is_conflict_reg_graph(g, rid, x, x + g->reglen, g->regovl)){ nskip ++; continue; }
				if(rid < nqry){
					thread_wait_one(mdbg);
				} else {
					thread_wait_next(mdbg);
				}
				if(mdbg->reg.closed == 0){
					if(hzm_debug){
						fprintf(hzm_debug_out, "QUERY: %s\t+\t%d\t%d\n", g->dmo->reads->buffer[mdbg->reg.rid].tag, mdbg->reg.beg, mdbg->reg.end);
						for(i=0;i<mdbg->hits->size;i++){
							hit = ref_wtovlv(mdbg->hits, i);
							fprintf(hzm_debug_out, "\t%s\t%c\t%d\t%d\t%d\t%d\t%d\n", g->dmo->reads->buffer[hit->pb2].tag, "+-"[hit->dir2], g->dmo->reads->buffer[hit->pb2].rdlen, hit->qb, hit->qe, hit->aln, hit->mat);
						}
					}
					if(alno){
						for(i=0;i<mdbg->hits->size;i++){
							hit = ref_wtovlv(mdbg->hits, i);
							fprintf(alno, "%s\t+\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%d\n", g->dmo->reads->buffer[mdbg->reg.rid].tag, mdbg->reg.beg, mdbg->reg.end,
								g->dmo->reads->buffer[hit->pb2].tag, "+-"[hit->dir2], g->dmo->reads->buffer[hit->pb2].rdlen, hit->qb, hit->qe, hit->aln, hit->mat);
						}
					}
					if(mdbg->hits->size >= 2){
						clear_rdregv(regs);
						for(i=0;i<mdbg->hits->size;i++){
							hit = ref_wtovlv(mdbg->hits, i);
							if(hit->score == WT_OVL_NULL_SCORE) continue;
							if(hit->dir2){
								beg = g->dmo->reads->buffer[hit->pb2].rdlen - hit->qe;
								end = g->dmo->reads->buffer[hit->pb2].rdlen - hit->qb;
								hit->qb = beg; hit->qe = end;
							}
							push_rdregv(regs, (rd_reg_t){mdbg->reg.node, hit->pb2, hit->dir2, hit->qb, hit->qe, 0});
						}
						if(regs->size >= 2){
							nd = &ND;
							nd->idx = 0;
							nd->cnt = regs->size;
							nd->rank = mdbg->reg.node;
							nd->fix = mdbg->reg.rid < g->n_fix;
							update_regs_core_graph(g, regs, nd, ins);
							dep = num_min(255, regs->size);
							kcnts[dep] ++;
						}
					}
					mdbg->reg.closed = 1;
				}
				if(rid < nqry){
					mdbg->reg = (reg_t){rank, rid, 0, x, x + g->reglen, 0, 0};
					clear_wtovlv(mdbg->hits);
					thread_wake(mdbg);
				}
			}
			if(rid < g->n_fix && round == 1) fix_node = rank;
		}
		if(ib < mi && !hzm_debug) fprintf(hzm_debug_out, "\r");
	}
	thread_beg_close(mdbg);
	free_wtovlv(mdbg->hits);
	thread_end_close(mdbg);
	if(!hzm_debug) fprintf(hzm_debug_out, "\r%u reads | %llu intervals | skipped %llu | total hits %llu\n", nqry, (u8i)rank, nskip, (u8i)g->regs->size);
	reset_index_dmo(g->dmo);
	if(0){
		for(i=1;i<51;i++){
			fprintf(hzm_debug_out, "%10llu", (long long unsigned int)kcnts[i]);
			if(((i - 1) % 10) == 9) fprintf(hzm_debug_out, "\n");
		}
		if(((i - 1) % 10) != 0) fprintf(hzm_debug_out, "\n");
	}
	free_u8v(ins);
	free_rdregv(regs);
	encap_regv(g->regs, 1);
	memset(g->regs->buffer + g->regs->size, 0xFF, sizeof(reg_t));
	print_proc_stat_info(0);
}

//#define DEBUG_REFINE_NODE

void refine_node_core_graph(Graph *g, DMOPar *par, u4i node, u32hash *rdh, uihash *ndh, hzmhash *hash, hzmrv *seeds, rdregv *regs, wtovlv *hits, rdregv *rets, DMOAux *aux){
	node_t *nd, *n;
	read_t *rd;
	reg_t  *rg, *r, *reg;
	rd_reg_t *rr;
	uihash_t *u;
	wt_ovl_t *hit;
	u8i idx, pre;
	u4i i, j;
	int off, dir1, dir2, k, x, y, max, tmp, beg, end, no_refine;
	clear_u32hash(rdh);
	clear_uihash(ndh);
	clear_rdregv(rets);
	nd = ref_nodev(g->nodes, node);
	max = 0; reg = NULL;
#ifdef DEBUG_REFINE_NODE
	fprintf(hzm_debug_out, "%s: N%u cnt=%d\n", __FUNCTION__, node, nd->regs.cnt);
#endif
	no_refine = (nd->regs.cnt >= g->max_node_cov);
	for(i=0;i<nd->regs.cnt;i++){
		rg = ref_regv(g->regs, nd->regs.idx + i);
		if(rg->closed) continue;
#ifdef DEBUG_REFINE_NODE
		if(no_refine == 0) fprintf(hzm_debug_out, " REG R%u %s_%c_%d_%d\n", rg->rid, g->dmo->reads->buffer[rg->rid].tag,  "FR"[rg->dir], rg->beg, rg->end - rg->beg);
#endif
		push_rdregv(rets, (rd_reg_t){node, rg->rid, rg->dir, rg->beg, rg->end, 0});
		if(no_refine) continue;
		if(rg->end - rg->beg > max){ max = rg->end - rg->beg; reg = rg; }
		put_u32hash(rdh, rg->rid);
		rd = ref_readv(g->reads, rg->rid);
		// find the previous node on read
		idx = rd->regs.idx;
		pre = 0;
		while(idx){
			r = ref_regv(g->regs, idx);
			if(r->node == node) break;
			pre = idx;
			idx = r->read_link;
		}
		if(pre){
			r = ref_regv(g->regs, pre);
			kv_put_uihash(ndh, r->node, ((((int)(rg->beg + rg->end) / 2) - ((int)(r->beg + r->end) / 2)) << 2) | ((!rg->dir) << 1) | (r->dir ^ rg->dir));
		}
		if(rg->read_link){
			r = ref_regv(g->regs, rg->read_link);
			kv_put_uihash(ndh, r->node, ((((int)(r->beg + r->end) / 2) - ((int)(rg->beg + rg->end) / 2)) << 2) | ((rg->dir) << 1) | (r->dir ^ rg->dir));
		}
	}
	clear_rdregv(regs);
	reset_iter_uihash(ndh);
	while((u = ref_iter_uihash(ndh))){
		off  = u->val >> 2;
		dir1 = (u->val >> 1) & 0x01;
		dir2 = u->val & 0x01;
		n = ref_nodev(g->nodes, u->key);
		for(i=0;i<n->regs.cnt;i++){
			r = ref_regv(g->regs, n->regs.idx + i);
			k = dir2 ^ r->dir;
			if((r->dir ^ dir1 ^ dir2)){
				x = ((int)(r->beg + r->end)) / 2 + off;
			} else {
				x = ((int)(r->beg + r->end)) / 2 - off;
			}
			y = x + 1.5 * g->reglen;
			x = x - 1.5 * g->reglen;
			x = num_max(x, 0);
			y = num_min(y, (int)g->dmo->reads->buffer[r->rid].rdlen);
			if(x + (int)g->reglen > y) continue;
			if(exists_u32hash(rdh, r->rid)) continue;
			push_rdregv(regs, (rd_reg_t){node, r->rid, k, x, y, 0});
		}
	}
	if(regs->size == 0) return;
	index_reg_dmo(g->dmo->rdseqs, g->dmo->reads->buffer[reg->rid].rdoff + reg->beg, reg->end - reg->beg, par, hash, seeds);
	for(i=0;i<regs->size;i++){
		rr = ref_rdregv(regs, i);
#ifdef DEBUG_REFINE_NODE
		fprintf(hzm_debug_out, " TED R%u %s_%c_%d_%d\n", rr->rid, g->dmo->reads->buffer[rr->rid].tag,  "FR"[rr->dir], rr->beg, rr->end - rr->beg);
#endif
		clear_wtovlv(hits);
		query_reg_dmo(g->dmo->rdseqs, rr->rid, g->dmo->reads->buffer[rr->rid].rdoff + rr->beg, rr->end - rr->beg, reg->rid, reg->end - reg->beg, hash, seeds, hits, NULL, par, aux);
		for(j=0;j<hits->size;j++){
			hit = ref_wtovlv(hits, j);
			tmp = hit->dir2? g->dmo->reads->buffer[rr->rid].rdlen - rr->end : rr->beg;
			hit->qb += tmp; hit->qe += tmp;
			if(hit->dir2){
				beg = g->dmo->reads->buffer[hit->pb2].rdlen - hit->qe;
				end = g->dmo->reads->buffer[hit->pb2].rdlen - hit->qb;
				hit->qb = beg; hit->qe = end;
			}
			hit->dir2 = hit->dir2 ^ reg->dir;
			push_rdregv(rets, (rd_reg_t){node, hit->pb2, hit->dir2, hit->qb, hit->qe, 0});
#ifdef DEBUG_REFINE_NODE
			fprintf(hzm_debug_out, " FND R%u %s_%c_%d_%d %d %d %d %d\n", hit->pb2, g->dmo->reads->buffer[hit->pb2].tag,  "FR"[hit->dir2], hit->qb, hit->qe - hit->qb, hit->tb, hit->te, hit->mat, hit->aln);
#endif
			if(rets->size >= g->max_node_cov) break;
		}
		if(rets->size >= g->max_node_cov) break;
	}
	sort_array(rets->buffer, rets->size, rd_reg_t, num_cmpgt(a.rid, b.rid));
	for(i=1;i<rets->size;i++){
		if(rets->buffer[i-1].rid == rets->buffer[i].rid){
			rets->buffer[i-1].closed = 1;
			rets->buffer[i].closed = 1;
		}
	}
}

thread_beg_def(mreg);
Graph *g;
DMOPar *par;
rdregv *rets;
u4i node;
thread_end_def(mreg);

thread_beg_func(mreg);
u32hash *rdh;
uihash *ndh;
DMOAux *aux;
hzmhash *hash;
hzmrv *seeds;
rdregv *regs;
wtovlv *hits;
rdh = init_u32hash(1023);
ndh = init_uihash(1023);
aux = init_dmoaux();
hash = init_hzmhash(1023);
seeds = init_hzmrv(32);
regs = init_rdregv(1024);
hits = init_wtovlv(32);
thread_beg_loop(mreg);
if(mreg->node == 0xFFFFFFFFU) continue;
refine_node_core_graph(mreg->g, mreg->par, mreg->node, rdh, ndh, hash, seeds, regs, hits, mreg->rets, aux);
thread_end_loop(mreg);
free_u32hash(rdh);
free_uihash(ndh);
free_dmoaux(aux);
free_hzmhash(hash);
free_hzmrv(seeds);
free_rdregv(regs);
free_wtovlv(hits);
thread_end_func(mreg);

void refine_nodes_graph(Graph *g, DMOPar *par, int ncpu, FILE *log){
	rdregv *regs;
	rnkrefv *nds;
	rnk_ref_t *nd;
	u8v *ins;
	u8i idx, rank, kcnts[256], upds[3];
	u4i node, i, ret, dep;
#ifdef DEBUG_REFINE_NODE
	ncpu = 1;
#endif
	thread_preprocess(mreg);
	regs = init_rdregv(1024);
	nds = init_rnkrefv(1024);
	memset(kcnts, 0, sizeof(uint64_t) * 256);
	thread_beg_init(mreg, ncpu);
	mreg->g = g;
	mreg->par = par;
	mreg->rets = init_rdregv(1024);
	mreg->node = 0xFFFFFFFFU;
	thread_end_init(mreg);
	fprintf(hzm_debug_out, "[%s] refining node alignments (%u,%llu)\n", date(), (u4i)g->nodes->size, (u8i)g->regs->size);
	for(node=0;node<g->nodes->size+(u4i)ncpu;node++){
		if(node < g->nodes->size){
			if(!hzm_debug && (node % 1000) == 0){ fprintf(hzm_debug_out, "\r%u|%llu", node, (u8i)regs->size); fflush(hzm_debug_out); }
			thread_wait_one(mreg);
		} else {
			thread_wait_next(mreg);
		}
		if(mreg->node != 0xFFFFFFFFU){
			for(i=ret=0;i<mreg->rets->size;i++){
				if(mreg->rets->buffer[i].closed) continue;
				ret ++;
				push_rdregv(regs, mreg->rets->buffer[i]);
			}
			if(log){
				fprintf(log, "update n%u: %u -> %u\n", node, g->nodes->buffer[mreg->node].regs.cnt, ret);
			}
		}
		if(node < g->nodes->size){
			mreg->node = node;
			thread_wake(mreg);
		} else {
			mreg->node = 0xFFFFFFFFU;
		}
	}
	thread_beg_close(mreg);
	free_rdregv(mreg->rets);
	thread_end_close(mreg);
	if(!hzm_debug) fprintf(hzm_debug_out, "\r%u nodes | hits %llu -> %llu\n", (u4i)g->nodes->size, (u8i)g->regs->size, (u8i)regs->size);
	clear_regv(g->regs);
	next_ref_regv(g->regs);
	clear_nodev(g->nodes);
	for(i=0;i<g->reads->size;i++){
		g->reads->buffer[i].regs.idx = 0;
		g->reads->buffer[i].regs.cnt = 0;
	}
	fprintf(hzm_debug_out, "[%s] sorting alignments\n", date());
	psort_array(regs->buffer, regs->size, rd_reg_t, ncpu, num_cmpgtx(a.node, b.node, a.rid, b.rid));
	rank = 0xFFFFFFFFFFFFFFFFLLU;
	nd = NULL;
	for(idx=0;idx<regs->size;idx++){
		if(regs->buffer[idx].node != rank){
			nd = next_ref_rnkrefv(nds);
			nd->idx = idx;
			nd->rank = regs->buffer[idx].node;
			nd->fix = 0;
			nd->cnt = 1;
			rank = regs->buffer[idx].node;
		} else {
			nd->cnt ++;
		}
	}
	psort_array(nds->buffer, nds->size, rnk_ref_t, ncpu, num_cmpgtx(b.cnt, a.cnt, a.rank, b.rank));
	if(0){
		for(idx=0;idx<nds->size;idx++){
			dep = num_min(nds->buffer[idx].cnt, 255);
			kcnts[dep] ++;
		}
		for(i=1;i<51;i++){
			fprintf(hzm_debug_out, "%10llu", (long long unsigned int)kcnts[i]);
			if(((i - 1) % 10) == 9) fprintf(hzm_debug_out, "\n");
		}
		if(((i - 1) % 10) != 0) fprintf(hzm_debug_out, "\n");
	}
	fprintf(hzm_debug_out, "[%s] selecing important intervals\n", date());
	upds[0] = upds[1] = upds[2] = 0;
	ins = init_u8v(32);
	for(idx=0;idx<nds->size;idx++){
		nd = ref_rnkrefv(nds, idx);
		upds[update_regs_core_graph(g, regs, nd, ins)] ++;
	}
	free_u8v(ins);
	free_rdregv(regs);
	free_rnkrefv(nds);
	encap_regv(g->regs, 1);
	memset(g->regs->buffer + g->regs->size, 0xFF, sizeof(reg_t));
}

void remove_all_edges_graph(Graph *g){
	node_t *n;
	uint64_t nid;
	free_edgev(g->edges);
	g->edges = init_edgev(32);
	free_edgehash(g->ehash);
	g->ehash = init_edgehash(1023);
	set_userdata_edgehash(g->ehash, g->edges);
	free_edgerefv(g->erefs);
	g->erefs = init_edgerefv(32);
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		n->edges[0] = n->edges[1] = PTR_REF_NULL;
	}
}

typedef struct {
	uint64_t idx:46; int off:18;
} edge_off_t;
define_list(edgeoffv, edge_off_t);

int estimate_edge_length(edge_off_t *ps, uint32_t size, uint32_t idxs[2]){
	int64_t tot, var;
	uint32_t i, b, e, mi;
	int v, max, len, avg, std;
	 idxs[0] = 0; idxs[1] = size;
	if(size == 0){ return 0; }
	if(size <= 2){ return ps[size/2].off; }
	b = 0; e = size;
	for(i=b,tot=0;i<e;i++) {tot += ps[i].off;} len = tot / (e - b);
	//fprintf(stdout, "b=%d\te=%d\tlen=%d\n", b, e, len);
	while(b + 2 < e){
		max = 0; mi = 0;
		for(i=b+1;i<e;i++){
			if(ps[i].off - ps[i-1].off > max){
				max = ps[i].off - ps[i-1].off;
				mi = i;
			}
		}
		if(max < len * 0.5) break;
		else if(max < 100) break;
		if(mi - b > e - mi) e = mi;
		else b = mi;
		for(i=b,tot=0;i<e;i++) {tot += ps[i].off;} avg = tot / (e - b);
		if(num_diff(avg, len) < num_max(avg * 0.2, 50)) break;
		len = avg;
	}
	//fprintf(stdout, "b=%d\te=%d\tlen=%d\n", b, e, len);
	if(0){
		if(b + 1 < e){
			for(i=b,var=0;i<e;i++){
				v = ((int)ps[i].off) - ((int)len);
				var += v * v;
			}
			std = sqrt(var / (e - b));
			//fprintf(stdout, "std=%d\n", std);
			b = 0; e = size;
			while(b < e && num_diff(ps[b].off, len) > 3 * std) b ++;
			while(b < e && num_diff(ps[e - 1].off, len) > 3 * std) e --;
		}
		idxs[0] = b;
		idxs[1] = e;
	}
	return len;
}

int calculate_edge_cov_off_graph(Graph *g, edge_t *e, edgeoffv *offs){
	node_t *n1, *n2;
	reg_t *r1, *r2;
	uint32_t dir1, dir2, cov, idxs[2];
	int off;
	n1 = ref_nodev(g->nodes, e->node1);
	n2 = ref_nodev(g->nodes, e->node2);
	r1 = ref_regv(g->regs, n1->regs.idx);
	r2 = ref_regv(g->regs, n2->regs.idx);
	cov = e->cov;
	clear_edgeoffv(offs);
	while(1){
		if(r1->rid > r2->rid){
			r2 ++;
			if(r2->node != e->node2) break;
		} else if(r1->rid < r2->rid){
			r1 ++;
			if(r1->node != e->node1) break;
		} else {
			if(r1->beg < r2->beg){
				if(r1->node < r2->node){
					dir1 = r1->dir;
					dir2 = r2->dir;
				} else {
					dir1 = !r2->dir;
					dir2 = !r1->dir;
				}
				off = r2->beg - r1->end;
			} else {
				if(r2->node < r1->node){
					dir1 = r2->dir;
					dir2 = r1->dir;
				} else {
					dir1 = !r1->dir;
					dir2 = !r2->dir;
				}
				off = ((int)r1->beg) - r2->end;
			}
			if(dir1 == e->dir1 && dir2 == e->dir2){
				push_edgeoffv(offs, (edge_off_t){e - g->edges->buffer, off});
			}
			r1 ++;
			if(r1->node != e->node1) break;
			r2 ++;
			if(r2->node != e->node2) break;
		}
	}
	e->off = estimate_edge_length(offs->buffer, offs->size, idxs);
	e->cov = idxs[1] - idxs[0];
	if(cov != e->cov) return 1;
	else return 0;
}

uint32_t estimate_genome_size(Graph *g, unsigned long long tot_bp, FILE *out){
	uint64_t kcnts[256];
	node_t *n;
	uint64_t nid, sum, ncnt, pmax;
	uint32_t i, dep, peak, mid;
	float avg;
	ncnt = g->nodes->size;
	memset(kcnts, 0, sizeof(uint64_t) * 256);
	sum = 0;
	for(nid=0;nid<ncnt;nid++){
		n = ref_nodev(g->nodes, nid);
		dep = num_min(n->regs.cnt, 255);
		sum += n->regs.cnt;
		kcnts[dep] ++;
	}
	mid = pmax = 0;
	while(mid < 255){
		pmax += kcnts[mid];
		if(pmax >= ncnt / 2) break;
		mid ++;
	}
	avg = 1.0 * sum / (ncnt + 1);
	fprintf(out, "[%s] median node cov = %d\n", date(), mid);
	return mid;
	//TODO: calculate the genome coverage
	for(i=1;i<51;i++){
		fprintf(out, "%10llu", (long long unsigned int)kcnts[i]);
		if(((i - 1) % 10) == 9) fprintf(out, "\n");
	}
	if(((i - 1) % 10) != 0) fprintf(out, "\n");
	return avg;
	pmax = 0; peak = avg;
	for(i=g->min_node_cov+1;i<256;i++){
		if(kcnts[i] > pmax){ peak = i; pmax = kcnts[i]; }
		else if(i > avg && kcnts[i] < 0.8 * pmax) break;
	}
	fprintf(out, "[%s] cov peak = %d\n", date(), peak);
	fprintf(out, "[%s] genome size = %llu\n", date(), tot_bp / peak);
	return peak;
}

void build_edges_graph(Graph *g, int ncpu){
	read_t *rd;
	node_t *n;
	edge_t *E;
	vplist *regs;
	reg_t *r1, *r2;
	edge_ref_t f1, f2;
	edgeoffv *offs;
	uint64_t idx, lst, cnt, x, *u;
	cnt = 0;
	uint64_t wushigang = cnt;
	uint64_t tmp = wushigang;
	wushigang = tmp;
	uint32_t rid, i, idxs[2];
	int exists;
	clear_edgev(g->edges);
	E = next_ref_edgev(g->edges);
	memset(E, 0, sizeof(edge_t));
	E->cov = 1;
	E->status = 0;
	clear_edgehash(g->ehash);
	offs = init_edgeoffv(32);
	regs = init_vplist(32);
	for(rid=0;rid<g->dmo->n_rd;rid++){
		rd = ref_readv(g->reads, rid);
		if(rd->regs.cnt < 2) continue;
		clear_vplist(regs);
		idx = rd->regs.idx;
		while(idx){
			r2 = ref_regv(g->regs, idx);
			idx = r2->read_link;
			//if(g->nodes->buffer[r2->node].closed) continue;
			for(i=0;i<regs->size;i++){
				r1 = (reg_t*)get_vplist(regs, i);
				E = ref_edgev(g->edges, 0);
				if(r1->node < r2->node){
					E->node1 = r1->node;
					E->node2 = r2->node;
					E->dir1  = r1->dir;
					E->dir2  = r2->dir;
				} else {
					E->node1 = r2->node;
					E->node2 = r1->node;
					E->dir1  = !r2->dir;
					E->dir2  = !r1->dir;
				}
				E->off = ((int)r2->beg) - r1->end;
				u = prepare_edgehash(g->ehash, 0, &exists);
				if(exists){
					if(g->edges->buffer[*u].cov < WT_MAX_EDGE_COV) g->edges->buffer[*u].cov ++;
				} else {
					*u = g->edges->size;
					push_edgev(g->edges, *E);
				}
				push_edgeoffv(offs, (edge_off_t){*u, ((int)r2->beg) - r1->end});
			}
			push_vplist(regs, r2);
		}
	}
	free_vplist(regs);
	if(g->edges->size == 0){
		free_edgeoffv(offs);
		return;
	}
	psort_array(offs->buffer, offs->size, edge_off_t, ncpu, num_cmpgtx(a.idx, b.idx, a.off, b.off));
	lst = 0;
	for(idx=1;idx<=offs->size;idx++){
		if(idx < offs->size && offs->buffer[idx].idx == offs->buffer[lst].idx) continue;
		//g->edges->buffer[offs->buffer[lst].idx].off = offs->buffer[(lst+idx)/2].cnt;
		g->edges->buffer[offs->buffer[lst].idx].off = estimate_edge_length(offs->buffer + lst, idx - lst, idxs);
		g->edges->buffer[offs->buffer[lst].idx].cov = idxs[1] - idxs[0];
		if(0){
			uint64_t m;
			for(m=lst;m<idx;m++) fprintf(stdout, "%u\t", offs->buffer[m].off);
			fprintf(stdout, "\n");
			for(m=lst+idxs[0];m<lst+idxs[1];m++) fprintf(stdout, "%u\t", offs->buffer[m].off);
			fprintf(stdout, "\n");
		}
		lst = idx;
	}
	free_edgeoffv(offs);
	clear_edgerefv(g->erefs);
	push_edgerefv(g->erefs, EDGE_REF_NULL);
	for(idx=1;idx<g->edges->size;idx++){
		if(g->edges->buffer[idx].cov < g->min_edge_cov){
			if(g->store_low_cov_edge) g->edges->buffer[idx].closed = WT_EDGE_CLOSED_LESS;
			else continue;
		}
		if(g->nodes->buffer[g->edges->buffer[idx].node1].closed || g->nodes->buffer[g->edges->buffer[idx].node2].closed){
			g->edges->buffer[idx].closed = WT_EDGE_CLOSED_HARD;
		}
		push_edgerefv(g->erefs, (edge_ref_t){idx, 0, 0});
		push_edgerefv(g->erefs, (edge_ref_t){idx, 1, 0});
	}
	psort_array(g->erefs->buffer + 1, g->erefs->size - 1, edge_ref_t, ncpu, num_cmpgtx(
		(a.flg? ((g->edges->buffer[a.idx].node2 << 1) | !g->edges->buffer[a.idx].dir2) : ((g->edges->buffer[a.idx].node1 << 1) | g->edges->buffer[a.idx].dir1)),
		(b.flg? ((g->edges->buffer[b.idx].node2 << 1) | !g->edges->buffer[b.idx].dir2) : ((g->edges->buffer[b.idx].node1 << 1) | g->edges->buffer[b.idx].dir1)),
		g->edges->buffer[a.idx].off, g->edges->buffer[b.idx].off));
	f1.idx = g->nodes->size; f1.flg = 0;
	cnt = 0;
	for(lst=idx=1;idx<g->erefs->size;idx++){
		if(g->erefs->buffer[idx].flg){
			f2.idx =  g->edges->buffer[g->erefs->buffer[idx].idx].node2;
			f2.flg = !g->edges->buffer[g->erefs->buffer[idx].idx].dir2;
		} else {
			f2.idx =  g->edges->buffer[g->erefs->buffer[idx].idx].node1;
			f2.flg =  g->edges->buffer[g->erefs->buffer[idx].idx].dir1;
		}
		if(f1.idx == f2.idx && f1.flg == f2.flg) continue;
		if(lst < idx){
			n = ref_nodev(g->nodes, f1.idx);
			n->edges[f1.flg].idx = lst;
			n->edges[f1.flg].cnt = 0;
			for(x=lst;x+1<idx;x++){
				g->erefs->buffer[x].next = x + 1;
				if(g->edges->buffer[g->erefs->buffer[x].idx].closed == WT_EDGE_CLOSED_NULL) n->edges[f1.flg].cnt ++;
			}
			if(g->edges->buffer[g->erefs->buffer[x].idx].closed == WT_EDGE_CLOSED_NULL) n->edges[f1.flg].cnt ++;
		}
		lst = idx;
		f1 = f2;
	}
	if(lst < idx){
		n = ref_nodev(g->nodes, f1.idx);
		n->edges[f1.flg].idx = lst;
		n->edges[f1.flg].cnt = 0;
		for(x=lst;x+1<idx;x++){
			g->erefs->buffer[x].next = x + 1;
			if(g->edges->buffer[g->erefs->buffer[x].idx].closed == WT_EDGE_CLOSED_NULL) n->edges[f1.flg].cnt ++;
		}
		if(g->edges->buffer[g->erefs->buffer[x].idx].closed == WT_EDGE_CLOSED_NULL) n->edges[f1.flg].cnt ++;
	}
}

void load_nodes_graph(Graph *g, FileReader *fr){
	read_t *rd;
	node_t *n;
	reg_t *reg, *r;
	uint64_t nid;
	uint32_t i, nreg;
	char *str, *tok;
	int ncol;
	clear_nodev(g->nodes);
	clear_regv(g->regs);
	next_ref_regv(g->regs);
	while((ncol = fread_table(fr)) != -1){
		if(fr->line->string[0] == '#') continue;
		if(ncol < 2) continue;
		nreg = atoi(get_col_str(fr, 1));
		//if(nreg == 0) continue;
		if(nreg < g->min_node_cov) continue;
		nid = atol(get_col_str(fr, 0) + 1);
		while(g->nodes->size < nid){
			n = next_ref_nodev(g->nodes);
			n->closed = 1;
		}
		n = next_ref_nodev(g->nodes);
		n->unvisit = 0;
		n->closed = 0;
		n->single_in = 0;
		n->bt_visit = 0;
		n->bt_idx = 0;
		n->regs.idx = g->regs->size;
		n->regs.cnt = nreg;
		n->edges[0] = PTR_REF_NULL;
		n->edges[1] = PTR_REF_NULL;
		for(i=0;i<nreg;i++){
			reg = next_ref_regv(g->regs);
			reg->node = nid;
			reg->closed = 0;
			reg->read_link = 0;
			str = get_col_str(fr, 2 + i);
			if(0){
				tok = index(str, '_'); *tok = '\0';
				reg->rid = kv_get_cuhash(g->dmo->tag2idx, str);
				str = tok + 1; *tok = '_'; tok = index(str, '_'); *tok = '\0';
				reg->dir = (str[0] == 'R');
				str = tok + 1; *tok = '_'; tok = index(str, '_'); *tok = '\0';
				reg->beg = atoi(str);
				str = tok + 1;
				reg->end = atoi(str);
				reg->end += reg->beg;
			} else {
				tok = rindex(str, '_'); *tok = '\0'; tok ++;
				reg->end = atoi(tok);
				tok = rindex(str, '_'); *tok = '\0'; tok ++;
				reg->beg = atoi(tok);
				reg->end += reg->beg;
				tok = rindex(str, '_'); *tok = '\0'; tok ++;
				reg->dir = (tok[0] == 'R');
				reg->rid = kv_get_cuhash(g->dmo->tag2idx, str);
			}
			rd = ref_readv(g->reads, reg->rid);
			if(rd->regs.idx){
				r = ref_regv(g->regs, rd->regs.idx);
				if(r->beg > reg->beg){
					reg->read_link = rd->regs.idx;
					rd->regs.idx = g->regs->size - 1;
				} else {
					while(1){
						if(r->read_link == 0) break;
						if(g->regs->buffer[r->read_link].beg > reg->beg) break;
						r = ref_regv(g->regs, r->read_link);
					}
					reg->read_link = r->read_link;
					r->read_link = g->regs->size - 1;
				}
			} else {
				rd->regs.idx = g->regs->size - 1;
			}
			rd->regs.cnt ++;
		}
	}
	encap_regv(g->regs, 1);
	g->regs->buffer[g->regs->size].node = WT_MAX_NODE;
}

// MUST be called before build_edges
uint64_t mask_nodes_by_cov_graph(Graph *g, FILE *out){
	node_t *n;
	uint64_t ret, i;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->regs.cnt > g->max_node_cov || n->regs.cnt < g->min_node_cov){
		//if(n->regs.cnt > g->max_node_cov){
			n->closed = 1;
			ret ++;
			if(out) fprintf(out, "MASK_COV\tN%llu\t%u\n", (u8i)i, (u4i)n->regs.cnt);
		}
	}
	return ret;
}

#define MAX_BT_NIDX	0xFFFFFFU
typedef struct {
	u8i node:46, visit:1, closed:1, cov:11, fix:1, sub_dir:1;
	u8i flag;
	u8i bt_dir:1, bt_open:10, bt_nidx:24, bt_score:20, bt_step:8, bt_hit:1;
	ptr_ref_t edges[2];
} subnode_t;
#define subnode_hashcode(E) u64hashcode((E).node)
#define subnode_hashequals(E1, E2) (E1).node == (E2).node
define_hashset(subnodehash, subnode_t, subnode_hashcode, subnode_hashequals);

typedef struct {
	subnode_t *node;
	uint32_t cov:28, visit:1, fwd:1, dir:1, closed:1;
	uint32_t next;
} subedge_t;
define_list(subedgev, subedge_t);

int evaluate_node_connectivity_graph(Graph *g, uint64_t nid, u4v *rds, subnodehash *nodes, subedgev *edges, ptrrefv *stack){
	node_t *nd;
	read_t *rd;
	reg_t  *rg;
	subnode_t N, *n, *n1, *n2;
	subedge_t *e;
	ptr_ref_t *p;
	u8i idx, edx, aim;
	uint32_t i, k, k1, k2, cnt;
	int exists;
	// collect reads containing nid
	clear_u4v(rds);
	nd = ref_nodev(g->nodes, nid);
	for(i=0;i<nd->regs.cnt;i++){
		rg = ref_regv(g->regs, nd->regs.idx + i);
		push_u4v(rds, rg->rid);
	}
	// prepare nodes in subgraph
	clear_subnodehash(nodes);
	clear_subedgev(edges);
	next_ref_subedgev(edges);
	memset(&N, 0, sizeof(subnode_t));
	N.cov = 1;
	for(i=0;i<rds->size;i++){
		rd = ref_readv(g->reads, rds->buffer[i]);
		rg = NULL;
		idx = rd->regs.idx;
		while(idx){
			rg = ref_regv(g->regs, idx);
			idx = rg->read_link;
			N.node = rg->node;
			n = prepare_subnodehash(nodes, N, &exists);
			if(exists){
				n->cov ++;
			} else {
				*n = N;
			}
		}
	}
	// mask low cov nodes
	reset_iter_subnodehash(nodes);
	while((n = ref_iter_subnodehash(nodes))){
		if(n->cov < g->min_node_cov) n->closed = 1;
	}
	// build edges
	for(i=0;i<rds->size;i++){
		rd = ref_readv(g->reads, rds->buffer[i]);
		n1 = NULL;
		k1 = 0;
		idx = rd->regs.idx;
		while(idx){
			rg = ref_regv(g->regs, idx);
			idx = rg->read_link;
			N.node = rg->node;
			n2 = get_subnodehash(nodes, N);
			k2 = rg->dir;
			if(n2->closed) continue;
			if(n1){
				// link n1 to n2
				edx = n1->edges[k1].idx;
				while(edx){
					e = ref_subedgev(edges, edx);
					if(e->node == n2 && e->dir == k2){
						e->cov ++;
						break;
					}
					edx = e->next;
				}
				if(edx == 0){
					edx = edges->size;
					e = next_ref_subedgev(edges);
					e->node = n2;
					e->dir = k2;
					e->cov = 1;
					e->next = n1->edges[k1].idx;
					n1->edges[k1].idx = edx;
					n1->edges[k1].cnt ++;
				}
				// link rev n2 to rev n1
				edx = n2->edges[!k2].idx;
				while(edx){
					e = ref_subedgev(edges, edx);
					if(e->node == n1 && e->dir == !k1){
						e->cov ++;
						break;
					}
					edx = e->next;
				}
				if(edx == 0){
					edx = edges->size;
					e = next_ref_subedgev(edges);
					e->node = n1;
					e->dir = !k1;
					e->cov = 1;
					e->next = n2->edges[!k2].idx;
					n2->edges[!k2].idx = edx;
					n2->edges[!k2].cnt ++;
				}
			}
			n1 = n2;
			k1 = k2;
		}
	}
	// find the nid node
	N.node = nid;
	n = get_subnodehash(nodes, N);
	n->visit = 1;
	// checking whether its out-edges collapse into one node
	for(k=0;k<2;k++){
		if(n->edges[k].cnt > 64) return 0;
		if(n->edges[k].cnt < 2) continue;
		idx = n->edges[k].idx;
		cnt = 0;
		while(idx){
			e = ref_subedgev(edges, idx);
			idx = e->next;
			if(e->cov == 1) continue; // don't track low cov out-edges
			cnt ++;
		}
		aim = 0xFFFFFFFFFFFFFFFFLLU >> (64 - cnt);
		cnt = 0;
		exists = 0;
		if(k){
			reset_iter_subnodehash(nodes);
			while((n1 = ref_iter_subnodehash(nodes))){ n1->flag = 0; }
		}
		idx = n->edges[k].idx;
		while(idx){
			e = ref_subedgev(edges, idx);
			idx = e->next;
			if(e->cov == 1) continue; // don't track low cov out-edges
			e->node->flag |= 1LLU << cnt;
			cnt ++;
			reset_iter_subnodehash(nodes);
			while((n1 = ref_iter_subnodehash(nodes))){ n1->visit = 0; }
			n->visit = 1;
			clear_ptrrefv(stack);
			push_ptrrefv(stack, (ptr_ref_t){offset_subnodehash(nodes, e->node), e->dir});
			while(stack->size){
				p = peer_ptrrefv(stack);
				n1 = nodes->array + p->idx;
				k1 = p->cnt;
				stack->size --;
				if(n1->flag == aim){ exists = 1; break; }
				if(n1->visit) continue;
				n1->visit = 1;
				edx = n1->edges[k1].idx;
				while(edx){
					e = ref_subedgev(edges, edx);
					edx = e->next;
					if(e->node->visit) continue;
					e->node->flag |= n1->flag;
					push_ptrrefv(stack, (ptr_ref_t){offset_subnodehash(nodes, e->node), e->dir});
				}
				if(exists) break;
			}
			if(exists) break;
		}
		if(exists == 0) return 0;
	}
	return 1;
}

void print_subgraph_dot(Graph *g, u8i id, subnodehash *nodes, subedgev *edges, FILE *out){
	subnode_t *n;
	subedge_t *e;
	u8i idx;
	int k;
	fprintf(out, "digraph N%llu {\n", id);
	fprintf(out, " N%llu [style=filled fillcolor=yellow]\n", id);
	reset_iter_subnodehash(nodes);
	while((n = ref_iter_subnodehash(nodes))){
		if(n->closed) continue;
		fprintf(out, "N%llu [label=\"N%llu(%llu)\"]\n", (u8i)n->node, (u8i)n->node, (u8i)g->nodes->buffer[n->node].rep_idx);
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			while(idx){
				e = ref_subedgev(edges, idx);
				idx = e->next;
				fprintf(out, " N%llu -> N%llu [label=\"%c%c:%d\"]\n", (u8i)n->node, (u8i)e->node->node, "+-"[k], "+-"[e->dir], e->cov);
			}
		}
	}
	fprintf(out, "}\n");
}

thread_beg_def(mrep);
Graph *g;
u8i ret;
u8v *reps;
thread_end_def(mrep);

thread_beg_func(mrep);
subnodehash *nodes;
subedgev *edges;
u4v *rds;
ptrrefv *stack;
u8i nid, tidx, ncpu;
nodes = init_subnodehash(1023);
edges = init_subedgev(32);
rds = init_u4v(32);
stack = init_ptrrefv(32);
tidx = mrep->t_idx;
ncpu = mrep->n_cpu;
thread_beg_loop(mrep);
for(nid=tidx;nid<mrep->g->nodes->size;nid+=ncpu){
	if(mrep->g->nodes->buffer[nid].closed) continue;
	if(evaluate_node_connectivity_graph(mrep->g, nid, rds, nodes, edges, stack) == 0){
		if(0){
			print_subgraph_dot(mrep->g, nid, nodes, edges, stdout);
		}
		mrep->g->nodes->buffer[nid].closed = 1;
		push_u8v(mrep->reps, nid);
		mrep->ret ++;
	}
}
thread_end_loop(mrep);
free_subnodehash(nodes);
free_subedgev(edges);
free_u4v(rds);
free_ptrrefv(stack);
thread_end_func(mrep);

u8i mask_nodes_by_connectivity_graph(Graph *g, int ncpu, FILE *out){
	node_t *n;
	u8i ret, i;
	thread_preprocess(mrep);
	ret = 0;
	if(0) ncpu = 1;
	thread_beg_init(mrep, ncpu);
	mrep->g   = g;
	mrep->ret = 0;
	mrep->reps = init_u8v(32);
	thread_end_init(mrep);
	thread_wake_all(mrep);
	thread_wait_all(mrep);
	thread_beg_close(mrep);
	if(out){
		for(i=0;i<mrep->reps->size;i++){
			n = ref_nodev(g->nodes, mrep->reps->buffer[i]);
			fprintf(out, "N%llu\t%u\tconn\n", (u8i)mrep->reps->buffer[i], (u4i)n->regs.cnt);
		}
	}
	ret += mrep->ret;
	free_u8v(mrep->reps);
	thread_end_close(mrep);
	return ret;
}

uint64_t mask_possible_tip_nodes_graph(Graph *g){
	node_t *n;
	reg_t *r;
	uint64_t ret, i;
	uint32_t j, cnt;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		cnt = 0;
		for(j=0;j<n->regs.cnt;j++){
			r = ref_regv(g->regs, n->regs.idx + j);
			if(r->read_link == 0) continue;
			if(n->regs.idx + j == g->reads->buffer[r->rid].regs.idx) continue;
			cnt ++;
		}
		if(cnt < g->min_edge_cov){
			n->closed = 1;
			ret ++;
		}
	}
	return ret;
}

void cut_edge_core_graph(Graph *g, edge_t *e, int closed_val){
	//if(e->closed == closed_val) return;
	if(e->closed) return;
	e->closed = closed_val;
	ref_nodev(g->nodes, e->node1)->edges[e->dir1].cnt --;
	ref_nodev(g->nodes, e->node2)->edges[!e->dir2].cnt --;
}

#define cut_edge_graph(g, e) cut_edge_core_graph(g, e, 1)

void cut_lnk_core_graph(Graph *g, lnk_t *e, int closed_val){
	if(e->closed) return;
	e->closed = closed_val;
	ref_frgv(g->frgs, e->frg1)->lnks[e->dir1].cnt --;
	ref_frgv(g->frgs, e->frg2)->lnks[!e->dir2].cnt --;
}

#define cut_lnk_graph(g, e) cut_lnk_core_graph(g, e, 1)

void revive_edge_graph(Graph *g, edge_t *e){
	if(e->closed == WT_EDGE_CLOSED_NULL) return;
	e->closed = WT_EDGE_CLOSED_NULL;
	ref_nodev(g->nodes, e->node1)->edges[e->dir1].cnt ++;
	ref_nodev(g->nodes, e->node2)->edges[!e->dir2].cnt ++;
}

void revive_lnk_graph(Graph *g, lnk_t *e){
	if(e->closed == WT_EDGE_CLOSED_NULL) return;
	e->closed = WT_EDGE_CLOSED_NULL;
	ref_frgv(g->frgs, e->frg1)->lnks[e->dir1].cnt ++;
	ref_frgv(g->frgs, e->frg2)->lnks[!e->dir2].cnt ++;
}

void print_node_edges_graph(Graph *g, u8i nid, int dir, FILE *out){
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	uint64_t idx;
	n = ref_nodev(g->nodes, nid);
	idx = n->edges[dir].idx;
	while(idx){
		f = ref_edgerefv(g->erefs, idx);
		idx = f->next;
		e = ref_edgev(g->edges, f->idx);
		if(f->flg){
			fprintf(out, "N%llu\t%c\tN%llu\t%c\t%d\t%d\n", (u8i)e->node2, "+-"[!e->dir1], (u8i)e->node1, "+-"[e->dir1], e->cov, e->off);
		} else {
			fprintf(out, "N%llu\t%c\tN%llu\t%c\t%d\t%d\n", (u8i)e->node1, "+-"[e->dir1], (u8i)e->node2, "+-"[e->dir2], e->cov, e->off);
		}
	}
}

edge_ref_t* first_living_edge_graph(Graph *g, node_t *n, int dir, int *info){
	edge_ref_t *f, *ret;
	uint64_t idx;
	ret = NULL;
	if(info){
		*info = WT_TRACE_MSG_ZERO;
		if(n->edges[dir].cnt == 0) return NULL;
		idx = n->edges[dir].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			if(g->edges->buffer[f->idx].closed) continue;
			if(ret){ *info = WT_TRACE_MSG_MORE; return NULL; }
			else { *info = WT_TRACE_MSG_ONE; ret = f; }
		}
	} else {
		if(n->edges[dir].cnt == 0) return NULL;
		idx = n->edges[dir].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			if(g->edges->buffer[f->idx].closed) continue;
			if(ret){ return NULL; }
			else { ret = f; }
		}
	}
	return ret;
}

edge_ref_t* first_living_lnk_graph(Graph *g, frg_t *n, int dir, int *info){
	edge_ref_t *f, *ret;
	uint64_t idx;
	ret = NULL;
	if(info){
		*info = WT_TRACE_MSG_ZERO;
		if(n->lnks[dir].cnt == 0) return NULL;
		idx = n->lnks[dir].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			if(g->lnks->buffer[f->idx].closed) continue;
			if(ret){ *info = WT_TRACE_MSG_MORE; return NULL; }
			else { *info = WT_TRACE_MSG_ONE; ret = f; }
		}
	} else {
		if(n->lnks[dir].cnt == 0) return NULL;
		idx = n->lnks[dir].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			if(g->lnks->buffer[f->idx].closed) continue;
			if(ret){ return NULL; }
			else { ret = f; }
		}
	}
	return ret;
}

#define count_living_edges_graph(g, n, dir) (n)->edges[dir].cnt

#define count_living_lnks_graph(g, n, dir) (n)->lnks[dir].cnt

// dir = 2 means either strand
edge_ref_t* edge_node2node_graph(Graph *g, u8i node1, int dir1, u8i node2, int dir2){
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	uint64_t idx;
	int dire;
	n = ref_nodev(g->nodes, node1);
	if(dir1 > 1){
		dir1 = 0; dire = 2;
	} else {
		dire = dir1 + 1;
	}
	while(dir1 < dire){
		idx = n->edges[dir1].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = ref_edgev(g->edges, f->idx);
			if(f->flg){
				if(e->node1 == node2 && (dir2 > 1? 1 : (dir2 == (!e->dir1)))) return f;
			} else {
				if(e->node2 == node2 && (dir2 > 1? 1 : (dir2 == e->dir2))) return f;
			}
		}
		dir1 ++;
	}
	return NULL;
}

uint64_t linear_trace_graph(Graph *g, tracev *path, uint64_t max_step, int *msg){
	trace_t *t;
	node_t *n;
	edge_t *e;
	edge_ref_t *f;
	uint64_t step;
	int dir, info;
	if(path->size == 0){
		if(msg) *msg = 3;
		return 0;
	}
	t = ref_tracev(path, path->size - 1);
	step = 0;
	while(step < max_step){
		f = first_living_edge_graph(g, ref_nodev(g->nodes, t->node), t->dir, &info);
		if(info != WT_TRACE_MSG_ONE){ if(msg) *msg = info; break; }
		e = g->edges->buffer + f->idx;
		n = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
		dir = f->flg? !e->dir1 : e->dir2;
		t->edges[t->dir] = *f;
		t = next_ref_tracev(path);
		t->node = n - g->nodes->buffer;
		t->dir = dir;
		t->edges[!t->dir] = (edge_ref_t){f->idx, !f->flg, 0};
		t->edges[t->dir] = EDGE_REF_NULL;
		f = first_living_edge_graph(g, n, !dir, &info);
		if(info != WT_TRACE_MSG_ONE){  path->size --; if(msg) *msg = -1 - info; break; }
		step ++;
	}
	return step;
}

uint64_t linear_path_graph(Graph *g, pathv *path, int max_len, int *msg){
	path_t *t;
	frg_t *n;
	lnk_t *e;
	edge_ref_t *f;
	int len;
	int dir, info;
	if(path->size == 0){
		if(msg) *msg = 3;
		return 0;
	}
	t = ref_pathv(path, path->size - 1);
	len = ref_frgv(g->frgs, t->frg)->len;
	while(len < max_len){
		f = first_living_lnk_graph(g, ref_frgv(g->frgs, t->frg), t->dir, &info);
		if(info != WT_TRACE_MSG_ONE){ if(msg) *msg = info; break; }
		e = g->lnks->buffer + f->idx;
		len += e->off;
		n = ref_frgv(g->frgs, f->flg? e->frg1 : e->frg2);
		dir = f->flg? !e->dir1 : e->dir2;
		t->lnks[t->dir] = *f;
		t = next_ref_pathv(path);
		t->frg = n - g->frgs->buffer;
		t->dir = dir;
		t->lnks[!t->dir] = (edge_ref_t){f->idx, !f->flg, 0};
		t->lnks[t->dir] = EDGE_REF_NULL;
		f = first_living_lnk_graph(g, n, !dir, &info);
		if(info != WT_TRACE_MSG_ONE){  path->size --; if(msg) *msg = -1 - info; break; }
		len += n->len;
	}
	return len;
}

int cal_offset_traces_graph(Graph *g, tracev *path, u8i beg, u8i end){
	trace_t *t;
	node_t *n;
	reg_t *r;
	edge_t *e;
	uint64_t i;
	int off;
	off = 0;
	for(i=beg;i<end;i++){
		t = ref_tracev(path, i);
		t->off = off;
		n = ref_nodev(g->nodes, t->node);
		r = ref_regv(g->regs, n->regs.idx);
		if(t->edges[t->dir].idx == EDGE_REF_NULL.idx){
			off += r->end - r->beg;
		} else {
			e = ref_edgev(g->edges, t->edges[t->dir].idx);
			off += r->end - r->beg + e->off;
		}
	}
	return off;
}

int cal_offset_paths_graph(Graph *g, pathv *path, u8i beg, u8i end){
	path_t *t;
	frg_t *n;
	lnk_t *e;
	uint64_t i;
	int off;
	off = 0;
	for(i=beg;i<end;i++){
		t = ref_pathv(path, i);
		t->off = off;
		n = ref_frgv(g->frgs, t->frg);
		if(t->lnks[t->dir].idx == EDGE_REF_NULL.idx){
			off += n->length;
		} else {
			e = ref_lnkv(g->lnks, t->lnks[t->dir].idx);
			off += ((i == beg)? n->length : n->len) + e->off;
		}
	}
	return off;
}

uint64_t true_linear_unique_trace_graph(Graph *g, tracev *path, uint64_t max_step, uint64_t visit, int *msg){
	trace_t *t;
	node_t *n;
	edge_t *e;
	edge_ref_t *f;
	uint64_t step;
	int dir, info;
	if(path->size == 0){
		if(msg) *msg = WT_TRACE_MSG_ZERO;
		return 0;
	}
	step = 0;
	while(step < max_step){
		t = ref_tracev(path, path->size - 1);
		f = first_living_edge_graph(g, ref_nodev(g->nodes, t->node), !t->dir, &info);
		if(info == WT_TRACE_MSG_MORE){ if(path->size > 1) path->size --; if(msg) *msg = -1 - info; break; }
		n = ref_nodev(g->nodes, t->node);
		n->bt_visit = visit;
		f = first_living_edge_graph(g, ref_nodev(g->nodes, t->node),  t->dir, &info);
		if(info == WT_TRACE_MSG_ZERO){ if(msg) *msg = info; break; }
		else if(info == WT_TRACE_MSG_MORE){ if(path->size > 1) path->size --; if(msg) *msg = info; break; }
		e = g->edges->buffer + f->idx;
		n = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
		if(n->bt_visit == visit){ if(msg) *msg = WT_TRACE_MSG_VISITED; break; }
		dir = f->flg? !e->dir1 : e->dir2;
		t->edges[t->dir] = *f;
		t = next_ref_tracev(path);
		t->node = n - g->nodes->buffer;
		t->dir = dir;
		t->edges[!t->dir] = (edge_ref_t){f->idx, !f->flg, 0};
		t->edges[t->dir] = EDGE_REF_NULL;
		step ++;
	}
	return step;
}

uint64_t true_linear_unique_path_graph(Graph *g, pathv *path, uint64_t max_step, uint64_t visit, int *msg){
	path_t *t;
	frg_t *n;
	lnk_t *e;
	edge_ref_t *f;
	uint64_t step;
	int dir, info;
	if(path->size == 0){
		if(msg) *msg = WT_TRACE_MSG_ZERO;
		return 0;
	}
	step = 0;
	while(step < max_step){
		t = ref_pathv(path, path->size - 1);
		f = first_living_lnk_graph(g, ref_frgv(g->frgs, t->frg), !t->dir, &info);
		if(info == WT_TRACE_MSG_MORE){ if(path->size > 1) path->size --; if(msg) *msg = -1 - info; break; }
		n = ref_frgv(g->frgs, t->frg);
		n->bt_visit = visit;
		f = first_living_lnk_graph(g, ref_frgv(g->frgs, t->frg),  t->dir, &info);
		if(info == WT_TRACE_MSG_ZERO){ if(msg) *msg = info; break; }
		else if(info == WT_TRACE_MSG_MORE){ if(path->size > 1) path->size --; if(msg) *msg = info; break; }
		e = g->lnks->buffer + f->idx;
		n = ref_frgv(g->frgs, f->flg? e->frg1 : e->frg2);
		if(n->bt_visit == visit){ if(msg) *msg = WT_TRACE_MSG_VISITED; break; }
		dir = f->flg? !e->dir1 : e->dir2;
		t->lnks[t->dir] = *f;
		t = next_ref_pathv(path);
		t->frg = n - g->frgs->buffer;
		t->dir = dir;
		t->lnks[!t->dir] = (edge_ref_t){f->idx, !f->flg, 0};
		t->lnks[t->dir] = EDGE_REF_NULL;
		step ++;
	}
	return step;
}

uint64_t count_linear_trace_graph(Graph *g, trace_t *t, uint64_t max_step, int *msg){
	node_t *n;
	edge_t *e;
	edge_ref_t *f;
	uint64_t step;
	int dir, info;
	step = 0;
	while(step < max_step){
		f = first_living_edge_graph(g, ref_nodev(g->nodes, t->node), t->dir, &info);
		if(info != WT_TRACE_MSG_ONE){ if(msg) *msg = info; break; }
		e = g->edges->buffer + f->idx;
		n = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
		dir = f->flg? !e->dir1 : e->dir2;
		t->node = n - g->nodes->buffer;
		t->dir = dir;
		f = first_living_edge_graph(g, n, !dir, &info);
		if(info != WT_TRACE_MSG_ONE){ if(msg) *msg = -1 - info; break; }
		step ++;
	}
	return step;
}

int count_linear_path_graph(Graph *g, path_t *t, int max_len, int *msg){
	frg_t *n;
	lnk_t *e;
	edge_ref_t *f;
	int len;
	int dir, info;
	len = ref_frgv(g->frgs, t->frg)->len;
	while(len < max_len){
		f = first_living_lnk_graph(g, ref_frgv(g->frgs, t->frg), t->dir, &info);
		if(info != WT_TRACE_MSG_ONE){ if(msg) *msg = info; break; }
		e = g->lnks->buffer + f->idx;
		len += e->off;
		n = ref_frgv(g->frgs, f->flg? e->frg1 : e->frg2);
		dir = f->flg? !e->dir1 : e->dir2;
		t->frg = n - g->frgs->buffer;
		t->dir = dir;
		f = first_living_lnk_graph(g, n, !dir, &info);
		if(info != WT_TRACE_MSG_ONE){ if(msg) *msg = -1 - info; break; }
		len += n->len;
	}
	return len;
}

void del_node_edges_graph(Graph *g, node_t *n){
	edge_ref_t *f;
	edge_t *e;
	uint64_t idx;
	uint32_t k;
	for(k=0;k<2;k++){
		idx = n->edges[k].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = g->edges->buffer + f->idx;
			cut_edge_core_graph(g, e, WT_EDGE_CLOSED_HARD);
		}
	}
}

void del_frg_lnks_graph(Graph *g, frg_t *n){
	edge_ref_t *f;
	lnk_t *e;
	uint64_t idx;
	uint32_t k;
	for(k=0;k<2;k++){
		idx = n->lnks[k].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			e = g->lnks->buffer + f->idx;
			cut_lnk_core_graph(g, e, WT_EDGE_CLOSED_HARD);
		}
	}
}

void del_node_graph(Graph *g, node_t *n){
	del_node_edges_graph(g, n);
	n->closed = 1;
}

uint64_t del_isolated_nodes_graph(Graph *g){
	node_t *n;
	uint64_t ret, i;
	int f, r;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		first_living_edge_graph(g, n, 0, &f);
		if(f != WT_TRACE_MSG_ZERO) continue;
		first_living_edge_graph(g, n, 1, &r);
		if(r != WT_TRACE_MSG_ZERO) continue;
		n->closed = 1;
		ret ++;
	}
	return ret;
}

uint64_t cut_binary_edges_graph(Graph *g){
	UUhash *hash;
	UUhash_t *u;
	node_t *n;
	edge_ref_t *f;
	edge_t *e, *p;
	uint64_t idx, nid, ret;
	ret = 0;
	hash = init_UUhash(15);
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		clear_UUhash(hash);
		idx = n->edges[0].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = ref_edgev(g->edges, f->idx);
			if(e->closed < WT_EDGE_CLOSED_LESS){
				kv_put_UUhash(hash, f->flg? e->node1 : e->node2, f->idx);
			}
		}
		idx = n->edges[1].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = ref_edgev(g->edges, f->idx);
			if(e->closed >= WT_EDGE_CLOSED_LESS) continue;
			if((u = get_UUhash(hash, (UUhash_t){f->flg? e->node1 : e->node2, 0})) == NULL) continue;
			p = ref_edgev(g->edges, u->val);
			if(0){
				if(p->cov > e->cov) cut_edge_core_graph(g, e, WT_EDGE_CLOSED_HARD);
				else cut_edge_core_graph(g, p, WT_EDGE_CLOSED_HARD);
				ret ++;
			} else {
				cut_edge_core_graph(g, e, WT_EDGE_CLOSED_HARD);
				cut_edge_core_graph(g, p, WT_EDGE_CLOSED_HARD);
				ret += 2;
			}
		}
	}
	free_UUhash(hash);
	return ret;
}

uint64_t cut_binary_lnks_graph(Graph *g){
	UUhash *hash;
	UUhash_t *u;
	frg_t *n;
	edge_ref_t *f;
	lnk_t *e, *p;
	uint64_t idx, nid, ret;
	ret = 0;
	hash = init_UUhash(15);
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		if(n->closed) continue;
		clear_UUhash(hash);
		idx = n->lnks[0].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			e = ref_lnkv(g->lnks, f->idx);
			if(e->closed < WT_EDGE_CLOSED_LESS){
				kv_put_UUhash(hash, f->flg? e->frg1 : e->frg2, f->idx);
			}
		}
		idx = n->lnks[1].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			e = ref_lnkv(g->lnks, f->idx);
			if(e->closed >= WT_EDGE_CLOSED_LESS) continue;
			if((u = get_UUhash(hash, (UUhash_t){f->flg? e->frg1 : e->frg2, 0})) == NULL) continue;
			p = ref_lnkv(g->lnks, u->val);
			if(1){
				if(p->cov > e->cov) cut_lnk_core_graph(g, e, WT_EDGE_CLOSED_HARD);
				else cut_lnk_core_graph(g, p, WT_EDGE_CLOSED_HARD);
				ret ++;
			} else {
				cut_lnk_core_graph(g, e, WT_EDGE_CLOSED_HARD);
				cut_lnk_core_graph(g, p, WT_EDGE_CLOSED_HARD);
				ret += 2;
			}
		}
	}
	free_UUhash(hash);
	return ret;
}

u8i cut_low_cov_lnks_graph(Graph *g, int low_cov){
	frg_t *n;
	edge_ref_t *f;
	lnk_t *e;
	u8i idx, nid, ret;
	u4i k;
	int max_cov;
	ret = 0;
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		if(n->closed) continue;
		for(k=0;k<2;k++){
			max_cov = 0;
			idx = n->lnks[k].idx;
			while(idx){
				f = ref_edgerefv(g->lrefs, idx);
				idx = f->next;
				e = ref_lnkv(g->lnks, f->idx);
				if(e->cov > max_cov) max_cov = e->cov;
			}
			if(max_cov <= low_cov || max_cov < (int)g->min_edge_cov) continue;
			idx = n->lnks[k].idx;
			while(idx){
				f = ref_edgerefv(g->lrefs, idx);
				idx = f->next;
				e = ref_lnkv(g->lnks, f->idx);
				if(e->cov <= low_cov){
					cut_lnk_core_graph(g, e, WT_EDGE_CLOSED_HARD);
					ret ++;
				}
			}
		}
	}
	return ret;
}

uint32_t rescue_low_cov_transitive_edges_graph(Graph *g, uint64_t nid, u8v *edges, UUhash *hash){
	node_t *n, *w, *v;
	edge_ref_t *f, *f2, *f3;
	edge_t *e, *e1, *e2;
	reg_t *r;
	r = NULL;
	reg_t *wushigang = r;
	reg_t *tmp = wushigang;
	wushigang = tmp;
	UUhash_t *u;
	uint64_t idx, nid2, nid3;
	uint32_t i, k, k2, k3, k4, ret;
	int off1, off2, yes;
	n = ref_nodev(g->nodes, nid);
	if(n->closed) return 0;
	ret = 0;
	for(k=0;k<2;k++){
		clear_UUhash(hash);
		clear_u8v(edges);
		idx = n->edges[k].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			e = ref_edgev(g->edges, f->idx);
			if(e->closed < WT_EDGE_CLOSED_LESS){
				kv_put_UUhash(hash, f->flg? e->node1: e->node2, idx);
			} else if(e->closed == WT_EDGE_CLOSED_LESS){
				push_u8v(edges, idx);
			}
			idx = f->next;
		}
		if(edges->size <= 1) continue;
		for(i=0;i<edges->size;i++){
			f = ref_edgerefv(g->erefs, edges->buffer[i]);
			e = ref_edgev(g->edges, f->idx);
			if(e->status) continue;
			if(e->closed != WT_EDGE_CLOSED_LESS) continue;
			if(f->flg){ nid2 = e->node1; k2 = !e->dir1; }
			else      { nid2 = e->node2; k2 =  e->dir2; }
			yes = 0;
			w = ref_nodev(g->nodes, nid2);
			idx = w->edges[k2].idx;
			while(idx){
				f2 = ref_edgerefv(g->erefs, idx);
				idx = f2->next;
				if(f2->flg){ nid3 = g->edges->buffer[f2->idx].node1; k3 = !g->edges->buffer[f2->idx].dir1; }
				else       { nid3 = g->edges->buffer[f2->idx].node2; k3 =  g->edges->buffer[f2->idx].dir2; }
				e2 = ref_edgev(g->edges, f2->idx);
				if(e2->closed >= WT_EDGE_CLOSED_LESS) continue;
				if((u = get_UUhash(hash, (UUhash_t){nid3, 0})) == NULL) continue;
				f3 = ref_edgerefv(g->erefs, u->val);
				e1 = ref_edgev(g->edges, f3->idx);
				k4 = f3->flg? !e1->dir1 : e1->dir2;
				if(k3 != k4) continue;
				v = ref_nodev(g->nodes, nid3);
				off1 = off2 = (g->regs->buffer[n->regs.idx].end - g->regs->buffer[n->regs.idx].beg) + (g->regs->buffer[v->regs.idx].end - g->regs->buffer[v->regs.idx].beg);
				off1 += e1->off;
				off2 += e->off + e2->off + (g->regs->buffer[w->regs.idx].end - g->regs->buffer[w->regs.idx].beg);
				r = ref_regv(g->regs, w->regs.idx);
				//if(e->off + e2->off + (r->end - r->beg) >= longest) continue;
				if(num_diff(off1, off2) >= num_min(off1, off2)) continue;
				yes = 1;
				//revive_edge_graph(g, e);
				e->status = 1;
				ret ++;
				break;
			}
			if(0){
				if(yes) continue;
				idx = w->edges[!k2].idx;
				while(idx){
					f2 = ref_edgerefv(g->erefs, idx);
					idx = f2->next;
					if(f2->flg){ nid3 = g->edges->buffer[f2->idx].node1; k3 = !g->edges->buffer[f2->idx].dir1; }
					else       { nid3 = g->edges->buffer[f2->idx].node2; k3 =  g->edges->buffer[f2->idx].dir2; }
					e2 = ref_edgev(g->edges, f2->idx);
					if(e2->closed >= WT_EDGE_CLOSED_LESS) continue;
					if((u = get_UUhash(hash, (UUhash_t){nid3, 0})) == NULL) continue;
					f3 = ref_edgerefv(g->erefs, u->val);
					e1 = ref_edgev(g->edges, f3->idx);
					k4 = f3->flg? !e1->dir1 : e1->dir2;
					if(k3 != k4) continue;
					v = ref_nodev(g->nodes, nid3);
					off1 = off2 = (g->regs->buffer[n->regs.idx].end - g->regs->buffer[n->regs.idx].beg) + (g->regs->buffer[w->regs.idx].end - g->regs->buffer[w->regs.idx].beg);
					off1 += e->off;
					off2 += e1->off + e2->off + (g->regs->buffer[v->regs.idx].end - g->regs->buffer[v->regs.idx].beg);
					r = ref_regv(g->regs, w->regs.idx);
					if(num_diff(off1, off2) >= num_min(off1, off2)) continue;
					yes = 1;
					//revive_edge_graph(g, e);
					e->status = 1;
					ret ++;
					break;
				}
			}
		}
	}
	return ret;
}

uint64_t rescue_low_cov_edges_graph(Graph *g){
	u8v *edges;
	UUhash *hash;
	uint64_t i, nid, ret;
	ret = 0;
	edges = init_u8v(32);
	hash = init_UUhash(13);
	for(i=0;i<g->edges->size;i++) g->edges->buffer[i].status = 0;
	for(nid=0;nid<g->nodes->size;nid++){
		rescue_low_cov_transitive_edges_graph(g, nid, edges, hash);
	}
	for(i=0;i<g->edges->size;i++){
		if(g->edges->buffer[i].status){
			revive_edge_graph(g, g->edges->buffer + i);
			g->edges->buffer[i].status = 0;
			ret ++;
		}
	}
	free_u8v(edges);
	free_UUhash(hash);
	return ret;
}

uint32_t rescue_low_cov_tip_edges_core(Graph *g, uint64_t nid){
	node_t *n, *w, *ww;
	edge_t *e, *ee;
	edge_ref_t *f;
	uint64_t idx, wid;
	uint32_t k, dir, ret;
	n = ref_nodev(g->nodes, nid);
	if(n->closed) return 0;
	if(n->edges[0].cnt == 0 && n->edges[1].cnt == 0) return 0;
	ret = 0;
	for(k=0;k<2;k++){
		if(n->edges[k].cnt) continue;
		idx = n->edges[k].idx;
		ee = NULL;
		ww = NULL;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = ref_edgev(g->edges, f->idx);
			if(f->flg){
				wid = e->node1; dir = !e->dir1;
			} else {
				wid = e->node2; dir = e->dir2;
			}
			w = ref_nodev(g->nodes, wid);
			//if(w->edges[!dir].cnt) continue;
			if(w->edges[dir].cnt == 0) continue;
			if(ee == NULL || e->cov > ee->cov || (e->cov == ee->cov && w->regs.cnt > ww->regs.cnt)){ ee = e; ww = w; }
		}
		if(ee){ revive_edge_graph(g, ee); ret ++; }
	}
	return ret;
}

uint64_t rescue_low_cov_tip_edges_graph(Graph *g){
	u8v *edges;
	UUhash *hash;
	uint64_t nid, ret;
	ret = 0;
	edges = init_u8v(32);
	hash = init_UUhash(13);
	for(nid=0;nid<g->nodes->size;nid++){
		ret += rescue_low_cov_tip_edges_core(g, nid);
	}
	free_u8v(edges);
	free_UUhash(hash);
	return ret;
}

uint32_t rescue_weak_tip_lnks_core(Graph *g, uint64_t nid){
	frg_t *n, *w, *ww;
	ww = NULL;
	frg_t *wushigang = ww;
	frg_t *tmp = wushigang;
	wushigang = tmp;
	lnk_t *e, *ee;
	edge_ref_t *f;
	uint64_t idx, wid;
	uint32_t k, dir, ret;
	n = ref_frgv(g->frgs, nid);
	if(n->closed) return 0;
	if(n->lnks[0].cnt == 0 && n->lnks[1].cnt == 0) return 0;
	ret = 0;
	for(k=0;k<2;k++){
		if(n->lnks[k].cnt) continue;
		idx = n->lnks[k].idx;
		ee = NULL;
		ww = NULL;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			e = ref_lnkv(g->lnks, f->idx);
			if(e->weak == 0) continue;
			if(f->flg){
				wid = e->frg1; dir = !e->dir1;
			} else {
				wid = e->frg2; dir = e->dir2;
			}
			w = ref_frgv(g->frgs, wid);
			if(w->lnks[!dir].cnt) continue;
			if(ee == NULL){ ee = e; ww = w; }
			else { ee = NULL; break; }
		}
		if(ee){ revive_lnk_graph(g, ee); ret ++; }
	}
	return ret;
}

uint64_t rescue_weak_tip_lnks2_graph(Graph *g){
	uint64_t nid, ret;
	ret = 0;
	for(nid=0;nid<g->frgs->size;nid++){
		ret += rescue_weak_tip_lnks_core(g, nid);
	}
	return ret;
}

u8i rescue_weak_tip_lnks_graph(Graph *g){
	u8v *weaks[2];
	frg_t *n;
	lnk_t *e;
	edge_ref_t *f;
	u8i nid, i, ret;
	u4i k, eidx, idx;
	weaks[0] = init_u8v(g->frgs->size);
	weaks[1] = init_u8v(g->frgs->size);
	ret = 0;
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		for(k=0;k<2;k++){
			if(n->lnks[k].cnt){
				eidx = 0;
			} else {
				idx = n->lnks[k].idx;
				eidx = 0;
				while(idx){
					f = ref_edgerefv(g->lrefs, idx);
					e = ref_lnkv(g->lnks, f->idx);
					if(e->weak){
						if(eidx == 0) eidx = f->idx;
						else eidx = MAX_VALUE_U4;
					}
					idx = f->next;
				}
			}
			push_u8v(weaks[k], eidx);
		}
	}
	for(k=0;k<2;k++){
		for(i=0;i<weaks[k]->size;i++){
			if(weaks[k]->buffer[i] == 0 || weaks[k]->buffer[i] == MAX_VALUE_U4) continue;
			e = ref_lnkv(g->lnks, weaks[k]->buffer[i]);
			if(i != e->frg1) continue;
			if(weaks[0]->buffer[e->frg2] == weaks[k]->buffer[i] || weaks[1]->buffer[e->frg2] == weaks[k]->buffer[i]){
				ret ++;
				revive_lnk_graph(g, e);
			}
		}
	}
	free_u8v(weaks[0]);
	free_u8v(weaks[1]);
	return ret;
}

int _scoring_edge_orders(Graph *g, uint64_t fidx){
	edge_ref_t *f;
	edge_t *e;
	node_t *n;
	int score;
	f = ref_edgerefv(g->erefs, fidx);
	e = ref_edgev(g->edges, f->idx);
	n = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
	score = e->off + (n->regs.cnt * -5) + (e->cov * -5);
	return score;
}

uint32_t reduce_transitive_edges_core_graph(Graph *g, uint64_t nid, u8v *edges, UUhash *hash, uint32_t closed_val){
	node_t *n, *w, *v;
	edge_ref_t *f, *f2, *f3;
	edge_t *e, *e1, *e2;
	UUhash_t *u;
	uint64_t idx, nid2, nid3;
	uint32_t i, k, k2, k3, k4, ret;
	int off1, off2;
	n = ref_nodev(g->nodes, nid);
	if(n->closed) return 0;
	ret = 0;
	for(k=0;k<2;k++){
		clear_UUhash(hash);
		clear_u8v(edges);
		idx = n->edges[k].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			e = ref_edgev(g->edges, f->idx);
			if(e->closed < closed_val){
				kv_put_UUhash(hash, f->flg? e->node1: e->node2, (idx << 1) | e->closed);
				push_u8v(edges, idx);
			}
			idx = f->next;
		}
		if(edges->size <= 1) continue;
		// sort the edges by composition of e->off and e->cov
		//sort_array(edges->buffer, edges->size, uint64_t, num_cmpgt(_scoring_edge_orders(g, a), _scoring_edge_orders(g, b)));
		for(i=0;i<edges->size;i++){
			f = ref_edgerefv(g->erefs, edges->buffer[i]);
			e = ref_edgev(g->edges, f->idx);
			if(f->flg){ nid2 = e->node1; k2 = !e->dir1; }
			else      { nid2 = e->node2; k2 =  e->dir2; }
			w = ref_nodev(g->nodes, nid2);
			idx = w->edges[k2].idx;
			while(idx){
				f2 = ref_edgerefv(g->erefs, idx);
				idx = f2->next;
				if(f2->flg){ nid3 = g->edges->buffer[f2->idx].node1; k3 = !g->edges->buffer[f2->idx].dir1; }
				else       { nid3 = g->edges->buffer[f2->idx].node2; k3 =  g->edges->buffer[f2->idx].dir2; }
				e2 = ref_edgev(g->edges, f2->idx);
				//if(e2->closed) continue;
				if(e2->closed >= closed_val) continue;
				if((u = get_UUhash(hash, (UUhash_t){nid3, 0})) == NULL) continue;
				if(u->val & 0x01) continue; // already deleted
				f3 = ref_edgerefv(g->erefs, u->val >> 1);
				e1 = ref_edgev(g->edges, f3->idx);
				k4 = f3->flg? !e1->dir1 : e1->dir2;
				if(k3 != k4) continue;
				v = ref_nodev(g->nodes, nid3);
				off1 = off2 = (g->regs->buffer[n->regs.idx].end - g->regs->buffer[n->regs.idx].beg) + (g->regs->buffer[v->regs.idx].end - g->regs->buffer[v->regs.idx].beg);
				off1 += e1->off;
				off2 += e->off + e2->off + (g->regs->buffer[w->regs.idx].end - g->regs->buffer[w->regs.idx].beg);
				// check whether off1 and off2 diff too much
				if(num_diff(off1, off2) >= num_min(off1, off2)) continue;
				u->val |= 1;
			}
		}
		reset_iter_UUhash(hash);
		while((u = ref_iter_UUhash(hash))){
			if(u->val & 0x01){
				e = ref_edgev(g->edges, ref_edgerefv(g->erefs, u->val >> 1)->idx);
				if(e->closed == WT_EDGE_CLOSED_NULL){
					cut_edge_graph(g, e);
					ret ++;
				}
			}
		}
	}
	return ret;
}

void set_init_ends_graph(Graph *g){
	node_t *n;
	u8i nid;
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->edges[0].cnt == 0 || n->edges[1].cnt == 0){
			n->init_end = 1;
		}
	}
}

uint64_t reduce_transitive_edges_graph(Graph *g){
	u8v *edges;
	UUhash *hash;
	uint64_t nid, ret;
	ret = 0;
	edges = init_u8v(32);
	hash = init_UUhash(13);
	for(nid=0;nid<g->nodes->size;nid++){
		ret += reduce_transitive_edges_core_graph(g, nid, edges, hash, 2);
	}
	free_u8v(edges);
	free_UUhash(hash);
	return ret;
}

uint32_t reduce_transitive_lnks_core_graph(Graph *g, uint64_t nid, u8v *lnks, UUhash *hash, uint32_t closed_val){
	frg_t *n, *w, *v;
	edge_ref_t *f, *f2, *f3;
	lnk_t *e, *e1, *e2;
	UUhash_t *u;
	uint64_t idx, nid2, nid3;
	uint32_t i, k, k2, k3, k4, ret;
	int off1, off2;
	n = ref_frgv(g->frgs, nid);
	if(n->closed) return 0;
	ret = 0;
	for(k=0;k<2;k++){
		clear_UUhash(hash);
		clear_u8v(lnks);
		idx = n->lnks[k].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			e = ref_lnkv(g->lnks, f->idx);
			if(e->closed < closed_val){
				kv_put_UUhash(hash, f->flg? e->frg1: e->frg2, (idx << 1) | e->closed);
				push_u8v(lnks, idx);
			}
			idx = f->next;
		}
		if(lnks->size <= 1) continue;
		// sort the lnks by composition of e->off and e->cov
		//sort_array(lnks->buffer, lnks->size, uint64_t, num_cmpgt(_scoring_lnk_orders(g, a), _scoring_lnk_orders(g, b)));
		for(i=0;i<lnks->size;i++){
			f = ref_edgerefv(g->lrefs, lnks->buffer[i]);
			e = ref_lnkv(g->lnks, f->idx);
			if(f->flg){ nid2 = e->frg1; k2 = !e->dir1; }
			else      { nid2 = e->frg2; k2 =  e->dir2; }
			w = ref_frgv(g->frgs, nid2);
			idx = w->lnks[k2].idx;
			while(idx){
				f2 = ref_edgerefv(g->lrefs, idx);
				idx = f2->next;
				if(f2->flg){ nid3 = g->lnks->buffer[f2->idx].frg1; k3 = !g->lnks->buffer[f2->idx].dir1; }
				else       { nid3 = g->lnks->buffer[f2->idx].frg2; k3 =  g->lnks->buffer[f2->idx].dir2; }
				e2 = ref_lnkv(g->lnks, f2->idx);
				//if(e2->closed) continue;
				if(e2->closed >= closed_val) continue;
				if((u = get_UUhash(hash, (UUhash_t){nid3, 0})) == NULL) continue;
				if(u->val & 0x01) continue; // already deleted
				f3 = ref_edgerefv(g->lrefs, u->val >> 1);
				e1 = ref_lnkv(g->lnks, f3->idx);
				k4 = f3->flg? !e1->dir1 : e1->dir2;
				if(k3 != k4) continue;
				v = ref_frgv(g->frgs, nid3);
				off1 = off2 = n->len + v->len;
				off1 += e1->off;
				off2 += e->off + e2->off + w->len;
				// check whether off1 and off2 diff too much
				if(num_diff(off1, off2) >= num_min(off1, off2)) continue;
				u->val |= 1;
			}
		}
		reset_iter_UUhash(hash);
		while((u = ref_iter_UUhash(hash))){
			if(u->val & 0x01){
				e = ref_lnkv(g->lnks, ref_edgerefv(g->lrefs, u->val >> 1)->idx);
				if(e->closed == WT_EDGE_CLOSED_NULL){
					cut_lnk_graph(g, e);
					ret ++;
				}
			}
		}
	}
	return ret;
}

uint64_t reduce_transitive_lnks_graph(Graph *g){
	u8v *lnks;
	UUhash *hash;
	uint64_t nid, ret;
	ret = 0;
	lnks = init_u8v(32);
	hash = init_UUhash(13);
	for(nid=0;nid<g->frgs->size;nid++){
		ret += reduce_transitive_lnks_core_graph(g, nid, lnks, hash, 2);
	}
	free_u8v(lnks);
	free_UUhash(hash);
	return ret;
}

uint64_t trim_tip_core_graph(Graph *g, uint16_t max_step, tracev *path, uint64_t nid, int hard_trim){
	trace_t *t, T;
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	uint64_t ret, idx;
	uint32_t i, dir, step, step2, found, n_in;
	int msg1, msg2;
	if(g->cut_tip == 0) return 0;
	ret = 0;
	n = ref_nodev(g->nodes, nid);
	if(n->closed) return 0;
	first_living_edge_graph(g, n, 0, &msg1);
	first_living_edge_graph(g, n, 1, &msg2);
	if(msg1 != WT_TRACE_MSG_ZERO){
		if(msg2 != WT_TRACE_MSG_ZERO) return 0;
		dir = 0;
	} else if(msg2 != WT_TRACE_MSG_ZERO){
		dir = 1;
	} else return 0;
	clear_tracev(path);
	t = next_ref_tracev(path);
	t->node = nid;
	t->edges[0] = EDGE_REF_NULL;
	t->edges[1] = EDGE_REF_NULL;
	t->dir = dir;
	step = linear_trace_graph(g, path, max_step, &msg1) + 1;
	if(step > max_step) return 0;
	//if(msg1 != -1 - WT_TRACE_MSG_MORE && msg1 != WT_TRACE_MSG_MORE) return 0;
	if(msg1 == WT_TRACE_MSG_MORE){
		if(!hard_trim) return 0;
		path->size --;
	} else if(msg1 == -1 - WT_TRACE_MSG_MORE){
		t = ref_tracev(path, path->size); // please see linear_trace_graph
		n = ref_nodev(g->nodes, t->node);
		dir = !t->dir;
		n_in = 0;
		idx = n->edges[dir].idx;
		found = 0;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			if(f->idx == t->edges[dir].idx) continue;
			e = g->edges->buffer + f->idx;
			if(e->closed) continue;
			n_in ++;
			T.node = f->flg? e->node1 : e->node2;
			T.dir  = f->flg? !e->dir1 : e->dir2;
			step2 = count_linear_trace_graph(g, &T, step + 1, &msg2) + 1;
			if(msg2 != WT_TRACE_MSG_ZERO) step2 ++;
			if(step2 >= step){ found = 1; break; }
			//if(step2 + 1 >= step && msg2 != WT_TRACE_MSG_ZERO){ found = 1; break; }
		}
		if(!found) return 0;
	} else return 0;
	for(i=0;i<path->size;i++){
		del_node_graph(g, ref_nodev(g->nodes, path->buffer[i].node));
		//del_node_edges_graph(g, ref_nodev(g->nodes, path->buffer[i].node));
		ret ++;
	}
	return ret;
}

uint64_t trim_tips_graph(Graph *g, uint16_t max_step, int hard_trim){
	tracev *path;
	node_t *n;
	uint64_t ret, nid;
	ret = 0;
	path = init_tracev(32);
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		if(trim_tip_core_graph(g, max_step, path, nid, hard_trim)) ret ++;
	}
	free_tracev(path);
	return ret;
}

// careful sharp_tip -> blunt_tip -> sharp_tip and so on
u4i trim_blunt_tip_core_graph(Graph *g, u8i nid){
	node_t *n;
	int k;
	if(g->cut_tip == 0) return 0;
	n = ref_nodev(g->nodes, nid);
	if(n->edges[0].cnt && n->edges[1].cnt) return 0;
	k = (n->edges[0].cnt == 0);
	if(n->edges[k].cnt < 2) return 0;
	del_node_graph(g, n);
	return 1;
}

uint64_t trim_blunt_tips_graph(Graph *g){
	node_t *n;
	uint64_t ret, nid;
	ret = 0;
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		if(trim_blunt_tip_core_graph(g, nid)) ret ++;
	}
	return ret;
}

uint64_t trim_frgtip_core_graph(Graph *g, int max_len, pathv *path, uint64_t nid){
	path_t *t, T;
	frg_t *n;
	edge_ref_t *f;
	lnk_t *e;
	uint64_t ret, idx;
	uint32_t i, dir, found, n_in;
	int len, len2;
	int msg1, msg2;
	if(g->cut_tip == 0) return 0;
	ret = 0;
	n = ref_frgv(g->frgs, nid);
	if(n->closed) return 0;
	first_living_lnk_graph(g, n, 0, &msg1);
	first_living_lnk_graph(g, n, 1, &msg2);
	if(msg1 != WT_TRACE_MSG_ZERO){
		if(msg2 != WT_TRACE_MSG_ZERO) return 0;
		dir = 0;
	} else if(msg2 != WT_TRACE_MSG_ZERO){
		dir = 1;
	} else return 0;
	clear_pathv(path);
	t = next_ref_pathv(path);
	t->frg = nid;
	t->lnks[0] = EDGE_REF_NULL;
	t->lnks[1] = EDGE_REF_NULL;
	t->dir = dir;
	len = linear_path_graph(g, path, max_len, &msg1) + 1;
	if(len > max_len) return 0;
	//if(msg1 != -1 - WT_TRACE_MSG_MORE && msg1 != WT_TRACE_MSG_MORE) return 0;
	if(msg1 == WT_TRACE_MSG_MORE){
		path->size --;
	} else if(msg1 == -1 - WT_TRACE_MSG_MORE){
		t = ref_pathv(path, path->size); // please see linear_path_graph
		n = ref_frgv(g->frgs, t->frg);
		dir = !t->dir;
		n_in = 0;
		idx = n->lnks[dir].idx;
		found = 0;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			if(f->idx == t->lnks[dir].idx) continue;
			e = g->lnks->buffer + f->idx;
			if(e->closed) continue;
			n_in ++;
			T.frg = f->flg? e->frg1 : e->frg2;
			T.dir  = f->flg? !e->dir1 : e->dir2;
			len2 = count_linear_path_graph(g, &T, len + 1, &msg2) + 1;
			if(msg2 != WT_TRACE_MSG_ZERO) len2 ++;
			if(len2 >= len){ found = 1; break; }
		}
		if(!found) return 0;
	} else return 0;
	for(i=0;i<path->size;i++){
		del_frg_lnks_graph(g, ref_frgv(g->frgs, path->buffer[i].frg));
		ret ++;
	}
	return ret;
}

uint64_t trim_frgtips_graph(Graph *g, int max_len){
	pathv *path;
	frg_t *n;
	uint64_t ret, nid;
	ret = 0;
	path = init_pathv(32);
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		if(n->closed) continue;
		if(trim_frgtip_core_graph(g, max_len, path, nid)) ret ++;
	}
	free_pathv(path);
	return ret;
}

typedef struct {
	node_t *n;
	edge_t *e; // incoming edge
	uint64_t dir:1, ind:1, step:8, bt:16, ending:16, score:20, keep:2;
} bt_t;
define_list(btv, bt_t);
#define WT_MAX_BTIDX	0xFFFF

uint32_t pop_bubble_backtrace_graph(Graph *g, btv *bts, uint32_t idx){
	bt_t *bt;
	uint32_t ret;
	while(idx){
		bt = ref_btv(bts, idx);
		bt->step = 0;
		idx = bt->bt;
	}
	ret = 0;
	for(idx=1;idx<bts->size;idx++){
		bt = ref_btv(bts, idx);
		if(bt->step == 0) continue;
		cut_edge_graph(g, bt->e);
		ret ++;
	}
	return ret;
}

typedef struct {
	frg_t *n;
	lnk_t *e; // incoming edge
	uint64_t dir:1, ind:1, step:8, bt:16, ending:16, score:22;
} frg_bt_t;
define_list(frgbtv, frg_bt_t);
uint32_t pop_frg_bubble_backtrace_graph(Graph *g, frgbtv *bts, uint32_t idx){
	frg_bt_t *bt;
	uint32_t ret;
	while(idx){
		bt = ref_frgbtv(bts, idx);
		bt->step = 0;
		idx = bt->bt;
	}
	ret = 0;
	for(idx=1;idx<bts->size;idx++){
		bt = ref_frgbtv(bts, idx);
		if(bt->step == 0) continue;
		cut_lnk_graph(g, bt->e);
		ret ++;
	}
	return ret;
}

uint32_t pop_bubble2_backtrace_graph(Graph *g, btv *bts, uint32_t _idx){
	bt_t *bt;
	uint32_t ret, i, idx;
	for(i=1;i<bts->size;i++){
		bt = ref_btv(bts, i);
		bt->keep = 2;
	}
	idx = _idx;
	while(idx){
		bt = ref_btv(bts, idx);
		bt->keep = 1;
		idx = bt->bt;
	}
	for(i=1;i<bts->size;i++){
		bt = ref_btv(bts, i);
		if(bt->keep == 1) continue;
		if(bt->ending == 0) continue;
		if(bts->buffer[bt->ending].keep == 1){
			idx = i;
			while(idx){
				bt = ref_btv(bts, idx);
				if(bt->keep != 2) break;
				idx = bt->bt;
				bt->keep = 0;
			}
		}
	}
	ret = 0;
	for(idx=1;idx<bts->size;idx++){
		bt = ref_btv(bts, idx);
		if(bt->keep) continue;
		cut_edge_graph(g, bt->e);
		ret ++;
	}
	return ret;
}

uint32_t pop_frg_bubble2_backtrace_graph(Graph *g, frgbtv *bts, uint32_t _idx){
	frg_bt_t *bt;
	uint32_t ret, i, idx;
	for(i=1;i<bts->size;i++){
		bt = ref_frgbtv(bts, i);
		if(bt->ending == _idx){
			idx = i;
			while(idx){
				bt = ref_frgbtv(bts, idx);
				bt->step = 0;
				idx = bt->bt;
			}
		}
	}
	idx = _idx;
	while(idx){
		bt = ref_frgbtv(bts, idx);
		bt->step = 1;
		idx = bt->bt;
	}
	ret = 0;
	for(idx=1;idx<bts->size;idx++){
		bt = ref_frgbtv(bts, idx);
		if(bt->step != 0) continue;
		cut_lnk_graph(g, bt->e);
		ret ++;
	}
	return ret;
}

uint32_t safe_cut_redundant_edges_graph(Graph *g, btv *bts, bt_t *b1, bt_t *b2){
	uint32_t ret;
	ret = 0;
	if(0){
		ret = 1;
		cut_edge_graph(g, b1->e);
		b1 = ref_btv(bts, b1->bt);
	}
	while(1){
		if(b1->step >= b2->step){
			if(b1 == b2) break;
			ret ++;
			cut_edge_graph(g, b1->e);
			b1 = ref_btv(bts, b1->bt);
		} else {
			b2 = ref_btv(bts, b2->bt);
		}
	}
	return ret;
}

uint32_t safe_cut_redundant_lnks_graph(Graph *g, frgbtv *bts, frg_bt_t *b1, frg_bt_t *b2){
	uint32_t ret;
	ret = 0;
	if(0){
		ret = 1;
		cut_lnk_graph(g, b1->e);
		b1 = ref_frgbtv(bts, b1->bt);
	}
	while(1){
		if(b1->step >= b2->step){
			if(b1 == b2) break;
			ret ++;
			cut_lnk_graph(g, b1->e);
			b1 = ref_frgbtv(bts, b1->bt);
		} else {
			b2 = ref_frgbtv(bts, b2->bt);
		}
	}
	return ret;
}

uint32_t pop_bubble_core_graph(Graph *g, uint16_t max_step, btv *bts, u4v *heap, uint64_t nid, uint32_t dir, uint64_t visit){
	bt_t *bt, *tb;
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	uint64_t ret, idx;
	ret = 0;
	uint64_t wushigang = ret;
	uint64_t tmp = wushigang;
	wushigang = tmp;
	uint32_t bidx, i, lst, unclosed;
	ret = 0;
	n = ref_nodev(g->nodes, nid);
	if(n->closed) return 0;
	if(count_living_edges_graph(g, n, dir) < 2) return 0;
	clear_btv(bts);
	next_ref_btv(bts);
	bt = next_ref_btv(bts);
	bt->n = n;
	bt->e = NULL;
	bt->dir = dir;
	bt->ind = 1;
	//bt->ind = 0;
	bt->step = 0;
	bt->bt = 0;
	bt->score = 0;
	bt->ending = 0;
	clear_u4v(heap);
	array_heap_push(heap->buffer, heap->size, heap->cap, uint32_t, bts->size - 1, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
	n->bt_visit = visit;
	n->bt_idx = bts->size - 1;
	n->single_in = 1;
	unclosed = 0;
	while(heap->size && heap->size < WT_MAX_BTIDX){
		bidx = array_heap_pop(heap->buffer, heap->size, heap->cap, uint32_t, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
		encap_btv(bts, bts->buffer[bidx].n->edges[bts->buffer[bidx].dir].cnt);
		bt = ref_btv(bts, bidx);
		if(bt->step >= max_step) return 0;
		if(bt->ind && bt->n->single_in == 0) bt->ind = 0;
		lst = bts->size;
		idx = bt->n->edges[bt->dir].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = g->edges->buffer + f->idx;
			if(e->closed) continue;
			tb = next_ref_btv(bts);
			tb->n = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
			if(tb->n == bts->buffer[1].n) return 0;
			tb->e = e;
			tb->dir = f->flg? !e->dir1 : e->dir2;
			tb->step = bt->step + 1;
			tb->bt   = bidx;
			tb->ind = 0;
			tb->score = bt->score + e->cov;
			tb->ending = 0;
		}
		if(bt->ind && (bt->bt == 0 || lst + 1 == bts->size)){
			for(i=lst;i<bts->size;i++) bts->buffer[i].ind = 1;
		}
		for(i=lst;i<bts->size;i++){
			tb = ref_btv(bts, i);
			if(tb->n->bt_visit != visit){
				tb->n->bt_visit = visit;
				tb->n->unvisit = count_living_edges_graph(g, tb->n, !tb->dir);
				if(tb->n->unvisit == 1) tb->n->single_in = 1;
				else tb->n->single_in = 0;
				tb->n->bt_idx = i;
				unclosed ++;
			} else {
				tb->ending = tb->n->bt_idx;
				if(tb->dir == bts->buffer[tb->n->bt_idx].dir){
					if(tb->ind && bts->buffer[tb->n->bt_idx].ind){
						if(tb->step == bts->buffer[tb->n->bt_idx].step){
							return safe_cut_redundant_edges_graph(g, bts, tb, ref_btv(bts, tb->n->bt_idx));
						} else {
							return safe_cut_redundant_edges_graph(g, bts, ref_btv(bts, tb->n->bt_idx), tb);
						}
					} else if(tb->ind){
						return safe_cut_redundant_edges_graph(g, bts, tb, ref_btv(bts, tb->n->bt_idx));
					} else if(bts->buffer[tb->n->bt_idx].ind){
						return safe_cut_redundant_edges_graph(g, bts, ref_btv(bts, tb->n->bt_idx), tb);
					}
				} else {
					// circle
					return 0;
				}
			}
			tb->n->unvisit --;
			if(tb->n->unvisit == 0){
				if(tb->step > bts->buffer[tb->n->bt_idx].step) tb->n->bt_idx = i;
				if(count_living_edges_graph(g, tb->n, tb->dir)){
					array_heap_push(heap->buffer, heap->size, heap->cap, uint32_t, tb->n->bt_idx, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
				}
				unclosed --;
			}
		}
		if(heap->size == 1 && unclosed == 0){
			return pop_bubble2_backtrace_graph(g, bts, heap->buffer[0]);
		}
	}
	return 0;
}

uint64_t pop_bubbles_graph(Graph *g, uint16_t max_step){
	btv *bts;
	u4v *heap;
	node_t *n;
	uint64_t nid, visit, ret, _ret;
	int dir;
	ret = 0;
	for(nid=0;nid<g->nodes->size;nid++) g->nodes->buffer[nid].bt_visit = 0;
	bts = init_btv(32);
	heap = init_u4v(32);
	visit = 0;
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		for(dir=0;dir<2;dir++){
			_ret = pop_bubble_core_graph(g, max_step, bts, heap, nid, dir, ++visit);
			if(_ret) ret ++;
		}
	}
	free_btv(bts);
	free_u4v(heap);
	return ret;
}

uint32_t pop_frg_bubble_core_graph(Graph *g, uint16_t max_step, frgbtv *bts, u4v *heap, uint64_t nid, uint32_t dir, uint64_t visit){
	frg_bt_t *bt, *tb;
	frg_t *n;
	edge_ref_t *f;
	lnk_t *e;
	uint64_t ret, idx;
	ret = 0;
	uint64_t wushigang = ret;
	uint64_t tmp = wushigang;
	wushigang = tmp;
	uint32_t bidx, i, lst, unclosed;
	ret = 0;
	n = ref_frgv(g->frgs, nid);
	if(n->closed) return 0;
	if(count_living_lnks_graph(g, n, dir) < 2) return 0;
	clear_frgbtv(bts);
	next_ref_frgbtv(bts);
	bt = next_ref_frgbtv(bts);
	bt->n = n;
	bt->e = NULL;
	bt->dir = dir;
	bt->ind = 1;
	//bt->ind = 0;
	bt->step = 0;
	bt->bt = 0;
	bt->score = 0;
	bt->ending = 0;
	clear_u4v(heap);
	array_heap_push(heap->buffer, heap->size, heap->cap, uint32_t, bts->size - 1, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
	n->bt_visit = visit;
	n->bt_idx = bts->size - 1;
	n->single_in = 1;
	unclosed = 0;
	while(heap->size && heap->size < WT_MAX_BTIDX){
		bidx = array_heap_pop(heap->buffer, heap->size, heap->cap, uint32_t, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
		encap_frgbtv(bts, bts->buffer[bidx].n->lnks[bts->buffer[bidx].dir].cnt);
		bt = ref_frgbtv(bts, bidx);
		if(bt->step >= max_step) return 0;
		if(bt->ind && bt->n->single_in == 0) bt->ind = 0;
		lst = bts->size;
		idx = bt->n->lnks[bt->dir].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			e = g->lnks->buffer + f->idx;
			if(e->closed) continue;
			tb = next_ref_frgbtv(bts);
			tb->n = ref_frgv(g->frgs, f->flg? e->frg1 : e->frg2);
			if(tb->n == bts->buffer[1].n) return 0;
			tb->e = e;
			tb->dir = f->flg? !e->dir1 : e->dir2;
			tb->step = bt->step + 1;
			tb->bt   = bidx;
			tb->ind = 0;
			tb->score = bt->score + e->cov;
			tb->ending = 0;
		}
		if(bt->ind && (bt->bt == 0 || lst + 1 == bts->size)){
			for(i=lst;i<bts->size;i++) bts->buffer[i].ind = 1;
		}
		for(i=lst;i<bts->size;i++){
			tb = ref_frgbtv(bts, i);
			if(tb->n->bt_visit != visit){
				tb->n->bt_visit = visit;
				tb->n->unvisit = count_living_lnks_graph(g, tb->n, !tb->dir);
				if(tb->n->unvisit == 1) tb->n->single_in = 1;
				else tb->n->single_in = 0;
				tb->n->bt_idx = i;
				unclosed ++;
			} else {
				tb->ending = tb->n->bt_idx;
				if(tb->dir == bts->buffer[tb->n->bt_idx].dir){
					if(tb->ind && bts->buffer[tb->n->bt_idx].ind){
						if(tb->step == bts->buffer[tb->n->bt_idx].step){
							return safe_cut_redundant_lnks_graph(g, bts, tb, ref_frgbtv(bts, tb->n->bt_idx));
						} else {
							return safe_cut_redundant_lnks_graph(g, bts, ref_frgbtv(bts, tb->n->bt_idx), tb);
						}
					} else if(tb->ind){
						return safe_cut_redundant_lnks_graph(g, bts, tb, ref_frgbtv(bts, tb->n->bt_idx));
					} else if(bts->buffer[tb->n->bt_idx].ind){
						return safe_cut_redundant_lnks_graph(g, bts, ref_frgbtv(bts, tb->n->bt_idx), tb);
					}
				}
			}
			tb->n->unvisit --;
			if(tb->n->unvisit == 0){
				if(tb->step > bts->buffer[tb->n->bt_idx].step) tb->n->bt_idx = i;
				if(count_living_lnks_graph(g, tb->n, tb->dir)){
					array_heap_push(heap->buffer, heap->size, heap->cap, uint32_t, tb->n->bt_idx, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
				}
				unclosed --;
			}
		}
		if(heap->size == 1 && unclosed == 0){
			return pop_frg_bubble2_backtrace_graph(g, bts, heap->buffer[0]);
		}
	}
	return 0;
}

uint64_t pop_frg_bubbles_graph(Graph *g, uint16_t max_step){
	frgbtv *bts;
	u4v *heap;
	frg_t *n;
	uint64_t nid, visit, ret, _ret;
	int dir;
	ret = 0;
	for(nid=0;nid<g->frgs->size;nid++) g->frgs->buffer[nid].bt_visit = 0;
	bts = init_frgbtv(32);
	heap = init_u4v(32);
	visit = 0;
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		if(n->closed) continue;
		for(dir=0;dir<2;dir++){
			_ret = pop_frg_bubble_core_graph(g, max_step, bts, heap, nid, dir, ++visit);
			if(_ret) ret ++;
		}
	}
	free_frgbtv(bts);
	free_u4v(heap);
	return ret;
}

u4i resolve_yarn_core_graph(Graph *g, u4i max_step, btv *bts, u4v *heap, u8i nid, u4i dir, u8i visit){
	bt_t *bt, *tb;
	node_t *n, *m;
	edge_ref_t *f;
	edge_t *e;
	uint64_t ret, idx, tip_idx;
	ret = 0;
	uint64_t wushigang = ret;
	uint64_t tmp = wushigang;
	wushigang = tmp;
	uint32_t bidx, i, lst, tip;
	ret = 0;
	n = ref_nodev(g->nodes, nid);
	if(n->closed) return 0;
	if(count_living_edges_graph(g, n, dir) < 2) return 0;
	clear_btv(bts);
	next_ref_btv(bts);
	bt = next_ref_btv(bts);
	bt->n = n;
	bt->e = NULL;
	bt->dir = dir;
	bt->step = 0;
	bt->bt = 0;
	bt->score = 0;
	bt->ending = 0;
	clear_u4v(heap);
	array_heap_push(heap->buffer, heap->size, heap->cap, uint32_t, bts->size - 1, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
	n->bt_visit = visit;
	n->bt_idx = bts->size - 1;
	n->single_in = 1;
	tip = 0; tip_idx = WT_MAX_BTIDX;
	while(heap->size && bts->size < WT_MAX_BTIDX){
		bidx = array_heap_pop(heap->buffer, heap->size, heap->cap, uint32_t, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
		encap_btv(bts, bts->buffer[bidx].n->edges[bts->buffer[bidx].dir].cnt);
		bt = ref_btv(bts, bidx);
		if(bt->step >= max_step) return 0;
		lst = bts->size;
		idx = bt->n->edges[bt->dir].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = g->edges->buffer + f->idx;
			if(e->closed) continue;
			tb = next_ref_btv(bts);
			tb->n = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
			//if(tb->n == bts->buffer[1].n) return 0;
			tb->e = e;
			tb->dir = f->flg? !e->dir1 : e->dir2;
			tb->step = bt->step + 1;
			tb->bt   = bidx;
			tb->score = bt->score + e->cov;
			tb->ending = 0;
		}
		for(i=lst;i<bts->size;i++){
			tb = ref_btv(bts, i);
			if(tb->n->bt_visit != visit){
				tb->n->bt_visit = visit;
				tb->n->bt_idx = i;
				if(count_living_edges_graph(g, tb->n, tb->dir)){
					array_heap_push(heap->buffer, heap->size, heap->cap, uint32_t, tb->n->bt_idx, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
				} else if(tip == 0){
					tip = 1;
					tip_idx = i;
				} else {
					return 0;
				}
			} else {
				if(tb->step > bts->buffer[tb->n->bt_idx].step){
					tb->n->bt_idx = i;
				}
			}
		}
		if(heap->size == 0 && tip == 1){
			return pop_bubble_backtrace_graph(g, bts, tip? tip_idx : heap->buffer[0]);
		} else if(heap->size == 1 && tip == 0){
			tb = ref_btv(bts, heap->buffer[0]);
			idx = bt->n->edges[bt->dir].idx;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				idx = f->next;
				e = g->edges->buffer + f->idx;
				if(e->closed) continue;
				m = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
				if(m->bt_visit == visit) continue;
				else tip ++;
			}
			if(tip == 0) return 0;
			return pop_bubble_backtrace_graph(g, bts, tip? tip_idx : heap->buffer[0]);
		}
	}
	return 0;
}

// very complicated local region, like yarn, but with Single In edge and Single Out edges
u8i resolve_yarns_graph(Graph *g, u4i max_step){
	btv *bts;
	u4v *heap;
	node_t *n;
	uint64_t nid, visit, ret, _ret;
	int dir;
	ret = 0;
	for(nid=0;nid<g->nodes->size;nid++) g->nodes->buffer[nid].bt_visit = 0;
	bts = init_btv(32);
	heap = init_u4v(32);
	visit = 0;
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		if(n->edges[0].cnt <= 1 && n->edges[1].cnt > 1){
			dir = 1;
		} else if(n->edges[1].cnt <= 1 && n->edges[0].cnt > 1){
			dir = 0;
		} else continue;
		_ret = resolve_yarn_core_graph(g, max_step, bts, heap, nid, dir, ++visit);
		if(_ret) ret ++;
	}
	free_btv(bts);
	free_u4v(heap);
	return ret;
}

u8i mask_all_branching_nodes_graph(Graph *g){
	node_t *n;
	u8i node, ret;
	ret = 0;
	for(node=0;node<g->nodes->size;node++){
		n = ref_nodev(g->nodes, node);
		if(n->closed) continue;
		if(n->edges[0].cnt > 1 || n->edges[1].cnt > 1){
			n->rep_idx = 1;
			ret ++;
		} else {
			n->rep_idx = 0;
		}
	}
	for(node=0;node<g->nodes->size;node++){
		n = ref_nodev(g->nodes, node);
		if(n->closed) continue;
		if(n->rep_idx == 0) continue;
		del_node_graph(g, n);
	}
	return ret;
}

uint64_t gen_unitigs_graph(Graph *g){
	tracev *path;
	u4v *lens;
	trace_t *t;
	node_t *n;
	uint64_t nid, nutg, i;
	for(i=0;i<g->utgs->size;i++) free_tracev(g->utgs->buffer[i]);
	clear_vplist(g->utgs);
	lens = init_u4v(1024);
	nutg = 0;
	for(nid=0;nid<g->nodes->size;nid++){
		g->nodes->buffer[nid].bt_visit = 0;
		g->nodes->buffer[nid].rep_idx  = MAX_REP_IDX;
	}
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		if(n->bt_visit) continue;
		path = init_tracev(4);
		nutg ++;
		t = next_ref_tracev(path);
		t->node = nid;
		t->edges[0] = EDGE_REF_NULL;
		t->edges[1] = EDGE_REF_NULL;
		t->dir = 0;
		true_linear_unique_trace_graph(g, path, 0xFFFFFFFFFFFFFFFFLLU, nutg, NULL);
		reverse_tracev(path);
		for(i=0;i<path->size;i++) path->buffer[i].dir = !path->buffer[i].dir;
		true_linear_unique_trace_graph(g, path, 0xFFFFFFFFFFFFFFFFLLU, nutg, NULL);
		push_u4v(lens, cal_offset_traces_graph(g, path, 0, path->size));
		for(i=0;i<path->size;i++){
			ref_nodev(g->nodes, path->buffer[i].node)->rep_idx = g->utgs->size;
		}
		push_vplist(g->utgs, path);
	}
	fprintf(hzm_debug_out, "[%s] ", date()); num_n50(lens, hzm_debug_out); fprintf(hzm_debug_out, "\n");
	free_u4v(lens);
	return nutg;
}

tracev* path2traces_graph(Graph *g, pathv *path){
	tracev *ts;
	trace_t *t1, *t2;
	path_t *ps[2];
	frg_t *frgs[2];
	edge_ref_t *f;
	lnk_t *l;
	u4i i, k, r;
	int z, x, y, d, found;
	ps[0] = ps[1] = NULL;
	frgs[0] = frgs[1] = NULL;
	// linking
	for(i=0;i<path->size;i++){
		ps[1] = ref_pathv(path, i);
		frgs[1] = ref_frgv(g->frgs, ps[1]->frg);
		ps[1]->tx = 0;
		ps[1]->ty = frgs[1]->tcnt - 1;
		if(ps[0]){
			f = ps[0]->lnks + ps[0]->dir;
			l = ref_lnkv(g->lnks, f->idx);
			k = f->flg ^ l->flag;
			r = !k;
			//TODO: F1(+) -> F2(+):10, F3(-) -> F2(-):9, may cause F2.tx > F2.ty
			if(ps[k]->dir ^ k){
				ps[k]->tx = l->tidx;
				t1 = ref_tracev(g->traces, frgs[k]->toff + ps[k]->tx);
			} else {
				ps[k]->ty = frgs[k]->ty + l->tidx;
				t1 = ref_tracev(g->traces, frgs[k]->toff + ps[k]->ty);
			}
			if(ps[r]->dir ^ r){
				x = 0; y = frgs[r]->ty; d = 1; y ++;
			} else {
				x = frgs[r]->tcnt - 1; y = frgs[r]->tx; d = -1; y --;
			}
			found = 0;
			for(z=x;z!=y;z+=d){
				t2 = ref_tracev(g->traces, frgs[r]->toff + z);
				if(t2->node == t1->node){
					found = 1;
					break;
				}
			}
			if(found == 0){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				return NULL;
			}
			//TODO: F1(+) -> F2(+):10, F3(-) -> F2(-):9, may cause F2.tx > F2.ty
			if(ps[r]->dir ^ r){
				ps[r]->tx = z;
			} else {
				ps[r]->ty = z;
			}
		}
		ps[0] = ps[1];
		frgs[0] = frgs[1];
	}
	ts = init_tracev(32);
	ps[0] = ps[1] = NULL;
	frgs[0] = frgs[1] = NULL;
	for(i=0;i<path->size;i++){
		ps[1] = ref_pathv(path, i);
		frgs[1] = ref_frgv(g->frgs, ps[1]->frg);
		x = ps[1]->tx; y = ps[1]->ty;
		if(ps[0]){
			t1 = ref_tracev(ts, ts->size - 1);
			if(ps[1]->dir){
				t2 = ref_tracev(g->traces, frgs[1]->toff + ps[1]->ty);
				t1->edges[t1->dir] = t2->edges[!t2->dir];
				y --;
			} else {
				t2 = ref_tracev(g->traces, frgs[1]->toff + ps[1]->tx);
				t1->edges[t1->dir] = t2->edges[t2->dir];
				x ++;
			}
		}
		if(ps[1]->dir){
			for(z=y;z>=x;z--){
				push_tracev(ts, g->traces->buffer[frgs[1]->toff + z]);
				ts->buffer[ts->size - 1].dir = !ts->buffer[ts->size - 1].dir;
			}
		} else {
			for(z=x;z<=y;z++) push_tracev(ts, g->traces->buffer[frgs[1]->toff + z]);
		}
		ps[0] = ps[1];
		frgs[0] = frgs[1];
	}
	return ts;
}

uint64_t gen_contigs_graph(Graph *g, FILE *out){
	pathv *path;
	tracev *ts;
	path_t *t;
	frg_t *n;
	uint64_t nid, nctg, i;
	for(i=0;i<g->ctgs->size;i++) free_tracev(g->ctgs->buffer[i]);
	clear_vplist(g->ctgs);
	nctg = 0;
	for(nid=0;nid<g->frgs->size;nid++) g->frgs->buffer[nid].bt_visit = 0;
	path = init_pathv(4);
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		if(n->closed) continue;
		if(n->bt_visit) continue;
		nctg ++;
		clear_pathv(path);
		t = next_ref_pathv(path);
		t->frg = nid;
		t->lnks[0] = EDGE_REF_NULL;
		t->lnks[1] = EDGE_REF_NULL;
		t->dir = 0;
		true_linear_unique_path_graph(g, path, 0xFFFFFFFFFFFFFFFFLLU, nctg, NULL);
		reverse_pathv(path);
		for(i=0;i<path->size;i++) path->buffer[i].dir = !path->buffer[i].dir;
		true_linear_unique_path_graph(g, path, 0xFFFFFFFFFFFFFFFFLLU, nctg, NULL);
		if((ts = path2traces_graph(g, path)) == NULL){
			continue;
		}
		cal_offset_paths_graph(g, path, 0, path->size);
		for(i=0;i<path->size;i++){
			t = ref_pathv(path, i);
			fprintf(out, "ctg%d\tF%d\t%c\t%d\n", (int)g->ctgs->size, t->frg, "+-*@"[t->dir], t->off);
		}
		push_vplist(g->ctgs, ts);
	}
	free_pathv(path);
	g->major_nctg = g->ctgs->size;
	return nctg;
}

u8i gen_complex_contigs_graph(Graph *g){
	tracev *ts;
	trace_t *t;
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	u8i i, idx, mi, cnt;
	u4i j, k, mk, mc;
	int mf;
	for(i=0;i<g->nodes->size;i++) g->nodes->buffer[i].unvisit = 1;
	cnt = 0;
	for(i=0;i<g->ctgs->size;i++){
		ts = (tracev*)get_vplist(g->ctgs, i);
		for(j=0;j<ts->size;j++){
			t = ref_tracev(ts, j);
			n = ref_nodev(g->nodes, t->node);
			n->unvisit = 0;
		}
	}
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->unvisit == 0) continue;
		if(n->regs.cnt < g->min_node_cov){
			n->unvisit = 0;
			continue;
		}
	}
	cnt = 0;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->unvisit == 0) continue;
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			mi = mk = mc = 0; mf = MAX_VALUE_B4;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				idx = f->next;
				e = g->edges->buffer + f->idx;
				if(e->node1 == e->node2) continue;
				if(e->cov < g->min_edge_cov) continue;
				if(e->off < mf){
					mi = f - g->erefs->buffer; mk = k; mc = e->cov; mf = e->off;
				}
			}
			if(mf == MAX_VALUE_B4) continue;
			f = ref_edgerefv(g->erefs, mi);
			if(f->flg) continue;
			e = ref_edgev(g->edges, f->idx);
			if(g->nodes->buffer[e->node2].unvisit == 0) continue;
			ts = init_tracev(2);
			t = next_ref_tracev(ts);
			t->node = i;
			t->dir = k;
			t->cov = 0;
			t->off = 0;
			t->edges[!k] = EDGE_REF_NULL;
			t->edges[k] = (edge_ref_t){f->idx, 0, 0};
			t = next_ref_tracev(ts);
			t->node = e->node2;
			t->dir  = e->dir2;
			t->cov = 0;
			t->off = 0;
			t->edges[e->dir2] = EDGE_REF_NULL;
			t->edges[!e->dir2] = (edge_ref_t){f->idx, 1, 0};
			push_vplist(g->ctgs, ts);
			cnt ++;
		}
	}
	return cnt;
}

void n50_stat_contigs_graph(Graph *g){
	tracev *ts;
	u4v *lens;
	int len;
	u8i i;
	lens = init_u4v(1024);
	for(i=0;i<g->major_nctg;i++){
		ts = (tracev*)get_vplist(g->ctgs, i);
		len = cal_offset_traces_graph(g, ts, 0, ts->size);
		push_u4v(lens, len);
	}
	fprintf(hzm_debug_out, "[%s] Estimated: ", date()); num_n50(lens, hzm_debug_out); fprintf(hzm_debug_out, "\n");
	free_u4v(lens);
}

// after gen_contigs
uint64_t print_isolated_nodes_dot_graph(Graph *g, FILE *out){
	tracev *ts;
	trace_t *t;
	node_t *n;
	reg_t *r, *rr;
	edge_ref_t *f;
	edge_t *e;
	u8i i, idx, cnt;
	u4i j, k, max;
	for(i=0;i<g->nodes->size;i++) g->nodes->buffer[i].unvisit = 1;
	cnt = 0;
	for(i=0;i<g->ctgs->size;i++){
		ts = (tracev*)get_vplist(g->ctgs, i);
		for(j=0;j<ts->size;j++){
			t = ref_tracev(ts, j);
			n = ref_nodev(g->nodes, t->node);
			n->unvisit = 0;
		}
	}
	fprintf(out, "digraph {\n");
	fprintf(out, "node [shape=record]\n");
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->unvisit == 0) continue;
		if(n->regs.cnt < g->min_node_cov){
			n->unvisit = 0;
			continue;
		}
		r = NULL; max = 0;
		for(j=0;j<n->regs.cnt;j++){
			rr = ref_regv(g->regs, n->regs.idx + j);
			if(g->reads->buffer[rr->rid].regs.cnt > max){
				r = rr;
				max = g->reads->buffer[rr->rid].regs.cnt;
			}
		}
		if(r == NULL) continue;
		fprintf(out, "N%llu [label=\"{N%llu %d | %s | %c_%d_%d}\"]\n", i, i, n->regs.cnt, g->dmo->reads->buffer[r->rid].tag, "FR"[r->dir], r->beg, r->end - r->beg);
		cnt ++;
	}
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->unvisit == 0) continue;
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				idx = f->next;
				e = g->edges->buffer + f->idx;
				if(e->cov < g->min_edge_cov) continue;
				if(f->flg){
					if(g->nodes->buffer[e->node1].unvisit == 0) continue;
					fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s%s]\n", i, (unsigned long long)e->node1, "+-"[k], "+-"[!e->dir1], e->cov, e->off, colors[k][!e->dir1], e->closed? " style=dashed" : "");
				} else {
					if(g->nodes->buffer[e->node2].unvisit == 0) continue;
					fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s%s]\n", i, (unsigned long long)e->node2, "+-"[k], "+-"[e->dir2], e->cov, e->off, colors[k][e->dir2], e->closed? " style=dashed" : "");
				}
			}
		}
	}
	fprintf(out, "}\n");
	return cnt;
}

// after gen_contigs
u4i count_isolated_reads_graph(Graph *g){
	tracev *ts;
	trace_t *t;
	node_t *n;
	read_t *rd;
	reg_t *r;
	u8i i, j, cnt, idx;
	int fnd;
	for(i=0;i<g->nodes->size;i++) g->nodes->buffer[i].unvisit = 1;
	cnt = 0;
	for(i=0;i<g->ctgs->size;i++){
		ts = (tracev*)get_vplist(g->ctgs, i);
		for(j=0;j<ts->size;j++){
			t = ref_tracev(ts, j);
			n = ref_nodev(g->nodes, t->node);
			n->unvisit = 0;
		}
	}
	for(i=0;i<g->reads->size;i++){
		rd = ref_readv(g->reads, i);
		fnd = 0;
		idx = rd->regs.idx;
		while(idx){
			r = ref_regv(g->regs, idx);
			idx = r->read_link;
			n = ref_nodev(g->nodes, r->node);
			if(n->unvisit == 0){ fnd = 1; break; }
		}
		if(fnd == 0) cnt ++;
	}
	return cnt;
}

void print_dot_subgraph(Graph *g, subnodehash *nodes, subedgev *edges, FILE *out){
	subnode_t *n1, *n2;
	subedge_t *e;
	u4i k, idx;
	fprintf(out, "digraph {\n");
	fprintf(out, " rankdir=LR\n");
	reset_iter_subnodehash(nodes);
	while((n1 = ref_iter_subnodehash(nodes))){
		if(n1->closed) continue;
		fprintf(out, "N%llu [label=\"N%llu(%llu) %d:%d:%d\" style=filled fillcolor=\"%s\" color=\"%s\"]\n", (u8i)n1->node, (u8i)n1->node, (u8i)g->nodes->buffer[n1->node].rep_idx, n1->cov, n1->visit, n1->bt_open, n1->fix? "yellow" : "white", n1->visit? "green" : (n1->cov > 2? "blue" : "black"));
	}
	reset_iter_subnodehash(nodes);
	while((n1 = ref_iter_subnodehash(nodes))){
		if(n1->closed) continue;
		for(k=0;k<2;k++){
			idx = n1->edges[k].idx;
			while(idx){
				e = ref_subedgev(edges, idx);
				idx = e->next;
				if(e->fwd == 0) continue;
				if(e->closed) continue;
				n2 = e->node;
				fprintf(out, " N%llu -> N%llu [label=\"%c%c:%d:%d\" color=\"%s\" %s]\n", (u8i)n1->node, (u8i)n2->node, "+-"[k], "+-"[e->dir], e->cov, e->visit, e->cov > 1? "blue" : "black", e->visit? "style=dashed":"");
			}
		}
	}
	fprintf(out, "}\n");
	fflush(out);
}

void fprintf_dot_subgraph(Graph *g, subnodehash *nodes, subedgev *edges, char *name_prefix, char *name_suffix){
	FILE *out;
	out = open_file_for_write(name_prefix, name_suffix, 1);
	print_dot_subgraph(g, nodes, edges, out);
	fclose(out);
}

typedef struct {
	u4i node:31, dir:1;
	u4i flag;
	u4i prev;
	u4i step;
	int score;
} sg_heap_t;
define_list(sgheapv, sg_heap_t);

typedef struct {
	u4i node:31, dir:1;
	u4i group:30, solid:1, closed:1;
} sg_tip_t;
define_list(sgtipv, sg_tip_t);

subedge_t* find_edge_subgraph(subnodehash *nodes, subedgev *edges, u4i node1, int dir1, u4i node2, int dir2){
	subnode_t *n;
	subedge_t *e;
	u8i idx;
	n = ref_subnodehash(nodes, node1);
	idx = n->edges[dir1].idx;
	e = ref_subedgev(edges, idx);
	if(offset_subnodehash(nodes, e->node) == node2 && e->dir == dir2){
		return e;
	}
	while((idx = e->next)){
		e = ref_subedgev(edges, idx);
		if(offset_subnodehash(nodes, e->node) == node2 && e->dir == dir2){
			return e;
		}
	}
	return NULL;
}

int cut_edge_core_subgraph(subnodehash *nodes, subedgev *edges, u4i node1, int dir1, u4i node2, int dir2){
	subnode_t *n;
	subedge_t *e, *p;
	u8i idx;
	n = ref_subnodehash(nodes, node1);
	idx = n->edges[dir1].idx;
	e = ref_subedgev(edges, idx);
	if(offset_subnodehash(nodes, e->node) == node2 && e->dir == dir2){
		e->closed = 1;
		n->edges[dir1].idx = e->next;
		n->edges[dir1].cnt --;
		return 1;
	}
	while((idx = e->next)){
		p = e;
		e = ref_subedgev(edges, idx);
		if(offset_subnodehash(nodes, e->node) == node2 && e->dir == dir2){
			e->closed = 1;
			p->next = e->next;
			n->edges[dir1].cnt --;
			return 1;
		}
	}
	return 0;
}

int cut_edge_subgraph(subnodehash *nodes, subedgev *edges, u4i node1, int dir1, u4i node2, int dir2){
	return cut_edge_core_subgraph(nodes, edges, node1, dir1, node2, dir2)
		+ cut_edge_core_subgraph(nodes, edges, node2, !dir2, node1, !dir1);
}

u4i cut_loopback_edges_subgraph(tracev *traces, u4i tb, subnodehash *nodes, subedgev *edges, sgheapv *heap){
	subnode_t *n, *n1, *n2, N;
	subedge_t *e, *e2;
	sg_heap_t *p, P;
	trace_t *t1;
	u8i idx;
	u4i k, i, ret, open, round;
	ret = 0;
	round = 0;
	while(1){
		round ++;
		clear_sgheapv(heap);
		reset_iter_subnodehash(nodes);
		while((n = ref_iter_subnodehash(nodes))){
			n->flag = 0xFFFFFFFFFFFFFFFLLU;
			n->visit = 0;
			n->bt_nidx = MAX_BT_NIDX;
			n->bt_dir = 0;
			n->bt_open = 0;
			n->bt_step = 0;
			n->bt_score = 0;
		}
		for(i=0;i<edges->size;i++) edges->buffer[i].visit = 0;
		t1 = ref_tracev(traces, tb);
		memset(&N, 0, sizeof(subnode_t));
		N.node = t1->node;
		n1 = get_subnodehash(nodes, N);
		array_heap_push(heap->buffer, heap->size, heap->cap, sg_heap_t, ((sg_heap_t){offset_subnodehash(nodes, n1), (n1->edges[0].cnt == 0), 0, 0xFFFFFFFFU, 0, 0}), num_cmp(b.score, a.score));
		p = &P;
		open = 0;
		while(heap->size){
			P = heap->buffer[0];
			array_heap_remove(heap->buffer, heap->size, heap->cap, sg_heap_t, 0, num_cmp(b.score, a.score));
			n = ref_subnodehash(nodes, p->node);
			k = p->dir;
			if(n->visit == 0){
				n->visit = 1;
				n->bt_open = n->edges[!k].cnt;
				n->bt_dir  = k;
				n->bt_nidx = p->prev;
				n->bt_step = p->step;
				n->bt_score = p->score;
				open ++;
			} else {
				if(k != n->bt_dir){
					// found circle
					//cut_edge_subgraph(nodes, edges, p->prev, nodes->array[p->prev].bt_dir, p->node, p->dir);
					//ret ++;
				}
			}
			if(n->bt_open) n->bt_open --;
			if(n->bt_open == 0){
				open --;
				idx = n->edges[k].idx;
				while(idx){
					e = ref_subedgev(edges, idx);
					idx = e->next;
					if(e->closed) continue;
					if(e->visit){
						static int warn_cnt = 0;
						warn_cnt ++;
						if(warn_cnt == 10){
							fprintf(stderr, " -- too many warnings in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						} else if(warn_cnt < 10){
							fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						} else {
							warn_cnt = 11;
						}
					}
					e->visit = 1; // will be visited
					array_heap_push(heap->buffer, heap->size, heap->cap, sg_heap_t, ((sg_heap_t){offset_subnodehash(nodes, e->node), e->dir, n->flag, p->node, p->step + 1, e->cov + n->bt_score}), num_cmp(b.score, a.score));
				}
			}
		}
		if(open == 0) break;
		//try to find the latest opened node
		n2 = NULL;
		reset_iter_subnodehash(nodes);
		while((n = ref_iter_subnodehash(nodes))){
			if(n->visit == 0) continue;
			if(n->bt_open == 0) continue;
			if(n2 == NULL) n2 = n;
			else if(n->bt_step < n2->bt_step){
				n2 = n;
			}
		}
		if(n2){
			n = n2;
			idx = n->edges[!n->bt_dir].idx;
			while(idx){
				e = ref_subedgev(edges, idx);
				idx = e->next;
				if(e->closed) continue;
				e2 = find_edge_subgraph(nodes, edges, offset_subnodehash(nodes, e->node), !e->dir, offset_subnodehash(nodes, n), n->bt_dir);
				if(e2->visit) continue;
				cut_edge_subgraph(nodes, edges, offset_subnodehash(nodes, n), !n->bt_dir, offset_subnodehash(nodes, e->node), e->dir);
				ret ++;
			}
		}
	}
	return ret;
}

//#define LOCAL_ASSEMBLY_NO_TRIMMING

int local_assembly_core_subgraph(Graph *g, tracev *traces, u4i tstar, int tdir, subnodehash *nodes, subedgev *edges, sgtipv *tips, sgheapv *heap, tracev *exts){
	subnode_t *n, *n1, *n2, N;
	subedge_t *e;
	sg_heap_t *p, P;
	edge_ref_t *f;
	trace_t *t1, *t2;
	u8i idx;
	u4i k, i, ntip, flg, ncov, ntot, bt_hit;
	int debug = 0;
	clear_sgtipv(tips);
	clear_tracev(exts);
	// scan tips
	clear_and_encap_sgheapv(heap, nodes->count);
	ntip = 0;
	t1 = ref_tracev(traces, tstar);
	reset_iter_subnodehash(nodes);
	while((n = ref_iter_subnodehash(nodes))){
		n->flag = 0xFFFFFFFFFFFFFFFLLU;
		n->visit = 0;
		n->bt_nidx = MAX_BT_NIDX;
		n->bt_dir = 0;
		n->bt_open = 0;
		n->bt_step = 0;
		n->bt_score = 0;
		n->bt_hit = (n->node == t1->node);
		if(n->fix) continue;
		if(n->edges[0].cnt && n->edges[1].cnt) continue;
		array_heap_push(heap->buffer, heap->size, heap->cap, sg_heap_t, ((sg_heap_t){offset_subnodehash(nodes, n), (n->edges[0].cnt == 0), ntip, 0xFFFFFFFFU, 0, 0}), num_cmpx(b.score, a.score, nodes->array[a.node].node, nodes->array[b.node].node)); //replace step into node_id
		push_sgtipv(tips, (sg_tip_t){offset_subnodehash(nodes, n), (n->edges[0].cnt == 0), ntip, 0, 0});
		++ ntip;
		if(debug) fprintf(stderr, "tip N%llu\n", (u8i)n->node);
	}
	if(ntip == 0) return 0; // MUST have circle
	// try to cut tips
	p = &P;
	bt_hit = 0;
	while(heap->size){
		P = heap->buffer[0];
		array_heap_remove(heap->buffer, heap->size, heap->cap, sg_heap_t, 0, num_cmpx(b.score, a.score, a.step, b.step));
		n = ref_subnodehash(nodes, p->node);
		flg = tips->buffer[p->flag].group;
		k = p->dir;
		if(debug) fprintf(stderr, "pop: N%llu tip=%u dir=%d visit=%d\n", (u8i)n->node, (int)flg, k, n->visit);
		if(n->visit == 0){
			n->visit = 1;
			n->bt_open = n->edges[!k].cnt;
			n->bt_dir  = k;
			n->bt_nidx = p->prev;
			n->bt_step = p->step;
			n->bt_score = p->score;
			n->flag    = flg;
			n->bt_hit  = bt_hit;
		} else {
			n->flag = tips->buffer[n->flag].group;
			if(k != n->bt_dir){
				if(debug) fprintf(stderr, "circle: N%llu\n", (u8i)n->node);
				// found circle
				return 0;
			} else if(n->flag != flg){
				if(debug) fprintf(stderr, "met: N%llu tips=%d,%d\n", (u8i)n->node, (int)flg, (int)n->flag);
				// found common ancester
#ifdef LOCAL_ASSEMBLY_NO_TRIMMING
				if(0){
#else
				if((tips->buffer[flg].solid)){
					if(tips->buffer[n->flag].solid){
						// cannot cut tip
						if(debug) fprintf(stderr, "crash: N%llu flag:%d,%d\n", (u8i)n->node, (int)flg, (int)n->flag);
						return 0;
					} else {
						// cut n->flag
						tips->buffer[n->flag].closed = 1;
						tips->buffer[n->flag].group  = flg;
						if(debug) fprintf(stderr, "cut: tip%d -> tip%d\n", (int)n->flag, flg);
						ntip --;
						n->flag = flg;
						n->bt_dir  = k;
						n->bt_nidx = p->prev;
					}
				} else if(tips->buffer[n->flag].solid){
					// cut flg
					tips->buffer[flg].closed = 1;
					tips->buffer[flg].group  = n->flag;
					if(debug) fprintf(stderr, "cut: tip%d -> tip%d\n", (int)flg, (int)n->flag);
					ntip --;
#endif
				} else if(p->score > n->bt_score){
					// cut n->flag
					tips->buffer[n->flag].closed = 1;
					tips->buffer[n->flag].group  = flg;
					if(debug) fprintf(stderr, "cut: tip%d -> tip%d\n", (int)n->flag, (int)n->flag);
					ntip --;
					n->flag = flg;
					n->bt_dir  = k;
					n->bt_nidx = p->prev;
					n->bt_score = p->score;
				} else {
					// cut flg
					tips->buffer[flg].closed = 1;
					tips->buffer[flg].group  = n->flag;
					if(debug) fprintf(stderr, "cut: tip%d -> tip%d\n", (int)flg, (int)n->flag);
					ntip --;
				}
			} else if(p->score > n->bt_score){
				n->bt_dir  = k;
				n->bt_nidx = p->prev;
				n->bt_score = p->score;
			}
		}
		if(n->bt_open) n->bt_open --;
		if(n->bt_open == 0){
			if(n->node == t1->node) bt_hit = 1;
			if(!n->fix && n->cov >= g->max_node_cov_sg) tips->buffer[n->flag].solid = 1;
			idx = n->edges[k].idx;
			while(idx){
				e = ref_subedgev(edges, idx);
				idx = e->next;
				if(e->closed) continue;
				//array_heap_push(heap->buffer, heap->size, heap->cap, sg_heap_t, ((sg_heap_t){offset_subnodehash(nodes, e->node), e->dir, n->flag, p->node, p->step + 1, e->node->cov + n->bt_score}), num_cmpx(b.score, a.score, a.step, b.step));
				array_heap_push(heap->buffer, heap->size, heap->cap, sg_heap_t, ((sg_heap_t){offset_subnodehash(nodes, e->node), e->dir, n->flag, p->node, p->step + 1, e->node->cov + e->cov + n->bt_score}), num_cmpx(b.score, a.score, a.step, b.step));
				if(debug) fprintf(stderr, "push: N%llu tip=%d dir=%d\n", (u8i)e->node->node, (int)n->flag, e->dir);
			}
		}
	}
#ifdef LOCAL_ASSEMBLY_NO_TRIMMING
#else
	if(ntip > 1){
		u4i fid;
		fid = MAX_VALUE_U4;
		ntip = 0;
		for(i=0;i<tips->size;i++){
			if(tips->buffer[i].closed) continue;
			if(tips->buffer[i].node == fid) continue;
			ntip ++;
		}
		if(ntip > 1){
			return 0;
		}
	}
#endif
	// whether the backtrace path has solid
	for(i=0;i<tips->size;i++){
		if(tips->buffer[i].closed) continue;
		if(tips->buffer[i].solid == 0){
			//fprintf(stderr, " -- N%llu in %s -- %s:%d --\n", (u8i)traces->buffer[tstar].node, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			return 0;
		}
	}
	// traceback from end trace
	t1 = ref_tracev(traces, tstar);
	memset(&N, 0, sizeof(subnode_t));
	N.node = t1->node;
	n1 = get_subnodehash(nodes, N);
	if(n1 == NULL){
		fprintf(stderr, " -- Something worng in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 0;
	}
	if(n1->visit == 0){
		return 0;
	}
	t1->dir = t1->dir ^ tdir;
	ntot = 0;
	u8i totn = 0;
	reset_iter_subnodehash(nodes);
	while((n = ref_iter_subnodehash(nodes))){
		if(n->fix) continue;
		if(n->bt_hit == 1) continue;
		ntot += n->cov;
		if(n->bt_hit == 1) continue;
		totn += n->cov;
	}
	ncov = 0;
	while(n1->bt_nidx != MAX_BT_NIDX){
		n2 = ref_subnodehash(nodes, n1->bt_nidx);
		ncov += n2->cov;
		f = edge_node2node_graph(g, n1->node, !n1->bt_dir, n2->node, !n2->bt_dir);
		if(f == NULL){
			// (store_low_cov is 0 and edge cov is less than g->min_edge_cov)
			break;
		}
		t1->edges[t1->dir] = (edge_ref_t){f->idx, f->flg, 0};
		t2 = next_ref_tracev(exts);
		t2->node = n2->node;
		t2->dir  = !n2->bt_dir;
		t2->cov  = n2->cov;
		t2->off  = 0;
		t2->edges[t2->dir] = EDGE_REF_NULL;
		t2->edges[(!t2->dir)] = (edge_ref_t){f->idx, !f->flg, 0};
		n1 = n2;
		t1 = t2;
	}
	t1 = ref_tracev(traces, tstar);
	t1->dir = t1->dir ^ tdir;
	if(tdir){
		reverse_tracev(exts);
		for(i=0;i<exts->size;i++) exts->buffer[i].dir = !exts->buffer[i].dir;
	}
	//if((ncov < ntot / 2) ^ (ncov < totn / 2)){
		//fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	//}
	N.node = t1->node;
	n1 = get_subnodehash(nodes, N);
	if(ncov < ntot / 2){
		return 0;
	}
	return 1;
}

typedef struct {
	u4i rid:31, rk:1;
	u8i beg_idx, end_idx;
} sub_rd_t;
define_list(subrdv, sub_rd_t);
#define subrd_hashcode(E) u32hashcode((E).rid)
#define subrd_hashequals(E1, E2) (E1).rid == (E2).rid
define_hashset(subrdhash, sub_rd_t, subrd_hashcode, subrd_hashequals);

u4i extending_unitig_core_graph(Graph *g, u4i tid, u4i max_end, u4i trgs[2], tracev *exts[2], subnodehash *nodes, subedgev *edges, subrdhash *hash, subrdv *srs, sgtipv *tips, sgheapv *heap){
	tracev *ts;
	trace_t *t;
	node_t *_n;
	reg_t *reg, *r;
	subnode_t N, *n, *n1, *n2;
	subedge_t *e;
	sub_rd_t *sr, SR;
	u8i idx, edx;
	u4i tbeg, tend, tk, k1, k2, rk, i, sri;
	int tb, te, ti, exists;
	ts = (tracev*)get_vplist(g->utgs, tid);
	memset(&N, 0, sizeof(subnode_t));
	memset(&SR, 0, sizeof(sub_rd_t));
	N.cov = 1;
	tbeg = 0;
	tend = ts->size;
	for(tk=0;tk<2;tk++){
		clear_tracev(exts[tk]);
		while(tbeg < tend){
			if(tk){
				te = tbeg - 1; tb = num_min(max_end + tbeg, tend) - 1;
			} else {
				te = tend; tb = (tbeg + max_end < tend)? tend - max_end : tbeg;
			}
			clear_subnodehash(nodes);
			clear_subedgev(edges);
			memset(next_ref_subedgev(edges), 0, sizeof(subedge_t));
			clear_subrdhash(hash);
			//collect reads
			for(ti=tb;ti!=te;tk? ti-- : ti++){
				t = ref_tracev(ts, ti);
				_n = ref_nodev(g->nodes, t->node);
				for(i=0;i<_n->regs.cnt;i++){
					reg = ref_regv(g->regs, _n->regs.idx + i);
					SR.rid = reg->rid;
					SR.rk = tk ^ t->dir ^ reg->dir;
					if(SR.rk){
						SR.beg_idx = g->reads->buffer[reg->rid].regs.idx;
						SR.end_idx = reg->read_link;
					} else {
						SR.beg_idx = _n->regs.idx + i;
						SR.end_idx = 0;
					}
					sr = prepare_subrdhash(hash, SR, &exists);
					if(exists) continue;
					*sr = SR;
				}
			}
			clear_subrdv(srs);
			reset_iter_subrdhash(hash);
			while((sr = ref_iter_subrdhash(hash))){
				push_subrdv(srs, *sr);
			}
			sort_array(srs->buffer, srs->size, sub_rd_t, num_cmpgtx(a.rid, b.rid, a.rk, b.rk));
			//build nodes
			for(sri=0;sri<srs->size;sri++){
				sr = ref_subrdv(srs, sri);
				rk = sr->rk;
				idx = sr->beg_idx;
				while(idx != sr->end_idx){
					r = ref_regv(g->regs, idx);
					idx = r->read_link;
					N.node = r->node;
					if(1){
						// try to ignore the nodes not in unitigs
						if(g->nodes->buffer[r->node].rep_idx == MAX_REP_IDX) continue;
					}
					n = prepare_subnodehash(nodes, N, &exists);
					if(exists){
						n->cov ++;
					} else {
						*n = N;
					}
					n->edges[rk ^ r->dir].cnt ++;
				}
			}
			// mask fixed nodes
			for(ti=tb;ti!=te;tk? ti-- : ti++){
				t = ref_tracev(ts, ti);
				N.node = t->node;
				n = get_subnodehash(nodes, N);
				if(n == NULL) continue;
				n->fix = 1;
			}
			// remove low occ strand of nodes
			reset_iter_subnodehash(nodes);
			while((n = ref_iter_subnodehash(nodes))){
				n->sub_dir = n->edges[0].cnt < n->edges[1].cnt;
				n->cov     = n->edges[n->sub_dir].cnt;
				n->edges[0].cnt = 0;
				n->edges[1].cnt = 0;
			}
			//build edges
			for(sri=0;sri<srs->size;sri++){
				sr = ref_subrdv(srs, sri);
				rk = sr->rk;
				idx = sr->beg_idx;
				n1 = n2 = NULL;
				k1 = k2 = 0;
				while(idx != sr->end_idx){
					r = ref_regv(g->regs, idx);
					idx = r->read_link;
					if(1){
						// try to ignore the nodes not in unitigs
						if(g->nodes->buffer[r->node].rep_idx == MAX_REP_IDX) continue;
					}
					N.node = r->node;
					k2 = r->dir;
					n2 = get_subnodehash(nodes, N);
					if((k2 ^ rk) != n2->sub_dir){
						continue;
					}
					if(n1){
						// link n1 to n2
						{
							edx = n1->edges[k1].idx;
							while(edx){
								e = ref_subedgev(edges, edx);
								if(e->node == n2 && e->dir == k2){
									e->cov ++; break;
								}
								edx = e->next;
							}
							if(edx == 0){
								edx = edges->size;
								e = next_ref_subedgev(edges);
								e->node = n2;
								e->fwd = !rk;
								e->dir = k2;
								e->cov = 1;
								e->closed = 0;
								e->visit = 0;
								e->next = n1->edges[k1].idx;
								n1->edges[k1].idx = edx;
								n1->edges[k1].cnt ++;
							}
						}
						// link n2 to n1
						{
							edx = n2->edges[!k2].idx;
							while(edx){
								e = ref_subedgev(edges, edx);
								if(e->node == n1 && e->dir == !k1){
									e->cov ++;
									break;
								}
								edx = e->next;
							}
							if(edx == 0){
								edx = edges->size;
								e = next_ref_subedgev(edges);
								e->node = n1;
								e->fwd = rk;
								e->dir = !k1;
								e->cov = 1;
								e->closed = 0;
								e->visit = 0;
								e->next = n2->edges[!k2].idx;
								n2->edges[!k2].idx = edx;
								n2->edges[!k2].cnt ++;
							}
						}
					}
					n1 = n2;
					k1 = k2;
				}
			}
			// try local assembly
			// TODO: if nodes used in this round are trimmed in next round, how about this assembly
			// Now, assuming the result of local assembly is right even trimming some nodes
			//if(ts->buffer[tk? tbeg : tend - 1].node == 9188){
				//fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			//}
			//if(ts->size > 100 && tk == 0 && tend == ts->size){ // never happen
				//fprintf_dot_subgraph(g, nodes, edges, "debug", ".dot");
			//}
			if(cut_loopback_edges_subgraph(ts, tb, nodes, edges, heap)){
				if(0) return 0;
			}
			if(local_assembly_core_subgraph(g, ts, tk? tbeg : tend - 1, tk, nodes, edges, tips, heap, exts[tk])){
				break;
			} else {
				//if(ref_nodev(g->nodes, ts->buffer[tk? tbeg : tend - 1].node)->init_end){
					//// Don't trim dead ends just after grapg construction
					//break;
				//}
				if(tk) tbeg ++;
				else tend --;
			}
		}
	}
	trgs[0] = tbeg;
	trgs[1] = tend;
	if(tbeg < tend) return tend - tbeg;
	else return 0;
}

thread_beg_def(uext);
Graph *g;
u4i max_end;
u4i tid;
tracev *exts[2];
u4i t_regs[2];
thread_end_def(uext);

thread_beg_func(uext);
subnodehash *nodes;
subedgev *edges;
sgtipv *tips;
sgheapv *heap;
subrdhash *hash;
subrdv *srs;
nodes = init_subnodehash(1023);
edges = init_subedgev(32);
tips = init_sgtipv(32);
heap = init_sgheapv(1024);
hash = init_subrdhash(13);
srs  = init_subrdv(32);
thread_beg_loop(uext);
if(uext->tid == 0xFFFFFFFFU) continue;
extending_unitig_core_graph(uext->g, uext->tid, uext->max_end, (u4i*)uext->t_regs, (tracev**)uext->exts, nodes, edges, hash, srs, tips, heap);
thread_end_loop(uext);
free_subnodehash(nodes);
free_subedgev(edges);
free_sgtipv(tips);
free_sgheapv(heap);
free_subrdhash(hash);
free_subrdv(srs);
thread_end_func(uext);

// MUST call gen_unitigs_graph first
u4i unitigs2frgs_graph(Graph *g, int ncpu){
	u4v *lens;
	frg_t *frg;
	node_t *n;
	tracev *ts;
	trace_t *t;
	u4i i, j, tid, ret;
	thread_preprocess(uext);
	thread_beg_init(uext, ncpu);
	uext->g = g;
	uext->max_end = g->max_sg_end;
	uext->tid = 0xFFFFFFFFU;
	uext->exts[0] = init_tracev(32);
	uext->exts[1] = init_tracev(32);
	uext->t_regs[0] = 0;
	uext->t_regs[1] = 0;
	thread_end_init(uext);
	{
		clear_frgv(g->frgs);
		clear_lnkv(g->lnks);
		clear_edgerefv(g->lrefs);
		clear_tracev(g->traces);
	}
	lens = init_u4v(32);
	ret = 0;
	for(tid=0;tid<=g->utgs->size+(u4i)ncpu;tid++){
		if(tid < g->utgs->size){
			if((tid % 1000) == 0 && !hzm_debug){ fprintf(hzm_debug_out, "\r%u", tid); fflush(hzm_debug_out); }
			thread_wait_one(uext);
		} else {
			thread_wait_next(uext);
		}
		if(uext->tid != 0xFFFFFFFFU){
			if(uext->t_regs[1] > uext->t_regs[0]){
				ret ++;
				frg = next_ref_frgv(g->frgs);
				frg->toff = g->traces->size;
				frg->lnks[0] = PTR_REF_NULL;
				frg->lnks[1] = PTR_REF_NULL;
				frg->closed = 0;
				ts = (tracev*)get_vplist(g->utgs, uext->tid);
				append_array_tracev(g->traces, uext->exts[1]->buffer, uext->exts[1]->size);
				frg->tx = uext->exts[1]->size;
				frg->ty = uext->exts[1]->size + uext->t_regs[1] - uext->t_regs[0];
				for(i=uext->t_regs[0];i<uext->t_regs[1];i++) ts->buffer[i].cov = 0;
				append_array_tracev(g->traces, ts->buffer + uext->t_regs[0], uext->t_regs[1] - uext->t_regs[0]);
				append_array_tracev(g->traces, uext->exts[0]->buffer, uext->exts[0]->size);
				frg->tcnt = g->traces->size - frg->toff;
				frg->len    = cal_offset_traces_graph(g, g->traces, frg->toff + frg->tx, frg->toff + frg->ty);
				frg->length = cal_offset_traces_graph(g, g->traces, frg->toff, frg->toff + frg->tcnt);
				push_u4v(lens, frg->length);
			}
		}
		if(tid < g->utgs->size){
			uext->tid = tid;
			thread_wake(uext);
		} else {
			uext->tid = 0xFFFFFFFFU;
		}
	}
	fprintf(hzm_debug_out, "\r%u untigs -> %u\n", (u4i)g->utgs->size, (u4i)g->frgs->size);
	fprintf(hzm_debug_out, "[%s] ", date()); num_n50(lens, hzm_debug_out); fprintf(hzm_debug_out, "\n"); fflush(hzm_debug_out);
	free_u4v(lens);
	thread_beg_close(uext);
	free_tracev(uext->exts[0]);
	free_tracev(uext->exts[1]);
	thread_end_close(uext);
	psort_array(g->frgs->buffer, g->frgs->size, frg_t, ncpu, num_cmpgt(b.length, a.length));
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		n->unvisit = 1;
		n->rep_idx = MAX_REP_IDX;
		n->bt_visit = 0;
	}
	for(i=0;i<g->frgs->size;i++){
		frg = ref_frgv(g->frgs, i);
		for(j=frg->tx;j<frg->ty;j++){
			t = ref_tracev(g->traces, frg->toff + j);
			n = ref_nodev(g->nodes, t->node);
			n->rep_idx = i;
			n->rep_dir = t->dir;
			n->bt_visit = j;
		}
	}
	return ret;
}

typedef struct {
	u8i node;
	u4i frgs[2];
	b4i offs[2];
	u1i dirs[2][2];
	u2i tidx[2];
	u2i nfrg, cov;
} weak_frg_node_t;
#define weakfn_hashcode(E) u64hashcode((E).node)
#define weakfn_hashequals(E1, E2) ((E1).node == (E2).node)
define_hashset(weakfnhash, weak_frg_node_t, weakfn_hashcode, weakfn_hashequals);

int _node_off2dist_sg(Graph *g, u4i f_idx, u4i t_idx, int fdir){
	frg_t *frg;
	trace_t *t0, *t1;
	node_t *n0, *n1;
	reg_t *rg0, *rg1;
	frg = ref_frgv(g->frgs, f_idx);
	t0 = ref_tracev(g->traces, frg->toff + (fdir? frg->tx : frg->ty - 1));
	n0 = ref_nodev(g->nodes, t0->node);
	rg0 = ref_regv(g->regs, n0->regs.idx);
	t1 = ref_tracev(g->traces, frg->toff + t_idx);
	n1 = ref_nodev(g->nodes, t1->node);
	rg1 = ref_regv(g->regs, n1->regs.idx);
	return fdir? t0->off - (t1->off + rg1->end - rg1->beg) : t1->off - (t0->off + rg0->end - rg0->beg);
}

u4i gen_lnks_graph(Graph *g, int ncpu){
	frg_t *frg;
	trace_t *t, *tt;
	node_t *n, *n1, *n2;
	reg_t *rg1, *rg2;
	lnkhash *hash;
	lnk_t *l, L;
	weakfnhash *weaks;
	weak_frg_node_t *w, W;
	edge_ref_t F1, F2;
	u8i lst, idx;
	u4i i, cnt, fdir;
	cnt = 0;
	u4i wushigang = cnt;
	u4i tmp = wushigang;
	wushigang = tmp;
	int exists, j, x, y, inc;
	clear_lnkv(g->lnks);
	memset(next_ref_lnkv(g->lnks), 0, sizeof(lnk_t));
	clear_edgerefv(g->lrefs);
	push_edgerefv(g->lrefs, EDGE_REF_NULL);
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		n->unvisit = 1;
		n->rep_idx = MAX_REP_IDX;
		n->bt_visit = 0;
	}
	for(i=0;i<g->frgs->size;i++){
		frg = ref_frgv(g->frgs, i);
		for(j=frg->tx;j<(int)frg->ty;j++){
			t = ref_tracev(g->traces, frg->toff + j);
			n = ref_nodev(g->nodes, t->node);
			n->rep_idx = i;
			n->rep_dir = t->dir;
			n->bt_visit = j;
		}
	}
	hash = init_lnkhash(1023);
	memset(&L, 0, sizeof(lnk_t));
	for(i=0;i<g->frgs->size;i++){
		frg = ref_frgv(g->frgs, i);
		for(fdir=0;fdir<2;fdir++){
			x = fdir? 0 : frg->ty;
			y = fdir? frg->tx : frg->tcnt;
			x = fdir? frg->tx : frg->ty;
			y = fdir? 0 : frg->tcnt;
			inc = fdir? -1 : 1;
			for(j=x;j!=y;j+=inc){
				t = ref_tracev(g->traces, frg->toff + j);
				//TODO: whether to use low coverage read path
				//if(t->cov < g->max_node_cov_sg) continue;
				n = ref_nodev(g->nodes, t->node);
				if(n->rep_idx == MAX_REP_IDX) continue;
				if(n->rep_idx == i) continue;
				{
					if(i < n->rep_idx){
						L.frg1 = i;
						L.frg2 = n->rep_idx;
						L.dir1 = fdir ^ 0;
						L.dir2 = fdir ^ (n->rep_dir ^ t->dir);
						L.off = _node_off2dist_sg(g, i, j, L.dir1) + _node_off2dist_sg(g, n->rep_idx, n->bt_visit, !L.dir2);
						L.flag = 0;
					} else {
						L.frg1 = n->rep_idx;
						L.frg2 = i;
						L.dir1 = fdir ^ (!n->rep_dir ^ t->dir);
						L.dir2 = fdir ^ 1;
						L.off = _node_off2dist_sg(g, i, j, !L.dir2) + _node_off2dist_sg(g, n->rep_idx, n->bt_visit, L.dir1);
						L.flag = 1;
					}
					if(L.off + (int)g->reglen < 0) continue;
					L.tidx = j - (fdir? 0 : x);
					L.cov = t->cov;
				}
				l = prepare_lnkhash(hash, L, &exists);
				if(exists){
					if(l->cov < L.cov){
						*l = L;
					}
				} else {
					*l = L;
				}
			}
		}
	}
	memset(&W, 0, sizeof(weak_frg_node_t));
	W.cov = 0xFFFFU;
	weaks = init_weakfnhash(1023);
	for(i=0;i<g->frgs->size;i++){
		frg = ref_frgv(g->frgs, i);
		for(fdir=0;fdir<2;fdir++){
			x = fdir? 0 : frg->ty;
			y = fdir? frg->tx : frg->tcnt;
			tt = ref_tracev(g->traces, frg->toff + (fdir? frg->tx : frg->ty - 1));
			n1  = ref_nodev(g->nodes, tt->node);
			rg1 = ref_regv(g->regs, n1->regs.idx);
			for(j=x;j<y;j++){
				t = ref_tracev(g->traces, frg->toff + j);
				n2 = ref_nodev(g->nodes, t->node);
				if(n2->rep_idx != MAX_REP_IDX) continue;
				if(t->cov == 1) continue; //TODO: singleton
				rg2 = ref_regv(g->regs, n2->regs.idx);
				W.node = t->node;
				w = prepare_weakfnhash(weaks, W, &exists);
				if(exists){
				} else {
					*w = W;
				}
				if(w->nfrg < 2){
					w->frgs[w->nfrg] = i;
					w->offs[w->nfrg] = fdir? tt->off - (t->off + rg2->end - rg2->beg) : t->off - (tt->off + rg1->end - rg1->beg);
					w->dirs[w->nfrg][0] = fdir;
					w->dirs[w->nfrg][1] = fdir ^ t->dir;
					w->tidx[w->nfrg] = j - x;
					if(w->cov > t->cov) w->cov = t->cov;
				}
				w->nfrg ++;
			}
		}
	}
	reset_iter_weakfnhash(weaks);
	L.weak = 1;
	L.closed = WT_EDGE_CLOSED_LESS;
	while((w = ref_iter_weakfnhash(weaks))){
		if(w->nfrg != 2) continue;
		if(w->frgs[0] == w->frgs[1]) continue;
		if((w->dirs[0][1]) == (w->dirs[1][1])) continue;
		if(w->cov < g->max_node_cov_sg) continue; // TODO
		rg1 = ref_regv(g->regs, ref_nodev(g->nodes, w->node)->regs.idx);
		L.cov = w->cov;
		L.off = w->offs[0] + w->offs[1] + rg1->end - rg1->beg;
		if(w->frgs[0] < w->frgs[1]){
			L.frg1 = w->frgs[0];
			L.frg2 = w->frgs[1];
			L.dir1 = w->dirs[0][0];
			L.dir2 = !w->dirs[1][0];
			L.flag = 0;
			L.tidx = w->tidx[0];
		} else {
			L.frg1 = w->frgs[1];
			L.frg2 = w->frgs[0];
			L.dir1 = !w->dirs[1][0];
			L.dir2 = w->dirs[0][0];
			L.flag = 0;
			L.tidx = w->tidx[1];
		}
		l = prepare_lnkhash(hash, L, &exists);
		if(exists){
			continue;
		} else {
			*l = L;
		}
	}
	free_weakfnhash(weaks);
	reset_iter_lnkhash(hash);
	while((l = ref_iter_lnkhash(hash))){
		push_edgerefv(g->lrefs, (edge_ref_t){g->lnks->size, 0, 0});
		push_edgerefv(g->lrefs, (edge_ref_t){g->lnks->size, 1, 0});
		push_lnkv(g->lnks, *l);
	}
	free_lnkhash(hash);
	// sort lrefs
	psort_array(g->lrefs->buffer + 1, g->lrefs->size - 1, edge_ref_t, ncpu, num_cmpgt(
		(a.flg? ((g->lnks->buffer[a.idx].frg2 << 1) | !g->lnks->buffer[a.idx].dir2) : ((g->lnks->buffer[a.idx].frg1 << 1) | g->lnks->buffer[a.idx].dir1)),
		(b.flg? ((g->lnks->buffer[b.idx].frg2 << 1) | !g->lnks->buffer[b.idx].dir2) : ((g->lnks->buffer[b.idx].frg1 << 1) | g->lnks->buffer[b.idx].dir1))
		));
	push_edgerefv(g->lrefs, (edge_ref_t){g->lnks->size, 0, 0}); memset(next_ref_lnkv(g->lnks), 0, sizeof(lnk_t)); g->lnks->size --;
	g->lrefs->size --;
	F1.idx = g->lnks->size; F1.flg = 0;
	cnt = 0;
	// update frg->lnks
	for(lst=idx=1;idx<=g->lrefs->size;idx++){
		if(g->lrefs->buffer[idx].flg){
			F2.idx =  g->lnks->buffer[g->lrefs->buffer[idx].idx].frg2;
			F2.flg = !g->lnks->buffer[g->lrefs->buffer[idx].idx].dir2;
		} else {
			F2.idx =  g->lnks->buffer[g->lrefs->buffer[idx].idx].frg1;
			F2.flg =  g->lnks->buffer[g->lrefs->buffer[idx].idx].dir1;
		}
		if(F1.idx == F2.idx && F1.flg == F2.flg) continue;
		if(lst < idx){
			frg = ref_frgv(g->frgs, F1.idx);
			frg->lnks[F1.flg].idx = lst;
			for(x=lst;x+1<(int)idx;x++){
				g->lrefs->buffer[x].next = x + 1;
				if(g->lnks->buffer[g->lrefs->buffer[x].idx].closed == WT_EDGE_CLOSED_NULL) frg->lnks[F1.flg].cnt ++;
			}
			if(g->lnks->buffer[g->lrefs->buffer[x].idx].closed == WT_EDGE_CLOSED_NULL) frg->lnks[F1.flg].cnt ++;
		}
		lst = idx;
		F1 = F2;
	}
	return g->lnks->size - 1;
}

int gen_seq_traces_graph(Graph *g, tracev *path, String *seq){
	trace_t *t1, *t2;
	reg_t *reg, *r1, *r2;
	edge_t *e;
	uint32_t i;
	int inc, found;
	clear_string(seq);
	t1 = NULL;
	for(i=0;i<path->size;i++){
		t2 = ref_tracev(path, i);
		if(t1){
			inc = 0;
			r1 = ref_regv(g->regs, ref_nodev(g->nodes, t1->node)->regs.idx);
			r2 = ref_regv(g->regs, ref_nodev(g->nodes, t2->node)->regs.idx);
			e = ref_edgev(g->edges, t1->edges[t1->dir].idx);
			do {
				inc = 0;
				found = 0;
				while(r1->node == t1->node && r2->node == t2->node){
					if(r1->rid > r2->rid){
						r2 ++;
					} else if(r1->rid < r2->rid){
						r1 ++;
					} else {
						if(r1->beg < r2->beg){
							if(t1->dir ^ r1->dir){ r1++; r2++; continue; }
							inc = r2->beg - r1->end;
							if(inc <= 0) break;
							encap_string(seq, inc);
							seq_basebank(g->dmo->rdseqs, g->dmo->reads->buffer[r1->rid].rdoff + r1->end, inc, seq->string + seq->size);
							seq->size += inc;
						} else {
							if(!(t1->dir ^ r1->dir)){ r1++; r2++; continue; }
							inc = r1->beg - r2->end;
							if(inc <= 0) break;
							encap_string(seq, inc);
							revseq_basebank(g->dmo->rdseqs, g->dmo->reads->buffer[r1->rid].rdoff + r2->end, inc, seq->string + seq->size);
							seq->size += inc;
						}
						inc = 0;
						found = 1; break;
					}
				}
				if(found == 0){ inc = e->off; break; }
			} while(0);
			if(inc > 0){ inc = 0; while(inc++ < e->off) add_char_string(seq, 'N'); }
			else if(inc < 0){
				if(seq->size + inc < 0) seq->size = 0;
				else seq->size += inc;
				seq->string[seq->size] = '\0';
			}
		}
		t1 = t2;
		reg = ref_regv(g->regs, ref_nodev(g->nodes, t1->node)->regs.idx);
		inc = reg->end - reg->beg;
		encap_string(seq, inc);
		if(t1->dir ^ reg->dir) revseq_basebank(g->dmo->rdseqs, g->dmo->reads->buffer[reg->rid].rdoff + reg->beg, inc, seq->string + seq->size);
		else                   seq_basebank(g->dmo->rdseqs, g->dmo->reads->buffer[reg->rid].rdoff + reg->beg, inc, seq->string + seq->size);
		seq->size += inc;
	}
	return seq->size;
}

typedef struct {
	uint64_t rid:26, dir:1, beg:18, end:18, view:1;
} lay_reg_t;
define_list(layregv, lay_reg_t);

typedef struct {
	u4i tidx;
	u8i roff:48, rcnt:16;
} lay_t;
define_list(layv, lay_t);

void gen_lay_regs_core_graph(Graph *g, tracev *path, u4i tidx, layregv *regs){
	trace_t *t1, *t2;
	node_t *n1, *n2;
	reg_t *r1, *r2;
	uint32_t rid, beg, end;
	t1 = ref_tracev(path, tidx);
	t2 = ref_tracev(path, tidx + 1);
	n1 = ref_nodev(g->nodes, t1->node);
	n2 = ref_nodev(g->nodes, t2->node);
	r1 = ref_regv(g->regs, n1->regs.idx);
	r2 = ref_regv(g->regs, n2->regs.idx);
	while(r1->node == t1->node && r2->node == t2->node){
		if(r1->rid > r2->rid){
			r2 ++;
		} else if(r1->rid < r2->rid){
			r1 ++;
		} else {
			rid = r1->rid;
			if(r1->beg < r2->beg){
				if(t1->dir ^ r1->dir){ r1 ++; r2 ++; continue; }
				beg = r1->beg; end = r2->end;
				push_layregv(regs, (lay_reg_t){rid, 0, beg, end, 0});
			} else {
				if(!(t1->dir ^ r1->dir)){ r1 ++; r2 ++; continue; }
				beg = r2->beg; end = r1->end;
				push_layregv(regs, (lay_reg_t){rid, 1, beg, end, 0});
			}
			r1 ++; r2 ++;
		}
	}
}

thread_beg_def(mlay);
Graph *g;
tracev *path;
u4i tidx;
layregv *regs;
thread_end_def(mlay);

thread_beg_func(mlay);
thread_beg_loop(mlay);
if(mlay->path && mlay->tidx != 0xFFFFFFFFU){
	gen_lay_regs_core_graph(mlay->g, mlay->path, mlay->tidx, mlay->regs);
	sort_array(mlay->regs->buffer, mlay->regs->size, lay_reg_t, num_cmpgt(b.end - b.beg, a.end - a.beg));
}
thread_end_loop(mlay);
thread_end_func(mlay);

u8i print_ctgs_graph(Graph *g, u8i uid, u8i beg, u8i end, char *prefix, char *utg_suffix, char *lay_suffix, u4i ncpu, FILE *log){
	FILE *o_utg, *o_lay;
	layv *lays;
	layregv *regs;
	tracev *path;
	String *seqs;
	u4v *lens;
	trace_t *t1, *t2;
	lay_t *lay;
	lay_reg_t *reg;
	u8i i;
	u4i j, c;
	int offset;
	thread_preprocess(mlay);
	o_utg = open_file_for_write(prefix, utg_suffix, 1);
	o_lay = open_file_for_write(prefix, lay_suffix, 1);
	seqs = init_string(1024);
	lens = init_u4v(1024);
	lays = init_layv(32);
	regs = init_layregv(32);
	thread_beg_init(mlay, ncpu);
	mlay->g = g;
	mlay->path = NULL;
	mlay->tidx = 0xFFFFFFFFU;
	mlay->regs = init_layregv(32);
	thread_end_init(mlay);
	for(i=beg;i<end;i++){
		path = (tracev*)get_vplist(g->ctgs, i);
		if(path->size < 2) continue;
		clear_layv(lays);
		clear_layregv(regs);
		clear_string(seqs);
		thread_beg_iter(mlay);
		mlay->path = path;
		mlay->tidx = 0xFFFFFFFFU;
		clear_layregv(mlay->regs);
		thread_end_iter(mlay);
		for(j=0;j<=path->size + ncpu;j++){
			if(j + 1 < path->size){
				thread_wait_one(mlay);
			} else {
				thread_wait_next(mlay);
			}
			if(mlay->tidx != 0xFFFFFFFFU){
				if(mlay->regs->size == 0){
					fprintf(stderr, " -- N%llu(%c) -> N%llu(%c) has no read path in %s -- %s:%d --\n", path->buffer[mlay->tidx].node, "+-"[path->buffer[mlay->tidx].dir], path->buffer[mlay->tidx + 1].node, "+-"[path->buffer[mlay->tidx + 1].dir], __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				}
				lay = next_ref_layv(lays);
				lay->tidx = mlay->tidx;
				lay->roff = regs->size;
				lay->rcnt = mlay->regs->size;
				append_layregv(regs, mlay->regs);
				clear_layregv(mlay->regs);
				mlay->tidx = 0xFFFFFFFFU;
			}
			if(j + 1 < path->size){
				mlay->tidx = j;
				thread_wake(mlay);
			}
		}
		if(lays->size == 0) continue;
		uid ++;
		sort_array(lays->buffer, lays->size, lay_t, num_cmpgt(a.tidx, b.tidx));
		fprintf(o_lay, ">ctg%llu nodes=%llu\n", uid, (u8i)path->size);
		fprintf(log, "OUTPUT\tctg%d -> ctg%d\n", (int)i, (int)uid);
		offset = 0;
		for(j=0;j<lays->size;j++){
			lay = ref_layv(lays, j);
			t1 = ref_tracev(path, lay->tidx);
			t2 = ref_tracev(path, lay->tidx + 1);
			if(j){
				if(seqs->size < (int)g->reglen){
					offset = 0;
				} else {
					offset = seqs->size - g->reglen;
				}
				trunc_string(seqs, offset);
			}
			fprintf(o_lay, "E\t%d\tN%llu\t%c\tN%llu\t%c\n", offset, t1->node, "+-"[t1->dir], t2->node, "+-"[t2->dir]);
			for(c=0;c<lay->rcnt;c++){
				reg = ref_layregv(regs, lay->roff + c);
				fprintf(o_lay, "%c\t%s\t%c\t%d\t%d\t", "Ss"[reg->view], g->dmo->reads->buffer[reg->rid].tag, "+-"[reg->dir], reg->beg, reg->end - reg->beg);
				if(reg->dir){
					print_revseq_basebank(g->dmo->rdseqs, g->dmo->reads->buffer[reg->rid].rdoff + reg->beg, reg->end - reg->beg, o_lay);
				} else {
					print_seq_basebank(g->dmo->rdseqs, g->dmo->reads->buffer[reg->rid].rdoff + reg->beg, reg->end - reg->beg, o_lay);
				}
				fprintf(o_lay, "\n");
			}
			c = lay->rcnt / 2;
			reg = ref_layregv(regs, lay->roff + c);
			encap_string(seqs, reg->end - reg->beg);
			if(reg->dir){
				revseq_basebank(g->dmo->rdseqs, g->dmo->reads->buffer[reg->rid].rdoff + reg->beg, reg->end - reg->beg, seqs->string + seqs->size);
			} else {
				seq_basebank(g->dmo->rdseqs, g->dmo->reads->buffer[reg->rid].rdoff + reg->beg, reg->end - reg->beg, seqs->string + seqs->size);
			}
			seqs->size += reg->end - reg->beg;
			seqs->string[seqs->size] = '\0';
		}
		if(seqs->size){
			fprintf(o_utg, ">ctg%llu len=%d nodes=%llu beg=N%llu end=N%llu\n", uid, seqs->size, (u8i)path->size,
				path->buffer[0].node, path->buffer[path->size - 1].node);
			print_pretty_seq(o_utg, seqs, 100);
			push_u4v(lens, seqs->size);
		}
	}
	fclose(o_lay);
	fclose(o_utg);
	thread_beg_close(mlay);
	free_layregv(mlay->regs);
	thread_end_close(mlay);
	fprintf(hzm_debug_out, "[%s] contigs: ", date()); num_n50(lens, hzm_debug_out); fprintf(hzm_debug_out, "\n");
	free_string(seqs);
	free_u4v(lens);
	free_layv(lays);
	free_layregv(regs);
	return uid;
}

uint32_t print_traces_graph(Graph *g, tracev *path, FILE *out){
	trace_t *t1, *t2;
	node_t *n1, *n2;
	reg_t *r1, *r2;
	edge_ref_t *f;
	edge_t *e;
	e = NULL;
	edge_t *wushigang = e;
	edge_t *tmp = wushigang;
	wushigang = tmp;
	int offset, fst;
	uint64_t j, beg, end;
	uint32_t i, rid;
	if(path->size < 2) return 0;
	offset = 0;
	t1 = ref_tracev(path, 0);
	for(i=1;i<path->size;i++){
		t2 = ref_tracev(path, i);
		{
			n1 = ref_nodev(g->nodes, t1->node);
			n2 = ref_nodev(g->nodes, t2->node);
			f  = t1->edges + t1->dir;
			e  = ref_edgev(g->edges, f->idx);
			r1 = ref_regv(g->regs, n1->regs.idx);
			r2 = ref_regv(g->regs, n2->regs.idx);
			fprintf(out, "E\t%d\tN%llu\t%c\tN%llu\t%c\n", offset, t1->node, "+-"[t1->dir], t2->node, "+-"[t2->dir]);
			fst = 1;
			while(r1->node == t1->node && r2->node == t2->node){
				if(r1->rid > r2->rid){
					r2 ++;
				} else if(r1->rid < r2->rid){
					r1 ++;
				} else {
					rid = r1->rid;
					if(r1->beg < r2->beg){
						if(t1->dir ^ r1->dir){ r1 ++; r2 ++; continue; }
						beg = r1->beg; end = r2->end;
						fprintf(out, "S\t%s\t", g->dmo->reads->buffer[rid].tag);
						fprintf(out, "+\t%d\t%d\t", (int)beg, (int)(end - beg));
						beg += g->dmo->reads->buffer[rid].rdoff;
						end += g->dmo->reads->buffer[rid].rdoff;
						for(j=beg;j<end;j++){
							fputc(bit_base_table[bits2bit(g->dmo->rdseqs->bits, j)], out);
						}
					} else {
						if(!(t1->dir ^ r1->dir)){ r1 ++; r2 ++; continue; }
						beg = r2->beg; end = r1->end;
						fprintf(out, "S\t%s\t", g->dmo->reads->buffer[rid].tag);
						fprintf(out, "-\t%d\t%d\t", (int)beg, (int)(end - beg));
						beg += g->dmo->reads->buffer[rid].rdoff;
						end += g->dmo->reads->buffer[rid].rdoff;
						for(j=end;j>beg;j--){
							fputc(bit_base_table[bits2revbit(g->dmo->rdseqs->bits, (j-1))], out);
						}
					}
					fprintf(out, "\n");
					if(fst){
						offset += end - beg; fst = 0;
					}
					r1 ++; r2 ++;
				}
			}
		}
		t1 = t2;
	}
	return offset;
}

uint64_t print_utgs_graph(Graph *g, char *prefix, char *utg, char *lay){
	FILE *o_utg, *o_lay, *files[4];
	tracev *path;
	String *seq;
	char *str;
	uint64_t i, uid, cnt, tot;
	char ch;
	int beg, end;
	files[0] = open_file_for_write(prefix, utg, 1);
	str = catstr(2, utg, ".filtered");
	files[1] = open_file_for_write(prefix, str, 1);
	free(str);
	files[2] = open_file_for_write(prefix, lay, 1);
	str = catstr(2, lay, ".filtered");
	files[3] = open_file_for_write(prefix, str, 1);
	free(str);
	seq = init_string(1024);
	tot = cnt = 0;
	for(i=uid=0;i<g->utgs->size;i++){
		path = (tracev*)get_vplist(g->utgs, i);
		if(gen_seq_traces_graph(g, path, seq) < g->min_ctg_len || (int)path->size < g->min_ctg_nds){
			o_utg = files[1];
			o_lay = files[3];
		} else {
			o_utg = files[0];
			o_lay = files[2];
			cnt ++;
			tot += seq->size;
		}
		uid ++;
		fprintf(o_utg, ">utg%llu len=%d nodes=%llu beg=N%llu end=N%llu\n", (unsigned long long)uid, seq->size, (unsigned long long)path->size,
			path->buffer[0].node, path->buffer[path->size - 1].node);
		for(beg=0;beg<seq->size;beg+=100){
			end = beg + 100;
			if(end > seq->size) end = seq->size;
			ch = seq->string[end];
			seq->string[end] = '\0';
			fprintf(o_utg, "%s\n", seq->string + beg);
			seq->string[end] = ch;
		}
		fprintf(o_lay, ">utg%llu len=%d nodes=%llu\n", (unsigned long long)uid, seq->size, (unsigned long long)path->size);
		print_traces_graph(g, path, o_lay);
	}
	free_string(seq);
	fprintf(hzm_debug_out, "[%s] %llu unitigs (>= %d bp), total %llu bp\n", date(), (unsigned long long)cnt, g->min_ctg_len, (unsigned long long)tot);
	fclose(files[0]);
	fclose(files[1]);
	fclose(files[2]);
	fclose(files[3]);
	return uid;
}

uint64_t _print_ctgs_graph(Graph *g, char *prefix, char *utg, char *lay){
	FILE *o_utg, *o_lay;
	tracev *path;
	String *seq;
	u4v *lens;
	uint64_t i, uid;
	char ch;
	int beg, end;
	o_utg = open_file_for_write(prefix, utg, 1);
	o_lay = open_file_for_write(prefix, lay, 1);
	seq = init_string(1024);
	lens = init_u4v(1024);
	for(i=uid=0;i<g->ctgs->size;i++){
		path = (tracev*)get_vplist(g->ctgs, i);
		if(path->size < 2) continue;
		gen_seq_traces_graph(g, path, seq);
		push_u4v(lens, seq->size);
		uid ++;
		fprintf(o_utg, ">ctg%llu len=%d nodes=%llu beg=N%llu end=N%llu\n", (unsigned long long)uid, seq->size, (unsigned long long)path->size,
			path->buffer[0].node, path->buffer[path->size - 1].node);
		for(beg=0;beg<seq->size;beg+=100){
			end = beg + 100;
			if(end > seq->size) end = seq->size;
			ch = seq->string[end];
			seq->string[end] = '\0';
			fprintf(o_utg, "%s\n", seq->string + beg);
			seq->string[end] = ch;
		}
		fprintf(o_lay, ">ctg%llu len=%d nodes=%llu\n", (unsigned long long)uid, seq->size, (unsigned long long)path->size);
		print_traces_graph(g, path, o_lay);
	}
	free_string(seq);
	fclose(o_utg);
	fclose(o_lay);
	fprintf(hzm_debug_out, "[%s] contigs: ", date()); num_n50(lens, hzm_debug_out); fprintf(hzm_debug_out, "\n");
	free_u4v(lens);
	return uid;
}

uint64_t print_dot_full_graph(Graph *g, FILE *out){
	node_t *n;
	reg_t *r, *rr;
	edge_ref_t *f;
	edge_t *e;
	unsigned long long i, idx;
	uint32_t j, k, max;
	fprintf(out, "digraph {\n");
	fprintf(out, "node [shape=record]\n");
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		//if(n->closed) continue;
		r = NULL; max = 0;
		for(j=0;j<n->regs.cnt;j++){
			rr = ref_regv(g->regs, n->regs.idx + j);
			if(g->reads->buffer[rr->rid].regs.cnt > max){
				r = rr;
				max = g->reads->buffer[rr->rid].regs.cnt;
			}
		}
		if(r == NULL) continue;
		fprintf(out, "N%llu [label=\"{N%llu %d | %s | %c_%d_%d}\"]\n", i, i, n->regs.cnt, g->dmo->reads->buffer[r->rid].tag, "FR"[r->dir], r->beg, r->end - r->beg);
	}
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		//if(n->closed) continue;
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				idx = f->next;
				e = g->edges->buffer + f->idx;
				//if(e->closed) continue;
				if(f->flg){
					fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s%s]\n", i, (unsigned long long)e->node1, "+-"[k], "+-"[!e->dir1], e->cov, e->off, colors[k][!e->dir1], e->closed? " style=dashed" : "");
				} else {
					fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s%s]\n", i, (unsigned long long)e->node2, "+-"[k], "+-"[e->dir2], e->cov, e->off, colors[k][e->dir2], e->closed? " style=dashed" : "");
				}
			}
		}
	}
	fprintf(out, "}\n");
	return 0;
}

uint64_t print_dot_graph(Graph *g, FILE *out){
	node_t *n;
	reg_t *r, *rr;
	edge_ref_t *f;
	edge_t *e;
	unsigned long long i, idx;
	uint32_t j, k, max;
	fprintf(out, "digraph {\n");
	fprintf(out, "node [shape=record]\n");
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		r = NULL; max = 0;
		for(j=0;j<n->regs.cnt;j++){
			rr = ref_regv(g->regs, n->regs.idx + j);
			if(g->reads->buffer[rr->rid].regs.cnt > max){
				r = rr;
				max = g->reads->buffer[rr->rid].regs.cnt;
			}
		}
		if(r == NULL) continue;
		fprintf(out, "N%llu [label=\"{N%llu %d | %s | %c_%d_%d}\"]\n", i, i, n->regs.cnt, g->dmo->reads->buffer[r->rid].tag, "FR"[r->dir], r->beg, r->end - r->beg);
	}
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				idx = f->next;
				e = g->edges->buffer + f->idx;
				if(e->closed) continue;
				if(f->flg){
					fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s]\n", i, (unsigned long long)e->node1, "+-"[k], "+-"[!e->dir1], e->cov, e->off, colors[k][!e->dir1]);
				} else {
					fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s]\n", i, (unsigned long long)e->node2, "+-"[k], "+-"[e->dir2], e->cov, e->off, colors[k][e->dir2]);
				}
			}
		}
	}
	fprintf(out, "}\n");
	return 0;
}

static uint64_t local_dot_node = 1;
static uint32_t local_dot_step = 10;

void get_subgraph_nodes_graph(Graph *g, ptrrefhash *nodes, u8v *stack, uint16_t max_step, uint32_t closed_val){
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	ptr_ref_t *p, *pp;
	uint64_t nid, idx;
	uint32_t k, cnt;
	int exists;
	clear_u8v(stack);
	reset_iter_ptrrefhash(nodes);
	while((p = ref_iter_ptrrefhash(nodes))){
		p->cnt = 0;
		push_u8v(stack, p->idx);
	}
	while(stack->size){
		p = get_ptrrefhash(nodes, (ptr_ref_t){stack->buffer[--stack->size], 0});
		if(p->cnt >> 16) continue;
		if((p->cnt & 0xFFFF) >= max_step) continue;
		n = ref_nodev(g->nodes, p->idx);
		cnt = p->cnt;
		p->cnt |= 1U << 16;
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				idx = f->next;
				e = ref_edgev(g->edges, f->idx);
				if(e->closed >= closed_val) continue;
				nid = f->flg? e->node1 : e->node2;
				pp = prepare_ptrrefhash(nodes, (ptr_ref_t){nid, 0}, &exists);
				if(exists) continue;
				pp->idx = nid; pp->cnt = cnt + 1;
				push_u8v(stack, nid);
			}
		}
	}
}

uint64_t print_local_dot_graph(Graph *g, FILE *out){
	ptrrefhash *hash;
	u8v *stack;
	ptr_ref_t *p;
	node_t *n;
	reg_t *r, *rr;
	edge_ref_t *f;
	edge_t *e;
	unsigned long long i, idx;
	uint32_t j, k, max;
	hash = init_ptrrefhash(1023);
	stack = init_u8v(32);
	put_ptrrefhash(hash, (ptr_ref_t){local_dot_node, 0});
	get_subgraph_nodes_graph(g, hash, stack, local_dot_step, 1);
	fprintf(out, "digraph {\n");
	fprintf(out, "node [shape=record]\n");
	reset_iter_ptrrefhash(hash);
	while((p = ref_iter_ptrrefhash(hash))){
		i = p->idx;
		n = ref_nodev(g->nodes, i);
		r = NULL; max = 0;
		for(j=0;j<n->regs.cnt;j++){
			rr = ref_regv(g->regs, n->regs.idx + j);
			if(g->reads->buffer[rr->rid].regs.cnt > max){
				r = rr;
				max = g->reads->buffer[rr->rid].regs.cnt;
			}
		}
		if(r == NULL) continue;
		fprintf(out, "N%llu [label=\"{N%llu %d | %s | %c_%d_%d}\"]\n", i, i, n->regs.cnt, g->dmo->reads->buffer[r->rid].tag, "FR"[r->dir], r->beg, r->end - r->beg);
	}
	reset_iter_ptrrefhash(hash);
	while((p = ref_iter_ptrrefhash(hash))){
		i = p->idx;
		n = ref_nodev(g->nodes, i);
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				idx = f->next;
				e = g->edges->buffer + f->idx;
				if(e->closed) continue;
				if(f->flg){
					//if(!exists_ptrrefhash(hash, (ptr_ref_t){e->node1, 0})) continue;
					fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s]\n", i, (unsigned long long)e->node1, "+-"[k], "+-"[!e->dir1], e->cov, e->off, colors[k][!e->dir1]);
				} else {
					//if(!exists_ptrrefhash(hash, (ptr_ref_t){e->node2, 0})) continue;
					fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s]\n", i, (unsigned long long)e->node2, "+-"[k], "+-"[e->dir2], e->cov, e->off, colors[k][e->dir2]);
				}
			}
		}
	}
	fprintf(out, "}\n");
	return 0;
}

uint64_t print_nodes_graph(Graph *g, FILE *out){
	node_t *n;
	reg_t *r;
	unsigned long long i;
	uint32_t j;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		fprintf(out, "N%llu\t%u", i, n->regs.cnt);
		for(j=0;j<n->regs.cnt;j++){
			r = ref_regv(g->regs, n->regs.idx + j);
			fprintf(out, "\t%s_%c_%d_%d", g->dmo->reads->buffer[r->rid].tag, "FR"[r->dir], r->beg, r->end - r->beg);
		}
		fprintf(out, "\n");
	}
	return i;
}

uint64_t print_reads_graph(Graph *g, FILE *out){
	read_t *rd;
	reg_t  *r;
	uint64_t idx;
	uint32_t i;
	for(i=0;i<g->dmo->n_rd;i++){
		rd = ref_readv(g->reads, i);
		fprintf(out, "%s\t%d\t%u", g->dmo->reads->buffer[i].tag, g->dmo->reads->buffer[i].rdlen, rd->regs.cnt);
		idx = rd->regs.idx;
		while(idx){
			r = ref_regv(g->regs, idx);
			fprintf(out, "\tN%llu%s:%c_%d_%d", (unsigned long long)r->node, g->nodes->buffer[r->node].closed? "*": "", "FR"[r->dir], r->beg, r->end - r->beg);
			idx = r->read_link;
		}
		fprintf(out, "\n");
	}
	return i;
}

uint64_t print_frgs_nodes_graph(Graph *g, FILE *out){
	frg_t *frg;
	trace_t *t;
	node_t *n;
	u4i i, j;
	for(i=0;i<g->frgs->size;i++){
		frg = ref_frgv(g->frgs, i);
		if(frg->closed) continue;
		fprintf(out, "F%u\t%d\t%d\t%u\t%u", i, frg->length, frg->len, frg->tcnt, frg->ty - frg->tx);
		for(j=0;j<frg->tcnt;j++){
			t = ref_tracev(g->traces, frg->toff + j);
			if(j < frg->tx || j >= frg->ty){
				n = ref_nodev(g->nodes, t->node);
				if(n->rep_idx == MAX_REP_IDX){
					fprintf(out, "\tn%llu:%c:%d:::%d", t->node, "+-"[t->dir], t->off, t->cov);
				} else {
					fprintf(out, "\tn%llu:%c:%d:F%llu:%c:%d", t->node, "+-"[t->dir], t->off, (u8i)n->rep_idx, "+-"[n->rep_dir], t->cov);
				}
			} else {
				fprintf(out, "\tN%llu:%c:%d", t->node, "+-"[t->dir], t->off);
			}
		}
		fprintf(out, "\n");
	}
	return i;
}

uint64_t print_frgs_dot_graph(Graph *g, FILE *out){
	frg_t *frg;
	trace_t *t1, *t2;
	node_t *n1, *n2;
	reg_t *r1, *r2, *rr;
	edge_ref_t *f;
	lnk_t *e;
	unsigned long long i, idx;
	uint32_t j, k, max;
	fprintf(out, "digraph {\n");
	fprintf(out, "node [shape=record]\n");
	for(i=0;i<g->frgs->size;i++){
		frg = ref_frgv(g->frgs, i);
		if(frg->closed) continue;
		//if(frg->ty - frg->tx < (u4i)g->min_ctg_nds) continue;
		t1 = ref_tracev(g->traces, frg->toff + frg->tx);
		t2 = ref_tracev(g->traces, frg->toff + frg->ty - 1);
		n1 = ref_nodev(g->nodes, t1->node);
		n2 = ref_nodev(g->nodes, t2->node);
		r1 = NULL; max = 0;
		for(j=0;j<n1->regs.cnt;j++){
			rr = ref_regv(g->regs, n1->regs.idx + j);
			if(g->reads->buffer[rr->rid].regs.cnt > max){
				r1 = rr;
				max = g->reads->buffer[rr->rid].regs.cnt;
			}
		}
		if(r1 == NULL) continue;
		r2 = NULL; max = 0;
		for(j=0;j<n2->regs.cnt;j++){
			rr = ref_regv(g->regs, n2->regs.idx + j);
			if(g->reads->buffer[rr->rid].regs.cnt > max){
				r2 = rr;
				max = g->reads->buffer[rr->rid].regs.cnt;
			}
		}
		if(r2 == NULL) continue;
		fprintf(out, "F%llu [label=\"{F%llu %u %u/%u | { {N%llu:%c | %s | %c_%d_%d} | {N%llu:%c | %s | %c_%d_%d}}}\"]\n", i, i, frg->ty - frg->tx, frg->len, frg->length,
			t1->node, "+-"[t1->dir], g->dmo->reads->buffer[r1->rid].tag, "FR"[r1->dir], r1->beg, r1->end - r1->beg,
			t2->node, "+-"[t2->dir], g->dmo->reads->buffer[r2->rid].tag, "FR"[r2->dir], r2->beg, r2->end - r2->beg);
	}
	for(i=0;i<g->frgs->size;i++){
		frg = ref_frgv(g->frgs, i);
		if(frg->closed) continue;
		//if(frg->ty - frg->tx < (u4i)g->min_ctg_nds) continue;
		for(k=0;k<2;k++){
			idx = frg->lnks[k].idx;
			while(idx){
				f = ref_edgerefv(g->lrefs, idx);
				idx = f->next;
				e = g->lnks->buffer + f->idx;
				if(e->closed) continue;
				if(f->flg){
					//if(g->frgs->buffer[e->frg1].ty - g->frgs->buffer[e->frg1].tx < (u4i)g->min_ctg_nds) continue;
					fprintf(out, "F%llu -> F%llu [label=\"%c%c:%d:%d\" color=%s style=%s]\n", i, (u8i)e->frg1, "+-"[k], "+-"[!e->dir1], e->cov, e->off, colors[k][!e->dir1], e->weak? "dashed" : "solid");
				} else {
					//if(g->frgs->buffer[e->frg2].ty - g->frgs->buffer[e->frg2].tx < (u4i)g->min_ctg_nds) continue;
					fprintf(out, "F%llu -> F%llu [label=\"%c%c:%d:%d\" color=%s style=%s]\n", i, (u8i)e->frg2, "+-"[k], "+-"[e->dir2], e->cov, e->off, colors[k][e->dir2], e->weak? "dashed" : "solid");
				}
			}
		}
	}
	fprintf(out, "}\n");
	return 0;
}

typedef uint64_t (*graph_print_func)(Graph *g, FILE *out);

uint64_t generic_print_graph(Graph *g, graph_print_func func, char *prefix, char *suffix){
	FILE *out;
	char *file;
	uint64_t cnt;
	{
		fprintf(hzm_debug_out, "[%s] output %s%s\n", date(), prefix, suffix? suffix : ""); fflush(hzm_debug_out);
		file = malloc(strlen(prefix) + (suffix? strlen(suffix) : 0) + 2);
		sprintf(file, "%s%s", prefix, suffix? suffix : "");
		out = fopen(file, "w");
		cnt = func(g, out);
		fclose(out);
		free(file);
	}
	return cnt;
}

int usage(){
	printf(
	"WTDBG: De novo assembler for long noisy sequences\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.1.006\n"
#ifdef TIMESTAMP
	"Compiled: %u\n"
#endif
	"Usage: wtdbg [options]\n"
	"Options:\n"
	" -t <int>    Number of threads, 0: all cores, [1]\n"
	" -i <string> Long reads sequences file, + *\n"
	" -I <string> Accurate contig sequences file, +\n"
	" -r          Search nodes on accurate contig sequences only\n"
	" -T <int>    Limitation of total base-pairs in Mb, usually for test [0]\n"
	//" -G <float>  Split the reads into parts of <-G> G bp for indexing, to save memory [0]\n"
	" -R          Rapid subsampling instead of aligning all intervals\n"
	" -j <int>    Expected length of node, [1000] bp\n"
	" -J <int>    Distance of next node start, 0: the same as <-j> [0] bp\n"
	" -X <float>  Max overlapped fraction of two adjacent nodes in a read, [0.2]\n"
	//" -Y <float>  Min fraction of overlapped length between two nodes' instances to be merged into one node, [0.9]\n"
	" -e <int>    Min cov of edges, [3]\n"
	" -n <int>    Min cov of nodes, 0: the same as <-e> [0]\n"
	" -N <int>    Max cov of nodes, remove too-high coverage nodes [200]\n"
	" -L <string> Load graph from \"<your_prefix>.nodes\" instead of build, [NULL]\n"
	" -b <int>    Max steps of bubble, [20]\n"
	" -d <int>    Max steps of tips, [5]\n"
	" -o <string> Prefix of output files, *\n"
	" -f          Force overwrite\n"
	//" -c <int>    Min length of contig, [10000]\n"
	//" -C <int>    Min nodes  of contig, [5]\n"
	" -v          Verbose\n"
	" -0          Output as less files as it can\n"
	" -----------parameters of node alignment-----------\n"
	" -H          Trun on homopolymer compression\n"
	" -k <int>    Kmer size, 5 <= value <= 32, [15]\n"
	" -W <int>    Max size of insertion in the middle of kmer when querying, [0]\n"
	"             PART1|ins|PART2, PART1 + PART2 = ksize, PART2 = ksize/2, ins <= max_kgap, max_kgap + ksize <= 32\n"
	" -K <float>  Filter high frequency kmers, maybe repetitive, [0.05]\n"
	"             if K >= 1, take the integer value as cutoff\n"
	"             else, mask the top fraction part high frequency kmers\n"
	" -E <int>    Min kmer frequency, [2]\n"
	" -S <float>  Subsampling kmers, 1/(int((<-S> * 100) %% 100) * int(<-S>)) kmers are indexed, [1.01]\n"
	"             Decimal part is used for minimizer window size\n"
	"             Fraction part is used for HASH based subsampling\n"
	" -x <int>    Intra-block: Max gap  [256]\n"
	" -y <int>    Intra-block: Max deviation [128]\n"
	" -z <int>    Intra-block: Min kmer [3]\n"
	" -w <int>    Inter-block: deviation penalty [1.0]\n"
	" -u <int>    Inter-block: gap penalty [0.1]\n"
	" -l <float>  Min fraction of alignment:aligned / <expected_node_length>, [0.7]\n"
	" -m <float>  Min fraction of alignment:matched / <expected_node_length>, [0.1]\n"
	" -s <float>  Max length variation of two aligned fragments, [0.2]\n"
	" ------------parameters of refine alignment-------\n"
	" -Z <int>    Try to align selected nodes against potential matched reads, [0]\n"
	"             5 <= value <= 27, suggested 10, like -k, set to 0 to disable refinement\n"
	"             All align parameters after -Z is set to realignment parameters\n"
	"             default: -Z 0 -K 0 -S 1.01 -x 256 -y 64 -z 3 -w 1.0 -u 0.1 -l 0.70 -m 0.20 -s 0.2\n"
	" -------------------------------------------------\n"
	" -M <int>    Test functions, +\n"
	"             1: don't skip to align a new node when it conflicts with previous repetitive nodes\n"
	"             2: output alignments to <prefix>.alignments\n"
	"             3: filter tip-like nodes\n"
	"             4: don't keep low coverage edges in graph\n"
	"             5: don't filter nodes are likely to be repetitive by local subgraph analysis\n"
	//"             6: use read path detachment in solving repeats\n"
	//"             7: don't cut tips\n"
	"\n"
#ifdef TIMESTAMP
	, TIMESTAMP
#endif
	);
	return 1;
}

int main(int argc, char **argv){
	rnk_ref_t wushigang1 = RNK_REF_NULL;
	rnk_ref_t tmp1 = wushigang1;
	wushigang1 = tmp1;
	vec_ref_t wushigang2 = VEC_REF_NULL;
	vec_ref_t tmp2 = wushigang2;
	wushigang2 = tmp2;
	obj_desc_t wushigang = rdregv_obj_desc;
	obj_desc_t tmp = wushigang;
	wushigang = tmp;
	wushigang = rdrepv_obj_desc;
	wushigang = regv_obj_desc;
	wushigang = xyv_obj_desc;
	wushigang = vecrefv_obj_desc;
	wushigang = ptrrefv_obj_desc;
	wushigang = rnkrefv_obj_desc;
	wushigang = ptrrefhash_obj_desc;
	wushigang = edgev_obj_desc;
	wushigang = edgehash_obj_desc;
	wushigang = edgerefv_obj_desc;
	wushigang = nodev_obj_desc;
	wushigang = readv_obj_desc;
	wushigang = tracev_obj_desc;
	wushigang = lnkv_obj_desc;
	wushigang = lnkhash_obj_desc;
	wushigang = frgv_obj_desc;
	wushigang = pathv_obj_desc;
	wushigang = edgeoffv_obj_desc;
	wushigang = subnodehash_obj_desc;
	wushigang = subedgev_obj_desc;
	wushigang = btv_obj_desc;
	wushigang = frgbtv_obj_desc;
	wushigang = sgheapv_obj_desc;
	wushigang = sgtipv_obj_desc;
	wushigang = subrdv_obj_desc;
	wushigang = subrdhash_obj_desc;
	wushigang = weakfnhash_obj_desc;
	wushigang = layregv_obj_desc;
	wushigang = layv_obj_desc;
	Graph *g;
	DMOPar pars[2];
	DMO *dmo;
	FileReader *fr;
	Sequence *seq;
	cplist *pbs, *ngs;
	FILE *evtlog;
	char *prefix, *load;
	int len, tag_size;
	unsigned long long tot_bp, cnt, bub, tip, rep, yarn, max_bp, max_idx_bp, nfix;
	uint32_t i;
	int c, val, ncpu, only_fix, node_cov, max_node_cov, edge_cov, store_low_cov_edge, reglen, regoff, bub_step, tip_step, rep_step;
	int frgtip_len;
	int less_out, tip_like, cut_tip, rep_filter, out_alns, cnn_filter, log_rep, rep_detach;
	int min_ctg_len, min_ctg_nds, max_trace_end, overwrite, len_order, node_order, align_mode, zalign, fast_mode;
	float regovl, node_mrg, fval;
	BEG_STAT_PROC_INFO(stderr, argc, argv);
	pbs = init_cplist(4);
	ngs = init_cplist(4);
	ncpu = 1;
	fast_mode = 0;
	max_bp = 0;
	max_idx_bp = 0LLU * 1000 * 1000 * 1000; // unlimited
	reglen = 1000;
	regoff = 0;
	regovl = 0.2;
	node_mrg = 0.9;
	only_fix = 0;
	node_cov = 0;
	max_node_cov = 200;
	edge_cov = 3;
	store_low_cov_edge = 1;
	load = NULL;
	bub_step = 20;
	tip_step = 5;
	rep_step = 0;
	max_trace_end = 5;
	frgtip_len = 50000;
	prefix = NULL;
	overwrite = 0;
	less_out = 0;
	rep_filter = 1;
	tip_like = 0;
	cut_tip = 1;
	cnn_filter = 1;
	log_rep = 1;
	rep_detach = 0;
	min_ctg_len = 10000;
	min_ctg_nds = 5;
	pars[0].min_mat = 100; // 0.1 * 1000
	pars[0].min_aln = 700; // 0.7 * 1000
	pars[0].min_sm  = 0.0;
	pars[0].aln_var = 0.2;
	pars[0].hk = 0;
	pars[0].ksize = 15;
	pars[0].kmax  = 0;
	pars[0].ktop  = 0.05;
	pars[0].kmin  = 2;
	pars[0].kgap  = 0;
	pars[0].max_hit = 0;
	pars[0].hzmh_kmer_mod = 1 * HZMH_KMER_MOD;
	pars[0].hzmh_kmer_win = 1;
	pars[0].xvar = 256;
	pars[0].yvar = 128;
	pars[0].zmin = 3;
	pars[0].max_overhang = 600;
	pars[0].deviation_penalty = 1.0;
	pars[0].gap_penalty = 0.1;
	pars[1].min_mat = 200;
	pars[1].min_aln = 700;
	pars[1].min_sm  = 0.0;
	pars[1].aln_var = 0.2;
	pars[1].hk = 0;
	pars[1].ksize = 0;
	pars[1].kmax  = 0;
	pars[1].ktop  = 0.00001;
	pars[1].kmin  = 1;
	pars[1].kgap  = 0;
	pars[1].max_hit = 0;
	pars[1].hzmh_kmer_mod = HZMH_KMER_MOD;
	pars[1].hzmh_kmer_win = 1;
	pars[1].xvar = 256;
	pars[1].yvar = 64;
	pars[1].zmin = 3;
	pars[1].max_overhang = 600;
	pars[1].deviation_penalty = 1.0;
	pars[1].gap_penalty = 0.1;
	zalign = 0;
	len_order = 1;
	node_order = 0;
	align_mode = 1;
	out_alns = 0;
	while((c = getopt(argc, argv, "ht:i:I:rT:G:Rj:J:X:Y:n:N:e:L:b:d:o:fc:C:0Hk:W:K:E:S:x:y:z:w:u:l:m:s:Z:F:vM:")) != -1){
		switch(c){
			case 't': ncpu = atoi(optarg); break;
			case 'i': push_cplist(pbs, optarg); break;
			case 'I': push_cplist(ngs, optarg); len_order = 0; break;
			case 'r': only_fix = 1; break;
			case 'R': fast_mode = 1; break;
			case 'T': max_bp = atol(optarg); break;
			case 'G': max_idx_bp = atof(optarg) * 1000LLU * 1000 * 1000; break;
			case 'j': reglen = atoi(optarg); break;
			case 'J': regoff = atoi(optarg); break;
			case 'X': regovl = atof(optarg); break;
			case 'Y': node_mrg = atof(optarg); break;
			case 'n': node_cov = atoi(optarg); break;
			case 'N': max_node_cov = atoi(optarg); break;
			case 'e': edge_cov = atoi(optarg); break;
			case 'L': load = optarg; break;
			case 'b': bub_step = atoi(optarg); break;
			case 'd': tip_step = atoi(optarg); break;
			//case 'r': rep_step = atoi(optarg); break;
			case 'o': prefix = optarg; break;
			case 'f': overwrite = 1; break;
			case 'c': min_ctg_len = atoi(optarg); break;
			case 'C': min_ctg_nds = atoi(optarg); break;
			case '0': less_out = 1; break;
			case 'H': pars[zalign].hk = 1; break;
			case 'k': pars[zalign].ksize = atoi(optarg); break;
			case 'W': pars[zalign].kgap  = atoi(optarg); break;
			case 'K': pars[zalign].kmax = atoi(optarg); break;
			case 'E': pars[zalign].kmin = atoi(optarg); break;
			case 'S': fval = atof(optarg);
					pars[zalign].hzmh_kmer_win = (int)fval;
					if(pars[zalign].hzmh_kmer_win < 1) pars[zalign].hzmh_kmer_win = 1;
					pars[zalign].hzmh_kmer_mod = (((int)(fval * 100)) % 100) * HZMH_KMER_MOD;
					if(pars[zalign].hzmh_kmer_mod < HZMH_KMER_MOD) pars[zalign].hzmh_kmer_mod = HZMH_KMER_MOD;
					break;
			case 'x': pars[zalign].xvar = atoi(optarg); break;
			case 'y': pars[zalign].yvar = atoi(optarg); break;
			case 'z': pars[zalign].zmin = atoi(optarg); break;
			case 'w': pars[zalign].deviation_penalty = atof(optarg); break;
			case 'u': pars[zalign].gap_penalty = atof(optarg); break;
			case 'l': pars[zalign].min_aln = atof(optarg) * 1000; break;
			case 'm': pars[zalign].min_mat = atof(optarg) * 1000; break;
			case 's': pars[zalign].aln_var = atof(optarg); break;
			case 'Z': zalign = 1; pars[zalign].ksize = atoi(optarg); break;
			case 'F': align_mode = atoi(optarg); break;
			case 'v': hzm_debug ++; break;
			case 'M': val = atoi(optarg);
				switch(val){
					case 1: rep_filter = 0; break;
					case 2: out_alns = 1; break;
					case 3: tip_like = 1; break;
					case 4: store_low_cov_edge = 0; break;
					case 5: cnn_filter = 0; break;
					case 6: rep_detach = 1; break;
					case 7: cut_tip = 0; break;
					default: fprintf(stderr, "Unknown mask number -M %d\n", val);
				}
				break;
			default: return usage();
		}
	}
	if(prefix == NULL) return usage();
	if(pbs->size + ngs->size == 0) return usage();
	if(pars[0].ksize > DMO_MAX_KSIZE || pars[0].ksize < 5) return usage();
	if(pars[1].ksize == 0) zalign = 0;
	else zalign = 1;
	if(zalign && (pars[1].ksize > DMO_MAX_KSIZE || pars[1].ksize < 5)) return usage();
	if(regoff == 0) regoff = reglen;
	if(node_cov == 0) node_cov = edge_cov;
	pars[0].max_overhang = 2 * pars[0].xvar;
	pars[1].max_overhang = 2 * pars[1].xvar;
	pars[0].min_aln = pars[0].min_aln * reglen / 1000;
	pars[1].min_aln = pars[1].min_aln * reglen / 1000;
	pars[0].min_mat = pars[0].min_mat * reglen / 1000;
	pars[1].min_mat = pars[1].min_mat * reglen / 1000;
	if(!overwrite && file_exists(prefix)){
		fprintf(hzm_debug_out, "File exists! '%s'\n\n", prefix);
		return usage();
	}
	if(max_idx_bp == 0) max_idx_bp = 0xFFFFFFFFFFFFFFFFLLU;
	if(ncpu <= 0 && _sig_proc_deamon) ncpu = _sig_proc_deamon->ncpu;
	dmo = init_dmo(pars + 0);
	dmo->len_order = len_order;
	fprintf(hzm_debug_out, "[%s] loading reads\n", date());
	tot_bp = 0;
	max_bp *= 1000000;
	nfix = 0;
	if(ngs->size){
		fr = fopen_m_filereader(ngs->size, ngs->buffer);
		seq = NULL;
		while(fread_seq(&seq, fr)){
			if(!hzm_debug && (dmo->n_rd % 10000) == 0){ fprintf(hzm_debug_out, "\r%u", (unsigned)dmo->n_rd); fflush(hzm_debug_out); }
			if(max_bp && tot_bp + (uint32_t)seq->seq.size > max_bp){ free_sequence(seq); break; }
			tag_size = seq->tag.size;
			for(i=0;(int)i<seq->seq.size;i+=WT_MAX_RDLEN){
				len = num_min(seq->seq.size - i, WT_MAX_RDLEN);
				if(i){
					append_string(&seq->tag, "_V", 2);
					add_int_string(&seq->tag, i / WT_MAX_RDLEN);
				}
				push_long_read_dmo(dmo, seq->name.string, seq->name.size, seq->seq.string + i, len);
				if(i){ seq->tag.size = tag_size; seq->tag.string[tag_size] = '\0'; }
				nfix ++;
			}
			tot_bp += seq->seq.size;
		}
		fclose_filereader(fr);
	}
	if(pbs->size){
		if((fr = fopen_m_filereader(pbs->size, pbs->buffer)) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", pbs->buffer[0], __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
		}
		seq = NULL;
		while(fread_seq(&seq, fr)){
			if(!hzm_debug && (dmo->n_rd % 10000) == 0){ fprintf(hzm_debug_out, "\r%u", (unsigned)dmo->n_rd); fflush(hzm_debug_out); }
			if(max_bp && tot_bp + (uint32_t)seq->seq.size > max_bp){ free_sequence(seq); break; }
			if(seq->seq.size > WT_MAX_RDLEN){ seq->seq.size = WT_MAX_RDLEN; seq->seq.string[seq->seq.size] = '\0'; }
			push_long_read_dmo(dmo, seq->name.string, seq->name.size, seq->seq.string, seq->seq.size);
			tot_bp += seq->seq.size;
		}
		fclose_filereader(fr);
	}
	if(!hzm_debug){ fprintf(hzm_debug_out, "\r%u reads", (unsigned)dmo->n_rd); fflush(hzm_debug_out); }
	fprintf(hzm_debug_out, "\n[%s] Done, %u reads, %llu bps\n", date(), (unsigned)dmo->n_rd, tot_bp); fflush(hzm_debug_out);
	ready_dmo(dmo);
	//if(tot_bp > max_idx_bp){
		//pars[0].kmin = 1;
	//}
	//print_proc_stat_info(0);
	g = init_graph(dmo);
	g->node_order = node_order;
	g->reglen = reglen;
	g->regoff = regoff;
	g->regovl = reglen * regovl;
	g->node_merge_cutoff = node_mrg;
	g->min_node_cov = node_cov;
	g->max_node_cov_sg = node_cov;
	g->max_node_cov = max_node_cov;
	g->min_edge_cov = edge_cov;
	g->max_sg_end = max_trace_end;
	g->store_low_cov_edge = store_low_cov_edge;
	g->bub_step = bub_step;
	g->tip_step = tip_step;
	g->rep_step = rep_step;
	g->min_ctg_len = min_ctg_len;
	g->min_ctg_nds = min_ctg_nds;
	g->n_fix = nfix;
	g->only_fix = only_fix;
	g->rep_filter = rep_filter;
	g->rep_detach = rep_detach;
	g->cut_tip = cut_tip;
	if(load){
		fprintf(hzm_debug_out, "[%s] load nodes from %s\n", date(), load);
		fr = fopen_filereader(load);
		load_nodes_graph(g, fr);
		fclose_filereader(fr);
		print_proc_stat_info(0);
		fprintf(hzm_debug_out, "[%s] Done, %llu nodes\n", date(), (unsigned long long)g->nodes->size);
	}
	if(log_rep){
		evtlog = open_file_for_write(prefix, ".events", 1);
	} else evtlog = NULL;
	if(load == NULL){
		FILE *alno;
		alno = out_alns? open_file_for_write(prefix, ".alignments", 1) : NULL;
		fprintf(hzm_debug_out, "[%s] generate nodes, %d threads\n", date(), ncpu);
		if(fast_mode){
			fast_build_nodes_graph(g, align_mode, ncpu, alno);
		} else {
			build_nodes_graph(g, align_mode, max_idx_bp, ncpu, alno);
		}
		if(pars[1].ksize){
			refine_nodes_graph(g, pars + 1, ncpu, evtlog);
		}
		if(alno) fclose(alno);
		fprintf(hzm_debug_out, "[%s] Done, %llu nodes\n", date(), (unsigned long long)g->nodes->size);
	}
	if(load == NULL){
		generic_print_graph(g, print_nodes_graph, prefix, ".1.nodes");
	}
	if(1){
		estimate_genome_size(g, tot_bp, hzm_debug_out);
		////uint32_t mid = estimate_genome_size(g, tot_bp, hzm_debug_out);
		cnt = mask_nodes_by_cov_graph(g, evtlog);
		fprintf(hzm_debug_out, "[%s] masked %llu high coverage nodes (>%d or <%d)\n", date(), (unsigned long long)cnt, max_node_cov, node_cov);
	}
	if(cnn_filter){
		cnt = mask_nodes_by_connectivity_graph(g, ncpu, evtlog);
		fprintf(hzm_debug_out, "[%s] masked %llu repeat-like nodes by local subgraph analysis\n", date(), (unsigned long long)cnt);
	}
	if(tip_like){
		cnt = mask_possible_tip_nodes_graph(g);
		fprintf(hzm_debug_out, "[%s] masked %llu tip-like nodes\n", date(), (unsigned long long)cnt);
	}
	fprintf(hzm_debug_out, "[%s] generate edges\n", date());
	build_edges_graph(g, ncpu);
	fprintf(hzm_debug_out, "[%s] Done, %llu edges\n", date(), (unsigned long long)g->edges->size);
	if(!less_out) generic_print_graph(g, print_reads_graph, prefix, ".1.reads");
	if(!less_out) generic_print_graph(g, print_dot_full_graph,   prefix, ".1.dot");
	fprintf(hzm_debug_out, "[%s] graph clean\n", date()); fflush(hzm_debug_out);
	//cnt = rescue_low_cov_tip_edges_graph(g);
	//fprintf(hzm_debug_out, "[%s] rescued %llu low cov edges\n", date(), (unsigned long long)cnt);
	cnt = cut_binary_edges_graph(g);
	fprintf(hzm_debug_out, "[%s] deleted %llu binary edges\n", date(), (unsigned long long)cnt);
	if(!g->rep_detach){
		cnt = del_isolated_nodes_graph(g);
		fprintf(hzm_debug_out, "[%s] deleted %llu isolated nodes\n", date(), (unsigned long long)cnt);
	}
	cnt = reduce_transitive_edges_graph(g);
	set_init_ends_graph(g);
	fprintf(hzm_debug_out, "[%s] cut %llu transitive edges\n", date(), (unsigned long long)cnt);
	if(!less_out) generic_print_graph(g, print_dot_graph,   prefix, ".2.dot");
	{
		bub = tip = rep = yarn = 0;
		do {
			c = 0;
			do {
				cnt = trim_tips_graph(g, tip_step, bub > 0);
				tip += cnt;
				if(cnt) c = 1;
			} while(cnt);
			do {
				cnt = pop_bubbles_graph(g, bub_step);
				bub += cnt;
				if(cnt) c = 1;
			} while(cnt);
			do {
				cnt = trim_blunt_tips_graph(g);
				tip += cnt;
				if(cnt) c = 1;
			} while(cnt);
			do {
				cnt = pop_bubbles_graph(g, bub_step);
				bub += cnt;
				if(cnt) c = 1;
			} while(cnt);
			do {
				cnt = resolve_yarns_graph(g, bub_step * 5);
				yarn += cnt;
				if(cnt) c = 1;
			} while(cnt);
		} while(c);
		if(bub + tip) {fprintf(hzm_debug_out, "[%s] %llu bubbles; %llu tips; %llu yarns;\n", date(), bub, tip, yarn);} fflush(hzm_debug_out);
	}
	{
		cnt = del_isolated_nodes_graph(g);
		fprintf(hzm_debug_out, "[%s] deleted %llu isolated nodes\n", date(), (unsigned long long)cnt);
	}
	if(!less_out) generic_print_graph(g, print_dot_graph,   prefix, ".3.dot");
	rep = mask_all_branching_nodes_graph(g);
	fprintf(hzm_debug_out, "[%s] cut %llu branching nodes\n", date(), rep);
	if(1){
		cnt = del_isolated_nodes_graph(g);
		fprintf(hzm_debug_out, "[%s] deleted %llu isolated nodes\n", date(), (unsigned long long)cnt);
	}
	fprintf(hzm_debug_out, "[%s] building unitigs\n", date());
	gen_unitigs_graph(g);
	fprintf(hzm_debug_out, "[%s] trimming and extending unitigs by local assembly, %d threads\n", date(), ncpu);
	unitigs2frgs_graph(g, ncpu);
	if(!less_out) generic_print_graph(g, print_frgs_nodes_graph, prefix, ".frg.nodes");
	fprintf(hzm_debug_out, "[%s] generating links\n", date());
	gen_lnks_graph(g, ncpu);
	if(!less_out) generic_print_graph(g, print_frgs_dot_graph, prefix, ".frg.dot");
	if(1){
		cnt = rescue_weak_tip_lnks_graph(g);
		fprintf(hzm_debug_out, "[%s] rescue %llu weak links\n", date(), (unsigned long long)cnt);
	}
	cnt = cut_binary_lnks_graph(g);
	fprintf(hzm_debug_out, "[%s] deleted %llu binary links\n", date(), (unsigned long long)cnt);
	cnt = reduce_transitive_lnks_graph(g);
	fprintf(hzm_debug_out, "[%s] cut %llu transitive links\n", date(), (unsigned long long)cnt);
	//cnt = cut_low_cov_lnks_graph(g, 1);
	//fprintf(hzm_debug_out, "[%s] deleted %llu low cov links\n", date(), (unsigned long long)cnt);
	//if(!less_out) generic_print_graph(g, print_frgs_dot_graph, prefix, ".frg.2.dot");
	cnt = trim_frgtips_graph(g, frgtip_len);
	fprintf(hzm_debug_out, "[%s] cut %llu tips\n", date(), (unsigned long long)cnt);
	//if(!less_out) generic_print_graph(g, print_frgs_dot_graph, prefix, ".frg.3.dot");
	cnt = pop_frg_bubbles_graph(g, bub_step);
	fprintf(hzm_debug_out, "[%s] pop %llu bubbles\n", date(), (unsigned long long)cnt);
	if(!less_out) generic_print_graph(g, print_frgs_dot_graph, prefix, ".ctg.dot");
	fprintf(hzm_debug_out, "[%s] building contigs\n", date());
	cnt = gen_contigs_graph(g, evtlog);
	fprintf(hzm_debug_out, "[%s] searched %llu contigs\n", date(), (unsigned long long)cnt);
	if(1){
		cnt = gen_complex_contigs_graph(g);
		u8i sum;
		sum = 0;
		for(i=g->major_nctg;i<g->ctgs->size;i++){
			sum += cal_offset_traces_graph(g, (tracev*)get_vplist(g->ctgs, i), 0, ((tracev*)get_vplist(g->ctgs, i))->size);
		}
		fprintf(hzm_debug_out, "[%s] added %llu unsolved repetitive contigs, %llu bp\n", date(), (unsigned long long)cnt, sum);
	}
	n50_stat_contigs_graph(g);
	//cnt = generic_print_graph(g, print_isolated_nodes_dot_graph, prefix, ".4.dot");
	//fprintf(hzm_debug_out, "[%s] %llu nodes not in contigs\n", date(), (unsigned long long)cnt);
	//cnt = count_isolated_reads_graph(g);
	//fprintf(hzm_debug_out, "[%s] %llu reads not in contigs\n", date(), (unsigned long long)cnt);
	fprintf(hzm_debug_out, "[%s] outputing contigs\n", date());
	cnt = print_ctgs_graph(g, 0, 0, g->major_nctg, prefix, ".ctg.fa", ".ctg.lay", ncpu, evtlog);
	fprintf(hzm_debug_out, "[%s] outputing reptigs\n", date());
	cnt = print_ctgs_graph(g, cnt, g->major_nctg, g->ctgs->size, prefix, ".rtg.fa", ".rtg.lay", ncpu, evtlog);
	if(evtlog) fclose(evtlog);
	free_cplist(pbs);
	free_cplist(ngs);
	free_dmo(dmo);
	free_graph(g);
	fprintf(hzm_debug_out, "[%s] Program Done\n", date());
	END_STAT_PROC_INFO(stderr);
	return 0;
}
