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

#include "kswx.h"
#include "dbgcns.h"
#include "dagcns.h"
#include "file_reader.h"

static int cns_debug = 0;

thread_beg_def(mcns);
uint32_t eid;
CNS *cns;
u1v *seq1, *seq2;
int reglen;
int W, M, X, I, D, E;
int twice_dbgcns, call_dagcns, candidate_mode;
f4i pM, pX, pI, pD;
kswx_t ret;
u32list *cigars;
int task;
thread_end_def(mcns);

thread_beg_func(mcns);
DAGCNS *dag;
GEGraph *g;
bdpnodev *bnodes;
bdpedgev *bedges;
bdplinkv *linkstack;
u8list *mem_cache;
u1v    *mem_buffer;
u32list *cigars[2];
kswx_t *xs[2];
u4i i, nbeg;
blk_t *blk;
int qb, qe, tb, te;
mem_cache = init_u8list(1024);
cigars[0] = mcns->cigars;
cigars[1] = NULL;
xs[0] = malloc(sizeof(kswx_t));
xs[1] = NULL;
dag = init_dagcns(mcns->W, mcns->M, mcns->X, mcns->I, mcns->D, mcns->E, mcns->pM, mcns->pX, mcns->pI, mcns->pD);
g = init_gegraph();
bnodes = init_bdpnodev(32);
bedges = init_bdpedgev(32);
linkstack = init_bdplinkv(32);
mem_buffer = init_u1v(1024);
thread_beg_loop(mcns);
if(mcns->task == 1){
	ready_cns(mcns->cns);
	run_cns(mcns->cns, mcns->candidate_mode);
	if(mcns->twice_dbgcns && mcns->cns->seq->size){
		if(cns_debug) fprintf(stderr, "DBG1_%d\t%d\t%s\n", mcns->eid, mcns->cns->seq->size, mcns->cns->seq->string);
		run_core_cns(mcns->cns, mcns->cns->cns->buffer, mcns->cns->cns->size);
		if(cns_debug) fprintf(stderr, "DBG2_%d\t%d\t%s\n", mcns->eid, mcns->cns->seq->size, mcns->cns->seq->string);
	} else if(mcns->call_dagcns && mcns->cns->seq->size){
		if(cns_debug) fprintf(stderr, "DBG1_%d\t%d\t%s\n", mcns->eid, mcns->cns->seq->size, mcns->cns->seq->string);
		clear_u8list(dag->cns);
		for(i=0;i<(u4i)mcns->cns->seq->size;i++){
			push_u8list(dag->cns, base_bit_table[(int)mcns->cns->seq->string[i]]);
		}
		gen_pregraph_dagcns(dag);
		for(i=0;i<mcns->cns->qblks->size;i++){
			blk = ref_blkv(mcns->cns->qblks, i);
			nbeg = branched_dynamic_programming_alignment(dag, mcns->cns->qseqs->buffer + blk->off, blk->len, g, bnodes, bedges, mem_buffer);
			if(nbeg == 0){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				continue;
			}
			bdpgraph2dagcns(dag, g, bnodes, bedges, nbeg, linkstack);
		}
		merge_nodes_dagcns(dag);
		gen_consensus_dagcns(dag, NULL);
		if(cns_debug){
			fprintf(stderr, "DAG2_%d\t%d\t", mcns->eid, (int)dag->cns->size);
			print_seq_dagcns(dag, stderr);
			fprintf(stderr, "\n");
		}
		clear_u8list(mcns->cns->cns);
		append_u8list(mcns->cns->cns, dag->cns);
		clear_string(mcns->cns->seq);
		for(i=0;i<dag->cns->size;i++) add_char_string(mcns->cns->seq, bit_base_table[dag->cns->buffer[i]]);
	}
} else if(mcns->task == 2){
	if(mcns->seq2->size == 0){
		mcns->ret = KSWX_NULL;
		continue;
	}
	qb = 0; qe = mcns->seq1->size;
	tb = 0; te = mcns->seq2->size;
	if(qe > mcns->reglen) qb = qe - mcns->reglen;
	if(te > mcns->reglen) te = mcns->reglen;
	kswx_overlap_align_core(xs, cigars, qe - qb, mcns->seq1->buffer + qb, te - tb, mcns->seq2->buffer + tb, 1, mcns->M, mcns->X, mcns->I, mcns->D, mcns->E, mem_cache);
	xs[0]->qb += qb;
	xs[0]->qe += qb;
	xs[0]->tb += tb;
	xs[0]->te += tb;
	mcns->ret = *xs[0];
}
thread_end_loop(mcns);
free(xs[0]);
free_u8list(mem_cache);
free_u1v(mem_buffer);
free_dagcns(dag);
free_gegraph(g);
free_bdpnodev(bnodes);
free_bdpedgev(bedges);
free_bdplinkv(linkstack);
thread_end_func(mcns);

int revise_joint_point(u32list *cigars, int *qe, int *te, int overhang){
	u4i i, op, ln;
	int qq, q, tt, t;
	q = t = 0;
	qq = tt = 0;
	for(i=1;i<=cigars->size;i++){
		op = cigars->buffer[cigars->size - i] & 0xF;
		ln = cigars->buffer[cigars->size - i] >> 4;
		switch(op){
			case 1: q += ln; break;
			case 2: t += ln; break;
			default: qq = q; tt = t; q += ln; t += ln;
		}
		if(q >= overhang && t >= overhang){
			if(cns_debug){
				fprintf(stderr, "qe = %d -> %d\n", *qe, (*qe) - qq);
				fprintf(stderr, "te = %d -> %d\n", *te, (*te) - tt);
				fflush(stderr);
			}
			*qe -= qq;
			*te -= tt;
			return 1;
		}
	}
	return 0;
}

uint32_t run(int reglen, int ksize, int Z, int W, int M, int X, int I, int D, int E, int XX, int OO, int EE, int twice_dbgcns, int call_dagcns, f4i pM, f4i pX, f4i pI, f4i pD, int candidate_mode, uint32_t ncpu, FileReader *fr, FILE *out){
	String *tag, *seq;
	u1v *cseqs;
	u4v *cxs, *cys, *tes, *qes;
	uint32_t i, m, eid, beg, end;
	int c, j, sl, b, e;
	char *ss;
	thread_preprocess(mcns);
	tag = init_string(32);
	seq = init_string(32);
	cseqs = init_u1v(32);
	cxs = init_u4v(32);
	cys = init_u4v(32);
	tes = init_u4v(32);
	qes = init_u4v(32);
	thread_beg_init(mcns, ncpu);
	mcns->eid = 0;
	mcns->cns = init_cns(ksize, Z, W, 0, X, I, D, E);
	mcns->seq1 = init_u1v(32);
	mcns->seq2 = init_u1v(32);
	mcns->reglen = reglen;
	mcns->W = 128;
	mcns->M = M;
	mcns->X = XX;
	mcns->I = OO;
	mcns->D = OO;
	mcns->E = EE;
	mcns->twice_dbgcns = twice_dbgcns;
	mcns->call_dagcns = call_dagcns;
	mcns->pM = pM;
	mcns->pX = pX;
	mcns->pI = pI;
	mcns->pD = pD;
	mcns->candidate_mode = candidate_mode;
	mcns->ret = KSWX_NULL;
	mcns->cigars = init_u32list(16);
	mcns->task = 0;
	thread_end_init(mcns);
	eid = 0;
	thread_wait_one(mcns);
	while(1){
		c = fread_table(fr);
		if(c == -1 || fr->line->string[0] == 'E' || fr->line->string[0] == '>'){
			thread_wake(mcns);
			thread_wait_one(mcns);
			if(mcns->task == 1 && mcns->cns->seq->size){
				if(cns_debug){
					fprintf(stderr, "%s_%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", tag->string, mcns->eid, mcns->cns->qlen, mcns->cns->seq->size, mcns->cns->max_score, mcns->cns->alns[0], mcns->cns->alns[1], mcns->cns->alns[2], mcns->cns->alns[3], mcns->cns->seq->string);
				}
				cxs->buffer[mcns->eid] = cseqs->size;
				for(j=0;j<mcns->cns->seq->size;j++) push_u1v(cseqs, base_bit_table[(int)mcns->cns->seq->string[j]]);
				cys->buffer[mcns->eid] = cseqs->size;
			}
			reset_cns(mcns->cns);
			clear_string(mcns->cns->seq);
			mcns->task = 1;
			mcns->eid = eid ++;
			push_u4v(cxs, 0);
			push_u4v(cys, 0);
			if(fr->line->string[0] == 'E') continue;
			if(tag->size){
				thread_beg_iter(mcns);
				thread_wait(mcns);
				if(mcns->task == 1 && mcns->cns->seq->size){
					if(cns_debug){
						fprintf(stderr, "%s_%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", tag->string, mcns->eid, mcns->cns->qlen, mcns->cns->seq->size, mcns->cns->max_score, mcns->cns->alns[0], mcns->cns->alns[1], mcns->cns->alns[2], mcns->cns->alns[3], mcns->cns->seq->string);
					}
					cxs->buffer[mcns->eid] = cseqs->size;
					for(j=0;j<mcns->cns->seq->size;j++) push_u1v(cseqs, base_bit_table[(int)mcns->cns->seq->string[j]]);
					cys->buffer[mcns->eid] = cseqs->size;
				}
				reset_cns(mcns->cns);
				clear_string(mcns->cns->seq);
				mcns->eid = 0;
				mcns->task = 0;
				thread_end_iter(mcns);
				push_u4v(qes, 0);
				push_u4v(tes, 0);
				for(i=1;i<eid;i++){
					thread_wait_one(mcns);
					if(mcns->task == 2 && mcns->ret.aln > 0){
						if(cns_debug){
							fprintf(stderr, "#%s_%d\t%d\t%d\t%d", tag->string, mcns->eid - 1, (int)mcns->seq1->size, mcns->ret.qb, mcns->ret.qe);
							fprintf(stderr, "\t%s_%d\t%d\t%d\t%d", tag->string, mcns->eid, (int)mcns->seq2->size, mcns->ret.tb, mcns->ret.te);
							fprintf(stderr, "\t%d\t%d\t%d\t%d\t%d\n", mcns->ret.aln, mcns->ret.mat, mcns->ret.mis, mcns->ret.ins, mcns->ret.del);
						}
						{
							b = mcns->ret.qe;
							e = mcns->ret.te;
							if(1){
								revise_joint_point(mcns->cigars, &b, &e, reglen / 5);
							}
							qes->buffer[mcns->eid] = b;
							tes->buffer[mcns->eid] = e;
						}
					}
					mcns->task = 2;
					mcns->eid = i;
					clear_u1v(mcns->seq1); append_array_u1v(mcns->seq1, cseqs->buffer + cxs->buffer[i-1], cys->buffer[i-1] - cxs->buffer[i-1]);
					clear_u1v(mcns->seq2); append_array_u1v(mcns->seq2, cseqs->buffer + cxs->buffer[i], cys->buffer[i] - cxs->buffer[i]);
					mcns->ret = KSWX_NULL;
					push_u4v(qes, mcns->seq1->size);
					push_u4v(tes, 0);
					thread_wake(mcns);
				}
				push_u4v(qes, cys->buffer[eid-1] - cxs->buffer[eid-1]);
				push_u4v(tes, 0);
				thread_beg_iter(mcns);
				thread_wait(mcns);
				if(mcns->task == 2 && mcns->ret.aln > 0){
					if(cns_debug){
						fprintf(stderr, "#%s_%d\t%d\t%d\t%d", tag->string, mcns->eid - 1, (int)mcns->seq1->size, mcns->ret.qb, mcns->ret.qe);
						fprintf(stderr, "\t%s_%d\t%d\t%d\t%d", tag->string, mcns->eid, (int)mcns->seq2->size, mcns->ret.tb, mcns->ret.te);
						fprintf(stderr, "\t%d\t%d\t%d\t%d\t%d\n", mcns->ret.aln, mcns->ret.mat, mcns->ret.mis, mcns->ret.ins, mcns->ret.del);
					}
					{
						b = mcns->ret.qe;
						e = mcns->ret.te;
						if(1){
							revise_joint_point(mcns->cigars, &b, &e, reglen / 5);
						}
						qes->buffer[mcns->eid] = b;
						tes->buffer[mcns->eid] = e;
					}
				}
				mcns->ret = KSWX_NULL;
				mcns->eid = 0;
				mcns->task = 0;
				thread_end_iter(mcns);
				// generate contig seq
				clear_string(seq);
				for(i=0;i<eid;i++){
					beg = cxs->buffer[i] + tes->buffer[i];
					end = cxs->buffer[i] + qes->buffer[i + 1];
					for(m=beg;m<end;m++){
						push_string(seq, bit_base_table[cseqs->buffer[m]]);
					}
					if(cns_debug){
						fprintf(stderr, "=%s_%d\t%d\t%d\n", tag->string, i, seq->size - (end - beg), seq->size);
					}
				}
				fprintf(out, ">%s len=%d\n", tag->string, seq->size);
				for(j=0;j<seq->size;j+=100){
					sl = num_min(j + 100, seq->size);
					char ch = seq->string[sl];
					seq->string[sl] = 0;
					fprintf(out, "%s\n", seq->string + j);
					seq->string[sl] = ch;
				}
			}
			if(c == -1) break;
			clear_string(tag);
			ss = get_col_str(fr, 0) + 1;
			sl = get_col_len(fr, 0) - 1;
			for(j=0;j<sl;j++){
				if(ss[j] == ' ') break;
				push_string(tag, ss[j]);
			}
			clear_string(seq);
			eid = 0;
			clear_u1v(cseqs);
			clear_u4v(cxs);
			clear_u4v(cys);
			clear_u4v(tes);
			clear_u4v(qes);
			thread_wait_one(mcns);
		} else if(fr->line->string[0] == 'S'){
			ss = get_col_str(fr, 5);
			sl = get_col_len(fr, 5);
			add_seq_cns(mcns->cns, ss, sl);
		}
	}
	thread_beg_close(mcns);
	free_cns(mcns->cns);
	free_u1v(mcns->seq1);
	free_u1v(mcns->seq2);
	free_u32list(mcns->cigars);
	thread_end_close(mcns);
	free_u1v(cseqs);
	free_u4v(cxs);
	free_u4v(cys);
	free_u4v(tes);
	free_u4v(qes);
	free_string(tag);
	free_string(seq);
	return 0;
}

int usage(){
	printf(
	"WTDBG-CNS: Consensuser for wtdbg\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtdbg-cns [options]\n"
	"Options:\n"
	" -t <int>    Number of threads, [1]\n"
	" -i <string> Input file(s) *.utg.cns from wtdbg, +, [STDIN]\n"
	" -o <string> Output files, [STDOUT]\n"
	" -f          Force overwrite\n"
	" -j <int>    Expected length of node, or say the overlap length of two adject units in layout file, [1000] bp\n"
	" -k <int>    Kmer size for long reads, [13]\n"
	" -Z <int>    Z-cutoff, drop the lower  (score / <-X>), [4]\n"
	" -W <int>    W-cutoff, drop the lagger (position), [48]\n"
	" -M <int>    Match score, [2]\n"
	" -X <int>    Mismatch score, [-5]\n"
	" -I <int>    Insertion score, [-3]\n"
	" -D <int>    Deletion score, [-4]\n"
	" -E <int>    Gap extension score, [-2]\n"
	" -m <int>    0: DBG correction; 1: DBG plus DBG; 2: DBG plus DAG, [0]\n"
	" -c <int>    Candidate strategy, 0: median length, 1: longest, 2: shortest, 3: first, [0]\n"
	" -v          Verbose\n"
	"\n");
	return 1;
}

int main(int argc, char **argv){
	FileReader *fr;
	cplist *infs;
	FILE *out;
	char *outf;
	int c, ncpu, overwrite, reglen, ksize, Z, W, C, M, X, I, D, E, XX, OO, EE, call_dag, twice_dbg;
	int candidate_mode;
	f4i pM, pX, pI, pD;
	BEG_STAT_PROC_INFO(stderr, argc, argv);
	ncpu = 1;
	reglen = 1000;
	ksize = 13;
	Z = 4;
	W = 48;
	C = 1;
	M = 2;
	X = -5;
	I = -3;
	D = -4;
	E = -2;
	XX = -4;
	OO = -2;
	EE = -1;
	twice_dbg = 0;
	call_dag = 0;
	candidate_mode = 0;
	pM = log(0.85);
	pX = log(0.10);
	pI = log(0.03);
	pD = log(0.02);
	infs = init_cplist(4);
	outf = NULL;
	overwrite = 0;
	while((c = getopt(argc, argv, "hvt:k:i:o:fj:Z:W:C:M:X:I:D:E:m:c:")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'k': ksize = atoi(optarg); break;
			case 'i': push_cplist(infs, optarg); break;
			case 'o': outf = optarg; break;
			case 'f': overwrite = 1; break;
			case 'j': reglen = atoi(optarg); break;
			case 'Z': Z = atoi(optarg); break;
			case 'W': W = atoi(optarg); break;
			case 'C': C = atoi(optarg); break;
			case 'M': M = atoi(optarg); break;
			case 'X': X = atoi(optarg); break;
			case 'I': I = atoi(optarg); break;
			case 'D': D = atoi(optarg); break;
			case 'E': E = atoi(optarg); break;
			case 'm': twice_dbg = atoi(optarg) == 1; call_dag = atoi(optarg) == 2; break;
			case 'c': candidate_mode = atoi(optarg); break;
			case 'v': cns_debug ++; break;
			default: return usage();
		}
	}
	if(outf && !overwrite && file_exists(outf)){
		fprintf(stderr, "File exists! '%s'\n\n", outf);
		return usage();
	}
	if(infs->size) fr = fopen_m_filereader(infs->size, infs->buffer);
	else fr = stdin_filereader();
	if(outf){
		out = open_file_for_write(outf, NULL, 1);
	} else out = stdout;
	run(reglen, ksize, Z, W, M, X, I, D, E, XX, OO, EE, twice_dbg, call_dag, pM, pX, pI, pD, candidate_mode, ncpu, fr, out);
	fclose_filereader(fr);
	if(outf) fclose(out);
	free_cplist(infs);
	END_STAT_PROC_INFO(stderr);
	return 0;
}

