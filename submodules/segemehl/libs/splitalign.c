/*
 *   segemehl - a read aligner
 *   Copyright (C) 2008-2017  Steve Hoffmann and Christian Otto
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */




/*
 *  splitalign.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 12/20/2012 13:25:36 CET
 *  
 */

#include "splitalign.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <float.h>
#include <math.h>

#include "basic-types.h"
#include "brendel.h"
#include "kdchain.h"
#include "mathematics.h"
#include "mapfrag.h"	//(skipped)
#include "sort.h"	//(skipped)
#include "sw.h"	//(skipped)
#include "sufarray.h"	//(skipped)
#include "multicharseq.h"	//(skipped)
#include "biofiles.h"	//(skipped)
#include "stringutils.h"	//(skipped)
#include "memory.h"	//(skipped)
#include "alignment.h"	//(skipped)
#include "locus.h"	//(skipped)




/*-------------------------- bl_initSplitAlignment ---------------------------
 *    
 * @brief initalize the split alignment
 * @author Steve Hoffmann 
 *   
 */
 
splitalignment_t*
bl_initSplitAlignment (splitalignment_t* aln, char* u, uint64_t ulen, uint64_t upos,
    uint64_t vpos, Uint vidx, char rc, uint64_t vdonpos, Uint vdonidx, char donrc, uint64_t vaccpos,
    Uint vaccidx, char accrc)
{
  aln->u =u;
  aln->ulen = ulen;
  aln->upos = upos;
  aln->vpos = vpos;
  aln->vidx = vidx;
  aln->rc = rc;

  aln->vdonpos = vdonpos;
  aln->vdonidx = vdonidx;
  aln->vdonrc = donrc;
  aln->noofsplits = 0;
  aln->splits = NULL;
  aln->vaccpos = vaccpos;
  aln->vaccidx = vaccidx;
  aln->vaccrc = accrc;

  return aln;
}


/*------------------------ bl_destructSplitAlignment -------------------------
 *    
 * @brief destruct the split alignment structure
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_destructSplitAlignment (splitalignment_t *aln)
{
  Uint i;

  for(i=0; i < aln->noofsplits; i++) {
    FREEMEMORY(NULL, aln->splits[i].cigar);
  }

  FREEMEMORY(NULL, aln->splits);
  aln->noofsplits = 0;
  return ;
}


/*-------------------------- bl_addSplitToAlignment --------------------------
 *    
 * @brief add a split to a split alignment with a partial alignment string
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_addSplitToAlignment (splitalignment_t *aln, uint64_t uoff, uint64_t ulen, 
    uint64_t voff, uint64_t vlen, char *cigar)
{

  aln->splits = ALLOCMEMORY(NULL, aln->splits, split_t, aln->noofsplits+1);
  aln->splits[aln->noofsplits].uoff = uoff;
  aln->splits[aln->noofsplits].voff = voff;
  aln->splits[aln->noofsplits].ulen = ulen;
  aln->splits[aln->noofsplits].vlen = vlen;
  aln->splits[aln->noofsplits].cigar = bl_strdup(cigar);

  aln->noofsplits++;

  return ;
}

/*-------------------------- bl_cigarGetSplitAlignment --------------------------
 *    
 * @brief decode a cigar string
 * @author Steve Hoffmann 
 *   
 */

splitalignment_t *
bl_cigarGetSplitAlignment(
    char* u, uint64_t ulen, uint64_t upos,
    uint64_t vpos, Uint idx, char rc, uint64_t vdonpos, Uint vdonidx, char donrc, 
    uint64_t vaccpos, Uint vaccidx, char accrc, char *cigar) {
  
  Uint i, len, allen=0, vallen=0, cur=0; 
  char *buffer, *string = NULL, op;
  splitalignment_t *aln = NULL;
  Uint val = 0;

  uint64_t myupos=upos, myvpos=vpos, lastupos=upos, lastvpos = vpos;

  len = strlen(cigar);
  buffer = calloc(len, sizeof(char));

  aln = ALLOCMEMORY(NULL, NULL, splitalignment_t, 1);
  bl_initSplitAlignment (aln, u, ulen, upos, vpos, idx, rc, vdonpos, 
      vdonidx, donrc, vaccpos, vaccidx, accrc);

  for(i=0; i < len; i++) {
    switch (cigar[i]) {
      case 'S':
        val = atoi(buffer);      
        vallen = snprintf(NULL, 0, "%dS", val);
        string = ALLOCMEMORY(NULL, string, char, allen+vallen+1);
        allen += snprintf(&string[allen], vallen+1, "%dS", atoi(buffer));
        myupos += val;
        memset (buffer, 0, len);
        cur = 0;
        break;
      case '=':
      case 'X':
      case 'M':
        op = cigar[i];
        val = atoi(buffer);    
        vallen = snprintf(NULL, 0, "%d%c", val, op);
        string = ALLOCMEMORY(NULL, string, char, allen+vallen+1);
        allen += snprintf(&string[allen], vallen+1, "%d%c", atoi(buffer), op);
        myupos += val;
        myvpos += val;
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'D':   
        val = atoi(buffer);      
        vallen = snprintf(NULL, 0, "%dD", val);
        string = ALLOCMEMORY(NULL, string, char, allen+vallen+1);
        allen += snprintf(&string[allen], vallen+1, "%dD", atoi(buffer));
        myvpos += val;
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'I':
        val = atoi(buffer);      
        vallen = snprintf(NULL, 0, "%dI", val);
        string = ALLOCMEMORY(NULL, string, char, allen+vallen+1);
        allen += snprintf(&string[allen], vallen+1, "%dI", atoi(buffer));
        myupos += val;
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'N':
        
        bl_addSplitToAlignment (aln, lastupos, myupos-lastupos+1, 
          lastvpos, myvpos-lastvpos+1, string);       
        
        myvpos += atoi(buffer)+1;
        lastvpos = myvpos;
        lastupos = myupos;
        allen = 0;
        FREEMEMORY(NULL, string);
        memset(buffer, 0, len);
        cur = 0;
        break;  //nomoves on skipping
      default :
        buffer[cur++] = cigar[i];
    }
  }
   
  bl_addSplitToAlignment (aln, lastupos, myupos-lastupos+1, 
          lastvpos, myvpos-lastvpos+1, string);
  
  FREEMEMORY(NULL, string);
  free(buffer);
  return aln;
}


/*--------------------------- bl_initSpliceEvents ----------------------------
 *    
 * @brief initalize splice events
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_initSpliceEvents (spliceevents_t *events)
{

  events->noofevents =0;
  events->event = NULL;

  return ;
}

/*------------------------- bl_destructSpliceEvents --------------------------
 *    
 * @brief destruct the splice events
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_destructSpliceEvents (spliceevents_t *events)
{

  FREEMEMORY(NULL, events->event);
  events->noofevents = 0;
  return ;
}

/*---------------------------- bl_addSpliceEvent -----------------------------
 *    
 * @brief add a splice event to the spliceevents_t
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_addSpliceEvent (spliceevents_t *events, uint64_t vdonpos, Uint vdonidx, 
    char vdonrc, uint64_t vaccpos, Uint vaccidx, char vaccrc)
{

  events->event = ALLOCMEMORY(NULL, events->event, spliceevent_t, 
      events->noofevents+1);

  events->event[events->noofevents].vdonpos = vdonpos;
  events->event[events->noofevents].vdonidx = vdonidx;
  events->event[events->noofevents].vdonrc = vdonrc;
  events->event[events->noofevents].donerr = 0;
  events->event[events->noofevents].donins = 0;

  events->event[events->noofevents].vaccpos = vaccpos;
  events->event[events->noofevents].vaccidx = vaccidx;
  events->event[events->noofevents].vaccrc = vaccrc;
  events->event[events->noofevents].accerr = 0;
  events->event[events->noofevents].accins = 0;

  events->noofevents++;

  return ;
}

/*-------------------------- bl_splitAlignmentDump ---------------------------
 *    
 * @brief dump split alignment
 * @author Steve Hoffmann 
 *   
 */

  spliceevents_t*
bl_splitAlignmentGetSpliceEvents (splitalignment_t *aln)
{
  Uint i;
  spliceevents_t *events;
  uint64_t lastvpos=0, lastvlen=0;
  Uint idx = aln->vidx;
  char rc = aln->rc;

  events = ALLOCMEMORY(NULL, NULL, spliceevents_t, 1);
  bl_initSpliceEvents(events);

  //make the right connection to the previous split
  if(aln->vdonpos != -1) { 
    //end to start
    if(!aln->vdonrc && rc) { 
      bl_addSpliceEvent(events, 
          aln->vdonpos, aln->vdonidx, aln->vdonrc,
          aln->splits[0].voff+aln->splits[0].vlen-1, idx, rc);
    //reverse splice on negative strand
    } else if(aln->vdonrc && rc){
      bl_addSpliceEvent(events, 
          aln->splits[0].voff, idx, rc, 
          aln->vdonpos, aln->vdonidx, aln->vdonrc);
    //all other cases: 
    } else { 
      bl_addSpliceEvent(events, 
          aln->vdonpos, aln->vdonidx, aln->vdonrc,
          aln->splits[0].voff, idx, rc);
    }
  }

  for(i=0; i < aln->noofsplits; i++) {
    if(i > 0) {
      //reverse splice on negative strand
      if(rc) {
        bl_addSpliceEvent(events,  
            aln->splits[i].voff, idx, rc,
            lastvpos+lastvlen-1, idx, rc);  
      } else { 
        bl_addSpliceEvent(events, 
            lastvpos+lastvlen-1, idx, rc, 
            aln->splits[i].voff, idx, rc);
      }
    }
    lastvpos = aln->splits[i].voff;
    lastvlen = aln->splits[i].vlen;
  }

  //make the right connection to the upcoming split
  if(aln->vaccpos != -1) { 
    if(rc && !aln->vaccrc) {
      bl_addSpliceEvent(events, lastvpos, idx, rc,
          aln->vaccpos, aln->vaccidx, aln->vaccrc);
    } else if (rc && aln->vaccrc){
      bl_addSpliceEvent(events, 
          aln->vaccpos, aln->vaccidx, aln->vaccrc,
          lastvpos+lastvlen-1, idx, rc);
    } else { 
      bl_addSpliceEvent(events, 
          lastvpos+lastvlen-1, idx, rc,
          aln->vaccpos, aln->vaccidx, aln->vaccrc);
    }
  }


  return events;
}


/*--------------------------- bl_dumpSpliceEvents ----------------------------
 *    
 * @brief print the splice events
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_dumpSpliceEvents (spliceevents_t *events)
{

  Uint i;
  char *arrow = "->\0";

  for(i=0; i < events->noofevents; i++) {
    if(events->event[i].vdonrc &&  events->event[i].vaccrc) arrow = "<-\0";
    else arrow = "->\0";

    fprintf(stderr, "%d: %d-%"PRIu64" (rc:%d) %s %d-%"PRIu64" (rc:%d)\n", i, 
        events->event[i].vdonidx, events->event[i].vdonpos, events->event[i].vdonrc, arrow,
        events->event[i].vaccidx, events->event[i].vaccpos, events->event[i].vaccrc);
  }

  return ;
}




/*------------------------- se_kdAlignEvalSplitAlign -------------------------
 *    
 * @brief post processing of the multisplitalignment
 * @author Steve Hoffmann 
 *   
 */
 
void
se_kdAlignEvalSplitAlign (MultiCharSeq *seq, MultiCharSeqAlignment *a,  
    mappingset_t* set, Uint *totalcover, int *totalscore, 
    unsigned char *trans, int *scores, int indel, unsigned int noofaligns, 
    char ismate, char **seqs, char **quals, segemehl_t *nfo)
{

  unsigned char laststrand=0;
  Uint i, k, edist, ulen, vlen, lastsubidx=0;
  MultiCharSeqAlignment *copy;
  int score; 
  char firstfragment = 1;
  Uint maxedist;

#ifdef DEBUGTRANSALIGN
  Uint sub_start, sub_end, ustart, vstart;
#endif
 
  //TODO: left and right clip!
  //now check if this is query or mate and whether 
  //this has already been
  //mapped. Flaw: only first alignment is checked 
 
  //bl_checkInterFragmentGaps (a, noofaligns, nfo->minfragmentalignlen, nfo->minfragmentalignscore, scores, indel);

  for(k=0, i=0; i < noofaligns; i++) {
 
    Uint previdx = -1, nextidx = -1, prevdist =-1, nextdist =-1;
    if(i > 0) {
      previdx = a[i-1].subidx;
      prevdist = (a[i].refstart > a[i-1].refstart)? a[i].refstart - a[i-1].refstart : a[i-1].refstart - a[i].refstart;
    }

    if(i+1 < noofaligns) {
      nextidx = a[i+1].subidx;
      nextdist = (a[i].refstart > a[i+1].refstart)? a[i].refstart - a[i+1].refstart : a[i+1].refstart - a[i].refstart;
    }

    edist = getEdist(a[i].al);
    ulen = getUalignlen(a[i].al);
    vlen = getValignlen(a[i].al);
    score = getAlignScore(a[i].al, scores, indel);

#ifdef DEBUGTRANSALIGN
    ustart = a[i].al->uoff;
    vstart = a[i].al->voff;

  
    getMultiCharSeqIdxBounds(seq, a[i].subidx, &sub_start, &sub_end);

    DBG("evaluating fragment align: i=%d, query=[%u,%u], uoff=%d, ref=[%u,%u], voff=%d, strand=%d, score=%d, edist=%d, ulen=%d, vlen=%d\n", 
        i, a[i].qrystart, a[i].qrystart+a[i].qrylen-1, ustart, 
        a[i].refstart-sub_start, a[i].refstart-sub_start+a[i].reflen-1, 
        vstart, a[i].strand,  score, edist, ulen, vlen);

    showAlign(a[i].al, DBGDEVICE);
#endif


    maxedist = getMaximumAlignmentEdist(ulen, nfo->accuracy);
    
    if(edist > maxedist) {  
#ifdef DEBUGTRANSALIGN
     DBG("purging fragment align: edist=%d > maxedist=%d\n", edist, maxedist);
#endif
    }
    
    if ((ulen >= nfo->minfragmentalignlen && 
        vlen >= nfo->minfragmentalignlen &&
        score >= nfo->minfragmentalignscore) || //(a[i].pass && ulen >= 10) ||
        (ulen >= 10 && edist <= 1 && ((previdx == -1 && nextidx==a[i].subidx && nextdist < 10000) 
                    ||  (nextidx == -1 && previdx==a[i].subidx && prevdist < 10000)
                    ||  (previdx == a[i].subidx && previdx == a[i].subidx && prevdist < 10000 && nextdist < 10000))
        )) { 

      *totalcover += ulen;
      *totalscore += score;
      k++;
     
   if(firstfragment) { 
     set->elem = ALLOCMEMORY(NULL, set->elem, mapping_t, set->n+1);
     bl_initMapping(&set->elem[set->n], seqs[0], quals[0], 0, 0); 
     set->n++;
     firstfragment = 0;
    }
 
#ifdef DEBUGTRANSALIGN
   DBG("assigned fragment align: i=%d", i);
#endif
      
      copy = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
      copy = bl_copyMCSA(copy, &a[i]);
      bl_addMapFrag(&set->elem[set->n-1], copy, NULL, ismate, 1);
 
      if (k > 1 && (laststrand != a[i].strand || lastsubidx != a[i].subidx)) {
        *trans = 1;
      }

      laststrand = a[i].strand;
      lastsubidx = a[i].subidx;
    } else {
#ifdef DEBUGTRANSALIGN
      DBG("not assigned fragment: i=%d", i);
#endif
      //purge =1;
    } 
  }
  
  return;
}



/*----------------------------- se_AlignSplitMap -----------------------------
 *    
 * @brief align split mapping
 * @author Steve Hoffmann 
 *   
 */

MultiCharSeqAlignment*
se_AlignSplitMap (Uint *mystarts, Uint *myends, uint64_t *mypos, char *myrc, Uint nooffrags, MultiCharSeq *seq, char *querydesc, 
    char **seqs, char **quals, Uint qrylen, int *scores, int indel, int transition)
{

  Uint i, j;
  MultiCharSeqAlignment *a;
  Uint maxedist = qrylen - floor(((double)0.8 * qrylen)/100.);
  unsigned int margin=50, maxmargin=150; //500;
  int ***M, **lmr, **lmv, **lmc;
  PairUint *bestscr;
  char ***K = NULL;

  a = ALLOCMEMORY(space, NULL, MultiCharSeqAlignment, nooffrags);
  Uint *reflens = ALLOCMEMORY(space, NULL, Uint, nooffrags);
  Uint *strands = ALLOCMEMORY(space, NULL, Uint, nooffrags);
  Uint *starts = ALLOCMEMORY(space, NULL, Uint, nooffrags);
  Uint *ends = ALLOCMEMORY(space, NULL, Uint, nooffrags);
  Uint *tstarts = ALLOCMEMORY(space, NULL, Uint, nooffrags);
  Uint *tends = ALLOCMEMORY(space, NULL, Uint, nooffrags);
  Uint *lengths = ALLOCMEMORY(space, NULL, Uint, nooffrags);
  char **refseqs = ALLOCMEMORY(space, NULL, char*, nooffrags);
  Alignment **aligns =  ALLOCMEMORY(space, NULL, Alignment*, nooffrags);
  PairUint *diag = ALLOCMEMORY(space, NULL, PairUint, nooffrags);

  
  for(i=0; i < nooffrags; i++) {

    Uint left = mystarts[i];
    Uint right = myends[i];

    Uint uloff = left; 
    Uint uroff = qrylen - right;
    Uint strand = myrc[i];
    Uint floff;
    Uint flen;
 
    if(i > 0) {
       //distance to previous fragment
       uloff = DISTREV(left, myends[i-1]);
    } 

    if(i < nooffrags-1) {
      //distance to next fragment
       uroff = DISTREV(mystarts[i+1], right);
    }

    if(strand) {
       //increase the length (reverse complement)
       floff = uroff + MIN(maxedist + margin, maxmargin);
       flen = floff + (right - left) + uloff + MIN(maxmargin, maxedist + margin);
    } else { 
       //increase the length 
       floff = maxedist + uloff + margin;
       flen = floff + (right - left) + uroff + MIN(maxmargin, maxedist + margin);
    }

    Uint start = mypos[i];
      
    initMultiCharSeqAlignmentOpt(NULL, &a[i], seq, start,
        querydesc, seqs[strand], quals[strand], left, right,
        qrylen, floff, flen, uloff, uroff, MIN(maxmargin, maxedist+margin), strand);

    a[i].pass = 1;
    strands[i] = strand;
    aligns[i]  = a[i].al;
    refseqs[i] = a[i].refseq;
    reflens[i] = a[i].reflen; 
    lengths[i] = a[i].maxalignlen;
    tstarts[i] = a[i].qrystart;
    tends[i] = tstarts[i]+lengths[i]-1;
    
    if(strands[i]==0) { 
      starts[i] = a[i].qrystart;
      ends[i] = starts[i]+lengths[i]-1;
    } else {
      starts[i] = qrylen - (a[i].qrystart + lengths[i]);
      assert(qrylen >= a[i].qrystart+lengths[i]);
      ends[i] = starts[i]+lengths[i]-1;
      assert(ends[i] <=  qrylen);
    }

    if(strands[i]==0) { 
    diag[i].b = a[i].floff;
    } else {
       if(a[i].reflen >= MIN(maxmargin, maxedist+margin) + uroff + (right - left) - 1){ 
         diag[i].b = a[i].reflen - MIN(maxmargin, maxedist+margin) - uroff - (right - left) + 1;
       } else {
         diag[i].b = 0;
       }
    }

    diag[i].a = left;
  }
  
  M = localmultisplicedmatrixopt(NULL, seqs[0], seqs[1], qrylen, lengths,
      refseqs, reflens, strands, starts, ends, tstarts, tends, nooffrags, indel, transition,
      constscr, scores, &lmv, &lmr, &lmc, &bestscr, &K, diag);

      
  localmultisplicedtracebackopt(NULL, M, seqs[0], seqs[1], qrylen, lengths, 
      refseqs, reflens, strands, starts, ends, tstarts, tends, 
      nooffrags, indel, transition, constscr, scores, 
      aligns, lmv, lmr, lmc, bestscr);

  for(i=0; i < nooffrags; i++) {
    FREEMEMORY(space, lmv[i]);
    FREEMEMORY(space, lmr[i]);
    FREEMEMORY(space, lmc[i]);

    for(j=0; j < lengths[i]+1; j++) {
      FREEMEMORY(space, M[i][j]);
    }

    FREEMEMORY(space, M[i]);
  }

  FREEMEMORY(space, reflens);
  FREEMEMORY(space, refseqs);
  FREEMEMORY(space, strands);
  FREEMEMORY(space, starts);
  FREEMEMORY(space, ends);
  FREEMEMORY(space, tstarts);
  FREEMEMORY(space, tends);
  FREEMEMORY(space, lengths);
  FREEMEMORY(space, diag);
  FREEMEMORY(space, bestscr);

  FREEMEMORY(space, lmv);
  FREEMEMORY(space, lmr);
  FREEMEMORY(space, lmc);
  FREEMEMORY(space, aligns);
  FREEMEMORY(space, M);
 
  return a;
}


/*--------------------------- se_kdAlignSplitChain ---------------------------
 *    
 * @brief align a chain of fragments using local multi spliced alignment
 * @author Steve Hoffmann 
 *   
 */

MultiCharSeqAlignment*
se_kdAlignSplitChain (void *space, branchChain_t *chains, Uint noofchains,
    Suffixarray *arr, MultiCharSeq *seq, char *querydesc, matchstem_t **stems,
    char **seqs, char **quals, Uint qrylen, int *scores, int indel, int transition, 
    unsigned char ismate, Uint *enctab, bitvector *D, Uint *noofaligns, segemehl_t *nfo) {

  Uint k, i, j, q, start, floff = 0, flen =0, 
    maxedist,  
    *strands, *starts, *ends, *tstarts, *tends, *lengths, *reflens,
    sub_start, sub_end;

  unsigned int margin, maxmargin; //50;
  branchChain_t *newchains;
  char **refseqs;
  char ***K = NULL;
  PairUint *bestscr;
  int ***M, **lmr, **lmv, **lmc;
  branchChain_t *cur;
  MultiCharSeqAlignment *a;
  Alignment **aligns; 
  PairUint *diag;
#ifdef DEBUGKBAND
  char ***B;
#endif

  margin = nfo->splitfragmargin;
  maxmargin = nfo->maxsplitfragmargin;
  maxedist = getMaximumAlignmentEdist(qrylen, nfo->accuracy);

  if(noofchains == 0) return NULL;

  qsort(chains, noofchains, sizeof(branchChain_t), cmp_chainscores);
  double maxchainscore = (double) chains[0].score;

  for(k=0; k < noofchains; k++) {
    if(chains[k].score < maxchainscore*nfo->chainscorescale) break;
  }

#ifdef DEBUGTRANSALIGN  
  DBG("call condenseChain: noofchains=%d, minscore=%d\n", 
      k, maxchainscore*nfo->chainscorescale);
  showChains(chains, k, arr, DBGDEVICE, seqs[1], qrylen);
#endif

  newchains = condenseChain(chains, k, seq, arr);
  
#ifdef DEBUGTRANSALIGN   
  DBG("call qsort: noofchains=%d\n", k);
  showChains(newchains, k, arr, DBGDEVICE, seqs[1], qrylen);
#endif

  qsort(newchains, k, sizeof(branchChain_t), cmp_chainlocality);

#ifdef DEBUGTRANSALIGN   
  DBG("after sort: noofchains=%d\n",k);
  showChains(newchains, k, arr, DBGDEVICE, seqs[1], qrylen);
#endif


  q = 0;
  floff = 0;
  flen = 0;

  cur = &newchains[q];

  if(cur->nooffragments <= 1) {  
    wrapChains(space, newchains, k);
    FREEMEMORY(space, newchains);
    return NULL;
  }


#ifdef DEBUBTRANSALIGN
  fprintf(stdout, "nooffrags: %d; scr1:%d scr:2:%d\n", cur->nooffragments, 
      chains[0].score, chains[1].score);
#endif

  a = ALLOCMEMORY(space, NULL, MultiCharSeqAlignment, cur->nooffragments);
  reflens = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  strands = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  starts = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  ends = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  tstarts = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  tends = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  lengths = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  refseqs = ALLOCMEMORY(space, NULL, char*, cur->nooffragments);
  aligns =  ALLOCMEMORY(space, NULL, Alignment*, cur->nooffragments);
  diag = ALLOCMEMORY(space, NULL, PairUint, cur->nooffragments);

  for(i=0; i < cur->nooffragments; i++) {

    Uint uloff = cur->f[i]->start;
    Uint uroff = qrylen - cur->f[i]->end;
 
    if(i > 0) {
       //distance to previous fragment
       uloff = DISTREV(cur->f[i]->start, cur->f[i-1]->end);
    } 

    if(i < cur->nooffragments-1) {
      //distance to next fragment
       uroff = DISTREV(cur->f[i+1]->start, cur->f[i]->end);
    }
    
    if(cur->f[i]->strand) {
       floff = uroff + MIN(maxedist + margin, maxmargin);
       flen = floff + (cur->f[i]->end - cur->f[i]->start) + uloff + MIN(maxmargin, maxedist + margin);
    } else { 
       floff = maxedist + uloff + margin;
       flen = floff + (cur->f[i]->end - cur->f[i]->start) + uroff + MIN(maxmargin, maxedist + margin);
    }

#ifdef DEBUGTRANSALIGN
    fprintf(stdout, "strand:%d floff:%d\tflen:%d\t[%d,%d]->%d\n", 
        cur->f[i]->strand, floff, flen, cur->f[i]->start, 
        cur->f[i]->end, (i < cur->nooffragments-1) ? cur->f[i+1]->start : qrylen);
#endif

    start = cur->f[i]->substart;
    initMultiCharSeqAlignmentOpt(space, &a[i], seq, start,
        querydesc, seqs[cur->f[i]->strand], quals[cur->f[i]->strand], cur->f[i]->start, cur->f[i]->end,
        qrylen, floff, flen, uloff, uroff, MIN(maxmargin, maxedist+margin), cur->f[i]->strand);

    a[i].pass = cur->f[i]->pass;
    aligns[i]  = a[i].al;
    refseqs[i] = a[i].refseq;
    reflens[i] = a[i].reflen; 
    strands[i] = cur->f[i]->strand;
    lengths[i] = a[i].maxalignlen;
    tstarts[i] = a[i].qrystart;
    tends[i] = tstarts[i]+lengths[i]-1;
    
    if(strands[i]==0) { 
      starts[i] = a[i].qrystart;
      ends[i] = starts[i]+lengths[i]-1;
    } else {
      starts[i] = qrylen - (a[i].qrystart + lengths[i]);
      assert(qrylen >= a[i].qrystart+lengths[i]);
      ends[i] = starts[i]+lengths[i]-1;
      assert(ends[i] <=  qrylen);
    }

    if(strands[i]==0) { 
    diag[i].b = a[i].floff;
    } else {
       if(a[i].reflen >= MIN(maxmargin, maxedist+margin) + uroff + (cur->f[i]->end - cur->f[i]->start) - 1){ 
         diag[i].b = a[i].reflen - MIN(maxmargin, maxedist+margin) - uroff - (cur->f[i]->end - cur->f[i]->start) + 1;
       } else {
         diag[i].b = 0;
       }
    }

    
    //diag[i].a = cur->f[i]->start - a[i].qrystart;
    diag[i].a = cur->f[i]->start;


// -DDEBUGMULTISPLICEOPT -DDEBUGTRANSALIGN  
#ifdef DEBUGTRANSALIGN 
    fprintf (stdout, "query sequence of fragment %d [%d->%d]\n", i, starts[i], ends[i]);
   
    Uint h=0;
    for(h=0; h < lengths[i]; h++) {
      if(h && (h%60) == 0) fprintf(stdout, "\n");
      fprintf(stdout, "%c",  a[i].query[starts[i]+h]);
    }
    fprintf(stdout,"\n");
   
    fprintf (stdout, "reference sequence of fragment %d\n", i);
    for(h=0; h < reflens[i]; h++) {
      if(h && (h%60) == 0) fprintf(stdout, "\n");
      fprintf(stdout, "%c",  refseqs[i][h]);
    }
    fprintf(stdout,"\n");

    fprintf(stdout, "%s\n qrylen:%d, fragment:%d, start:%d, strand:%d, curstart:%d, curend:%d, maxedist:%d mapping to [%d,%d]\n", 
            querydesc, qrylen, i, start, strands[i],  starts[i], ends[i], maxedist, a[i].refstart, a[i].refstart+a[i].reflen -1);
    fprintf(stdout, "\n");
#endif
  }

  M = localmultisplicedmatrixopt(space, seqs[0], seqs[1], qrylen, lengths,
      refseqs, reflens, strands, starts, ends, tstarts, tends, cur->nooffragments, indel, transition,
      constscr, scores, &lmv, &lmr, &lmc, &bestscr, &K, diag);


  if(M == NULL) {
    fprintf(stderr, "empty matrix returned for seqs: '%s'/'%s' (%d)\n", 
        seqs[0], seqs[1], qrylen);

    for(i=0; i < cur->nooffragments; i++) {

      getMultiCharSeqIdxBounds(seq, a[i].subidx, &sub_start, &sub_end);
      fprintf(stderr, "fragment %d: %d in %d[%d,%d] '", 
          i, 1 /*arr->suftab[bd[beststart][i]]*/, a[i].subidx, sub_start, sub_end);
      for(j=0; j< qrylen; j++) fprintf(stderr, "%c", refseqs[i][j]);
      fprintf(stderr, "'(%d) strand:%d\n", reflens[i], strands[i]);
    }
    return NULL;
  }
  
#ifdef DEBUGKBAND 
  B = 
#endif 
      localmultisplicedtracebackopt(space, M, seqs[0], seqs[1], qrylen, lengths, 
      refseqs, reflens, strands, starts, ends, tstarts, tends, 
      cur->nooffragments, indel, transition, constscr, scores, 
      aligns, lmv, lmr, lmc, bestscr);

  for(i=0; i < cur->nooffragments; i++) {
    FREEMEMORY(space, lmv[i]);
    FREEMEMORY(space, lmr[i]);
    FREEMEMORY(space, lmc[i]);


#ifdef DEBUGKBAND
    Uint uloff;
    Uint uroff;

    if(i > 0) {
       uloff = (cur->f[i-1]->end < cur->f[i]->start) ? cur->f[i]->start - cur->f[i-1]->end : 0; 
    } else {
       uloff = cur->f[i]->start;
    }

    if(i < cur->nooffragments-1) {
       uroff = (cur->f[i]->end < cur->f[i+1]->start) ? cur->f[i+1]->start - cur->f[i]->end : 0; 
    } else {
       uroff = qrylen-cur->f[i]->end;
    }
   
    fprintf(stderr, "matrix %d of %d\n", i, cur->nooffragments);
    fprintf (stderr, "query sequence of fragment %d (%d,%d)[%d->%d], starts:%u ends:%u fragstart:%d fragend:%d, uloff:%d, uroff:%d, floff:%d, reflen:%d\n", i, cur->f[i]->start, cur->f[i]->end, starts[i], ends[i], diag[i].a, diag[i].b, cur->f[i]->start, cur->f[i]->end, uloff, uroff, a[i].floff, a[i].reflen);

    dumprecursionmatrix2D(stderr, M[i], B[i], K[i], lengths[i], reflens[i], &diag[i]);
#endif

    for(j=0; j < lengths[i]+1; j++) {
#ifdef DEBUGKBAND
      FREEMEMORY(space, B[i][j]);
      FREEMEMORY(space, K[i][j]);
#endif
      FREEMEMORY(space, M[i][j]);
    }
#ifdef DEBUGKBAND
    FREEMEMORY(space, B[i]);
    FREEMEMORY(space, K[i]);
#endif
    FREEMEMORY(space, M[i]);
  }
#ifdef DEBUGKBAND
  FREEMEMORY(space, B);
  FREEMEMORY(space, K);
#endif

  FREEMEMORY(space, M);
  FREEMEMORY(space, bestscr);
  FREEMEMORY(space, diag);

  *noofaligns = cur->nooffragments;

  wrapChains(space, newchains, k);
  FREEMEMORY(space, newchains);

  FREEMEMORY(space, reflens);
  FREEMEMORY(space, refseqs);
  FREEMEMORY(space, strands);
  FREEMEMORY(space, starts);
  FREEMEMORY(space, ends);
  FREEMEMORY(space, tstarts);
  FREEMEMORY(space, tends);
  FREEMEMORY(space, lengths);

 
  FREEMEMORY(space, lmv);
  FREEMEMORY(space, lmr);
  FREEMEMORY(space, lmc);
  FREEMEMORY(space, aligns);

  return a;
}



/*------------------------------ se_kdSplitRead ------------------------------
 *    
 * @brief find the splits of a chimeric reads from matchstem data
 * @author Steve Hoffmann 
 *   
 */

mappingset_t*
bl_splitAlign(void *space, Suffixarray *arr, MultiCharSeq *seq, mappingset_t *set, 
    char *querydesc, matchstem_t **stems, char **seqs, char **quals, Uint len, 
    karlin_t *stats, Uint *enctab, bitvector *D, unsigned char ismate, segemehl_t *nfo) 
{

  unsigned char trans;
  Uint noofchains, totalcover=0;
  int totalscore=0;
  branchfragment_t* fragments;
  branchChain_t *chains;
  MultiCharSeqAlignment *a;
  Uint noofaligns=0, i ;

  chains = branchChain(space, arr, stems, seqs, len, stats, 
      &noofchains, &fragments, nfo->maxfragocc, nfo->maxsplitevalue, 
      nfo->minentropy);

#ifdef DEBUGTRANSALIGN   
  fprintf(stdout, "before kdAlignSplitChain\n");
  showChains(chains, noofchains, arr, stdout, seqs[1], len); 
#endif

  a = se_kdAlignSplitChain (space, chains, noofchains,
      arr, seq, querydesc, stems, seqs, quals, len, nfo->scores, nfo->scores[2], 
      nfo->transition, ismate, enctab, D, &noofaligns, nfo);

  se_kdAlignEvalSplitAlign (seq, a,  set, &totalcover, &totalscore, 
   &trans, nfo->scores, nfo->scores[2], noofaligns, ismate, seqs, quals, nfo);

  for(i=0; i < noofaligns; i++) { 
    wrapMultiCharSeqAlignment(space, &a[i]);
  }

  FREEMEMORY(space, a);
  wrapChains(space, chains, noofchains);
  FREEMEMORY(space, fragments);
  FREEMEMORY(space, chains);

  return set;
}



/*---------------------------- bl_checkSplitAlign ----------------------------
 *    
 * @brief check split align for dangling ends and query gaps
 * @author Steve Hoffmann 
 *   
 */

char
bl_checkSplitAlign (mapping_t *m, char ismate, int *scores, int indel, 
    int maxrgap, int maxlgap, int maxedist5, int maxedist3, int tailsize, int *curscore) {
  Uint j, ulen, lgap=0, rgap=0;
  MultiCharSeqAlignment *mcsa;
  int score=0, fragscore=0;
  int edist3, edist5;
  char checkaligns = 0, trans;
  mapfrag_t *f;


  for(j=0; j < m->n; j++) { 

    f = &m->f[j];

    if(bl_getMapFragIsMate(f) == ismate) { 

      edist3 = 0; 
      edist5 = 0;
      trans = 0;

      mcsa = bl_getMapFragMCSA(f);
      ulen = getUalignlen(mcsa->al);
      fragscore = getAlignScore(mcsa->al, scores, indel);
      score += fragscore;

      //dangling ends with errors
      if(ulen >= tailsize){ 
        edist5 = getSubstringEdist(mcsa->al, 0, tailsize); //15
        edist3 = getSubstringEdist(mcsa->al, ulen-tailsize, ulen);
      } else {
        edist5= getEdist(mcsa->al);
        edist3= edist5;
      }

      if(edist5 > maxedist5) checkaligns |= 1 << 5; //4
      if(edist3 > maxedist3) checkaligns |= 1 << 3; //4

      lgap = f->leftgap;
      rgap = f->rightgap;
      trans = f->nextnoncollinear;

      if(trans) {
        checkaligns |= 1 << 1;
      }

      //first fragment of respective mate
      if(j == 0 || bl_getMapFragIsMate(&m->f[j-1]) != ismate) {
        if((lgap > maxlgap || trans || (ulen < 5 && edist5 > 0))) { //&& !mcsa->strand) { //5
          if(!mcsa->strand) { 
            checkaligns |= 1 << 2;
          } else { 
            checkaligns |= 1 << 4;
          }
        }
      }

      //last fragment of respective mate
      if(j+1 == m->n || bl_getMapFragIsMate(&m->f[j+1]) != ismate) { 
        if((rgap > maxrgap || trans || (ulen < 5 && edist3 > 0) )) { // && !mcsa->strand) { //5
          if(!mcsa->strand) { 
            checkaligns |= 1 << 4;
          } else { 
            checkaligns |= 1 << 2;
          }
        }
      }

    }
  }

  *curscore = score;
  return checkaligns ;
}


/*-------------------------- bl_fixSplitAlignSplit ---------------------------
 *    
 * @brief fix split align with a second round
 * @author Steve Hoffmann 
 *   
 */

mappingset_t*
bl_fixSplitAlignHoffmann (Suffixarray *arr, mappingset_t *set, MultiCharSeq *mseq, char *qrydesc, 
    char **seqs, char **quals, Uint qrylen, char ismate, segemehl_t *nfo)
{
  Uint j, ret;
  int curscore;
  locus_t *mylocus;
  locuslist_t *mylist;
  mappingset_t *tmp; 
  Uint *leftseeds[2];
  Uint leftseedno[] = {0, 0};
  Uint *rightseeds[2];
  Uint rightseedno[] = {0, 0};
  char check;
  Uint maxedist;
  unsigned char trans;
  char checkcnt = 0;
  int transition = -2;
  char leftstrand, rightstrand;
  Uint totalcover=0;
  int totalscore=0;
  const Uint repairFragmentLen = nfo->minfragmentalignlen;

  MultiCharSeqAlignment *a;

  if(repairFragmentLen == 0 || qrylen < repairFragmentLen) return set;

  //fprintf(stdout, "fix split align split of mate:%d\n", ismate);
  //bl_dumpMappingSet(stdout, set);

  //for the first 2 mappings
  for(j=0; j < set->n; j++) { 
    
    check = bl_checkSplitAlign(&set->elem[j], ismate, nfo->scores, nfo->scores[2], 6, 6, 3, 3, 12,  &curscore);
    //fprintf(stdout, "check says %d\n", check);

    if(check) {

      checkcnt++;

      if((check & ( 1 << 2 | 1 << 5 | 1 << 1))) {
        //check gaps left 
        leftseeds[0] = searchSuffixList(NULL, arr, seqs[0], repairFragmentLen, nfo->maxfixsplitocc, &leftseedno[0]);
        leftseeds[1] = searchSuffixList(NULL, arr, seqs[1], repairFragmentLen, nfo->maxfixsplitocc, &leftseedno[1]);
      } 

      if((check & ( 1 << 4 | 1 << 3 | 1 << 1))) {
        //check gaps right
        rightseeds[0] = searchSuffixList(NULL, arr, &seqs[0][qrylen-repairFragmentLen], repairFragmentLen, nfo->maxfixsplitocc, &rightseedno[0]);
        rightseeds[1] = searchSuffixList(NULL, arr, &seqs[1][qrylen-repairFragmentLen], repairFragmentLen, nfo->maxfixsplitocc, &rightseedno[1]);
      }

      mylist = bl_getMappingLocusList(&set->elem[j], mseq, ismate); 

      uint64_t leftd=0, rightd=0;
      leftstrand = 0;
      rightstrand = 0;

      if((check & ( 1 << 2 | 1 << 5))) {
        mylocus = &mylist->loci[0];
        Uint rpos = bl_getLocusChromPosOffset(mylocus, 0);
        Uint strand = bl_getLocusStrand(mylocus);
        leftstrand = strand;

        if(leftseedno[strand]) {  
          ret = binarySearch_left(leftseeds[strand], leftseedno[strand], &rpos, cmp_Uint_bin, NULL);
          if(ret > 0 && leftseeds[strand][ret-1] + repairFragmentLen + 6 <  rpos && (rpos - leftseeds[strand][ret-1] + repairFragmentLen + 1) < nfo->fixsplitdist){
            leftd = leftseeds[strand][ret-1];
          }
        }
      }

      if((check & ( 1 << 4 | 1 << 3))) {
        mylocus = &mylist->loci[mylist->noofloci-1];
        Uint rpos = bl_getLocusChromPosOffset(mylocus, 0);
        Uint strand = bl_getLocusStrand(mylocus);
        Uint rlen = bl_getLocusLen(mylocus);
        rightstrand = strand;
        if(rightseedno[strand]) { 
          ret = binarySearch_right(rightseeds[strand], rightseedno[strand], &rpos, cmp_Uint_bin, NULL);
          Uint rend = rpos + rlen -1;
          if(ret < rightseedno[strand] && rightseeds[strand][ret] > rend+5 &&  (rightseeds[strand][ret] - rpos) < nfo->fixsplitdist) {
            rightd = rightseeds[strand][ret]; 
          }
        }
      }

      if(rightd || leftd) {

        Uint *mystart = bl_locusListGetReadStart(mylist);
        Uint *myends = bl_locusListGetReadEnd(mylist);
        uint64_t *mypos = bl_locusListGetStartPos(mylist);
        char *myrc = bl_locusListStrand(mylist);
        Uint mysize = mylist->noofloci;

        if(leftd) {
          //fprintf(stdout, "leftd : %u\n", leftd);

          mystart = ALLOCMEMORY(space, mystart, Uint, mysize+1);
          myends = ALLOCMEMORY(space, myends, Uint, mysize+1);
          myrc = ALLOCMEMORY(space, myrc, char, mysize+1);
          mypos = ALLOCMEMORY(space, mypos, uint64_t, mysize+1);


         if(leftstrand) {
            mystart[mysize] = qrylen-repairFragmentLen;
            myends[mysize] = qrylen-1;
            mypos[mysize] = leftd;
            myrc[mysize] = leftstrand;
          } else { 
            memmove(&mystart[1], &mystart[0], sizeof(Uint)*mysize);
            memmove(&myends[1], &myends[0], sizeof(Uint)*mysize);
            memmove(&myrc[1], &myrc[0], sizeof(char)*mysize);
            memmove(&mypos[1], &mypos[0], sizeof(uint64_t)*mysize);
            
            mystart[0] = 0;
            myends[0] = repairFragmentLen;
            mypos[0] = leftd;
            myrc[0] = leftstrand;
          }

          mysize++;
        }

        if(rightd) {

          //fprintf(stdout, "rightd : %u\n", rightd);
          mystart = ALLOCMEMORY(space, mystart, Uint, mysize+1);
          myends = ALLOCMEMORY(space, myends, Uint, mysize+1);
          myrc = ALLOCMEMORY(space, myrc, char, mysize+1);
          mypos = ALLOCMEMORY(space, mypos, uint64_t ,mysize+1);
        
          if(rightstrand) {
            memmove(&mystart[1], &mystart[0], sizeof(Uint)*mysize);
            memmove(&myends[1], &myends[0], sizeof(Uint)*mysize);
            memmove(&myrc[1], &myrc[0], sizeof(char)*mysize);     
            memmove(&mypos[1], &mypos[0], sizeof(uint64_t)*mysize);

            mystart[0] = 0;
            myends[0] = repairFragmentLen;
            mypos[0] = rightd;
            myrc[0] = rightstrand;
          } else { 
            
            mystart[mysize] = qrylen-repairFragmentLen;
            myends[mysize] = qrylen-1;
            mypos[mysize] = rightd;
            myrc[mysize] = rightstrand;
          }

          mysize++;
        }
    
      /*  for(Uint i=0; i < mysize; i++) {
            fprintf(stdout, "Hoffmann dumping alignments by AlignSplitMap\n");
            fprintf(stdout, "%d:[%d,%d]->%u(%d)\n", i, mystart[i], myends[i], mypos[i], myrc[i]);
          }
      */
          a= se_AlignSplitMap (mystart, myends, mypos, myrc, mysize,  mseq, qrydesc, seqs, quals, qrylen, 
              nfo->scores, nfo->scores[2], transition);

      /*    for(Uint i=0; i < mysize; i++) {
            fprintf(stdout, "Hoffmann dumping alignments by AlignSplitMap with transition score %d\n", transition);
            showAlign(a[i].al, stdout);
          }
      */
          tmp = ALLOCMEMORY(NULL, NULL, mappingset_t, 1);
          bl_initMappingSet(tmp);

          se_kdAlignEvalSplitAlign (mseq, a,  tmp, &totalcover, &totalscore, 
              &trans, nfo->scores, nfo->scores[2], mysize, ismate, seqs, quals, nfo);
          Uint mindist = -1;

          if(tmp->n ==1) { 
            mindist = bl_getMappingMinFragDist (&tmp->elem[0]);
          }

          //fprintf(stdout, "previous alignment follows\n");
          //bl_dumpMapping(stdout, &set->elem[j]);

          //fprintf(stdout, "new alignment follows oldscore:%d, newscore:%d\n", curscore, totalscore);
          //bl_dumpMappingSet(stdout, tmp);

          maxedist = getMaximumAlignmentEdist(qrylen, nfo->accuracy);
          if(tmp->n == 1 && curscore < totalscore && mindist > maxedist) {
            //remove mate   
            bl_removeMappingQM(&set->elem[j], ismate);
            //fprintf(stdout, "%d elements are added to replace prev align\n", tmp->elem[0].n);
            //add new mate
            for(Uint i=0; i < tmp->elem[0].n; i++) {

              MultiCharSeqAlignment *mcsa = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
              mcsa = bl_copyMCSA(mcsa, tmp->elem[0].f[i].mcsa);
              bl_addMapFrag(&set->elem[j], mcsa, NULL, ismate, 1);
            }
          }

          bl_removeMappingSet(tmp);

          FREEMEMORY(NULL, tmp); 
          FREEMEMORY(NULL , mystart);
          FREEMEMORY(NULL , myends);
          FREEMEMORY(NULL , mypos);
          FREEMEMORY(NULL , myrc);


          for(Uint i=0; i < mysize; i++) { 
            wrapMultiCharSeqAlignment(NULL, &a[i]);
          }

          FREEMEMORY(NULL, a);
      }
          
      wrapLocuslist(mylist);
      FREEMEMORY(space, mylist);

      if(leftseedno[0]) FREEMEMORY(space, leftseeds[0]);
      if(leftseedno[1]) FREEMEMORY(space, leftseeds[1]);
      if(rightseedno[0]) FREEMEMORY(space, rightseeds[0]);
      if(rightseedno[1]) FREEMEMORY(space, rightseeds[1]);

    }
  }

  return set;
}


/*----------------------------- bl_fixSplitAlign -----------------------------
 *    
 * @brief fix dangling ends and query gaps
 * @author Steve Hoffmann 
 *   
 */

void bl_fixSplitAlignBrendel (Suffixarray *arr, mappingset_t *set, MultiCharSeq *mseq, char *qrydesc, 
    char **seqs, char **quals, Uint qrylen, char ismate, segemehl_t *nfo) {

  Uint i, j, ret;
  int curscore, newscore;
  locus_t *mylocus;
  locuslist_t *mylist, *mylist2;
  mapping_t *mymapping; 
  Uint *leftseeds[2];
  Uint leftseedno[] = {0, 0};
  Uint *rightseeds[2];
  Uint rightseedno[] = {0, 0};
  char check;
  int indel = -2; //-3
  int scores[]={1, -2};
  Uint leftmargin, rightmargin;
  //Uint nooffrags;
  char checkcnt = 0;
  const Uint repairFragmentLen = nfo->minfragmentalignlen;

  if(repairFragmentLen == 0 || qrylen < repairFragmentLen) return;

  bl_sortMappingSetByScore(set, scores, indel);

  for(j=0; j < set->n ; j++) { 

    //if(bl_isSplitMappingQM (&set->elem[j], ismate) && checkcnt < 2) { 
    if(checkcnt < 1) { 

      //fprintf(stdout, " fixSplitAlign ----------------- mateno:%d mappingno %d\n", ismate, j);
      //check all fragments in mapping j that belong to ismate
      check = bl_checkSplitAlign(&set->elem[j], ismate, scores, indel, 6, 6, 3, 3, 12,  &curscore);

      if(check) {

        checkcnt++;

        if((check & ( 1 << 2 | 1 << 5 | 1 << 1))) {
          //check gaps left 
          leftseeds[0] = searchSuffixList(NULL, arr, seqs[0], repairFragmentLen, 30000, &leftseedno[0]);
          leftseeds[1] = searchSuffixList(NULL, arr, seqs[1], repairFragmentLen, 30000, &leftseedno[1]);
        } 

        if((check & ( 1 << 4 | 1 << 3 | 1 << 1))) {
          //check gaps right
          rightseeds[0] = searchSuffixList(NULL, arr, &seqs[0][qrylen-repairFragmentLen], repairFragmentLen, 30000, &rightseedno[0]);
          rightseeds[1] = searchSuffixList(NULL, arr, &seqs[1][qrylen-repairFragmentLen], repairFragmentLen, 30000, &rightseedno[1]);
        }


        //extract the split loci
        mylist = bl_getMappingLocusList(&set->elem[j], mseq, ismate); 
        //bl_sortLocuslistByLength(mylist);

        uint64_t leftd=0, rightd=0;

        //myfrag = bl_getMapFragsQM(&set->elem[j], &nooffrags, ismate);
        mylist2 = bl_mergeLocibyDistance (mylist, 10000);
        bl_sortLocuslistByLength(mylist2);
        wrapLocuslist(mylist);
        FREEMEMORY(space, mylist);

        mylist = mylist2;

        //realign the read to all loci  
        for(i=0; i < mylist->noofloci; i++) {
          if(bl_getLocusLen(&mylist->loci[i]) >= 30) { 

            //        fprintf(stdout, "------------------ locusno %d\n", i);

            leftmargin = 100; 
            rightmargin = 100;


            leftd=0, rightd=0;
            mylocus = &mylist->loci[i];
            Uint rpos = bl_getLocusChromPosOffset(mylocus, 0);
            Uint rlen = bl_getLocusLen(mylocus);
            Uint strand = bl_getLocusStrand(mylocus);

            if((check & ( 1 << 2 | 1 << 5))) {
              leftmargin = 1000;
              if(leftseedno[strand]) {  
                //          fprintf(stdout, "left search for closest match to rpos: %u\n", rpos);
                //          fprintf(stdout, "searching in string (%d): %s\n", strand, seqs[strand]);
                ret = binarySearch_left(leftseeds[strand], leftseedno[strand], &rpos, cmp_Uint_bin, NULL);
                //          fprintf(stdout, "returned %d of %d\n", ret, leftseedno[strand]);
                if(ret > 0 && rpos > leftseeds[strand][ret-1] + repairFragmentLen + 1 + 5) {

                  leftd = rpos - leftseeds[strand][ret-1] + repairFragmentLen + 1;

                  //                 fprintf(stdout, "this is %u ; d1:%" PRIu64 "\n", leftseeds[strand][ret-1], leftd);
                }
              }
            }


            if((check & ( 1 << 4 | 1 << 3))) {
              rightmargin = 1000;
              if(rightseedno[strand]) { 
                //          fprintf(stdout, "right search for closest match to rpos: %u\n", rpos);
                //          fprintf(stdout, "searching in string (%d): %s\n", strand, seqs[strand]);
                ret = binarySearch_right(rightseeds[strand], rightseedno[strand], &rpos, cmp_Uint_bin, NULL);
                //          fprintf(stdout, "returned %d of %d\n", ret, rightseedno[strand]);
                if(ret < rightseedno[strand] && rpos + rlen - 1 + 5 < rightseeds[strand][ret] ) {

                  rightd = rightseeds[strand][ret] - rpos + repairFragmentLen + 1;
                  /*            
                                if(ret >0) fprintf(stdout, "prev %" PRId64 " (%" PRId64 ")\n", (int64_t) rightseeds[strand][ret-1], ((int64_t)rpos - rightseeds[strand][ret-1]));
                                fprintf(stdout,"cur  %" PRId64 " (%" PRId64 ")\n", (int64_t) rightseeds[strand][ret], ((int64_t)rpos - rightseeds[strand][ret]));
                                if(ret+1 < rightseedno[strand]) fprintf(stdout, "next %" PRId64 " (%" PRId64 ")\n", (int64_t) rightseeds[strand][ret+1], ((int64_t)rpos - rightseeds[strand][ret+1]));
                                fprintf(stdout, "this is %u ; d1:%" PRIu64 "\n", rightseeds[strand][ret], rightd);
                                */
                }
              }
            }

            if(leftd > 0 && leftd < 15000) leftmargin = leftd; //leftd;
            if(rightd > 0 && rightd < 15000) rightmargin = rightd; //rightd;


            //original pos
            //Uint opos = bl_getLocusChromPosOffset(mylocus, 0); 
            //original length of locus
            //Uint ll = bl_getLocusLenOffset(mseq, mylocus, 0, 0);
            //with left offset
            Uint vpos = bl_getLocusChromPosOffset(mylocus, leftmargin);
            //offset of locus wrt to vpos
            //Uint off = opos - vpos;
            //reference sequence with left offset
            char *genome = bl_getLocusSeqOffset(mseq, mylocus, leftmargin);
            //length of locus with left and right margin
            Uint n = bl_getLocusLenOffset(mseq, mylocus, leftmargin, rightmargin) ;

            //read coords that gave rise to locus
            //Uint a = bl_getMapFragLeft(&myfrag[i]);
            //Uint b = bl_getMapFragRight(&myfrag[i]);
            /*     Uint l = opos;
                   Uint r = opos + ll -1;

                   Uint virtualstart = a;
                   Uint virtualend = b;

            //the following is the case when (1) there is a seed left or right 
            //and (2) [a,b] is not middle
            if(leftd > rightd) {
            virtualstart = ((a - leftd)/2); //lower?
            //virtualend =   ((b ))

            } */

            //fprintf(stdout, "initializing seed %d-[%u,%u]-%d\n", off, a, b, ll-opos+1);
            //        fprintf(stdout, "initializing alignment for %d-[%u,%u]-%d\n", leftmargin, vpos, vpos+n, rightmargin);

            gene_t *model;
            Alignment *myal;
            char *query, *qual;

            if(!strand) {
              //myal = splicedaligndpopt(seqs[0], qrylen, genome, n, &model, l, r, a, b);

              myal = splicedaligndp(seqs[0], qrylen, genome, n, &model);
              query = seqs[0];
              qual = quals[0];
            } else {
              //myal = splicedaligndpopt(seqs[1], qrylen, genome, n, &model, l, r, a, b);

              myal = splicedaligndp(seqs[1], qrylen, genome, n, &model);
              query = seqs[1];
              qual = quals[1];
            }

            //showAlignModel(myal, stdout, model);

            mymapping = bl_dpsplicealign2map(myal, model, mseq, vpos, n, strand, qrydesc,
                query, qual, qrylen, ismate);

            newscore = bl_getMappingScore(mymapping, scores, indel);
            char check = bl_checkSpliceAlign(mymapping);

            //fprintf(stdout, "new score: %d, old score: %d, check: %d\n", newscore, curscore, check);
            //replace mapping frags of query or mate
            if(newscore > curscore && check) { 
              check = 1;

              for(Uint u=0; u < mymapping->n; u++) {
                //   showAlign(new[k], stdout);
                Alignment *new = mymapping->f[u].mcsa->al;
                if((getUalignlen(new) >= nfo->minfragmentalignlen && 
                      getValignlen(new) >= nfo->minfragmentalignlen &&
                      getAlignScore(new, scores, indel) >= nfo->minfragmentalignscore) ||
                    ((getUalignlen(new) >= 6 && getEdist(new) == 0) ||
                     (getUalignlen(new) >= 10 && getEdist(new) <= 1))) {
                } else {
                  check = 0;
                }
              }

              if(check) {   
                //fprintf(stdout, "replacing mapping %d of mate %d with %d fragments\n", j, ismate, mymapping->n);
                //remove old
                bl_removeMappingQM(&set->elem[j], ismate);
                //add new
                for(Uint u=0; u < mymapping->n; u++) {

                  MultiCharSeqAlignment *mcsa = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
                  mcsa = bl_copyMCSA(mcsa, mymapping->f[u].mcsa);
                  bl_addMapFrag(&set->elem[j], mcsa, NULL, ismate, 1);
                }
              }
            } 

            bl_removeMapping(mymapping);
            FREEMEMORY(NULL, mymapping);

            //cleanup
            wrapAlignment(myal);
            FREEMEMORY(NULL, myal);
            bl_wrapGene(model);
            FREEMEMORY(NULL, model);
          }
        }

        wrapLocuslist(mylist);
        FREEMEMORY(space, mylist);

        if(leftseedno[0]) FREEMEMORY(space, leftseeds[0]);
        if(leftseedno[1]) FREEMEMORY(space, leftseeds[1]);
        if(rightseedno[0]) FREEMEMORY(space, rightseeds[0]);
        if(rightseedno[1]) FREEMEMORY(space, rightseeds[1]);

      }
    }
  }

  return ;
  }

char
bl_crosscorrection(Suffixarray *arr, mappingset_t *set, MultiCharSeq *mseq, char *qrydesc, 
    char **seqs, char **quals, Uint qrylen, Uint matelen, Uint checkedist, char ismate, segemehl_t *nfo) {

  Uint j, k, loop, corrected=0;
  char qcheck;
  locuslist_t *mylist, *merged;
  int indel = -3; //-2
  int scores[]={1, -2}, curscore=0;
  Alignment *al, **new;
  Uint mstrand, curstrand;
  Uint edist[] = {0,0}; 
  int score[] = {0,0};
  Uint midx, qidx;
  char correctcnt = 0;

  MultiCharSeqAlignment *mcsa;
  char coll = 1;

  //the refmate is the other one
  Uint refmate = (ismate) ? 0 : 1;


  //extract splice sites for all mappings
  for(j=0; j < set->n; j++) { 

    if(bl_isSplitMappingQM(&set->elem[j], refmate) && 
        bl_isPairedMapping(&set->elem[j]) && correctcnt < 2) {

      loop = 0;
      corrected = 0;
      correctcnt++;

      qcheck = bl_checkSplitAlign(&set->elem[j], ismate, scores, indel, 2, 2, 1, 2, 12,  &curscore);
      mstrand = bl_getMappingStrandQM(&set->elem[j], refmate);
      coll = bl_isCollinearMapping(&set->elem[j], ismate);
      qidx = bl_getQueryStartIdx(&set->elem[j]);
      midx = bl_getMateStartIdx(&set->elem[j]);
      if(qidx < (Uint) -1)
        qidx = bl_getMapFragChrIdx(&set->elem[j].f[qidx]);
      if(midx < (Uint) -1)
        midx = bl_getMapFragChrIdx(&set->elem[j].f[midx]);

      bl_getMappingEdistQM(&set->elem[j], &edist[0], &edist[1]);
      bl_getMappingScoreQM(&set->elem[j], scores, indel, &score[0], &score[1]);

      //        fprintf(stdout, "cross correction mapping %d: isCollinear:%d, qcheck:%d, refissplit:%d, thisedist:%d (threshold:%d), refedist:%d, matelen:%d, refidx:%d thisidx:%d, mstrand:%d\n", 
      //           j, coll, qcheck,  bl_isSplitMappingQM(&set->elem[j], refmate), edist[(Uint)ismate], checkedist, edist[refmate], matelen, midx, qidx, mstrand);

      if(!coll || ((qcheck || edist[(Uint)ismate] > checkedist || qidx != midx) &&
            edist[refmate] < (((double)matelen)*(0.1)))) { 

        //loop the strands start with rc of refmat
        curstrand = (mstrand) ? 0 : 1;
        do {    
          //             fprintf(stdout, "checking %s.\n", seqs[curstrand]);
          loop++;
          //sort and merge all loci (or only the mate locus) wrt the genomic coordiantes
          //only if the reference is collinear. Otherwise, the mapping of the
          //reference mate "decides" the order of fragments to be concatenated

          if(loop == 2) {
            //join the genomic loci of ref mate aligns (i.e. junctions are re-joined)
            mylist = bl_getMappingLocusList(&set->elem[j], mseq, 2); //refmate
            //         bl_showLocusList(stdout, mylist);
            qsort(mylist->loci, mylist->noofloci, sizeof(locus_t), bl_cmpLocusPosNoStrand);
            merged = bl_mergeLocibyDistanceNoStrand (mylist, 0);

            wrapLocuslist(mylist); 
            FREEMEMORY(space, mylist);
            mylist = merged;
          } else {          
            mylist = bl_getMappingLocusList(&set->elem[j], mseq, refmate);
            qsort(mylist->loci, mylist->noofloci, sizeof(locus_t), bl_cmpLocusPosNoStrand);
            //bl_showLocusList(stdout, mylist);
          }

          if(bl_getLocusListCheck(mylist, 250000)) {  

            //add offsets to the left and right site *
            bl_locusListAddOffset(mseq, mylist, qrylen, qrylen);
            //align to locus list *
            al = bl_locusListAlign(mseq, mylist, seqs[curstrand], qrylen, scores, indel);
            //fprintf(stdout, "locuslistalign\n");
            //showAlign(al, stdout);

            //calculated score
            if( (getEdist(al) < edist[(Uint)ismate]) || 
                (getEdist(al) == edist[(Uint)ismate] && (getAlignScore(al, scores, indel)) > score[(Uint)ismate]) 
                || (getEdist(al) == edist[(Uint)ismate] && !coll)
              ) { 


              //expand the alignment *
              new = expandAlignmentDisjoint(al, mylist, mseq);
              corrected  = 1;
             // fprintf(stdout, "after expansion:\n");
              for(k=0; k < mylist->noofloci; k++) {
                if(new[k]) {
               //   showAlign(new[k], stdout);
                  if((getUalignlen(new[k]) >= nfo->minfragmentalignlen && 
                        getValignlen(new[k]) >= nfo->minfragmentalignlen &&
                        getAlignScore(new[k], scores, indel) >= nfo->minfragmentalignscore) ||
                      ((getUalignlen(new[k]) >= 6 && getEdist(new[k]) == 0) ||
                       (getUalignlen(new[k]) >= 10 && getEdist(new[k]) <= 1))) {
                  } else {
                    corrected = 0;
                  }
                }
              }

              if(corrected) { 
                              
                //fprintf(stdout, "cross correction accepted in loop %d\n", loop);
                //remove query part of the mapping *
                bl_removeMappingQM(&set->elem[j], ismate);
                //new mapping is a split mapping if it has more fragments
                char hassplit = (mylist->noofloci > 1) ?  1: 0;
                //iter all loci and add associated alignment pieces where avail
                for(k=0; k < mylist->noofloci; k++) {      
                  if(new[k]) {

                    mcsa = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);

                    initMultiCharSeqAlignment(NULL, mcsa, mseq, 
                        mylist->loci[k].pos, 0, new[k]->vlen, curstrand, 
                        qrydesc, seqs[curstrand], quals[curstrand], qrylen);

                    wrapAlignment(mcsa->al);
                    FREEMEMORY(NULL, mcsa->al);

                    mcsa->al = new[k];
                    bl_addMapFrag(&set->elem[j], mcsa, NULL, ismate, hassplit);
                  }
                }
                //inform about correction
                FREEMEMORY(space, new);
              } else {
                for(k=0; k < mylist->noofloci; k++) {      
                  if(new[k]) {
                    wrapAlignment(new[k]);
                  }
                }
                FREEMEMORY(space, new);
              }
            }

            FREEMEMORY(NULL, al->v);
            wrapAlignment(al);
            FREEMEMORY(NULL, al);
          }

          wrapLocuslist(mylist);
          FREEMEMORY(NULL, mylist);

        } while (!corrected && loop < 2);
      }
    }
  }

  return corrected;
}

