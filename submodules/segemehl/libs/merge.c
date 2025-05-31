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
 * merge.c
 * functions to merge matches
 *
 */

// #ifndef _DEFAULT_SOURCE
//   #define _DEFAULT_SOURCE
// #endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "debug.h"
#include "info.h"
#include "basic-types.h"
#include "stringutils.h"
#include "mathematics.h"
#include "biofiles.h"
#include "fileBins.h"
#include "matchfilesfields.h"
#include "merge.h"
#include "filebuffer.h"
#include "sort.h"
#include "segemehl.h"
#include "bamio.h"

/*------------------------- bl_mergefilesInit ---------------------------------
 *
 * @brief init container for multiple merge files
 * @author Christian Otto
 *
 */
void bl_mergefilesInit(bl_mergefiles_t *files, Uint nooffiles){
  files->f = ALLOCMEMORY(NULL, NULL, bl_mergefile_t, nooffiles);
  files->nooffiles = nooffiles;
}


/*----------------------- bl_mergefilesDestruct -------------------------------
 *
 * @brief destruct container for multiple merge files
 * @author Christian Otto
 *
 */
void bl_mergefilesDestruct(bl_mergefiles_t *files){
  if (files->f != NULL){
    FREEMEMORY(NULL, files->f);
    files->f = NULL;
  }
  files->nooffiles = 0;
}

/*--------------------------- bl_mergefileInit --------------------------------
 *
 * @brief init merge file container
 * @author Christian Otto
 *
 */
void bl_mergefileInit(bl_mergefile_t *file, FILE *fp){
  file->fp = fp;

#ifdef FILEBUFFEREDMERGE
  bl_circBufferInit(&file->cb, 10000000, fp, NULL);
  int n = fread(file->cb.buffer, sizeof(char), file->cb.size, file->cb.dev);
  file->cb.end = n-1;
#endif
  file->eof = 0;
  file->buffer = NULL;
  file->complete = 0;
  file->entry = ALLOCMEMORY(NULL, NULL, bl_mergefilematch_t, 1);
  bl_mergefilematchInit(file->entry);
  file->mtx = ALLOCMEMORY(space, NULL, pthread_mutex_t, 1);
  pthread_mutex_init(file->mtx, NULL);
}

/*------------------------ bl_mergefileDestruct -------------------------------
 *
 * @brief destruct merge file container
 * @author Christian Otto
 *
 */
void bl_mergefileDestruct(bl_mergefile_t *file){
  file->fp = NULL;
  file->eof = 0;
  file->complete = 0;
  if (file->buffer != NULL){
    FREEMEMORY(NULL, file->buffer);
    file->buffer = NULL;
  }
  if (file->entry != NULL){
    bl_mergefilematchDestruct(file->entry);
    FREEMEMORY(NULL, file->entry);
    file->entry = NULL;
  }

#ifdef FILEBUFFEREDMERGE
  bl_circBufferDestruct(&file->cb);
#endif
    FREEMEMORY(space, file->mtx);
}

/*------------------------ bl_mergefilematchInit ------------------------------
 *
 * @brief init merge file entry
 * @author Christian Otto
 *
 */
void bl_mergefilematchInit(bl_mergefilematch_t *entry){
  entry->qname = NULL;
  entry->matchid = -1;
  entry->read = NULL;
  entry->mate = NULL;
}


/*---------------------- bl_mergefilematchDestruct ----------------------------
 *
 * @brief destruct merge file entry
 * @author Christian Otto
 *
 */
void bl_mergefilematchDestruct(bl_mergefilematch_t *entry){
  entry->qname = NULL;
  entry->matchid = -1;
  if (entry->read != NULL){
    bl_samDestruct(entry->read);
    FREEMEMORY(NULL, entry->read);
    entry->read = NULL;
  }
  if (entry->mate != NULL){
    bl_samDestruct(entry->mate);
    FREEMEMORY(NULL, entry->mate);
    entry->mate = NULL;
  }
}


/*--------------------- bl_mergefilematchCompareUnmapSeed ---------------------
 *
 * @brief compare two sam record entries regarding best found seed
 *        regarding maximum sum of seed length (stored in SAM tag ZL). If
 *        there is no unmapped alignment or ZL is not found, the first entry
 *        is deemed better. Returns -1 if first is better than or equal to
 *        second entry, 1 otherwise.
 *        Christian Otto
 *
 */
int bl_mergeCompareUnmapped(samrec_t *i, samrec_t *j){
  int tmpi, tmpj;
  samtag_t *tag = NULL;

  if (i == NULL || j == NULL){
    return -1;
  }

  tag = bl_samgetTag(i, "ZL");
  if (!tag){
    return -1;
  }
  tmpi = atoi(tag->val);

  tag = bl_samgetTag(j, "ZL");
  if (!tag){
    return -1;
  }
  tmpj = atoi(tag->val);

  if (tmpi >= tmpj){
    return -1;
  }
  else {
    return 1;
  }
}


/*------------------ bl_mergefilematchRemoveRedundantUnmapped -----------------
 *
 * @brief remove redundant unmapped alignment records in mergefilematch list.
 *        These are reads/mates for which either aligned records exist or
 *        "better" unmapped records (regarding the found seed length).
 *        Christian Otto
 *
 */
void bl_mergefilematchRemoveRedundantUnmapped(bl_mergefilematch_t **list, Uint n){
  int i, cmp, tmpiread, tmpimate, bestread, bestmate, tmpbestread, tmpbestmate;

  bestread = -1, bestmate = -1;
  for (i = 0; i < n; i++){
    if (list[i]->read != NULL){
      /* compare against previous best read */
      if (bestread != -1){
        tmpbestread = (list[bestread]->read->flag & 0x4);
        tmpiread = (list[i]->read->flag & 0x4);
        /* previous was aligned */
        if (!tmpbestread){
          /* remove current if unmapped */
          if (tmpiread){
            bl_samDestruct(list[i]->read);
            FREEMEMORY(NULL, list[i]->read);
            list[i]->read = NULL;
          }
        }
        else {
          /* remove best if current is aligned */
          if (!tmpiread){
            bl_samDestruct(list[bestread]->read);
            FREEMEMORY(NULL, list[bestread]->read);
            list[bestread]->read = NULL;
            bestread = i;
          }
          /* otherwise: take best unmapped */
          else {
            cmp = bl_mergeCompareUnmapped(list[bestread]->read, list[i]->read);
            if (cmp < 0){
              bl_samDestruct(list[i]->read);
              FREEMEMORY(NULL, list[i]->read);
              list[i]->read = NULL;
            }
            else {
              bl_samDestruct(list[bestread]->read);
              FREEMEMORY(NULL, list[bestread]->read);
              list[bestread]->read = NULL;
              bestread = i;
            }
          }
        }
      }
      else {
        bestread = i;
      }
    }
    if (list[i]->mate != NULL){
      /* compare against previous best mate */
      if (bestmate != -1){
        tmpbestmate = (list[bestmate]->mate->flag & 0x4);
        tmpimate = (list[i]->mate->flag & 0x4);
        /* previous was aligned */
        if (!tmpbestmate){
          /* remove current if unmapped */
          if (tmpimate){
            bl_samDestruct(list[i]->mate);
            FREEMEMORY(NULL, list[i]->mate);
            list[i]->mate = NULL;
          }
        }
        else {
          /* remove best if current is aligned */
          if (!tmpimate){
            bl_samDestruct(list[bestmate]->mate);
            FREEMEMORY(NULL, list[bestmate]->mate);
            list[bestread]->mate = NULL;
            bestmate = i;
          }
          /* otherwise: take best unmapped */
          else {
            cmp = bl_mergeCompareUnmapped(list[bestmate]->mate, list[i]->mate);
            if (cmp < 0){
              bl_samDestruct(list[i]->mate);
              FREEMEMORY(NULL, list[i]->mate);
              list[i]->mate = NULL;
            }
            else {
              bl_samDestruct(list[bestmate]->mate);
              FREEMEMORY(NULL, list[bestmate]->mate);
              list[bestmate]->mate = NULL;
              bestmate = i;
            }
          }
        }
      }
      else {
        bestmate = i;
      }
    }
  }
  return;
}
/*------------------- bl_mergefilematchRemoveNextInformation ------------------
 *
 * @brief remove/unset next information in terms of RNEXT, PNEXT, and TLEN.
 *        Christian Otto
 *
 */
void bl_mergefilematchRemoveNextInformation(bl_mergefilematch_t **list, Uint n){
  int i;
  for (i = 0; i < n; i++){
    if (list[i]->read != NULL){
      FREEMEMORY(NULL, list[i]->read->rnext);
      list[i]->read->rnext = bl_strdup("*\0");
      list[i]->read->pnext = 0;
      list[i]->read->tlen = 0;
    }
    if (list[i]->mate != NULL){
      FREEMEMORY(NULL, list[i]->mate->rnext);
      list[i]->mate->rnext = bl_strdup("*\0");
      list[i]->mate->pnext = 0;
      list[i]->mate->tlen = 0;
    }
  }
  return;
}
/*------------------- bl_mergefilematchRemoveSuboptimalPairs ------------------
 *
 * @brief remove suboptimal paired alignments in a mergefilematch list which
 *        are alignments that exceeded the minimal pair edit distance.
 *        Christian Otto
 *
 */
void bl_mergefilematchRemoveSuboptimalPairs(bl_mergefilematch_t **list, Uint n){
  int i, edist, bestedist = -1;

  for (i = 0; i < n; i++){
    if (list[i]->read != NULL && list[i]->mate != NULL){
      edist = atoi(bl_samgetTag(list[i]->read, "NM")->val);
      edist += atoi(bl_samgetTag(list[i]->mate, "NM")->val);
      if (bestedist == -1 || edist < bestedist){
        bestedist = edist;
      }
    }
  }

  for (i = 0; i < n; i++){
    if (list[i]->read != NULL && list[i]->mate != NULL){
      edist = atoi(bl_samgetTag(list[i]->read, "NM")->val);
      edist += atoi(bl_samgetTag(list[i]->mate, "NM")->val);
      if (edist > bestedist) {
        bl_samDestruct(list[i]->read);
        FREEMEMORY(NULL, list[i]->read);
        list[i]->read = NULL;
        bl_samDestruct(list[i]->mate);
        FREEMEMORY(NULL, list[i]->mate);
        list[i]->mate = NULL;
      }
    }
  }
  return;
}
/*----------------- bl_mergefilematchRemoveSuboptimalSingletons ---------------
 *
 * @brief remove suboptimal singleton alignments in a mergefilematch list which
 *        are alignments that exceeded the minimal read/mate edit distance.
 *        Note that the minimal edit distance may be different for read and mate
 *        alignments.
 *        Christian Otto
 *
 */
void bl_mergefilematchRemoveSuboptimalSingletons(bl_mergefilematch_t **list, Uint n){
  int i, edist, bestreadedist = -1, bestmateedist = -1;

  for (i = 0; i < n; i++){
    if (list[i]->read != NULL && !(list[i]->read->flag & 0x4)){
      edist = atoi(bl_samgetTag(list[i]->read, "NM")->val);
      if (bestreadedist == -1 || edist < bestreadedist){
        bestreadedist = edist;
      }
    }
    if (list[i]->mate != NULL && !(list[i]->mate->flag & 0x4)){
      edist = atoi(bl_samgetTag(list[i]->mate, "NM")->val);
      if (bestmateedist == -1 || edist < bestmateedist){
        bestmateedist = edist;
      }
    }
  }
  for (i = 0; i < n; i++){
    if (list[i]->read != NULL && !(list[i]->read->flag & 0x4)){
      edist = atoi(bl_samgetTag(list[i]->read, "NM")->val);
      if (edist > bestreadedist){
        bl_samDestruct(list[i]->read);
        FREEMEMORY(NULL, list[i]->read);
        list[i]->read = NULL;
      }
    }
    if (list[i]->mate != NULL && !(list[i]->mate->flag & 0x4)){
      edist = atoi(bl_samgetTag(list[i]->mate, "NM")->val);
      if (edist > bestmateedist){
        bl_samDestruct(list[i]->mate);
        FREEMEMORY(NULL, list[i]->mate);
        list[i]->mate = NULL;
      }
    }
  }
  return;
}
/*--------------------- bl_mergefileComparePairingState -----------------------
 *
 * @brief compare two merge file entries regarding SAM flag in case of
 *        paired-end reads (i.e. best pairing state), returns -1 if
 *        first is better, 1 if second is better, 0 otherwise
 *        @author Christian Otto
 *
 */
int bl_mergefilematchComparePairingState(bl_mergefilematch_t *i, bl_mergefilematch_t *j,
    Uint *pairingstate){
  int tmpi, tmpj, tmpiread, tmpimate, tmpjread, tmpjmate;

  tmpiread = (i->read != NULL) && !(i->read->flag & 0x4);
  tmpimate = (i->mate != NULL) && !(i->mate->flag & 0x4);
  tmpjread = (j->read != NULL) && !(j->read->flag & 0x4);
  tmpjmate = (j->mate != NULL) && !(j->mate->flag & 0x4);

  tmpi = tmpiread & tmpimate;
  tmpj = tmpjread & tmpjmate;

  /* check if either one is fully mapped and properly paired */
  if (!(tmpi | tmpj)){
    *pairingstate = NOT_PAIRED;
  }
  else {
    *pairingstate = PAIRED;
    if (tmpi != tmpj){
      return tmpj - tmpi;
    }
    else {
      assert((i->read->flag & 0x2) == (i->mate->flag & 0x2));
      assert((j->read->flag & 0x2) == (j->mate->flag & 0x2));
      tmpi = (i->read->flag & 0x2);
      tmpj = (j->read->flag & 0x2);
      if (tmpi != tmpj){
        return tmpj - tmpi;
      }
    }
  }
  return 0;
}
/*------------------------- bl_mergefileFastaIDCompare -------------------------------
 *
 * @brief compare two fasta descriptions if they contain the same fasta ID,
 *        in case of paired-end data, it disregards /1 or /2 at the end or
 *        any  differences after the first white space
 *        returns 1 if both descriptions contain the same ID,  0 otherwise
 * @author Christian Otto
 *
 */
unsigned char
bl_mergefileFastaIDCompare(char *desc1, Uint desc1len, char *desc2, Uint desc2len) {

  char *id, *id2, *tok1, *tok2;
  unsigned char res;

  id = ALLOCMEMORY(NULL, NULL, char, desc1len+2);
  id2 = ALLOCMEMORY(NULL, NULL, char, desc2len+2);

  strcpy(id, desc1);
  strcpy(id2, desc2);

  tok1 = strtok(id, "/");
  tok2 = strtok(id2, "/");
  res = (strcmp(tok1, tok2)==0);

  if(!res) {
    FREEMEMORY(NULL, id);
    FREEMEMORY(NULL, id2);

    id = ALLOCMEMORY(NULL, NULL, char, desc1len+2);
    id2 = ALLOCMEMORY(NULL, NULL, char, desc2len+2);

    strcpy(id, desc1);
    strcpy(id2, desc2);

    tok1 = strtok(id, " ");
    tok2 = strtok(id2, " ");
    res = (strcmp(tok1, tok2)==0);
  }

  FREEMEMORY(NULL, id);
  FREEMEMORY(NULL, id2);
  return res;
}
/*------------------------- bl_mergeParseLine ---------------------------------
 *
 * @brief parses a SAM-formatted line (single or paired-end) and
 *        inserts the data in the given container
 *        NOTE: split reads not supported up to now
 * @author Christian Otto
 *
 */
unsigned char bl_mergeParseLine(samheader_t *head, bl_mergefilematch_t *entry, char *line, Uint *len){
  samrec_t *rec = NULL;
  samtag_t *tag = NULL;
  int matchid = -1;

  /* parse record */
  rec = bl_samline2rec(line, *len, head);

  /* get matchid (in case of mapped) */
  tag = bl_samgetTag(rec, "HI");
  if (tag == NULL || strcmp(tag->type, "i") || atoi(tag->val) < 0){
    DBG("Error in reading HI tag for SAM entry: %sExit forced.\n",
        bl_samprintSamrec2Buffer(rec, '\n'));
    exit(-1);
  }
  matchid = atoi(tag->val);

  if (rec->flag & 0x800){
    DBG("Split reads not supported yet. Exit forced.\n", NULL);
    exit(-1);
  }

  /* abort if non-valid flags if paired (either first or second in pair) */
  if (rec->flag & 0x1){
    if (!((rec->flag & 0x40) ^
          (rec->flag & 0x80))){
      DBG("Invalid SAM flag for entry: %sExit forced.\n",
          bl_samprintSamrec2Buffer(rec, '\n'));
      exit(-1);
    }
  } else {
    if ((rec->flag & 0x2) ||
        (rec->flag & 0x40) ||
        (rec->flag & 0x80)){
      DBG("Invalid SAM flag for entry: %sExit forced.\n",
          bl_samprintSamrec2Buffer(rec, '\n'));
      exit(-1);
    }
  }

  if (entry->qname == NULL){
    entry->qname = rec->qname;
    entry->matchid = matchid;
  } else {
    if (! bl_mergefileFastaIDCompare(entry->qname, strlen(entry->qname),
          rec->qname, strlen(rec->qname)) ||
        entry->matchid != matchid){
      bl_samDestruct(rec);
      FREEMEMORY(NULL, rec);
      return 1;
    }
  }

  /* assign as read if unpaired or mate 1 */
  if (!(rec->flag & 0x1) || (rec->flag & 0x40)){
    if (entry->read != NULL){
      DBG("Multiple alignments for read %s with same HI tag value found. Exit forced.\n",
          entry->qname);
      exit(-1);
    }
    entry->read = rec;
  } else {
    /* assign as mate otherwise */
    if (entry->mate != NULL){
      DBG("Multiple alignments for read %s with same HI tag value found. Exit forced.\n",
          entry->qname);
      exit(-1);
    }
    entry->mate = rec;
  }
  /* set input line as processed via len */
  *len = 0;

  if (rec->flag & 0x8){
    return 1;
  } else {
    return 0;
  }
}


/*-------------------------- bl_mergeReadNext ---------------------------------
 *
 * @brief read next entry in mergefile
 * @author Christian Otto
 *
 */
void bl_mergeReadNext(samheader_t *head, bl_mergefile_t *file){
  Uint len = 0;
  char *buffer = NULL;

  if (file->buffer != NULL){
    //a string for the next iteration was stored
    len = strlen(file->buffer);
    file->complete = bl_mergeParseLine(head, file->entry, file->buffer, &len);
    assert(len == 0);
    FREEMEMORY(NULL, file->buffer);
    file->buffer = NULL;
  }

  if (!file->complete && !file->eof){

#ifndef FILEBUFFEREDMERGE
    char ch;
    Uint buffersize = 1024;
    buffer = ALLOCMEMORY(NULL, NULL, char, buffersize);
    len = 0;

    while((ch = getc(file->fp)) != EOF) {
      /* extend buffer */
      if(len == buffersize-1) {
        buffersize = 2*buffersize+1;
        buffer = ALLOCMEMORY(NULL, buffer, char, buffersize);
      }
      /* process buffer */
      if(ch == '\n' && len > 0) {
        buffer[len] = '\0';

        file->complete = bl_mergeParseLine(head, file->entry, buffer, &len);

        if (file->complete){
          if (len > 0){
            file->buffer = buffer;
          }
          break;
        } else {
          len = 0;
          continue;
        }
      } else {
        if(ch != '\n') buffer[len++] = ch;
      }
    }
#else
    while((buffer = bl_circBufferReadLine(&file->cb, &len))) {
      assert(len == strlen(buffer));
      file->complete = bl_mergeParseLine(head, file->entry, buffer, &len);

      if (file->complete){
        if (len > 0){
          //the next sam for next iteration
          file->buffer = buffer;
        }
        break;
      } else {
        len = 0;
        FREEMEMORY(NULL, buffer);
        continue;
      }
    }

#endif
#ifndef FILEBUFFEREDMERGE
    /* set end of file */
    if (!file->eof && ch == EOF){
      file->eof = 1;
      if (file->entry->qname != NULL){
        file->complete = 1;
      }
    }

    if (file->buffer == NULL){
      FREEMEMORY(NULL, buffer);
    }
#else

    if (!file->eof && file->cb.feof){
      file->eof = 1;
      if (file->entry->qname != NULL){
        file->complete = 1;
      }
    }

    if (file->buffer == NULL){
      FREEMEMORY(NULL, buffer);
    }
#endif
  }
}


/*------------------------- bl_mergeUpdateTag ---------------------------------
 *
 * @brief update SAM tags NH and HI for given record
 * @author Christian Otto
 *
 */
void bl_mergeUpdateTag(samrec_t *rec, Uint matchid, Uint noofmatches){
  samtag_t *tag;

  /* get and update HI tag (if necessary) */
  tag = bl_samgetTag(rec, "HI");
  if (tag == NULL || strcmp(tag->type, "i")){
    DBG("HI tag is missing or invalid in SAM entry: %s",
        bl_samprintSamrec2Buffer(rec, '\n'));
    exit(-1);
  }
  FREEMEMORY(NULL, tag->tag);
  FREEMEMORY(NULL, tag->val);
  bl_asprintf(&tag->tag, "HI:i:%d", matchid);
  bl_asprintf(&tag->val, "%d", matchid);

  /* get and update NH tag (if necessary) */
  tag = bl_samgetTag(rec, "NH");
  if (tag == NULL || strcmp(tag->type, "i")){
    DBG("NH tag is missing or invalid in SAM entry: %s",
        bl_samprintSamrec2Buffer(rec, '\n'));
    exit(-1);
  }
  FREEMEMORY(NULL, tag->tag);
  FREEMEMORY(NULL, tag->val);
  bl_asprintf(&tag->tag, "NH:i:%d", noofmatches);
  bl_asprintf(&tag->val, "%d", noofmatches);
}
/*-------------------- bl_mergeCountAlignedMappings ---------------------------
 *
 * @brief count aligned mappings (i.e. excluding unmapped) for query and mate
 *         in list of matches
 * @author Christian Otto
 *
 */
void bl_mergeCountAlignedMappings(bl_mergefilematch_t **list, int n,
    Uint *nqueries, Uint *nmates){
  Uint i, noofmates = 0, noofqueries = 0;
  for (i=0; i < n; i++){
    if (list[i]->read != NULL && !(list[i]->read->flag & 0x4)){
      noofqueries++;
    }
    if (list[i]->mate != NULL && !(list[i]->mate->flag & 0x4)){
      noofmates++;
    }
  }
  *nqueries = noofqueries;
  *nmates = noofmates;

  return;
}
/*----------------------- bl_mergeCountMappings -------------------------------
 *
 * @brief count mappings for query and mate in list of matches
 * @author Christian Otto
 *
 */
void bl_mergeCountMappings(bl_mergefilematch_t **list, int n,
    Uint *nqueries, Uint *nmates){
  Uint i, noofmates = 0, noofqueries = 0;
  for (i=0; i < n; i++){
    if (list[i]->read != NULL){
      noofqueries++;
    }
    if (list[i]->mate != NULL){
      noofmates++;
    }
  }
  *nqueries = noofqueries;
  *nmates = noofmates;

  return;
}
/*-------------------------- cmp_mergemap_qsort  -----------------------------
 *
 * @brief qsort function to sort the merge map according to the read intervals
 * @author Steve Hoffmann
 *
 */


int cmp_mergemap_qsort(const void *a, const void *b) {
  mergemap_t  *first = (mergemap_t*) a;
  mergemap_t  *secnd = (mergemap_t*) b;

  if (first->beg > secnd->beg) return 1;
  if (first->beg < secnd->beg) return -1;
  if (first->end > secnd->end) return 1;
  if (first->end < secnd->end) return -1;

  return 0;

}
/*--------------------------- se_mergeStatsUpdate ------------------------------
 *
 * @brief update the structure for read mapping statistics
 * @author Steve Hoffmann
 *
 */
void
se_mergeStatsUpdate( bl_mergefilematch_t **list, Uint n, Uint pairingstate,
    char hasMate, mappingstats_t *stats, pthread_mutex_t *stats_mtx) {
  Uint noofqueries;
  Uint noofmates;

  bl_mergeCountAlignedMappings(list, n, &noofqueries, &noofmates);
  /* LOCK THE OUTPUT MTX9*/
  pthread_mutex_lock(stats_mtx);

  stats->total +=1;

  if (hasMate){
    stats->total +=1;
  }

  if(pairingstate == PAIRED) {
    //mapped in pair (+2 mapped reads)
    stats->mapped+=2;
    stats->paired+=1;

    if (noofqueries > 1){
      //multiple pair (+2 multiple mapped reads)
      stats->multiplemapped+=2;
      stats->multiplepaired+=1;
    } else {
      //unique pair (+2 uniquely mapped reads)
      stats->uniquemapped+=2;
      stats->uniquepaired+=1;
    }
  } else {
    //unmapped in pair
    if(noofqueries > 0){
      //query mapped (+1 mapped reads)
      stats->mapped+=1;
      if (hasMate){
        stats->singlequerymapped+=1;
      }
      if(noofqueries > 1) {
        stats->multiplemapped+=1;
      } else {
        stats->uniquemapped+=1;
      }
    } else {
      stats->unmapped+=1;
    }
    if(noofmates > 0){
      //mate mapped (+1 mapped reads)
      stats->mapped+=1;
      if (hasMate){
        stats->singlematemapped+=1;
      }
      if(noofmates > 1) {
        stats->multiplemapped+=1;
      } else {
        stats->uniquemapped+=1;
      }
    } else if(hasMate) {
      stats->unmapped+=1;
    }
  }

  pthread_mutex_unlock(stats_mtx);

  return;
}
/*--------------------------- se_mergeWriteOutput ------------------------------
 *
 * @brief write the output to the selected output device
 * @author Steve Hoffmann
 *
 */

void
se_mergeWriteOutput( bl_mergefilematch_t **best, Uint noofbest,
    mappingstats_t *stats, bl_fileBinDomains_t *odms, segemehl_t *nfo) {

  pthread_mutex_t *write_mtx;
  Uint noofqueries;
  Uint noofmates;
  Uint j, q, l;
  FILE *fp;

  write_mtx = nfo->mtx9;

  q = 0, l = 0;
  bl_mergeCountMappings(best, noofbest, &noofqueries, &noofmates);
  for (j = 0; j < noofbest; j++){
    /* updating SAM flag and tags */
    if (best[j]->read != NULL){
      if (noofqueries > 1){
        best[j]->read->flag |= 0x100;
      }
      bl_mergeUpdateTag(best[j]->read, q++, noofqueries);
    }
    if (best[j]->mate != NULL){
      if (noofmates > 1){
        best[j]->mate->flag |= 0x100;
      }
      bl_mergeUpdateTag(best[j]->mate, l++, noofmates);
    }

    /* report output */

    if(best[j]->read != NULL){


      if(!nfo->bam) {
        pthread_mutex_lock(write_mtx);
        fp = nfo->dev;
        /* select output device in case of chrdomains */
        if (odms != NULL) {
          fp = bl_fileBinsOpen(NULL, bl_fileBinsDomainGetBin(odms,
                best[j]->read->rname, best[j]->read->pos), "w");
        }
        bl_samprintSamrec(fp, best[j]->read, '\n');
        pthread_mutex_unlock(write_mtx);
      } else {
        bl_bamPrintBamrec (nfo->bamdev, best[j]->read, nfo->bamhdr, write_mtx);
      }
    }

    if(best[j]->mate != NULL){

      if(!nfo->bam) {
        pthread_mutex_lock(write_mtx);
        fp = nfo->dev;
        /* select output device in case of chrdomains */
        if (odms != NULL){
          fp = bl_fileBinsOpen(NULL, bl_fileBinsDomainGetBin(odms,
                best[j]->mate->rname, best[j]->mate->pos), "w");
        }
        bl_samprintSamrec(fp, best[j]->mate, '\n');
        pthread_mutex_unlock(write_mtx);
      } else {
        bl_bamPrintBamrec (nfo->bamdev, best[j]->mate, nfo->bamhdr, write_mtx);
      }

    }
  }
  return;
}
/*----------------------------- se_initMergeInfo --------------------------------
 *
 * @brief update the structure for read mapping statistics
 * @author Steve Hoffmann
 *
 */
mergeinfo_t *
se_initMergeInfo(bl_fileBinDomains_t *dms, fasta_t *reads,
    samheader_t *head, FILE *dev, bl_fileBinDomains_t *odms, Uchar remove,
    Uint bestonly, mappingstats_t *stats, segemehl_t *nfo, mergemap_t *map,
    Uint mapsz, bl_mergefiles_t *files, pthread_mutex_t *mtx) {

  mergeinfo_t *mi = ALLOCMEMORY(NULL, NULL, mergeinfo_t, 1);

  mi->bestonly = bestonly;
  mi->dev = dev;
  mi->dms = dms;
  mi->head = head;
  mi->map = map;
  mi->mapsz = mapsz;
  mi->nfo = nfo;
  mi->odms = odms;
  mi->reads = reads;
  mi->remove = remove;
  mi->stats = stats;
  mi->files = files;
  mi->cur = 0;
  mi->mastermtx = mtx;
  mi->processed = 0;

  return mi;
}
/*--------------------------- se_mergeComplexMaster ---------------------------------
 *
 * @brief return the next interval in map data structure for threaded processing,
 *        instead of waiting for the subsequent interval, this master looks for the
 *        next set of files that is ready for processing of a new interval
 * @author Steve Hoffmann
 *
 */

Uint se_mergeComplexMaster(mergeinfo_t *mi) {
  Uint k = -1, v, u;
  Uint processed = 0, avail = 0;
  //int s;

  pthread_mutex_lock(mi->mastermtx);

  /*
   * iter until an interval can be processed
   * or all intervals have been processed
   */
  while(1) {
    processed = 0;
    for(v=0; v < mi->mapsz; v++) {

      //unprocessed interval
      if(!mi->map[v].processed) {

        //check all necessary ressoruces for availabliltiy
        avail = 0;
        for (u=0; u < mi->map[v].n; u++) {

          Uint j = mi->map[v].binids[u];

          if (mi->map[v].queuepos[u] == mi->queuebeg[j]) {
            avail++;
          }
        }
        //all ressources available
        if (avail == mi->map[v].n) {
          break;
        }
      } else {
        //increase processed counter
        processed++;
      }
    }
    /*
     * in case that there is an unprocessed interval
     * for which all ressoruces are available
     * break and report k
     */
    if(v != mi->mapsz) {
      k=v;
      break;
    }
    //context sensitive: this point should never
    //be reached if for loop did not make a full
    //cycle
    assert(processed == mi->processed);
    /*
     * break here in case all intervals have been
     * processed, merging has finished
     *
     */
    if(processed == mi->mapsz) {
      break;
    }
  }
  /*
   * free interval found
   *
   */
  if(k != -1) {
    //mark itnerval as processed
    mi->map[k].processed = 1;
    mi->processed += 1;
    for(Uint u=0; u < mi->map[k].n; u++) {
      Uint j = mi->map[k].binids[u];
      //increase the queuebeg for the ressources
      /*
       * using trylock here as a safeguard in
       * case something went wrong. All ressources
       * must be free for this interval by construction
       *
       */
      int s = pthread_mutex_trylock(mi->files->f[j].mtx);
      if (s != 0) {
        NFO("failed lock bin %d [%d,%d]\n", j, mi->map[k].beg, mi->map[k].end);
        handle_error_en(s, "pthread_mutex_trylock");
        exit(-1);
      }
    }
  }

  pthread_mutex_unlock(mi->mastermtx);

  return k;
}
/*------------------------------- bl_mergeWorker ---------------------------------
 *
 * @brief the worker processing and merging the interval delivered by master
 * @author Steve Hoffmann
 *
 */

void* se_mergeWorker(void *x) {

  Uint u, j, k, i, q;
  mergeinfo_t *mi = (mergeinfo_t*)x;
  mergemap_t *map = mi->map;
  bl_mergefiles_t *files = mi->files;


  bl_mergefilematch_t **best = NULL;
  Uint noofbest = 0;
  Uint allocbest = 1000;
  int cmp = 0;
  Uint pairingstate = 0;
  fasta_t *reads=NULL;

  while((k = se_mergeComplexMaster(mi)) != -1) {


    best = ALLOCMEMORY(NULL, NULL, bl_mergefilematch_t**, allocbest);

    if(mi->reads) {
      reads = bl_fastxCopyIndex (NULL, mi->reads, map[k].chk, 1);
    }

    Uint len=map[k].end+1-map[k].beg;

    char *curkey = NULL;
    Uint curlen = 0;
    //for(i=map[k].beg; i <= map[k].end; i++) [
    for(i=0; i < len; i++) {

      noofbest = 0;
      if(reads) {
        curkey = bl_fastaGetDescription(reads, i);
        curlen = bl_fastaGetDescriptionLength(reads, i);
      } else {
        if(curkey) FREEMEMORY(NULL, curkey);
        curkey = NULL;
        curlen = 0;
      }

      /*
       *all bins, i.e. files that have processed this read
       *
       * */
      for(u=0; u < map[k].n; u++) {
        j = map[k].binids[u];

        while (1){

          /* read next entry */
          if (!files->f[j].complete){
            bl_mergeReadNext(mi->head, &files->f[j]);
            if(!reads && !curkey && files->f[j].entry->qname) {
              curkey = strdup(files->f[j].entry->qname);
              curlen = strlen(curkey);
            }
          }
          /* no match left */
          if (!files->f[j].complete){

            break;
          }

          if(!reads && !curkey && files->f[j].entry->qname) {
            curkey = strdup(files->f[j].entry->qname);
            curlen = strlen(curkey);
          }

          if (curkey && ! bl_mergefileFastaIDCompare(curkey, curlen,
                files->f[j].entry->qname, strlen(files->f[j].entry->qname))){
            break;
          }

          /* compare current with previous
           * best match based on pairing
           * state (no pair, pair, proper pair)
           * */

          if (noofbest > 0){
            cmp = bl_mergefilematchComparePairingState(best[0],
                files->f[j].entry, &pairingstate);

            if (cmp > 0){
              /* new best match found => destruct previous ones */
              for (q = 0; q < noofbest; q++){
                bl_mergefilematchDestruct(best[q]);
                FREEMEMORY(NULL, best[q]);
              }
              noofbest = 0;
            }
          }
          else {
            cmp = 0;
          }

          if (cmp >= 0){
            /* extend best buffer */
            if (noofbest == allocbest - 1){
              allocbest *= 2;
              best = ALLOCMEMORY(NULL, best, bl_mergefilematch_t **, allocbest);
            }
            /* append current match to best */
            best[noofbest++] = files->f[j].entry;

            files->f[j].entry = ALLOCMEMORY(NULL, NULL, bl_mergefilematch_t, 1);
            bl_mergefilematchInit(files->f[j].entry);
          }
          /* better match already found => clear data */
          else {
            bl_mergefilematchDestruct(files->f[j].entry);
          }

          /* remove redundant unmapped alignments and next information */
          if (noofbest > 1 && pairingstate == NOT_PAIRED){
            bl_mergefilematchRemoveRedundantUnmapped(best, noofbest);
            bl_mergefilematchRemoveNextInformation(best, noofbest);
          }

          /* apply best-only */
          if (noofbest > 1 && mi->bestonly){
            /* compare pair edit distance in case of PAIR or PROPER_PAIR */
            if (pairingstate != NOT_PAIRED)
            {
              bl_mergefilematchRemoveSuboptimalPairs(best, noofbest);
            }
            /* otherwise: compare edist distance for reads and mates separately*/
            else {
              bl_mergefilematchRemoveSuboptimalSingletons(best, noofbest);
            }
          }

          files->f[j].complete = 0;
        }

      }
      /* STATS */
      char hasMate = 1;
      if(reads) {
        hasMate = bl_fastaHasMate(reads);
      }
      se_mergeStatsUpdate(best, noofbest, pairingstate, hasMate, mi->stats,
          mi->nfo->mtx9);

      /* REPORTING */
      se_mergeWriteOutput(best, noofbest, mi->stats, mi->odms, mi->nfo) ;


      for (j = 0; j < noofbest; j++){
        /* clear data */
        bl_mergefilematchDestruct(best[j]);
        FREEMEMORY(NULL, best[j]);
      }

    } /*end of read iteration*/

    if(!reads && curkey) FREEMEMORY(NULL, curkey);
    if(reads) {
      bl_fastxDestructSequence(NULL, reads);
      bl_fastxDestructChunkIndex(NULL, reads);
      FREEMEMORY(NULL, reads);
    }
    FREEMEMORY(NULL, best);


    for(u=0; u < map[k].n; u++) {
      Uint j = map[k].binids[u];
      int s = pthread_mutex_unlock(files->f[j].mtx);
      if (s != 0) {
        handle_error_en(s, "pthread_mutex_unlock");
        exit(-1);
      }
      mi->queuebeg[j] += 1;
    }

  }

  return NULL;
}
/*------------------------------ se_mergeBisulfiteBins ---------------------------------
 *
 * @brief function to start the threaded merging of bisulfite bins
 * @author Steve Hoffmann
 *
 */
void se_mergeBisulfiteBinsNew(bl_fileBinDomains_t *dms, fasta_t *reads,
    samheader_t *head, FILE *dev, bl_fileBinDomains_t *odms, Uchar remove,
    Uint bestonly, mappingstats_t *stats, segemehl_t *nfo) {

  Uint i, j, k, noofbins, noofdms, mapsz=0;
  mergemap_t *map=NULL;
  bl_mergefiles_t files;
  pthread_t *threads;
  /*
   * ensure the integrity of domains as each has to have the
   * same number of bins, i.e. one bin for each alignment thread
   */
  assert(dms->noofdomains > 0);
  noofdms = dms->noofdomains;
  noofbins = dms->domain[0].bins.noofbins;
  for (i = 0; i < noofdms; i++){
    if (dms->domain[i].bins.noofbins != dms->domain[0].bins.noofbins){
      DBG("Inconsistent noofbins in domains. Exit forced.\n", NULL);
      exit(-1);
    }
  }
  /*
   * merging all the bins from all domains to a consecutive list
   *
   */
  NFO("Merging bisulfite bins now.\n", NULL);
  /* init and open files */
  bl_mergefilesInit(&files, noofbins*noofdms);
  k = 0;
  for (i = 0; i < noofdms; i++){
    for (j = 0; j < dms->domain[i].bins.noofbins; j++){
     // fprintf(stdout, "%d\t%d\t%d\t%d\t%s\n", i, dms->domain[i].domainsize, j, k,
     //     dms->domain[i].bins.b[j].fname);

      bl_mergefileInit(&files.f[k],
          bl_fileBinsOpen(NULL, &dms->domain[i].bins.b[j], "r"));

      k++;
    }
  }
  /*
   * get the sorted map of intervals and necessary files
   *
   */
  map = nfo->merge->map;
  mapsz = nfo->merge->mapsz;
  qsort(map, mapsz, sizeof(mergemap_t), cmp_mergemap_qsort);

  //map = se_mergeGetMap(nfo, noofbins, &mapsz);
  /*
   * sort bins and output test
   *
   */
  /*fprintf(stderr, "mapsize: %d\n", mapsz);
  for(k=0; k < mapsz; k++) {
    qsort(map[k].binids, map[k].n, sizeof(Uint), cmp_Uint_qsort);
    fprintf(stdout, "%d\t%d\t ", map[k].beg, map[k].end);
    for(Uint u=0; u < map[k].n; u++) {
      if(u>0) fprintf(stdout, ";");
      fprintf(stdout, "%d", map[k].binids[u]);
    }
    fprintf(stdout, "\n");
  }*/
  /*
   * determine the chunk position of each interval
   * within the respective file ressource
   *
   */
  Uint *queuebeg = calloc(2*nfo->threadno, sizeof(Uint));
  for(k=0; k < mapsz; k++) {
    for(j=0; j < 2; j++) {
      map[k].queuepos[j] = queuebeg[map[k].binids[j]];
      queuebeg[map[k].binids[j]]++;
    }
  }
  /*
   *  reset queuebeg as it is the pointer to the
   *  beginning of the queue in each file
   *
   */
  memset(queuebeg, 0, sizeof(Uint)*2*nfo->threadno);
  /*
   * init info structure
   *
   */
  pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;
  mergeinfo_t *w = se_initMergeInfo(dms, reads, head, dev, odms,
      remove, bestonly, stats, nfo, map, mapsz, &files, &mtx);

  w->queuebeg = queuebeg;

  pthread_attr_t attr;
  threads = ALLOCMEMORY(NULL, NULL, pthread_t, nfo->threadno);

  int s = pthread_attr_init(&attr);

  if (s != 0) {
    handle_error_en(s, "pthread_attr_init");
  }

  s = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  if (s != 0) {
    handle_error_en(s, "pthread_attr_setdetachstate");
  }

  NFO("merging with %u threads\n", nfo->threadno);
  for(i=0; i < nfo->threadno; i++) {
    s = pthread_create(&threads[i], NULL, se_mergeWorker, w);
    if(s != 0) {
      handle_error_en(s, "pthread_create");
    }
  }
  /*
   * join threads
   *
   */
  for(i=0; i < nfo->threadno; i++) {
    pthread_join(threads[i], NULL);
  }
  /*
   * check if files are entirely processed
   * and destruct the structrures
   *
   */
  for (j = 0; j < files.nooffiles; j++){
    /* check whether match file is entirely processed */
    bl_mergeReadNext(head, &files.f[j]);
    if (!files.f[j].eof || files.f[j].complete){
      DBG("File %d not yet entirely processed. Exit forced.\n", j);
      exit(-1);
    }
    /* destruct */
    bl_mergefileDestruct(&files.f[j]);
  }
  bl_mergefilesDestruct(&files);

  /*
   * close files and unlink
   *
   */

  for (i = 0; i < noofdms; i++){
    for (j = 0; j < noofbins; j++){
      bl_fileBinsClose(&dms->domain[i].bins.b[j]);
      if (remove){
        bl_rm(NULL, dms->domain[i].bins.b[j].fname);
        dms->domain[i].bins.b[j].unlinked = 1;
      }
    }
  }

  /*
   * clean up map
   *
   */
  for(k=0; k < mapsz; k++) {
    FREEMEMORY(NULL, map[k].binids);
    FREEMEMORY(NULL, map[k].queuepos);
  }

  FREEMEMORY(NULL, threads);
  FREEMEMORY(NULL, queuebeg);
  FREEMEMORY(NULL, w);
  FREEMEMORY(NULL, map);
  return;
}
