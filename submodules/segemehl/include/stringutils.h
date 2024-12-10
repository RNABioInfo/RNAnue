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

#ifndef STRINGUTILS_H
#define STRINGUTILS_H

#include <stdio.h>

#include "basic-types.h"

#ifndef ALLOCMEMORY
#include "memory.h"
#endif

#define TAB '\t' /*0x09*/
#define LF '\n'  /*0x0A*/
#define VT '\v'
#define FVAL '\f'
#define CR '\r' /*0x0D*/
#define SP ' '
#define DQT '\"'
#define SQT '\''

#define COPYSTR(M, D, S, L)             \
    D = ALLOCMEMORY(M, D, char, L + 1); \
    strncpy(D, S, L);                   \
    D[L] = '\0'

#define INDENT(s, x, c)                              \
    {                                                \
        int p;                                       \
        for (p = 0; p < x; p++) fprintf(s, "%c", c); \
    }
#define ISWHITESPACE(C) (C == SP || C == TAB || C == LF || C == VT || C == CR || C == FVAL)
#define ISQUOTE(C) (C == DQT || C == SQT)
#define SETSTR(X, I) (X)->strings[(I)].str
#define SETSTRLEN(X, I) (X)->strings[(I)].len

#define APPENDCHAR(S, X, L, Y)            \
    (X) = ALLOCMEMORY(S, X, char, L + 1); \
    X[L - 1] = Y;                         \
    X[L] = 0

typedef struct {
    char *str;
    Uint len;

} string_t;

typedef struct {
    string_t *strings;
    Uint noofstrings;

} stringset_t;

char *sprintflt(char **str, double flt);
char *sprintstr(char **str, char *src, Uint len);
char *sprintUint(char **str, Uint n);
char *sprintint(char **str, int n);
char *sprintchar(char **str, char chr);
char *strtok_bl(char *, char *, char **);
stringset_t *tokensToStringset(void *, char *, char *, Uint);
stringset_t *initStringset(void *);
char *strrev(char *str, Uint len);
char *strtrim(void *, char *, Uint *);
char *strtrimquote(void *, char *, Uint *);
char *strclip(void *, char *, Uint *);
void strconvert(char *, Uint, char, char);
void addString(void *, stringset_t *, char *, Uint);
void destructStringset(void *, stringset_t *);
char *concat(void *spacetab, char *strA, char *strB, int lenA, int lenB);
char *concatdelim(void *spacetab, char *strA, char *strB, int lenA, int lenB, char delim);
char *strreverse(char *, Uint);
char *my_itoa(int, char *, Uint);
char *bl_strdup(const char *s);
char *attachext(void *, char *, Uint, char *, Uint);
int checkmd5(unsigned char *a, unsigned char *b);
void fprintStringset(FILE *dev, stringset_t *set);
int bl_asprintf(char **strp, const char *fmt, ...);
int bl_bsprintf(char **strp, const char *fmt, ...);
int bl_strnfill(char **strp, int len, char ch);
#endif
