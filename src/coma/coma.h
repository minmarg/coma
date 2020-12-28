/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __coma_h__
#define __coma_h__


#define DEFAULT_SCORING_SCHEME          ( AbstractScoreMatrix::ProfileSpecific )
#define DEFAULT_STATISTICAL_BEHAVIOUR   ( AbstractScoreMatrix::ComputeStatistics )
#define DEFAULT_PRECISION               ( AbstractScoreMatrix::AutoScalling )
#define DEFAULT_MASKING                 ( MaskToIgnore )


// Version history:
//
// 0.2  . using integer schoring scheme; this enables fast and effective lambda and entropy computation
//      . computation of scaling parameter lambda has been changed; now it is more robust
//      . relative entropy computation has been introduced
// 0.3  . program has been added the statistical significance estimates derived using version 0.2
//          (integer scores, robust lambda computation)
// 0.4  . computation of the statistical significance parameters is now based on the effective sequence
//          profile lengths rather than on the raw lengths
// 0.5  . statistical significance calculations are fully changed and now analytical computations as in
//          psi-blast are used; searching against the database is included in this version
// 0.6  . minor adjustments including alignment of profile sequences and using fixed gap costs
// 0.7  . changed weighted frequency adjustment when all of them are zero in the column; now background
//          frequencies are used
// 0.8  . each PSSM is supposed to be scaled (or is scaled) before performing alignment procedure
// 0.9  . universal frequency vectors-based scoring system has been implemented here; the universal
//          scoring system uses two-level hash system to solve instant access problem
// 0.10 . thickness control added to the program usage of which enables to ignore positions containing
//          insufficient number of sequences
// 0.11 . configuration of parameters maintained through config files added
// 0.12 . bug of score computing cleared; there wesn't explicit type conversion to double within round function
// 0.13 . hashing of frequencies is changed; the key is now used as the combination of frequencies and
//          scores to overcome key-vale (frequencies as keys, scores as values) ambiguity problem which
//          arises when individual profiles are scaled before compiling a database and beginning of searching;
//          Also adjusted scoring method has been added: now each score matrix between profile pair is
//          scaled to attain lambda specific to the whole (big) profile-profile score matrix.
// 0.14 . gapped lambda patrameter is left unchanged after ungapped parameters are found;
//        fixed: real scores are used within ProfileMatrix instead of integers;
//        calculus of increased precision for statistical parameters and for alignment procedure has been
//        implemented
// 0.15 . internal structure of the program's scoring systems changed in order to employ more conveniently
//        integer score system with increased precision
// 0.16 . option of information content that controls alignement searching has been appended
// 0.17 . options of scoring schemes have been added
// 0.18 . SEG implementation and design for profiles has been added;
// 0.19 . Expectation option to compute expected mean alignment length has been added;
// 0.20 . Automatic gap open costs implemented by using autocorrelation of the max scores
// 0.21 . Scaling of masked positions has been introduced
// 0.22 . increased precision of gaps for any scoring scheme
// 0.23 . frequency vectors in database now with effective thickness added
// 0.24 . effective thickness of a profile column corresponds to the number of sequences in the extent
//          at the position
// 0.25 . fixed: following some version, original scores were left in input sequence profile despite its
//          scaling, so unscaled sequence profile version was always used
// 0.26 . full information included in sequence profile construction procedure, i.e. delete states are
//          now processed
// 0.27 . profile format changed; Alignment furnished with delete states 
// 0.28 . computation of weighted frequencies adjusted for profiles: now each extent encompasses at least a
//          given minimum number of positions
// 0.29 . processing of extents of one sequence has been changed: computing anyway
// 0.30 . new functionality and options added; 2nd pass incorporation, autocorrelation autocorrection
//          computations, automatic information content adjustment on the 2nd pass, etc.
// 0.31 . input alignment processing fixed so that beginning and trailing gaps in query sequence are properly processed
// 0.32 . gap probability factor depending on the effective thickness of profiles introduced
// 0.33 . relative entropy calculations for each position of score system derived and applied for gap cost computations
// 0.34 . option of positional relative entropy calculations added
// 0.35 . adjustments by context are now accomplished using evalue pair alignment
// 0.40 . naming changed; correction of information content threshold is determined by two curves;
//          statistical parameters of gapped alignments estimated
// 0.41 . usual SEG for each sequence in input multiple alignment added;
//          analytically computed entropies disabled in global score system
// 0.42 . output of statistical parameters slightly changed in cases of non-negative expected scores
// 0.43 . Expected value of Expect has been used to change `n/a' in the coma output
// ----
//
// 1.01 . heuristics of high-scoring pairs implemented
// 1.02 . file of options added
// 1.05 . format changed: new section of background probabilities
// 1.10 . optimization of profile's target frequencies implemented


static const char*  version = "1.10";
static const char*  verdate = "";



static const char*  instructions = "\n\
<>\n\
\n\
A protein sequence profile comparison and alignment searching tool.\n\
(C)2010 Mindaugas Margelevicius,Institute of Biotechnology,Vilnius\n\
\n\
Usage:\n\
<> -i <query> -d <database> [-o <output>] [-p <options>]\n\
\n\
Parameters:\n\
\n\
-i <query>      [Filename]  Either multiple alignment file in fasta or\n\
                            profile made by makepro.\n\
-d <database>   [Filename]  Either name of database made by makedb or\n\
                            another multiple alignment file.\n\
-o <output>     [Filename]  Output file of alignments.\n\
-p <options>    [Filename]  Input file of options;\n\
                            By default, file in configuration\n\
                            directory of this package is searched.\n\
\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
";

#endif//__coma_h__
