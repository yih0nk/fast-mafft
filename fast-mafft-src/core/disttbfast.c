#ifndef _GNU_SOURCE
#define _GNU_SOURCE // fast-mafft: glibc hides pthread_rwlock_t under -std=c99 without a feature macro
#endif
#include "mltaln.h"
#include "progtiming.h"
#include "lazygap.h"
#include "dphelp.h"
#ifdef PROGTIMING
typedef struct { int l, thread, clus1, clus2, len1, len2, c0, c1; double t[8]; double acc[PT_NACC]; } ptrec_t;
static ptrec_t *ptlog = NULL;
static TLS double pt_growt = 0.0; // work/v8s
static double pt_base = 0.0;
static void ptlog_write( int njob, int pass )
{
	char *fn = getenv( "MAFFT_PROGLOG" );
	FILE *fp;
	int i, k;
	if( !fn || !ptlog ) return;
	fp = fopen( fn, "a" );
	if( !fp ) return;
	for( i=0; i<njob-1; i++ )
	{
		ptrec_t *r = ptlog+i;
		fprintf( fp, "%d %d %d %d %d %d %d %d %d", pass, r->l, r->thread, r->c0, r->c1, r->clus1, r->clus2, r->len1, r->len2 );
		for( k=0; k<8; k++ ) fprintf( fp, " %.6f", r->t[k] );
		for( k=0; k<PT_NACC; k++ ) fprintf( fp, " %.6f", r->acc[k] );
		fprintf( fp, "\n" );
	}
	fclose( fp );
}
#define PTT(k) do { if( ptlog ) ptlog[l].t[k] = pt_now() - pt_base; } while( 0 )
#else
#define PTT(k)
#endif


#define REPORTCOSTS 0

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0
#define SKIP 1

#define ITERATIVECYCLE 2

#define END_OF_VEC -1

static int treein;
static int topin;
static int treeout;
static int noalign;
static int tuplesize;
static int subalignment;
static int subalignmentoffset;
static int nguidetree;
static int sparsepickup;
static int keeplength;
static int ndeleted;
static int mapout;
static int smoothing;
static double maxdistmtxsize;
static int nthreadtb;
static int nthreadtbq = -1;       // work/v7: threads for the progressive passes when the step order cannot change the result (-C a-b-c: c); -1 = nthreadtb
static int tb_origserial = 0;     // work/largen: the original would run the progressive passes in the serial treebase() (njob > 20000)
static size_t tb_keepbytes = 0;   // work/largen: > 0: after each step a thread frees its row and traceback arenas larger than this
static double tb_heavymin = 64.0 * 1048576.0; // work/largen: a step whose full DP traceback needs this many bytes is heavy (see treebasethread)
static int useexternalanchors;
static int oneiteration;

#if 0
#define PLENFACA 0.0123
#define PLENFACB 10252
#define PLENFACC 10822
#define PLENFACD 0.5
#define DLENFACA 0.01
#define DLENFACB 2445
#define DLENFACC 2412
#define DLENFACD 0.1
#else
#define PLENFACA 0.01
#define PLENFACB 10000
#define PLENFACC 10000
#define PLENFACD 0.1
#define D6LENFACA 0.01
#define D6LENFACB 2500
#define D6LENFACC 2500
#define D6LENFACD 0.1
#define D10LENFACA 0.01
#define D10LENFACB 1000000
#define D10LENFACC 1000000
#define D10LENFACD 0.0
#endif

typedef struct _jobtable
{
    int i;  
    int j;  
} Jobtable;

typedef struct _msacompactdistmtxthread_arg
{
	int njob;
	int thread_no;
	int *selfscore;
	double **partmtx;
	char **seq;
	int **skiptable;
	double *mindist;
	int *mindistfrom;
 	int *jobpospt;
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} msacompactdistmtxthread_arg_t;

typedef struct _compactdistmtxthread_arg
{
	int njob;
	int thread_no;
	int *nogaplen;
	int **pointt;
	int *selfscore;
	double **partmtx;
	int *jobpospt;
	double *mindist;
	int *mindistfrom;
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} compactdistmtxthread_arg_t;

typedef struct _msadistmtxthread_arg
{
	int njob;
	int thread_no;
	int *selfscore;
	double **iscore;
	double **partmtx;
	char **seq;
	int **skiptable;
	Jobtable *jobpospt;
	void *fast; // msadistfast_t *: exact fast scorer; NULL -> original naivepairscorefast()
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} msadistmtxthread_arg_t;

#ifdef enablemultithread
// ue futatsu ha singlethread demo tsukau
typedef struct _treebasethread_arg
{
	int thread_no;
	int njob;
	int *nrunpt;
	int *nlen;
	int *jobpospt;
	int ***topol;
	Treedep *dep;
	double ***cpmxhist;
	int **memhist;
	char **aseq;
	double *effarr;
	int *alloclenpt;
	int *fftlog;
	char *mergeoralign;
	double **newdistmtx;
	int *selfscore;
	ExtAnch *extanch;
	int **anchindex;
	pthread_mutex_t *mutex;
	pthread_cond_t *treecond;
	// work/prog: ready-queue scheduling (see tbsched_* below)
	struct _tbsched *sched;
	// work/prog: guards the rows of aseq against ReallocateCharMtx (write lock) while member rows
	// are copied in and out (read lock), so those copies no longer hold the global mutex
	pthread_rwlock_t *aseqlock;
	// work/prog3: the steps work on the rows of aseq in place (no local copies); each row grows on
	// its own when a step that owns it needs more room (rowcap: allocated bytes of each row)
	int *rowcap;
	int nocopy;
	struct _lzstate *lz; // work/prog3: deferred gap insertion (lazygap.h), NULL = off
	int nrunmax; // work/v7: threads of this attempt (the limit on steps in flight; was nthreadtb)
	int helperonly; // work/v10w: in-order attempts: a thread beyond the original's nthreadtb, which never takes a step and only helps
} treebasethread_arg_t;

typedef struct _distancematrixthread_arg
{
	int thread_no;
	int njob;
	int *jobpospt;
	int **pointt;
	double **mtx;
	pthread_mutex_t *mutex;
} distancematrixthread_arg_t;
#endif


void arguments( int argc, char *argv[] )
{
    int c;

	nthread = 1;
	nthreadpair = 1;
	nthreadtb = 1;
	outnumber = 0;
	topin = 0;
	treein = 0;
	treeout = 0;
	distout = 0;
	noalign = 0;
	nevermemsave = 0;
	inputfile = NULL;
	nadd = 0;
	addprofile = 1;
	fftkeika = 0;
	constraint = 0;
	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	devide = 0;
	use_fft = 0;
	useexternalanchors = 0;
	oneiteration = 0;
	force_fft = 0;
	fftscore = 1;
	fftRepeatStop = 0;
	fftNoAnchStop = 0;
    weight = 3;
    utree = 1;
	tbutree = 1;
    refine = 0;
    check = 1;
    cut = 0.0;
    disp = 0;
    outgap = 1;
    alg = 'A';
    mix = 0;
	tbitr = 0;
	scmtd = 5;
	tbweight = 0;
	tbrweight = 3;
	checkC = 0;
	treemethod = 'X';
	sueff_global = 0.1;
	contin = 0;
	scoremtx = 1;
	kobetsubunkatsu = 0;
	dorp = NOTSPECIFIED;
	ppenalty_dist = NOTSPECIFIED;
	ppenalty = -1530;
	ppenalty_ex = NOTSPECIFIED;
	penalty_shift_factor = 1000.0;
	poffset = -123;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	TMorJTT = JTT;
	scoreout = 0;
	spscoreout = 0;
	tuplesize = 6;
	subalignment = 0;
	subalignmentoffset = 0;
	legacygapcost = 0;
	specificityconsideration = 0.0;
	nguidetree = 1;
	sparsepickup = 0;
	keeplength = 0;
	mapout = 0;
	smoothing = 0;
	nwildcard = 0;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					reporterr(       "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'I':
					nadd = myatoi( *++argv );
					reporterr(       "nadd = %d\n", nadd );
					--argc;
					goto nextoption;
				case 'V':
					ppenalty_dist = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
//					reporterr(       "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'Q':
					penalty_shift_factor = atof( *++argv );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					reporterr(       "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
//					reporterr(       "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
					reporterr(       "kappa = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = myatoi( *++argv );
					scoremtx = 1;
//					reporterr(       "blosum %d / kimura 200 \n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					reporterr(       "jtt/kimura %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
					reporterr(       "tm %d\n", pamN );
					--argc;
					goto nextoption;
				case 'C':
					nthreadpair = nthread = myatoi( *++argv );
					reporterr(       "nthread = %d\n", nthread );
					reporterr(       "nthreadpair = %d\n", nthread );
					if( strchr( *argv, '-' ) )
					{
						nthreadtb = myatoi( strchr( *argv, '-' )+1 );
						/* work/v7: -C a-b-c.  c = threads for the progressive passes whenever their result
						   cannot depend on the step order (see nthreadtbq in disttbfast()); b stays what the
						   original uses, and the order-dependent (in-order) dispatch keeps it.  mafft.tmpl
						   passes c = --thread when it caps b at 16. */
						if( strchr( strchr( *argv, '-' )+1, '-' ) )
							nthreadtbq = myatoi( strchr( strchr( *argv, '-' )+1, '-' )+1 );
					}
					else
						nthreadtb = nthread;
					reporterr(       "nthreadtb = %d\n", nthreadtb );
					--argc; 
					goto nextoption;
				case 's':
					specificityconsideration = (double)myatof( *++argv );
//					reporterr(       "specificityconsideration = %f\n", specificityconsideration );
					--argc; 
					goto nextoption;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'K':
					addprofile = 0;
					break;
				case 'y':
                                        distout = *(*++argv);
					reporterr(       "distout=%c\n", distout );
                                        if( distout != 'c' && distout != 'h' )
                                        {
                                            reporterr(       "Set -y c or -y h in v>=7.521.\n" );
                                            exit( 1 );
                                        }
					--argc; 
					goto nextoption;
				case 't':
					treeout = 1;
					break;
				case '^':
					treeout = 2;
					break;
				case 'T':
					noalign = 1;
					break;
				case 'r':
					oneiteration = 1;
					break;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
				case 'L':
					legacygapcost = 1;
					break;
				case 'e':
					fftscore = 0;
					break;
				case 'x':
					terminalmargin = myatoi( *++argv );
					--argc; 
					goto nextoption;
				case 'H':
					subalignment = 1;
					subalignmentoffset = myatoi( *++argv );
					--argc;
					goto nextoption;
#if 0
				case 'R':
					fftRepeatStop = 1;
					break;
#endif
				case 'n' :
					outnumber = 1;
					break;
#if 0
				case 's':
					treemethod = 's';
					break;
				case 'q':
					treemethod = 'q'; // minimum
					break;
#endif
				case 'q':
					sparsepickup = myatoi( *++argv );
//					reporterr(       "sparsepickup = %d\n", sparsepickup );
					--argc; 
					goto nextoption;
				case 'X':
					treemethod = 'X';
					sueff_global = atof( *++argv );
//					fprintf( stderr, "sueff_global = %f\n", sueff_global );
					--argc;
					goto nextoption;
				case 'E':
					nguidetree = myatoi( *++argv );
//					reporterr(       "nguidetree = %d\n", nguidetree );
					--argc; 
					goto nextoption;
#if 0
				case 'a':
					alg = 'a';
					break;
				case 'H':
					alg = 'H';
					break;
				case 'R':
					alg = 'R';
					break;
#endif
				case 'A':
					alg = 'A';
					break;
				case '&':
					alg = 'a';
					break;
				case '@':
					alg = 'd';
					break;
				case 'N':
					nevermemsave = 1;
					break;
				case 'M':
					alg = 'M';
					break;
#if 0
				case 'S' :
					scoreout = 1; // for checking parallel calculation
					break;
#else
				case 'S' :
					spscoreout = 1; // 2014/Dec/30, sp score
					break;
#endif
				case 'B': // hitsuyou! memopt -M -B no tame
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'l':
					useexternalanchors = 1;
				case 'G':
					use_fft = 1;
					force_fft = 1;
					break;
#if 0
				case 'V':
					topin = 1;
					break;
#endif
				case 'U':
					treein = 1;
					break;
				case 'u':
					weight = 0;
					tbrweight = 0;
					break;
				case 'v':
					tbrweight = 3;
					break;
#if 1
				case 'd':
					disp = 1;
					break;
#endif
#if 1
				case 'O':
					outgap = 0;
					break;
#else
				case 'O':
					fftNoAnchStop = 1;
					break;
#endif
				case 'J':
					tbutree = 0;
					break;
#if 0
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
#endif
				case 'w':
					fftWinSize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'W':
					tuplesize = myatoi( *++argv );
					--argc;
					goto nextoption;
#if 0
				case 'Z':
					checkC = 1;
					break;
#endif
				case 'Y':
					keeplength = 1;
					break;
				case 'z':
					mapout = 2;
					break;
				case 'Z':
					mapout = 1;
					break;
				case 'p':
					smoothing = 1;
					break;
				case ':':
					nwildcard = 1;
					break;
                default:
                    reporterr(       "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
		}
		nextoption:
			;
	}
    if( argc == 1 )
    {
        cut = atof( (*argv) );
        argc--;
    }
    if( argc != 0 ) 
    {
        reporterr(       "options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		reporterr(       "conflicting options : o, m or u\n" );
		exit( 1 );
	}
}

static int varpairscore( int nseq, int npick, int nlenmax, char **seq, int seed )
{
	int i, j, npair;
	int *slist;
	char **pickseq;
	double score;
	double scoreav;
	double scoreav2;
	double scorestd;
	double scorevar;
	slist = calloc( nseq, sizeof( int ) );
	pickseq = AllocateCharMtx( npick, nlenmax );
	reporterr( "nseq = %d, nlenmax=%d, seed=%d\n", nseq, nlenmax, seed );

	srand( seed );

	for( i=0; i<nseq; i++ ) slist[i] = i;
//	for( i=0; i<nseq; i++ ) reporterr( "slist[%d] = %d\n", i, slist[i] );

	stringshuffle( slist, nseq );
	for( i=0; i<npick; i++ ) gappick0( pickseq[i], seq[slist[i]] );

	scoreav = 0.0;
	scoreav2 = 0.0;
	npair = npick * (npick-1) / 2;
	for( i=1; i<npick; i++ ) 
	{
		reporterr( "%d / %d\r", i, npick );
		for( j=0; j<i; j++ ) 
		{
			score = G__align11_noalign( n_dis_consweight_multi, -1200, -60, pickseq+i, pickseq+j, nlenmax );
			scoreav += score;
			scoreav2 += score * score;
			printf( "score = %d\n", (int)score );
		}
	}

	scoreav /= (double)npair;
	scoreav2 /= (double)npair;
	scorevar = ( scoreav2 - scoreav * scoreav )*npair/(npair-1);
	scorestd = sqrt( scorevar );
	printf( "av = %f\n", scoreav );
	printf( "stddev = %f\n", scorestd );
	printf( "cv = %f\n", scorestd/scoreav );

	FreeCharMtx( pickseq );

	if( scorestd/scoreav < 0.2 ) return( 's' );
	else return( 't' );
}

static void pickup( int n, int *seqlen, int ***topol, char **name, char **seq ) // memsave ni mitaiou
{
	int i, j, k, m;
	int **longestseq;
	int **longestlen;
	int *select;
	char **nameout, **seqout;
	int *nlenout;
	char **namenotused, **seqnotused;
	int *nlennotused;
	FILE *notusedfp;

	longestseq = AllocateIntMtx( n-1, 2 );
	longestlen = AllocateIntMtx( n-1, 2 );
	select = AllocateIntVec( n );
	for( i=0; i<n; i++ ) select[i] = 0;
	nameout = AllocateCharMtx( n, 0 );
	seqout = AllocateCharMtx( n, 0 );
	nlenout = AllocateIntVec( n );
	namenotused = AllocateCharMtx( n, 0 );
	seqnotused = AllocateCharMtx( n, 0 );
	nlennotused = AllocateIntVec( n );

	for( i=0; i<n-1; i++ )
	{
//		reporterr( "STEP %d\n", i );
		longestlen[i][0] = -1;
		longestseq[i][0] = -1;
		for( j=0; (m=topol[i][0][j])!=-1; j++ ) // sukoshi muda
		{
			if( seqlen[m] > longestlen[i][0] )
			{
				longestlen[i][0] = seqlen[m];
				longestseq[i][0] = m;
			}
//			reporterr( "%d ", topol[i][0][j] );
		}
//		reporterr( "longest = %d (%d)\n", longestlen[i][0], longestseq[i][0] );


		longestlen[i][1] = -1;
		longestseq[i][1] = -1;
		for( j=0; (m=topol[i][1][j])!=-1; j++ ) // sukoshi muda
		{
			if( seqlen[m] > longestlen[i][1] )
			{
				longestlen[i][1] = seqlen[m];
				longestseq[i][1] = m;
			}
//			reporterr( "%d ", topol[i][1][j] );
		}
//		reporterr( "longest = %d (%d)\n", longestlen[i][1], longestseq[i][1] );
	}

	m = 1;
	for( i=n-2; i>-1; i-- )
	{
//		reporterr( "longest[%d][0] = %d (%d)\n", i, longestlen[i][0], longestseq[i][0] );
//		reporterr( "longest[%d][1] = %d (%d)\n", i, longestlen[i][1], longestseq[i][1] );
		select[longestseq[i][0]] = 1;
		select[longestseq[i][1]] = 1;
		m += 1;
		if( m >= sparsepickup ) break;
	}
	for( i=0, k=0, j=0; i<n; i++ ) 
	{
		if( select[i] )
		{
			nameout[k] = name[i];
			seqout[k] = seq[i];
			nlenout[k] = strlen( seqout[k] );
			k++;
		}
		else
		{
			namenotused[j] = name[i];
			seqnotused[j] = seq[i];
			nlennotused[j] = strlen( seqnotused[j] );
			j++;
		}
	}
	writeData_pointer( stdout, m, nameout, nlenout, seqout );

	notusedfp = fopen( "notused", "w" );
	writeData_pointer( notusedfp, n-m, namenotused, nlennotused, seqnotused );
	fclose( notusedfp );


	free( nameout );
	free( nlenout );
	free( seqout );
	free( namenotused );
	free( nlennotused );
	free( seqnotused );
	FreeIntMtx( longestseq );
	FreeIntMtx( longestlen );
	free( select );
}


static int nunknown = 0;

void seq_grp_nuc( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++];
		if( tmp < 4 )
			*grp++ = tmp;
		else
			nunknown++;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < tuplesize )
	{
//		reporterr(       "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

void seq_grp( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(unsigned char)*seq++];
		if( tmp < 6 )
			*grp++ = tmp;
		else
			nunknown++;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < 6 )
	{
//		reporterr(       "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

void makecompositiontable_p( int *table, int *pointt )
{
	int point;

	while( ( point = *pointt++ ) != END_OF_VEC )
	{
#if 1
		table[point]++;
#else // kakunin shinai
		if( (unsigned int)table[point]++ >= INT_MAX )
		{
			reporterr( "Overflow. table[point]=%d>INT_MAX(%d).\n", table[point], INT_MAX );
			exit( 1 );
		}
#endif
	}
}


void makepointtable_nuc_dectet( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *262144;
	point += *n++ * 65536;
	point += *n++ * 16384;
	point += *n++ *  4096;
	point += *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ *262144;
		point *= 4;
		point += *n++;
		*pointt++ = point;

	}
	*pointt = END_OF_VEC;
}

void makepointtable_nuc_octet( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ * 16384;
	point += *n++ *  4096;
	point += *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 16384;
		point *= 4;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

void makepointtable_nuc( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 1024;
		point *= 4;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

void makepointtable( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *  7776;
	point += *n++ *  1296;
	point += *n++ *   216;
	point += *n++ *    36;
	point += *n++ *     6;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 7776;
		point *= 6;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

static double preferenceval( int ori, int pos, int max ) // for debug
{
	pos -= ori;
	if( pos < 0 ) pos += max;
	return( 0.00000000000001 * pos );
}

static void *compactdisthalfmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	compactdistmtxthread_arg_t *targ = (compactdistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **partmtx = targ->partmtx;
	int *nogaplen = targ->nogaplen;
	int **pointt = targ->pointt;
 	int *jobpospt = targ->jobpospt;
	double *mindist = targ->mindist;
	int *mindistfrom = targ->mindistfrom;
	int i, j;
	double tmpdist, preference, tmpdistx; //, tmpdisty;
	int *table1;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthreadpair ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = *jobpospt;
			if( i == -1 )
			{
				pthread_mutex_unlock( targ->mutex );
				commonsextet_p( NULL, NULL );
				return( NULL );
			}
			*jobpospt = i-1;
			pthread_mutex_unlock( targ->mutex );
		}
		else
#endif 
		{
			i = *jobpospt;
			if( i == -1 )
			{
				commonsextet_p( NULL, NULL );
				return( NULL );
			}
			*jobpospt = i-1;
		}

		table1 = (int *)calloc( tsize, sizeof( int ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		if( i % 100 == 0 )
		{
			if( nthreadpair )
				reporterr(       "\r% 5d / %d (thread %4d)", njob-i, njob, thread_no );
			else
				reporterr(       "\r% 5d / %d", njob-i, njob );
		}
		makecompositiontable_p( table1, pointt[i] );

//		for( j=i+1; j<njob; j++ ) 
		for( j=i-1; j>-1; j-- ) 
		{

			tmpdist = distcompact( nogaplen[i], nogaplen[j], table1, pointt[j], selfscore[i], selfscore[j] );
			preference = preferenceval( i, j, njob );
			tmpdistx = tmpdist + preference;
			if( tmpdistx < mindist[i] )
			{
				mindist[i] = tmpdistx;
				mindistfrom[i] = j;
			}

//			preference = preferenceval( j, i, njob );
//			tmpdisty = tmpdist + preference;
//			if( tmpdisty < mindist[j] )
//			{
//				mindist[j] = tmpdisty;
//				mindistfrom[j] = i;
//			}

			if( partmtx[i] ) partmtx[i][j] = tmpdist;
			if( partmtx[j] ) partmtx[j][i] = tmpdist;
		} 
		free( table1 );
	}
}

static void *ylcompactdisthalfmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	compactdistmtxthread_arg_t *targ = (compactdistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **partmtx = targ->partmtx;
	int *nogaplen = targ->nogaplen;
	int **pointt = targ->pointt;
 	int *jobpospt = targ->jobpospt;
	double *mindist = targ->mindist;
	int *mindistfrom = targ->mindistfrom;
	int i, j;
	double tmpdist, preference, tmpdistx, tmpdisty;
	int *table1;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthreadpair ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = *jobpospt;
			if( i == njob-1 )
			{
				pthread_mutex_unlock( targ->mutex );
				commonsextet_p( NULL, NULL );
				return( NULL );
			}
			*jobpospt = i+1;
			pthread_mutex_unlock( targ->mutex );
		}
		else
#endif 
		{
			i = *jobpospt;
			if( i == njob-1 )
			{
				commonsextet_p( NULL, NULL );
				return( NULL );
			}
			*jobpospt = i+1;
		}

		table1 = (int *)calloc( tsize, sizeof( int ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		if( i % 100 == 0 )
		{
			if( nthreadpair )
				reporterr(       "\r% 5d / %d (thread %4d)", i+1, njob, thread_no );
			else
				reporterr(       "\r% 5d / %d", i+1, njob );
		}
		makecompositiontable_p( table1, pointt[i] );

		for( j=i+1; j<njob; j++ ) 
//		for( j=i-1; j>-1; j-- ) 
		{

			tmpdist = distcompact( nogaplen[i], nogaplen[j], table1, pointt[j], selfscore[i], selfscore[j] );
			preference = preferenceval( i, j, njob );
			tmpdistx = tmpdist + preference;
			if( tmpdistx < mindist[i] )
			{
				mindist[i] = tmpdistx;
				mindistfrom[i] = j;
			}

			preference = preferenceval( j, i, njob );
			tmpdisty = tmpdist + preference;
			if( tmpdisty < mindist[j] )
			{
				mindist[j] = tmpdisty;
				mindistfrom[j] = i;
			}

			if( partmtx[i] ) partmtx[i][j] = tmpdist;
			if( partmtx[j] ) partmtx[j][i] = tmpdist;
		} 
		free( table1 );
	}
}


static void *msacompactdisthalfmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	msacompactdistmtxthread_arg_t *targ = (msacompactdistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **partmtx = targ->partmtx;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	double *mindist = targ->mindist;
	int *mindistfrom = targ->mindistfrom;
 	int *jobpospt = targ->jobpospt;
	double tmpdist, preference, tmpdistx; //, tmpdisty;
	int i, j;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthreadpair ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = *jobpospt;
			if( i == -1 )
			{
				pthread_mutex_unlock( targ->mutex );
				return( NULL );
			}
			*jobpospt = i-1;
			pthread_mutex_unlock( targ->mutex );
		}
		else
#endif
		{
			i = *jobpospt;
			if( i == -1 )
			{
				return( NULL );
			}
			*jobpospt = i-1;
		}

		if( i % 100 == 0 ) 
		{
			if( nthreadpair )
				fprintf( stderr, "\r% 5d / %d (thread %4d)", njob-i, njob, thread_no );
			else
				fprintf( stderr, "\r% 5d / %d", i, njob );
		}

		for( j=i-1; j>-1; j-- ) 
//		for( j=i+1; j<njob; j++ ) 
		{
			tmpdist = distcompact_msa( seq[i], seq[j], skiptable[i], skiptable[j], selfscore[i], selfscore[j] ); // osoikedo,

			preference = preferenceval( i, j, njob );
			tmpdistx = tmpdist + preference;
			if( tmpdistx < mindist[i] )
			{
				mindist[i] = tmpdistx;
				mindistfrom[i] = j;
			}

//			preference = preferenceval( j, i, njob );
//			tmpdisty = tmpdist + preference;
//			if( tmpdisty < mindist[j] )
//			{
//				mindist[j] = tmpdisty;
//				mindistfrom[j] = i;
//			}
			if( partmtx[i] ) partmtx[i][j] = tmpdist;
			if( partmtx[j] ) partmtx[j][i] = tmpdist;
		}
	}
}

static void *ylmsacompactdisthalfmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	msacompactdistmtxthread_arg_t *targ = (msacompactdistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **partmtx = targ->partmtx;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	double *mindist = targ->mindist;
	int *mindistfrom = targ->mindistfrom;
 	int *jobpospt = targ->jobpospt;
	double tmpdist, preference, tmpdistx, tmpdisty;
	int i, j;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthreadpair ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = *jobpospt;
			if( i == njob-1 )
			{
				pthread_mutex_unlock( targ->mutex );
				return( NULL );
			}
			*jobpospt = i+1;
			pthread_mutex_unlock( targ->mutex );
		}
		else
#endif
		{
			i = *jobpospt;
			if( i == njob-1 )
			{
				return( NULL );
			}
			*jobpospt = i+1;
		}

		if( i % 100 == 0 ) 
		{
			if( nthreadpair )
				fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
			else
				fprintf( stderr, "\r% 5d / %d", i, njob );
		}

//		for( j=i-1; j>-1; j-- ) 
		for( j=i+1; j<njob; j++ ) 
		{
			tmpdist = distcompact_msa( seq[i], seq[j], skiptable[i], skiptable[j], selfscore[i], selfscore[j] ); // osoikedo,

			preference = preferenceval( i, j, njob );
			tmpdistx = tmpdist + preference;
			if( tmpdistx < mindist[i] )
			{
				mindist[i] = tmpdistx;
				mindistfrom[i] = j;
			}

			preference = preferenceval( j, i, njob );
			tmpdisty = tmpdist + preference;
			if( tmpdisty < mindist[j] )
			{
				mindist[j] = tmpdisty;
				mindistfrom[j] = i;
			}
			if( partmtx[i] ) partmtx[i][j] = tmpdist;
			if( partmtx[j] ) partmtx[j][i] = tmpdist;
		}
	}
}

/*
 * fast-mafft (work/x86simd): instruction-set level of the two distance kernels
 * below (msadistfast, sextetfast).
 *   level 0  portable C
 *   level 1  128-bit SIMD: NEON on arm64; on x86-64 SSE2 (sextet, always
 *            available) and POPCNT (msadist; portable C if the CPU lacks it)
 *   level 2  x86-64 AVX2 (msadist also needs POPCNT and BMI1)
 * On x86-64 (gcc or clang) the SSE2/POPCNT/AVX2 kernels are compiled with
 * __attribute__((target(...))), so a default build (no -m flags) has them, and
 * the level is chosen at run time with __builtin_cpu_supports().  Only integer
 * kernels carry a target attribute; every double operation stays in default
 * code, so no FMA contraction or other rounding change can occur.  The level is
 * fixed before any worker thread starts.
 *   MAFFT_SIMD=scalar|sse2|neon|popcnt|avx2   lower the level (debug / A-B);
 *   a level the CPU lacks is not used (warning)
 */
/* BEGIN fmsimd dispatch */
#include <stdint.h>
#if defined(__x86_64__) && ( defined(__GNUC__) || defined(__clang__) )
// mltaln.h defines one-letter macros (M N B C D G ...) that some compilers'
// intrinsics headers use as parameter names (Apple clang 15: B)
#pragma push_macro( "M" )
#pragma push_macro( "N" )
#pragma push_macro( "B" )
#pragma push_macro( "C" )
#pragma push_macro( "D" )
#pragma push_macro( "G" )
#pragma push_macro( "TM" )
#pragma push_macro( "WIN" )
#pragma push_macro( "MAX" )
#pragma push_macro( "MIN" )
#undef M
#undef N
#undef B
#undef C
#undef D
#undef G
#undef TM
#undef WIN
#undef MAX
#undef MIN
#include <immintrin.h>
#pragma pop_macro( "MIN" )
#pragma pop_macro( "MAX" )
#pragma pop_macro( "WIN" )
#pragma pop_macro( "TM" )
#pragma pop_macro( "G" )
#pragma pop_macro( "D" )
#pragma pop_macro( "C" )
#pragma pop_macro( "B" )
#pragma pop_macro( "N" )
#pragma pop_macro( "M" )
#define FMSIMD_X86 1
#else
#define FMSIMD_X86 0
#endif
#define FMSIMD_SCALAR 0
#define FMSIMD_BASE 1
#define FMSIMD_AVX2 2

// the highest level this CPU runs; *popcnt = x86 POPCNT present, *avx2ok = the
// msadist AVX2 kernel can run (AVX2, POPCNT, BMI1); either pointer may be NULL
static int fmsimd_cpu( int *popcnt, int *avx2ok )
{
	int pc = 0, lv = FMSIMD_SCALAR, a2 = 0;
#if FMSIMD_X86
	__builtin_cpu_init();
	pc = __builtin_cpu_supports( "popcnt" ) ? 1 : 0;
	a2 = __builtin_cpu_supports( "avx2" ) ? 1 : 0;
	lv = a2 ? FMSIMD_AVX2 : FMSIMD_BASE; // SSE2 is part of x86-64
	if( avx2ok ) *avx2ok = a2 && pc && __builtin_cpu_supports( "bmi" ); // msadist AVX2 kernel: avx2,popcnt,bmi
#elif defined(__aarch64__) && defined(__ARM_NEON)
	lv = FMSIMD_BASE;
	if( avx2ok ) *avx2ok = 0;
#else
	if( avx2ok ) *avx2ok = 0;
#endif
	if( popcnt ) *popcnt = pc;
	(void)a2;
	return( lv );
}

// level to use: the CPU's, lowered by MAFFT_SIMD
static int fmsimd_level( void )
{
	int cpu = fmsimd_cpu( NULL, NULL ), want = cpu;
	const char *e = getenv( "MAFFT_SIMD" );
	if( e && *e )
	{
		if( !strcmp( e, "scalar" ) ) want = FMSIMD_SCALAR;
		else if( !strcmp( e, "sse2" ) || !strcmp( e, "neon" ) || !strcmp( e, "popcnt" ) ) want = FMSIMD_BASE;
		else if( !strcmp( e, "avx2" ) ) want = FMSIMD_AVX2;
		else reporterr( "MAFFT_SIMD=%s not recognised (scalar|sse2|neon|popcnt|avx2)\n", e );
		if( want > cpu )
		{
			reporterr( "MAFFT_SIMD=%s: not available on this CPU, using level %d\n", e, cpu );
			want = cpu;
		}
	}
	return( want );
}
/* END fmsimd dispatch */

/*
 * fast-mafft: exact fast replacement for naivepairscorefast() in the
 * "Making a distance matrix from msa" stage.
 *
 * naivepairscorefast( s1, s2, skip1, skip2, penal ) drops the columns where
 * both rows are gaps and returns
 *     sum over columns with residues in both rows of amino_dis[c1][c2]
 *   + penal * ( number of maximal runs of (gap,residue) columns
 *             + number of maximal runs of (residue,gap) columns )
 * accumulated in a double.  Every summand is an int and every partial sum is
 * an integer far below 2^53, so the double sum is exact; it equals the same
 * sum taken in int64 in any order and converted to double at the end.
 *
 * Gap runs: once gap-gap columns are removed, a run of (gap in i, residue
 * in j) columns is exactly one maximal gap block of row i (in the aligned
 * row) that holds at least one residue of row j.  With G the gap bitmask of
 * row i and X = G & ~G_j, the multiword sum G + X carries out of each such
 * block into the column just after it, a column where G is 0, so the count
 * is popcount( (G + X) & ~G ).  The i->j and j->i result bits never share a
 * column, so one popcount of their OR counts both.
 *
 * Residue scores, bit-plane path: when the alignment has at most 4 residue
 * kinds and some 2-bit coding makes amino_dis[x][y] a function Td of
 * code(x)^code(y) (true for MAFFT's nucleotide matrix: match, transition,
 * transversion), each row is stored as gap mask + 2 code-bit planes and
 *     sum = Td[3]*n_rr + (Td[0]-Td[3])*n_0 + (Td[1]-Td[3])*n_1 + (Td[2]-Td[3])*n_2
 * where n_rr counts both-residue columns and n_d those with code xor d;
 * all counts are popcounts (NEON on arm64; POPCNT or AVX2 on x86-64, see
 * fmsimd above).  With AVX2 the gap runs of row i are counted for 4 rows j at
 * a time in 64-bit lanes (msadist_gaprow_avx2) and the pair loop only adds the
 * residue popcounts.
 * Bit-plane path with exceptions (fast-mafft work/x86simd): with more than 4
 * residue kinds, the 4 most frequent go in the planes if an xor code exists
 * for them and the rest are rare (IUPAC codes in DNA: N, R, Y, ...).  An
 * exception is stored as code 0 and, for every column where one row holds an
 * exception and the other a residue, amino_dis[a][b] - Td[code(a)^code(b)]
 * is added from per-row exception lists.
 * Column path (any other residue bytes < 0x80, e.g. protein): for row i the
 * residue sums with all later rows are accumulated at once by walking, for
 * each residue column of i, the list of later rows with a residue in that
 * column (msadist_colrow); the gap runs come from gap masks as above.  With
 * AVX2, a column where >= 1/8 of the rows have a residue is instead added for
 * 32 later rows at a time from column-major residue codes (vpshufb lookups of
 * amino_dis as int16; msadist_densecol_avx2).  The older generic path
 * (both-residue columns visited bit by bit) remains as a fallback.
 * Rows of unequal length or bytes >= 0x80 fall back to the original code.
 *
 * Debug (environment):
 *   MAFFT_MSADIST=ref      original function for every pair
 *   MAFFT_MSADIST=check    fast and original for every pair; exit(1) on any
 *                          bitwise difference ("generic,check" also works)
 *   MAFFT_MSADIST=stat     print the path chosen (check does too)
 *   MAFFT_MSADIST=generic  skip the bit-plane path
 *   MAFFT_MSADIST=noexc    no exceptions (bit-plane path only for <= 4 kinds)
 *   MAFFT_MSADIST=ctz      generic path instead of the column path
 *   MAFFT_MSADIST=nolane   AVX2 without the 4-lane gap runs
 *   MAFFT_MSADIST=nodense  AVX2 column path without dense columns
 *   MAFFT_MSADIST_DUMP=f   write the finished matrix (row i, j>i, doubles)
 *   MAFFT_MSADIST_SAVEMSA=f  write the aligned rows + amino_dis (bench input)
 */
/* BEGIN msadistfast kernel */
#include <stdint.h>
#if defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#define MSADIST_NEON 1
#else
#define MSADIST_NEON 0
#endif

#define MSADIST_CHUNK 4000 // 128-column blocks per u16 accumulation chunk (16*4000 < 65536)
#define MSADIST_EXC 8      // lut flag: a residue outside the 4 plane kinds (bit-plane path with exceptions)

typedef struct _msadistfast
{
	int nseq;
	int L;        // aligned length
	int W;        // 64-bit words per row mask; even, and W*64 > L so the last carry has room
	int penal;
	int check;    // 1: also run naivepairscorefast() and compare bitwise
	int timing;   // work/v9p, MAFFT_MSADIST=time: thread-seconds of the row passes and the pair loops
	double tcol, tgap, tpair;
	int xorpath;  // 1: bit-plane path, 0: generic path
	int kern;     // kernel: 0 portable C, 1 NEON (arm64) or POPCNT (x86-64), 2 AVX2 (x86-64)
	char **seq;   // the aligned rows (generic path and exceptions read residues from them)
	/* bit-plane path: <= 4 plane kinds and amino_dis[x][y] depends only on code(x)^code(y) */
	uint64_t *planes; // per row, per 128 columns: G[2] b0[2] b1[2] (gap mask, code bit 0, code bit 1)
	long long Td[4];  // amino_dis value for code xor d
	/* exceptions: residues outside the plane kinds (IUPAC codes in DNA); they are
	   stored in the planes as code 0 and the difference is added per column */
	int nexckinds;    // residue kinds outside the plane kinds (0: none)
	int *exn;         // per row: number of exception columns
	uint32_t **excol; // per row: ( column << 8 ) | byte for each exception column, ascending (NULL if none)
	int *CXa, *CXb;   // ( exception byte e << 3 ) | k, p = the plane kind of code k: amino_dis[e][p] - Td[k] (CXa),
	                  // amino_dis[p][e] - Td[k] (CXb); k = 4 (a gap): 0
	int *CB;          // ( a << 7 ) | b, both exceptions: amino_dis[a][b] - Td[0]
	long long nexc;   // total exception columns
	int nexcrows;     // rows with at least one
	/* generic path: any residue bytes < 0x80 */
	uint64_t *masks;  // per row: W words of gap mask G, then W words of residue mask R
	int *T;           // 0x80*0x80 copy of amino_dis (the gap row and column are never used)
	/* column path (generic alphabets): the residue sum of row i with every later
	   row is accumulated at once from per-column lists (msadist_colrow); the gap
	   runs come from the gap masks */
	int colpath;
	uint32_t *centry; // per column, ascending rows: ( row << 8 ) | residue byte
	size_t *coff;     // L+1 offsets into centry
	uint32_t *rpos;   // per row, per residue (ascending column): index of its own entry in centry
	size_t *roff;     // nseq+1 offsets into rpos
	size_t *tcnt;     // preparation: per thread, per column: count, then fill cursor
	/* dense columns (column path, AVX2): a column where >= nseq/8 rows have a residue is
	   added for 32 later rows at a time from column-major residue codes instead of
	   walking its list (msadist_densecol_avx2) */
	int ndense;       // number of dense columns (0: none)
	int *dcol;        // per column: its dense index, or -1
	int np32;         // rows rounded up to a multiple of 32
	uint8_t *dcode;   // ( dense index ) * np32 + row: residue code 0..30, 31 for a gap (and padding rows)
	uint8_t *dtab;    // per residue byte a, 64 bytes: amino_dis[a][kind of code k] as int16 split into
	                  // low bytes (codes 0-15, 16-31) and high bytes (codes 0-15, 16-31); code 31 -> 0
	unsigned char kcode[0x80]; // residue byte -> code
	/* AVX2 (kern 2): the gap runs of row i with 4 rows at a time, in 64-bit lanes
	   (msadist_gaprow_avx2); the gap masks are stored transposed in groups of 4 rows */
	int ne4;          // 1: g4 is built and the gap runs come from msadist_rowprep()
	uint64_t *g4;     // ( ( row >> 2 ) * W + w ) * 4 + ( row & 3 ): gap word w of row
	/* work/v9p: a gap-mask word that is all gaps in both rows adds nothing to the gap runs and
	   passes both carries on unchanged (a = b = ~0: a & ~b = 0, s = a + carry), so the gap-run
	   loops visit only the words where either row (any of the 4 rows of a group) has a residue */
	int nmw;          // 64-bit words of an activity mask: ( W + 63 ) / 64
	uint64_t *ract;   // per row, nmw words: bit w = gap word w of the row is not all gaps
	uint64_t *gact;   // ne4: per group of 4 rows, the OR of their masks
	int noskip;       // MAFFT_MSADIST=noskip: visit every word
	int skip1;        // bit-plane path: Td[1] == Td[3] (DNA), so n1 is not needed
	/* preparation */
	int nthreads;
	int phase;        // 0: lengths and bytes seen, 1: masks, 2: byte counts, 3: exception lists,
	                  // 4: residues per column, 5: column lists, 6: transposed gap masks, 7: dense codes
	unsigned char lut[0x80];
} msadistfast_t;

typedef struct _msadistprep
{
	msadistfast_t *f;
	int t;                 // thread number
	int r0, r1;
	int bad;               // phase 0: a row of another length, or a byte >= 0x80
	unsigned char seen[0x80];
	int nknown;            // phase 0: bytes seen so far, in known[] (NEON fast scan)
	unsigned char known[0x80];
	long long cnt[0x80];   // phase 2: byte counts
} msadistprep_t;

static void msadistfast_free( msadistfast_t *f )
{
	int i;
	if( f == NULL ) return;
	if( f->excol ) for( i=0; i<f->nseq; i++ ) free( f->excol[i] );
	free( f->excol );
	free( f->exn );
	free( f->CXa );
	free( f->CXb );
	free( f->CB );
	free( f->planes );
	free( f->masks );
	free( f->T );
	free( f->centry ); free( f->coff ); free( f->rpos ); free( f->roff ); free( f->tcnt );
	free( f->g4 );
	free( f->ract ); free( f->gact ); // work/v9p
	free( f->dcol ); free( f->dcode ); free( f->dtab );
	free( f );
}

#if MSADIST_NEON
// bit k of the result = top bit of byte k of m0|m1|m2|m3 (each byte 0x00 or 0xff)
static inline uint64_t msadist_movemask64( uint8x16_t m0, uint8x16_t m1, uint8x16_t m2, uint8x16_t m3 )
{
	static const uint8_t wt[16] = { 1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128 };
	const uint8x16_t bw = vld1q_u8( wt );
	uint8x16_t s0 = vpaddq_u8( vandq_u8( m0, bw ), vandq_u8( m1, bw ) );
	uint8x16_t s1 = vpaddq_u8( vandq_u8( m2, bw ), vandq_u8( m3, bw ) );
	s0 = vpaddq_u8( s0, s1 );
	s0 = vpaddq_u8( s0, s0 );
	return( vgetq_lane_u64( vreinterpretq_u64_u8( s0 ), 0 ) );
}

// 64 bytes (all < 0x80) through f->lut
static inline void msadist_lut64( const msadistfast_t *f, const unsigned char *s, uint8x16_t *c0, uint8x16_t *c1, uint8x16_t *c2, uint8x16_t *c3 )
{
	const uint8x16x4_t lo = vld1q_u8_x4( f->lut ), hi = vld1q_u8_x4( f->lut + 64 );
	const uint8x16_t k64 = vdupq_n_u8( 64 );
	uint8x16_t v;
	v = vld1q_u8( s );    *c0 = vqtbx4q_u8( vqtbl4q_u8( lo, v ), hi, vsubq_u8( v, k64 ) );
	v = vld1q_u8( s+16 ); *c1 = vqtbx4q_u8( vqtbl4q_u8( lo, v ), hi, vsubq_u8( v, k64 ) );
	v = vld1q_u8( s+32 ); *c2 = vqtbx4q_u8( vqtbl4q_u8( lo, v ), hi, vsubq_u8( v, k64 ) );
	v = vld1q_u8( s+48 ); *c3 = vqtbx4q_u8( vqtbl4q_u8( lo, v ), hi, vsubq_u8( v, k64 ) );
}
#endif

#if FMSIMD_X86
// bit k = ( s[k] == c ) for 64 bytes (SSE2, part of x86-64)
static inline uint64_t msadist_eqmask64( const unsigned char *s, unsigned char c )
{
	const __m128i v = _mm_set1_epi8( (char)c );
	uint64_t m0 = (uint32_t)_mm_movemask_epi8( _mm_cmpeq_epi8( _mm_loadu_si128( (const __m128i *)s ), v ) );
	uint64_t m1 = (uint32_t)_mm_movemask_epi8( _mm_cmpeq_epi8( _mm_loadu_si128( (const __m128i *)( s + 16 ) ), v ) );
	uint64_t m2 = (uint32_t)_mm_movemask_epi8( _mm_cmpeq_epi8( _mm_loadu_si128( (const __m128i *)( s + 32 ) ), v ) );
	uint64_t m3 = (uint32_t)_mm_movemask_epi8( _mm_cmpeq_epi8( _mm_loadu_si128( (const __m128i *)( s + 48 ) ), v ) );
	return( m0 | ( m1 << 16 ) | ( m2 << 32 ) | ( m3 << 48 ) );
}
#endif

// phase 0: check lengths and collect the residue bytes; phase 1: build the row bit masks;
// phase 2: count the bytes; phase 3: list the exception columns of each row
static void *msadist_prepthread( void *arg )
{
	msadistprep_t *p = (msadistprep_t *)arg;
	msadistfast_t *f = p->f;
	int i, c, k, w, cmax, L = f->L, W = f->W;
	uint32_t *buf = NULL;
	if( f->phase == 3 ) buf = malloc( ( L + 1 ) * sizeof( uint32_t ) );
	for( i=p->r0; i<p->r1; i++ )
	{
		const unsigned char *s = (const unsigned char *)f->seq[i];
		if( f->phase == 0 )
		{
			unsigned char hi = 0;
			if( (int)strlen( f->seq[i] ) != L ) { p->bad = 1; continue; }
			c = 0;
#if MSADIST_NEON
			// 16 bytes at a time while they are all bytes already seen (<= 8 kinds)
			while( p->nknown <= 8 && c + 16 <= L )
			{
				uint8x16_t v = vld1q_u8( s+c ), m = vdupq_n_u8( 0 );
				for( k=0; k<p->nknown; k++ ) m = vorrq_u8( m, vceqq_u8( v, vdupq_n_u8( p->known[k] ) ) );
				if( vminvq_u8( m ) != 0xff )
				{
					for( k=c; k<c+16; k++ )
					{
						hi |= s[k];
						if( !p->seen[s[k]&0x7f] ) { p->seen[s[k]&0x7f] = 1; p->known[p->nknown++] = s[k]&0x7f; }
					}
				}
				c += 16;
			}
#elif FMSIMD_X86
			// the same with SSE2
			while( p->nknown <= 8 && c + 16 <= L )
			{
				__m128i v = _mm_loadu_si128( (const __m128i *)( s+c ) ), m = _mm_setzero_si128();
				for( k=0; k<p->nknown; k++ ) m = _mm_or_si128( m, _mm_cmpeq_epi8( v, _mm_set1_epi8( (char)p->known[k] ) ) );
				if( _mm_movemask_epi8( m ) != 0xffff )
				{
					for( k=c; k<c+16; k++ )
					{
						hi |= s[k];
						if( !p->seen[s[k]&0x7f] ) { p->seen[s[k]&0x7f] = 1; p->known[p->nknown++] = s[k]&0x7f; }
					}
				}
				c += 16;
			}
#endif
			for( ; c<L; c++ ) { hi |= s[c]; p->seen[s[c]&0x7f] = 1; }
			if( hi & 0x80 ) p->bad = 1;
		}
		else if( f->phase == 2 )
		{
			for( c=0; c<L; c++ ) p->cnt[s[c]]++;
		}
		else if( f->phase == 7 )
		{
			for( c=0; c<L; c++ ) if( f->dcol[c] >= 0 ) f->dcode[(size_t)f->dcol[c] * f->np32 + i] = f->kcode[s[c]];
		}
		else if( f->phase == 8 ) // work/v9p: activity masks of the gap words
		{
			uint64_t *ra = f->ract + (size_t)i * f->nmw;
			for( w=0; w<W; w++ )
			{
				uint64_t g = f->xorpath ? f->planes[(size_t)i * 3 * W + ( w >> 1 ) * 6 + ( w & 1 )] : f->masks[(size_t)i * 2 * W + w];
				if( ~g ) ra[w>>6] |= 1ULL << ( w & 63 );
			}
		}
		else if( f->phase == 6 )
		{
			uint64_t *g4 = f->g4 + (size_t)( i >> 2 ) * W * 4 + ( i & 3 );
			if( f->xorpath ) for( w=0; w<W; w++ ) g4[4*w] = f->planes[(size_t)i * 3 * W + ( w >> 1 ) * 6 + ( w & 1 )];
			else for( w=0; w<W; w++ ) g4[4*w] = f->masks[(size_t)i * 2 * W + w];
		}
		else if( f->phase == 4 || f->phase == 5 )
		{
			// residues of row i in column order, from its residue mask
			const uint64_t *mk = f->masks + (size_t)i * 2 * W;
			size_t *tc = f->tcnt + (size_t)p->t * L;
			uint32_t *rp = f->rpos ? f->rpos + f->roff[i] : NULL;
			int m = 0;
			for( w=0; w<W; w++ )
			{
				uint64_t r = mk[W+w];
				while( r )
				{
					c = w*64 + __builtin_ctzll( r );
					if( f->phase == 4 ) tc[c]++;
					else
					{
						size_t q = tc[c]++;
						f->centry[q] = ( (uint32_t)i << 8 ) | s[c];
						rp[m] = (uint32_t)q;
					}
					m++;
					r &= r - 1;
				}
			}
			if( f->phase == 4 ) f->roff[i+1] = m;
		}
		else if( f->phase == 3 )
		{
			int m = 0;
			for( c=0; c<L; c++ ) if( f->lut[s[c]] & MSADIST_EXC ) buf[m++] = ( (uint32_t)c << 8 ) | s[c];
			f->exn[i] = m;
			if( m )
			{
				f->excol[i] = malloc( m * sizeof( uint32_t ) );
				memcpy( f->excol[i], buf, m * sizeof( uint32_t ) );
			}
		}
		else if( f->xorpath )
		{
			uint64_t *pl = f->planes + (size_t)i * 3 * W;
			for( w=0; w*64<L; w++ )
			{
				uint64_t g = 0, b0 = 0, b1 = 0, v;
				uint64_t *blk = pl + ( w >> 1 ) * 6 + ( w & 1 );
				cmax = L - w*64 < 64 ? L - w*64 : 64;
#if MSADIST_NEON
				if( cmax == 64 )
				{
					uint8x16_t c0, c1, c2, c3;
					msadist_lut64( f, s + w*64, &c0, &c1, &c2, &c3 );
					blk[0] = msadist_movemask64( vtstq_u8( c0, vdupq_n_u8( 4 ) ), vtstq_u8( c1, vdupq_n_u8( 4 ) ), vtstq_u8( c2, vdupq_n_u8( 4 ) ), vtstq_u8( c3, vdupq_n_u8( 4 ) ) );
					blk[2] = msadist_movemask64( vtstq_u8( c0, vdupq_n_u8( 1 ) ), vtstq_u8( c1, vdupq_n_u8( 1 ) ), vtstq_u8( c2, vdupq_n_u8( 1 ) ), vtstq_u8( c3, vdupq_n_u8( 1 ) ) );
					blk[4] = msadist_movemask64( vtstq_u8( c0, vdupq_n_u8( 2 ) ), vtstq_u8( c1, vdupq_n_u8( 2 ) ), vtstq_u8( c2, vdupq_n_u8( 2 ) ), vtstq_u8( c3, vdupq_n_u8( 2 ) ) );
					continue;
				}
#endif
				for( k=0; k<cmax; k++ )
				{
					v = f->lut[s[w*64+k]];
					g |= ( ( v >> 2 ) & 1 ) << k;
					b0 |= ( v & 1 ) << k;
					b1 |= ( ( v >> 1 ) & 1 ) << k;
				}
				blk[0] = g; blk[2] = b0; blk[4] = b1;
			}
		}
		else
		{
			uint64_t *mk = f->masks + (size_t)i * 2 * W;
			for( w=0; w*64<L; w++ )
			{
				uint64_t g = 0, r = 0, v;
				cmax = L - w*64 < 64 ? L - w*64 : 64;
#if MSADIST_NEON
				if( cmax == 64 )
				{
					const uint8x16_t gap = vdupq_n_u8( '-' );
					g = msadist_movemask64( vceqq_u8( vld1q_u8( s+w*64 ), gap ), vceqq_u8( vld1q_u8( s+w*64+16 ), gap ), vceqq_u8( vld1q_u8( s+w*64+32 ), gap ), vceqq_u8( vld1q_u8( s+w*64+48 ), gap ) );
					mk[w] = g; mk[W+w] = ~g;
					continue;
				}
#elif FMSIMD_X86
				if( cmax == 64 )
				{
					g = msadist_eqmask64( s+w*64, '-' );
					mk[w] = g; mk[W+w] = ~g;
					continue;
				}
#endif
				for( k=0; k<cmax; k++ )
				{
					v = ( s[w*64+k] == '-' );
					g |= v << k;
					r |= ( v ^ 1 ) << k;
				}
				mk[w] = g; mk[W+w] = r;
			}
		}
	}
	free( buf );
	return( NULL );
}

static void msadist_prepphase( msadistfast_t *f, msadistprep_t *p, int phase )
{
	int t, nth = f->nthreads;
	f->phase = phase;
	for( t=0; t<nth; t++ )
	{
		p[t].f = f;
		p[t].t = t;
		p[t].r0 = (int)( (long long)f->nseq * t / nth );
		p[t].r1 = (int)( (long long)f->nseq * ( t + 1 ) / nth );
	}
#ifdef enablemultithread
	if( nth > 1 )
	{
		pthread_t *handle = calloc( nth, sizeof( pthread_t ) );
		for( t=0; t<nth; t++ ) pthread_create( handle+t, NULL, msadist_prepthread, (void *)(p+t) );
		for( t=0; t<nth; t++ ) pthread_join( handle[t], NULL );
		free( handle );
	}
	else
#endif
		for( t=0; t<nth; t++ ) msadist_prepthread( (void *)(p+t) );
}

// find a 2-bit code for each of the K<=4 residue kinds such that
// amino_dis[chr[a]][chr[b]] is a function of code[a]^code[b]
static int msadist_findxorcode( int K, const unsigned char *chr, int *code, long long *Td )
{
	int p[4], a, b, d, ok, set[4];
	long long v, f[4];
	for( p[0]=0; p[0]<4; p[0]++ ) for( p[1]=0; p[1]<4; p[1]++ ) for( p[2]=0; p[2]<4; p[2]++ ) for( p[3]=0; p[3]<4; p[3]++ )
	{
		if( p[0]==p[1] || p[0]==p[2] || p[0]==p[3] || p[1]==p[2] || p[1]==p[3] || p[2]==p[3] ) continue;
		ok = 1;
		for( d=0; d<4; d++ ) { set[d] = 0; f[d] = 0; }
		for( a=0; a<K && ok; a++ ) for( b=0; b<K; b++ )
		{
			d = p[a] ^ p[b];
			v = amino_dis[chr[a]][chr[b]];
			if( !set[d] ) { set[d] = 1; f[d] = v; }
			else if( f[d] != v ) { ok = 0; break; }
		}
		if( ok )
		{
			for( a=0; a<4; a++ ) code[a] = p[a];
			for( d=0; d<4; d++ ) Td[d] = f[d]; // an unset entry never occurs (no such pair of kinds)
			return( 1 );
		}
	}
	return( 0 );
}

static const char *msadist_kname( const msadistfast_t *f )
{
	if( f->kern == 2 ) return( "avx2" );
#if MSADIST_NEON
	if( f->kern == 1 ) return( "neon" );
#endif
	if( f->kern == 1 ) return( "popcnt" );
	return( "portable" );
}

// Returns NULL when the original function must be used.
// MAFFT_MSADIST may contain "ref" (use the original), "check" (compare every pair
// with the original), "generic" (skip the bit-plane path), "noexc" (bit-plane
// path only for <= 4 residue kinds) and "stat" (print the path chosen).
static msadistfast_t *msadistfast_prepare( int nseq, char **seq, int penal, int nthreads )
{
	msadistfast_t *f;
	msadistprep_t *p;
	int code[4];
	unsigned char chr[4], all[0x80];
	int t, a, b, c, K, Kp, bad, lv, pc, a2;
	long long Td[4], cnt[0x80], nres;
	unsigned char seen[0x80];
	const char *mode = getenv( "MAFFT_MSADIST" );
	int stat;

	if( mode && strstr( mode, "ref" ) ) return( NULL );
	if( nseq < 2 ) return( NULL );
	if( nthreads < 1 ) nthreads = 1;
	if( nthreads > nseq ) nthreads = nseq;
	stat = ( mode && ( strstr( mode, "stat" ) || strstr( mode, "check" ) ) );

	f = calloc( 1, sizeof( msadistfast_t ) );
	p = calloc( nthreads, sizeof( msadistprep_t ) );
	f->nseq = nseq;
	f->seq = seq;
	f->L = strlen( seq[0] );
	f->W = ( f->L / 64 + 2 ) & ~1;
	f->penal = penal;
	f->check = ( mode && strstr( mode, "check" ) );
	f->timing = ( mode && strstr( mode, "time" ) );
	f->nthreads = nthreads;

	lv = fmsimd_level();
	fmsimd_cpu( &pc, &a2 );
	f->kern = 0;
#if MSADIST_NEON
	if( lv >= FMSIMD_BASE ) f->kern = 1;
#elif FMSIMD_X86
	if( lv >= FMSIMD_AVX2 && a2 ) f->kern = 2;
	else if( lv >= FMSIMD_BASE && pc ) f->kern = 1;
#endif
	(void)pc; (void)a2;

	msadist_prepphase( f, p, 0 );
	bad = 0;
	memset( seen, 0, sizeof( seen ) );
	for( t=0; t<nthreads; t++ )
	{
		bad |= p[t].bad;
		for( c=0; c<0x80; c++ ) seen[c] |= p[t].seen[c];
	}
	if( bad || seen[0] )
	{
		if( stat ) reporterr( "msadist fast: not used (rows of unequal length or bytes >= 0x80): original code for every pair\n" );
		free( p );
		msadistfast_free( f );
		return( NULL );
	}

	K = 0;
	for( c=1; c<0x80; c++ ) if( seen[c] && c != '-' )
	{
		if( K < 4 ) chr[K] = c;
		all[K] = c;
		K++;
	}
	Kp = K; // kinds in the planes
	f->xorpath = 0;
	if( !( mode && strstr( mode, "generic" ) ) )
	{
		if( K <= 4 ) f->xorpath = msadist_findxorcode( K, chr, code, Td );
		else if( !( mode && strstr( mode, "noexc" ) ) )
		{
			// more than 4 kinds: the 4 most frequent go in the planes if an xor code
			// exists for them and the others are rare (<= L/32 per row on average)
			msadist_prepphase( f, p, 2 );
			memset( cnt, 0, sizeof( cnt ) );
			for( t=0; t<nthreads; t++ ) for( c=0; c<0x80; c++ ) cnt[c] += p[t].cnt[c];
			for( a=0; a<4; a++ ) // selection of the top 4 (ties: smaller byte first)
			{
				int best = -1;
				for( b=0; b<K; b++ )
				{
					int used = 0;
					for( c=0; c<a; c++ ) if( chr[c] == all[b] ) used = 1;
					if( !used && ( best < 0 || cnt[all[b]] > cnt[all[best]] ) ) best = b;
				}
				chr[a] = all[best];
			}
			nres = 0;
			for( b=0; b<K; b++ ) nres += cnt[all[b]];
			f->nexc = nres - cnt[chr[0]] - cnt[chr[1]] - cnt[chr[2]] - cnt[chr[3]];
			Kp = 4;
			if( f->nexc * 32 <= (long long)nseq * f->L && msadist_findxorcode( 4, chr, code, Td ) )
			{
				f->xorpath = 1;
				f->nexckinds = K - 4;
			}
			else if( stat ) reporterr( "msadist fast: %d residue kinds, %lld of %lld residues outside the top 4; generic path\n", K, f->nexc, nres );
		}
	}
	if( f->xorpath )
	{
		for( a=0; a<4; a++ ) f->Td[a] = Td[a];
		memset( f->lut, 0, sizeof( f->lut ) );
		f->lut['-'] = 4;
		if( f->nexckinds ) for( b=0; b<K; b++ ) f->lut[all[b]] = MSADIST_EXC; // code 0 in the planes
		for( a=0; a<Kp; a++ ) f->lut[chr[a]] = code[a];
		f->planes = calloc( (size_t)nseq * 3 * f->W, sizeof( uint64_t ) );
		if( f->nexckinds )
		{
			f->exn = calloc( nseq, sizeof( int ) );
			f->excol = calloc( nseq, sizeof( uint32_t * ) );
			f->CXa = calloc( 0x80 * 8, sizeof( int ) );
			f->CXb = calloc( 0x80 * 8, sizeof( int ) );
			f->CB = calloc( 0x80 * 0x80, sizeof( int ) );
			// an exception is code 0 in the planes, so the plane sum counted Td[0^k] = Td[k]
			// against the plane kind chr[b] of code k = code[b], and Td[0] against another exception
			for( a=1; a<0x80; a++ ) if( f->lut[a] & MSADIST_EXC )
			{
				for( b=0; b<4; b++ )
				{
					f->CXa[(a<<3)|code[b]] = amino_dis[a][chr[b]] - (int)Td[code[b]];
					f->CXb[(a<<3)|code[b]] = amino_dis[chr[b]][a] - (int)Td[code[b]];
				}
				for( b=1; b<0x80; b++ ) if( f->lut[b] & MSADIST_EXC ) f->CB[(a<<7)|b] = amino_dis[a][b] - (int)Td[0];
			}
			msadist_prepphase( f, p, 3 );
			f->nexc = 0;
			for( a=0; a<nseq; a++ )
			{
				f->nexc += f->exn[a];
				if( f->exn[a] ) f->nexcrows++;
			}
		}
		if( stat ) reporterr( "msadist fast (bit-plane path, %s%s): L=%d K=%d Td=%lld,%lld,%lld,%lld penal=%d; exceptions: %d kinds, %lld columns in %d of %d rows\n", msadist_kname( f ), ( f->kern == 2 && !( mode && strstr( mode, "nolane" ) ) ) ? ", 4-lane gap runs" : "", f->L, K, Td[0], Td[1], Td[2], Td[3], penal, f->nexckinds, f->nexc, f->nexcrows, nseq );
	}
	else
	{
		f->T = calloc( 0x80 * 0x80, sizeof( int ) );
		for( a=0; a<0x80; a++ ) for( b=0; b<0x80; b++ ) f->T[(a<<7)|b] = ( a == '-' || b == '-' ) ? 0 : amino_dis[a][b];
		f->masks = calloc( (size_t)nseq * 2 * f->W, sizeof( uint64_t ) );
		// column path unless switched off or an int sum could overflow
		{
			long long tmax = 0;
			for( a=0; a<0x80*0x80; a++ ) if( llabs( (long long)f->T[a] ) > tmax ) tmax = llabs( (long long)f->T[a] );
			f->colpath = !( mode && strstr( mode, "ctz" ) ) && nseq < ( 1 << 24 ) && tmax * f->L < 0x7fffffffLL;
		}
	}
	msadist_prepphase( f, p, 1 );
	if( f->colpath )
	{
		size_t tot, q;
		int L = f->L;
		f->tcnt = calloc( (size_t)nthreads * L, sizeof( size_t ) );
		f->roff = calloc( nseq + 1, sizeof( size_t ) );
		f->coff = calloc( L + 1, sizeof( size_t ) );
		msadist_prepphase( f, p, 4 );
		for( a=0; a<nseq; a++ ) f->roff[a+1] += f->roff[a];
		for( c=0; c<L; c++ ) for( t=0; t<nthreads; t++ ) f->coff[c+1] += f->tcnt[(size_t)t*L+c];
		for( c=0; c<L; c++ ) f->coff[c+1] += f->coff[c];
		tot = f->coff[L];
		if( tot < 0xffffffffUL ) // entry indices are 32-bit
		{
			for( c=0; c<L; c++ ) // fill cursors: thread t's rows follow those of threads < t
			{
				q = f->coff[c];
				for( t=0; t<nthreads; t++ ) { size_t m = f->tcnt[(size_t)t*L+c]; f->tcnt[(size_t)t*L+c] = q; q += m; }
			}
			f->centry = malloc( ( tot + 1 ) * sizeof( uint32_t ) );
			f->rpos = malloc( ( tot + 1 ) * sizeof( uint32_t ) );
		}
		if( f->centry && f->rpos ) msadist_prepphase( f, p, 5 );
		else // too many residues or no memory: the generic path
		{
			free( f->centry ); free( f->rpos ); free( f->coff ); free( f->roff );
			f->centry = f->rpos = NULL; f->coff = f->roff = NULL;
			f->colpath = 0;
		}
		free( f->tcnt ); f->tcnt = NULL;
	}
	f->skip1 = f->xorpath && f->Td[1] == f->Td[3];
#if FMSIMD_X86
	if( f->colpath && f->kern == 2 && !( mode && strstr( mode, "nodense" ) ) )
	{
		int K2 = 0, nd = 0, ok = 1, kb[31];
		memset( f->kcode, 31, sizeof( f->kcode ) );
		for( c=1; c<0x80; c++ ) if( seen[c] && c != '-' )
		{
			if( K2 == 31 ) { ok = 0; break; }
			kb[K2] = c;
			f->kcode[c] = K2++;
		}
		for( a=0; a<K2 && ok; a++ ) for( b=0; b<K2; b++ ) if( amino_dis[kb[a]][kb[b]] < -32768 || amino_dis[kb[a]][kb[b]] > 32767 ) ok = 0;
		if( ok )
		{
			f->dcol = malloc( f->L * sizeof( int ) );
			for( c=0; c<f->L; c++ )
			{
				f->dcol[c] = -1;
				if( (long long)( f->coff[c+1] - f->coff[c] ) * 8 >= nseq && nseq >= 64 ) f->dcol[c] = nd++;
			}
			f->np32 = ( nseq + 31 ) & ~31;
			if( nd && !posix_memalign( (void **)&f->dcode, 32, (size_t)nd * f->np32 ) && ( f->dtab = calloc( 0x80, 64 ) ) )
			{
				memset( f->dcode, 31, (size_t)nd * f->np32 );
				for( a=0; a<K2; a++ ) for( b=0; b<K2; b++ )
				{
					int v = amino_dis[kb[a]][kb[b]];
					uint8_t *tb = f->dtab + kb[a] * 64;
					tb[( b >> 4 ) * 16 + ( b & 15 )] = (uint8_t)( v & 0xff );
					tb[32 + ( b >> 4 ) * 16 + ( b & 15 )] = (uint8_t)( ( v >> 8 ) & 0xff );
				}
				msadist_prepphase( f, p, 7 );
				f->ndense = nd;
			}
			else
			{
				free( f->dcol ); free( f->dcode ); free( f->dtab );
				f->dcol = NULL; f->dcode = NULL; f->dtab = NULL;
			}
		}
	}

	if( f->kern == 2 && ( f->xorpath || f->colpath ) && !( mode && strstr( mode, "nolane" ) ) )
	{
		size_t ng = ( (size_t)nseq + 3 ) / 4;
		if( !posix_memalign( (void **)&f->g4, 32, ng * f->W * 4 * sizeof( uint64_t ) ) )
		{
			memset( f->g4, 0, ng * f->W * 4 * sizeof( uint64_t ) );
			msadist_prepphase( f, p, 6 );
			f->ne4 = 1;
		}
		else f->g4 = NULL;
	}
#endif
	// work/v9p: activity masks for skipping the words that are all gaps in both rows
	f->noskip = ( mode && strstr( mode, "noskip" ) );
	f->nmw = ( f->W + 63 ) / 64;
	if( !f->noskip && ( f->ne4 || f->colpath ) )
	{
		f->ract = calloc( (size_t)nseq * f->nmw, sizeof( uint64_t ) );
		if( f->ract ) msadist_prepphase( f, p, 8 );
		if( f->ract && f->ne4 && ( f->gact = calloc( (size_t)( ( nseq + 3 ) / 4 ) * f->nmw, sizeof( uint64_t ) ) ) )
		{
			for( a=0; a<nseq; a++ ) for( c=0; c<f->nmw; c++ ) f->gact[(size_t)( a >> 2 ) * f->nmw + c] |= f->ract[(size_t)a * f->nmw + c];
		}
		if( f->ne4 && !f->gact ) { free( f->ract ); f->ract = NULL; } // without the group masks: no skipping
		// skipping pays only when many words can be skipped (the bit loop costs a few operations per
		// visited word): groups of 4 (AVX2) or single rows (a pair visits the union of two) below 3/4 / 1/2
		if( f->ract && !( mode && strstr( mode, "forceskip" ) ) )
		{
			long long nr = 0, ng = 0;
			if( f->gact ) for( a=0; a<(nseq+3)/4; a++ ) for( c=0; c<f->nmw; c++ ) ng += __builtin_popcountll( f->gact[(size_t)a*f->nmw+c] );
			else for( a=0; a<nseq; a++ ) for( c=0; c<f->nmw; c++ ) nr += __builtin_popcountll( f->ract[(size_t)a*f->nmw+c] );
			if( f->gact ? ng * 4 > 3LL * ( ( nseq + 3 ) / 4 ) * f->W : nr * 2 > (long long)nseq * f->W )
			{
				free( f->ract ); free( f->gact );
				f->ract = f->gact = NULL;
			}
		}
	}
	if( !f->xorpath && stat ) reporterr( "msadist fast (%s path, %s%s): L=%d K=%d penal=%d; %d dense columns\n", f->colpath ? "column" : "generic", msadist_kname( f ), f->ne4 ? ", 4-lane gap runs" : "", f->L, K, penal, f->ndense );
	if( stat && f->ract ) // work/v9p
	{
		long long nr = 0, ng = 0;
		for( a=0; a<nseq; a++ ) for( c=0; c<f->nmw; c++ ) nr += __builtin_popcountll( f->ract[(size_t)a*f->nmw+c] );
		if( f->gact ) for( a=0; a<(nseq+3)/4; a++ ) for( c=0; c<f->nmw; c++ ) ng += __builtin_popcountll( f->gact[(size_t)a*f->nmw+c] );
		reporterr( "msadist fast: gap-run words with a residue: %.1f of %d per row%s%.1f per group of 4\n", (double)nr / nseq, f->W, f->gact ? ", " : "", f->gact ? (double)ng / ( ( nseq + 3 ) / 4 ) : 0.0 );
	}
	free( p );
	return( f );
}

// one 64-column word of the gap-run count: blocks of a (gaps of one row) that
// contain a residue of the other row end where the carry of a + (a & ~b) lands
#define MSADIST_RUNWORD( a, b, ca, cb, e ) \
	do { \
		uint64_t a_ = (a), b_ = (b), ta_, tb_, sa_, sb_; \
		ta_ = a_ + ( a_ & ~b_ ); sa_ = ta_ + (ca); \
		tb_ = b_ + ( b_ & ~a_ ); sb_ = tb_ + (cb); \
		(ca) = ( ta_ < a_ ) | ( sa_ < ta_ ); \
		(cb) = ( tb_ < b_ ) | ( sb_ < tb_ ); \
		(e) = ( sa_ & ~a_ ) | ( sb_ & ~b_ ); \
	} while( 0 )

// nA counts columns with a gap in either row, nB those that are not identical
// residue pairs, n1/n2 the both-residue columns with code xor 1/2, nE the gap
// runs.  Padding columns (>= len) have no gap bit and code 0 in both rows, so
// they add to neither nA nor nB.  Both-residue columns: len - nA; with code
// xor 0: len - nB; xor 1: n1; xor 2: n2; xor 3: the rest.
static inline long long msadist_xorsum( const msadistfast_t *f, long long nA, long long nB, long long n1, long long n2, long long nE )
{
	long long len = f->L;
	return( f->Td[3] * ( len - nA ) + ( f->Td[0] - f->Td[3] ) * ( len - nB ) + ( f->Td[1] - f->Td[3] ) * n1 + ( f->Td[2] - f->Td[3] ) * n2
	      + (long long)f->penal * nE );
}

// portable C
static long long msadist_xorscore_c( const msadistfast_t *f, const uint64_t *pi, const uint64_t *pj )
{
	long long nA = 0, nB = 0, n1 = 0, n2 = 0, nE = 0;
	int nblk = f->W / 2, k, h;
	uint64_t ci = 0, cj = 0, e0, e1;
	for( k=0; k<nblk; k++ )
	{
		const uint64_t *qi = pi + 6*k, *qj = pj + 6*k;
		for( h=0; h<2; h++ )
		{
			uint64_t gor = qi[h] | qj[h], x0 = qi[2+h] ^ qj[2+h], x1 = qi[4+h] ^ qj[4+h];
			nA += __builtin_popcountll( gor );
			nB += __builtin_popcountll( gor | x0 | x1 );
			n1 += __builtin_popcountll( x0 & ~( gor | x1 ) );
			n2 += __builtin_popcountll( x1 & ~( gor | x0 ) );
		}
		MSADIST_RUNWORD( qi[0], qj[0], ci, cj, e0 );
		MSADIST_RUNWORD( qi[1], qj[1], ci, cj, e1 );
		nE += __builtin_popcountll( e0 ) + __builtin_popcountll( e1 );
	}
	return( msadist_xorsum( f, nA, nB, n1, n2, nE ) );
}

#if MSADIST_NEON
static long long msadist_xorscore_neon( const msadistfast_t *f, const uint64_t *pi, const uint64_t *pj )
{
	long long nA = 0, nB = 0, n1 = 0, n2 = 0, nE = 0;
	int nblk = f->W / 2, k0, k1, k;
	uint64_t ci = 0, cj = 0, e0, e1;
	for( k0=0; k0<nblk; k0=k1 )
	{
		uint16x8_t aA = vdupq_n_u16( 0 ), aB = vdupq_n_u16( 0 ), a1 = vdupq_n_u16( 0 ), a2 = vdupq_n_u16( 0 ), aE = vdupq_n_u16( 0 );
		k1 = k0 + MSADIST_CHUNK < nblk ? k0 + MSADIST_CHUNK : nblk;
		for( k=k0; k<k1; k++ )
		{
			const uint64_t *qi = pi + 6*k, *qj = pj + 6*k;
			uint64x2_t gi = vld1q_u64( qi ), i0 = vld1q_u64( qi+2 ), i1 = vld1q_u64( qi+4 );
			uint64x2_t gj = vld1q_u64( qj ), j0 = vld1q_u64( qj+2 ), j1 = vld1q_u64( qj+4 );
			uint64x2_t gor = vorrq_u64( gi, gj );
			uint64x2_t x0 = veorq_u64( i0, j0 ), x1 = veorq_u64( i1, j1 );
			uint64x2_t t0 = vorrq_u64( gor, x0 ), t1 = vorrq_u64( gor, x1 );
			aA = vpadalq_u8( aA, vcntq_u8( vreinterpretq_u8_u64( gor ) ) );
			aB = vpadalq_u8( aB, vcntq_u8( vreinterpretq_u8_u64( vorrq_u64( t0, x1 ) ) ) );
			a1 = vpadalq_u8( a1, vcntq_u8( vreinterpretq_u8_u64( vbicq_u64( x0, t1 ) ) ) );
			a2 = vpadalq_u8( a2, vcntq_u8( vreinterpretq_u8_u64( vbicq_u64( x1, t0 ) ) ) );
			MSADIST_RUNWORD( qi[0], qj[0], ci, cj, e0 );
			MSADIST_RUNWORD( qi[1], qj[1], ci, cj, e1 );
			aE = vpadalq_u8( aE, vcntq_u8( vreinterpretq_u8_u64( vcombine_u64( vcreate_u64( e0 ), vcreate_u64( e1 ) ) ) ) );
		}
		nA += vaddlvq_u16( aA ); nB += vaddlvq_u16( aB ); n1 += vaddlvq_u16( a1 ); n2 += vaddlvq_u16( a2 ); nE += vaddlvq_u16( aE );
	}
	return( msadist_xorsum( f, nA, nB, n1, n2, nE ) );
}
#endif

#if FMSIMD_X86
// The gap-run words of one 128-column block with x86 add-with-carry (adc):
// a0,a1 the gap words of row i, b0,b1 of row j; ca, cb the carries of the
// two chains (unsigned char).  Same result as two MSADIST_RUNWORDs.
#define MSADIST_ADC2( a0, a1, b0, b1, ca, cb, e0, e1 ) \
	do { \
		unsigned long long ga0_ = (a0), ga1_ = (a1), gb0_ = (b0), gb1_ = (b1), sa0_, sa1_, sb0_, sb1_; \
		unsigned long long xa0_ = ga0_ & ~gb0_, xa1_ = ga1_ & ~gb1_, xb0_ = gb0_ & ~ga0_, xb1_ = gb1_ & ~ga1_; \
		(ca) = _addcarry_u64( (ca), ga0_, xa0_, &sa0_ ); \
		(ca) = _addcarry_u64( (ca), ga1_, xa1_, &sa1_ ); \
		(cb) = _addcarry_u64( (cb), gb0_, xb0_, &sb0_ ); \
		(cb) = _addcarry_u64( (cb), gb1_, xb1_, &sb1_ ); \
		(e0) = ( sa0_ & ~ga0_ ) | ( sb0_ & ~gb0_ ); \
		(e1) = ( sa1_ & ~ga1_ ) | ( sb1_ & ~gb1_ ); \
	} while( 0 )

// one 128-column block with scalar popcounts (inside a target("popcnt") function)
#define MSADIST_POPBLOCK( qi, qj ) \
	do { \
		int h_; \
		uint64_t e0_, e1_; \
		for( h_=0; h_<2; h_++ ) \
		{ \
			uint64_t gor_ = (qi)[h_] | (qj)[h_], x0_ = (qi)[2+h_] ^ (qj)[2+h_], x1_ = (qi)[4+h_] ^ (qj)[4+h_]; \
			uint64_t t0_ = gor_ | x0_, t1_ = gor_ | x1_; \
			nA += __builtin_popcountll( gor_ ); \
			nB += __builtin_popcountll( t0_ | x1_ ); \
			n1 += __builtin_popcountll( x0_ & ~t1_ ); \
			n2 += __builtin_popcountll( x1_ & ~t0_ ); \
		} \
		MSADIST_ADC2( (qi)[0], (qi)[1], (qj)[0], (qj)[1], ci, cj, e0_, e1_ ); \
		nE += __builtin_popcountll( e0_ ) + __builtin_popcountll( e1_ ); \
	} while( 0 )

__attribute__((target("popcnt")))
static long long msadist_xorscore_popcnt( const msadistfast_t *f, const uint64_t *pi, const uint64_t *pj )
{
	long long nA = 0, nB = 0, n1 = 0, n2 = 0, nE = 0;
	int nblk = f->W / 2, k;
	unsigned char ci = 0, cj = 0;
	for( k=0; k<nblk; k++ ) MSADIST_POPBLOCK( pi + 6*k, pj + 6*k );
	return( msadist_xorsum( f, nA, nB, n1, n2, nE ) );
}

// byte-wise popcount of a 256-bit vector (nibble table, pshufb)
#define MSADIST_POPB( v ) _mm256_add_epi8( _mm256_shuffle_epi8( lut4, _mm256_and_si256( (v), m4 ) ), _mm256_shuffle_epi8( lut4, _mm256_and_si256( _mm256_srli_epi16( (v), 4 ), m4 ) ) )
#define MSADIST_HSUM( v, out ) do { __m128i h_ = _mm_add_epi64( _mm256_castsi256_si128( v ), _mm256_extracti128_si256( v, 1 ) ); (out) = (long long)_mm_cvtsi128_si64( h_ ) + (long long)_mm_cvtsi128_si64( _mm_unpackhi_epi64( h_, h_ ) ); } while( 0 )

// AVX2: two 128-column blocks per step.  The four residue popcounts run in
// 256-bit vectors (byte counts, widened with vpsadbw every 31 steps); the gap-run
// chain runs in scalar registers with adc and popcnt alongside.
__attribute__((target("avx2,popcnt,bmi")))
static long long msadist_xorscore_avx2( const msadistfast_t *f, const uint64_t *pi, const uint64_t *pj )
{
	const __m256i lut4 = _mm256_setr_epi8( 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 );
	const __m256i m4 = _mm256_set1_epi8( 0x0f ), z = _mm256_setzero_si256();
	__m256i sA = z, sB = z, s1 = z, s2 = z;
	long long nA = 0, nB = 0, n1 = 0, n2 = 0, nE = 0;
	int nblk = f->W / 2, npr = nblk / 2, k, k0, k1;
	unsigned char ci = 0, cj = 0;
	for( k0=0; k0<npr; k0=k1 )
	{
		__m256i bA = z, bB = z, b1 = z, b2 = z; // byte counts, at most 8 per step
		k1 = k0 + 31 < npr ? k0 + 31 : npr;
		for( k=k0; k<k1; k++ )
		{
			const uint64_t *qi = pi + 12*k, *qj = pj + 12*k;
			uint64_t e0, e1, e2, e3;
			// two blocks = 12 words, as 3 vectors: [G b0] [b1 G'] [b0' b1'] (128 bits each)
			__m256i ai = _mm256_loadu_si256( (const __m256i *)qi ), bi = _mm256_loadu_si256( (const __m256i *)( qi + 4 ) ), di = _mm256_loadu_si256( (const __m256i *)( qi + 8 ) );
			__m256i aj = _mm256_loadu_si256( (const __m256i *)qj ), bj = _mm256_loadu_si256( (const __m256i *)( qj + 4 ) ), dj = _mm256_loadu_si256( (const __m256i *)( qj + 8 ) );
			__m256i xd = _mm256_xor_si256( di, dj );
			__m256i gor = _mm256_blend_epi32( _mm256_or_si256( ai, aj ), _mm256_or_si256( bi, bj ), 0xF0 ); // [G|G, G'|G']
			__m256i x0 = _mm256_permute2x128_si256( _mm256_xor_si256( ai, aj ), xd, 0x21 );          // [b0^b0, b0'^b0']
			__m256i x1 = _mm256_blend_epi32( _mm256_xor_si256( bi, bj ), xd, 0xF0 );                 // [b1^b1, b1'^b1']
			__m256i t0 = _mm256_or_si256( gor, x0 ), t1 = _mm256_or_si256( gor, x1 );
			bA = _mm256_add_epi8( bA, MSADIST_POPB( gor ) );
			bB = _mm256_add_epi8( bB, MSADIST_POPB( _mm256_or_si256( t0, x1 ) ) );
			b1 = _mm256_add_epi8( b1, MSADIST_POPB( _mm256_andnot_si256( t1, x0 ) ) );
			b2 = _mm256_add_epi8( b2, MSADIST_POPB( _mm256_andnot_si256( t0, x1 ) ) );
			MSADIST_ADC2( qi[0], qi[1], qj[0], qj[1], ci, cj, e0, e1 );
			MSADIST_ADC2( qi[6], qi[7], qj[6], qj[7], ci, cj, e2, e3 );
			nE += __builtin_popcountll( e0 ) + __builtin_popcountll( e1 ) + __builtin_popcountll( e2 ) + __builtin_popcountll( e3 );
		}
		sA = _mm256_add_epi64( sA, _mm256_sad_epu8( bA, z ) );
		sB = _mm256_add_epi64( sB, _mm256_sad_epu8( bB, z ) );
		s1 = _mm256_add_epi64( s1, _mm256_sad_epu8( b1, z ) );
		s2 = _mm256_add_epi64( s2, _mm256_sad_epu8( b2, z ) );
	}
	MSADIST_HSUM( sA, nA ); MSADIST_HSUM( sB, nB ); MSADIST_HSUM( s1, n1 ); MSADIST_HSUM( s2, n2 );
	if( nblk & 1 ) MSADIST_POPBLOCK( pi + 6*(nblk-1), pj + 6*(nblk-1) );
	return( msadist_xorsum( f, nA, nB, n1, n2, nE ) );
}

// AVX2 with the gap runs nE already known (msadist_gaprow_avx2): the residue
// popcounts only, without n1 when skip1 (its coefficient Td[1]-Td[3] is 0)
static inline __attribute__((always_inline,target("avx2,popcnt"))) long long msadist_xorresbody( const msadistfast_t *f, const uint64_t *pi, const uint64_t *pj, long long nE, const int skip1 )
{
	const __m256i lut4 = _mm256_setr_epi8( 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 );
	const __m256i m4 = _mm256_set1_epi8( 0x0f ), z = _mm256_setzero_si256();
	__m256i sA = z, sB = z, s1 = z, s2 = z;
	long long nA = 0, nB = 0, n1 = 0, n2 = 0;
	int nblk = f->W / 2, npr = nblk / 2, h;
	const uint64_t *qe = pi + 12 * npr, *qe1;
	while( pi < qe )
	{
		__m256i bA = z, bB = z, b1 = z, b2 = z; // byte counts, at most 8 per step
		qe1 = qe - pi > 12*31 ? pi + 12*31 : qe;
		for( ; pi<qe1; pi+=12, pj+=12 )
		{
			// two blocks = 12 words, as 3 vectors: [G b0] [b1 G'] [b0' b1'] (128 bits each)
			__m256i ai = _mm256_loadu_si256( (const __m256i *)pi ), bi = _mm256_loadu_si256( (const __m256i *)( pi + 4 ) ), di = _mm256_loadu_si256( (const __m256i *)( pi + 8 ) );
			__m256i aj = _mm256_loadu_si256( (const __m256i *)pj ), bj = _mm256_loadu_si256( (const __m256i *)( pj + 4 ) ), dj = _mm256_loadu_si256( (const __m256i *)( pj + 8 ) );
			__m256i xd = _mm256_xor_si256( di, dj );
			__m256i gor = _mm256_blend_epi32( _mm256_or_si256( ai, aj ), _mm256_or_si256( bi, bj ), 0xF0 ); // [G|G, G'|G']
			__m256i x0 = _mm256_permute2x128_si256( _mm256_xor_si256( ai, aj ), xd, 0x21 );          // [b0^b0, b0'^b0']
			__m256i x1 = _mm256_blend_epi32( _mm256_xor_si256( bi, bj ), xd, 0xF0 );                 // [b1^b1, b1'^b1']
			__m256i t0 = _mm256_or_si256( gor, x0 ), t1 = _mm256_or_si256( gor, x1 );
			bA = _mm256_add_epi8( bA, MSADIST_POPB( gor ) );
			bB = _mm256_add_epi8( bB, MSADIST_POPB( _mm256_or_si256( t0, x1 ) ) );
			if( !skip1 ) b1 = _mm256_add_epi8( b1, MSADIST_POPB( _mm256_andnot_si256( t1, x0 ) ) );
			b2 = _mm256_add_epi8( b2, MSADIST_POPB( _mm256_andnot_si256( t0, x1 ) ) );
		}
		sA = _mm256_add_epi64( sA, _mm256_sad_epu8( bA, z ) );
		sB = _mm256_add_epi64( sB, _mm256_sad_epu8( bB, z ) );
		if( !skip1 ) s1 = _mm256_add_epi64( s1, _mm256_sad_epu8( b1, z ) );
		s2 = _mm256_add_epi64( s2, _mm256_sad_epu8( b2, z ) );
	}
	MSADIST_HSUM( sA, nA ); MSADIST_HSUM( sB, nB ); MSADIST_HSUM( s1, n1 ); MSADIST_HSUM( s2, n2 );
	if( nblk & 1 ) for( h=0; h<2; h++ ) // last block
	{
		uint64_t gor = pi[h] | pj[h], x0 = pi[2+h] ^ pj[2+h], x1 = pi[4+h] ^ pj[4+h];
		uint64_t t0 = gor | x0, t1 = gor | x1;
		nA += __builtin_popcountll( gor );
		nB += __builtin_popcountll( t0 | x1 );
		n1 += __builtin_popcountll( x0 & ~t1 );
		n2 += __builtin_popcountll( x1 & ~t0 );
	}
	if( skip1 ) n1 = 0;
	return( msadist_xorsum( f, nA, nB, n1, n2, nE ) );
}

__attribute__((target("avx2,popcnt")))
static long long msadist_xorres_avx2_3( const msadistfast_t *f, const uint64_t *pi, const uint64_t *pj, long long nE )
{
	return( msadist_xorresbody( f, pi, pj, nE, 1 ) );
}

__attribute__((target("avx2,popcnt")))
static long long msadist_xorres_avx2_4( const msadistfast_t *f, const uint64_t *pi, const uint64_t *pj, long long nE )
{
	return( msadist_xorresbody( f, pi, pj, nE, 0 ) );
}

// AVX2: gap runs of row i with every later row, 4 rows per vector (64-bit lanes,
// rows 4g..4g+3 from the transposed masks g4).  Per lane the same multiword sum
// as MSADIST_RUNWORD: s = a + (a & ~b) + carry, carry-out = top bit of
// ( (a & ~b) | ( a & ~s ) ) (the majority of the top bits, as a & ~b is inside a).
// ne[j] = gap runs of pair (i,j) for i < j < nseq.
__attribute__((target("avx2")))
static void msadist_gaprow_avx2( const msadistfast_t *f, int i, int *ne )
{
	const __m256i lut4 = _mm256_setr_epi8( 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 );
	const __m256i m4 = _mm256_set1_epi8( 0x0f ), z = _mm256_setzero_si256();
	const int W = f->W, n = f->nseq;
	const uint64_t *gi = f->g4 + (size_t)( i >> 2 ) * W * 4 + ( i & 3 ); // word w of row i: gi[4*w]
	int g, w, w0, w1, l;
	long long v[4];
	for( g=( i + 1 ) >> 2; g*4<n; g++ )
	{
		const uint64_t *gj = f->g4 + (size_t)g * W * 4;
		__m256i ci = z, cj = z, sum = z;
		if( f->gact ) // work/v9p: only the words where row i or a row of the group has a residue
		{
			const uint64_t *ai = f->ract + (size_t)i * f->nmw, *ag = f->gact + (size_t)g * f->nmw;
			__m256i b8 = z; // byte counts, at most 8 per word
			int m, nb = 0;
			for( m=0; m<f->nmw; m++ )
			{
				uint64_t bits = ai[m] | ag[m];
				while( bits )
				{
					w = m * 64 + __builtin_ctzll( bits );
					bits &= bits - 1;
					{
						const __m256i a = _mm256_set1_epi64x( (long long)gi[4*w] );
						const __m256i b = _mm256_load_si256( (const __m256i *)( gj + 4*w ) );
						const __m256i xa = _mm256_andnot_si256( b, a ), xb = _mm256_andnot_si256( a, b );
						const __m256i sa = _mm256_add_epi64( _mm256_add_epi64( a, xa ), ci );
						const __m256i sb = _mm256_add_epi64( _mm256_add_epi64( b, xb ), cj );
						ci = _mm256_srli_epi64( _mm256_or_si256( xa, _mm256_andnot_si256( sa, a ) ), 63 );
						cj = _mm256_srli_epi64( _mm256_or_si256( xb, _mm256_andnot_si256( sb, b ) ), 63 );
						b8 = _mm256_add_epi8( b8, MSADIST_POPB( _mm256_or_si256( _mm256_andnot_si256( a, sa ), _mm256_andnot_si256( b, sb ) ) ) );
					}
					if( ++nb == 31 )
					{
						sum = _mm256_add_epi64( sum, _mm256_sad_epu8( b8, z ) );
						b8 = z;
						nb = 0;
					}
				}
			}
			sum = _mm256_add_epi64( sum, _mm256_sad_epu8( b8, z ) );
			_mm256_storeu_si256( (__m256i *)v, sum );
			for( l=0; l<4; l++ ) if( 4*g + l > i && 4*g + l < n ) ne[4*g+l] = (int)v[l];
			continue;
		}
		for( w0=0; w0<W; w0=w1 )
		{
			__m256i b8 = z; // byte counts, at most 8 per word
			w1 = w0 + 31 < W ? w0 + 31 : W;
			for( w=w0; w<w1; w++ )
			{
				const __m256i a = _mm256_set1_epi64x( (long long)gi[4*w] );
				const __m256i b = _mm256_load_si256( (const __m256i *)( gj + 4*w ) );
				const __m256i xa = _mm256_andnot_si256( b, a ), xb = _mm256_andnot_si256( a, b );
				const __m256i sa = _mm256_add_epi64( _mm256_add_epi64( a, xa ), ci );
				const __m256i sb = _mm256_add_epi64( _mm256_add_epi64( b, xb ), cj );
				ci = _mm256_srli_epi64( _mm256_or_si256( xa, _mm256_andnot_si256( sa, a ) ), 63 );
				cj = _mm256_srli_epi64( _mm256_or_si256( xb, _mm256_andnot_si256( sb, b ) ), 63 );
				b8 = _mm256_add_epi8( b8, MSADIST_POPB( _mm256_or_si256( _mm256_andnot_si256( a, sa ), _mm256_andnot_si256( b, sb ) ) ) );
			}
			sum = _mm256_add_epi64( sum, _mm256_sad_epu8( b8, z ) );
		}
		_mm256_storeu_si256( (__m256i *)v, sum );
		for( l=0; l<4; l++ ) if( 4*g + l > i && 4*g + l < n ) ne[4*g+l] = (int)v[l];
	}
}
#endif // FMSIMD_X86

// class of column c of a row from its planes: 4 for a gap, else the 2-bit code
static inline int msadist_class( const uint64_t *pl, uint32_t c )
{
	const uint64_t *q = pl + ( c >> 7 ) * 6 + ( ( c >> 6 ) & 1 );
	int sh = c & 63;
	if( ( q[0] >> sh ) & 1 ) return( 4 );
	return( (int)( ( ( q[2] >> sh ) & 1 ) | ( ( ( q[4] >> sh ) & 1 ) << 1 ) ) );
}

// Exceptions: an exception is code 0 in the planes, so the plane score counted
// Td[code(a)^code(b)] at its column; add amino_dis[a][b] - Td[code(a)^code(b)] for
// every column where one row holds an exception and the other a residue.  The
// two sorted lists carry the exception bytes; the other row's gap/code comes
// from its planes (just read by the kernel), so no sequence byte is loaded.
static long long msadist_exccorr( const msadistfast_t *f, int i, int j )
{
	const uint32_t *pi = f->excol[i], *pj = f->excol[j], *ei = pi + f->exn[i], *ej = pj + f->exn[j];
	const uint64_t *pli = f->planes + (size_t)i * 3 * f->W, *plj = f->planes + (size_t)j * 3 * f->W;
	long long s = 0;
	while( pi < ei && pj < ej )
	{
		uint32_t ci = *pi >> 8, cj = *pj >> 8;
		if( ci == cj ) { s += f->CB[( ( *pi & 0xff ) << 7 ) | ( *pj & 0xff )]; pi++; pj++; }
		else if( ci < cj ) { s += f->CXa[( ( *pi & 0xff ) << 3 ) | msadist_class( plj, ci )]; pi++; }
		else { s += f->CXb[( ( *pj & 0xff ) << 3 ) | msadist_class( pli, cj )]; pj++; }
	}
	for( ; pi<ei; pi++ ) s += f->CXa[( ( *pi & 0xff ) << 3 ) | msadist_class( plj, *pi >> 8 )];
	for( ; pj<ej; pj++ ) s += f->CXb[( ( *pj & 0xff ) << 3 ) | msadist_class( pli, *pj >> 8 )];
	return( s );
}

static inline __attribute__((always_inline)) long long msadist_genericbody( const msadistfast_t *f, int i, int j )
{
	const uint64_t *mi = f->masks + (size_t)i * 2 * f->W, *mj = f->masks + (size_t)j * 2 * f->W;
	const unsigned char *si = (const unsigned char *)f->seq[i], *sj = (const unsigned char *)f->seq[j];
	const int *T = f->T;
	uint64_t ci = 0, cj = 0, e, rr;
	long long s = 0, nE = 0;
	int w, W = f->W;
	for( w=0; w<W; w++ )
	{
		const unsigned char *a = si + w*64, *b = sj + w*64;
		rr = mi[W+w] & mj[W+w];
		while( rr )
		{
			int k = __builtin_ctzll( rr );
			s += T[(a[k]<<7)|b[k]];
			rr &= rr - 1;
		}
		MSADIST_RUNWORD( mi[w], mj[w], ci, cj, e );
		nE += __builtin_popcountll( e );
	}
	return( s + (long long)f->penal * nE );
}

static long long msadist_genericscore_c( const msadistfast_t *f, int i, int j )
{
	return( msadist_genericbody( f, i, j ) );
}

#if FMSIMD_X86
__attribute__((target("popcnt")))
static long long msadist_genericscore_popcnt( const msadistfast_t *f, int i, int j )
{
	return( msadist_genericbody( f, i, j ) );
}
#endif

// Column path, row i: acc[j] += amino_dis[ seq[i][c] ][ seq[j][c] ] for every
// later row j and every column c where both have a residue (acc[j] starts at 0;
// the caller clears acc[i+1..] after use).
#if FMSIMD_X86
// Dense column (AVX2): acc[j] += amino_dis[a][ residue of row j ] for the rows j of
// the 32-row blocks from j0 & ~31 to np32 (codes of gaps and padding rows give 0).
// Rows j < j0 of the first block and the padding rows also receive values; the
// caller clears them (msadist_colrow).
__attribute__((target("avx2")))
static void msadist_densecol_avx2( const uint8_t *code, const uint8_t *tab, int j0, int np32, int *acc )
{
	const __m256i t0l = _mm256_broadcastsi128_si256( _mm_loadu_si128( (const __m128i *)tab ) );
	const __m256i t1l = _mm256_broadcastsi128_si256( _mm_loadu_si128( (const __m128i *)( tab + 16 ) ) );
	const __m256i t0h = _mm256_broadcastsi128_si256( _mm_loadu_si128( (const __m128i *)( tab + 32 ) ) );
	const __m256i t1h = _mm256_broadcastsi128_si256( _mm_loadu_si128( (const __m128i *)( tab + 48 ) ) );
	const __m256i m15 = _mm256_set1_epi8( 15 ), m16 = _mm256_set1_epi8( 16 );
	int j;
	for( j=j0&~31; j<np32; j+=32 )
	{
		const __m256i c = _mm256_loadu_si256( (const __m256i *)( code + j ) );
		const __m256i lo4 = _mm256_and_si256( c, m15 );
		const __m256i hi = _mm256_cmpeq_epi8( _mm256_and_si256( c, m16 ), m16 ); // codes 16..31
		const __m256i vl = _mm256_blendv_epi8( _mm256_shuffle_epi8( t0l, lo4 ), _mm256_shuffle_epi8( t1l, lo4 ), hi );
		const __m256i vh = _mm256_blendv_epi8( _mm256_shuffle_epi8( t0h, lo4 ), _mm256_shuffle_epi8( t1h, lo4 ), hi );
		const __m256i w0 = _mm256_unpacklo_epi8( vl, vh ), w1 = _mm256_unpackhi_epi8( vl, vh ); // int16: rows j+0..7,16..23 / j+8..15,24..31
		__m256i *a = (__m256i *)( acc + j );
		_mm256_storeu_si256( a,     _mm256_add_epi32( _mm256_loadu_si256( a ),     _mm256_cvtepi16_epi32( _mm256_castsi256_si128( w0 ) ) ) );
		_mm256_storeu_si256( a + 1, _mm256_add_epi32( _mm256_loadu_si256( a + 1 ), _mm256_cvtepi16_epi32( _mm256_castsi256_si128( w1 ) ) ) );
		_mm256_storeu_si256( a + 2, _mm256_add_epi32( _mm256_loadu_si256( a + 2 ), _mm256_cvtepi16_epi32( _mm256_extracti128_si256( w0, 1 ) ) ) );
		_mm256_storeu_si256( a + 3, _mm256_add_epi32( _mm256_loadu_si256( a + 3 ), _mm256_cvtepi16_epi32( _mm256_extracti128_si256( w1, 1 ) ) ) );
	}
}
#endif

// acc must hold np32 + 32 ints when f->ndense
static void msadist_colrow( const msadistfast_t *f, int i, int *acc )
{
	const unsigned char *s = (const unsigned char *)f->seq[i];
	const uint64_t *mi = f->masks + (size_t)i * 2 * f->W;
	const uint32_t *rp = f->rpos + f->roff[i], *ce = f->centry;
	int w, c, k = 0, j;
	for( w=0; w<f->W; w++ )
	{
		uint64_t r = mi[f->W+w];
		while( r )
		{
			const int *Tr;
			const uint32_t *e, *end;
			c = w*64 + __builtin_ctzll( r );
#if FMSIMD_X86
			if( f->ndense && f->dcol[c] >= 0 )
			{
				msadist_densecol_avx2( f->dcode + (size_t)f->dcol[c] * f->np32, f->dtab + s[c] * 64, i + 1, f->np32, acc );
				k++;
				r &= r - 1;
				continue;
			}
#endif
			Tr = f->T + ( s[c] << 7 );
			e = ce + rp[k++] + 1;
			end = ce + f->coff[c+1];
			for( ; e<end; e++ ) acc[*e >> 8] += Tr[*e & 0xff];
			r &= r - 1;
		}
	}
	if( f->ndense ) // the dense columns also wrote rows ( i+1 ) & ~31 .. i and the padding rows
	{
		for( j=( i + 1 ) & ~31; j<=i; j++ ) acc[j] = 0;
		for( j=f->nseq; j<f->np32; j++ ) acc[j] = 0;
	}
}

// gap runs of pair (i,j) from the gap masks (W contiguous words per row)
static inline __attribute__((always_inline)) long long msadist_gaprunsbody( const msadistfast_t *f, int i, int j )
{
	const uint64_t *mi = f->masks + (size_t)i * 2 * f->W, *mj = f->masks + (size_t)j * 2 * f->W, *me = mi + f->W;
	long long nE = 0;
	if( f->ract ) // work/v9p: only the words where either row has a residue
	{
		const uint64_t *ai = f->ract + (size_t)i * f->nmw, *aj = f->ract + (size_t)j * f->nmw;
		int m, w;
#if FMSIMD_X86
		unsigned char ci = 0, cj = 0;
		unsigned long long sa, sb;
#else
		uint64_t ci = 0, cj = 0;
#endif
		uint64_t e;
		for( m=0; m<f->nmw; m++ )
		{
			uint64_t bits = ai[m] | aj[m];
			while( bits )
			{
				w = m * 64 + __builtin_ctzll( bits );
				bits &= bits - 1;
#if FMSIMD_X86
				{
					unsigned long long ga = mi[w], gb = mj[w];
					ci = _addcarry_u64( ci, ga, ga & ~gb, &sa );
					cj = _addcarry_u64( cj, gb, gb & ~ga, &sb );
					e = ( sa & ~ga ) | ( sb & ~gb );
				}
#else
				MSADIST_RUNWORD( mi[w], mj[w], ci, cj, e );
#endif
				nE += __builtin_popcountll( e );
			}
		}
		return( nE );
	}
#if FMSIMD_X86
	unsigned char ci = 0, cj = 0;
	uint64_t e0, e1;
	for( ; mi<me; mi+=2, mj+=2 ) // W is even
	{
		MSADIST_ADC2( mi[0], mi[1], mj[0], mj[1], ci, cj, e0, e1 );
		nE += __builtin_popcountll( e0 ) + __builtin_popcountll( e1 );
	}
#else
	uint64_t ci = 0, cj = 0, e;
	for( ; mi<me; mi++, mj++ )
	{
		MSADIST_RUNWORD( mi[0], mj[0], ci, cj, e );
		nE += __builtin_popcountll( e );
	}
#endif
	return( nE );
}

static long long msadist_gapruns_c( const msadistfast_t *f, int i, int j )
{
	return( msadist_gaprunsbody( f, i, j ) );
}

#if FMSIMD_X86
__attribute__((target("popcnt")))
static long long msadist_gapruns_popcnt( const msadistfast_t *f, int i, int j )
{
	return( msadist_gaprunsbody( f, i, j ) );
}
#endif

// Per-row work before the pairs (i,j), j > i: the column path's residue sums
// (acc, must start zeroed; the caller clears acc[i+1..] afterwards) and, with
// ne4, the gap runs of all those pairs (ne).  Either pointer may be NULL when
// unused (f->colpath == 0, f->ne4 == 0).
static void msadist_rowprep( const msadistfast_t *f, int i, int *acc, int *ne )
{
	if( f->colpath ) msadist_colrow( f, i, acc );
#if FMSIMD_X86
	if( f->ne4 ) msadist_gaprow_avx2( f, i, ne );
#endif
	(void)ne;
}

// exact value of naivepairscorefast( seq[i], seq[j], skip[i], skip[j], penal ) for this
// alignment; acc and ne as left by msadist_rowprep( f, i, acc, ne )
static inline double msadistfast_score( const msadistfast_t *f, int i, int j, const int *acc, const int *ne )
{
	long long s;
	if( f->colpath )
	{
		if( f->ne4 ) s = ne[j];
		else
#if FMSIMD_X86
		if( f->kern >= 1 ) s = msadist_gapruns_popcnt( f, i, j );
		else
#endif
			s = msadist_gapruns_c( f, i, j );
		s = acc[j] + (long long)f->penal * s;
	}
	else if( f->xorpath )
	{
		const uint64_t *pi = f->planes + (size_t)i * 3 * f->W, *pj = f->planes + (size_t)j * 3 * f->W;
#if MSADIST_NEON
		if( f->kern >= 1 ) s = msadist_xorscore_neon( f, pi, pj );
		else
#elif FMSIMD_X86
		if( f->ne4 ) s = f->skip1 ? msadist_xorres_avx2_3( f, pi, pj, ne[j] ) : msadist_xorres_avx2_4( f, pi, pj, ne[j] );
		else if( f->kern == 2 ) s = msadist_xorscore_avx2( f, pi, pj );
		else if( f->kern == 1 ) s = msadist_xorscore_popcnt( f, pi, pj );
		else
#endif
			s = msadist_xorscore_c( f, pi, pj );
		if( f->nexckinds && ( f->exn[i] | f->exn[j] ) ) s += msadist_exccorr( f, i, j );
	}
	else
	{
#if FMSIMD_X86
		if( f->kern >= 1 ) s = msadist_genericscore_popcnt( f, i, j );
		else
#endif
			s = msadist_genericscore_c( f, i, j );
	}
	return( (double)s );
}
/* END msadistfast kernel */

static double msadist_pairscore( msadistmtxthread_arg_t *targ, int i, int j, const int *acc, const int *ne )
{
	msadistfast_t *f = (msadistfast_t *)targ->fast;
	double v;
	if( f == NULL ) return( naivepairscorefast( targ->seq[i], targ->seq[j], targ->skiptable[i], targ->skiptable[j], penalty_dist ) );
	v = msadistfast_score( f, i, j, acc, ne );
	if( f->check )
	{
		double r = naivepairscorefast( targ->seq[i], targ->seq[j], targ->skiptable[i], targ->skiptable[j], penalty_dist );
		if( memcmp( &r, &v, sizeof( double ) ) )
		{
			reporterr( "msadist check FAILED at i=%d j=%d: fast=%.17g ref=%.17g\n", i, j, v, r );
			exit( 1 );
		}
	}
	return( v );
}

static void msadist_savemsa( const char *fn, int njob, char **seq, int *selfscore ) // debug: input for a standalone kernel benchmark
{
	FILE *fp = fopen( fn, "wb" );
	int i, L = strlen( seq[0] );
	if( fp == NULL ) { reporterr( "Cannot open %s\n", fn ); return; }
	fwrite( &penalty_dist, sizeof( int ), 1, fp );
	for( i=0; i<0x80; i++ ) fwrite( amino_dis[i], sizeof( int ), 0x80, fp );
	fwrite( &njob, sizeof( int ), 1, fp );
	fwrite( &L, sizeof( int ), 1, fp );
	fwrite( selfscore, sizeof( int ), njob, fp );
	for( i=0; i<njob; i++ ) fwrite( seq[i], 1, L, fp );
	fclose( fp );
}

static void msadist_dump( const char *fn, int njob, double **iscore )
{
	FILE *fp = fopen( fn, "wb" );
	int i;
	if( fp == NULL ) { reporterr( "Cannot open %s\n", fn ); return; }
	for( i=0; i<njob-1; i++ ) fwrite( iscore[i]+1, sizeof( double ), njob-1-i, fp );
	fclose( fp );
}

static double sextet_now( void ); // work/v9p (sextetfast block)
#if 1
static void *msadistmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	msadistmtxthread_arg_t *targ = (msadistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **iscore = targ->iscore;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	Jobtable *jobpospt = targ->jobpospt;


	double ssi, ssj, bunbo, iscoretmp;
	int i, j;
	int nlim = njob-1;
	msadistfast_t *mf = (msadistfast_t *)targ->fast; // fast-mafft
	int *acc = ( mf && mf->colpath ) ? calloc( njob + 64, sizeof( int ) ) : NULL; // column path: residue sums of row i with every later row (+64: dense-column blocks)
	int *ne = ( mf && mf->ne4 ) ? calloc( njob, sizeof( int ) ) : NULL; // AVX2: gap runs of row i with every later row
	double v9tcol = 0.0, v9tgap = 0.0, v9tpair = 0.0, v9tp0 = 0.0; // work/v9p: MAFFT_MSADIST=time
#if FMSIMD_X86
	int v9nt = mf && njob >= 4096 && !( getenv( "MAFFT_MSADIST" ) && strstr( getenv( "MAFFT_MSADIST" ), "nont" ) ); // work/v9p: non-temporal matrix stores (large N: the rows are not in any cache)
#endif

	while( 1 )
	{
#ifdef enablemultithread
		if( nthreadpair ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = jobpospt->i; // (jobpospt-i)++ dato, shuuryou hantei no mae ni ++ surunode, tomaranakunaru.

			if( i == nlim )
			{
				if( mf && mf->timing ) { mf->tcol += v9tcol; mf->tgap += v9tgap; mf->tpair += v9tpair; } // under the mutex
				pthread_mutex_unlock( targ->mutex );
				free( acc ); free( ne );
				return( NULL );
			}
			jobpospt->i += 1;
			pthread_mutex_unlock( targ->mutex );
			if( i % 100 == 0 ) fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
		}
		else
#endif
		{
			i = (jobpospt->i)++;
			if( i == nlim ) { if( mf && mf->timing ) { mf->tcol += v9tcol; mf->tgap += v9tgap; mf->tpair += v9tpair; } free( acc ); free( ne ); return( NULL ); }
			if( i % 100 == 0 ) fprintf( stderr, "\r% 5d / %d", i, njob );
		}

		ssi = selfscore[i];
		if( mf && mf->timing ) // work/v9p
		{
			double t0_ = sextet_now(), t1_;
			if( mf->colpath ) msadist_colrow( mf, i, acc );
			t1_ = sextet_now(); v9tcol += t1_ - t0_;
#if FMSIMD_X86
			if( mf->ne4 ) msadist_gaprow_avx2( mf, i, ne );
#endif
			v9tgap += sextet_now() - t1_;
		}
		else
		if( acc || ne ) msadist_rowprep( mf, i, acc, ne );
		if( mf && mf->timing ) v9tp0 = sextet_now();
		for( j=i+1; j<njob; j++ )
		{
			ssj = selfscore[j];
			bunbo = MIN( ssi, ssj );
//fprintf( stderr, "bunbo = %f\n", bunbo );
//fprintf( stderr, "naivepairscorefast() = %f\n", naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) );
			if( bunbo == 0.0 )
				iscoretmp = 2.0; // 2013/Oct/17
			else
			{
				iscoretmp = ( 1.0 - msadist_pairscore( targ, i, j, acc, ne ) / bunbo ) * 2.0; // 2014/Aug/15 fast; fast-mafft: exact replacement for naivepairscorefast()
				if( iscoretmp > 10 ) iscoretmp = 10.0;  // 2015/Mar/17
	
			}
			if( iscoretmp < 0.0 ) 
			{
				reporterr( "WARNING: negative distance, iscoretmp = %f\n", iscoretmp );
				iscoretmp = 0.0;
			}
#if FMSIMD_X86
			if( v9nt ) // work/v9p: the same 8 bytes, written around the caches (no read of the old line)
			{
				long long b_;
				memcpy( &b_, &iscoretmp, sizeof( double ) );
				_mm_stream_si64( (long long *)( iscore[i] + ( j - i ) ), b_ );
			}
			else
#endif
			iscore[i][j-i] = iscoretmp;
//			printf( "i,j=%d,%d, iscoretmp=%f\n", i, j, iscoretmp );

		}
#if FMSIMD_X86
		if( v9nt ) _mm_sfence();
#endif
		if( mf && mf->timing ) v9tpair += sextet_now() - v9tp0;
		if( acc ) memset( acc + i + 1, 0, ( njob - i - 1 ) * sizeof( int ) );
	}
}
#else
static void *msadistmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	msadistmtxthread_arg_t *targ = (msadistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **iscore = targ->iscore;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	Jobtable *jobpospt = targ->jobpospt;


	double ssi, ssj, bunbo, iscoretmp;
	int i, j;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthreadpair ) pthread_mutex_lock( targ->mutex );
#endif
		j = jobpospt->j;
		i = jobpospt->i;
		j++;
		if( j == njob )
		{
			i++;
			j = i + 1;
			if( i == njob-1 )
			{
#ifdef enablemultithread
				if( nthreadpair ) pthread_mutex_unlock( targ->mutex );
#endif
				return( NULL );
			}
		}
		jobpospt->j = j;
		jobpospt->i = i;
#ifdef enablemultithread
		if( nthreadpair ) pthread_mutex_unlock( targ->mutex );
#endif


		if( nthreadpair )
		{
			if( j==i+1 && i % 10 == 0 ) fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
		}
		else
		{
			if( j==i+1 && i % 10 == 0 ) fprintf( stderr, "\r% 5d / %d", i, njob );
		}
		ssi = selfscore[i];
		ssj = selfscore[j];
		bunbo = MIN( ssi, ssj );
//fprintf( stderr, "bunbo = %f\n", bunbo );
//fprintf( stderr, "naivepairscorefast() = %f\n", naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) );
		if( bunbo == 0.0 )
			iscoretmp = 2.0; // 2013/Oct/17
		else
		{
			iscoretmp = ( 1.0 - naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) / bunbo ) * 2.0; // 2014/Aug/15 fast 
			if( iscoretmp > 10 ) iscoretmp = 10.0;  // 2015/Mar/17

		}
		iscore[i][j-i] = iscoretmp;


	}
}
#endif

#ifdef enablemultithread
static void *distancematrixthread( void *arg )
{
	distancematrixthread_arg_t *targ = (distancematrixthread_arg_t *)arg;
	int thread_no = targ->thread_no;
	int njob = targ->njob;
	int *jobpospt = targ->jobpospt;
	int **pointt = targ->pointt;
	double **mtx = targ->mtx;

	int *table1;
	int i, j;

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		i = *jobpospt;
		if( i == njob )
		{
			pthread_mutex_unlock( targ->mutex );
			commonsextet_p( NULL, NULL );
			return( NULL );
		}
		*jobpospt = i+1;
		pthread_mutex_unlock( targ->mutex );

		table1 = (int *)calloc( tsize, sizeof( int ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		if( i % 100 == 0 )
		{
			reporterr(       "\r% 5d / %d (thread %4d)", i+1, njob, thread_no );
		}
		makecompositiontable_p( table1, pointt[i] );

		for( j=i; j<njob; j++ ) 
		{
			mtx[i][j-i] = (double)commonsextet_p( table1, pointt[j] );
		} 
		free( table1 );
	}
}
#endif // enablemultithread

/*
 * fast-mafft (work/sextet): exact fast replacement for the pair loop of the
 * "Making a distance matrix .." (6-mer) stage.
 *
 * The original loop stores, for every pair j >= i,
 *     mtx[i][j-i] = (double)commonsextet_p( table_i, pointt[j] )
 * where table_i[p] = c_i(p) is the number of times tuple p occurs in pointt[i]
 * and commonsextet_p() returns the int  sum over p of min( c_i(p), c_j(p) )
 * (0 when pointt[j] is the short-sequence marker -1; table_i is all zero when
 * pointt[i] is).  Only this integer is computed differently here.  The serial
 * conversion loop that follows (count -> distance with lenfac and
 * bunbo = MIN( mtx[i][0], mtx[j][0] )) is done per pair by the worker threads
 * right after the count, with the same statements on the same double values
 * (so the same fmadd contraction and rounding), and the original loop is
 * skipped.  The diagonal mtx[i][0] is the number of tuples of i.
 *
 * DNA (tsize 4096): each row is a dense table t_i[p] = min( c_i(p), 255 ) of
 * 4096 bytes.  With T_i = sum_p t_i(p),
 *     sum_p min( t_i, t_j ) = ( T_i + T_j - sum_p |t_i - t_j| ) / 2
 * and the absolute-difference sum is one NEON uaba per 16 tuples per pair,
 * computed for 4x4 blocks of rows so each loaded vector is used 4 times.  The
 * uaba bytes are widened into u16 lanes every K = 255/mc 16-tuple chunks,
 * mc = largest count in the 8 rows, so no byte overflows; the u16 lanes hold
 * at most 2*256*mc, so this path needs mc <= 127.  Blocks with a larger count
 * use a plain loop over min( t_i, t_j ).  Saturation: min( t_i, t_j ) equals
 * min( c_i, c_j ) unless both counts are >= 255; for those tuples the exact
 * counts are kept in a per-row list and  min( c_i, c_j ) - 255  is added.
 * Rows are handed out 16 at a time (a "super-row"); for each group of 4
 * columns the 4 column tables are loaded once for all 4 row blocks.
 *
 * x86-64 (work/x86simd): the same identity with psadbw (SSE2) / vpsadbw (AVX2),
 * which sums |a-b| of 8 bytes into a 64-bit lane, so there is no limit on the
 * counts and every block uses it (plus the saturation correction).
 *
 * Amino acids (tsize 46656): each row keeps its list of distinct tuples with
 * their counts (u16).  For row i a dense u16 table of c_i is filled once and
 * each j sums min( table[p], c_j(p) ) over j's list.  If any count exceeds
 * 65535 the original code is used.  Other tuple sizes use the original code.
 * Inverted index (kind 3, work/x86simd): for every tuple the list of rows that
 * hold it, ascending; row i adds min( c_i(p), c_j(p) ) to acc[j] for every
 * tuple p of i and every later row j in p's list.  That costs the sum over
 * tuples of (rows holding it)^2/2 instead of n^2/2 times tuples per row
 * (15-70x fewer steps on homfam sets); it is used when the estimate says so.
 *
 * Debug (environment):
 *   MAFFT_SEXTET=ref     original loop
 *   MAFFT_SEXTET=check   also run commonsextet_p() for every pair and exit(1)
 *                        on any bitwise difference of the stored double
 *   MAFFT_SEXTET=scalar  DNA: skip the SIMD block path (plain loop + overflow lists);
 *                        MAFFT_SIMD=scalar does the same
 *   MAFFT_SEXTET=list    DNA: use the list path
 *   MAFFT_SEXTET=inv     list path: always use the inverted index (if it fits)
 *   MAFFT_SEXTET=noinv   list path: never use the inverted index
 *   MAFFT_SEXTET=counts  store counts only; the original conversion loop runs
 *   MAFFT_SEXTET=nolftab compute lenfac for every pair (no length table)
 *   MAFFT_SEXTET=time    print the time of each phase
 *   MAFFT_SEXTET_DUMP=f  write the finished distance matrix (row i, j>=i, doubles)
 *   (these combine, e.g. "scalar,check")
 */
/* BEGIN sextetfast kernel */
#include <stdint.h>
#include <sys/time.h>
#if defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#define SEXTET_NEON 1
#else
#define SEXTET_NEON 0
#endif

#define SEXTET_DW 4096 // DNA 6-mer table width (4^6)
#define SEXTET_SR 16   // rows per job on the DNA path (a multiple of 4)
#define SEXTET_LR 4    // rows per job on the list path
#define SEXTET_LFMAX 512 // largest span of lengths for the lenfac table

typedef struct _sextetfast
{
	int n, npad;
	int kind;        // 1: dense DNA tables, 2: tuple lists, 3: inverted index of the tuple lists, 4: dense element masks + index of the rest (work/v9p)
	int check, scalar, timing;
	int simd;        // fmsimd level (scalar, MAFFT_SIMD and the CPU)
	long nsimd, nslow; // dense path, with timing: 4x4 blocks by the SIMD kernel / by the plain loop
	int **pointt;
	double **mtx;
	int *nogaplen;
	double *self;    // per row: number of tuples = the diagonal count mtx[i][0]
	int convert;     // 1: store distances (the conversion of the original code), 0: counts
	double *lftab;   // lftab[(li-lflo)*lfn + (lj-lflo)] = sextet_lenfac( li, lj ), or NULL
	int lflo, lfn;
	int tw;          // table width (tsize)
	/* dense path */
	uint8_t *tab;    // npad rows of SEXTET_DW bytes, min( count, 255 ); rows >= n stay zero
	int *tot;        // per row: sum of the saturated counts
	unsigned char *mc; // per row: largest saturated count
	int *ovn;        // per row: number of tuples with count >= 255 (almost always 0)
	int **ovp, **ovc;  // their tuples (ascending) and exact counts
	/* list path */
	int *ln;         // per row: number of distinct tuples
	size_t *off;     // per row: offset into lp/lc
	uint16_t *lp, *lc; // tuple, count
	/* inverted index (kind 3): for every tuple the rows that hold it, ascending */
	uint32_t *post;  // ( row << cb ) | count
	size_t *poff;    // tw+1 offsets into post
	uint32_t *ppos;  // per list entry (same layout as lp/lc): index of the row's own entry in post
	int cb;          // count bits
	double listwork, invwork; // estimated operations of the two pair loops
	/* kind 4 (work/v9p) */
	int dk, dw;      // dense elements, 64-bit words per row mask
	uint64_t *dmask; // npad rows of dw words
	size_t *rstart;  // per row: its entries of the index of the rest, rstart[i] .. rstart[i+1]
	uint32_t *rpp;   // per such entry: its position in post (whose count is c - kd)
	uint16_t *rtup;  // and its tuple (post[poff[p] .. poff[p+1]) are p's rows)
	size_t nres;
	int popcnt;      // x86: the POPCNT version of the pair loop
	double densework;
	/* jobs */
	int phase;       // 0: build tables/lists, 1: pairs
	int jobpos;
	int bad;         // list path: a count does not fit in u16
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} sextetfast_t;

static double sextet_now( void )
{
	struct timeval tv;
	gettimeofday( &tv, NULL );
	return( tv.tv_sec + 1e-6 * tv.tv_usec );
}

static int sextet_claim( sextetfast_t *s, int step )
{
	int i;
#ifdef enablemultithread
	if( s->mutex ) pthread_mutex_lock( s->mutex );
#endif
	i = s->jobpos;
	if( i < s->n ) s->jobpos = i + step;
#ifdef enablemultithread
	if( s->mutex ) pthread_mutex_unlock( s->mutex );
#endif
	return( i );
}

static void sextet_buildrow( sextetfast_t *s, int i, int *cnt )
{
	int *p = s->pointt[i];
	int k, m, point, tot, mc;

	m = 0;
	if( *p != -1 ) while( p[m] != END_OF_VEC ) m++;
	s->self[i] = (double)m; // commonsextet_p( table_i, pointt[i] ) counts every tuple of i

	if( s->kind == 1 )
	{
		uint8_t *t = s->tab + (size_t)i * SEXTET_DW;
		if( *p != -1 ) for( ; ( point = *p ) != END_OF_VEC; p++ ) if( t[point] != 255 ) t[point]++;
		tot = mc = 0;
		for( k=0; k<SEXTET_DW; k++ )
		{
			tot += t[k];
			if( t[k] > mc ) mc = t[k];
		}
		s->tot[i] = tot;
		s->mc[i] = mc;
		if( mc == 255 ) // exact counts of the saturated tuples
		{
			for( p=s->pointt[i]; ( point = *p ) != END_OF_VEC; p++ ) cnt[point]++;
			for( m=0, k=0; k<SEXTET_DW; k++ ) if( cnt[k] >= 255 ) m++;
			s->ovn[i] = m;
			s->ovp[i] = calloc( m, sizeof( int ) );
			s->ovc[i] = calloc( m, sizeof( int ) );
			for( m=0, k=0; k<SEXTET_DW; k++ )
			{
				if( cnt[k] >= 255 )
				{
					s->ovp[i][m] = k;
					s->ovc[i][m] = cnt[k];
					m++;
				}
				cnt[k] = 0;
			}
		}
	}
	else
	{
		uint16_t *lp = s->lp + s->off[i], *lc = s->lc + s->off[i];
		m = 0;
		if( *p != -1 ) for( ; ( point = *p ) != END_OF_VEC; p++ ) if( cnt[point]++ == 0 ) lp[m++] = point;
		for( k=0; k<m; k++ )
		{
			if( cnt[lp[k]] > 65535 ) s->bad = 1;
			lc[k] = cnt[lp[k]];
			cnt[lp[k]] = 0;
		}
		s->ln[i] = m;
	}
}

static int sextet_ovcorr( sextetfast_t *s, int i, int j )
{
	int a = 0, b = 0, v = 0;
	int na = s->ovn[i], nb = s->ovn[j];
	int *pa = s->ovp[i], *pb = s->ovp[j], *ca = s->ovc[i], *cb = s->ovc[j];
	while( a < na && b < nb )
	{
		if( pa[a] < pb[b] ) a++;
		else if( pa[a] > pb[b] ) b++;
		else
		{
			v += MIN( ca[a], cb[b] ) - 255;
			a++; b++;
		}
	}
	return( v );
}

static int sextet_pairslow( sextetfast_t *s, int i, int j )
{
	const uint8_t *a = s->tab + (size_t)i * SEXTET_DW;
	const uint8_t *b = s->tab + (size_t)j * SEXTET_DW;
	unsigned int v = 0;
	int k;
	for( k=0; k<SEXTET_DW; k++ ) v += ( a[k] < b[k] ) ? a[k] : b[k];
	if( s->ovn[i] && s->ovn[j] ) v += sextet_ovcorr( s, i, j );
	return( (int)v );
}

#if SEXTET_NEON
// d[4r+c] = sum_p | ti[r][p] - tj[c][p] |  for 4 rows ti (consecutive tables) and 4 columns tj;
// every count is <= mc <= 127
static void sextet_absdiff4x4( const uint8_t *ti, const uint8_t *tj, int mc, uint32_t *d )
{
	const uint8_t *i0 = ti, *i1 = ti + SEXTET_DW, *i2 = ti + 2 * SEXTET_DW, *i3 = ti + 3 * SEXTET_DW;
	const uint8_t *j0 = tj, *j1 = tj + SEXTET_DW, *j2 = tj + 2 * SEXTET_DW, *j3 = tj + 3 * SEXTET_DW;
	int K = ( mc > 0 ) ? 255 / mc : SEXTET_DW / 16; // 16-byte chunks per u8 accumulation
	int g, e, k;
	uint16x8_t w00 = vdupq_n_u16( 0 ), w01 = w00, w02 = w00, w03 = w00, w10 = w00, w11 = w00, w12 = w00, w13 = w00,
	           w20 = w00, w21 = w00, w22 = w00, w23 = w00, w30 = w00, w31 = w00, w32 = w00, w33 = w00;
	for( g=0; g<SEXTET_DW/16; g+=K )
	{
		uint8x16_t a00 = vdupq_n_u8( 0 ), a01 = a00, a02 = a00, a03 = a00, a10 = a00, a11 = a00, a12 = a00, a13 = a00,
		           a20 = a00, a21 = a00, a22 = a00, a23 = a00, a30 = a00, a31 = a00, a32 = a00, a33 = a00;
		e = g + K;
		if( e > SEXTET_DW/16 ) e = SEXTET_DW/16;
		for( k=g*16; k<e*16; k+=16 )
		{
			uint8x16_t x0 = vld1q_u8( i0 + k ), x1 = vld1q_u8( i1 + k ), x2 = vld1q_u8( i2 + k ), x3 = vld1q_u8( i3 + k );
			uint8x16_t y0 = vld1q_u8( j0 + k ), y1 = vld1q_u8( j1 + k ), y2 = vld1q_u8( j2 + k ), y3 = vld1q_u8( j3 + k );
			a00 = vabaq_u8( a00, x0, y0 ); a01 = vabaq_u8( a01, x0, y1 ); a02 = vabaq_u8( a02, x0, y2 ); a03 = vabaq_u8( a03, x0, y3 );
			a10 = vabaq_u8( a10, x1, y0 ); a11 = vabaq_u8( a11, x1, y1 ); a12 = vabaq_u8( a12, x1, y2 ); a13 = vabaq_u8( a13, x1, y3 );
			a20 = vabaq_u8( a20, x2, y0 ); a21 = vabaq_u8( a21, x2, y1 ); a22 = vabaq_u8( a22, x2, y2 ); a23 = vabaq_u8( a23, x2, y3 );
			a30 = vabaq_u8( a30, x3, y0 ); a31 = vabaq_u8( a31, x3, y1 ); a32 = vabaq_u8( a32, x3, y2 ); a33 = vabaq_u8( a33, x3, y3 );
		}
		w00 = vpadalq_u8( w00, a00 ); w01 = vpadalq_u8( w01, a01 ); w02 = vpadalq_u8( w02, a02 ); w03 = vpadalq_u8( w03, a03 );
		w10 = vpadalq_u8( w10, a10 ); w11 = vpadalq_u8( w11, a11 ); w12 = vpadalq_u8( w12, a12 ); w13 = vpadalq_u8( w13, a13 );
		w20 = vpadalq_u8( w20, a20 ); w21 = vpadalq_u8( w21, a21 ); w22 = vpadalq_u8( w22, a22 ); w23 = vpadalq_u8( w23, a23 );
		w30 = vpadalq_u8( w30, a30 ); w31 = vpadalq_u8( w31, a31 ); w32 = vpadalq_u8( w32, a32 ); w33 = vpadalq_u8( w33, a33 );
	}
	d[0]  = vaddlvq_u16( w00 ); d[1]  = vaddlvq_u16( w01 ); d[2]  = vaddlvq_u16( w02 ); d[3]  = vaddlvq_u16( w03 );
	d[4]  = vaddlvq_u16( w10 ); d[5]  = vaddlvq_u16( w11 ); d[6]  = vaddlvq_u16( w12 ); d[7]  = vaddlvq_u16( w13 );
	d[8]  = vaddlvq_u16( w20 ); d[9]  = vaddlvq_u16( w21 ); d[10] = vaddlvq_u16( w22 ); d[11] = vaddlvq_u16( w23 );
	d[12] = vaddlvq_u16( w30 ); d[13] = vaddlvq_u16( w31 ); d[14] = vaddlvq_u16( w32 ); d[15] = vaddlvq_u16( w33 );
}
#endif

#if FMSIMD_X86
// d[4r+c] = sum_p | ti[r][p] - tj[c][p] | for 4 rows ti and 4 columns tj, any
// byte values (psadbw sums 8 absolute byte differences into a 64-bit lane;
// a pair's total is at most 4096*255).  Two passes of 4 rows x 2 columns keep
// the 8 accumulators and 6 loaded vectors in the 16 vector registers.
static void sextet_absdiff4x4_sse2( const uint8_t *ti, const uint8_t *tj, uint32_t *d )
{
	const uint8_t *i0 = ti, *i1 = ti + SEXTET_DW, *i2 = ti + 2 * SEXTET_DW, *i3 = ti + 3 * SEXTET_DW;
	int h, k;
	for( h=0; h<4; h+=2 )
	{
		const uint8_t *j0 = tj + h * SEXTET_DW, *j1 = j0 + SEXTET_DW;
		__m128i a00 = _mm_setzero_si128(), a01 = a00, a10 = a00, a11 = a00, a20 = a00, a21 = a00, a30 = a00, a31 = a00, s;
		for( k=0; k<SEXTET_DW; k+=16 )
		{
			__m128i y0 = _mm_load_si128( (const __m128i *)( j0 + k ) ), y1 = _mm_load_si128( (const __m128i *)( j1 + k ) );
			__m128i x0 = _mm_load_si128( (const __m128i *)( i0 + k ) ), x1 = _mm_load_si128( (const __m128i *)( i1 + k ) );
			__m128i x2 = _mm_load_si128( (const __m128i *)( i2 + k ) ), x3 = _mm_load_si128( (const __m128i *)( i3 + k ) );
			a00 = _mm_add_epi64( a00, _mm_sad_epu8( x0, y0 ) ); a01 = _mm_add_epi64( a01, _mm_sad_epu8( x0, y1 ) );
			a10 = _mm_add_epi64( a10, _mm_sad_epu8( x1, y0 ) ); a11 = _mm_add_epi64( a11, _mm_sad_epu8( x1, y1 ) );
			a20 = _mm_add_epi64( a20, _mm_sad_epu8( x2, y0 ) ); a21 = _mm_add_epi64( a21, _mm_sad_epu8( x2, y1 ) );
			a30 = _mm_add_epi64( a30, _mm_sad_epu8( x3, y0 ) ); a31 = _mm_add_epi64( a31, _mm_sad_epu8( x3, y1 ) );
		}
#define SEXTET_HS128( a ) ( s = _mm_add_epi64( (a), _mm_unpackhi_epi64( (a), (a) ) ), (uint32_t)_mm_cvtsi128_si32( s ) )
		d[0+h] = SEXTET_HS128( a00 );  d[1+h] = SEXTET_HS128( a01 );
		d[4+h] = SEXTET_HS128( a10 );  d[5+h] = SEXTET_HS128( a11 );
		d[8+h] = SEXTET_HS128( a20 );  d[9+h] = SEXTET_HS128( a21 );
		d[12+h] = SEXTET_HS128( a30 ); d[13+h] = SEXTET_HS128( a31 );
#undef SEXTET_HS128
	}
}

// the same with 256-bit vectors (vpsadbw)
__attribute__((target("avx2")))
static void sextet_absdiff4x4_avx2( const uint8_t *ti, const uint8_t *tj, uint32_t *d )
{
	const uint8_t *i0 = ti, *i1 = ti + SEXTET_DW, *i2 = ti + 2 * SEXTET_DW, *i3 = ti + 3 * SEXTET_DW;
	int h, k;
	for( h=0; h<4; h+=2 )
	{
		const uint8_t *j0 = tj + h * SEXTET_DW, *j1 = j0 + SEXTET_DW;
		__m256i a00 = _mm256_setzero_si256(), a01 = a00, a10 = a00, a11 = a00, a20 = a00, a21 = a00, a30 = a00, a31 = a00;
		__m128i s;
		for( k=0; k<SEXTET_DW; k+=32 )
		{
			__m256i y0 = _mm256_load_si256( (const __m256i *)( j0 + k ) ), y1 = _mm256_load_si256( (const __m256i *)( j1 + k ) );
			__m256i x0 = _mm256_load_si256( (const __m256i *)( i0 + k ) ), x1 = _mm256_load_si256( (const __m256i *)( i1 + k ) );
			__m256i x2 = _mm256_load_si256( (const __m256i *)( i2 + k ) ), x3 = _mm256_load_si256( (const __m256i *)( i3 + k ) );
			a00 = _mm256_add_epi64( a00, _mm256_sad_epu8( x0, y0 ) ); a01 = _mm256_add_epi64( a01, _mm256_sad_epu8( x0, y1 ) );
			a10 = _mm256_add_epi64( a10, _mm256_sad_epu8( x1, y0 ) ); a11 = _mm256_add_epi64( a11, _mm256_sad_epu8( x1, y1 ) );
			a20 = _mm256_add_epi64( a20, _mm256_sad_epu8( x2, y0 ) ); a21 = _mm256_add_epi64( a21, _mm256_sad_epu8( x2, y1 ) );
			a30 = _mm256_add_epi64( a30, _mm256_sad_epu8( x3, y0 ) ); a31 = _mm256_add_epi64( a31, _mm256_sad_epu8( x3, y1 ) );
		}
#define SEXTET_HS256( a ) ( s = _mm_add_epi64( _mm256_castsi256_si128( a ), _mm256_extracti128_si256( (a), 1 ) ), s = _mm_add_epi64( s, _mm_unpackhi_epi64( s, s ) ), (uint32_t)_mm_cvtsi128_si32( s ) )
		d[0+h] = SEXTET_HS256( a00 );  d[1+h] = SEXTET_HS256( a01 );
		d[4+h] = SEXTET_HS256( a10 );  d[5+h] = SEXTET_HS256( a11 );
		d[8+h] = SEXTET_HS256( a20 );  d[9+h] = SEXTET_HS256( a21 );
		d[12+h] = SEXTET_HS256( a30 ); d[13+h] = SEXTET_HS256( a31 );
#undef SEXTET_HS256
	}
}
#endif

// lenfac of the conversion loop of the original code (below the pair loop in
// disttbfast()) for nogaplen[i] = li, nogaplen[j] = lj; same statements
static double sextet_lenfac( int li, int lj )
{
	double longer, shorter;
	double lenfac;

	if( li > lj )
	{
		longer=(double)li;
		shorter=(double)lj;
	}
	else
	{
		longer=(double)lj;
		shorter=(double)li;
	}
	lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
	return( lenfac );
}

// Stores the count of pair (i,j), j >= i.  With s->convert, j > i gets the
// distance, computed by the same statements as the conversion loop of the
// original code, where mtx[i][j-i] holds (double)cnt and mtx[i][0], mtx[j][0]
// hold the diagonal counts.  lenfac depends only on the two lengths; when they
// span at most SEXTET_LFMAX values it is read from a table filled by
// sextet_lenfac() (the same double), which saves 3 of the 4 divisions.
static void sextet_store( sextetfast_t *s, int i, int j, int cnt )
{
	double lenfac;
	double bunbo;
	double **mtx = s->mtx;

	mtx[i][j-i] = (double)cnt;
	if( j == i || !s->convert ) return;

	if( s->lftab )
		lenfac = s->lftab[( s->nogaplen[i] - s->lflo ) * s->lfn + ( s->nogaplen[j] - s->lflo )];
	else
		lenfac = sextet_lenfac( s->nogaplen[i], s->nogaplen[j] );
	bunbo = MIN( s->self[i], s->self[j] );
	if( bunbo == 0.0 )
		mtx[i][j-i] = 2.0; // 2013/Oct/17 -> 2bai
	else
		mtx[i][j-i] = ( 1.0 - mtx[i][j-i] / bunbo ) * lenfac * 2.0; // 2013/Oct/17 -> 2bai
}

// rows i0 .. i0+SEXTET_SR-1, every j >= i
static void sextet_densejob( sextetfast_t *s, int i0 )
{
	int n = s->n;
	int i1 = MIN( i0 + SEXTET_SR, n );
	int ib, jb, r, c, i, j, mcb, mcj;
	long nsimd = 0, nslow = 0;
#if SEXTET_NEON || FMSIMD_X86
	uint32_t d[16];
#endif

	for( jb=i0; jb<n; jb+=4 )
	{
		mcj = MAX( MAX( s->mc[jb], s->mc[jb+1] ), MAX( s->mc[jb+2], s->mc[jb+3] ) );
		for( ib=i0; ib<i1 && ib<=jb; ib+=4 )
		{
			mcb = MAX( mcj, MAX( MAX( s->mc[ib], s->mc[ib+1] ), MAX( s->mc[ib+2], s->mc[ib+3] ) ) );
#if SEXTET_NEON
			if( mcb <= 127 && s->simd >= FMSIMD_BASE )
			{
				sextet_absdiff4x4( s->tab + (size_t)ib * SEXTET_DW, s->tab + (size_t)jb * SEXTET_DW, mcb, d );
				for( r=0; r<4; r++ ) for( c=0; c<4; c++ )
				{
					i = ib + r; j = jb + c;
					if( i < n && j < n && j >= i )
						sextet_store( s, i, j, ( s->tot[i] + s->tot[j] - (int)d[r*4+c] ) / 2 );
				}
				nsimd++;
				continue;
			}
#elif FMSIMD_X86
			// psadbw takes any byte values, so there is no limit on the counts; tuples
			// where both rows have a saturated count (>= 255) get the overflow correction
			if( s->simd >= FMSIMD_BASE )
			{
				if( s->simd >= FMSIMD_AVX2 ) sextet_absdiff4x4_avx2( s->tab + (size_t)ib * SEXTET_DW, s->tab + (size_t)jb * SEXTET_DW, d );
				else sextet_absdiff4x4_sse2( s->tab + (size_t)ib * SEXTET_DW, s->tab + (size_t)jb * SEXTET_DW, d );
				for( r=0; r<4; r++ ) for( c=0; c<4; c++ )
				{
					i = ib + r; j = jb + c;
					if( i < n && j < n && j >= i )
						sextet_store( s, i, j, ( s->tot[i] + s->tot[j] - (int)d[r*4+c] ) / 2 + ( ( s->ovn[i] && s->ovn[j] ) ? sextet_ovcorr( s, i, j ) : 0 ) );
				}
				nsimd++;
				continue;
			}
#endif
			(void)mcb;
			nslow++;
			for( r=0; r<4; r++ ) for( c=0; c<4; c++ )
			{
				i = ib + r; j = jb + c;
				if( i < n && j < n && j >= i )
					sextet_store( s, i, j, sextet_pairslow( s, i, j ) );
			}
		}
	}
	if( s->timing ) { __sync_fetch_and_add( &s->nsimd, nsimd ); __sync_fetch_and_add( &s->nslow, nslow ); }
}

static void sextet_listrow( sextetfast_t *s, int i, uint16_t *tab )
{
	int n = s->n;
	int j, k, m;
	uint16_t *lp = s->lp + s->off[i], *lc = s->lc + s->off[i];

	m = s->ln[i];
	for( k=0; k<m; k++ ) tab[lp[k]] = lc[k];
	for( j=i; j<n; j++ )
	{
		const uint16_t *p = s->lp + s->off[j], *c = s->lc + s->off[j];
		unsigned int v0 = 0, v1 = 0, a, b;
		int mj = s->ln[j];
		for( k=0; k+1<mj; k+=2 )
		{
			a = tab[p[k]];   b = c[k];   v0 += ( a < b ) ? a : b;
			a = tab[p[k+1]]; b = c[k+1]; v1 += ( a < b ) ? a : b;
		}
		if( k < mj )
		{
			a = tab[p[k]]; b = c[k]; v0 += ( a < b ) ? a : b;
		}
		sextet_store( s, i, j, (int)( v0 + v1 ) );
	}
	for( k=0; k<m; k++ ) tab[lp[k]] = 0;
}

// Inverted index: for every tuple p of row i, every later row j that holds p
// gets min( c_i(p), c_j(p) ) added to acc[j]; acc[j] is then the count of
// pair (i,j), the same integer as the list loop computes.  The work is the
// sum over tuples of (rows holding it)^2/2 instead of n^2/2 * tuples per row.
static void sextet_invrow( sextetfast_t *s, int i, int *acc )
{
	int n = s->n, j, k, m = s->ln[i], cb = s->cb, self = 0;
	const uint16_t *lp = s->lp + s->off[i], *lc = s->lc + s->off[i];
	const uint32_t *pp = s->ppos + s->off[i];
	const uint32_t mask = ( 1u << cb ) - 1;
	for( k=0; k<m; k++ )
	{
		const uint32_t ci = lc[k];
		const uint32_t *e = s->post + pp[k] + 1, *end = s->post + s->poff[lp[k]+1];
		self += ci;
		for( ; e<end; e++ )
		{
			uint32_t cj = *e & mask;
			acc[*e >> cb] += ( ci < cj ) ? ci : cj;
		}
	}
	sextet_store( s, i, i, self );
	for( j=i+1; j<n; j++ )
	{
		sextet_store( s, i, j, acc[j] );
		acc[j] = 0;
	}
}

// Builds the inverted index from the tuple lists and decides between it and
// the list loop by their operation counts.  Returns 1 to use it.
static int sextet_buildinv( sextetfast_t *s, int force )
{
	int n = s->n, tw = s->tw, i, k, rb;
	size_t tot = s->off[n], q, *cur;
	double inv = 0.0, list = 0.0;
	uint32_t maxc = 0;

	for( rb=1; ( 1L << rb ) < n; rb++ ) ;
	s->cb = 32 - rb;
	if( s->cb < 2 || tot >= 0xffffffffUL ) return( 0 );
	s->poff = calloc( tw + 1, sizeof( size_t ) );
	if( !s->poff ) return( 0 );
	for( i=0; i<n; i++ )
	{
		const uint16_t *lp = s->lp + s->off[i], *lc = s->lc + s->off[i];
		for( k=0; k<s->ln[i]; k++ )
		{
			s->poff[lp[k]+1]++;
			if( lc[k] > maxc ) maxc = lc[k];
		}
		list += (double)s->ln[i] * ( i + 1 ); // row i's list is read for every row <= i
	}
	for( k=0; k<tw; k++ ) inv += (double)s->poff[k+1] * ( s->poff[k+1] + 1 ) / 2;
	s->listwork = list;
	s->invwork = inv;
	// a scattered add costs about twice a list step
	if( ( !force && 2.0 * inv >= list ) || (double)maxc >= (double)( 1u << s->cb ) )
	{
		free( s->poff ); s->poff = NULL;
		return( 0 );
	}
	for( k=0; k<tw; k++ ) s->poff[k+1] += s->poff[k];
	s->post = malloc( ( tot + 1 ) * sizeof( uint32_t ) );
	s->ppos = malloc( ( tot + 1 ) * sizeof( uint32_t ) );
	cur = malloc( tw * sizeof( size_t ) );
	if( !s->post || !s->ppos || !cur )
	{
		free( cur );
		free( s->post ); free( s->ppos ); free( s->poff );
		s->post = s->ppos = NULL; s->poff = NULL;
		return( 0 );
	}
	memcpy( cur, s->poff, tw * sizeof( size_t ) );
	for( i=0; i<n; i++ ) // ascending rows, so every posting list is in row order
	{
		const uint16_t *lp = s->lp + s->off[i], *lc = s->lc + s->off[i];
		uint32_t *pp = s->ppos + s->off[i];
		for( k=0; k<s->ln[i]; k++ )
		{
			q = cur[lp[k]]++;
			s->post[q] = ( (uint32_t)i << s->cb ) | lc[k];
			pp[k] = (uint32_t)q;
		}
	}
	free( cur );
	return( 1 );
}

// Reference for the check: the original count (commonsextet_p) and, with
// s->convert, the original conversion statements applied to it.
static void sextet_checkrow( sextetfast_t *s, int i )
{
	int j;
	double v;
	double longer, shorter;
	double lenfac;
	double bunbo;
	int *nogaplen = s->nogaplen;
	int *table1 = (int *)calloc( s->tw, sizeof( int ) );
	if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
	makecompositiontable_p( table1, s->pointt[i] );
	for( j=i; j<s->n; j++ )
	{
		v = (double)commonsextet_p( table1, s->pointt[j] );
		if( j == i && memcmp( &v, &s->self[i], sizeof( double ) ) )
		{
			reporterr( "\nsextetfast check FAILED: self count of %d: original=%.17g fast=%.17g\n", i, v, s->self[i] );
			exit( 1 );
		}
		if( j > i && s->convert )
		{
			if( nogaplen[i] > nogaplen[j] )
			{
				longer=(double)nogaplen[i];
				shorter=(double)nogaplen[j];
			}
			else
			{
				longer=(double)nogaplen[j];
				shorter=(double)nogaplen[i];
			}
			lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
			bunbo = MIN( s->self[i], s->self[j] );
			if( bunbo == 0.0 )
				v = 2.0;
			else
				v = ( 1.0 - v / bunbo ) * lenfac * 2.0;
		}
		if( memcmp( &v, &s->mtx[i][j-i], sizeof( double ) ) )
		{
			reporterr( "\nsextetfast check FAILED: i=%d j=%d original=%.17g fast=%.17g (kind %d, %s)\n", i, j, v, s->mtx[i][j-i], s->kind, s->convert ? "distance" : "count" );
			exit( 1 );
		}
	}
	free( table1 );
}

/* fast-mafft work/v9p, kind 4 (protein): the count of a pair is the size of the intersection of
   the two rows' sets of "elements" (p,k), k = 1..c(p) (min( c_i(p), c_j(p) ) of them per tuple p).
   The SEXTET_DMAX at most most frequent elements (each a bit of a dense mask per row, AND +
   popcount per pair) are chosen by an operation-count estimate; per tuple they are a prefix
   k <= kd(p) of its elements (an element's frequency does not grow with k; ties: smaller k first),
   so what is left of a tuple is  min( c_i - kd, c_j - kd )  for the rows with c > kd(p), added
   from an inverted index of those rows (as kind 3).  On homfam_rvp 128 elements leave 0.09 index
   steps per pair (kind 3: 67, the list loop: 88).  Exact integers throughout. */
#define SEXTET_KD 8      // occurrences per tuple that can be dense elements (k < SEXTET_KD)
#define SEXTET_DMAX 1024 // largest number of dense elements

static int sextet_u64cmp( const void *a, const void *b )
{
	uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
	return( x < y ? -1 : x > y ? 1 : 0 );
}

// Chooses the dense elements and builds the masks and the index of the rest.  Returns 1 to use
// kind 4 (force: whenever there is something to make dense).
static int sextet_builddense( sextetfast_t *s, int force )
{
	int n = s->n, tw = s->tw, i, k, p, c, e, m, K, bestK, rb, dw;
	int *f = NULL, *kd = NULL, *bitof = NULL;
	uint64_t *el = NULL;
	size_t q, nres, *cur = NULL;
	double pairs = (double)n * ( n - 1 ) / 2, sp, sp0, best, cost, list = 0.0;
	uint32_t maxc = 0;

	f = calloc( (size_t)tw * SEXTET_KD, sizeof( int ) );
	kd = calloc( tw, sizeof( int ) );
	if( !f || !kd ) goto fail;
	for( i=0; i<n; i++ )
	{
		const uint16_t *lp = s->lp + s->off[i], *lc = s->lc + s->off[i];
		for( k=0; k<s->ln[i]; k++ )
		{
			c = lc[k] < SEXTET_KD ? lc[k] : SEXTET_KD;
			for( e=0; e<c; e++ ) f[lp[k]*SEXTET_KD+e]++;
		}
		list += (double)s->ln[i] * ( i + 1 );
	}
	for( m=0, q=0; q<(size_t)tw*SEXTET_KD; q++ ) if( q % SEXTET_KD < SEXTET_KD - 1 && f[q] >= 2 ) m++;
	el = malloc( ( m + 1 ) * sizeof( uint64_t ) );
	if( !el ) goto fail;
	for( m=0, q=0; q<(size_t)tw*SEXTET_KD; q++ ) if( q % SEXTET_KD < SEXTET_KD - 1 && f[q] >= 2 ) el[m++] = ( (uint64_t)( 0x7fffffff - f[q] ) << 32 ) | q; // frequency descending, then (p,k) ascending
	qsort( el, m, sizeof( uint64_t ), sextet_u64cmp );
	// operation counts: a mask word per pair per 64 elements; two per index step (as kind 3)
	for( sp=0.0, p=0; p<tw; p++ ) sp += (double)f[p*SEXTET_KD] * f[p*SEXTET_KD] / 2;
	sp0 = sp;
	best = 2.0 * sp;
	bestK = 0;
	for( K=1; K<=m && K<=SEXTET_DMAX; K++ )
	{
		q = el[K-1] & 0xffffffff;
		p = (int)( q / SEXTET_KD ); e = (int)( q % SEXTET_KD ); // e == kd[p]
		sp += ( (double)f[q+1] * f[q+1] - (double)f[q] * f[q] ) / 2;
		kd[p]++;
		if( K % 64 == 0 && ( cost = pairs * ( K / 64 ) + 2.0 * sp ) < best )
		{
			best = cost;
			bestK = K;
		}
	}
	s->listwork = list;
	s->invwork = sp0;
	s->densework = best;
	if( bestK == 0 || ( !force && best >= 0.8 * MIN( list, 2.0 * sp0 ) ) ) goto fail;

	// the dense elements: bit b = their rank
	s->dk = bestK;
	s->dw = dw = bestK / 64;
	bitof = malloc( (size_t)tw * SEXTET_KD * sizeof( int ) );
	if( !bitof ) goto fail;
	memset( kd, 0, tw * sizeof( int ) );
	for( K=0; K<bestK; K++ )
	{
		q = el[K] & 0xffffffff;
		p = (int)( q / SEXTET_KD );
		bitof[q] = K;
		kd[p]++;
	}
	if( posix_memalign( (void **)&s->dmask, 64, (size_t)s->npad * dw * sizeof( uint64_t ) ) ) { s->dmask = NULL; goto fail; }
	memset( s->dmask, 0, (size_t)s->npad * dw * sizeof( uint64_t ) );
	// the rest: per tuple the rows with c > kd(p), count c - kd(p)
	s->poff = calloc( tw + 1, sizeof( size_t ) );
	s->rstart = calloc( n + 1, sizeof( size_t ) );
	if( !s->poff || !s->rstart ) goto fail;
	for( nres=0, i=0; i<n; i++ )
	{
		const uint16_t *lp = s->lp + s->off[i], *lc = s->lc + s->off[i];
		uint64_t *mk = s->dmask + (size_t)i * dw;
		for( k=0; k<s->ln[i]; k++ )
		{
			p = lp[k];
			c = lc[k];
			for( e=0; e<c && e<kd[p]; e++ ) mk[bitof[p*SEXTET_KD+e]>>6] |= 1ULL << ( bitof[p*SEXTET_KD+e] & 63 );
			if( c > kd[p] )
			{
				s->poff[p+1]++;
				if( (uint32_t)( c - kd[p] ) > maxc ) maxc = c - kd[p];
				nres++;
			}
		}
		s->rstart[i+1] = nres;
	}
	for( rb=1; ( 1L << rb ) < n; rb++ ) ;
	s->cb = 32 - rb;
	if( s->cb < 2 || nres >= 0xffffffffUL || (double)maxc >= (double)( 1u << s->cb ) ) goto fail;
	for( p=0; p<tw; p++ ) s->poff[p+1] += s->poff[p];
	s->post = malloc( ( nres + 1 ) * sizeof( uint32_t ) );
	s->rpp = malloc( ( nres + 1 ) * sizeof( uint32_t ) );
	s->rtup = malloc( ( nres + 1 ) * sizeof( uint16_t ) );
	cur = malloc( tw * sizeof( size_t ) );
	if( !s->post || !s->rpp || !s->rtup || !cur ) goto fail;
	memcpy( cur, s->poff, tw * sizeof( size_t ) );
	for( nres=0, i=0; i<n; i++ ) // ascending rows, so every posting list is in row order
	{
		const uint16_t *lp = s->lp + s->off[i], *lc = s->lc + s->off[i];
		for( k=0; k<s->ln[i]; k++ )
		{
			p = lp[k];
			c = lc[k];
			if( c > kd[p] )
			{
				q = cur[p]++;
				s->post[q] = ( (uint32_t)i << s->cb ) | (uint32_t)( c - kd[p] );
				s->rpp[nres] = (uint32_t)q;
				s->rtup[nres] = (uint16_t)p;
				nres++;
			}
		}
	}
	s->nres = nres;
	free( cur ); free( f ); free( kd ); free( bitof ); free( el );
	return( 1 );
fail:
	free( cur ); free( f ); free( kd ); free( bitof ); free( el );
	free( s->dmask ); s->dmask = NULL;
	free( s->poff ); s->poff = NULL;
	free( s->rstart ); s->rstart = NULL;
	free( s->post ); s->post = NULL;
	free( s->rpp ); s->rpp = NULL;
	free( s->rtup ); s->rtup = NULL;
	s->dk = s->dw = 0;
	return( 0 );
}

// kind 4, row i: the rest from the index into acc[], then every pair's count and distance
static inline __attribute__((always_inline)) void sextet_denserowbody( sextetfast_t *s, int i, int *acc )
{
	int n = s->n, j, w, dw = s->dw, cb = s->cb, c;
	const uint32_t mask = ( 1u << cb ) - 1;
	const uint64_t *bi = s->dmask + (size_t)i * dw, *bj;
	size_t k;
	for( k=s->rstart[i]; k<s->rstart[i+1]; k++ )
	{
		const uint32_t *e = s->post + s->rpp[k], *end = s->post + s->poff[s->rtup[k]+1];
		const uint32_t ci = *e & mask;
		for( e++; e<end; e++ )
		{
			uint32_t cj = *e & mask;
			acc[*e >> cb] += ( ci < cj ) ? ci : cj;
		}
	}
	sextet_store( s, i, i, (int)s->self[i] ); // the number of tuples, as commonsextet_p( table_i, pointt[i] )
	bj = bi + dw;
	if( dw == 2 )
	{
		const uint64_t b0 = bi[0], b1 = bi[1];
		for( j=i+1; j<n; j++, bj+=2 )
		{
			c = acc[j] + __builtin_popcountll( b0 & bj[0] ) + __builtin_popcountll( b1 & bj[1] );
			acc[j] = 0;
			sextet_store( s, i, j, c );
		}
		return;
	}
	for( j=i+1; j<n; j++, bj+=dw )
	{
		c = acc[j];
		acc[j] = 0;
		for( w=0; w<dw; w++ ) c += __builtin_popcountll( bi[w] & bj[w] );
		sextet_store( s, i, j, c );
	}
}

static void sextet_denserow( sextetfast_t *s, int i, int *acc )
{
	sextet_denserowbody( s, i, acc );
}

#if FMSIMD_X86
__attribute__((target("popcnt")))
static void sextet_denserow_popcnt( sextetfast_t *s, int i, int *acc )
{
	sextet_denserowbody( s, i, acc );
}
#endif


static void *sextet_thread( void *arg )
{
	sextetfast_t *s = (sextetfast_t *)arg;
	int i, i0, i1, step;
	int *cnt, *acc = NULL;
	uint16_t *tab = NULL;

	if( s->phase == 0 )
	{
		cnt = (int *)calloc( s->tw, sizeof( int ) );
		if( !cnt ) ErrorExit( "Cannot allocate cnt\n" );
		while( ( i0 = sextet_claim( s, 64 ) ) < s->n )
		{
			i1 = MIN( i0 + 64, s->n );
			for( i=i0; i<i1; i++ ) sextet_buildrow( s, i, cnt );
		}
		free( cnt );
		return( NULL );
	}

	step = ( s->kind == 1 ) ? SEXTET_SR : SEXTET_LR;
	if( s->kind == 2 )
	{
		tab = (uint16_t *)calloc( s->tw, sizeof( uint16_t ) );
		if( !tab ) ErrorExit( "Cannot allocate tab\n" );
	}
	if( s->kind == 3 || s->kind == 4 )
	{
		acc = (int *)calloc( s->n, sizeof( int ) );
		if( !acc ) ErrorExit( "Cannot allocate acc\n" );
	}
	while( ( i0 = sextet_claim( s, step ) ) < s->n )
	{
		i1 = MIN( i0 + step, s->n );
		if( i0 % 100 == 0 || i0 / 100 != ( i1 - 1 ) / 100 ) reporterr( "\r% 5d / %d", i1, s->n );
		if( s->kind == 1 ) sextet_densejob( s, i0 );
		else if( s->kind == 3 ) for( i=i0; i<i1; i++ ) sextet_invrow( s, i, acc );
#if FMSIMD_X86
		else if( s->kind == 4 && s->popcnt ) for( i=i0; i<i1; i++ ) sextet_denserow_popcnt( s, i, acc );
#endif
		else if( s->kind == 4 ) for( i=i0; i<i1; i++ ) sextet_denserow( s, i, acc );
		else for( i=i0; i<i1; i++ ) sextet_listrow( s, i, tab );
		if( s->check ) for( i=i0; i<i1; i++ ) sextet_checkrow( s, i );
	}
	if( s->check ) commonsextet_p( NULL, NULL );
	free( tab );
	free( acc );
	return( NULL );
}


static void sextetfast_free( sextetfast_t *s )
{
	int i;
	if( s->ovp ) for( i=0; i<s->npad; i++ )
	{
		free( s->ovp[i] );
		free( s->ovc[i] );
	}
	free( s->tab ); free( s->tot ); free( s->mc ); free( s->ovn ); free( s->ovp ); free( s->ovc );
	free( s->ln ); free( s->off ); free( s->lp ); free( s->lc );
	free( s->post ); free( s->poff ); free( s->ppos );
	free( s->dmask ); free( s->rstart ); free( s->rpp ); free( s->rtup ); // work/v9p
	free( s->self );
	free( s->lftab );
}

// Fills mtx[i][j-i] (j >= i) with the common tuple counts that the original
// pair loop stores and returns 1, or also converts them to distances, exactly
// as the conversion loop that follows it, and returns 2 (then that loop must
// be skipped).  Returns 0, with mtx untouched, when the original code must be used.
static int sextetfast_run( int njob, int **pointt, double **mtx, int *nogaplen )
{
	sextetfast_t s;
	const char *mode = getenv( "MAFFT_SEXTET" );
	int i, m, t, nthr = 1;
	int *p;
	double t0, t1, t2;
#ifdef enablemultithread
	pthread_t *handle = NULL;
	pthread_mutex_t mutex;
#endif

	if( mode && strstr( mode, "ref" ) ) return( 0 );
	if( njob < 1 ) return( 0 );
	memset( &s, 0, sizeof( s ) );
	s.n = njob;
	s.npad = ( njob + 3 ) & ~3;
	s.pointt = pointt;
	s.mtx = mtx;
	s.nogaplen = nogaplen;
	s.convert = !( mode && strstr( mode, "counts" ) );
	s.tw = tsize;
	s.check = ( mode && strstr( mode, "check" ) );
	s.scalar = ( mode && strstr( mode, "scalar" ) );
	s.timing = ( mode && strstr( mode, "time" ) );
	s.simd = s.scalar ? FMSIMD_SCALAR : fmsimd_level();
	if( tsize == SEXTET_DW && !( mode && strstr( mode, "list" ) ) ) s.kind = 1;
	else if( tsize == SEXTET_DW || tsize == 46656 ) s.kind = 2;
	else return( 0 );

	t0 = sextet_now();
	s.self = calloc( s.npad, sizeof( double ) );
	if( !s.self ) return( 0 );
	if( s.convert && !( mode && strstr( mode, "nolftab" ) ) )
	{
		int lo = nogaplen[0], hi = nogaplen[0], a, b;
		for( i=1; i<njob; i++ )
		{
			if( nogaplen[i] < lo ) lo = nogaplen[i];
			if( nogaplen[i] > hi ) hi = nogaplen[i];
		}
		if( hi - lo < SEXTET_LFMAX && ( s.lftab = malloc( sizeof( double ) * ( hi - lo + 1 ) * ( hi - lo + 1 ) ) ) )
		{
			s.lflo = lo;
			s.lfn = hi - lo + 1;
			for( a=0; a<s.lfn; a++ ) for( b=0; b<s.lfn; b++ ) s.lftab[a*s.lfn+b] = sextet_lenfac( lo + a, lo + b );
		}
	}
	if( s.kind == 1 )
	{
		if( posix_memalign( (void **)&s.tab, 128, (size_t)s.npad * SEXTET_DW ) )
		{
			s.tab = NULL;
			sextetfast_free( &s );
			return( 0 );
		}
		memset( s.tab, 0, (size_t)s.npad * SEXTET_DW );
		s.tot = calloc( s.npad, sizeof( int ) );
		s.mc = calloc( s.npad, 1 );
		s.ovn = calloc( s.npad, sizeof( int ) );
		s.ovp = calloc( s.npad, sizeof( int * ) );
		s.ovc = calloc( s.npad, sizeof( int * ) );
		if( !s.tot || !s.mc || !s.ovn || !s.ovp || !s.ovc )
		{
			sextetfast_free( &s );
			return( 0 );
		}
	}
	else
	{
		s.ln = calloc( njob, sizeof( int ) );
		s.off = calloc( njob + 1, sizeof( size_t ) );
		if( !s.ln || !s.off )
		{
			sextetfast_free( &s );
			return( 0 );
		}
		for( i=0; i<njob; i++ )
		{
			p = pointt[i];
			m = 0;
			if( *p != -1 ) while( *p++ != END_OF_VEC ) m++;
			s.off[i+1] = s.off[i] + m;
		}
		s.lp = malloc( ( s.off[njob] + 1 ) * sizeof( uint16_t ) );
		s.lc = malloc( ( s.off[njob] + 1 ) * sizeof( uint16_t ) );
		if( !s.lp || !s.lc )
		{
			sextetfast_free( &s );
			return( 0 );
		}
	}

#ifdef enablemultithread
	if( nthreadpair > 1 )
	{
		nthr = nthreadpair;
		pthread_mutex_init( &mutex, NULL );
		s.mutex = &mutex;
		handle = calloc( nthr, sizeof( pthread_t ) );
	}
#endif
	t1 = t0;
	for( s.phase=0; s.phase<2; s.phase++ )
	{
		s.jobpos = 0;
#ifdef enablemultithread
		if( nthr > 1 )
		{
			for( t=0; t<nthr; t++ ) pthread_create( handle+t, NULL, sextet_thread, (void *)&s );
			for( t=0; t<nthr; t++ ) pthread_join( handle[t], NULL );
		}
		else
#endif
			sextet_thread( (void *)&s );
		if( s.phase == 0 )
		{
			if( s.bad ) break;
			if( s.kind == 2 && tsize == 46656 && !( mode && ( strstr( mode, "nodense" ) || strstr( mode, "noinv" ) ) ) && sextet_builddense( &s, mode && strstr( mode, "dense" ) && !strstr( mode, "nodense" ) ) ) // work/v9p
			{
				int pc = 0;
				s.kind = 4;
				fmsimd_cpu( &pc, NULL );
				s.popcnt = pc && s.simd >= FMSIMD_BASE;
			}
			else
			if( s.kind == 2 && !( mode && strstr( mode, "noinv" ) ) && sextet_buildinv( &s, mode && strstr( mode, "inv" ) && !strstr( mode, "noinv" ) ) ) s.kind = 3;
			t1 = sextet_now();
		}
	}
	t2 = sextet_now();
#ifdef enablemultithread
	if( nthr > 1 )
	{
		pthread_mutex_destroy( &mutex );
		free( handle );
	}
#endif
	sextetfast_free( &s );
	if( s.bad ) return( 0 ); // a list count did not fit; mtx was not written
	if( s.timing ) reporterr( "\nsextetfast: kind %d, simd level %d, %d threads, tables %.3f s, pairs %.3f s; dense 4x4 blocks: %ld simd, %ld plain; list work %.3g, inverted work %.3g, dense+rest work %.3g (%d dense elements, %zu rest entries)\n", s.kind, s.simd, nthr, t1 - t0, t2 - t1, s.nsimd, s.nslow, s.listwork, s.invwork, s.densework, s.dk, s.nres );
	if( s.check ) reporterr( "\nsextetfast check passed: kind %d, simd level %d, %d rows, every pair's %s bitwise identical\n", s.kind, s.simd, njob, s.convert ? "distance" : "count" );
	return( s.convert ? 2 : 1 );
}
// MAFFT_SEXTET_DUMP: the 6-mer distance matrix, row i = mtx[i][0..n-1-i] as raw doubles
static void sextet_dump( const char *fn, int n, double **mtx )
{
	int i;
	FILE *fp = fopen( fn, "wb" );
	if( !fp ) return;
	for( i=0; i<n; i++ ) fwrite( mtx[i], sizeof( double ), n - i, fp );
	fclose( fp );
}
/* END sextetfast kernel */

/* BEGIN v9p halfmtx */
/*
 * fast-mafft work/v9p: the half distance matrix in one anonymous mapping, rows pointing into it
 * (the original allocates every row with calloc, serially, which on glibc above the dynamic mmap
 * threshold means zeroing all of it in one thread: 1.8 s on homfam_sdr, 3.7 s on homfam_rvp, and
 * 4 KB pages, one TLB miss per strided read in the guide-tree step).  The mapping is zero-filled
 * by the kernel on first touch (by the threads that write the distances), with transparent huge
 * pages requested.  Only for the default guide-tree path (fixed_musclesupg_..._v9/_mt know the
 * slab and do not free its rows one by one); everything else uses AllocateFloatHalfMtx.
 *   MAFFT_V9SLAB=0   original allocation;  =1  never MADV_HUGEPAGE;  =2  always MADV_HUGEPAGE
 */
#include <sys/mman.h>
extern char *v9_slab_lo, *v9_slab_hi;
static double **v9_slab_mtx = NULL;
static size_t v9_slab_len = 0;

static double **v9_halfmtx_alloc( int n )
{
	const char *e = getenv( "MAFFT_V9SLAB" );
	size_t tot, off;
	double **m;
	char *base;
	int i;
	if( ( e && *e == '0' ) || v9_slab_mtx || n < 1 ) return( AllocateFloatHalfMtx( n ) );
	tot = (size_t)n * ( (size_t)n + 1 ) / 2 * sizeof( double );
	base = (char *)mmap( NULL, tot, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0 );
	if( base == MAP_FAILED ) return( AllocateFloatHalfMtx( n ) );
#ifdef MADV_HUGEPAGE
	// transparent huge pages: with the system mode "always" they come without asking (and without
	// the synchronous compaction that MADV_HUGEPAGE requests under the usual defrag=madvise); with
	// "madvise" they must be asked for
	if( !( e && *e == '1' ) )
	{
		char thp[128] = "";
		FILE *fp = fopen( "/sys/kernel/mm/transparent_hugepage/enabled", "r" );
		if( fp )
		{
			if( !fgets( thp, sizeof( thp ), fp ) ) thp[0] = 0;
			fclose( fp );
		}
		if( strstr( thp, "[madvise]" ) || ( e && *e == '2' ) ) madvise( base, tot, MADV_HUGEPAGE );
	}
#endif
	m = (double **)calloc( (size_t)n + 1, sizeof( double * ) );
	if( m == NULL )
	{
		munmap( base, tot );
		return( AllocateFloatHalfMtx( n ) );
	}
	for( off=0, i=0; i<n; i++ )
	{
		m[i] = (double *)( base + off );
		off += (size_t)( n - i ) * sizeof( double );
	}
	m[n] = NULL;
	v9_slab_lo = base;
	v9_slab_hi = base + tot;
	v9_slab_len = tot;
	v9_slab_mtx = m;
	return( m );
}

static void v9_halfmtx_free( double **m, int n )
{
	int i;
	if( m == NULL ) return;
	if( m != v9_slab_mtx )
	{
		FreeFloatHalfMtx( m, n );
		return;
	}
	for( i=0; i<n; i++ ) if( m[i] && !( (char *)m[i] >= v9_slab_lo && (char *)m[i] < v9_slab_hi ) ) free( m[i] );
	free( m );
	munmap( v9_slab_lo, v9_slab_len );
	v9_slab_lo = v9_slab_hi = NULL;
	v9_slab_len = 0;
	v9_slab_mtx = NULL;
}

/* The guide tree from the half matrix: fixed_musclesupg_..._v9 (mltaln9.c, the same tree as the
   _mt version and the original), or with MAFFT_UPG=old the _mt version.  MAFFT_UPG=check also
   builds the tree with the _mt version from a copy of the matrix and exits on any difference
   in topol, len, dep or the matrix left behind (needs memory for a second matrix). */
static void v9_upgma( int n, double **mtx, int ***topol, double **len, Treedep *dep, int progressout, int efffree, int nthread )
{
	const char *mode = getenv( "MAFFT_UPG" );
	double **cp;
	int ***topol2;
	double **len2;
	Treedep *dep2;
	int i, k, x, bad = 0;
	if( mode && strstr( mode, "old" ) )
	{
		fixed_musclesupg_double_realloc_nobk_halfmtx_memsave_mt( n, mtx, topol, len, dep, progressout, efffree, nthread );
		return;
	}
	if( !( mode && strstr( mode, "check" ) ) || n < 2 )
	{
		fixed_musclesupg_double_realloc_nobk_halfmtx_memsave_v9( n, mtx, topol, len, dep, progressout, efffree, nthread );
		return;
	}
	cp = (double **)calloc( n + 1, sizeof( double * ) );
	for( i=0; i<n; i++ )
	{
		cp[i] = (double *)malloc( ( n - i ) * sizeof( double ) );
		if( cp[i] == NULL ) ErrorExit( "MAFFT_UPG=check: no memory for a copy of the matrix\n" );
		memcpy( cp[i], mtx[i], ( n - i ) * sizeof( double ) );
	}
	topol2 = AllocateIntCub( n, 2, 0 );
	len2 = AllocateFloatMtx( n, 2 );
	dep2 = (Treedep *)calloc( n, sizeof( Treedep ) );
	fixed_musclesupg_double_realloc_nobk_halfmtx_memsave_mt( n, cp, topol2, len2, dep ? dep2 : NULL, 0, efffree, nthread );
	fixed_musclesupg_double_realloc_nobk_halfmtx_memsave_v9( n, mtx, topol, len, dep, progressout, efffree, nthread );
	for( k=0; k<n-1 && !bad; k++ )
	{
		for( x=0; x<2; x++ )
		{
			if( topol[k][x][0] != topol2[k][x][0] || topol[k][x][1] != topol2[k][x][1] ) bad = 1;
			if( memcmp( &len[k][x], &len2[k][x], sizeof( double ) ) ) bad = 2;
		}
		if( dep && ( dep[k].child0 != dep2[k].child0 || dep[k].child1 != dep2[k].child1 || memcmp( &dep[k].distfromtip, &dep2[k].distfromtip, sizeof( double ) ) ) ) bad = 3;
		if( bad ) reporterr( "\nUPGMA check FAILED at step %d (%s): v9 %d,%d len %.17g,%.17g; _mt %d,%d len %.17g,%.17g\n", k, bad == 1 ? "topol" : bad == 2 ? "len" : "dep", topol[k][0][0], topol[k][1][0], len[k][0], len[k][1], topol2[k][0][0], topol2[k][1][0], len2[k][0], len2[k][1] );
	}
	for( i=0; i<n && !bad; i++ )
	{
		if( ( mtx[i] == NULL ) != ( cp[i] == NULL ) ) bad = 4;
		else if( mtx[i] && memcmp( mtx[i], cp[i], ( n - i ) * sizeof( double ) ) ) bad = 5;
		if( bad ) reporterr( "\nUPGMA check FAILED: matrix row %d differs (%d)\n", i, bad );
	}
	if( bad ) exit( 1 );
	reporterr( "\nUPGMA check passed: n=%d efffree=%d, %d merges, topol/len/dep and the matrix identical to the _mt version\n", n, efffree, n - 1 );
	for( i=0; i<n; i++ ) if( cp[i] ) free( cp[i] );
	free( cp );
	FreeIntCub( topol2 );
	FreeFloatMtx( len2 );
	free( dep2 );
}
/* END v9p halfmtx */

static void recountpositions( ExtAnch *pairanch, int n1, int n2, char **seq1, char **seq2 ) // loop no junban kentou
{
	int i, j, k, len, pos;
	int *map;

	len = strlen( seq1[0] )+1;
	map = calloc( sizeof( int ), len );

	for( k=0; k<n1; k++ )
	{
		pos = 0;
		for( i=0; i<len; i++ )
		{
			if( seq1[k][i] != '-' )
			{
				map[pos] = i;
				pos++;
			}
		}

		for( j=0; pairanch[j].i>-1; j++ )
		{
			if( pairanch[j].i == k ) 
			{
//				reporterr( "pairanch[%d].endi: %d->%d\n", j, pairanch[j].endi, map[pairanch[j].endi] );
				pairanch[j].starti = map[pairanch[j].starti];
				pairanch[j].endi = map[pairanch[j].endi];
			}
		}
	}
	free( map );

	len = strlen( seq2[0] )+1;
	map = calloc( sizeof( int ), len );
	for( k=0; k<n2; k++ )
	{
		pos = 0;
		for( i=0; i<len; i++ )
		{
			if( seq2[k][i] != '-' )
			{
				map[pos] = i;
				pos++;
			}
		}
		for( j=0; pairanch[j].i>-1; j++ )
		{
			if( pairanch[j].j == k )
			{
//				reporterr( "pairanch[%d].endj: %d->%d\n", j, pairanch[j].endj, map[pairanch[j].endj] );
				pairanch[j].startj = map[pairanch[j].startj];
				pairanch[j].endj = map[pairanch[j].endj];
			}
		}
	}
	free( map );
}

static int anchidcomp( const void *p, const void *q )
{
	if ( ((ExtAnch *)q)->i != ((ExtAnch *)p)->i )
		return ((ExtAnch *)p)->i - ((ExtAnch *)q)->i;
	return ((ExtAnch *)p)->j - ((ExtAnch *)q)->j;
}

static int anchcomp( const void *p, const void *q )
{
	if ( ((ExtAnch *)q)->starti != ((ExtAnch *)p)->starti )
		return ((ExtAnch *)p)->starti - ((ExtAnch *)q)->starti;
	return (int)((void *)p - (void *)q);
}

static int anchscorecomp( const void *p, const void *q )
{
	if ( ((ExtAnch *)q)->score != ((ExtAnch *)p)->score )
		return ((ExtAnch *)q)->score - ((ExtAnch *)p)->score;
	return (int)((void *)q - (void *)p);
}


static void indexanchors( ExtAnch *a, int **idx )
{
	int n;
	for( n=0; a[n].i>-1; n++ )
		;

	qsort( a, n, sizeof( ExtAnch ), anchidcomp );


	for( n=0; a[n].i>-1; n++ )
	{
//		reporterr( "%d, %dx%d, %d-%d x %d-%d\n", n, a[n].i, a[n].j, a[n].starti, a[n].endi, a[n].startj, a[n].endj );
		if( idx[a[n].i][a[n].j] == -1 ) idx[a[n].i][a[n].j] = n;
	}
#if 0
	int m;
	for( n=0; n<njob; n++ ) for( m=n+1; m<njob; m++ )
		reporterr( "%dx%d -> %d\n", n, m, idx[n][m] );
	exit( 1 );
#endif
}

static void checkanchors_strongestfirst( ExtAnch *a, int s )
{
	int p, q;
	qsort( a, s, sizeof( ExtAnch ), anchscorecomp );

	for( p=0; a[p].i>-1; p++ )
	{
		if( a[p].starti == -1 ) continue;

//		reporterr( "P score=%d, %d-%d, %d-%d\n", a[p].score, a[p].starti, a[p].endi, a[p].startj, a[p].endj );
		for( q=p+1; a[q].i>-1; q++ )
		{
			if( a[q].starti == -1 ) continue;
//			reporterr( "Q score=%d, %d-%d, %d-%d\n", a[q].score, a[q].starti, a[q].endi, a[q].startj, a[q].endj );


			if( a[p].endi < a[q].starti && a[p].endj < a[q].startj ) 
			{
//				reporterr( "consistent\n" );
				;
			}
			else if( a[p].endi == a[q].starti && a[p].endj < a[q].startj && a[q].starti<a[q].endi ) 
			{
				a[q].starti += 1; // 1 zai overlap
			}
			else if( a[p].endi < a[q].starti && a[p].endj == a[q].startj && a[q].startj<a[q].endj ) 
			{
				a[q].startj += 1; // 1 zai overlap
			}
			else if( a[q].endi < a[p].starti && a[q].endj < a[p].startj )
			{
//				reporterr( "consistent\n" );
				;
			}
			else if( a[q].endi == a[p].starti && a[q].endj < a[p].startj && a[q].starti<a[q].endi ) // bug in v7.442
			{
				a[q].endi -= 1; // 1 zai overlap
			}
			else if( a[q].endi < a[p].starti && a[q].endj == a[p].startj && a[q].startj<a[q].endj )
			{
				a[q].endj -= 1; // 1 zai overlap
			}
			else 
			{
//				reporterr( "INconsistent\n" );
				a[q].starti = a[q].startj = a[q].startj = a[q].endj = -1;
			}
		}
		if( p % 1000 == 0 ) reporterr( "%d/%d\r", p, s );
	}

	qsort( a, s, sizeof( ExtAnch ), anchcomp );
#if 0
	reporterr( "after filtering and sorting\n" );
	for( p=0; a[p].i>-1; p++ )
	{
		reporterr( "a[%d].starti,j=%d,%d, score=%d\n", p, a[p].starti, a[p].startj, a[p].score );
	}
#endif
}

static void	pickpairanch( ExtAnch **pairanch, ExtAnch *extanch, int **anchindex, int n1, int n2, int *m1, int *m2, char **seq1, char **seq2 ) // loop no junban wo kaeta hou ga iikamo
{
	int i, j, k, s;
	s = 0;
#if 0
	reporterr( "m1,m2=\n" );
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		reporterr( "%d,%d\n", m1[i], m2[j] );
	}
#endif
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
#if 1
		if( m1[i] < m2[j] ) 
		{
//			reporterr( "%dx%d, %dx%d -> jump to %d\n", i, j, m1[i], m2[j], anchindex[m1[i]][m2[j]] );
			k = anchindex[m1[i]][m2[j]];
			while( ( k!=-1 ) && ( extanch[k].i == m1[i] && extanch[k].j == m2[j] ) )
			{
				s++;
				k++;
			}
		}
		else
		{
//			reporterr( "%dx%d, %dx%d -> jump to %d\n", j, i, m1[i], m2[j], anchindex[m2[j]][m1[i]] );
			k = anchindex[m2[j]][m1[i]];
			while( ( k!=-1 ) && ( extanch[k].i == m2[j] && extanch[k].j == m1[i] ) )
			{
				s++;
				k++;
			}
		}
#else
		k = 0;
		while( extanch[k].i > -1 ) // kanari muda
		{
			//reporterr( "m1[i],m2[j]=%d,%d ? extanch[k].i,j=%d,%d k=%d\n", m1[i], m2[j], extanch[k].i, extanch[k].j, k );
			if( ( extanch[k].i == m1[i] && extanch[k].j == m2[j] ) || ( extanch[k].i == m2[j] && extanch[k].j == m1[i] ) ) 
			{
				//reporterr( "hit, extanch[k].startj=%d\n", extanch[k].startj );
				s++;
			}
			k++;
		}
#endif
	}
	*pairanch = calloc( sizeof( ExtAnch ), s+1 );

	s = 0;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
#if 1
		if( m1[i] < m2[j] ) 
		{
			k = anchindex[m1[i]][m2[j]];
			while( ( k!=-1 ) && ( extanch[k].i == m1[i] && extanch[k].j == m2[j] ) )
			{
//				if( extanch[k].starti + 1 < extanch[k].endi )
				{
					(*pairanch)[s].i = i;
					(*pairanch)[s].j = j;
					(*pairanch)[s].starti = extanch[k].starti; // map mae
					(*pairanch)[s].endi = extanch[k].endi; // map mae
					(*pairanch)[s].startj = extanch[k].startj; // map mae 
					(*pairanch)[s].endj = extanch[k].endj; // map mae
					(*pairanch)[s].score = extanch[k].score;
					s++;
				}
				k++;
			}
		}
		else
		{
			k = anchindex[m2[j]][m1[i]];
			while( ( k!=-1 ) && ( extanch[k].i == m2[j] && extanch[k].j == m1[i] ) )
			{
//				if( extanch[k].starti + 1 < extanch[k].endi )
				{
					(*pairanch)[s].i = i;
					(*pairanch)[s].j = j;
					(*pairanch)[s].starti = extanch[k].startj; // map mae
					(*pairanch)[s].endi = extanch[k].endj; // map mae
					(*pairanch)[s].startj = extanch[k].starti; // map mae 
					(*pairanch)[s].endj = extanch[k].endi; // map mae
					(*pairanch)[s].score = extanch[k].score;
					s++;
				}
				k++;
			}
		}
#else
		k = 0;
		while( extanch[k].i > -1 ) // kanari muda
		{
			if( extanch[k].i == m1[i] && extanch[k].j == m2[j] )
			{
				(*pairanch)[s].i = i;
				(*pairanch)[s].j = j;
				(*pairanch)[s].starti = extanch[k].starti; // map mae
				(*pairanch)[s].endi = extanch[k].endi; // map mae
				(*pairanch)[s].startj = extanch[k].startj; // map mae 
				(*pairanch)[s].endj = extanch[k].endj; // map mae
				(*pairanch)[s].score = extanch[k].score;
				s++;
			}
			if( extanch[k].j == m1[i] && extanch[k].i == m2[j] )
			{
				(*pairanch)[s].i = i;
				(*pairanch)[s].j = j;
				(*pairanch)[s].starti = extanch[k].startj; // map mae
				(*pairanch)[s].endi = extanch[k].endj; // map mae
				(*pairanch)[s].startj = extanch[k].starti; // map mae 
				(*pairanch)[s].endj = extanch[k].endi; // map mae
				(*pairanch)[s].score = extanch[k].score;
				s++;
			}
			k++;
		}
#endif
	}
	(*pairanch)[s].i = (*pairanch)[s].j = -1;

	recountpositions( *pairanch, n1, n2, seq1, seq2 );
//	truncateseq_group( *pairanch, seq1, seq2, n1, n2 );
//	copybackanchors( *pairanch, ddn1, n2, seq1, seq2 ); // tabun dame

#if 0
	reporterr( "Before check\n" );
	for( k=0; (*pairanch)[k].i>-1; k++ )
	{
		if( (*pairanch)[k].starti!=-1)
			reporterr( "seq1-%d,seq2-%d %d-%d,%d-%d\n", (*pairanch)[k].i, (*pairanch)[k].j, (*pairanch)[k].starti, (*pairanch)[k].endi, (*pairanch)[k].startj, (*pairanch)[k].endj );
	}
#endif

#if 0
	reporterr( "\ngroup1=\n" );
	for( i=0; m1[i]>-1; i++ )
		reporterr( "%d ", m1[i] );
	reporterr( "\n" );
	reporterr( "\ngroup2=\n" );
	for( i=0; m2[i]>-1; i++ )
		reporterr( "%d ", m2[i] );
	reporterr( "\n" );
#endif

	reporterr( "Checking external anchors\n" );
	checkanchors_strongestfirst( *pairanch, s );


#if 0
	reporterr( "After check\n" );
	for( k=0; (*pairanch)[k].i>-1; k++ )
	{
		if( (*pairanch)[k].starti!=-1)
			reporterr( "seq1-%d,seq2-%d %d-%d,%d-%d\n", (*pairanch)[k].i, (*pairanch)[k].j, (*pairanch)[k].starti, (*pairanch)[k].endi, (*pairanch)[k].startj, (*pairanch)[k].endj );
	}
#endif
}

#ifdef enablemultithread
/* work/prog: ready-queue scheduling of the guide-tree merges.
   The original loop hands out merge steps strictly in index order and a thread that takes step l
   blocks until both children of l are done, so independent subtrees overlap only within a window
   of nthreadtb consecutive steps.  Here a step is queued as soon as both children are done, and the
   ready step with the largest estimated remaining critical path (sum over the step and its
   ancestors of size^0.75, size = number of sequences under the node) runs first.  The alignment
   computed at each step depends only on its two children's results, so the output is unchanged. */
typedef struct _tbsched
{
	int *heap;     // ready steps, max-heap on prio (ties: smaller index first)
	int nheap;
	int ntaken;    // steps handed out so far
	int ndone;     // steps finished (work/prog2: idle threads stay until all are done, to help)
	int *parent;   // parent step of each step, -1 for the root
	double *prio;
	int inorder;   // 1: original in-order dispatch
	int optimistic; // work/prog2: ready queue while the memsave switch has not fired; abort on it (see below)
	int abort;     // work/prog2: set when a step finds a child longer than 30000 columns
	char *started; // work/prog2: step handed out (for the clean-up after an abort)
	int limit;     // work/largen: only steps with a smaller index run in this attempt
	int detect;    // work/largen: a step that finds a child longer than 30000 columns lowers limit to
	               // its own index and is put back, instead of abandoning the pass (see disttbfast())
	double heavyfactor;   // work/largen: > 0: heavy steps (below) may hold at most heavyfactor * njob * alloclen
	double heavyinflight; // bytes of traceback in heavy steps now running
	int *deferred;        // heavy steps put back for lack of budget, pushed again when a heavy step ends
	int ndeferred;
	int ndeferrals;       // (count, for the test builds' report)
} tbsched_t;

static int tbsched_before( tbsched_t *sc, int a, int b )
{
	if( sc->prio[a] != sc->prio[b] ) return( sc->prio[a] > sc->prio[b] );
	return( a < b );
}

static void tbsched_readyprio( tbsched_t *sc ) // work/v8s: the best ready step's priority for critical helping (dphelp.h)
{
	if( dphelppool ) dphelppool->readyprio = sc->nheap ? sc->prio[sc->heap[0]] : -1.0;
}

static void tbsched_push( tbsched_t *sc, int l )
{
	int k = sc->nheap++, up;
	while( k > 0 )
	{
		up = ( k - 1 ) / 2;
		if( !tbsched_before( sc, l, sc->heap[up] ) ) break;
		sc->heap[k] = sc->heap[up];
		k = up;
	}
	sc->heap[k] = l;
	tbsched_readyprio( sc );
}

static int tbsched_pop( tbsched_t *sc )
{
	int top = sc->heap[0], last = sc->heap[--sc->nheap], k = 0, c;
	while( ( c = 2 * k + 1 ) < sc->nheap )
	{
		if( c + 1 < sc->nheap && tbsched_before( sc, sc->heap[c+1], sc->heap[c] ) ) c++;
		if( !tbsched_before( sc, sc->heap[c], last ) ) break;
		sc->heap[k] = sc->heap[c];
		k = c;
	}
	if( sc->nheap ) sc->heap[k] = last;
	tbsched_readyprio( sc );
	return( top );
}

static int tbsched_isready( Treedep *dep, int l )
{
	return( ( dep[l].child0 == -1 || dep[dep[l].child0].done ) && ( dep[l].child1 == -1 || dep[dep[l].child1].done ) );
}

static void tbsched_init( tbsched_t *sc, int njob, Treedep *dep, int inorder, int optimistic, int limit, int ***topol, int *seqlen )
{
	int l, nstep = njob - 1;
	double *size = calloc( nstep, sizeof( double ) );
	double *lmax = calloc( nstep, sizeof( double ) ); // work/v8s: longest member sequence of each step's group
	sc->heap = calloc( nstep, sizeof( int ) );
	sc->parent = calloc( nstep, sizeof( int ) );
	sc->prio = calloc( nstep, sizeof( double ) );
	sc->started = calloc( nstep, sizeof( char ) );
	sc->deferred = calloc( nstep, sizeof( int ) );
	sc->ndeferred = 0;
	sc->ndeferrals = 0;
	sc->heavyfactor = 0.0;
	sc->heavyinflight = 0.0;
	sc->nheap = 0;
	sc->ntaken = 0;
	sc->ndone = 0;
	sc->inorder = inorder;
	sc->optimistic = optimistic && !inorder;
	sc->abort = 0;
	sc->limit = limit;  // work/largen: njob-1 = every step
	sc->detect = 0;
#ifdef PROGTIMING
	if( getenv( "MAFFT_INORDER" ) ) sc->inorder = 1;
#endif
	for( l=0; l<nstep; l++ ) sc->parent[l] = -1;
	for( l=0; l<nstep; l++ ) // children always have smaller indices
	{
		size[l] = ( dep[l].child0 == -1 ? 1.0 : size[dep[l].child0] ) + ( dep[l].child1 == -1 ? 1.0 : size[dep[l].child1] );
		if( dep[l].child0 != -1 ) sc->parent[dep[l].child0] = l;
		if( dep[l].child1 != -1 ) sc->parent[dep[l].child1] = l;
	}
	/* work/v8s: a step's cost estimate is now size x Lmax + 3 x L0 x L1 (characters its passes read, plus
	   the DP cells; Lmax the longest member sequence of the merged group, L0 / L1 those of the two
	   children), summed to the root as before (prog2: size^0.75).  The old estimate ignored the
	   lengths: on 16S.T (lengths 506-4,066) it ranked the real critical path ~1,000th of ~2,000 ready
	   steps, which started 0.2-0.3 s late in each pass at 32 threads (a replay of the measured step
	   times, tools/priosim.py: 1.83 -> 1.51 s and 1.51 -> 1.29 s; rnasim10k unchanged).  The order of
	   the ready queue never changes a result (see treebasethread).  MAFFT_SCHEDPRIO=0 (timing builds)
	   restores the old estimate. */
	{
		int old = 0;
#ifdef PROGTIMING
		old = ( getenv( "MAFFT_SCHEDPRIO" ) && atoi( getenv( "MAFFT_SCHEDPRIO" ) ) == 0 );
#endif
		for( l=0; l<nstep; l++ ) // children first
		{
			double l0 = ( dep[l].child0 == -1 ) ? ( seqlen ? (double)seqlen[topol[l][0][0]] : 1.0 ) : lmax[dep[l].child0];
			double l1 = ( dep[l].child1 == -1 ) ? ( seqlen ? (double)seqlen[topol[l][1][0]] : 1.0 ) : lmax[dep[l].child1];
			lmax[l] = MAX( l0, l1 );
			sc->prio[l] = old ? pow( size[l], 0.75 ) : size[l] * lmax[l] + 3.0 * l0 * l1; // (the step's own cost, for now)
		}
		for( l=nstep-1; l>=0; l-- ) // parents always have larger indices
			sc->prio[l] += ( sc->parent[l] == -1 ? 0.0 : sc->prio[sc->parent[l]] );
	}
	for( l=0; l<nstep; l++ ) // work/largen: steps done in an earlier phase of the pass are skipped
	{
		if( dep[l].done ) sc->ndone++;
		else if( l < sc->limit && tbsched_isready( dep, l ) ) tbsched_push( sc, l );
	}
	free( size );
	free( lmax );
}

static void tbsched_free( tbsched_t *sc )
{
	free( sc->started );
	free( sc->deferred );
	free( sc->heap );
	free( sc->parent );
	free( sc->prio );
}

static void tbsched_done( tbsched_t *sc, Treedep *dep, int l ) // call with the mutex held, after dep[l].done = 1
{
	int p = sc->parent[l];
	sc->ndone++;
	if( p != -1 && p < sc->limit && tbsched_isready( dep, p ) ) tbsched_push( sc, p );
}

/* work/prog3: deferred gap insertion (see lazygap.h).
   op[c]    gap table that the parent of step c returned instead of inserting it into the rows of
            group c (to be applied after everything pending inside c's subtree), or NULL
   stale[c] some table is pending inside c's subtree (the rows of group c are not current)
   len[c]   logical length of the rows of group c (what strlen would give if they were current)
   Rows of single sequences are always current: their side's gaps are inserted directly. */
typedef struct { char *gt; int lin, lout; } lzop_t;
typedef struct _lzstate
{
	lzop_t **op;
	char *stale;
	int *len;
	int nstep;
	char **shadow; // MAFFT_LAZYCHECK (timing builds): every row with all of its gaps applied at once
	int *shadowcap;
} lzstate_t;

static lzstate_t *lz_new( int njob, char **aseq, int check )
{
	lzstate_t *z = calloc( 1, sizeof( lzstate_t ) );
	int i;
	z->nstep = njob - 1;
	z->op = calloc( njob, sizeof( lzop_t * ) );
	z->stale = calloc( njob, sizeof( char ) );
	z->len = calloc( njob, sizeof( int ) );
	if( check )
	{
		z->shadow = calloc( njob, sizeof( char * ) );
		z->shadowcap = calloc( njob, sizeof( int ) );
		for( i=0; i<njob; i++ ) { z->shadowcap[i] = strlen( aseq[i] ) + 1; z->shadow[i] = malloc( z->shadowcap[i] ); strcpy( z->shadow[i], aseq[i] ); }
	}
	return( z );
}

static void lz_reset( lzstate_t *z, int njob, char **aseq ) // after an abandoned attempt (rows restored)
{
	int i;
	for( i=0; i<z->nstep; i++ )
	{
		if( z->op[i] ) { free( z->op[i]->gt ); free( z->op[i] ); z->op[i] = NULL; }
		z->stale[i] = 0; z->len[i] = 0;
	}
	if( z->shadow )
		for( i=0; i<njob; i++ )
		{
			int n = strlen( aseq[i] ) + 1;
			if( n > z->shadowcap[i] ) { z->shadow[i] = realloc( z->shadow[i], n ); z->shadowcap[i] = n; }
			strcpy( z->shadow[i], aseq[i] );
		}
}

static void lz_free( lzstate_t *z, int njob )
{
	int i;
	if( z == NULL ) return;
	for( i=0; i<z->nstep; i++ ) if( z->op[i] ) { free( z->op[i]->gt ); free( z->op[i] ); }
	free( z->op ); free( z->stale ); free( z->len );
	if( z->shadow ) { for( i=0; i<njob; i++ ) free( z->shadow[i] ); free( z->shadow ); free( z->shadowcap ); }
	free( z );
}

/* insert the gaps of gt (lnew entries, 'o' = next residue, '-' = gap) into s, which holds lori chars,
   in place from the end (the same as gapireru_inplace in Salignmm.c) */
static void lz_gapins( char *s, int lori, const char *gt, int lnew )
{
	const char gapchar = '-';
	int p, q = lori - 1;
	s[lnew] = 0;
	for( p=lnew-1; p>=0; p-- )
		s[p] = ( gt[p] == '-' ) ? gapchar : s[q--];
}

/* the table that applies op first and then comp (comp: lcomp entries over op->lout chars) */
static char *lz_compose( const char *comp, int lcomp, const lzop_t *op )
{
	char *out = malloc( lcomp + 1 );
	int p, q = 0;
	for( p=0; p<lcomp; p++ )
		out[p] = ( comp[p] == '-' ) ? '-' : op->gt[q++];
	out[lcomp] = 0;
	if( q != op->lout ) { reporterr( "\nlz_compose: %d != %d\n", q, op->lout ); exit( 1 ); }
	return( out );
}

typedef struct { int row; const char *gt; int lin, lout; char *nw, *old; } lzapp_t; // work/v8s: nw, old (lz_materialize)
/* work/v8s: the in-place rows of a step are grown (realloc to localcopylen + 11 bytes, prog3's rule) by
   the step's thread and idle helpers together.  Near the root, after alloclen has grown, that is a
   realloc (copy) of every row of the alignment, ~0.1 GB on 16S.T, on the critical chain.  Each row
   belongs to this step alone; realloc keeps its bytes, whichever thread calls it. */
typedef struct { char **aseq; int *rowcap; const int *mem; int need; } rowgrowctx_t;
static void rowgrow_rows( void *vctx, int lo, int hi )
{
	rowgrowctx_t *c = (rowgrowctx_t *)vctx;
	int i, j;
	for( i=lo; i<hi; i++ )
	{
		j = c->mem[i];
		if( c->rowcap[j] < c->need )
		{
			c->aseq[j] = realloc( c->aseq[j], c->need );
			if( c->aseq[j] == NULL ) ErrorExit( "Cannot reallocate a row in treebasethread.\n" );
			c->rowcap[j] = c->need;
		}
	}
}
static void rowgrow_all( char **aseq, int *rowcap, const int *mem, int need )
{
	rowgrowctx_t c;
	int n, ngrow = 0;
	for( n=0; mem[n]!=-1; n++ ) if( rowcap[mem[n]] < need ) ngrow++;
	if( ngrow == 0 ) return;
	c.aseq = aseq; c.rowcap = rowcap; c.mem = mem; c.need = need;
#ifdef enablemultithread
	dphelp_for( n, MAX( 16, n / 64 ), (double)ngrow * need, 4.0e6, rowgrow_rows, &c );
#else
	rowgrow_rows( &c, 0, n );
#endif
}

typedef struct { char **aseq; int *rowcap; int need; lzapp_t *app; } lzappctx_t;
static void lz_apply_rows( void *vctx, int lo, int hi )
{
	lzappctx_t *c = (lzappctx_t *)vctx;
	int k;
	for( k=lo; k<hi; k++ )
	{
		lzapp_t *a = c->app + k;
#ifdef PROGTIMING
		if( (int)strlen( c->aseq[a->row] ) != a->lin ) { reporterr( "\nlazy gaps: row %d has %d characters, its table expects %d\n", a->row, (int)strlen( c->aseq[a->row] ), a->lin ); exit( 1 ); }
#endif
		if( c->rowcap[a->row] < c->need )
		{
			/* work/v8s: a row that must grow gets its gaps on the way into the new buffer (one pass:
			   the same characters in the same places as realloc + lz_gapins, forward instead of
			   backward), instead of realloc copying it first.  The new buffers are allocated and the
			   old ones freed by the poster, outside the parallel part (lz_materialize): on glibc,
			   threads growing their arenas at once serialize on the address-space lock, which the
			   page faults of the parallel part also need */
			const char *old = c->aseq[a->row], *gt = a->gt;
			char *nw = a->nw ? a->nw : malloc( c->need );
			int p, q = 0;
			if( nw == NULL ) ErrorExit( "Cannot reallocate a row (lazy gaps).\n" );
			for( p=0; p<a->lout; p++ ) nw[p] = ( gt[p] == '-' ) ? '-' : old[q++];
			nw[a->lout] = 0;
			if( a->nw ) a->old = (char *)old; else free( (char *)old );
			c->aseq[a->row] = nw;
			c->rowcap[a->row] = c->need;
		}
		else
			lz_gapins( c->aseq[a->row], a->lin, a->gt, a->lout );
	}
}

static int lz_twophase( void ) // work/v8s: MAFFT_LZTWOPHASE=0 (timing builds) allocates in the parallel part (A/B)
{
#ifdef PROGTIMING
	static int v = -1;
	if( v < 0 ) { char *e = getenv( "MAFFT_LZTWOPHASE" ); v = !( e && atoi( e ) == 0 ); }
	return( v );
#else
	return( 1 );
#endif
}

/* bring the rows of group x (step x) up to date: every table pending inside x's subtree, composed
   per row (deepest first) and inserted once.  need: row size to guarantee (> the new length). */
static void lz_materialize( lzstate_t *z, int x, Treedep *dep, int ***topol, char **aseq, int *rowcap, int need )
{
	typedef struct { int node; char *comp; int lcomp, lin; } fr_t; // lin: the input length of comp
	fr_t *stk;
	int nstk = 0, capstk = 64, napp = 0, capapp = 64, nown = 0, capown = 64, k, s;
	lzapp_t *app;
	char **own; // composed tables to free at the end
	int *usedop, nused = 0, capused = 64;
	if( !z->stale[x] ) return;
	stk = malloc( sizeof( fr_t ) * capstk );
	app = malloc( sizeof( lzapp_t ) * capapp );
	own = malloc( sizeof( char * ) * capown );
	usedop = malloc( sizeof( int ) * capused );
	stk[nstk].node = x; stk[nstk].comp = NULL; stk[nstk].lcomp = 0; stk[nstk].lin = 0; nstk++;
	while( nstk )
	{
		fr_t f = stk[--nstk];
		z->stale[f.node] = 0;
		for( s=0; s<2; s++ )
		{
			int c = s ? dep[f.node].child1 : dep[f.node].child0;
			char *comp = f.comp;
			int lcomp = f.lcomp, lin = f.lin;
			if( c != -1 && z->op[c] )
			{
				lin = z->op[c]->lin; // a composite's input is the input of its innermost table
				if( comp == NULL ) { comp = z->op[c]->gt; lcomp = z->op[c]->lout; }
				else
				{
					comp = lz_compose( comp, lcomp, z->op[c] );
					if( nown == capown ) own = realloc( own, sizeof( char * ) * ( capown *= 2 ) );
					own[nown++] = comp;
				}
				if( nused == capused ) usedop = realloc( usedop, sizeof( int ) * ( capused *= 2 ) );
				usedop[nused++] = c;
			}
			if( c == -1 )
			{
				if( comp )
				{
					if( napp == capapp ) app = realloc( app, sizeof( lzapp_t ) * ( capapp *= 2 ) );
					app[napp].row = topol[f.node][s][0];
					app[napp].gt = comp;
					app[napp].lout = lcomp;
					app[napp].lin = lin;
					app[napp].nw = app[napp].old = NULL;
					napp++;
				}
			}
			else if( comp || z->stale[c] )
			{
				if( nstk == capstk ) stk = realloc( stk, sizeof( fr_t ) * ( capstk *= 2 ) );
				stk[nstk].node = c; stk[nstk].comp = comp; stk[nstk].lcomp = lcomp; stk[nstk].lin = lin; nstk++;
			}
		}
	}
	{
		lzappctx_t ctx;
		int ngrow = 0;
		ctx.aseq = aseq; ctx.rowcap = rowcap; ctx.need = need; ctx.app = app;
		if( lz_twophase() ) // work/v8s: the new buffers of the rows that grow, by this thread (see lz_apply_rows)
			for( k=0; k<napp; k++ )
				if( rowcap[app[k].row] < need )
				{
					if( ( app[k].nw = malloc( need ) ) == NULL ) ErrorExit( "Cannot reallocate a row (lazy gaps).\n" );
					ngrow++;
				}
#ifdef enablemultithread
		dphelp_for( napp, MAX( 16, napp / 64 ), (double)napp * need, 4.0e6, lz_apply_rows, &ctx );
#else
		lz_apply_rows( &ctx, 0, napp );
#endif
		if( ngrow ) for( k=0; k<napp; k++ ) free( app[k].old ); // (NULL for the rows that did not grow)
	}
	for( k=0; k<nown; k++ ) free( own[k] );
	for( k=0; k<nused; k++ ) { free( z->op[usedop[k]]->gt ); free( z->op[usedop[k]] ); z->op[usedop[k]] = NULL; }
	free( stk ); free( app ); free( own ); free( usedop );
}

#ifdef PROGTIMING
/* MAFFT_LAZYCHECK=1: shadow rows get every gap table the moment it is made; the real rows must equal
   them whenever they are current (checked before a step uses current rows and after every
   materialization). */
static void lz_check_rows( lzstate_t *z, char **aseq, int *mem, const char *what, int l )
{
	int i, j;
	for( i=0; (j=mem[i])!=-1; i++ )
		if( strcmp( z->shadow[j], aseq[j] ) )
		{
			reporterr( "\nLAZYCHECK MISMATCH (%s) at step %d, row %d\n", what, l, j );
			exit( 1 );
		}
}
static void lz_shadow_set( lzstate_t *z, char **aseq, int *mem ) // the rows are current: copy them
{
	int i, j;
	for( i=0; (j=mem[i])!=-1; i++ )
	{
		int n = strlen( aseq[j] ) + 1;
		if( n > z->shadowcap[j] ) { z->shadow[j] = realloc( z->shadow[j], n ); z->shadowcap[j] = n; }
		strcpy( z->shadow[j], aseq[j] );
	}
}
static void lz_shadow_ins( lzstate_t *z, int *mem, const char *gt, int lin, int lout ) // a deferred table
{
	int i, j;
	for( i=0; (j=mem[i])!=-1; i++ )
	{
		if( (int)strlen( z->shadow[j] ) != lin ) { reporterr( "\nLAZYCHECK: shadow row %d has length %d, table expects %d\n", j, (int)strlen( z->shadow[j] ), lin ); exit( 1 ); }
		if( lout + 1 > z->shadowcap[j] ) { z->shadow[j] = realloc( z->shadow[j], lout + 1 ); z->shadowcap[j] = lout + 1; }
		lz_gapins( z->shadow[j], lin, gt, lout );
	}
}
static long lz_nmat = 0, lz_ndefer = 0;
static void lz_report( void ) __attribute__(( destructor ));
static void lz_report( void ) { if( getenv( "MAFFT_LAZYCHECK" ) ) reporterr( "\nLAZYCHECK: %ld tables deferred, %ld materializations, no mismatch\n", lz_ndefer, lz_nmat ); }
#endif

/* work/prog2: copying the member rows in and out of the local copies, split among idle helpers
   for large groups (each row is copied by exactly one thread; strcpy as before). */
/* work/v7: with prog3's / largen's in-place rows (targ->nocopy, always on in normal builds) none of
   the copy code below, aseqlock or the ReallocateCharMtx of all rows runs; they remain for the
   MAFFT_COPYROWS=1 A/B switch of the timing builds (-DPROGTIMING) only, with gccpar's deadlock fix. */
typedef struct { char **dst, **src; int *mem; } rowcopy_t;
static void rowcopy_body( void *vctx, int lo, int hi )
{
	rowcopy_t *c = (rowcopy_t *)vctx;
	int i;
	for( i=lo; i<hi; i++ ) strcpy( c->dst[c->mem[i]], c->src[c->mem[i]] );
}
/* work/gccpar: DEADLOCK FIX.  rowcopy_all() is called with targ->aseqlock held for reading, and
   dphelp_for() takes the pool mutex, which is the pass's global mutex (dphelp_post/dphelp_finish).
   A thread that reallocates aseq holds that global mutex while it waits in
   pthread_rwlock_wrlock( aseqlock ) for the readers to leave, so a reader that reaches
   dphelp_post/dphelp_finish waits for the mutex forever, and so does everyone else (observed:
   rnasim10k, --thread 10, gcc build, all threads blocked, stacks in
   treebasethread -> dphelp_for -> dphelp_finish and treebasethread -> pthread_rwlock_wrlock).
   The copies are now done by the calling thread alone (strcpy of the same rows, same result);
   helpers keep doing the DP tiles, profiles and gap insertion, which run outside the lock. */
static void rowcopy_all( char **dst, char **src, int *mem, int n, int rowlen )
{
	rowcopy_t c;
	c.dst = dst; c.src = src; c.mem = mem;
	(void)rowlen;
	if( n > 0 ) rowcopy_body( &c, 0, n );
}

#ifdef PROGTIMING
#define PT_LZMAT( mem ) do { if( targ->lz->shadow ) { lz_check_rows( targ->lz, aseq, (mem), "materialized", l ); __atomic_add_fetch( &lz_nmat, 1, __ATOMIC_RELAXED ); } } while( 0 )
#else
#define PT_LZMAT( mem ) do { } while( 0 )
#endif

static void *treebasethread( void *arg )
{
	treebasethread_arg_t *targ = (treebasethread_arg_t *)arg;
	int thread_no = targ->thread_no;
	int *nrunpt = targ->nrunpt;
	int njob = targ->njob;
	int *nlen = targ->nlen;
	int *jobpospt = targ->jobpospt;
	int ***topol = targ->topol;
	Treedep *dep = targ->dep;
	double ***cpmxhist = targ->cpmxhist;
	int **memhist = targ->memhist;
	char **aseq = targ->aseq;
	double *effarr = targ->effarr;
	int *alloclen = targ->alloclenpt;
	int *fftlog = targ->fftlog;
	char *mergeoralign = targ->mergeoralign;
	double **newdistmtx = targ->newdistmtx;
	int *selfscore = targ->selfscore;
	ExtAnch *extanch = targ->extanch;
	int **anchindex = targ->anchindex;

	char **mseq1, **mseq2;
	char **localcopy;
	int i, m, j, l;
	int immin, immax;
	int len1, len2;
	int clus1, clus2;
	double pscore, tscore;
	char *indication1, *indication2;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
//	double dumfl = 0.0;
	double dumdb = 0.0;
	int ffttry;
	int m1, m2;
	double **dynamicmtx;
	int ssi, ssm, bunbo;
	int tm, ti;
	int **localmem = NULL;
	double ***cpmxchild0, ***cpmxchild1;
	double orieff1, orieff2;
	ExtAnch *pairanch = NULL;
	int localcopylen = 0;
	double heavybytes = 0.0; // work/largen: traceback bytes this step counts against the heavy-step budget
	int lz_on = 0, lz_stale0 = 0, lz_stale1 = 0, lz_defer0 = 0, lz_defer1 = 0; // work/prog3
	int gskip[2] = { 0, 0 }; // work/v8s: side left untouched by this step (deferred, not materialized): its rows are not grown
	char *lcblk = NULL; // work/prog2: per-thread block holding the local copies of the member rows
	size_t lccap = 0;
#if SKIP
	int **skiptable1 = NULL, **skiptable2 = NULL;
#endif
#if 0
	int i, j;
#endif



	tscore = 0;
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );
	localcopy = calloc( njob, sizeof( char * ) );
	for( i=0; i<njob; i++ ) localcopy[i] = NULL;
	effarr1 = AllocateDoubleVec( njob );
	effarr2 = AllocateDoubleVec( njob );
	indication1 = AllocateCharVec( 150 );
	indication2 = AllocateCharVec( 150 );
	if( specificityconsideration )
		dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );
	localmem = (int **)calloc( sizeof( int * ), 2 );



#if 0
	reporterr(       "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		reporterr(       "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

//	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );


//	writePre( njob, name, nlen, aseq, 0 );

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		if( targ->sched && !targ->sched->inorder )
		{
			tbsched_t *sc = targ->sched;
			/* work/prog2: a thread with no ready step helps with posted jobs (large alignments of
			   other threads, see dphelp.h) and stays until every step is finished, so that the
			   last merges at the root still find helpers.  Without a pool it leaves as soon as
			   every step has been handed out, as before. */
			while( 1 )
			{
				if( sc->abort ) break; // work/prog2: the pass is being restarted in order (see below)
				while( sc->nheap > 0 && sc->heap[0] >= sc->limit ) tbsched_pop( sc ); // work/largen: left for a later phase
				if( sc->nheap > 0 )
				{
					/* work/v8s: a posted job of a step with a higher priority than every ready step (the
					   critical path) first (dphelp.h) */
					if( dphelppool && dphelppool->jobs && dphelp_crit_on() && dphelp_serve_crit_locked( sc->prio[sc->heap[0]] ) ) continue;
					break;
				}
				if( dphelppool && dphelp_serve_locked() ) continue;
				if( dphelppool ? ( sc->ndone == njob-1 ) : ( sc->ntaken == njob-1 ) ) break;
				if( *nrunpt == 0 ) break; // work/largen: nothing ready, nothing running: the rest belongs to a later phase
				if( dphelppool ) dphelppool->nidle++;
				pthread_cond_wait( targ->treecond, targ->mutex );
				if( dphelppool ) dphelppool->nidle--;
			}
			if( sc->nheap == 0 || sc->abort )
			{
				l = njob-1; // nothing left: exit below
				pthread_cond_broadcast( targ->treecond ); // wake idle threads so they can exit
			}
			else
			{
				l = tbsched_pop( sc );
				sc->ntaken++;
				sc->started[l] = 1;
				dphelp_myprio = sc->prio[l]; // work/v8s: this thread's step (critical helping, dphelp.h)
			}
		}
		else
		{
			/* work/v10w: the original in-order dispatch.  A thread that finds every step handed out
			   used to leave; with a helper pool it now stays and helps the steps still running (the
			   last merges at the root, which run one at a time) until the root is done.  Threads beyond
			   the original's count (helperonly: --thread above the 16 that mafft.tmpl gives the pass)
			   never take a step; they only help.  Which thread runs which step, and in which order steps
			   are handed out and check the memsave switch, is the original's. */
			l = targ->helperonly ? njob-1 : *jobpospt;
			if( l == njob-1 && dphelppool && njob > 1 )
			{
				/* (memory: the step code's per-thread buffers are given back first, as at the exit below) */
				pthread_mutex_unlock( targ->mutex );
				if( commonIP ) FreeIntMtx( commonIP );
				commonIP = NULL;
				Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
				Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
				A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
				G__align11( NULL, NULL, NULL, 0, 0, 0 );
				pthread_mutex_lock( targ->mutex );
				while( !dep[njob-2].done )
				{
					if( dphelp_serve_locked() ) continue;
					dphelppool->nidle++;
					pthread_cond_wait( targ->treecond, targ->mutex );
					dphelppool->nidle--;
				}
			}
		}
		if( l == njob-1 )
		{
			pthread_mutex_unlock( targ->mutex );
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
			Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
			Falign_givenanchors( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
			A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
			D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
			G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
			free( mseq1 );
			free( mseq2 );
			free( localcopy );
			free( effarr1 );
			free( effarr2 );
			free( indication1 );
			free( indication2 );
			if( specificityconsideration )
				FreeDoubleMtx( dynamicmtx );
			free( localmem );
			free( lcblk );
			return( NULL );
		}
		if( !targ->sched || targ->sched->inorder ) *jobpospt = l+1;
		PTT(0);

//		reporterr( "l=%d, child0=%d, child1=%d\n", l, dep[l].child0, dep[l].child1 );

		if( dep[l].child0 != -1 )
		{
			while( dep[dep[l].child0].done == 0 )
			{
				if( dphelppool && dphelp_serve_locked() ) continue; // work/prog2: help while waiting
				if( dphelppool ) dphelppool->nidle++;
				pthread_cond_wait( targ->treecond, targ->mutex );
				if( dphelppool ) dphelppool->nidle--;
			}
		}
		if( dep[l].child1 != -1 ) 
		{
			while( dep[dep[l].child1].done == 0 )
			{
				if( dphelppool && dphelp_serve_locked() ) continue;
				if( dphelppool ) dphelppool->nidle++;
				pthread_cond_wait( targ->treecond, targ->mutex );
				if( dphelppool ) dphelppool->nidle--;
			}
		}




//		while( *nrunpt >= nthread ) // bug
		while( *nrunpt >= targ->nrunmax ) // tabun iranai  (work/v7: the threads of this attempt, nthreadtb in the original)
			pthread_cond_wait( targ->treecond, targ->mutex ); // tabun iranai
		(*nrunpt)++;
		PTT(1);

		m1 = topol[l][0][0];
		m2 = topol[l][1][0];
#if 0
		localmem[0][0] = -1;
		posinmem=topolorderz( localmem[0], topol, dep, l, 0 ) - localmem[0];
		localmem[1][0] = -1;
		posinmem=topolorderz( localmem[1], topol, dep, l, 1 ) - localmem[1];
#else
		if( dep[l].child0 == -1 ) 
		{
			localmem[0] = calloc( sizeof( int ), 2 );
			localmem[0][0] = m1;
			localmem[0][1] = -1;
			clus1 = 1;
		}
		else
		{
			localmem[0] = memhist[dep[l].child0];
			clus1 = intlen( localmem[0] );
		}
		if( dep[l].child1 == -1 ) 
		{
			localmem[1] = calloc( sizeof( int ), 2 );
			localmem[1][0] = m2;
			localmem[1][1] = -1;
			clus2 = 1;
		}
		else
		{
			localmem[1] = memhist[dep[l].child1];
			clus2 = intlen( localmem[1] );
		}

		if( l != njob-2 )
		{
			memhist[l] = calloc( sizeof( int ), clus1+clus2+1 );
			intcpy( memhist[l], localmem[0] );
			intcpy( memhist[l]+clus1, localmem[1] );
			memhist[l][clus1+clus2] = -1;
		}
#endif

// moved, 2018/Mar/10. Must be after changing memhist[l]
		if( mergeoralign[l] == 'n' )
		{
//			reporterr(       "SKIP!\n" );
			dep[l].done = 1;
			if( targ->sched ) tbsched_done( targ->sched, dep, l );
			(*nrunpt)--;
			pthread_cond_broadcast( targ->treecond );
//			free( topol[l][0] ); topol[l][0] = NULL;
//			free( topol[l][1] ); topol[l][1] = NULL;
//			free( topol[l] ); topol[l] = NULL;
			pthread_mutex_unlock( targ->mutex );
			free( localmem[0] );
			free( localmem[1] );
			continue;
		}


//		reporterr( "l=%d, dep[l].child0=%d, dep[l].child1=%d\n", l, dep[l].child0, dep[l].child1 );
		if( dep[l].child0 == -1 ) cpmxchild0 = NULL; else cpmxchild0 = cpmxhist+dep[l].child0;
		if( dep[l].child1 == -1 ) cpmxchild1 = NULL; else cpmxchild1 = cpmxhist+dep[l].child1;
//		reporterr( "cpmxchild0=%p, cpmxchild1=%p\n", cpmxchild0, cpmxchild1 );


//		reporterr(       "\ndistfromtip = %f\n", dep[l].distfromtip );
		if( specificityconsideration )
			makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip );
		else
			dynamicmtx = n_dis_consweight_multi;
//		reporterr( "dynamicmtx[0][1] = %f\n", dynamicmtx[0][1] );

		len1 = ( targ->lz && dep[l].child0 != -1 ) ? targ->lz->len[dep[l].child0] : strlen( aseq[m1] ); // work/prog3: rows may be stale (lazygap.h)
		len2 = ( targ->lz && dep[l].child1 != -1 ) ? targ->lz->len[dep[l].child1] : strlen( aseq[m2] );
		if( *alloclen <= len1 + len2 )
		{
			reporterr(       "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			if( !targ->nocopy ) // work/prog3: with nocopy, each step grows only its own rows (below)
			{
			pthread_rwlock_wrlock( targ->aseqlock ); // wait until no thread is copying rows of aseq
			ReallocateCharMtx( aseq, njob, *alloclen + 10  ); 
			pthread_rwlock_unlock( targ->aseqlock );
			}
			reporterr(       "done. *alloclen = %d\n", *alloclen );
		}
		localcopylen = *alloclen; // work/prog: the member rows are copied after the mutex is released

		/* work/prog2: with the ready queue, steps run out of index order, and the one-way switch
		   below would then take effect at a different step than in the original in-order dispatch.
		   So the first step that would flip it abandons the pass instead; disttbfast() then runs the
		   pass again from the start with the original dispatch (bitwise the same as before as long
		   as the switch never fires, since until then every step uses the same aligner). */
		if( targ->sched && targ->sched->optimistic && targ->sched->detect && !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000  ) ) )
		{
			/* work/largen (more than 20000 sequences, see disttbfast()): the original, serial,
			   switches to the memsave algorithm at the first step in index order that finds a
			   child longer than 30000 columns.  Put this step back as if it had never been taken
			   and run no step with an index >= l in this attempt; the pass is finished in later
			   phases once the smallest such index is known.  (work/v7: nothing of the step has
			   been touched yet: no row grown, no lazy-gap table materialized or requested.) */
			tbsched_t *sc = targ->sched;
			if( l < sc->limit ) sc->limit = l;
			sc->started[l] = 0;
			sc->ntaken--;
			if( l != njob-2 ) { free( memhist[l] ); memhist[l] = NULL; }
			if( dep[l].child0 == -1 ) free( localmem[0] ); // the children's member lists stay in memhist
			if( dep[l].child1 == -1 ) free( localmem[1] );
			(*nrunpt)--;
			pthread_cond_broadcast( targ->treecond );
			pthread_mutex_unlock( targ->mutex );
			continue;
		}
		if( targ->sched && targ->sched->optimistic && !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000  ) ) )
		{
			targ->sched->abort = 1;
			(*nrunpt)--;
			pthread_cond_broadcast( targ->treecond );
			pthread_mutex_unlock( targ->mutex );
			free( localmem[0] ); // the children's member lists, consumed as a finished step would
			free( localmem[1] );
			continue;
		}
		if( !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000  ) ) )
		{
			reporterr(       "\nlen1=%d, len2=%d, Switching to the memsave mode\n", len1, len2 );
			alg = 'M';
		}

		/* work/largen: memory (above 20000 sequences, see disttbfast()).  A step that will run the
		   full A__align DP (no FFT, not 1-vs-1) with a traceback of at least tb_heavymin bytes is
		   heavy.  Heavy steps run together only while their tracebacks fit in heavyfactor * njob *
		   alloclen bytes (the original's A__align allocates a whole njob x alloclen copy of the
		   alignment, mseq, at the root step) and no heavy step is waiting; otherwise the step is put
		   back, as if never taken, until a heavy step ends.  A heavy step always runs when no other
		   one does.  Idle threads meanwhile help the running steps.  (work/v7: len1/len2 are the
		   logical lengths with deferred gaps; the put-back happens before anything of the step is
		   touched.) */
		heavybytes = 0.0;
		if( targ->sched && targ->sched->heavyfactor > 0.0 && !targ->sched->inorder && alg != 'M' && !useexternalanchors
			&& !( clus1 == 1 && clus2 == 1 ) && !( force_fft || ( use_fft && nlen[m1] > clus1 && nlen[m2] > clus2 ) ) )
		{
			tbsched_t *sc = targ->sched;
			double hb = (double)( len1 + 2 ) * (double)( len2 + 2 ) * sizeof( int );
			if( hb >= tb_heavymin )
			{
				if( sc->heavyinflight > 0.0 && ( sc->ndeferred > 0 || sc->heavyinflight + hb > sc->heavyfactor * (double)njob * (double)*alloclen ) ) // (no overtaking a deferred step)
				{
					sc->deferred[sc->ndeferred++] = l;
					sc->ndeferrals++;
					sc->started[l] = 0;
					sc->ntaken--;
					if( l != njob-2 ) { free( memhist[l] ); memhist[l] = NULL; }
					if( dep[l].child0 == -1 ) free( localmem[0] );
					if( dep[l].child1 == -1 ) free( localmem[1] );
					(*nrunpt)--;
					pthread_cond_broadcast( targ->treecond );
					pthread_mutex_unlock( targ->mutex );
					continue;
				}
				sc->heavyinflight += hb;
				heavybytes = hb;
			}
		}

		/* work/v10w: the step's aligner is taken from alg here, under the mutex, right after the switch
		   check, instead of after the unlock (the original reads it a few microseconds later, unlocked).
		   That is one of the original's own timings: whatever the original's step could read, it can
		   read this value (alg only goes from 'A' to 'M', under this mutex).  It fixes the race at the
		   start of an in-order pass, where the original's row copies inside the mutex let the steps
		   claimed before the switching one read 'A' nearly always, and v9's shorter critical section
		   (rows grown after the unlock) made them read 'M' as often as not (notes/v10w-lever.md).
		   Reads inside Falign (per segment) stay as in the original. */
		char algstep = alg;
		if( alg == 'M' ) // hoka no thread ga M ni shitakamo shirenainode
		{
//			reporterr(       "Freeing commonIP (thread %d)\n", thread_no );
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}

		pthread_mutex_unlock( targ->mutex );

		if( targ->nocopy )
		{
			/* work/prog3: no local copies.  The rows of the two children belong to this step alone
			   until it is marked done (the parent starts only after that, through the mutex), and no
			   step reallocates rows it does not own any more, so the aligners can work on the rows of
			   aseq directly.  Each row only has to be as large as the local copy was (localcopylen,
			   i.e. alloclen at the length check, > len1+len2); rows are grown here, one by one, to
			   localcopylen + 11 bytes (what ReallocateCharMtx( aseq, njob, alloclen + 10 ) gave every
			   row in the original). */
			int need = localcopylen + 11, s;
#ifdef PROGTIMING
			double _g0 = pt_now();
#endif
			/* work/v8s: a side this step will hand to A__align as deferred (lazy gaps) without
			   materializing it is neither read beyond its current length nor written (lazygap.h), so
			   its rows need no room for the new length now; they are grown when they are materialized
			   (lz_apply_rows grows each row as it inserts the gaps, with helpers).  Near the root that
			   moves the realloc of ~7,000 rows per alloclen increase (30-80 ms per step on CARC, on the
			   critical chain) into the one materialization at the root.  The prediction repeats the
			   decisions of the lazy-gap block below (checked there). */
			gskip[0] = gskip[1] = 0;
			if( targ->lz )
			{
				int ff = ( nlen[m1] > clus1 && nlen[m2] > clus2 ); // = ffttry below
				int useA = !useexternalanchors && !( force_fft || ( use_fft && ff ) ) && alg == 'A' && !( clus1 == 1 && clus2 == 1 ) && len1 > 0 && len2 > 0;
				if( useA )
				{
					int big = ( clus1 + clus2 > 20 && l != njob-2 && specificityconsideration == 0.0 );
					for( s=0; s<2; s++ )
					{
						int c = s ? dep[l].child1 : dep[l].child0;
						double ***cc = s ? cpmxchild1 : cpmxchild0;
						int st = ( c != -1 && targ->lz->stale[c] );
						if( st ) gskip[s] = ( cc && *cc ); // kept stale (else materialized: grown first, see below)
						else gskip[s] = ( c != -1 && big ); // current rows, gaps deferred
					}
				}
			}
			for( s=0; s<2; s++ ) // work/v8s: with idle helpers for large groups (rowgrow_all)
				if( !gskip[s] ) rowgrow_all( aseq, targ->rowcap, localmem[s], need );
#ifdef PROGTIMING
			pt_growt = pt_now() - _g0; // work/v8s: row growth time of this step (logged as acc[41])
#endif
			for( i=0; (j=localmem[0][i])!=-1; i++ ) localcopy[j] = aseq[j];
			for( i=0; (j=localmem[1][i])!=-1; i++ ) localcopy[j] = aseq[j];
		}
		else
		{
		/* work/prog: copy the member rows outside the global mutex.  The rows belong to finished
		   children, so nobody writes them now; the read lock only excludes a concurrent
		   ReallocateCharMtx() by another thread.  Buffer size as in the original (alloclen at the
		   time of the length check, > len1+len2). */
		/* work/prog2: the copies live in one per-thread block (rows of localcopylen chars) that is
		   reused across steps, instead of a calloc() per row per step. */
		{
			size_t need = (size_t)localcopylen * ( clus1 + clus2 ), k = 0;
			if( need > lccap )
			{
				free( lcblk );
				lccap = need + need / 2;
				lcblk = malloc( lccap );
				if( lcblk == NULL ) ErrorExit( "Cannot allocate the local copies in treebasethread.\n" );
			}
			for( i=0; (j=localmem[0][i])!=-1; i++ ) { localcopy[j] = lcblk + k; k += localcopylen; }
			for( i=0; (j=localmem[1][i])!=-1; i++ ) { localcopy[j] = lcblk + k; k += localcopylen; }
		}
		pthread_rwlock_rdlock( targ->aseqlock );
		if( clus1 + clus2 > 256 ) // work/prog2: large groups: row copies shared with idle helpers
		{
			rowcopy_all( localcopy, aseq, localmem[0], clus1, len1 );
			rowcopy_all( localcopy, aseq, localmem[1], clus2, len2 );
		}
		else
		{
		for( i=0; (j=localmem[0][i])!=-1; i++ ) strcpy( localcopy[j], aseq[j] );
		for( i=0; (j=localmem[1][i])!=-1; i++ ) strcpy( localcopy[j], aseq[j] );
		}
		pthread_rwlock_unlock( targ->aseqlock );
		}

#if 1 // CHUUI@@@@
#ifdef PROGTIMING
		PTT(2);
		if( ptlog ) { int k; for( k=0; k<PT_NACC; k++ ) pt_acc[k] = 0.0; }
#endif
		clus1 = fastconjuction_noname( localmem[0], localcopy, mseq1, effarr1, effarr, indication1, 0.0, &orieff1 );
		clus2 = fastconjuction_noname( localmem[1], localcopy, mseq2, effarr2, effarr, indication2, 0.0, &orieff2 );
#else
		clus1 = fastconjuction_noweight( topol[l][0], localcopy, mseq1, effarr1,  indication1 );
		clus2 = fastconjuction_noweight( topol[l][1], localcopy, mseq2, effarr2,  indication2 );
#endif


#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
#endif

#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
#endif


		PTT(3);
//		fprintf( trap_g, "\nSTEP-%d\n", l );
//		fprintf( trap_g, "group1 = %s\n", indication1 );
//		fprintf( trap_g, "group2 = %s\n", indication2 );

//		reporterr(       "\rSTEP % 5d / %d %d-%d", l+1, njob-1, clus1, clus2 );
		if( l < 500 || l % 100 == 0 ) reporterr(       "\rSTEP % 5d / %d (thread %4d)", l+1, njob-1, thread_no );
#if 0
		reporterr( "\nclus1=%d, clus2=%d\n", clus1, clus2 );
#endif

#if 0
		reporterr(       "STEP %d /%d\n", l+1, njob-1 );
		reporterr(       "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
		reporterr(       "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
#endif

/*
		reporterr(       "before align all\n" );
		display( aseq, njob );
		reporterr(       "\n" );
		reporterr(       "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		reporterr(       "\n" );
		reporterr(       "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		reporterr(       "\n" );
*/



//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		else						   ffttry = 0;
		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000); // v6.708
//		reporterr(       "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (double)len1/fftlog[m1], clus1, (double)len2/fftlog[m2], clus2 );
//		reporterr( "fftlog=%d,%d, ffttry=%d\n", fftlog[m1], fftlog[m2], ffttry );

		/* work/prog3: deferred gap insertion (lazygap.h).  Rows a step reads must be current; a group
		   with a cached profile that goes into A__align may stay stale, and its new gaps (and those of
		   the other group, when the result will have a cached profile too) can be deferred again. */
		lz_on = lz_stale0 = lz_stale1 = 0;
		if( targ->lz )
		{
			lzstate_t *z = targ->lz;
			int c0 = dep[l].child0, c1 = dep[l].child1;
			int useA = !useexternalanchors && !( force_fft || ( use_fft && ffttry ) ) && alg == 'A' && !( clus1 == 1 && clus2 == 1 ) && len1 > 0 && len2 > 0;
			lz_stale0 = ( c0 != -1 && z->stale[c0] );
			lz_stale1 = ( c1 != -1 && z->stale[c1] );
			if( ( gskip[0] || gskip[1] ) && !useA ) ErrorExit( "work/v8s: row growth skipped for a step that does not use A__align.\n" ); // (cannot happen: see the prediction)
			if( lz_stale0 && !( useA && cpmxchild0 && *cpmxchild0 ) )
			{
				if( gskip[0] ) ErrorExit( "work/v8s: row growth skipped for a side that is materialized.\n" );
				lz_materialize( z, c0, dep, topol, aseq, targ->rowcap, localcopylen + 11 );
				lz_stale0 = 0;
				PT_LZMAT( localmem[0] );
			}
			if( lz_stale1 && !( useA && cpmxchild1 && *cpmxchild1 ) )
			{
				if( gskip[1] ) ErrorExit( "work/v8s: row growth skipped for a side that is materialized.\n" );
				lz_materialize( z, c1, dep, topol, aseq, targ->rowcap, localcopylen + 11 );
				lz_stale1 = 0;
				PT_LZMAT( localmem[1] );
			}
#ifdef PROGTIMING
			if( z->shadow ) // the current rows must equal the shadow rows
			{
				if( !lz_stale0 ) lz_check_rows( z, aseq, localmem[0], "current rows, side 0", l );
				if( !lz_stale1 ) lz_check_rows( z, aseq, localmem[1], "current rows, side 1", l );
			}
#endif
			if( useA )
			{
				int big = ( clus1 + clus2 > 20 && l != njob-2 && specificityconsideration == 0.0 );
				lz_defer0 = lz_stale0 || ( c0 != -1 && big );
				lz_defer1 = lz_stale1 || ( c1 != -1 && big );
				if( ( gskip[0] && !lz_defer0 ) || ( gskip[1] && !lz_defer1 ) ) ErrorExit( "work/v8s: row growth skipped for a side that is not deferred.\n" ); // (cannot happen)
				if( lz_defer0 || lz_defer1 )
				{
					A_lz.on = 1; A_lz.used = 0;
					A_lz.stale[0] = lz_stale0; A_lz.stale[1] = lz_stale1;
					A_lz.len[0] = len1; A_lz.len[1] = len2;
					A_lz.defer[0] = lz_defer0; A_lz.defer[1] = lz_defer1;
					lz_on = 1;
				}
			}
		}

		if( useexternalanchors )
		{
//			reporterr( "%%%% %d vs %d\n", m1, m2 );
			pickpairanch( &pairanch, extanch, anchindex, clus1, clus2, localmem[0], localmem[1], mseq1, mseq2 );
//			reporterr( "pairanch: %d:%d\n", pairanch[0].starti, pairanch[0].startj );
			pscore = Falign_givenanchors( pairanch, NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			free( pairanch );
			pairanch = NULL;
		}
		else if( force_fft || ( use_fft && ffttry ) )
		{
			if( l < 500 || l % 100 == 0 ) reporterr(       " f\b\b" );
			if( algstep == 'M' ) // work/v10w: see algstep
			{
				if( l < 500 || l % 100 == 0 ) reporterr(       "m" );
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
			{
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
			}
		}
		else
		{
			if( l < 500 || l % 100 == 0 ) reporterr(       " d\b\b" );
			fftlog[m1] = 0;
			switch( algstep ) // work/v10w: see algstep
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					if( l < 500 || l % 100 == 0 ) reporterr(       "m" );
					if( l < 500 || l % 100 == 0 ) if( ( cpmxchild1 && *cpmxchild1 ) || ( cpmxchild0 && *cpmxchild0 ) ) reporterr(       " h" );
//					reporterr(       "%d-%d", clus1, clus2 );
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					break;
				case( 'd' ):
					if( 1 && clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = D__align_ls( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					}
					break;
				case( 'A' ):
					if( clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						if( l < 500 || l % 100 == 0 ) if( ( cpmxchild1 && *cpmxchild1 ) || ( cpmxchild0 && *cpmxchild0 ) ) reporterr(       " h" );
						pscore = A__align( dynamicmtx, penalty, penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, -1, -1, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					}
					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}
		if( targ->lz ) // work/prog3: record what the step left pending (lazygap.h)
		{
			lzstate_t *z = targ->lz;
			int c0 = dep[l].child0, c1 = dep[l].child1, s;
			if( lz_on )
			{
				if( !A_lz.used ) ErrorExit( "A__align did not take the lazy-gap request.\n" );
				for( s=0; s<2; s++ )
				{
					int c = s ? c1 : c0;
					if( A_lz.outgt[s] )
					{
						lzop_t *op = malloc( sizeof( lzop_t ) );
						op->gt = A_lz.outgt[s];
						op->lin = s ? len2 : len1;
						op->lout = A_lz.outlen;
						if( c == -1 || z->op[c] ) ErrorExit( "lazy gaps: bad table owner.\n" );
						z->op[c] = op;
					}
				}
#ifdef PROGTIMING
				if( z->shadow )
				{
					__atomic_add_fetch( &lz_ndefer, ( A_lz.outgt[0] != NULL ) + ( A_lz.outgt[1] != NULL ), __ATOMIC_RELAXED );
					if( lz_defer0 ) { if( A_lz.outgt[0] ) lz_shadow_ins( z, localmem[0], A_lz.outgt[0], len1, A_lz.outlen ); }
					else lz_shadow_set( z, aseq, localmem[0] );
					if( lz_defer1 ) { if( A_lz.outgt[1] ) lz_shadow_ins( z, localmem[1], A_lz.outgt[1], len2, A_lz.outlen ); }
					else lz_shadow_set( z, aseq, localmem[1] );
				}
#endif
				z->len[l] = A_lz.outlen;
				z->stale[l] = lz_stale0 || lz_stale1 || ( c0 != -1 && z->op[c0] ) || ( c1 != -1 && z->op[c1] );
			}
			else
			{
				z->len[l] = strlen( mseq1[0] );
				z->stale[l] = 0;
#ifdef PROGTIMING
				if( z->shadow ) { lz_shadow_set( z, aseq, localmem[0] ); lz_shadow_set( z, aseq, localmem[1] ); }
#endif
			}
			if( l == njob-2 && z->stale[l] ) // the root: every pending table goes in now
			{
				lz_materialize( z, l, dep, topol, aseq, targ->rowcap, localcopylen + 11 );
				PT_LZMAT( localmem[0] );
				PT_LZMAT( localmem[1] );
			}
		}
#ifdef PROGTIMING
		PTT(4);
		if( ptlog ) { int k; for( k=0; k<PT_NACC; k++ ) ptlog[l].acc[k] = pt_acc[k]; ptlog[l].acc[41] = pt_growt; pt_growt = 0.0; ptlog[l].l = l; ptlog[l].thread = thread_no; ptlog[l].clus1 = clus1; ptlog[l].clus2 = clus2; ptlog[l].len1 = len1; ptlog[l].len2 = len2; ptlog[l].c0 = dep[l].child0; ptlog[l].c1 = dep[l].child1; }
#endif
#if SCOREOUT
		reporterr(       "score = %10.2f\n", pscore );
#endif
		tscore += pscore;
		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

		if( disp ) display( localcopy, njob );

		if( newdistmtx ) // tsukawanai
		{
#if 0
			reporterr( "group1 = " );
			for( i=0; i<clus1; i++ ) reporterr( "%d ", topol[l][0][i] );
			reporterr( "\n" );
			reporterr( "group2 = " );
			for( m=0; m<clus2; m++ ) reporterr( "%d ", topol[l][1][m] );
			reporterr( "\n" );
#endif
#if SKIP
			skiptable1 = AllocateIntMtx( clus1, 0 );
			skiptable2 = AllocateIntMtx( clus2, 0 );
			makeskiptable( clus1, skiptable1, mseq1 ); // allocate suru.
			makeskiptable( clus2, skiptable2, mseq2 ); // allocate suru.
#endif
			for( i=0; i<clus1; i++ ) 
			{
				ti = localmem[0][i];
				ssi = selfscore[localmem[0][i]];
				for( m=0; m<clus2; m++ )
				{
					ssm = selfscore[localmem[1][m]];
					tm = localmem[1][m];
					if( ti<tm )
					{
						immin = ti;
						immax = tm;
					}
					else
					{
						immin = tm;
						immax = ti;
					}
					bunbo = MIN( ssi, ssm );
					if( bunbo == 0 )
						newdistmtx[immin][immax-immin] = 2.0; // 2013/Oct/17
					else
#if SKIP
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscorefast( mseq1[i], mseq2[m], skiptable1[i], skiptable2[m], penalty_dist ) / bunbo ) * 2.0;
#else
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscore11( mseq1[i], mseq2[m], penalty_dist ) / bunbo ) * 2.0;
#endif
				}
			}
#if SKIP
			FreeIntMtx( skiptable1 ); skiptable1 = NULL;
			FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif
		}






		/* work/prog: write the aligned rows back before the step is marked done (the parent reads
		   them only after that), under the aseq read lock instead of the global mutex. */
		if( !targ->nocopy ) // work/prog3: the rows were aligned in place
		{
		pthread_rwlock_rdlock( targ->aseqlock );
		if( clus1 + clus2 > 256 ) // work/prog2: see the copy-in above
		{
			int newlen = strlen( localcopy[localmem[0][0]] );
			rowcopy_all( aseq, localcopy, localmem[0], clus1, newlen );
			rowcopy_all( aseq, localcopy, localmem[1], clus2, newlen );
		}
		else
		{
		for( i=0; (j=localmem[0][i])!=-1; i++ )
			strcpy( aseq[j], localcopy[j] );
		for( i=0; (j=localmem[1][i])!=-1; i++ )
			strcpy( aseq[j], localcopy[j] );
		}
		pthread_rwlock_unlock( targ->aseqlock );
		}

		pthread_mutex_lock( targ->mutex );
		dep[l].done = 1;
		if( targ->sched ) tbsched_done( targ->sched, dep, l );
		if( heavybytes > 0.0 ) // work/largen: budget back, deferred heavy steps are ready again
		{
			targ->sched->heavyinflight -= heavybytes;
			while( targ->sched->ndeferred > 0 ) tbsched_push( targ->sched, targ->sched->deferred[--targ->sched->ndeferred] );
		}
		PTT(5);
		(*nrunpt)--;
		pthread_cond_broadcast( targ->treecond );

//		reporterr( "at step %d\n", l );
//		use_getrusage();

		pthread_mutex_unlock( targ->mutex );
		PTT(6);
		if( tb_keepbytes ) // work/largen: memory as with one step at a time (see disttbfast())
		{
			A__align_trim( tb_keepbytes );
			Falign_trim( tb_keepbytes );
		}



		for( i=0; (j=localmem[0][i])!=-1; i++ ) localcopy[j] = NULL; // work/prog2: rows live in lcblk
		for( i=0; (j=localmem[1][i])!=-1; i++ ) localcopy[j] = NULL;


//		if( topol[l][0] ) free( topol[l][0] );
//		topol[l][0] = NULL;
//		if( topol[l][1] ) free( topol[l][1] );
//		topol[l][1] = NULL;
//		if( topol[l] ) free( topol[l] );
//		topol[l] = NULL;


//		reporterr(       "\n" );

		PTT(7);
		free( localmem[0] );
		free( localmem[1] );
	}
#if SCOREOUT
	reporterr(       "totalscore = %10.2f\n\n", tscore );
#endif
}
#endif // enablemultithread

static int dooneiteration( int *nlen, char **aseq, int nadd, char *mergeoralign, char **mseq1, char **mseq2, int ***topol, Treedep *dep, int **memhist, double ***cpmxhist, double *effarr, double **newdistmtx, int *selfscore, ExtAnch *extanch, int **anchindex, int *alloclen, int (*callback)(int, int, char*) )
{
	int l, ll, len1, len2, i, j;
	int clus1, clus2;
	double pscore;
	char *indication1 = NULL, *indication2 = NULL;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
	int *fftlog = NULL; // fixed at 2006/07/26
//	double dumfl = 0.0;
	double dumdb = 0.0;
	int ffttry;
	int m1, m2;
	int *alreadyaligned = NULL;
	double **dynamicmtx = NULL;
	int **localmem = NULL;
	double ***cpmxchild0, ***cpmxchild1;
	double orieff1, orieff2;
	double oscore, nscore;
	ExtAnch *pairanch;
	char **oseq1, **oseq2;
#if SKIP
	int **skiptable1 = NULL, **skiptable2 = NULL;
#endif
#if 0
	int i, j;
#endif


	if( effarr1 == NULL ) 
	{
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
		fftlog = AllocateIntVec( njob );
		alreadyaligned = AllocateIntVec( njob );
		if( specificityconsideration )
			dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );
		localmem = calloc( sizeof( int * ), 2 );
	}
	for( i=0; i<njob-nadd; i++ ) alreadyaligned[i] = 1;
	for( i=njob-nadd; i<njob; i++ ) alreadyaligned[i] = 0;

	if( callback && callback( 0, 50, "Progressive alignment" ) ) goto chudan_tbfast;

	for( l=0; l<njob; l++ ) fftlog[l] = 1;

#if 0 // chain you
	localmem[0][0] = -1;
	localmem[1][0] = -1;
	clus1 = 1;// chain ni hitsuyou
#endif

#if 0
	reporterr(       "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		reporterr(       "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

//	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );


//	writePre( njob, name, nlen, aseq, 0 );

	for( ll=0; ll<njob*ITERATIVECYCLE; ll++ )
	{
		l = ll % njob;
		cpmxchild0 = NULL;
		cpmxchild1 = NULL;

		localmem[0] = calloc( sizeof( int ), 2 );
		localmem[0][0] = l;
		localmem[0][1] = -1;
		clus1 = 1;
		m1 = localmem[0][0];

		localmem[1] = calloc( sizeof( int ), njob );
		for( i=0,j=0; i<njob; i++ )
			if( i != l ) localmem[1][j++] = i;
		localmem[1][j] = -1;
		clus2 = njob-1;
		m2 = localmem[1][0];

//		reporterr(       "\ndistfromtip = %f\n", dep[l].distfromtip );
		if( specificityconsideration )
			makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip );
		else
			dynamicmtx = n_dis_consweight_multi;
//		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, ( dep[l].distfromtip - 0.2 ) * 3 );

		len1 = strlen( aseq[m1] );
		len2 = strlen( aseq[m2] );
		if( *alloclen < len1 + len2 )
		{
			reporterr(       "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  ); 
		}

#if 1 // CHUUI@@@@
		clus1 = fastconjuction_noname( localmem[0], aseq, mseq1, effarr1, effarr, indication1, 0.0, &orieff1 );
		clus2 = fastconjuction_noname( localmem[1], aseq, mseq2, effarr2, effarr, indication2, 0.0, &orieff2 );
#else
		clus1 = fastconjuction_noname( topol[l][0], aseq, mseq1, effarr1, effarr, indication1, 0.0 );
		clus2 = fastconjuction_noname( topol[l][1], aseq, mseq2, effarr2, effarr, indication2, 0.0 );
//		clus1 = fastconjuction_noweight( topol[l][0], aseq, mseq1, effarr1,  indication1 );
//		clus2 = fastconjuction_noweight( topol[l][1], aseq, mseq2, effarr2,  indication2 );
#endif

		intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, strlen( mseq1[0] ), &oscore );

		oseq1 = AllocateCharMtx( 1, len1+1 );
		oseq2 = AllocateCharMtx( njob-1, len2+1 );
		for( i=0; i<clus1; i++ ) strcpy( oseq1[i], mseq1[i] );
		for( i=0; i<clus2; i++ ) strcpy( oseq2[i], mseq2[i] );

		newgapstr = "-";
		commongappick( clus2, mseq2 );
		commongappick( clus1, mseq1 );

		if( l < 500 || l % 100 == 0 ) reporterr(       "\rIteration % 5d / %d ", ll+1, njob*ITERATIVECYCLE );
		if( callback && callback( 0, 50+50*l/(njob-1), "Progressive alignment" ) ) goto chudan_tbfast;
#if 0
		reporterr( "\nclus1=%d, clus2=%d\n", clus1, clus2 );
		for( i=0; i<clus1; i++ ) reporterr( "effarr1[%d]=%f\n", i, effarr1[i] );
		for( i=0; i<clus2; i++ ) reporterr( "effarr2[%d]=%f\n", i, effarr2[i] );
#endif

#if 0
		reporterr(       "STEP %d /%d\n", l+1, njob-1 );
		reporterr(       "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
		reporterr(       "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
#endif

/*
		reporterr(       "before align all\n" );
		display( aseq, njob );
		reporterr(       "\n" );
		reporterr(       "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		reporterr(       "\n" );
		reporterr(       "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		reporterr(       "\n" );
*/

		if( !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000  ) ) )
		{
			reporterr(       "\nlen1=%d, len2=%d, Switching to the memsave mode\n", len1, len2 );
			alg = 'M';
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}

//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000); // v6.708
//		reporterr(       "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (double)len1/fftlog[m1], clus1, (double)len2/fftlog[m2], clus2 );

		if( useexternalanchors )
		{
			pickpairanch( &pairanch, extanch, anchindex, clus1, clus2, localmem[0], localmem[1], mseq1, mseq2 );
//			reporterr( "pairanch: %d:%d\n", pairanch[0].starti, pairanch[0].startj );
			pscore = Falign_givenanchors( pairanch, NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			free( pairanch );
			pairanch = NULL;
		}
		else if( force_fft || ( use_fft && ffttry ) )
		{
			if( l < 500 || l % 100 == 0 ) reporterr(       " f\b\b" );
			if( alg == 'M' )
			{
				if( l < 500 || l % 100 == 0 ) reporterr(       "m" );
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
			{
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
//				reporterr(       "######### mseq1[0] = %s\n", mseq1[0] );
			}
		}
		else
		{
			if( l < 500 || l % 100 == 0 ) reporterr(       " d\b\b" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					if( l < 500 || l % 100 == 0 ) reporterr(       "m" );
					if( l < 500 || l % 100 == 0 ) if( ( cpmxchild1 && *cpmxchild1 ) || ( cpmxchild0 && *cpmxchild0 ) ) reporterr(       " h" );
//					reporterr(       "%d-%d", clus1, clus2 );
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					break;
				case( 'd' ):
					if( 1 && clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = D__align_ls( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					}
					break;
				case( 'A' ):
					if( clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
						if( l < 500 || l % 100 == 0 ) if( ( cpmxchild1 && *cpmxchild1 ) || ( cpmxchild0 && *cpmxchild0 ) ) reporterr(       " h" );
//						reporterr(       "\n\n %d - %d (%d x %d) : \n", topol[l][0][0], topol[l][1][0], clus1, clus2 );
						pscore = A__align( dynamicmtx, penalty, penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, localmem[0][0], 1, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					}

					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}
		intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, strlen( mseq1[0] ), &nscore );
#if SCOREOUT
		reporterr(       "score = %10.2f\n", pscore );
#endif
		if( nscore < oscore )
		{
			for( i=0; i<clus1; i++ ) strcpy( mseq1[i], oseq1[i] );
			for( i=0; i<clus2; i++ ) strcpy( mseq2[i], oseq2[i] );
		}
		FreeCharMtx( oseq1 );
		FreeCharMtx( oseq2 );

		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

//		writePre( njob, name, nlen, aseq, 0 );

		if( disp ) display( aseq, njob );
//		reporterr(       "\n" );



#if 0
		if( localmem[1][0] == 13 ) 
		{
			reporterr( "OUTPUT!\n" );
			for( i=0; i<clus1; i++ ) reporterr( ">g1\n%s\n", mseq1[i] );
			for( i=0; i<clus2; i++ ) reporterr( ">g2\n%s\n", mseq2[i] );
			exit( 1 );
		}
#endif

//		free( topol[l][0] ); topol[l][0] = NULL;
//		free( topol[l][1] ); topol[l][1] = NULL;
//		free( topol[l] ); topol[l] = NULL;


//		reporterr(       ">514\n%s\n", aseq[514] );
		free( localmem[0] );
		free( localmem[1] );

	}

	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	Falign_givenanchors( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
	D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
	free( effarr1 );
	free( effarr2 );
	free( indication1 );
	free( indication2 );
	free( fftlog );
	if( specificityconsideration )
		FreeDoubleMtx( dynamicmtx );
	free( alreadyaligned );
	free( localmem );
	effarr1 = NULL;
	return( 0 );

	chudan_tbfast:

	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	Falign_givenanchors( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
	D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
	if( effarr1 ) free( effarr1 ); effarr1 = NULL;
	if( effarr2 ) free( effarr2 ); effarr2 = NULL;
	if( indication1 ) free( indication1 ); indication1 = NULL;
	if( indication2 ) free( indication2 ); indication2 = NULL;
	if( fftlog ) free( fftlog ); fftlog = NULL;
	if( alreadyaligned ) free( alreadyaligned ); alreadyaligned = NULL;
	if( specificityconsideration )
	{
		if( dynamicmtx ) FreeDoubleMtx( dynamicmtx ); dynamicmtx = NULL;
	}
	if( localmem ) free( localmem ); localmem = NULL;
#if SKIP
	if( skiptable1 ) FreeIntMtx( skiptable1 ); skiptable1 = NULL;
	if( skiptable2 ) FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif

	return( 1 );
}
static int treebase( int *nlen, char **aseq, int nadd, char *mergeoralign, char **mseq1, char **mseq2, int ***topol, Treedep *dep, int **memhist, double ***cpmxhist, double *effarr, double **newdistmtx, int *selfscore, ExtAnch *extanch, int **anchindex, int *alloclen, int (*callback)(int, int, char*) )
{
	int l, len1, len2, i, m, immin, immax;
	int len1nocommongap, len2nocommongap;
	int clus1, clus2;
	double pscore, tscore;
	char *indication1 = NULL, *indication2 = NULL;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
	int *fftlog = NULL; // fixed at 2006/07/26
//	double dumfl = 0.0;
	double dumdb = 0.0;
	int ffttry;
	int m1, m2;
	int *gaplen = NULL;
	int *gapmap = NULL;
	int *alreadyaligned = NULL;
	double **dynamicmtx = NULL;
	double ssi, ssm, bunbo;
	int tm, ti;
	int gapmaplen;
	int **localmem = NULL;
	double ***cpmxchild0, ***cpmxchild1;
	double orieff1, orieff2;
	ExtAnch *pairanch;
#if SKIP
	int **skiptable1 = NULL, **skiptable2 = NULL;
#endif
#if 0
	int i, j;
#endif


	if( effarr1 == NULL ) 
	{
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
		fftlog = AllocateIntVec( njob );
		gaplen = AllocateIntVec( *alloclen+10 );
		gapmap = AllocateIntVec( *alloclen+10 );
		alreadyaligned = AllocateIntVec( njob );
		if( specificityconsideration )
			dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );
		localmem = calloc( sizeof( int * ), 2 );
	}
	for( i=0; i<njob-nadd; i++ ) alreadyaligned[i] = 1;
	for( i=njob-nadd; i<njob; i++ ) alreadyaligned[i] = 0;

	if( callback && callback( 0, 50, "Progressive alignment" ) ) goto chudan_tbfast;

	for( l=0; l<njob; l++ ) fftlog[l] = 1;

#if 0 // chain you
	localmem[0][0] = -1;
	localmem[1][0] = -1;
	clus1 = 1;// chain ni hitsuyou
#endif

#if 0
	reporterr(       "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		reporterr(       "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

//	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );


//	writePre( njob, name, nlen, aseq, 0 );

	tscore = 0.0;
	for( l=0; l<njob-1; l++ )
	{
		m1 = topol[l][0][0];
		m2 = topol[l][1][0];
//		reporterr( " at the beginning of the loop, clus1,clus2=%d,%d\n", clus1, clus2 );

//		reporterr( "l=%d, dep[l].child0=%d, dep[l].child1=%d\n", l, dep[l].child0, dep[l].child1 );
		if( dep[l].child0 == -1 ) cpmxchild0 = NULL; else cpmxchild0 = cpmxhist+dep[l].child0;
		if( dep[l].child1 == -1 ) cpmxchild1 = NULL; else cpmxchild1 = cpmxhist+dep[l].child1;
//		reporterr( "cpmxchild0=%p, cpmxchild1=%p\n", cpmxchild0, cpmxchild1 );

#if 0
		if(  l > 0 && dep[l].child0 == l-1 && dep[l].child1 == -1 && dep[dep[l].child0].child1 == -1 )
		{
			localmem[0][clus1] = topol[l-1][1][0];
			localmem[0][clus1+1] = -1;

			localmem[1][0] = topol[l][1][0];
			localmem[1][1] = -1;
		}
		else
		{
			localmem[0][0] = -1;
			posinmem = topolorderz( localmem[0], topol, dep, l, 0 ) - localmem[0];
			localmem[1][0] = -1;
			posinmem = topolorderz( localmem[1], topol, dep, l, 1 ) - localmem[1];
		}
#else
		if( dep[l].child0 == -1 ) 
		{
			localmem[0] = calloc( sizeof( int ), 2 );
			localmem[0][0] = m1;
			localmem[0][1] = -1;
			clus1 = 1;
		}
		else
		{
			localmem[0] = memhist[dep[l].child0];
			clus1 = intlen( localmem[0] );
		}
		if( dep[l].child1 == -1 ) 
		{
			localmem[1] = calloc( sizeof( int ), 2 );
			localmem[1][0] = m2;
			localmem[1][1] = -1;
			clus2 = 1;
		}
		else
		{
			localmem[1] = memhist[dep[l].child1];
			clus2 = intlen( localmem[1] );
		}

		if( l != njob-2 )
		{
			memhist[l] = calloc( sizeof( int ), clus1+clus2+1 );
			intcpy( memhist[l], localmem[0] );
			intcpy( memhist[l]+clus1, localmem[1] );
			memhist[l][clus1+clus2] = -1;
		}
#endif

		if( mergeoralign[l] == 'n' )
		{
//			reporterr(       "SKIP!\n" );
//			free( topol[l][0] ); topol[l][0] = NULL;
//			free( topol[l][1] ); topol[l][1] = NULL;
//			free( topol[l] ); topol[l] = NULL;
			free( localmem[0] );
			free( localmem[1] );
			continue;
		}

//		reporterr(       "\ndistfromtip = %f\n", dep[l].distfromtip );
		if( specificityconsideration )
			makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip );
		else
			dynamicmtx = n_dis_consweight_multi;
//		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, ( dep[l].distfromtip - 0.2 ) * 3 );


		len1 = strlen( aseq[m1] );
		len2 = strlen( aseq[m2] );
		if( *alloclen < len1 + len2 )
		{
			reporterr(       "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  ); 
			gaplen = realloc( gaplen, ( *alloclen + 10 ) * sizeof( int ) );
			if( gaplen == NULL )
			{
				reporterr(       "Cannot realloc gaplen\n" );
				exit( 1 );
			}
			gapmap = realloc( gapmap, ( *alloclen + 10 ) * sizeof( int ) );
			if( gapmap == NULL )
			{
				reporterr(       "Cannot realloc gapmap\n" );
				exit( 1 );
			}
			reporterr(       "done. *alloclen = %d\n", *alloclen );
		}

#if 1 // CHUUI@@@@
		clus1 = fastconjuction_noname( localmem[0], aseq, mseq1, effarr1, effarr, indication1, 0.0, &orieff1 );
		clus2 = fastconjuction_noname( localmem[1], aseq, mseq2, effarr2, effarr, indication2, 0.0, &orieff2 );
#else
		clus1 = fastconjuction_noname( topol[l][0], aseq, mseq1, effarr1, effarr, indication1, 0.0 );
		clus2 = fastconjuction_noname( topol[l][1], aseq, mseq2, effarr2, effarr, indication2, 0.0 );
//		clus1 = fastconjuction_noweight( topol[l][0], aseq, mseq1, effarr1,  indication1 );
//		clus2 = fastconjuction_noweight( topol[l][1], aseq, mseq2, effarr2,  indication2 );
#endif











		if( mergeoralign[l] == '1' || mergeoralign[l] == '2' )
		{
			newgapstr = "=";
		}
		else
			newgapstr = "-";

		len1nocommongap = len1;
		len2nocommongap = len2;
		if( mergeoralign[l] == '1' ) // nai
		{
			findcommongaps( clus2, mseq2, gapmap );
			commongappick( clus2, mseq2 );
			len2nocommongap = strlen( mseq2[0] );
		}
		else if( mergeoralign[l] == '2' )
		{
			findcommongaps( clus1, mseq1, gapmap );
			commongappick( clus1, mseq1 );
			len1nocommongap = strlen( mseq1[0] );
		}

#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
#endif


#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
#endif


//		fprintf( trap_g, "\nSTEP-%d\n", l );
//		fprintf( trap_g, "group1 = %s\n", indication1 );
//		fprintf( trap_g, "group2 = %s\n", indication2 );

//		reporterr(       "\rSTEP % 5d / %d %d-%d", l+1, njob-1, clus1, clus2 );
		if( l < 500 || l % 100 == 0 ) reporterr(       "\rSTEP % 5d / %d ", l+1, njob-1 );
		if( callback && callback( 0, 50+50*l/(njob-1), "Progressive alignment" ) ) goto chudan_tbfast;
#if 0
		reporterr( "\nclus1=%d, clus2=%d\n", clus1, clus2 );
#endif

#if 0
		reporterr(       "STEP %d /%d\n", l+1, njob-1 );
		reporterr(       "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
		reporterr(       "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
#endif

/*
		reporterr(       "before align all\n" );
		display( aseq, njob );
		reporterr(       "\n" );
		reporterr(       "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		reporterr(       "\n" );
		reporterr(       "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		reporterr(       "\n" );
*/


		if( !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000  ) ) )
		{
			reporterr(       "\nlen1=%d, len2=%d, Switching to the memsave mode\n", len1, len2 );
			alg = 'M';
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}

//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		else						   ffttry = 0;
		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		reporterr(       "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (double)len1/fftlog[m1], clus1, (double)len2/fftlog[m2], clus2 );
//		reporterr( "fftlog=%d,%d, ffttry=%d\n", fftlog[m1], fftlog[m2], ffttry );

		if( useexternalanchors )
		{
			pickpairanch( &pairanch, extanch, anchindex, clus1, clus2, localmem[0], localmem[1], mseq1, mseq2 );
//			reporterr( "pairanch: %d:%d\n", pairanch[0].starti, pairanch[0].startj );
			pscore = Falign_givenanchors( pairanch, NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			free( pairanch );
			pairanch = NULL;
		}
		else if( force_fft || ( use_fft && ffttry ) )
		{
			if( l < 500 || l % 100 == 0 ) reporterr(       " f\b\b" );
			if( alg == 'M' )
			{
				if( l < 500 || l % 100 == 0 ) reporterr(       "m" );
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
			{
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
//				reporterr(       "######### mseq1[0] = %s\n", mseq1[0] );
			}
		}
		else
		{
			if( l < 500 || l % 100 == 0 ) reporterr(       " d\b\b" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					if( l < 500 || l % 100 == 0 ) reporterr(       "m" );
					if( l < 500 || l % 100 == 0 ) if( ( cpmxchild1 && *cpmxchild1 ) || ( cpmxchild0 && *cpmxchild0 ) ) reporterr(       " h" );
//					reporterr(       "%d-%d", clus1, clus2 );
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					break;
				case( 'd' ):
					if( 1 && clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = D__align_ls( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					}
					break;
				case( 'A' ):
					if( clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
						if( l < 500 || l % 100 == 0 ) if( ( cpmxchild1 && *cpmxchild1 ) || ( cpmxchild0 && *cpmxchild0 ) ) reporterr(       " h" );
//						reporterr(       "\n\n %d - %d (%d x %d) : \n", topol[l][0][0], topol[l][1][0], clus1, clus2 );
						pscore = A__align( dynamicmtx, penalty, penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, localmem[0][0], 1, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					}

					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}
#if SCOREOUT
		reporterr(       "score = %10.2f\n", pscore );
#endif
		tscore += pscore;
		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

//		writePre( njob, name, nlen, aseq, 0 );

		if( disp ) display( aseq, njob );
//		reporterr(       "\n" );

		if( mergeoralign[l] == '1' ) // jissainiha nai. atarashii hairetsu ha saigo dakara.
		{
			reporterr( "Check source!!!\n" );
			exit( 1 );
		}
		if( mergeoralign[l] == '2' )
		{
//			if( localkeeplength ) ndeleted += deletenewinsertions( clus1, clus2, mseq1, mseq2, NULL );
//			for( i=0; i<clus1; i++ ) reporterr(       ">STEP0 mseq1[%d] = \n%s\n", i, mseq1[i] );
//			for( i=0; i<clus2; i++ ) reporterr(       ">STEP0 mseq2[%d] = \n%s\n", i, mseq2[i] );
			gapmaplen = strlen( mseq1[0] )-len1nocommongap+len1;
			adjustgapmap( gapmaplen, gapmap, mseq1[0] );
#if 0
			reporterr( "\n" );
			for( i=0; i<clus1; i++ ) reporterr(       ">STEP1 mseq1[%d] = \n%s\n", i, mseq1[i] );
			for( i=0; i<clus2; i++ ) reporterr(       ">STEP1 mseq2[%d] = \n%s\n", i, mseq2[i] );
#endif
//			if( clus1 + clus2 < njob ) restorecommongaps( njob, aseq, topol[l][0], topol[l][1], gapmap, *alloclen, '-' );
			if( smoothing )
			{
				restorecommongapssmoothly( njob, njob-(clus1+clus2), aseq, localmem[0], localmem[1], gapmap, *alloclen, '-' );
				findnewgaps( clus1, 0, mseq1, gaplen );
				insertnewgaps_bothorders( njob, alreadyaligned, aseq, localmem[0], localmem[1], gaplen, gapmap, gapmaplen, *alloclen, alg, '-' );
			}
			else
			{
				restorecommongaps( njob, njob-(clus1+clus2), aseq, localmem[0], localmem[1], gapmap, *alloclen, '-' );
				findnewgaps( clus1, 0, mseq1, gaplen );
				insertnewgaps( njob, alreadyaligned, aseq, localmem[0], localmem[1], gaplen, gapmap, *alloclen, alg, '-' );
			}

#if 0
			reporterr( "\n" );
			for( i=0; i<clus1; i++ ) reporterr(       ">STEP3 mseq1[%d] = \n%s\n", i, mseq1[i] );
			for( i=0; i<clus2; i++ ) reporterr(       ">STEP3 mseq2[%d] = \n%s\n", i, mseq2[i] );
#endif

#if 0
			for( i=0; i<njob; i++ ) eq2dash( aseq[i] );
			for( i=0; i<clus1; i++ ) 
			{
				reporterr( "mseq1[%d] bef change = %s\n", i, mseq1[i] );
				eq2dash( mseq1[i] );
				reporterr( "mseq1[%d] aft change = %s\n", i, mseq1[i] );
			}
			for( i=0; i<clus2; i++ ) 
			{
				reporterr( "mseq2[%d] bef change = %s\n", i, mseq2[i] );
				eq2dash( mseq2[i] );
				reporterr( "mseq2[%d] aft change = %s\n", i, mseq2[i] );
			}
			for( i=0; i<clus1; i++ ) eq2dash( mseq1[i] );
			for( i=0; i<clus2; i++ ) eq2dash( mseq2[i] );
#endif


			eq2dashmatometehayaku( mseq1, clus1 );
			eq2dashmatometehayaku( mseq2, clus2 );

			for( i=0; (m=localmem[1][i])>-1; i++ ) alreadyaligned[m] = 1;
		}

		if( newdistmtx ) // tsukawanai
		{
#if 0
			reporterr( "group1 = " );
			for( i=0; i<clus1; i++ ) reporterr( "%d ", topol[l][0][i] );
			reporterr( "\n" );
			reporterr( "group2 = " );
			for( m=0; m<clus2; m++ ) reporterr( "%d ", topol[l][1][m] );
			reporterr( "\n" );
#endif
#if SKIP
			skiptable1 = AllocateIntMtx( clus1, 0 );
			skiptable2 = AllocateIntMtx( clus2, 0 );
			makeskiptable( clus1, skiptable1, mseq1 ); // allocate suru.
			makeskiptable( clus2, skiptable2, mseq2 ); // allocate suru.
#endif
			for( i=0; i<clus1; i++ ) 
			{
#if SKIP
//				makeskiptable( 1, skiptable1, mseq1+i ); // allocate suru.
#endif
				ti = localmem[0][i];
				ssi = selfscore[localmem[0][i]];
				for( m=0; m<clus2; m++ )
				{
					ssm = selfscore[localmem[1][m]];
					tm = localmem[1][m];
					if( ti<tm )
					{
						immin = ti;
						immax = tm;
					}
					else
					{
						immin = tm;
						immax = ti;
					}
					bunbo = MIN( ssi, ssm );
					if( bunbo == 0.0 )
						newdistmtx[immin][immax-immin] = 2.0; // 2013/Oct/17
					else
#if SKIP
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscorefast( mseq1[i], mseq2[m], skiptable1[i], skiptable2[m], penalty_dist ) / bunbo ) * 2.0;
#else
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscore11( mseq1[i], mseq2[m], penalty_dist ) / bunbo ) * 2.0;
#endif
				}
			}
#if SKIP
			FreeIntMtx( skiptable1 ); skiptable1 = NULL;
			FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif
		}

#if 0
		if( localmem[1][0] == 13 ) 
		{
			reporterr( "OUTPUT!\n" );
			for( i=0; i<clus1; i++ ) reporterr( ">g1\n%s\n", mseq1[i] );
			for( i=0; i<clus2; i++ ) reporterr( ">g2\n%s\n", mseq2[i] );
			exit( 1 );
		}
#endif

//		free( topol[l][0] ); topol[l][0] = NULL;
//		free( topol[l][1] ); topol[l][1] = NULL;
//		free( topol[l] ); topol[l] = NULL;


//		reporterr(       ">514\n%s\n", aseq[514] );
		free( localmem[0] );
		free( localmem[1] );
	}

#if SCOREOUT
	reporterr(       "totalscore = %10.2f\n\n", tscore );
#endif
	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	Falign_givenanchors( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
	D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
	free( effarr1 );
	free( effarr2 );
	free( indication1 );
	free( indication2 );
	free( fftlog );
	free( gaplen );
	free( gapmap );
	if( specificityconsideration )
		FreeDoubleMtx( dynamicmtx );
	free( alreadyaligned );
	free( localmem );
	effarr1 = NULL;
	return( 0 );

	chudan_tbfast:

	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	Falign_givenanchors( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
	D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
	if( effarr1 ) free( effarr1 ); effarr1 = NULL;
	if( effarr2 ) free( effarr2 ); effarr2 = NULL;
	if( indication1 ) free( indication1 ); indication1 = NULL;
	if( indication2 ) free( indication2 ); indication2 = NULL;
	if( fftlog ) free( fftlog ); fftlog = NULL;
	if( gaplen ) free( gaplen ); gaplen = NULL;
	if( gapmap ) free( gapmap ); gapmap = NULL;
	if( alreadyaligned ) free( alreadyaligned ); alreadyaligned = NULL;
	if( specificityconsideration )
	{
		if( dynamicmtx ) FreeDoubleMtx( dynamicmtx ); dynamicmtx = NULL;
	}
	if( localmem ) free( localmem ); localmem = NULL;
#if SKIP
	if( skiptable1 ) FreeIntMtx( skiptable1 ); skiptable1 = NULL;
	if( skiptable2 ) FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif

	return( 1 );
}

static void WriteOptions( FILE *fp )
{

	if( dorp == 'd' ) fprintf( fp, "DNA\n" );
	else
	{
		if     ( scoremtx ==  0 ) fprintf( fp, "JTT %dPAM\n", pamN );
		else if( scoremtx ==  1 ) fprintf( fp, "BLOSUM %d\n", nblosum );
		else if( scoremtx ==  2 ) fprintf( fp, "M-Y\n" );
	}
    reporterr(       "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );
    if( use_fft ) fprintf( fp, "FFT on\n" );

	fprintf( fp, "tree-base method\n" );
	if( tbrweight == 0 ) fprintf( fp, "unweighted\n" );
	else if( tbrweight == 3 ) fprintf( fp, "clustalw-like weighting\n" );
	if( tbitr || tbweight ) 
	{
		fprintf( fp, "iterate at each step\n" );
		if( tbitr && tbrweight == 0 ) fprintf( fp, "  unweighted\n" ); 
		if( tbitr && tbrweight == 3 ) fprintf( fp, "  reversely weighted\n" ); 
		if( tbweight ) fprintf( fp, "  weighted\n" ); 
		fprintf( fp, "\n" );
	}

   	 fprintf( fp, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );

	if( alg == 'a' )
		fprintf( fp, "Algorithm A\n" );
	else if( alg == 'A' ) 
		fprintf( fp, "Algorithm A+\n" );
	else
		fprintf( fp, "Unknown algorithm\n" );

	if( treemethod == 'X' )
		fprintf( fp, "Tree = UPGMA (mix).\n" );
	else if( treemethod == 'E' )
		fprintf( fp, "Tree = UPGMA (average).\n" );
	else if( treemethod == 'q' )
		fprintf( fp, "Tree = Minimum linkage.\n" );
	else
		fprintf( fp, "Unknown tree.\n" );

    if( use_fft )
    {
        fprintf( fp, "FFT on\n" );
        if( dorp == 'd' )
            fprintf( fp, "Basis : 4 nucleotides\n" );
        else
        {
            if( fftscore )
                fprintf( fp, "Basis : Polarity and Volume\n" );
            else
                fprintf( fp, "Basis : 20 amino acids\n" );
        }
        fprintf( fp, "Threshold   of anchors = %d%%\n", fftThreshold );
        fprintf( fp, "window size of anchors = %dsites\n", fftWinSize );
    }
	else
        fprintf( fp, "FFT off\n" );
	fflush( fp );
}

static double **preparepartmtx( int nseq )
{
	int i;
	double **val;
	double size;

	val = (double **)calloc( nseq, sizeof( double *) );;
	size = 0;

	if( compacttree == 1 )
	{
		for( i=0; i<nseq; i++ )
		{
			size += (double)sizeof( double ) * nseq;
			if( size > maxdistmtxsize ) 
			{
				reporterr( "\n\nThe size of full distance matrix is estimated to exceed %.2fGB.\n", maxdistmtxsize / 1000 / 1000 /1000 );
				reporterr( "Will try the calculation using a %d x %d matrix.\n", nseq, i );
				reporterr( "This calculation will be slow due to the limited RAM space.\n", i, nseq );
				reporterr( "To avoid the slowdown, please try '--initialramusage xGB' (x>>%.2f),\n", maxdistmtxsize / 1000 / 1000 /1000 );
				reporterr( "if larger RAM space is available.\n" );
				reporterr( "Note that xGB is NOT the upper limit of RAM usage.\n" );
				reporterr( "Two to three times larger space may be used for building a guide tree.\n" );
				reporterr( "Memory usage of the MSA stage depends on similarity of input sequences.\n\n" );
//				reporterr( "If the RAM is small, try '--initialramusage xGB' with a smaller x value.\n" );
				reporterr( "The '--memsavetree' option uses smaller RAM space.\n" );
				reporterr( "If tree-like relationship can be ignored, try '--pileup' or '--randomchain'.\n\n" );
				reporterr( "The result of --initialramusage xGB is almost identical to the default, except for rounding differences.\n" );

				reporterr( "In the cases of --memsavetree, --pileup and --randomchain, the result differs from the default.\n\n" );
				break;
			}
			val[i] = (double *)calloc( nseq, sizeof( double ) );
		}
		if( i == nseq ) reporterr( "The full matrix will be used.\n" );

		for( ;i<nseq; i++ ) val[i] = NULL; // nen no tame
	}
	else
	{
		for( i=0; i<nseq; i++ ) val[i] = NULL; // nen no tame
	}
	return( val );
}

/* fast-mafft (overhead lever): the 6-mer count -> distance conversion, one row per call.
   The statements are those of the original serial loop in disttbfast(), unchanged (the lenfac
   statement is compiled to the same fmadd); every mtx[i][j-i] depends only on its own count
   and the diagonal entries mtx[i][0], mtx[j][0], which the loop never writes, so rows can be
   converted in any order, or at the same time. */
static void dist6conv_row( double **mtx, int *nogaplen, int njob, int i )
{
	int j;
	double longer, shorter;
	double lenfac;
	double bunbo;
	for( j=i+1; j<njob; j++ ) 
	{
		if( nogaplen[i] > nogaplen[j] )
		{
			longer=(double)nogaplen[i];
			shorter=(double)nogaplen[j];
		}
		else
		{
			longer=(double)nogaplen[j];
			shorter=(double)nogaplen[i];
		}
		lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
		bunbo = MIN( mtx[i][0], mtx[j][0] );
		if( bunbo == 0.0 )
			mtx[i][j-i] = 2.0; // 2013/Oct/17 -> 2bai
		else
			mtx[i][j-i] = ( 1.0 - mtx[i][j-i] / bunbo ) * lenfac * 2.0; // 2013/Oct/17 -> 2bai
	}
}

#ifdef enablemultithread
typedef struct
{
	double **mtx;
	int *nogaplen;
	int njob;
	int *jobpospt;
	pthread_mutex_t *mutex;
} dist6convthread_arg_t;

static void *dist6convthread( void *arg )
{
	dist6convthread_arg_t *targ = (dist6convthread_arg_t *)arg;
	int i, i0, i1;
	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		i0 = *targ->jobpospt;
		i1 = MIN( i0 + 32, targ->njob-1 );
		*targ->jobpospt = i1;
		pthread_mutex_unlock( targ->mutex );
		if( i0 >= targ->njob-1 ) return( NULL );
		for( i=i0; i<i1; i++ ) dist6conv_row( targ->mtx, targ->nogaplen, targ->njob, i );
	}
}
#endif

int disttbfast( int ngui, int lgui, char **namegui, char **seqgui, int argc, char **argv, int (*callback)(int, int, char*))
{
	int  *nlen = NULL;	
	int  *nogaplen = NULL;	
	char **name = NULL, **seq = NULL;
	char **mseq1 = NULL, **mseq2 = NULL;
	char **bseq = NULL;
	double *eff = NULL;
	int i, j;
	int ***topol = NULL;
	int *addmem = NULL;
	Treedep *dep = NULL;
	int **memhist = NULL;
	double ***cpmxhist = NULL;
	double **len = NULL;
	FILE *infp = NULL;
//	FILE *adfp;
	char c;
	int alloclen;
	double longer, shorter;
	double lenfac;
	double bunbo;
	int sextetdone = 0; // fast-mafft (work/sextet)

	FILE *orderfp = NULL, *hat2p = NULL;
	int *grpseq = NULL;
	char *tmpseq = NULL;
	int  **pointt = NULL;
	double **mtx = NULL; // by D. Mathog
	int *table1 = NULL;
	char b[B];
	int ien, nlim;
	int includememberres0, includememberres1;
	double unweightedspscore;
	int alignmentlength;
	char *mergeoralign = NULL;
	int foundthebranch;
	int nsubalignments = 0, maxmem;
	int **subtable = NULL;
	int *insubtable = NULL;
	int *preservegaps = NULL;
	char ***subalnpt = NULL;
	int val;
	char **tmpargv = NULL;
	int iguidetree;
	int *selfscore = NULL;
	int calcpairdists;
	int	**skiptable = NULL;
	char algbackup;
	char *originalgaps = NULL;
	char **addbk = NULL;
	GapPos **deletelist = NULL;
	FILE *dlf = NULL;
	int randomseed;
	int **localmem = NULL;
	int posinmem;
// for compacttree
	int *mindistfrom = NULL;
	double *mindist = NULL;
	double **partmtx = NULL;
// for compacttree
	ExtAnch *extanch = NULL;
	int **anchindex = NULL;


	if( ngui )
	{
		initglobalvariables();
		njob = ngui;
		nlenmax = 0;
		for( i=0; i<njob; i++ )
		{
			ien = strlen( seqgui[i] );
			if( ien > nlenmax ) nlenmax = ien;
		}
		infp = NULL;
//		stderr = fopen( "/dev/null", "a" ); // Windows????
		tmpargv = AllocateCharMtx( argc, 0 );
		for( i=0; i<argc; i++ ) tmpargv[i] = argv[i];
		gmsg = 1;
	}
	else
		gmsg = 0; // iranai

	arguments( argc, argv );
	algbackup = alg; // tbfast wo disttbfast ni ketsugou shitatame.
#ifndef enablemultithread
	nthreadpair = nthread = 0;
#endif


	if( ngui )
	{
		for( i=0; i<argc; i++ ) 
		{
//			free( tmpargv[i] );
			argv[i] = tmpargv[i];
		}
		free( tmpargv );
	}
	else
	{
		if( inputfile )
		{
			infp = fopen( inputfile, "rb" );
			if( !infp )
			{
				reporterr(       "Cannot open %s\n", inputfile );
				exit( 1 );
			}
		}
		else
			infp = stdin;
	
		getnumlen( infp );
		rewind( infp );
	}
	
	if( njob > 10000000 )
	{
		reporterr(       "The number of sequences must be < %d\n", 10000000 );
//		reporterr(       "Please try the --parttree option for such large data.\n" );
		exit( 1 );
	}

	if( njob < 2 )
	{
		seq = AllocateCharMtx( 2, nlenmax*1+1 );
    	name = AllocateCharMtx( 2, B+1 );
	    nlen = AllocateIntVec( 2 ); 
		readData_pointer( infp, name, nlen, seq );
		fclose( infp );
		gappick0( seq[1], seq[0] );
		writeData_pointer( stdout, njob, name, nlen, seq+1 );
		reporterr( "Warning: Only %d sequence found.\n", njob ); 
		FreeCharMtx( seq );
		FreeCharMtx( name );
		free( nlen );
		exit( 0 );
	}

	if( specificityconsideration != 0.0 && nlenmax)
	{
		if( nlenmax > 100000 )
		{
			reporterr( "\n" );
			reporterr( "Too long to apply --allowshift or --unalignlevel>0\n" );
			reporterr( "Please use the normal mode.\n" );
			reporterr( "Please also note that MAFFT does not assume genomic rearrangements.\n" );
			reporterr( "\n" );
			exit( 1 );
		}
	}


#if !defined(mingw) && !defined(_MSC_VER)
	setstacksize( 200 * njob ); // topolorder() de ookime no stack wo shiyou.
#endif

	if( subalignment )
	{
		readsubalignmentstable( njob, NULL, NULL, &nsubalignments, &maxmem );
		reporterr(       "nsubalignments = %d\n", nsubalignments );
		reporterr(       "maxmem = %d\n", maxmem );
		subtable = AllocateIntMtx( nsubalignments, maxmem+1 );
		insubtable = AllocateIntVec( njob );
		preservegaps = AllocateIntVec( njob );
		for( i=0; i<njob; i++ ) insubtable[i] = 0;
		for( i=0; i<njob; i++ ) preservegaps[i] = 0;
		subalnpt = AllocateCharCub( nsubalignments, maxmem, 0 );
		readsubalignmentstable( njob, subtable, preservegaps, NULL, NULL );
	}



	seq = AllocateCharMtx( njob, nlenmax*1+1 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );

	eff = AllocateDoubleVec( njob );
	mergeoralign = AllocateCharVec( njob );

	dep = (Treedep *)calloc( njob, sizeof( Treedep ) );

	if( nadd ) addmem = AllocateIntVec( nadd+1 );

	localmem = AllocateIntMtx( 2, njob+1 );

#if 0
	Read( name, nlen, seq );
	readData( infp, name, nlen, seq );
#else
    name = AllocateCharMtx( njob, B+1 );
    nlen = AllocateIntVec( njob ); 
    nogaplen = AllocateIntVec( njob ); 
	if( ngui )
	{
		if( copydatafromgui( namegui, seqgui, name, nlen, seq ) )
			exit( 1 );
	}
	else
	{
		readData_pointer( infp, name, nlen, seq );
		fclose( infp );
	}
#endif

	if( useexternalanchors ) // nogaplen ha ato de uwagaki sareru kamo
	{
		char *tmpseq = calloc( nlenmax+5, sizeof( char ) );
		for( i=0; i<njob; i++ )
		{
//			reporterr( "i=%d, nlenmax=%d, len=%d\n", i, nlenmax, strlen( seq[i] ) );
			gappick0( tmpseq, seq[i] );
			nogaplen[i] = strlen( tmpseq );
		}
		extanch = calloc( sizeof( ExtAnch ), 1 );
		extanch[0].i=-1;
		extanch[0].j=-1;
		reporterr( "reading anchors\n" );
		readexternalanchors( &extanch, njob, nogaplen ); // allocate sareru
		anchindex = AllocateIntMtx( njob, njob ); // sukoshi muda
		for( i=0; i<njob; i++ ) for( j=0; j<njob; j++ ) anchindex[i][j] = -1;
		reporterr( "sorting anchors\n" );
		indexanchors( extanch, anchindex );
		//checkanchors_internal( extanch ); // comment out -> equivalent to v7.448
		free( tmpseq );
	}

	constants( njob, seq );


#if 0
	reporterr(       "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM();

	initFiles();

	WriteOptions( trap_g );

	c = seqcheck( seq );
	if( c )
	{
		reporterr(       "Illegal character %c\n", c );
		exit( 1 );
	}

	reporterr(       "\n" );

//	reporterr(       "tuplesize = %d, dorp = %c\n", tuplesize, dorp );
	if( dorp == 'p' && tuplesize != 6 )
	{
		reporterr(       "tuplesize must be 6 for aa sequence\n" );
		exit( 1 );
	}
	if( dorp == 'd' && tuplesize != 6 && tuplesize != 10 )
	{
		reporterr(       "tuplesize must be 6 or 10 for dna sequence\n" );
		exit( 1 );
	}

	if( treein )
	{
		int npickx;
		treein = check_guidetreefile( &randomseed, &npickx, &maxdistmtxsize );
		if( treein == 't' )
		{
			varpairscore( njob, npickx, nlenmax, seq, randomseed );
			exit( 1 );
		}
		else if( treein == 'c' )
		{
			compacttree = 1;
			treein = 0;
//			use_fft = 0; // kankeinai?
//			maxdistmtxsize = 5 * 1000 * 1000; // 5GB. ato de kahen ni suru.
//			maxdistmtxsize =  1.0 * 1000 * 1000 * 1000; // 5GB. ato de kahen ni suru.
		}
		else if( treein == 'Y' )
		{
			compacttree = 4; // youngest linkage, 3 ha tbfast de tsukaunode ichiou sakeru
			treein = 0;
//			use_fft = 0; // kankeinai?
		}
		else if( treein == 'S' || treein == 'C' )
		{
			compacttree = 2; // 3 ha tbfast de tsukaunode ichiou sakeru
			treein = 0;
//			use_fft = 0; // kankeinai?
		}
		else if( treein == 'a' )
		{
//			reporterr( "Compute pairwise scores\n" );
			if( njob > 200000 )
			{
				reporterr( "Chain?\n" );
				treein = 's';
				nguidetree = 1;
			}
			else if( njob < 100 || 't' == varpairscore( njob, npickx, nlenmax, seq, randomseed ) )
			{
				if( treein == 'c' ) exit( 1 );
				reporterr( "Tree!\n" );
				treein = 0;
				nguidetree = 2;
			}
			else
			{
				reporterr( "Chain!\n" );
				treein = 's';
				nguidetree = 1;
			}
		}
		else if ( treein != 0 ) // auto no toki arieru
			nguidetree = 1;
	}

# if 0 // tameshini
	if( sueff_global < 0.0001 || compacttree == 2 )
	{
		nthread = 0;
		nthreadtb = 0;
	}
#endif
//	if( njob > 10000 ) nthreadtb = 0;
//	if( njob > 20000 ) nthreadtb = 0;
	/* work/largen: the original runs both progressive passes in the serial treebase() above 20000
	   sequences, whatever --thread says, for the two reasons given below (2018/Jan).  Neither holds
	   in this version:
	   1. treebasethread() builds the member lists from memhist (2018/Mar) and never calls
	      topolorderz()/topolorder_lessargs(), the recursion over the tree whose depth can reach
	      njob; the worker threads also get a stack of max( 8 MB, 200 * njob ) like the main thread
	      (pthread_create below).
	   2. there are no local copies any more: each step aligns its member rows in place (see the
	      note above treebasethread()).
	   With --thread N, nthread stays N on both paths, so A__align never reuses profiles (that needs
	   nthread == 0), and each step computes the same function of its two children in treebase()
	   and in treebasethread().  The one order-dependent piece of state is the memsave switch
	   (alg = 'M' once a child exceeds 30000 columns): above 20000 sequences the threads reproduce
	   where the serial pass switches, exactly (the phases in the driver below; at or below 20000
	   the original's threaded, in-order behaviour is kept, as in work/prog2).  With one thread,
	   treebase() as before. */
	tb_origserial = ( njob > 20000 );
#ifdef GCCPAR_TBTEST
	if( getenv( "MAFFT_LARGEN_ORIGSERIAL" ) ) tb_origserial = 1; // test builds only: behave as above 20000 sequences at any size
#endif
#ifdef PROGTIMING
	if( tb_origserial && nthreadtb <= 1 && !getenv( "MAFFT_LARGEN_TB1" ) ) nthreadtb = 0; // timing builds: MAFFT_LARGEN_TB1 keeps treebasethread() at one thread (per-step logs)
#else
	if( tb_origserial && nthreadtb <= 1 ) nthreadtb = 0;
#endif
#ifdef LARGEN_CUTOFF
	if( njob > LARGEN_CUTOFF ) nthreadtb = 0; // build option: restore a size cutoff
#endif
	/* work/v7: nthreadtbq, the threads for an attempt of a progressive pass whose result cannot
	   depend on the step order: the optimistic ready queue (it gives the original's serial result
	   or abandons the attempt), the phases above 20000 sequences, and a pass that runs in 'M' from
	   its start.  mafft.tmpl caps nthreadtb at 16 for --thread >= 16 (the original's
	   "multithread: 16 made"); it passes the uncapped count as c in -C a-b-c.  The in-order
	   dispatch, whose result with a memsave switch depends on the dispatch (the original's own
	   racy behaviour), keeps nthreadtb threads, i.e. what the original uses. */
	if( nthreadtbq < nthreadtb || nthreadtb <= 1 ) nthreadtbq = nthreadtb;
#ifdef GCCPAR_TBTEST
	if( getenv( "MAFFT_SERIALTB" ) ) nthreadtb = 0; // test builds only: force the serial treebase() path
#endif
// 2018/Jan.  Hairetsu ga ooi toki
// 1. topolorder_lessargs no stack ga tarinakunaru
// 2. localcopy no tame kouritsu warui

	if( compacttree == 1 )
	{
		if( maxdistmtxsize > (double)njob * (njob-1) * sizeof( double ) / 2 ) 
		{
			reporterr( "Use conventional tree.\n" );
			compacttree = 0;
		}
	}

	if( !treein )
	{
		reporterr(       "\n\nMaking a distance matrix ..\n" );
		if( callback && callback( 0, 0, "Distance matrix" ) ) goto chudan;

		double v9pt0 = getenv( "MAFFT_V9PTIME" ) ? sextet_now() : 0.0; // fast-mafft work/v9p: stage timers
		tmpseq = AllocateCharVec( nlenmax+1 );
		grpseq = AllocateIntVec( nlenmax+1 );
		pointt = AllocateIntMtx( njob, nlenmax+1 );
		if( v9pt0 > 0.0 ) reporterr( "\nv9p time: pointt %.3f s\n", sextet_now() - v9pt0 );
		if( !compacttree ) mtx = ( !subalignment && !treeout ) ? v9_halfmtx_alloc( njob ) : AllocateFloatHalfMtx( njob ); // fast-mafft work/v9p: one mapping for the default guide-tree path
		if( v9pt0 > 0.0 ) reporterr( "\nv9p time: + halfmtx %.3f s\n", sextet_now() - v9pt0 );
		if( dorp == 'd' ) tsize = (int)pow( 4, tuplesize );
		else              tsize = (int)pow( 6, 6 );

		if( dorp == 'd' && tuplesize == 6 )
		{
			lenfaca = D6LENFACA;
			lenfacb = D6LENFACB;
			lenfacc = D6LENFACC;
			lenfacd = D6LENFACD;
		}
		else if( dorp == 'd' && tuplesize == 10 )
		{
			lenfaca = D10LENFACA;
			lenfacb = D10LENFACB;
			lenfacc = D10LENFACC;
			lenfacd = D10LENFACD;
		}
		else    
		{
			lenfaca = PLENFACA;
			lenfacb = PLENFACB;
			lenfacc = PLENFACC;
			lenfacd = PLENFACD;
		}

		maxl = 0;
		for( i=0; i<njob; i++ ) 
		{
			gappick0( tmpseq, seq[i] );
			nogaplen[i] = strlen( tmpseq );
			if( nogaplen[i] < 6 )
			{
//				reporterr(       "Seq %d, too short, %d characters\n", i+1, nogaplen[i] );
//				reporterr(       "Please use mafft-ginsi, mafft-linsi or mafft-ginsi\n\n\n" );
//				exit( 1 );
			}
			if( nogaplen[i] > maxl ) maxl = nogaplen[i];
			if( dorp == 'd' ) /* nuc */
			{
				seq_grp_nuc( grpseq, tmpseq );
//				makepointtable_nuc( pointt[i], grpseq );
//				makepointtable_nuc_octet( pointt[i], grpseq );
				if( tuplesize == 10 )
					makepointtable_nuc_dectet( pointt[i], grpseq );
				else if( tuplesize == 6 )
					makepointtable_nuc( pointt[i], grpseq );
				else
				{
					reporterr(       "tuplesize=%d: not supported\n", tuplesize );
					exit( 1 );
				}
			}
			else                 /* amino */
			{
				seq_grp( grpseq, tmpseq );
				makepointtable( pointt[i], grpseq );
			}
		}
		if( nunknown ) reporterr(       "\nThere are %d ambiguous characters.\n", nunknown );
		if( v9pt0 > 0.0 ) reporterr( "\nv9p time: + prep %.3f s\n", sextet_now() - v9pt0 );


		if( compacttree )
		{

			reporterr( "Compact tree, step 1\n" );
			mindistfrom = (int *)calloc( njob, sizeof( int ) );
			mindist = (double *)calloc( njob, sizeof( double ) );
			selfscore = (int *)calloc( njob, sizeof( int ) );
			partmtx = preparepartmtx( njob );


			for( i=0; i<njob; i++ )
			{
				table1 = (int *)calloc( tsize, sizeof( int ) );
				if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
				makecompositiontable_p( table1, pointt[i] );
				selfscore[i] = commonsextet_p( table1, pointt[i] );
				free( table1 );
				table1 = NULL;
			}
			commonsextet_p( NULL, NULL );

#ifdef enablemultithread
			if( nthreadpair > 0 )
			{
				compactdistmtxthread_arg_t *targ;
				int jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
				double **mindistthread;
				int **mindistfromthread;
				if( compacttree == 4 )
					jobpos = 0;
				else
					jobpos = njob-1;
				targ = calloc( nthreadpair, sizeof( compactdistmtxthread_arg_t ) );
				handle = calloc( nthreadpair, sizeof( pthread_t ) );
				mindistthread = AllocateDoubleMtx( nthreadpair, njob );
				mindistfromthread = AllocateIntMtx( nthreadpair, njob );
				pthread_mutex_init( &mutex, NULL );

		
				for( j=0; j<nthreadpair; j++ )
				{
					for( i=0; i<njob; i++ )
					{
						mindistthread[j][i] = 999.9;
						mindistfromthread[j][i] = -1;
					}
					targ[j].thread_no = j;
					targ[j].nogaplen = nogaplen;
					targ[j].pointt = pointt;
					targ[j].selfscore = selfscore;
					targ[j].partmtx = partmtx;
					targ[j].njob = njob;
					targ[j].mindist = mindistthread[j];
					targ[j].mindistfrom = mindistfromthread[j];
					targ[j].jobpospt = &jobpos;
					targ[j].mutex = &mutex;
	
					if( compacttree == 4 )
						pthread_create( handle+j, NULL, ylcompactdisthalfmtxthread, (void *)(targ+j) );
					else
						pthread_create( handle+j, NULL, compactdisthalfmtxthread, (void *)(targ+j) );
				}
		
				for( j=0; j<nthreadpair; j++ ) pthread_join( handle[j], NULL );

				for( i=0; i<njob; i++ )
				{
					mindist[i] = 999.9;
					mindistfrom[i] = -1;
					for( j=0; j<nthreadpair; j++ )
					{
						if( mindistthread[j][i] < mindist[i] )
						{
							mindist[i] = mindistthread[j][i];
							mindistfrom[i] = mindistfromthread[j][i];
						}
					}
				}
				for( i=0; i<njob; i++ ) mindist[i] -= preferenceval( i, mindistfrom[i], njob ); // for debug


				pthread_mutex_destroy( &mutex );
				FreeDoubleMtx( mindistthread );
				FreeIntMtx( mindistfromthread );
				free( handle );
				free( targ );
	
			}
			else
#endif
			{
				compactdistmtxthread_arg_t *targ;
				int jobpos;
		
				if( compacttree == 4 )
					jobpos = 0;
				else
					jobpos = njob-1;
				targ = calloc( 1, sizeof( compactdistmtxthread_arg_t ) );
		
				{
					for( i=0; i<njob; i++ )
					{
						mindist[i] = 999.9;
						mindistfrom[i] = -1;
					}
					targ[0].thread_no = 0;
					targ[0].nogaplen = nogaplen;
					targ[0].pointt = pointt;
					targ[0].selfscore = selfscore;
					targ[0].partmtx = partmtx;
					targ[0].njob = njob;
					targ[0].mindist = mindist;
					targ[0].mindistfrom = mindistfrom;
					targ[0].jobpospt = &jobpos;
	
					if( compacttree == 4 )
						ylcompactdisthalfmtxthread( targ );
					else
						compactdisthalfmtxthread( targ );
				}

				free( targ );

				for( i=0; i<njob; i++ ) mindist[i] -= preferenceval( i, mindistfrom[i], njob ); // for debug
			}

//			reporterr( "\n" );
//			for( i=0; i<njob; i++ ) reporterr( "mindist[%d] = %f, mindistfrom[%d] = %d\n", i, mindist[i], i, mindistfrom[i] );
			reporterr( "\ndone.\n" );

#if 0
			reporterr( "\npartmtx = .\n" );
			for( i=0; i<njob; i++ )
			{
				reporterr( "i=%d\n", i );
				if( partmtx[i] ) for( j=0; j<njob; j++ ) reporterr( "%f ", partmtx[i][j]);
				else reporterr( "nil" );
				reporterr( "\n", i );
			}
#endif
		}
		else
		{
			sextetdone = callback ? 0 : sextetfast_run( njob, pointt, mtx, nogaplen ); // fast-mafft (work/sextet)
			if( !sextetdone ) {
#ifdef enablemultithread
			if( nthreadpair > 0 )
			{
				distancematrixthread_arg_t *targ; 
				int jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
	
				jobpos = 0; 
				targ = calloc( nthreadpair, sizeof( distancematrixthread_arg_t ) ); 
				handle = calloc( nthreadpair, sizeof( pthread_t ) ); 
				pthread_mutex_init( &mutex, NULL );
	
				for( i=0; i<nthreadpair; i++ )
				{
					targ[i].thread_no = i;
					targ[i].njob = njob;
					targ[i].jobpospt = &jobpos;
					targ[i].pointt = pointt;
					targ[i].mtx = mtx;
					targ[i].mutex = &mutex;
	
					pthread_create( handle+i, NULL, distancematrixthread, (void *)(targ+i) );
				}
			
				for( i=0; i<nthreadpair; i++ )
				{
					pthread_join( handle[i], NULL );
				}
				pthread_mutex_destroy( &mutex );
				free( handle );
				free( targ );
			}
			else
#endif
			{
				for( i=0; i<njob; i++ )
				{
					table1 = (int *)calloc( tsize, sizeof( int ) );
					if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
					if( i % 100 == 0 )
					{
						reporterr(       "\r% 5d / %d", i+1, njob );
						if( callback && callback( 0, i*25/njob, "Distance matrix" ) ) goto chudan;
					}
					makecompositiontable_p( table1, pointt[i] );
			
					for( j=i; j<njob; j++ ) 
					{
						mtx[i][j-i] = (double)commonsextet_p( table1, pointt[j] );
					} 
					free( table1 ); table1 = NULL;
				}
			}
			} // fast-mafft (work/sextet)
			if( v9pt0 > 0.0 )
			{
				FILE *smf = fopen( "/proc/self/smaps_rollup", "r" ); // Linux: transparent huge pages in use
				char sml[256];
				reporterr( "\nv9p time: + pairs %.3f s\n", sextet_now() - v9pt0 );
				if( smf )
				{
					while( fgets( sml, sizeof( sml ), smf ) ) if( !strncmp( sml, "Rss:", 4 ) || !strncmp( sml, "AnonHugePages:", 14 ) ) reporterr( "v9p mem: %s", sml );
					fclose( smf );
				}
			}
			reporterr(       "\ndone.\n\n" );
			ien = njob-1;

			// fast-mafft (overhead lever): the original loop, row by row (dist6conv_row), on
			// nthreadpair threads when there are enough rows.
			// fast-mafft (merge): skipped when work/sextet already stored final distances (sextetdone == 2)
			if( sextetdone != 2 )
			{
#ifdef enablemultithread
			if( nthreadpair > 1 && njob >= 1000 )
			{
				dist6convthread_arg_t *ctarg;
				pthread_t *chandle;
				pthread_mutex_t cmutex;
				int cjobpos = 0;
				ctarg = calloc( nthreadpair, sizeof( dist6convthread_arg_t ) );
				chandle = calloc( nthreadpair, sizeof( pthread_t ) );
				pthread_mutex_init( &cmutex, NULL );
				for( i=0; i<nthreadpair; i++ )
				{
					ctarg[i].mtx = mtx;
					ctarg[i].nogaplen = nogaplen;
					ctarg[i].njob = njob;
					ctarg[i].jobpospt = &cjobpos;
					ctarg[i].mutex = &cmutex;
					pthread_create( chandle+i, NULL, dist6convthread, (void *)(ctarg+i) );
				}
				for( i=0; i<nthreadpair; i++ ) pthread_join( chandle[i], NULL );
				pthread_mutex_destroy( &cmutex );
				free( chandle );
				free( ctarg );
			}
			else
#endif
			for( i=0; i<ien; i++ ) dist6conv_row( mtx, nogaplen, njob, i );
			} // fast-mafft (merge): end sextetdone != 2
			if( disopt )
			{
				for( i=0; i<njob; i++ ) 
				{
					sprintf( b, "=lgth = %04d", nogaplen[i] );
					strins( b, name[i] );
				}
			}
			if( getenv( "MAFFT_SEXTET_DUMP" ) ) sextet_dump( getenv( "MAFFT_SEXTET_DUMP" ), njob, mtx ); // fast-mafft (work/sextet)
			FreeIntMtx( pointt ); pointt = NULL;
			commonsextet_p( NULL, NULL );
		}
		free( grpseq ); grpseq = NULL;
		free( tmpseq ); tmpseq = NULL;

#if 0 // writehat2 wo kakinaosu -> iguidetree loop nai ni idou
		if( distout )
		{
			hat2p = fopen( "hat2", "w" );
			WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, mtx );
			fclose( hat2p );
		}
#endif

	}
#if 0 
	else 
	{
		reporterr(       "Loading 'hat2' ... " );
		prep = fopen( "hat2", "r" );
		if( prep == NULL ) ErrorExit( "Make hat2." );
		readhat2_double( prep, njob, name, mtx ); // name chuui
		fclose( prep );
		reporterr(       "done.\n" );
	}
#endif

//	reporterr( "after computing distance matrix," );
//	use_getrusage();

	if( nadd && keeplength )
	{
		originalgaps = (char *)calloc( nlenmax+1, sizeof( char) );
		recordoriginalgaps( originalgaps, njob-nadd, seq );

		if( mapout )
		{
			addbk = (char **)calloc( nadd+1, sizeof( char * ) );
			for( i=0; i<nadd; i++ )
			{
				ien = strlen( seq[njob-nadd+i] );
				addbk[i] = (char *)calloc( ien + 1, sizeof( char ) );
				gappick0( addbk[i], seq[njob-nadd+i] );
			}
			addbk[nadd] = NULL;
		}
		else
			addbk = NULL;
	}
	else
	{
		originalgaps = NULL;
		addbk = NULL;
	}
	


	for( iguidetree=0; iguidetree<nguidetree; iguidetree++ )
//	for( iguidetree=0; ; iguidetree++ )
	{

		alg = algbackup; // tbfast wo disttbfast ni ketsugou shitatame.


		topol = AllocateIntCub( njob, 2, 0 );
		len = AllocateFloatMtx( njob, 2 );

		if( iguidetree == nguidetree - 1 ) calcpairdists = 0;
		else                               calcpairdists = 1;
	
		if( treein )
		{
			nguidetree = 1; //  iranai
			calcpairdists = 0; // iranai
			if( treein == (int)'l' )
			{
				loadtree( njob, topol, len, name, nogaplen, dep, treeout );
			}
			else if( treein == (int)'s' )
			{
				createchain( njob, topol, len, name, nogaplen, dep, treeout, 1, randomseed );
				nthreadpair = nthread = 0;
				weight = 0; // mafft.tmpl kara idou
				tbrweight = 0; // mafft.tmpl kara idou
				nthreadtb = 0; // 2017/Nov
			}
			else if( treein == (int)'p' )
			{
				createchain( njob, topol, len, name, nogaplen, dep, treeout, 0, randomseed );
				nthreadpair = nthread = 0;
				weight = 0; // mafft.tmpl kara idou
				tbrweight = 0; // mafft.tmpl kara idou
				nthreadtb = 0; // 2017/Nov
			}
			else
			{
				reporterr( "Error. treein = %d or %c\n", treein, treein );
				exit( 1 );
			}
		}
		else if( topin )
		{
			reporterr(       "Loading a topology ... " );
			reporterr(       "--topin has been disabled\n" );
			exit( 1 );
//			loadtop( njob, mtx, topol, len );
//			FreeFloatHalfMtx( mtx, njob );
		}
		else
		{
			if( distout )
			{
				hat2p = fopen( "hat2", "w" );
				WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, mtx );
				// writehat2 wo kakinaosu
				fclose( hat2p );

				if( !treeout && noalign )  // 2016Jul31
				{
					writeData_pointer( stdout, njob, name, nlen, seq );
					reporterr(       "\n" );
					SHOWVERSION;
					goto chudan;
//					return( 0 );
				}
			}

			if( subalignment ) // merge ha localmem ni mitaiou
			{
				reporterr(       "Constructing a UPGMA tree ... " );
				fixed_supg_double_realloc_nobk_halfmtx_treeout_constrained( njob, mtx, topol, len, name, nlen, dep, nsubalignments, subtable, !calcpairdists );
				if( !calcpairdists ) 
				{
					FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
				}
			}
			else if( compacttree )
			{
				reporterr(       "Constructing a tree ... " );
// 2017Oct26
				if( compacttree == 4 )
					compacttree_memsaveselectable( njob, partmtx, mindistfrom, mindist, pointt, selfscore, bseq, skiptable, topol, len, name, nogaplen, dep, treeout, 2, 1 ); // howcompact == 2
				else
					compacttreegivendist( njob, mindist, mindistfrom, topol, len, name, dep, treeout );
				if( mindistfrom ) free( mindistfrom ); mindistfrom = NULL;
				if( mindist ) free( mindist );; mindist = NULL;
				if( selfscore ) free( selfscore ); selfscore = NULL;
				if( bseq ) FreeCharMtx( bseq ); bseq = NULL; // nikaime dake
				if( skiptable) FreeIntMtx( skiptable ); skiptable = NULL; // nikaime dake
				if( pointt ) FreeIntMtx( pointt ); pointt = NULL; // ikkaime dake.
				free( partmtx );
			}
			else if( treeout )
			{
				reporterr(       "Constructing a UPGMA tree (treeout, efffree=%d) ... ", !calcpairdists );
				fixed_musclesupg_double_realloc_nobk_halfmtx_treeout_memsave( njob, mtx, topol, len, name, nogaplen, dep, !calcpairdists, treeout );
				if( !calcpairdists )
				{
					FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
				}
			}
			else
			{
				reporterr(       "Constructing a UPGMA tree (efffree=%d) ... ", !calcpairdists );
#ifdef enablemultithread
				v9_upgma( njob, mtx, topol, len, dep, 1, !calcpairdists, nthreadpair ); // fast-mafft: same tree, faster (mltaln9.c; work/v9p)
#else
				v9_upgma( njob, mtx, topol, len, dep, 1, !calcpairdists, 1 ); // fast-mafft: same tree, faster (mltaln9.c; work/v9p)
#endif
				if( !calcpairdists )
				{
					v9_halfmtx_free( mtx, njob ); mtx = NULL;
				}
			}
		}
//		else 
//			ErrorExit( "Unknown tree method\n" );




		if( calcpairdists ) selfscore = AllocateIntVec( njob );


		if( callback && callback( 0, 25, "Guide tree" ) ) goto chudan;
		reporterr(       "\ndone.\n\n" );
		if( callback && callback( 0, 50, "Guide tree" ) ) goto chudan;

		if( sparsepickup && iguidetree == nguidetree-1 )
		{
			reporterr(       "Sparsepickup! \n" );
			pickup( njob, nogaplen, topol, name, seq );
			reporterr(       "done. \n" );
			SHOWVERSION;
			goto chudan;
		}
//		reporterr( "after tree building" );
//		use_getrusage();


		if( treein == 's' || treein == 'p' )
		{
			localmem[0][0] = topol[0][0][0];
			for( i=1; i<njob; i++ )
				localmem[0][i] = topol[i-1][1][0];
		}
		else
		{
			localmem[0][0] = -1;
			posinmem = topolorderz( localmem[0], topol, dep, njob-2, 2 ) - localmem[0];
		}
	
		orderfp = fopen( "order", "w" );
		if( !orderfp )
		{
			reporterr(       "Cannot open 'order'\n" );
			exit( 1 );
		}
#if 0
		for( i=0; (j=topol[njob-2][0][i])!=-1; i++ )
		{
			fprintf( orderfp, "%d\n", j );
		}
		for( i=0; (j=topol[njob-2][1][i])!=-1; i++ )
		{
			fprintf( orderfp, "%d\n", j );
		}
#else
		for( i=0; i<njob; i++ )
			fprintf( orderfp, "%d\n", localmem[0][i] );
#endif
		fclose( orderfp );

	
		if( ( treeout || distout )  && noalign ) 
		{
			writeData_pointer( stdout, njob, name, nlen, seq );
			reporterr(       "\n" );
			SHOWVERSION;
			goto chudan;
//			return( 0 );
		}
		
		if( tbrweight )
		{
			weight = 3; 
#if 0
			utree = 0; counteff( njob, topol, len, eff ); utree = 1;
#else
			counteff_simple_double_nostatic_memsave( njob, topol, len, dep, eff );
//			counteff_simple_double_nostatic( njob, topol, len, eff );
#endif
		}
		else
		{
			for( i=0; i<njob; i++ ) eff[i] = 1.0;
		}
	
#if 0
		for( i=0; i<njob; i++ )
			reporterr(       "eff[%d] = %20.16f\n", i, eff[i] );
		exit( 1 );
#endif
	
	
		FreeFloatMtx( len ); len = NULL;
	
		bseq = AllocateCharMtx( njob, nlenmax*2+1 );
		alloclen = nlenmax*2+1;


	
		if( nadd )
		{
			alignmentlength = strlen( seq[0] );
			for( i=0; i<njob-nadd; i++ )
			{
				if( alignmentlength != strlen( seq[i] ) )
				{
					reporterr(       "#################################################################################\n" );
					reporterr(       "# ERROR!\n" );
					reporterr(       "# The original %d sequences must be aligned\n", njob-nadd );
					reporterr(       "# alignmentlength = %d, but strlen(seq[%d])=%d\n", alignmentlength, i, (int)strlen( seq[i] ) );
					reporterr(       "#################################################################################\n" );
					goto chudan; // TEST!!
					//exit( 1 );
				}
			}
			if( addprofile )
			{
				alignmentlength = strlen( seq[njob-nadd] );
				for( i=njob-nadd; i<njob; i++ )
				{
					if( alignmentlength != strlen( seq[i] ) )
					{
						reporterr(       "###############################################################################\n" );
						reporterr(       "# ERROR!\n" );
						reporterr(       "# The %d additional sequences must be aligned\n", nadd );
						reporterr(       "# Otherwise, try the '--add' option, instead of '--addprofile' option.\n" );
						reporterr(       "###############################################################################\n" );
						exit( 1 );
					}
				}
				for( i=0; i<nadd; i++ ) addmem[i] = njob-nadd+i;
				addmem[nadd] = -1;
				foundthebranch = 0;
				for( i=0; i<njob-1; i++ )
				{
					localmem[0][0] = -1;
					posinmem = topolorderz( localmem[0], topol, dep, i, 0 ) - localmem[0];
					localmem[1][0] = -1;
					posinmem = topolorderz( localmem[1], topol, dep, i, 1 ) - localmem[1];

					if( samemember( localmem[0], addmem ) ) // jissainiha nai
					{
						mergeoralign[i] = '1';
						foundthebranch = 1;
					}
					else if( samemember( localmem[1], addmem ) ) // samemembern ni henkou kanou
					{
						mergeoralign[i] = '2';
						foundthebranch = 1;
					}
					else
					{
						mergeoralign[i] = 'n';
					}
				}
				if( !foundthebranch )
				{
					reporterr(       "###############################################################################\n" );
					reporterr(       "# ERROR!\n" );
					reporterr(       "# There is no appropriate position to add the %d sequences in the guide tree.\n", nadd );
					reporterr(       "# Check whether the %d sequences form a monophyletic cluster.\n", nadd );
					reporterr(       "# If not, try the '--add' option, instead of the '--addprofile' option.\n" );
					reporterr(       "############################################################################### \n" );
					exit( 1 );
				}
				commongappick( nadd, seq+njob-nadd );
				for( i=njob-nadd; i<njob; i++ ) strcpy( bseq[i], seq[i] );
			}
			else
			{
				for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'n';
#if 0
				for( j=njob-nadd; j<njob; j++ )
				{
					addmem[0] = j;
					addmem[1] = -1;
					for( i=0; i<njob-1; i++ )
					{
						reporterr( "Looking for samemember, %d-%d/%d\n", j, i, njob );
						localmem[0][0] = -1;
						posinmem = 0;
						topolorder( localmem[0], &posinmem, topol, dep, i, 0 );
						localmem[1][0] = -1;
						posinmem = 0;
						topolorder( localmem[1], &posinmem, topol, dep, i, 1 );

						if( samemembern( localmem[0], addmem, 1 ) ) // arieru
						{
//							reporterr(       "HIT!\n" );
							if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
							else mergeoralign[i] = '1';
						}
						else if( samemembern( localmem[1], addmem, 1 ) )
						{
//							reporterr(       "HIT!\n" );
							if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
							else mergeoralign[i] = '2';
						}
					}
				}
#else
				for( i=0; i<njob-1; i++ )
				{
//					reporterr( "Looking for samemember, %d-%d/%d\n", j, i, njob );
					localmem[0][0] = -1;
					posinmem = topolorderz( localmem[0], topol, dep, i, 0 ) - localmem[0];
					localmem[1][0] = -1;
					posinmem = topolorderz( localmem[1], topol, dep, i, 1 ) - localmem[1];

					for( j=njob-nadd; j<njob; j++ )
					{
						addmem[0] = j;
						addmem[1] = -1;

						if( samemembern( localmem[0], addmem, 1 ) ) // arieru
						{
//							reporterr(       "HIT!\n" );
							if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
							else mergeoralign[i] = '1';
						}
						else if( samemembern( localmem[1], addmem, 1 ) )
						{
//							reporterr(       "HIT!\n" );
							if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
							else mergeoralign[i] = '2';
						}
					}
				}
#endif
		
				for( i=0; i<nadd; i++ ) addmem[i] = njob-nadd+i;
				addmem[nadd] = -1;
				nlim = njob-1;
//				for( i=0; i<njob-1; i++ )
				for( i=0; i<nlim; i++ )
				{
					localmem[0][0] = -1;
					posinmem = topolorderz( localmem[0], topol, dep, i, 0 ) - localmem[0];
					localmem[1][0] = -1;
					posinmem = topolorderz( localmem[1], topol, dep, i, 1 ) - localmem[1];

					includememberres0 = includemember( localmem[0], addmem );
					includememberres1 = includemember( localmem[1], addmem );
//					if( includemember( topol[i][0], addmem ) && includemember( topol[i][1], addmem ) )
					if( includememberres0 && includememberres1 )
					{
						mergeoralign[i] = 'w';
					}
					else if( includememberres0 )
					{
						mergeoralign[i] = '1';
					}
					else if( includememberres1 )
					{
						mergeoralign[i] = '2';
					}
				}
#if 0
				for( i=0; i<njob-1; i++ )
				{
					reporterr(       "mem0 = " );
					for( j=0; topol[i][0][j]>-1; j++ )	reporterr(       "%d ", topol[i][0][j] );
					reporterr(       "\n" );
					reporterr(       "mem1 = " );
					for( j=0; topol[i][1][j]>-1; j++ )	reporterr(       "%d ", topol[i][1][j] );
					reporterr(       "\n" );
					reporterr(       "i=%d, mergeoralign[] = %c\n", i, mergeoralign[i] );
				}
#endif
				for( i=njob-nadd; i<njob; i++ ) gappick0( bseq[i], seq[i] );
			}
	
//			if( !keeplength ) commongappick( njob-nadd, seq );
			commongappick( njob-nadd, seq );

			for( i=0; i<njob-nadd; i++ ) strcpy( bseq[i], seq[i] );

		}
//--------------- kokokara ----
		else if( subalignment )
		{
			for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'a';
			for( i=0; i<nsubalignments; i++ )
			{
				reporterr(       "Checking subalignment %d:\n", i+1 );
				alignmentlength = strlen( seq[subtable[i][0]] );
//				for( j=0; subtable[i][j]!=-1; j++ )
//					reporterr(       " %d. %-30.30s\n", subtable[i][j]+1, name[subtable[i][j]]+1 );
				for( j=0; subtable[i][j]!=-1; j++ )
				{
					if( subtable[i][j] >= njob ) // check sumi
					{
						reporterr(       "No such sequence, %d.\n", subtable[i][j]+1 );
						exit( 1 );
					}
					if( alignmentlength != strlen( seq[subtable[i][j]] ) )
					{
						reporterr(       "\n" );
						reporterr(       "###############################################################################\n" );
						reporterr(       "# ERROR!\n" );
						reporterr(       "# Subalignment %d must be aligned.\n", i+1 );
						reporterr(       "# Please check the alignment lengths of following sequences.\n" );
						reporterr(       "#\n" );
						reporterr(       "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][0]+1, name[subtable[i][0]]+1, alignmentlength );
						reporterr(       "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][j]+1, name[subtable[i][j]]+1, (int)strlen( seq[subtable[i][j]] ) );
						reporterr(       "#\n" );
						reporterr(       "# See http://mafft.cbrc.jp/alignment/software/merge.html for details.\n" );
						if( subalignmentoffset )
						{
							reporterr(       "#\n" );
							reporterr(       "# You specified seed alignment(s) consisting of %d sequences.\n", subalignmentoffset );
							reporterr(       "# In this case, the rule of numbering is:\n" );
							reporterr(       "#   The aligned seed sequences are numbered as 1 .. %d\n", subalignmentoffset );
							reporterr(       "#   The input sequences to be aligned are numbered as %d .. %d\n", subalignmentoffset+1, subalignmentoffset+njob );
						}
						reporterr(       "###############################################################################\n" );
						reporterr(       "\n" );
						goto chudan; // TEST!!
						//exit( 1 );
					}
					insubtable[subtable[i][j]] = 1;
				}
				for( j=0; j<njob-1; j++ )
				{
					if( includemember( topol[j][0], subtable[i] ) && includemember( topol[j][1], subtable[i] ) )
					{
						mergeoralign[j] = 'n';
					}
				}
				foundthebranch = 0;
				for( j=0; j<njob-1; j++ )
				{
					if( samemember( topol[j][0], subtable[i] ) || samemember( topol[j][1], subtable[i] ) )
					{
						foundthebranch = 1;
						reporterr(       " -> OK\n" );
						break;
					}
				}
				if( !foundthebranch )
				{
					system( "cp infile.tree GuideTree" ); // tekitou
					reporterr(       "\n" );
					reporterr(       "###############################################################################\n" );
					reporterr(       "# ERROR!\n" );
					reporterr(       "# Subalignment %d does not seem to form a monophyletic cluster\n", i+1 );
					reporterr(       "# in the guide tree ('GuideTree' in this directory) internally computed.\n" );
					reporterr(       "# If you really want to use this subalignment, pelase give a tree with --treein \n" );
					reporterr(       "# http://mafft.cbrc.jp/alignment/software/treein.html\n" );
					reporterr(       "# http://mafft.cbrc.jp/alignment/software/merge.html\n" );
					if( subalignmentoffset )
					{
						reporterr(       "#\n" );
						reporterr(       "# You specified seed alignment(s) consisting of %d sequences.\n", subalignmentoffset );
						reporterr(       "# In this case, the rule of numbering is:\n" );
						reporterr(       "#   The aligned seed sequences are numbered as 1 .. %d\n", subalignmentoffset );
						reporterr(       "#   The input sequences to be aligned are numbered as %d .. %d\n", subalignmentoffset+1, subalignmentoffset+njob );
					}
					reporterr(       "############################################################################### \n" );
					reporterr(       "\n" );
					goto chudan; // TEST!!
					//exit( 1 );
				}
//				commongappick( seq[subtable[i]], subalignment[i] ); // irukamo
			}
#if 0
			for( i=0; i<njob-1; i++ )
			{
				reporterr(       "STEP %d\n", i+1 );
				reporterr(       "group1 = " );
				for( j=0; topol[i][0][j] != -1; j++ )
					reporterr(       "%d ", topol[i][0][j]+1 );
				reporterr(       "\n" );
				reporterr(       "group2 = " );
				for( j=0; topol[i][1][j] != -1; j++ )
					reporterr(       "%d ", topol[i][1][j]+1 );
				reporterr(       "\n" );
				reporterr(       "%d -> %c\n\n", i, mergeoralign[i] );
			}
#endif
	
			for( i=0; i<njob; i++ ) 
			{
				if( insubtable[i] ) strcpy( bseq[i], seq[i] );
				else gappick0( bseq[i], seq[i] );
			}
	
			for( i=0; i<nsubalignments; i++ ) 
			{
				for( j=0; subtable[i][j]!=-1; j++ ) subalnpt[i][j] = bseq[subtable[i][j]];
				if( !preservegaps[i] ) commongappick( j, subalnpt[i] );
			}
	
#if 0 // --> iguidetree loop no soto he
			FreeIntMtx( subtable );
			free( insubtable );
			for( i=0; i<nsubalignments; i++ ) free( subalnpt[i] );
			free( subalnpt );
			free( preservegaps );
#endif
		}
//--------------- kokomade ----
		else
		{
			for( i=0; i<njob; i++ ) gappick0( bseq[i], seq[i] );
			for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'a';
		}

		if( calcpairdists ) for( i=0; i<njob; i++ ) selfscore[i] = naivepairscore11( seq[i], seq[i], penalty_dist ); // (int)?
	
		reporterr(       "Progressive alignment %d/%d... \n", iguidetree+1, nguidetree );

//		reporterr( "\nbefore treebase" );
//		use_getrusage();

		cpmxhist = (double ***)calloc( njob-1, sizeof( double ** ) );
		for( i=0; i<njob-1; i++ ) cpmxhist[i] = NULL; 

		memhist = (int **)calloc( njob-1, sizeof( int * ) );
		for( i=0; i<njob-1; i++ ) memhist[i] = NULL; 

#if REPORTCOSTS
		time_t starttime, startclock;
		starttime = time(NULL);
		startclock = clock();
#endif
	
		int tb_runserial = 1; // work/largen: 0 once the threaded pass below has completed
#ifdef enablemultithread
		if( nthreadtb > 0 && nadd == 0 ) // nthreadpair ha minai
		{
			treebasethread_arg_t *targ;
			int jobpos;
			pthread_t *handle;
			pthread_mutex_t mutex;
			pthread_cond_t treecond;
			int *fftlog;
			int nrun;
			int nthread_yoyu;
			int nthread_max, nthread_q; // work/v7: threads of the order-independent attempts (nthreadtbq), allocated for the larger count
			tbsched_t tbsched;
			pthread_rwlock_t aseqlock;
			dphelppool_t tbpool;
			pthread_attr_t tbattr;

			tb_runserial = 0;
			nthread_q = ( nthreadtbq > nthreadtb ) ? nthreadtbq : nthreadtb;
			nthread_max = nthread_q;
			nthread_yoyu = nthreadtb * 1;
			/* work/largen: give the worker threads a stack as large as the one the main thread gets
			   for the tree recursions (setstacksize( 200 * njob ) above), at least 8 MB, instead of the
			   platform default (512 KB on macOS, RLIMIT_STACK or 2 MB with glibc).  The step code does
			   not recurse over the tree (see the note at tb_origserial); a margin, and only address space. */
			pthread_attr_init( &tbattr );
			{
				size_t ss = 0, want = (size_t)200 * njob;
				if( want < ( (size_t)8 << 20 ) ) want = (size_t)8 << 20;
				pthread_attr_getstacksize( &tbattr, &ss );
				if( ss < want ) pthread_attr_setstacksize( &tbattr, want );
			}
			/* work/largen: memory.  Above 20000 sequences the original holds the buffers of one step
			   at a time (plus the distance matrix, during pass 1).  Here every thread keeps its row
			   and traceback arenas between steps, and big traceback blocks go to a shared pool, so
			   several threads' large buffers would stay allocated together.  For these passes the
			   pool is off and a thread frees any arena larger than tb_keepbytes when its step ends:
			   what stays allocated is then only what the running steps use, which are disjoint
			   groups of rows. */
			tb_keepbytes = tb_origserial ? ( (size_t)16 << 20 ) : 0;
#ifdef GCCPAR_TBTEST
			if( getenv( "MAFFT_LARGEN_KEEP" ) ) tb_keepbytes = (size_t)atof( getenv( "MAFFT_LARGEN_KEEP" ) ); // test builds: 0 = keep everything (and pool)
#endif
			A__align_nopool = ( tb_keepbytes > 0 );
			nrun = 0;
			jobpos = 0;
			targ = calloc( nthread_max, sizeof( treebasethread_arg_t ) );
			fftlog = AllocateIntVec( njob );
			handle = calloc( nthread_max, sizeof( pthread_t ) );
			pthread_mutex_init( &mutex, NULL );
			pthread_cond_init( &treecond, NULL );
			pthread_rwlock_init( &aseqlock, NULL );

			/* The only order-dependent state in the progressive pass is the one-way switch to the
			   memsave algorithm (alg = 'M') at the first step (in dispatch order) with a child longer
			   than 30000 columns.  work/prog2:
			   - one thread: always the original in-order dispatch (bitwise the order of the original;
			     the ready queue gains nothing there);
			   - more threads, switch still possible: the ready queue runs "optimistically"; the first
			     step that would flip the switch abandons the pass, and the pass is run again from the
			     start with the original in-order dispatch (tbsched_t.abort, the attempt loop below).
			     Until the switch fires every step uses the same aligner whatever the order, so a
			     pass that is not abandoned is exactly what the original computes; an abandoned one is
			     redone with the original's own dispatch (and its own run-to-run behaviour);
			   - switch impossible (--nomemsave, or already 'M' from pass 1): ready queue. */
			static int tb_abandonedbefore = 0; // an earlier pass had to be redone: expect the same here
			/* work/v10w: a pass that starts with alg == 'M' (--memsave; note that alg is reset to its
			   initial value at the start of every pass, alg = algbackup above, so a switch in pass 1 does
			   not carry over into pass 2, whatever the comment above says) has no order-dependent state:
			   every step uses the memsave aligners whatever the order, so the ready queue gives the
			   in-order pass's result (as with --nomemsave) */
			int tb_inorder = ( nthreadtb <= 1 || ( tb_abandonedbefore && alg != 'M' ) );
			/* work/v10w: when a sequence itself is longer than 30000 columns the switch is certain (at the
			   first step that has it as a child, or earlier): the optimistic attempt can only be abandoned,
			   so the pass starts in MAFFT's in-order dispatch at once, with fresh threads, as MAFFT does
			   (the abandoned attempt made the redo's first claims race differently; see notes) */
			if( !tb_inorder && !nevermemsave && alg != 'M' && !tb_origserial )
			{
				for( i=0; i<njob; i++ ) if( strlen( bseq[i] ) > 30000 ) { tb_inorder = 1; break; }
			}
#ifdef PROGTIMING
			if( nthreadtb <= 1 && getenv( "MAFFT_READYQ1" ) && !tb_abandonedbefore ) tb_inorder = 0; // A/B: optimistic ready queue at one thread
#endif
			int tb_optimistic = ( !nevermemsave && alg != 'M' );
			int tb_attempt;
			/* work/prog3: rows are aligned in place (see treebasethread); MAFFT_COPYROWS=1 in timing
			   builds restores the local copies (A/B) */
#ifdef PROGTIMING
			int tb_nocopy = ( getenv( "MAFFT_COPYROWS" ) == NULL );
#else
			int tb_nocopy = 1;
#endif
			int *tb_rowcap = calloc( njob, sizeof( int ) );
			lzstate_t *tb_lz = NULL;
#ifdef PROGTIMING
			int tb_lazyok = ( getenv( "MAFFT_NOLAZY" ) == NULL ), tb_lazycheck = ( getenv( "MAFFT_LAZYCHECK" ) != NULL );
#else
			int tb_lazyok = 1, tb_lazycheck = 0;
#endif
			for( i=0; i<njob; i++ ) tb_rowcap[i] = alloclen; // AllocateCharMtx( njob, alloclen ): alloclen bytes per row
			char **tb_saveseq = NULL;
			int *tb_savenlen = NULL;
			if( !tb_inorder && tb_optimistic ) // what an abandoned attempt changes, to put it back
			{
				tb_saveseq = calloc( njob, sizeof( char * ) );
				for( i=0; i<njob; i++ ) { tb_saveseq[i] = malloc( strlen( bseq[i] ) + 1 ); strcpy( tb_saveseq[i], bseq[i] ); }
				tb_savenlen = calloc( njob, sizeof( int ) );
				for( i=0; i<njob; i++ ) tb_savenlen[i] = nlen[i];
			}
			/* work/largen: above 20000 sequences the original runs the pass serially, so the memsave
			   switch takes effect exactly at s*, the first step in index order with a child longer
			   than 30000 columns: steps < s* use alg, steps >= s* use 'M'.  Here that result is
			   reproduced with the threads, in phases:
			   0  every step with alg, from the ready queue; a step that finds a child longer than
			      30000 columns is put back and stops every step with a larger or equal index from
			      starting (tbsched_t.limit).  The set of such steps is closed under taking parents
			      (alignments only grow towards the root) and parents have larger indices, so every
			      step < s* still runs (with alg, as in the original) and the smallest step put back
			      is s* itself: at the end, limit = s*.
			   1  only if some step >= s* had already run with alg before limit came down to s*: those
			      results are wrong.  Each topmost such step's subtree is reset (rows and nlen of its
			      members from the start of the pass, done flags, profile caches, work/v7: and its
			      deferred gap tables) and its steps < s* are run again (limit = s*), which gives the
			      same results as in phase 0.
			   2  alg = 'M', then every remaining step (all >= s*), in any order.
			   No switch: phase 0 is the whole pass, as before. */
			int tb_phase = 0, tb_limit = njob-1, tb_sstar = -1;
			for( tb_attempt=0; ; tb_attempt++ )
			{
			nrun = 0;
			jobpos = 0;
			if( tb_phase == 0 )
			{
				for( i=0; i<njob; i++ ) dep[i].done = 0;
				for( i=0; i<njob; i++ ) fftlog[i] = 1;
			}
			tbsched_init( &tbsched, njob, dep, tb_inorder, tb_optimistic, tb_limit, topol, nlen ); // work/v8s: lengths for the priorities
			tbsched.detect = ( tb_origserial && tbsched.optimistic );
			tbsched.heavyfactor = ( tb_keepbytes > 0 ) ? 1.0 : 0.0; // work/largen: see treebasethread()
#ifdef GCCPAR_TBTEST
			if( getenv( "MAFFT_LARGEN_HEAVY" ) ) tbsched.heavyfactor = atof( getenv( "MAFFT_LARGEN_HEAVY" ) ); // test builds: 0 = no budget
			if( getenv( "MAFFT_LARGEN_HEAVYMIN" ) ) tb_heavymin = atof( getenv( "MAFFT_LARGEN_HEAVYMIN" ) ); // test builds: bytes
#endif
			/* work/v7: the order-independent attempts (ready queue) run nthread_q threads, the in-order
			   dispatch nthreadtb, as the original (see nthreadtbq) */
			nthread_yoyu = tbsched.inorder ? nthreadtb : nthread_q;
			/* work/v10w: an in-order attempt also starts nthread_q - nthreadtb helper-only threads (see treebasethread) */
			int nthread_create = ( tbsched.inorder && nthread_yoyu > 1 && nthread_q > nthread_yoyu ) ? nthread_q : nthread_yoyu;
			/* work/prog3: deferred gap insertion needs the rows in place, and a fixed aligner per step:
			   either the optimistic ready queue (a memsave switch abandons the attempt) or one thread
			   (the switch happens in the same thread, before the step looks at alg) */
			/* work/v7: and in the phases above 20000 sequences: phase 0 and 1 are the optimistic ready
			   queue (a step that would switch is put back before it touches anything), phase 2 runs
			   every remaining step with 'M' (nothing is deferred; pending tables are materialized
			   before a step reads its rows); a phase-1 reset drops the tables of the reset subtrees */
			if( tb_nocopy && ( !tbsched.inorder || nthread_yoyu == 1 ) && tb_lazyok )
			{
				if( tb_lz == NULL ) tb_lz = lz_new( njob, bseq, tb_lazycheck );
			}
			else if( tb_lz ) { lz_free( tb_lz, njob ); tb_lz = NULL; }
			/* work/prog2: idle threads help with large alignments (with the in-order dispatch, the
			   threads that wait for the children of their step help as well) */
#ifdef PROGTIMING
			if( nthread_yoyu > 1 && !getenv( "MAFFT_NOHELP" ) ) // A/B switch in timing builds
#else
			if( nthread_yoyu > 1 )
#endif
			{
				tbpool.mutex = &mutex;
				tbpool.cond = &treecond;
				tbpool.nidle = 0;
				tbpool.jobs = NULL;
				tbpool.pending = NULL;
				tbpool.readyprio = -1.0; // work/v8s
				dphelppool = &tbpool;
				if( !tbsched.inorder ) tbsched_readyprio( &tbsched ); // work/v8s: the leaves are already in the heap
			}
#ifdef PROGTIMING
			if( getenv( "MAFFT_PROGLOG" ) ) { ptlog = calloc( njob, sizeof( ptrec_t ) ); pt_base = pt_now(); }
			struct rusage pt_ru0; double pt_w0 = pt_now(); getrusage( RUSAGE_SELF, &pt_ru0 );
#endif
	
			for( i=0; i<nthread_create; i++ )
			{
				targ[i].helperonly = ( i >= nthread_yoyu ); // work/v10w
				targ[i].thread_no = i;
				targ[i].njob = njob;
				targ[i].nrunpt = &nrun;
				targ[i].nlen = nlen;
				targ[i].jobpospt = &jobpos;
				targ[i].topol = topol;
				targ[i].dep = dep;
				targ[i].cpmxhist = cpmxhist;
				targ[i].memhist = memhist;
				targ[i].aseq = bseq;
				targ[i].effarr = eff;
				targ[i].alloclenpt = &alloclen;
				targ[i].fftlog = fftlog;
				targ[i].mergeoralign = mergeoralign;
				targ[i].extanch = extanch;
				targ[i].anchindex = anchindex;
#if 1 // tsuneni SEPARATELYCALCPAIRDISTS
				targ[i].newdistmtx = NULL;
				targ[i].selfscore = NULL;
#else
				if( calcpairdists ) // except for last cycle
				{
					targ[i].newdistmtx = mtx;
					targ[i].selfscore = selfscore;
				}
				else
				{
					targ[i].newdistmtx = NULL;
					targ[i].selfscore = NULL;
				}
#endif
				targ[i].mutex = &mutex;
				targ[i].treecond = &treecond;
				targ[i].sched = &tbsched;
				targ[i].aseqlock = &aseqlock;
				targ[i].rowcap = tb_rowcap;
				targ[i].nocopy = tb_nocopy;
				targ[i].lz = tb_lz;
				targ[i].nrunmax = nthread_yoyu;

				if( pthread_create( handle+i, &tbattr, treebasethread, (void *)(targ+i) ) ) ErrorExit( "Cannot create a thread for the progressive alignment.\n" ); // work/largen: stack size
			}

			for( i=0; i<nthread_create; i++ ) // work/v10w
			{
				pthread_join( handle[i], NULL );
			}
#ifdef PROGTIMING
			{
				struct rusage pt_ru1; double pt_w1 = pt_now(); getrusage( RUSAGE_SELF, &pt_ru1 );
				reporterr( "\nPTRUSAGE pass=%d wall=%.3f user=%.3f sys=%.3f minflt=%ld majflt=%ld nivcsw=%ld nvcsw=%ld\n", iguidetree, pt_w1 - pt_w0,
					( pt_ru1.ru_utime.tv_sec - pt_ru0.ru_utime.tv_sec ) + 1e-6 * ( pt_ru1.ru_utime.tv_usec - pt_ru0.ru_utime.tv_usec ),
					( pt_ru1.ru_stime.tv_sec - pt_ru0.ru_stime.tv_sec ) + 1e-6 * ( pt_ru1.ru_stime.tv_usec - pt_ru0.ru_stime.tv_usec ),
					pt_ru1.ru_minflt - pt_ru0.ru_minflt, pt_ru1.ru_majflt - pt_ru0.ru_majflt, pt_ru1.ru_nivcsw - pt_ru0.ru_nivcsw, pt_ru1.ru_nvcsw - pt_ru0.ru_nvcsw );
			}
			ptlog_write( njob, iguidetree );
			free( ptlog ); ptlog = NULL;
#endif
			dphelppool = NULL;
#ifdef GCCPAR_TBTEST
			if( getenv( "MAFFT_LARGEN_HEAVY" ) ) reporterr( "\nheavy steps put back in this attempt: %d\n", tbsched.ndeferrals );
#endif
			if( !tbsched.abort )
			{
				int goswitch = 0;
				if( tbsched.detect && tbsched.limit < tb_limit ) // work/largen: phase 0 found the switch (see above)
				{
					int nbad = 0, p, l, m;
					char *insub;
					if( tb_phase != 0 ) ErrorExit( "A step below the memsave switch found a child longer than 30000 columns.\n" ); // cannot happen
					tb_sstar = tbsched.limit;
					insub = calloc( njob, sizeof( char ) );
					for( l=njob-2; l>=0; l-- ) // parents have larger indices than their children
					{
						p = tbsched.parent[l];
						if( dep[l].done && l >= tb_sstar && ( p == -1 || ( !insub[p] && !dep[p].done ) ) ) // a topmost step that ran with the wrong algorithm
						{
							if( memhist[l] == NULL ) ErrorExit( "No member list for a step to redo.\n" ); // only the root has none, and it cannot be done here
							for( i=0; (m=memhist[l][i])!=-1; i++ )
							{
								strcpy( bseq[m], tb_saveseq[m] ); nlen[m] = tb_savenlen[m];
#ifdef PROGTIMING
								if( tb_lz && tb_lz->shadow ) // work/v7: the shadow row starts over too
								{
									int n = strlen( bseq[m] ) + 1;
									if( n > tb_lz->shadowcap[m] ) { tb_lz->shadow[m] = realloc( tb_lz->shadow[m], n ); tb_lz->shadowcap[m] = n; }
									strcpy( tb_lz->shadow[m], bseq[m] );
								}
#endif
							}
							free( memhist[l] ); memhist[l] = NULL;
							insub[l] = 1;
							nbad++;
						}
						else if( p != -1 && insub[p] ) insub[l] = 1;
						if( insub[l] )
						{
							dep[l].done = 0; // (the member lists of the steps below it were freed by their parents)
							if( cpmxhist[l] ) { FreeDoubleMtx( cpmxhist[l] ); cpmxhist[l] = NULL; }
							if( tb_lz ) // work/v7: deferred gap tables made inside the subtree are void (its rows start over);
							{           // the topmost step's own table would come from its parent, which has not run
								if( tb_lz->op[l] ) { free( tb_lz->op[l]->gt ); free( tb_lz->op[l] ); tb_lz->op[l] = NULL; }
								tb_lz->stale[l] = 0; tb_lz->len[l] = 0;
							}
						}
					}
					free( insub );
					tbsched_free( &tbsched );
					reporterr( "\nAn alignment exceeds 30000 columns: the memsave mode starts at step %d of %d, as in the original's serial pass", tb_sstar+1, njob-1 );
					if( nbad ) reporterr( "; %d subtree(s) that ran past it are redone", nbad );
					reporterr( ".\n" );
					if( nbad ) // phase 1
					{
						tb_phase = 1;
						tb_limit = tb_sstar;
						continue;
					}
					goswitch = 1;
				}
				else
				{
					tbsched_free( &tbsched );
					if( tb_phase != 1 ) // the pass is complete
					{
						if( tb_lz && tb_lz->stale[njob-2] ) ErrorExit( "lazy gaps: the root was not materialized.\n" );
#ifdef PROGTIMING
						if( tb_lz && tb_lz->shadow ) // every row against its shadow
							for( i=0; i<njob; i++ )
								if( strcmp( tb_lz->shadow[i], bseq[i] ) ) { reporterr( "\nLAZYCHECK MISMATCH at the end of the pass, row %d\n", i ); exit( 1 ); }
#endif
						break;
					}
					goswitch = 1; // phase 1 done: every step < s* is done
				}
				if( goswitch ) // phase 2
				{
					reporterr( "\nSwitching to the memsave mode\n" );
					alg = 'M';
					tb_phase = 2;
					tb_limit = njob-1;
					tb_optimistic = 0;
					continue;
				}
			}
			/* work/prog2: abandoned (see above): undo the attempt and run the pass again in order.
			   The aligned rows and nlen go back to their saved values; profile caches and member
			   lists are dropped (a member list is still owned by its step when that step was handed
			   out but its parent was not); done flags and fftlog are reset at the top. */
			reporterr( "\nAn alignment exceeds 30000 columns: redoing this pass with the original step order.\n" );
			for( i=0; i<njob; i++ ) strcpy( bseq[i], tb_saveseq[i] );
			for( i=0; i<njob; i++ ) nlen[i] = tb_savenlen[i];
			if( tb_lz ) lz_reset( tb_lz, njob, bseq ); // work/prog3: pending tables belong to the abandoned attempt
			for( i=0; i<njob-1; i++ )
			{
				if( cpmxhist[i] ) { FreeDoubleMtx( cpmxhist[i] ); cpmxhist[i] = NULL; }
				if( memhist[i] && tbsched.started[i] && tbsched.parent[i] != -1 && !tbsched.started[tbsched.parent[i]] ) free( memhist[i] );
				memhist[i] = NULL;
			}
			tbsched_free( &tbsched );
			tb_inorder = 1;
			tb_optimistic = 0;
			tb_abandonedbefore = 1;
			}
			if( tb_saveseq )
			{
				for( i=0; i<njob; i++ ) free( tb_saveseq[i] );
				free( tb_saveseq );
				free( tb_savenlen );
			}
			free( tb_rowcap );
			lz_free( tb_lz, njob );
			pthread_mutex_destroy( &mutex );
			pthread_cond_destroy( &treecond );
			pthread_rwlock_destroy( &aseqlock );
			pthread_attr_destroy( &tbattr );
			tb_keepbytes = 0;
			A__align_nopool = 0;
			free( handle );
			free( targ );
			free( fftlog );

		}
#endif
		if( tb_runserial )
		{
#if 0
			if( calcpairdists ) // except for last
			{
				if( treebase( nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, eff, mtx, selfscore, &alloclen, callback ) ) goto chudan;
			}
			else
#endif
			{
//				if( treebase( keeplength && (iguidetree==nguidetree-1), nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, eff, NULL, NULL, deletemap, deletelag, &alloclen, callback ) ) goto chudan;
				if( treebase( nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, memhist, cpmxhist, eff, NULL, NULL, extanch, anchindex, &alloclen, callback ) ) goto chudan;
			}
		}


#if REPORTCOSTS
//		use_getrusage();
		reporterr( "\ntb %d, real = %f min\n", iguidetree, (float)(time(NULL) - starttime)/60.0 );
		reporterr( "tb %d, user = %f min\n", iguidetree, (float)(clock()-startclock)/CLOCKS_PER_SEC/60);
#endif
//		reporterr( "after treebase, " );
//		use_getrusage();
		reporterr(       "\ndone.\n\n" );
		if( callback && callback( 0, 100, "Progressive alignment" ) ) goto chudan;
//		free( topol[njob-1][0] ); topol[njob-1][0]=NULL;
//		free( topol[njob-1][1] ); topol[njob-1][1]=NULL;
//		free( topol[njob-1] ); topol[njob-1]=NULL;
//		free( topol ); topol=NULL;
		FreeIntCub( topol ); topol = NULL;
#if 1 // 2021/Jun/24
		if( cpmxhist ) // nakutemo yoi
		{
			for( i=0; i<njob-1; i++ )
			{
				if( cpmxhist[i] ) 
				{
//					reporterr( "freeing cpmxhist[%d]\n", i );
					FreeDoubleMtx( cpmxhist[i] ); cpmxhist[i] = NULL;
				}
			}
			free( cpmxhist ); cpmxhist = NULL;
		}
#else
		if( cpmxhist[njob-2] ) 
		{
//			reporterr( "freeing cpmxhist[njob-2]\n" );
			FreeDoubleMtx( cpmxhist[njob-2] ); cpmxhist[njob-2] = NULL;
		}
		free( cpmxhist ); cpmxhist = NULL;
#endif

		free( memhist ); memhist = NULL;
//		reporterr( "after freeing topol, " );
//		use_getrusage();


//		reporterr( "compacttree = %d, calcpairdist = %d\n", compacttree, calcpairdists );


//		reporterr( "\nbseq[njob-3] = %s\n", bseq[njob-3] );
//		reporterr( "bseq[njob-2] = %s\n", bseq[njob-2] );
//		reporterr( "bseq[njob-1] = %s\n", bseq[njob-1] );



// Distance matrix from MSA SEPARATELYCALCPAIRDISTS
//		if( iguidetree < nguidetree-1 )
#ifdef enablemultithread
//		if( nthread>0 && nadd==0 ) if( calcpairdists )
		if( calcpairdists && !compacttree )
#else
//		if( 0 && nadd==0 ) if( calcpairdists ) // zettai nai
		if( calcpairdists && !compacttree )
#endif
		{
			reporterr( "Making a distance matrix from msa.. \n" );
			double v9pm0 = getenv( "MAFFT_V9PTIME" ) ? sextet_now() : 0.0; // fast-mafft work/v9p: stage timers
			if( getenv( "MAFFT_MSADIST_SAVEMSA" ) ) msadist_savemsa( getenv( "MAFFT_MSADIST_SAVEMSA" ), njob, bseq, selfscore );
#ifdef enablemultithread
			msadistfast_t *msadistfast = msadistfast_prepare( njob, bseq, penalty_dist, nthreadpair ); // fast-mafft; NULL -> original scorer
#else
			msadistfast_t *msadistfast = msadistfast_prepare( njob, bseq, penalty_dist, 1 ); // fast-mafft; NULL -> original scorer
#endif
			skiptable = AllocateIntMtx( njob, 0 );
			if( msadistfast == NULL || msadistfast->check ) makeskiptable( njob, skiptable, bseq ); // allocate suru. (only the original scorer uses it)
			if( v9pm0 > 0.0 ) reporterr( "\nv9p time: msadist prepare %.3f s\n", sextet_now() - v9pm0 );
#ifdef enablemultithread
			if( nthreadpair > 0 )
			{
				msadistmtxthread_arg_t *targ;
				Jobtable jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
	
				jobpos.i = 0;
				jobpos.j = 0;
	
				targ = calloc( nthreadpair, sizeof( msadistmtxthread_arg_t ) );
				handle = calloc( nthreadpair, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex, NULL );
	
				for( i=0; i<nthreadpair; i++ )
				{
					targ[i].thread_no = i;
					targ[i].njob = njob;
					targ[i].selfscore = selfscore;
					targ[i].iscore = mtx;
					targ[i].seq = bseq;
					targ[i].skiptable = skiptable;
					targ[i].jobpospt = &jobpos;
					targ[i].fast = msadistfast;
					targ[i].mutex = &mutex;
	
					pthread_create( handle+i, NULL, msadistmtxthread, (void *)(targ+i) );
				}
	
				for( i=0; i<nthreadpair; i++ )
				{
					pthread_join( handle[i], NULL );
				}
				pthread_mutex_destroy( &mutex );
				free( handle );
				free( targ );
			}
			else
#endif
			{
//				reporterr( "Check source!\n" );
//				exit( 1 );

#if 1
				msadistmtxthread_arg_t *targ;
				Jobtable jobpos;

				jobpos.i = 0;
				jobpos.j = 0;
	
				targ = calloc( 1, sizeof( msadistmtxthread_arg_t ) );
	
				{
					targ[0].thread_no = 0;
					targ[0].njob = njob;
					targ[0].selfscore = selfscore;
					targ[0].iscore = mtx;
					targ[0].seq = bseq;
					targ[0].skiptable = skiptable;
					targ[0].jobpospt = &jobpos;
					targ[0].fast = msadistfast;
	
					msadistmtxthread( targ );
				}
	
				free( targ );
#endif
			}
			if( v9pm0 > 0.0 ) reporterr( "\nv9p time: msadist + pairs %.3f s\n", sextet_now() - v9pm0 );
			if( msadistfast && msadistfast->timing ) reporterr( "\nmsadist time (thread-s): residue rows %.3f, gap-run rows %.3f, pair loops %.3f\n", msadistfast->tcol, msadistfast->tgap, msadistfast->tpair );
			msadistfast_free( msadistfast ); msadistfast = NULL;
			if( getenv( "MAFFT_MSADIST_DUMP" ) ) msadist_dump( getenv( "MAFFT_MSADIST_DUMP" ), njob, mtx );
			if( skiptable) FreeIntMtx( skiptable ); skiptable = NULL;
			reporterr(       "\ndone.\n\n" );
			free( selfscore ); selfscore = NULL;
			FreeCharMtx( bseq ); bseq = NULL;
		}
		else if( calcpairdists && compacttree )
		{
			reporterr( "Making a compact tree from msa, step 1.. \n" );
			skiptable = AllocateIntMtx( njob, 0 );
			makeskiptable( njob, skiptable, bseq ); // allocate suru.
			mindistfrom = (int *)calloc( njob, sizeof( int ) );
			mindist = (double *)calloc( njob, sizeof( double ) );
			partmtx = preparepartmtx( njob );
#ifdef enablemultithread
			if( nthreadpair > 0 )
			{
				msacompactdistmtxthread_arg_t *targ;
				int jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
				double **mindistthread;
				int **mindistfromthread;

				mindistthread = AllocateDoubleMtx( nthreadpair, njob );
				mindistfromthread = AllocateIntMtx( nthreadpair, njob );
				targ = calloc( nthreadpair, sizeof( msacompactdistmtxthread_arg_t ) );
				handle = calloc( nthreadpair, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex, NULL );
				if( compacttree == 4 )
					jobpos = 0;
				else
					jobpos = njob-1;

				for( i=0; i<nthreadpair; i++ )
				{
					for( j=0; j<njob; j++ )
					{
						mindistthread[i][j] = 999.9;
						mindistfromthread[i][j] = -1;
					}
					targ[i].thread_no = i;
					targ[i].njob = njob;
					targ[i].selfscore = selfscore;
					targ[i].partmtx = partmtx;
					targ[i].seq = bseq;
					targ[i].skiptable = skiptable;
					targ[i].jobpospt = &jobpos;
					targ[i].mindistfrom = mindistfromthread[i];
					targ[i].mindist = mindistthread[i];
					targ[i].mutex = &mutex;
	
					if( compacttree == 4 )
						pthread_create( handle+i, NULL, ylmsacompactdisthalfmtxthread, (void *)(targ+i) );
					else
						pthread_create( handle+i, NULL, msacompactdisthalfmtxthread, (void *)(targ+i) );
				}
	
				for( i=0; i<nthreadpair; i++ ) pthread_join( handle[i], NULL );
				pthread_mutex_destroy( &mutex );

				for( i=0; i<njob; i++ )
				{
					mindist[i] = 999.9;
					mindistfrom[i] = -1;
					for( j=0; j<nthreadpair; j++ )
					{
						if( mindistthread[j][i] < mindist[i] )
						{
							mindist[i] = mindistthread[j][i];
							mindistfrom[i] = mindistfromthread[j][i];
						}
					}
				}
				for( i=0; i<njob; i++ ) mindist[i] -= preferenceval( i, mindistfrom[i], njob ); // for debug

				free( handle );
				free( targ );
				FreeDoubleMtx( mindistthread );
				FreeIntMtx( mindistfromthread );
			}
			else
#endif
			{
				msacompactdistmtxthread_arg_t *targ;
				int jobpos;
				if( compacttree == 4 )
					jobpos = 0;
				else
					jobpos = njob-1;
				targ = calloc( 1, sizeof( msacompactdistmtxthread_arg_t ) );

				{
					for( j=0; j<njob; j++ )
					{
						mindist[j] = 999.9;
						mindistfrom[j] = -1;
					}
					targ[0].thread_no = 0;
					targ[0].njob = njob;
					targ[0].selfscore = selfscore;
					targ[0].partmtx = partmtx;
					targ[0].seq = bseq;
					targ[0].skiptable = skiptable;
					targ[0].jobpospt = &jobpos;
					targ[0].mindistfrom = mindistfrom;
					targ[0].mindist = mindist;
	
					if( compacttree == 4 )
						ylmsacompactdisthalfmtxthread( targ );
					else
						msacompactdisthalfmtxthread( targ );
//					msacompactdistmtxthread( targ );
				}
				free( targ );
				for( i=0; i<njob; i++ ) mindist[i] -= preferenceval( i, mindistfrom[i], njob ); // for debug
			}
			reporterr( "\n" );
//			free( selfscore ); selfscore = NULL; // mada tsukau
//			FreeCharMtx( bseq ); bseq = NULL; // mada tsukau
//			if( skiptable) FreeIntMtx( skiptable ); skiptable = NULL;

//			reporterr( "\n" );
//			for( i=0; i<njob; i++ ) reporterr( "mindist[%d] = %f\n", i, mindist[i] );
//			exit( 1 );
		}
// Distance matrix from MSA end
//		reporterr( "at the end of guidetree loop, " );
//		use_getrusage();

	}


	if( oneiteration )
	{
		reporterr( "Iterative refinement (one vs others)\n" );
		dooneiteration( nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, memhist, cpmxhist, eff, NULL, NULL, extanch, anchindex, &alloclen, callback );
	}

#if DEBUG
	reporterr(       "closing trap_g\n" );
#endif
//	fclose( trap_g );
//	reporterr( "after guidetree loop, " );
//	use_getrusage();

	if( keeplength )
	{

		dlf = fopen( "_deletelist", "w" );
		deletelist = (GapPos **)calloc( nadd+1, sizeof( GapPos * ) );
		for( i=0; i<nadd; i++ )
		{
			deletelist[i] = calloc( 1, sizeof( GapPos ) );
			deletelist[i][0].pos = -1;
			deletelist[i][0].len = 0;
		}
		deletelist[nadd] = NULL;
		ndeleted = deletenewinsertions_whole( njob-nadd, nadd, bseq, bseq+njob-nadd, deletelist );

		for( i=0; i<nadd; i++ )
		{
			if( deletelist[i] )
				for( j=0; deletelist[i][j].pos!=-1; j++ )
//					fprintf( dlf, "%d %d\n", njob-nadd+i, deletelist[i][j] ); // 0origin
					fprintf( dlf, "%d %d %d\n", njob-nadd+i, deletelist[i][j].pos, deletelist[i][j].len ); // 0origin
		}
		fclose( dlf );

		restoreoriginalgaps( njob, bseq, originalgaps );
		free( originalgaps ); originalgaps = NULL; // 2017/Nov/15

		if( mapout )
		{
			dlf = fopen( "_deletemap", "w" );
			if( mapout == 1 )
				reconstructdeletemap( nadd, addbk, deletelist, bseq+njob-nadd, dlf, name+njob-nadd );
			else
				reconstructdeletemap_compact( nadd, addbk, deletelist, seq+njob-nadd, dlf, name+njob-nadd );
			FreeCharMtx( addbk );
			addbk = NULL;
			fclose( dlf );
		}

//		FreeIntMtx( deletelist );
//		deletelist = NULL;
		for( i=0; deletelist[i] != NULL; i++ ) free( deletelist[i] );
		free( deletelist );
		deletelist = NULL;
	}

	if( scoreout )
	{
		unweightedspscore = plainscore( njob, bseq );
		reporterr(       "\nSCORE %s = %.0f, ", "(treebase)", unweightedspscore );
		reporterr(       "SCORE / residue = %f", unweightedspscore / ( njob * strlen( bseq[0] ) ) );
		reporterr(       "\n\n" );
	}

#if DEBUG
	reporterr(       "writing alignment to stdout\n" );
#endif


	val = 0;
	if( ngui ) 
	{
		ien = strlen( bseq[0] );
		if( ien > lgui )
		{
			reporterr( "alignmentlength = %d, gui allocated %d", ien, lgui );
			val = GUI_LENGTHOVER;
		}
		else
		{
			for( i=0; i<njob; i++ ) 
			{
#if 1
				strcpy( seqgui[i], bseq[i] );
#else
				free( seqgui[i] );
				seqgui[i] =  bseq[i];
#endif
			}
		}
	}
	else
	{
		writeData_pointer( stdout, njob, name, nlen, bseq );
	} 

	if( spscoreout ) reporterr( "Unweighted sum-of-pairs score = %10.5f\n", sumofpairsscore( njob, bseq ) );
	nthread = nthreadtb; // toriaezu
	SHOWVERSION;
	if( ndeleted > 0 )
	{
		reporterr( "\nTo keep the alignment length, %d letters were DELETED.\n", ndeleted );
		if( mapout )
			reporterr( "The deleted letters are shown in the (filename).map file.\n" );
		else
			reporterr( "To know the positions of deleted letters, rerun the same command with the --mapout option.\n" );
	}



	if( subalignment )
	{
		FreeIntMtx( subtable );
		free( insubtable );
		for( i=0; i<nsubalignments; i++ ) free( subalnpt[i] );
		free( subalnpt );
		free( preservegaps );
	}
	if( useexternalanchors )
	{
		free( extanch );
		FreeIntMtx( anchindex );
	}


#if 1 // seqgui[i] =  bseq[i] no toki bseq ha free shinai
	FreeCharMtx( bseq );
#endif
	FreeCharMtx( name );
    free( nlen );

	free( mergeoralign );
	FreeCharMtx( seq );
    free( nogaplen );

	free( mseq1 );
	free( mseq2 );
//	FreeIntCub( topol ); // 
//	FreeFloatMtx( len ); //
//	free( mergeoralign ); //
	free( dep );

	if( nadd ) free( addmem );
	FreeIntMtx( localmem );
	free( eff );
	freeconstants();
	closeFiles();
	FreeCommonIP();
	if( originalgaps ) free( originalgaps ); originalgaps = NULL;
	if( deletelist ) 
	{
		for( i=0; deletelist[i] != NULL; i++ ) free( deletelist[i] );
		free( deletelist );
		deletelist = NULL;
	}

//	use_getrusage();

	return( val );

chudan:

	if( nlen ) free( nlen ); nlen = NULL;
	if( seq ) FreeCharMtx( seq ); seq = NULL;
	if( mseq1 ) free( mseq1 ); mseq1 = NULL;
	if( mseq2 ) free( mseq2 ); mseq2 = NULL;
//	if( topol ) 
//	{
//		for( i=0; i<njob; i++ )
//		{
//			if( topol[i] && topol[i][0] ) 
//			{
//				free( topol[i][0] ); topol[i][0] = NULL;
//			}
//			if( topol[i] && topol[i][1] ) 
//			{
//				free( topol[i][1] ); topol[i][1] = NULL;
//			}
//			if( topol[i] ) free( topol[i] ); topol[i] = NULL;
//		}
//		free( topol ); topol = NULL;
//	}
	if( topol ) FreeIntCub( topol ); topol = NULL;
	if( len ) FreeFloatMtx( len ); len = NULL;
	if( eff ) free( eff ); eff = NULL;
	if( mergeoralign ) free( mergeoralign ); mergeoralign = NULL;
	if( dep ) free( dep ); dep = NULL;
	if( addmem ) free( addmem ); addmem = NULL;
	if( localmem ) FreeIntMtx( localmem ); localmem = NULL;
	if( name ) FreeCharMtx( name ); name = NULL;
	if( nogaplen ) free( nogaplen ); nogaplen = NULL;

	if( tmpseq ) free( tmpseq ); tmpseq = NULL;
	if( grpseq ) free( grpseq ); grpseq = NULL;
	if( pointt ) FreeIntMtx( pointt ); pointt = NULL;
	if( mtx ) v9_halfmtx_free( mtx, njob ); mtx = NULL; // fast-mafft work/v9p
	if( table1 ) free( table1 ); table1 = NULL;

	if( bseq ) FreeCharMtx( bseq ); bseq = NULL;
	if( selfscore ) free( selfscore ); selfscore = NULL;
	if( skiptable ) FreeIntMtx( skiptable ); skiptable = NULL;
	if( originalgaps ) free( originalgaps ); originalgaps = NULL;
	if( deletelist ) 
	{
		for( i=0; deletelist[i] != NULL; i++ ) free( deletelist[i] );
		free( deletelist );
		deletelist = NULL;
	}


	if( subtable ) FreeIntMtx( subtable ); subtable = NULL;
	if( insubtable ) free( insubtable ); insubtable = NULL;
	for( i=0; i<nsubalignments; i++ ) 
	{
		if( subalnpt[i] ) free( subalnpt[i] ); subalnpt[i] = NULL;
	}
	if( subalnpt ) free( subalnpt ); subalnpt = NULL;
	if( preservegaps ) free( preservegaps ); preservegaps = NULL;


	if( mindistfrom ) free( mindistfrom ); mindistfrom = NULL;
	if( mindist ) free( mindist ); mindist = NULL;

	if( cpmxhist )
	{
		for( i=0; i<njob-1; i++ )
		{
			if( cpmxhist[i] ) FreeDoubleMtx( cpmxhist[i] ); cpmxhist[i] = NULL;
		}
		free( cpmxhist ); cpmxhist = NULL;
	}

	if( memhist )
	{
		for( i=0; i<njob-1; i++ )
		{
			if( memhist[i] ) free( memhist[i] ); memhist[i] = NULL;
		}
		free( memhist ); memhist = NULL;
	}

	freeconstants();
	closeFiles();
	FreeCommonIP();

	return( GUI_CANCEL );
}

int main( int argc, char **argv )
{
	int res;
#if !defined(mingw) && !defined(_MSC_VER)
	/* fast-mafft (overhead lever): read the input and write the alignment through 1 MB stdio
	   buffers instead of one read()/write() per 4 KB block; the bytes are the same. */
	if( !isatty( fileno( stdin ) ) ) setvbuf( stdin, NULL, _IOFBF, 1<<20 );
	if( !isatty( fileno( stdout ) ) ) setvbuf( stdout, NULL, _IOFBF, 1<<20 );
#endif
	res = disttbfast( 0, 0, NULL, NULL, argc, argv, NULL );
	if( res == GUI_CANCEL ) res = 0; // treeout de goto chudan wo riyousuru
	return res;
}
