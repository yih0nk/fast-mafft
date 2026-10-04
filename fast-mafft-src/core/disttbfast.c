#ifndef _GNU_SOURCE
#define _GNU_SOURCE // fast-mafft: glibc hides pthread_rwlock_t under -std=c99 without a feature macro
#endif
#include "mltaln.h"
#include "progtiming.h"
#include "dphelp.h"
#ifdef PROGTIMING
typedef struct { int l, thread, clus1, clus2, len1, len2, c0, c1; double t[8]; double acc[PT_NACC]; } ptrec_t;
static ptrec_t *ptlog = NULL;
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
						nthreadtb = myatoi( strchr( *argv, '-' )+1 );
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
 * all counts are popcounts (NEON when available).
 * Generic path (any other residue bytes < 0x80): the same gap masks plus a
 * residue mask; the both-residue columns are visited bit by bit and scored
 * from a flat copy of amino_dis.
 * Rows of unequal length or bytes >= 0x80 fall back to the original code.
 *
 * Debug (environment):
 *   MAFFT_MSADIST=ref      original function for every pair
 *   MAFFT_MSADIST=check    fast and original for every pair; exit(1) on any
 *                          bitwise difference ("generic,check" also works)
 *   MAFFT_MSADIST=generic  skip the bit-plane path
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

typedef struct _msadistfast
{
	int nseq;
	int L;        // aligned length
	int W;        // 64-bit words per row mask; even, and W*64 > L so the last carry has room
	int penal;
	int check;    // 1: also run naivepairscorefast() and compare bitwise
	int xorpath;  // 1: bit-plane path, 0: generic path
	char **seq;   // the aligned rows (generic path reads residues from them)
	/* bit-plane path: <= 4 residue kinds and amino_dis[x][y] depends only on code(x)^code(y) */
	uint64_t *planes; // per row, per 128 columns: G[2] b0[2] b1[2] (gap mask, code bit 0, code bit 1)
	long long Td[4];  // amino_dis value for code xor d
	/* generic path: any residue bytes < 0x80 */
	uint64_t *masks;  // per row, per 64 columns: G, R (gap mask, residue mask)
	int *T;           // 0x80*0x80 copy of amino_dis (the gap row and column are never used)
	/* preparation */
	int nthreads;
	int phase;
	unsigned char lut[0x80];
} msadistfast_t;

typedef struct _msadistprep
{
	msadistfast_t *f;
	int r0, r1;
	int bad;               // phase 0: a row of another length, or a byte >= 0x80
	unsigned char seen[0x80];
	int nknown;            // phase 0: bytes seen so far, in known[] (NEON fast scan)
	unsigned char known[0x80];
} msadistprep_t;

static void msadistfast_free( msadistfast_t *f )
{
	if( f == NULL ) return;
	free( f->planes );
	free( f->masks );
	free( f->T );
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

// phase 0: check lengths and collect the residue bytes; phase 1: build the row bit masks
static void *msadist_prepthread( void *arg )
{
	msadistprep_t *p = (msadistprep_t *)arg;
	msadistfast_t *f = p->f;
	int i, c, k, w, cmax, L = f->L, W = f->W;
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
#endif
			for( ; c<L; c++ ) { hi |= s[c]; p->seen[s[c]&0x7f] = 1; }
			if( hi & 0x80 ) p->bad = 1;
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
					mk[2*w] = g; mk[2*w+1] = ~g;
					continue;
				}
#endif
				for( k=0; k<cmax; k++ )
				{
					v = ( s[w*64+k] == '-' );
					g |= v << k;
					r |= ( v ^ 1 ) << k;
				}
				mk[2*w] = g; mk[2*w+1] = r;
			}
		}
	}
	return( NULL );
}

static void msadist_prepphase( msadistfast_t *f, msadistprep_t *p, int phase )
{
	int t, nth = f->nthreads;
	f->phase = phase;
	for( t=0; t<nth; t++ )
	{
		p[t].f = f;
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

// Returns NULL when the original function must be used.
// MAFFT_MSADIST may contain "ref" (use the original), "check" (compare every pair
// with the original) and "generic" (skip the bit-plane path).
static msadistfast_t *msadistfast_prepare( int nseq, char **seq, int penal, int nthreads )
{
	msadistfast_t *f;
	msadistprep_t *p;
	int code[4];
	unsigned char chr[4];
	int t, a, b, c, K, bad;
	long long Td[4];
	unsigned char seen[0x80];
	const char *mode = getenv( "MAFFT_MSADIST" );

	if( mode && strstr( mode, "ref" ) ) return( NULL );
	if( nseq < 2 ) return( NULL );
	if( nthreads < 1 ) nthreads = 1;
	if( nthreads > nseq ) nthreads = nseq;

	f = calloc( 1, sizeof( msadistfast_t ) );
	p = calloc( nthreads, sizeof( msadistprep_t ) );
	f->nseq = nseq;
	f->seq = seq;
	f->L = strlen( seq[0] );
	f->W = ( f->L / 64 + 2 ) & ~1;
	f->penal = penal;
	f->check = ( mode && strstr( mode, "check" ) );
	f->nthreads = nthreads;

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
		free( p );
		msadistfast_free( f );
		return( NULL );
	}

	K = 0;
	for( c=1; c<0x80; c++ ) if( seen[c] && c != '-' )
	{
		if( K < 4 ) chr[K] = c;
		K++;
	}
	f->xorpath = ( K <= 4 && !( mode && strstr( mode, "generic" ) ) && msadist_findxorcode( K, chr, code, Td ) );
	if( f->xorpath )
	{
		for( a=0; a<4; a++ ) f->Td[a] = Td[a];
		memset( f->lut, 0, sizeof( f->lut ) );
		f->lut['-'] = 4;
		for( a=0; a<K; a++ ) f->lut[chr[a]] = code[a];
		f->planes = calloc( (size_t)nseq * 3 * f->W, sizeof( uint64_t ) );
		if( f->check ) reporterr( "msadist fast (bit-plane path): L=%d K=%d Td=%lld,%lld,%lld,%lld penal=%d\n", f->L, K, Td[0], Td[1], Td[2], Td[3], penal );
	}
	else
	{
		f->T = calloc( 0x80 * 0x80, sizeof( int ) );
		for( a=0; a<0x80; a++ ) for( b=0; b<0x80; b++ ) f->T[(a<<7)|b] = ( a == '-' || b == '-' ) ? 0 : amino_dis[a][b];
		f->masks = calloc( (size_t)nseq * 2 * f->W, sizeof( uint64_t ) );
		if( f->check ) reporterr( "msadist fast (generic path): L=%d K=%d penal=%d\n", f->L, K, penal );
	}
	msadist_prepphase( f, p, 1 );
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

static inline long long msadist_xorscore( const msadistfast_t *f, const uint64_t *pi, const uint64_t *pj )
{
	long long len = f->L, nA = 0, nB = 0, n1 = 0, n2 = 0, nE = 0;
	int nblk = f->W / 2, k0, k1, k;
	uint64_t ci = 0, cj = 0, e0, e1;
#if MSADIST_NEON
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
#else
	(void)k0; (void)k1;
	for( k=0; k<nblk; k++ )
	{
		const uint64_t *qi = pi + 6*k, *qj = pj + 6*k;
		int h;
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
#endif
	// nA counts columns with a gap in either row, nB those that are not identical
	// residue pairs.  Padding columns (>= len) have no gap bit and code 0 in both
	// rows, so they add to neither.  Both-residue columns: len - nA; with code
	// xor 0: len - nB; xor 1: n1; xor 2: n2; xor 3: the rest.
	return( f->Td[3] * ( len - nA ) + ( f->Td[0] - f->Td[3] ) * ( len - nB ) + ( f->Td[1] - f->Td[3] ) * n1 + ( f->Td[2] - f->Td[3] ) * n2
	      + (long long)f->penal * nE );
}

static inline long long msadist_genericscore( const msadistfast_t *f, int i, int j )
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
		rr = mi[2*w+1] & mj[2*w+1];
		while( rr )
		{
			int k = __builtin_ctzll( rr );
			s += T[(a[k]<<7)|b[k]];
			rr &= rr - 1;
		}
		MSADIST_RUNWORD( mi[2*w], mj[2*w], ci, cj, e );
		nE += __builtin_popcountll( e );
	}
	return( s + (long long)f->penal * nE );
}

// exact value of naivepairscorefast( seq[i], seq[j], skip[i], skip[j], penal ) for this alignment
static inline double msadistfast_score( const msadistfast_t *f, int i, int j )
{
	long long s;
	if( f->xorpath )
		s = msadist_xorscore( f, f->planes + (size_t)i * 3 * f->W, f->planes + (size_t)j * 3 * f->W );
	else
		s = msadist_genericscore( f, i, j );
	return( (double)s );
}
/* END msadistfast kernel */

static double msadist_pairscore( msadistmtxthread_arg_t *targ, int i, int j )
{
	msadistfast_t *f = (msadistfast_t *)targ->fast;
	double v;
	if( f == NULL ) return( naivepairscorefast( targ->seq[i], targ->seq[j], targ->skiptable[i], targ->skiptable[j], penalty_dist ) );
	v = msadistfast_score( f, i, j );
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

	while( 1 )
	{
#ifdef enablemultithread
		if( nthreadpair ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = jobpospt->i; // (jobpospt-i)++ dato, shuuryou hantei no mae ni ++ surunode, tomaranakunaru.

			if( i == nlim )
			{
				pthread_mutex_unlock( targ->mutex );
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
			if( i == nlim ) return( NULL );
			if( i % 100 == 0 ) fprintf( stderr, "\r% 5d / %d", i, njob );
		}

		ssi = selfscore[i];
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
				iscoretmp = ( 1.0 - msadist_pairscore( targ, i, j ) / bunbo ) * 2.0; // 2014/Aug/15 fast; fast-mafft: exact replacement for naivepairscorefast()
				if( iscoretmp > 10 ) iscoretmp = 10.0;  // 2015/Mar/17
	
			}
			if( iscoretmp < 0.0 ) 
			{
				reporterr( "WARNING: negative distance, iscoretmp = %f\n", iscoretmp );
				iscoretmp = 0.0;
			}
			iscore[i][j-i] = iscoretmp;
//			printf( "i,j=%d,%d, iscoretmp=%f\n", i, j, iscoretmp );

		}
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
 * Amino acids (tsize 46656): each row keeps its list of distinct tuples with
 * their counts (u16).  For row i a dense u16 table of c_i is filled once and
 * each j sums min( table[p], c_j(p) ) over j's list.  If any count exceeds
 * 65535 the original code is used.  Other tuple sizes use the original code.
 *
 * Debug (environment):
 *   MAFFT_SEXTET=ref     original loop
 *   MAFFT_SEXTET=check   also run commonsextet_p() for every pair and exit(1)
 *                        on any bitwise difference of the stored double
 *   MAFFT_SEXTET=scalar  DNA: skip the NEON block path (plain loop + overflow lists)
 *   MAFFT_SEXTET=list    DNA: use the list path
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
	int kind;        // 1: dense DNA tables, 2: tuple lists
	int check, scalar, timing;
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
#if SEXTET_NEON
	uint32_t d[16];
#endif

	for( jb=i0; jb<n; jb+=4 )
	{
		mcj = MAX( MAX( s->mc[jb], s->mc[jb+1] ), MAX( s->mc[jb+2], s->mc[jb+3] ) );
		for( ib=i0; ib<i1 && ib<=jb; ib+=4 )
		{
			mcb = MAX( mcj, MAX( MAX( s->mc[ib], s->mc[ib+1] ), MAX( s->mc[ib+2], s->mc[ib+3] ) ) );
#if SEXTET_NEON
			if( mcb <= 127 && !s->scalar )
			{
				sextet_absdiff4x4( s->tab + (size_t)ib * SEXTET_DW, s->tab + (size_t)jb * SEXTET_DW, mcb, d );
				for( r=0; r<4; r++ ) for( c=0; c<4; c++ )
				{
					i = ib + r; j = jb + c;
					if( i < n && j < n && j >= i )
						sextet_store( s, i, j, ( s->tot[i] + s->tot[j] - (int)d[r*4+c] ) / 2 );
				}
				continue;
			}
#endif
			for( r=0; r<4; r++ ) for( c=0; c<4; c++ )
			{
				i = ib + r; j = jb + c;
				if( i < n && j < n && j >= i )
					sextet_store( s, i, j, sextet_pairslow( s, i, j ) );
			}
		}
	}
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

static void *sextet_thread( void *arg )
{
	sextetfast_t *s = (sextetfast_t *)arg;
	int i, i0, i1, step;
	int *cnt;
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
	while( ( i0 = sextet_claim( s, step ) ) < s->n )
	{
		i1 = MIN( i0 + step, s->n );
		if( i0 % 100 == 0 || i0 / 100 != ( i1 - 1 ) / 100 ) reporterr( "\r% 5d / %d", i1, s->n );
		if( s->kind == 1 ) sextet_densejob( s, i0 );
		else for( i=i0; i<i1; i++ ) sextet_listrow( s, i, tab );
		if( s->check ) for( i=i0; i<i1; i++ ) sextet_checkrow( s, i );
	}
	if( s->check ) commonsextet_p( NULL, NULL );
	free( tab );
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
			t1 = sextet_now();
			if( s.bad ) break;
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
	if( s.timing ) reporterr( "\nsextetfast: kind %d, %d threads, tables %.3f s, pairs %.3f s\n", s.kind, nthr, t1 - t0, t2 - t1 );
	if( s.check ) reporterr( "\nsextetfast check passed: kind %d, %d rows, every pair's %s bitwise identical\n", s.kind, njob, s.convert ? "distance" : "count" );
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
} tbsched_t;

static int tbsched_before( tbsched_t *sc, int a, int b )
{
	if( sc->prio[a] != sc->prio[b] ) return( sc->prio[a] > sc->prio[b] );
	return( a < b );
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
	return( top );
}

static int tbsched_isready( Treedep *dep, int l )
{
	return( ( dep[l].child0 == -1 || dep[dep[l].child0].done ) && ( dep[l].child1 == -1 || dep[dep[l].child1].done ) );
}

static void tbsched_init( tbsched_t *sc, int njob, Treedep *dep, int inorder, int optimistic )
{
	int l, nstep = njob - 1;
	double *size = calloc( nstep, sizeof( double ) );
	sc->heap = calloc( nstep, sizeof( int ) );
	sc->parent = calloc( nstep, sizeof( int ) );
	sc->prio = calloc( nstep, sizeof( double ) );
	sc->started = calloc( nstep, sizeof( char ) );
	sc->nheap = 0;
	sc->ntaken = 0;
	sc->ndone = 0;
	sc->inorder = inorder;
	sc->optimistic = optimistic && !inorder;
	sc->abort = 0;
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
	for( l=nstep-1; l>=0; l-- ) // parents always have larger indices
		sc->prio[l] = pow( size[l], 0.75 ) + ( sc->parent[l] == -1 ? 0.0 : sc->prio[sc->parent[l]] );
	for( l=0; l<nstep; l++ )
		if( tbsched_isready( dep, l ) ) tbsched_push( sc, l );
	free( size );
}

static void tbsched_free( tbsched_t *sc )
{
	free( sc->started );
	free( sc->heap );
	free( sc->parent );
	free( sc->prio );
}

static void tbsched_done( tbsched_t *sc, Treedep *dep, int l ) // call with the mutex held, after dep[l].done = 1
{
	int p = sc->parent[l];
	sc->ndone++;
	if( p != -1 && tbsched_isready( dep, p ) ) tbsched_push( sc, p );
}

/* work/prog2: copying the member rows in and out of the local copies, split among idle helpers
   for large groups (each row is copied by exactly one thread; strcpy as before). */
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
				if( sc->nheap > 0 ) break;
				if( dphelppool && dphelp_serve_locked() ) continue;
				if( dphelppool ? ( sc->ndone == njob-1 ) : ( sc->ntaken == njob-1 ) ) break;
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
			}
		}
		else
			l = *jobpospt;
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
		while( *nrunpt >= nthreadtb ) // tabun iranai
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

		len1 = strlen( aseq[m1] );
		len2 = strlen( aseq[m2] );
		if( *alloclen <= len1 + len2 )
		{
			reporterr(       "\nReallocating.." );
			pthread_rwlock_wrlock( targ->aseqlock ); // wait until no thread is copying rows of aseq
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  ); 
			pthread_rwlock_unlock( targ->aseqlock );
			reporterr(       "done. *alloclen = %d\n", *alloclen );
		}
		localcopylen = *alloclen; // work/prog: the member rows are copied after the mutex is released

		/* work/prog2: with the ready queue, steps run out of index order, and the one-way switch
		   below would then take effect at a different step than in the original in-order dispatch.
		   So the first step that would flip it abandons the pass instead; disttbfast() then runs the
		   pass again from the start with the original dispatch (bitwise the same as before as long
		   as the switch never fires, since until then every step uses the same aligner). */
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

		if( alg == 'M' ) // hoka no thread ga M ni shitakamo shirenainode
		{
//			reporterr(       "Freeing commonIP (thread %d)\n", thread_no );
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}

		pthread_mutex_unlock( targ->mutex );

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
			if( alg == 'M' )
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
//						reporterr(       "%d-%d", clus1, clus2 );
						if( l < 500 || l % 100 == 0 ) if( ( cpmxchild1 && *cpmxchild1 ) || ( cpmxchild0 && *cpmxchild0 ) ) reporterr(       " h" );
						pscore = A__align( dynamicmtx, penalty, penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, -1, -1, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					}
					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}
#ifdef PROGTIMING
		PTT(4);
		if( ptlog ) { int k; for( k=0; k<PT_NACC; k++ ) ptlog[l].acc[k] = pt_acc[k]; ptlog[l].l = l; ptlog[l].thread = thread_no; ptlog[l].clus1 = clus1; ptlog[l].clus2 = clus2; ptlog[l].len1 = len1; ptlog[l].len2 = len2; ptlog[l].c0 = dep[l].child0; ptlog[l].c1 = dep[l].child1; }
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

		pthread_mutex_lock( targ->mutex );
		dep[l].done = 1;
		if( targ->sched ) tbsched_done( targ->sched, dep, l );
		PTT(5);
		(*nrunpt)--;
		pthread_cond_broadcast( targ->treecond );

//		reporterr( "at step %d\n", l );
//		use_getrusage();

		pthread_mutex_unlock( targ->mutex );
		PTT(6);



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
	if( njob > 20000 ) nthreadtb = 0; 
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

		tmpseq = AllocateCharVec( nlenmax+1 );
		grpseq = AllocateIntVec( nlenmax+1 );
		pointt = AllocateIntMtx( njob, nlenmax+1 );
		if( !compacttree ) mtx = AllocateFloatHalfMtx( njob ); 
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
				fixed_musclesupg_double_realloc_nobk_halfmtx_memsave_mt( njob, mtx, topol, len, dep, 1, !calcpairdists, nthreadpair ); // fast-mafft: same tree, faster (mltaln9.c)
#else
				fixed_musclesupg_double_realloc_nobk_halfmtx_memsave_mt( njob, mtx, topol, len, dep, 1, !calcpairdists, 1 ); // fast-mafft: same tree, faster (mltaln9.c)
#endif
				if( !calcpairdists )
				{
					FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
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
			tbsched_t tbsched;
			pthread_rwlock_t aseqlock;
			dphelppool_t tbpool;
	
			nthread_yoyu = nthreadtb * 1;
			nrun = 0;
			jobpos = 0; 
			targ = calloc( nthread_yoyu, sizeof( treebasethread_arg_t ) ); 
			fftlog = AllocateIntVec( njob );
			handle = calloc( nthread_yoyu, sizeof( pthread_t ) ); 
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
			int tb_inorder = ( nthreadtb <= 1 || tb_abandonedbefore );
#ifdef PROGTIMING
			if( nthreadtb <= 1 && getenv( "MAFFT_READYQ1" ) && !tb_abandonedbefore ) tb_inorder = 0; // A/B: optimistic ready queue at one thread
#endif
			int tb_optimistic = ( !nevermemsave && alg != 'M' );
			int tb_attempt;
			char **tb_saveseq = NULL;
			int *tb_savenlen = NULL;
			if( !tb_inorder && tb_optimistic ) // what an abandoned attempt changes, to put it back
			{
				tb_saveseq = calloc( njob, sizeof( char * ) );
				for( i=0; i<njob; i++ ) { tb_saveseq[i] = malloc( strlen( bseq[i] ) + 1 ); strcpy( tb_saveseq[i], bseq[i] ); }
				tb_savenlen = calloc( njob, sizeof( int ) );
				for( i=0; i<njob; i++ ) tb_savenlen[i] = nlen[i];
			}
			for( tb_attempt=0; ; tb_attempt++ )
			{
			nrun = 0;
			jobpos = 0;
			for( i=0; i<njob; i++ ) dep[i].done = 0; 
			for( i=0; i<njob; i++ ) fftlog[i] = 1; 
			tbsched_init( &tbsched, njob, dep, tb_inorder, tb_optimistic );
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
				dphelppool = &tbpool;
			}
#ifdef PROGTIMING
			if( getenv( "MAFFT_PROGLOG" ) ) { ptlog = calloc( njob, sizeof( ptrec_t ) ); pt_base = pt_now(); }
			struct rusage pt_ru0; double pt_w0 = pt_now(); getrusage( RUSAGE_SELF, &pt_ru0 );
#endif
	
			for( i=0; i<nthread_yoyu; i++ )
			{
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
	
				pthread_create( handle+i, NULL, treebasethread, (void *)(targ+i) );
			}

			for( i=0; i<nthread_yoyu; i++ )
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
			if( !tbsched.abort )
			{
				tbsched_free( &tbsched );
				break;
			}
			/* work/prog2: abandoned (see above): undo the attempt and run the pass again in order.
			   The aligned rows and nlen go back to their saved values; profile caches and member
			   lists are dropped (a member list is still owned by its step when that step was handed
			   out but its parent was not); done flags and fftlog are reset at the top. */
			reporterr( "\nAn alignment exceeds 30000 columns: redoing this pass with the original step order.\n" );
			for( i=0; i<njob; i++ ) strcpy( bseq[i], tb_saveseq[i] );
			for( i=0; i<njob; i++ ) nlen[i] = tb_savenlen[i];
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
			pthread_mutex_destroy( &mutex );
			pthread_cond_destroy( &treecond );
			pthread_rwlock_destroy( &aseqlock );
			free( handle );
			free( targ );
			free( fftlog );

		}
		else
#endif
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
			if( getenv( "MAFFT_MSADIST_SAVEMSA" ) ) msadist_savemsa( getenv( "MAFFT_MSADIST_SAVEMSA" ), njob, bseq, selfscore );
#ifdef enablemultithread
			msadistfast_t *msadistfast = msadistfast_prepare( njob, bseq, penalty_dist, nthreadpair ); // fast-mafft; NULL -> original scorer
#else
			msadistfast_t *msadistfast = msadistfast_prepare( njob, bseq, penalty_dist, 1 ); // fast-mafft; NULL -> original scorer
#endif
			skiptable = AllocateIntMtx( njob, 0 );
			if( msadistfast == NULL || msadistfast->check ) makeskiptable( njob, skiptable, bseq ); // allocate suru. (only the original scorer uses it)
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
	if( mtx ) FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
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
