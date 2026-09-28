#include "mltaln.h"

#define DEBUG 0


static char *comment;
static char *orderfile;
static int format;
static int namelen;
static int excludedashseq;
static int extendedalphabet;

static void fillspace( char *seq, int lenmax )
{
	int len = strlen( seq );
	seq += len;
	lenmax -= len;
	while( lenmax-- ) *seq++ = ' ';
	*seq = 0;
}

void setmark_clustal( int nlen, int nseq, char **seq, char *mark )
{
	int i, j, k, nalpha;
	char firstletter;
	char *strong[9];
	char *weaker[11];
	int nstrong, nweaker;
	char s;

	if( dorp == 'd' ) 
	{
		strong[0] = "TU";
		nstrong = 1;
		weaker[0] = "AG";
		weaker[1] = "CT";
		weaker[2] = "CU";
		nweaker = 2;
		nalpha = 10;
	}
	else
	{
		strong[0] = "STA";
		strong[1] = "NEQK";
		strong[2] = "NHQK";
		strong[3] = "NDEQ";
		strong[4] = "QHRK";
		strong[5] = "MILV";
		strong[6] = "MILF";
		strong[7] = "HY";
		strong[8] = "FYW";
		nstrong = 9;
		weaker[0] = "CSA";
		weaker[1] = "ATV";
		weaker[2] = "SAG";
		weaker[3] = "STNK";
		weaker[4] = "STPA";
		weaker[5] = "SGND";
		weaker[6] = "SNDEQK";
		weaker[7] = "NDEQHK";
		weaker[8] = "NEQHRK";
		weaker[9] = "FVLIM";
		weaker[10] = "HFY";
		nweaker = 11;
		nalpha = 20;
	}

	for( i=0; i<nlen; i++ )
	{
		mark[i] = ' ';
		for( j=0; j<nseq; j++ )
		{
			s = seq[j][i];
			if( '-' == s || ' ' == s ) break;
		}
		if( j != nseq ) 
		{
			continue;
		}
		if( extendedalphabet )
		{
			firstletter = seq[0][i];
			if( amino_n[(unsigned char)firstletter] < 0 ) continue;
	
			for( j=0; j<nseq; j++ )
				if( seq[j][i] != firstletter ) break;
			if( j == nseq ) 
			{
				mark[i] = '*';
				continue;
			}
		}
		else 
		{
			firstletter = toupper( seq[0][i] );
			if( amino_n[(unsigned char)firstletter] >= nalpha || amino_n[(unsigned char)firstletter] < 0 ) continue;
	
			for( j=0; j<nseq; j++ )
				if( toupper( seq[j][i] ) != firstletter ) break;
			if( j == nseq ) 
			{
				mark[i] = '*';
				continue;
			}
			for( k=0; k<nstrong; k++ )
			{
				for( j=0; j<nseq; j++ )
				{
					if( !strchr( strong[k], toupper( seq[j][i] ) ) ) break;
				}
				if( j == nseq ) break;
			}
			if( k < nstrong )
			{
				mark[i] = ':';
				continue;
			}
			for( k=0; k<nweaker; k++ )
			{
				for( j=0; j<nseq; j++ )
				{
					if( !strchr( weaker[k], toupper( seq[j][i] ) ) ) break;
				}
				if( j == nseq ) break;
			}
			if( k < nweaker )
			{
				mark[i] = '.';
				continue;
			}
		}
	}
	mark[nlen] = 0;
}

void setmark( int nlen, int nseq, char **seq, char *mark )
{
	int i, j;

	for( i=0; i<nlen; i++ )
	{
		mark[i] = ' ';
		for( j=0; j<nseq; j++ )
			if( '-' == seq[j][i] ) break;
		if( j != nseq ) 
		{
			continue;
		}
		for( j=0; j<nseq; j++ )
			if( seq[0][i] != seq[j][i] ) break;
		if( j == nseq ) 
		{
			mark[i] = '*';
			continue;
		}
		for( j=0; j<nseq; j++ )
			if( amino_grp[(unsigned char)seq[0][i]] != amino_grp[(unsigned char)seq[j][i]] ) break;
		if( j == nseq ) 
		{
			mark[i] = '.';
			continue;
		}
	}
	mark[nlen] = 0;
}

void arguments( int argc, char *argv[] )
{
    int c;
	excludedashseq = 0;
	namelen = -1;
	scoremtx = 1;
	nblosum = 62;
	LineLengthInFASTA = -1;
	dorp = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	inputfile = NULL;
	comment = NULL;
	orderfile = NULL;
	format = 'c';
	extendedalphabet = 0;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'c':
					comment = *++argv;
					fprintf( stderr, "comment = %s\n", comment );
					--argc;
					goto nextoption;
				case 'r':
					orderfile = *++argv;
					fprintf( stderr, "orderfile = %s\n", orderfile );
					--argc;
					goto nextoption;
				case 'n':
					namelen = myatoi( *++argv );
					fprintf( stderr, "namelen = %d\n", namelen );
					--argc;
					goto nextoption;
				case 'l':
					LineLengthInFASTA = myatoi( *++argv );
					if( LineLengthInFASTA == 0 )
					{
						reporterr( "LineLengthInFASTA = 0 ??\n" );
						exit( 1 );
					}
					--argc;
					goto nextoption;
				case 'f':
					format = 'f';
					break;
				case 'd':
					excludedashseq = 1;
					break;
				case 'y':
					format = 'y';
					break;
				case 'E':
					extendedalphabet = 1;
					nblosum = -2;
					break;
				case 'N':
					extendedalphabet = 0;
					break;
                default:
                    fprintf( stderr, "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
		}
		nextoption:
			;
	}
    if( argc != 0 ) 
    {
        fprintf( stderr, "options: Check source file !\n" );
        exit( 1 );
    }
}

#if !defined(mingw) && !defined(_MSC_VER) && !defined(_WIN32)
#include <sys/stat.h>
/* fast-mafft (overhead lever): f2cl reads its input three times with getc()
   (countKUorWA, getnumlen_casepreserve, readData_pointer_casepreserve). For the plain FASTA
   that disttbfast writes, the whole file is read once into memory and split there instead.
   The result (njob, nlenmax, the dorp guess, names and sequences) is what the original
   functions produce for such input. Anything else (not a regular file, not at offset 0,
   junk before the first '>', a name line without '\n', NUL bytes, a sequence line with
   '=', '<' or '>', a sequence of N or more bytes) is rewound and read the original way. */
typedef struct
{
	char *buf;
	size_t len;
	int n;
	size_t *nameb, *namee, *bodyb, *bodye;
	int *flen;
} f2cl_fast_t;

/* Counts the bytes charfilter() removes (LF, space, CR) and tells whether there is one that
   makes it stop or exit (NUL, '=', '<', '>'). Branch-free so that the compiler vectorizes it. */
static int f2cl_scanbody( const unsigned char *p, size_t n, size_t *cntp )
{
	size_t k, cnt = 0;
	unsigned int bad = 0;
	for( k=0; k<n; k++ )
	{
		unsigned char c = p[k];
		bad |= ( c == 0 ) | ( c == '=' ) | ( c == '<' ) | ( c == '>' );
		cnt += ( c == 0x0a ) | ( c == 0x20 ) | ( c == 0x0d );
	}
	*cntp = cnt;
	return( bad != 0 );
}

static void f2cl_fastfree( f2cl_fast_t *f )
{
	if( !f ) return;
	free( f->buf ); free( f->nameb ); free( f->namee ); free( f->bodyb ); free( f->bodye ); free( f->flen );
	free( f );
}

static f2cl_fast_t *f2cl_fastscan( FILE *fp, int *nlenminpt )
{
	struct stat st;
	f2cl_fast_t *f;
	size_t cap, got, p, e, bs, be, q, k, cnt;
	int i, cap_n, total, nsite, atgcnum, len1, fl;
	char *tmp, *t, *nl;
	double atgcfreq;

	if( fstat( fileno( fp ), &st ) != 0 || !S_ISREG( st.st_mode ) ) return( NULL );
	if( ftell( fp ) != 0 ) return( NULL );

	f = calloc( 1, sizeof( f2cl_fast_t ) );
	if( !f ) return( NULL );
	cap = (size_t)st.st_size + 1;
	f->buf = malloc( cap );
	if( !f->buf ) { free( f ); return( NULL ); }
	f->len = 0;
	while( ( got = fread( f->buf + f->len, 1, cap - f->len, fp ) ) > 0 )
	{
		f->len += got;
		if( f->len == cap )
		{
			char *nb = realloc( f->buf, cap * 2 );
			if( !nb ) goto fallback;
			f->buf = nb;
			cap *= 2;
		}
	}
	if( ferror( fp ) ) goto fallback;

	if( f->len == 0 || f->buf[0] != '>' ) goto fallback;
	cap_n = 1024;
	f->nameb = malloc( cap_n * sizeof( size_t ) );
	f->namee = malloc( cap_n * sizeof( size_t ) );
	f->bodyb = malloc( cap_n * sizeof( size_t ) );
	f->bodye = malloc( cap_n * sizeof( size_t ) );
	f->flen = malloc( cap_n * sizeof( int ) );
	if( !f->nameb || !f->namee || !f->bodyb || !f->bodye || !f->flen ) goto fallback;
	f->n = 0;
	p = 0;
	while( p < f->len ) // f->buf[p] is a '>' at the start of a line
	{
		nl = memchr( f->buf + p + 1, '\n', f->len - p - 1 );
		if( !nl ) goto fallback; // the original myfgets() would not stop here
		e = nl - f->buf;
		if( memchr( f->buf + p + 1, 0, e - p - 1 ) ) goto fallback;
		bs = e + 1;
		if( bs < f->len && f->buf[bs] == '>' ) // empty sequence
			be = bs;
		else
		{
			be = f->len;
			for( q=bs; q<f->len; q++ )
			{
				nl = memchr( f->buf + q, '\n', f->len - q );
				if( !nl ) break;
				q = nl - f->buf;
				if( q + 1 < f->len && f->buf[q+1] == '>' ) { be = q + 1; break; }
			}
		}
		if( be - bs >= N ) goto fallback; // load1SeqWithoutName_realloc_casepreserve() would reallocate
		if( f2cl_scanbody( (unsigned char *)f->buf + bs, be - bs, &cnt ) ) goto fallback; // charfilter() stops or exits
		if( f->n == cap_n )
		{
			cap_n *= 2;
			if( !( f->nameb = realloc( f->nameb, cap_n * sizeof( size_t ) ) ) ) goto fallback;
			if( !( f->namee = realloc( f->namee, cap_n * sizeof( size_t ) ) ) ) goto fallback;
			if( !( f->bodyb = realloc( f->bodyb, cap_n * sizeof( size_t ) ) ) ) goto fallback;
			if( !( f->bodye = realloc( f->bodye, cap_n * sizeof( size_t ) ) ) ) goto fallback;
			if( !( f->flen = realloc( f->flen, cap_n * sizeof( int ) ) ) ) goto fallback;
		}
		f->nameb[f->n] = p + 1;
		f->namee[f->n] = e;
		f->bodyb[f->n] = bs;
		f->bodye[f->n] = be;
		f->flen[f->n] = (int)( be - bs - cnt );
		f->n++;
		p = be;
	}

	// getnumlen_casepreserve()
	njob = f->n;
	nlenmax = 0;
	*nlenminpt = 99999999;
	atgcnum = 0;
	total = 0;
	nsite = 0;
	tmp = NULL;
	for( i=0; i<njob; i++ )
	{
		fl = f->flen[i];
		if( fl > nlenmax ) nlenmax = fl;
		if( fl < *nlenminpt ) *nlenminpt = fl;
		if( total < 1000000 )
		{
			if( !( t = realloc( tmp, fl + 1 ) ) ) { free( tmp ); goto fallback; }
			tmp = t;
			len1 = 0;
			for( k=f->bodyb[i]; k<f->bodye[i]; k++ )
			{
				char c = f->buf[k];
				if( c != 0x0a && c != 0x20 && c != 0x0d ) tmp[len1++] = c;
			}
			tmp[len1] = 0;
			atgcnum += countATGC( tmp, &nsite );
			total += nsite;
		}
	}
	free( tmp );
	atgcfreq = (double)atgcnum / total;
	if( dorp == NOTSPECIFIED )
	{
		if( atgcfreq > 0.75 ) 	
			dorp = 'd';
		else                  
			dorp = 'p';
	}
	return( f );

fallback:
	f2cl_fastfree( f );
	clearerr( fp );
	rewind( fp );
	return( NULL );
}

/* readData_pointer_casepreserve() on the scanned input */
static void f2cl_fastfill( f2cl_fast_t *f, char **name, int *nlen, char **seq )
{
	int i;
	size_t k, l;
	char *d;
	for( i=0; i<f->n; i++ )
	{
		name[i][0] = '=';
		l = f->namee[i] - f->nameb[i];
		if( l > B-2 ) l = B-2; // myfgets( name[i]+1, B-2, fp )
		memcpy( name[i]+1, f->buf + f->nameb[i], l );
		name[i][1+l] = 0;
		d = seq[i];
		l = f->bodye[i] - f->bodyb[i];
		if( (size_t)f->flen[i] + 1 == l && f->buf[f->bodye[i]-1] == '\n' ) // one line: nothing to filter but its LF
		{
			memcpy( d, f->buf + f->bodyb[i], l - 1 );
			d += l - 1;
		}
		else
		{
			for( k=f->bodyb[i]; k<f->bodye[i]; k++ )
			{
				char c = f->buf[k];
				if( c != 0x0a && c != 0x20 && c != 0x0d ) *d++ = c;
			}
		}
		*d = 0;
		nlen[i] = (int)( d - seq[i] );
	}
}
#endif

int main( int argc, char *argv[] )
{
	static int  *nlen, *onlen;	
	static char **name, **oname, **seq, **oseq, *mark;
	static int *order, *oorder;
	int i, j;
	FILE *infp;
	FILE *orderfp;
	char gett[B];
	int nlenmin;
	int nout;
#if !defined(mingw) && !defined(_MSC_VER) && !defined(_WIN32)
	f2cl_fast_t *fastin = NULL;
#endif

#if !defined(mingw) && !defined(_MSC_VER) && !defined(_WIN32)
	/* fast-mafft: one write() per 1 MB of output instead of one per sequence; same bytes */
	if( !isatty( fileno( stdout ) ) ) setvbuf( stdout, NULL, _IOFBF, 1<<20 );
#endif
	arguments( argc, argv );


	if( inputfile )
	{
		infp = fopen( inputfile, "rb" );
		if( !infp )
		{
			fprintf( stderr, "Cannot open %s\n", inputfile );
			exit( 1 );
		}
	}
	else
		infp = stdin;

#if !defined(mingw) && !defined(_MSC_VER) && !defined(_WIN32)
	fastin = f2cl_fastscan( infp, &nlenmin );
	if( !fastin )
#endif
	{
	getnumlen_casepreserve( infp, &nlenmin );
	rewind( infp );
	}

#if !defined(mingw) && !defined(_MSC_VER) && !defined(_WIN32)
	if( fastin && format == 'f' && !excludedashseq )
	{
		/* fast-mafft: in this case every later use of seq[i] stops at its NUL (the writer, and
		   constants()), so a row of strlen+1 bytes holds all that is ever read. This skips
		   zeroing njob*(nlenmax*2+1) bytes. FreeCharMtx() frees these rows the same way. */
		seq = (char **)calloc( njob+1, sizeof( char * ) );
		if( seq == NULL ) { fprintf( stderr, "Cannot allocate seq\n" ); exit( 1 ); }
		for( i=0; i<njob; i++ )
		{
			seq[i] = (char *)malloc( fastin->flen[i] + 1 );
			if( seq[i] == NULL ) { fprintf( stderr, "Cannot allocate seq\n" ); exit( 1 ); }
		}
		seq[njob] = NULL;
	}
	else
#endif
	seq = AllocateCharMtx( njob, nlenmax*2+1 );
	mark = AllocateCharVec( nlenmax*2+1 );
	order = AllocateIntVec( njob );
	name = AllocateCharMtx( njob, B+1 );
    nlen = AllocateIntVec( njob );

	if( orderfile )
	{
		orderfp = fopen( orderfile, "r" );
		if( !orderfp )
		{
			fprintf( stderr, "Cannot open %s\n", orderfile );
			exit( 1 );
		}
		for( i=0; i<njob; i++ )
		{
			fgets( gett, B-1, orderfp );
			order[i] = atoi( gett );
		}
		fclose( orderfp );
	}
	else
	{
		for( i=0; i<njob; i++ ) order[i] = i;
	}


#if !defined(mingw) && !defined(_MSC_VER) && !defined(_WIN32)
	if( fastin )
	{
		f2cl_fastfill( fastin, name, nlen, seq );
		f2cl_fastfree( fastin ); fastin = NULL;
	}
	else
#endif
	readData_pointer_casepreserve( infp, name, nlen, seq );
	fclose( infp );

	if( format == 'c' || format == 'y' ) for( i=0; i<njob; i++ ) fillspace( seq[i], nlenmax );
	constants( njob, seq );

//	initSignalSM();

//	initFiles();

	if( excludedashseq )
	{
		nout = 0;
		for( i=0; i<njob; i++ )
		{
			if( !strncmp( name[order[i]]+1, "DASH|", 5 ) || strstr( name[order[i]]+1, "_numo_e_DASH|" ) ) continue;
			nout++;
		}

		reporterr( "nout=%d\n", nout );

		oseq = calloc( sizeof( char * ), nout );
		oname = calloc( sizeof( char * ), nout );
   		onlen = AllocateIntVec( nout );
   		oorder = AllocateIntVec( nout );


		for( i=0,j=0; i<njob; i++ )
		{
			if( !strncmp( name[order[i]]+1, "DASH|", 5 ) || strstr( name[order[i]]+1, "_numo_e_DASH|" ) ) continue;
			oname[j] = name[order[i]];
			oseq[j] = seq[order[i]];
			onlen[j] = nlen[order[i]];
   			oorder[j] = j;
			j++;
		}
		commongappick( nout, oseq );
	}
	else
	{
		oseq = seq;
		oname = name;
   		onlen = nlen;
   		oorder = order;
		nout = njob;
	}

	reporterr( "seq = %p, oseq=%p\n", seq, oseq );



//	setmark( nlenmax, njob, seq, mark );
	/* fast-mafft (overhead lever): the conservation marks are printed only in the CLUSTAL
	   format ('c'); the FASTA ('f') and PHYLIP ('y') writers never read them. */
	if( format == 'c' ) setmark_clustal( nlenmax, nout, oseq, mark );

#if mingw
	setmode( fileno( stdout ), O_TEXT ); // windows deha saishuu tekina output nomi text mode
#endif

	if( format == 'f' )
		writeData_reorder_pointer( stdout, nout, oname, onlen, oseq, oorder );
	else if( format == 'c' )
		clustalout_pointer( stdout, nout, nlenmax, oseq, oname, mark, comment, oorder, namelen );
	else if( format == 'y' )
		phylipout_pointer( stdout, nout, nlenmax, oseq, oname, oorder, namelen );
	else
		fprintf( stderr, "Unknown format\n" );

	FreeCharMtx( seq ); seq = NULL;
	FreeCharMtx( name ); name = NULL;
	free( nlen ); nlen = NULL;
	free( order ); order = NULL;

	if( excludedashseq )
	{
		free( oseq ); oseq = NULL;
		free( oname ); oname = NULL;
		free( onlen ); onlen = NULL;
		free( oorder ); oorder = NULL;
	}

	free( mark ); mark = NULL;
	freeconstants();

//	SHOWVERSION;
	return( 0 );
}
