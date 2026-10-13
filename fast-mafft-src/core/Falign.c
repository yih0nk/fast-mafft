#include "mltaln.h"
#include "progtiming.h"
#include "dphelp.h"
#ifdef PROGTIMING
TLS double pt_acc[PT_NACC];
#endif

#if 0
static FILE *fftfp;
#endif
static TLS int n20or4or2;

#define KEIKA 0
#define RND   0
#define DEBUG 0


#if RND // by D.Mathog
static void generateRndSeq( char *seq, int len )
{
	while( len-- )
#if 1
		*seq++ = (int)( rnd() * n20or4or2 );
#else
		*seq++ = (int)1;
#endif
}
#endif

static void vec_init( Fukusosuu *result, int nlen )
{
	while( nlen-- )
	{
		result->R = result->I = 0.0;
		result++;
	}
}

#if 0 // by D.Mathog
static void vec_init2( Fukusosuu **result, char *seq, double eff, int st, int ed )
{
	int i;
	for( i=st; i<ed; i++ )
		result[(int)*seq++][i].R += eff;
}
#endif

static void seq_vec_2( Fukusosuu *result, double *score, double incr, char *seq )
{
	static TLS int n;
	for( ; *seq; result++ )
	{
		n = amino_n[(int)*seq++];
		if( n < 20 && n >= 0 ) result->R += incr * score[n];
#if 0
		fprintf( stderr, "n=%d, score=%f, inc=%f R=%f\n",n,  score[n], incr * score[n], result->R );
#endif
	}
}

static void seq_vec_3( Fukusosuu **result, double incr, char *seq )
{
	int i;
	int n;
	for( i=0; *seq; i++ )
	{
		n = amino_n[(int)*seq++];
		if( n < n20or4or2 && n >= 0 ) result[n][i].R += incr;
	}
}

static void seq_vec_5( Fukusosuu *result, double *score1, double *score2, double incr, char *seq )
{
	int n;
	for( ; *seq; result++ )
	{
		n = amino_n[(int)*seq++];
		if( n > 20 ) continue;
		result->R += incr * score1[n];
		result->I += incr * score2[n];
#if 0
		fprintf( stderr, "n=%d, score=%f, inc=%f R=%f\n",n,  score[n], incr * score[n], result->R );
#endif
	}
}


static void seq_vec_4( Fukusosuu *result, double incr, char *seq )
{
	char s;
	for( ; *seq; result++ )
	{
		s = *seq++;
		if( s == 'a' )
			result->R += incr;
		else if( s == 't' )
			result->R -= incr;
		else if( s == 'g' )
			result->I += incr;
		else if( s == 'c' )
			result->I -= incr;
	}
}

#if 0 // by D.Mathog
static void seq_vec( Fukusosuu *result, char query, double incr, char *seq )
{
#if 0
	int bk = nlen;
#endif
	while( *seq )
	{
		if( *seq++ == query ) result->R += incr;
		result++;
#if 0
fprintf( stderr, "i = %d result->R = %f\n", bk-nlen, (result-1)->R );
#endif
	}
}

static int checkRepeat( int num, int *cutpos )
{
	int tmp, buf;

	buf = *cutpos;
	while( num-- )
	{
		if( ( tmp = *cutpos++ ) < buf ) return( 1 );
		buf = tmp;
	}
	return( 0 );
}

static int segcmp( void *ptr1, void *ptr2 )
{
	int diff;
	Segment **seg1 = (Segment **)ptr1;
	Segment **seg2 = (Segment **)ptr2;
#if 0
	return( (*seg1)->center - (*seg2)->center );
#else
	diff = (*seg1)->center - (*seg2)->center;
	if( diff ) return( diff );

	diff = (*seg1)->start - (*seg2)->start;
	if( diff ) return( diff );

	diff = (*seg1)->end - (*seg2)->end;
	if( diff ) return( diff );

	fprintf( stderr, "USE STABLE SORT !!\n" );
	exit( 1 );
	return( 0 );
#endif
}
#endif


static void mymergesort( int first, int last, Segment **seg )
{
	int middle;
	static TLS int i, j, k, p;
	static TLS int allo = 0;
	static TLS Segment **work = NULL;

	if( seg == NULL )
	{
		if( work ) free( work ); 
		work = NULL;
		allo = 0;
		return;
	}

	if( last > allo )
	{
		allo = last;
		if( work ) free( work );
		work = (Segment **)calloc( allo / 2 + 1, sizeof( Segment *) );
	}

	if( first < last )
	{
		middle = ( first + last ) / 2;
		mymergesort( first, middle, seg );
		mymergesort( middle+1, last, seg );
		p = 0;
		for( i=first; i<=middle; i++ ) work[p++] = seg[i];
		i = middle + 1; j = 0; k = first;
		while( i <= last && j < p )
		{
			if( work[j]->center <= seg[i]->center ) 
				seg[k++] = work[j++];
			else
				seg[k++] = seg[i++];
		}
		while( j < p ) seg[k++] = work[j++];
	}
}


double Fgetlag( 
				double **n_dynamicmtx, 
				char  **seq1, char  **seq2, 
			    double *eff1, double *eff2, 
			    int    clus1, int    clus2,
			    int alloclen )
{
	int i, j, k, l, m;
	int nlen, nlen2, nlen4;
	static TLS int crossscoresize = 0;
	static TLS char **tmpseq1 = NULL;
	static TLS char **tmpseq2 = NULL;
	static TLS char **tmpptr1 = NULL;
	static TLS char **tmpptr2 = NULL;
	static TLS char **tmpres1 = NULL;
	static TLS char **tmpres2 = NULL;
	static TLS char **result1 = NULL;
	static TLS char **result2 = NULL;
#if RND
	static TLS char **rndseq1 = NULL;
	static TLS char **rndseq2 = NULL;
#endif
	static TLS Fukusosuu **seqVector1 = NULL;
	static TLS Fukusosuu **seqVector2 = NULL;
	static TLS Fukusosuu **naiseki = NULL;   
	static TLS Fukusosuu *naisekiNoWa = NULL; 
	static TLS double *soukan = NULL;
	static TLS double **crossscore = NULL;
	int nlentmp;
	static TLS int *kouho = NULL;
	static TLS Segment *segment = NULL;
	static TLS Segment *segment1 = NULL;
	static TLS Segment *segment2 = NULL;
	static TLS Segment **sortedseg1 = NULL;
	static TLS Segment **sortedseg2 = NULL;
	static TLS int *cut1 = NULL;
	static TLS int *cut2 = NULL;
	static TLS int localalloclen = 0;
	int lag;
	int tmpint;
	int count, count0;
	int len1, len2;
	int totallen;
	double dumdb = 0.0;
	int headgp, tailgp;

	len1 = strlen( seq1[0] );
	len2 = strlen( seq2[0] );
	nlentmp = MAX( len1, len2 );

	nlen = 1;
	while( nlentmp >= nlen ) nlen <<= 1;
#if 0
	fprintf( stderr, "###   nlen    = %d\n", nlen );
#endif

	nlen2 = nlen/2; nlen4 = nlen2 / 2;

#if DEBUG
	fprintf( stderr, "len1 = %d, len2 = %d\n", len1, len2 );
	fprintf( stderr, "nlentmp = %d, nlen = %d\n", nlentmp, nlen );
#endif

	if( !localalloclen )
	{
		kouho = AllocateIntVec( NKOUHO );
		cut1 = AllocateIntVec( MAXSEG );
		cut2 = AllocateIntVec( MAXSEG );
		tmpptr1 = AllocateCharMtx( njob, 0 );
		tmpptr2 = AllocateCharMtx( njob, 0 );
		result1 = AllocateCharMtx( njob, alloclen );
		result2 = AllocateCharMtx( njob, alloclen );
		tmpres1 = AllocateCharMtx( njob, alloclen );
		tmpres2 = AllocateCharMtx( njob, alloclen );
//		crossscore = AllocateDoubleMtx( MAXSEG, MAXSEG );
		segment = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		segment1 = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		segment2 = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		sortedseg1 = (Segment **)calloc( MAXSEG, sizeof( Segment * ) );
		sortedseg2 = (Segment **)calloc( MAXSEG, sizeof( Segment * ) );
		if( !( segment && segment1 && segment2 && sortedseg1 && sortedseg2 ) )
			ErrorExit( "Allocation error\n" );

		if     ( scoremtx == -1 ) n20or4or2 = 4;
		else if( fftscore == 1  ) n20or4or2 = 2;
		else                      n20or4or2 = 20;
	}
	if( localalloclen < nlen )
	{
		if( localalloclen )
		{
#if 1
			FreeFukusosuuMtx ( seqVector1 );
			FreeFukusosuuMtx ( seqVector2 );
			FreeFukusosuuVec( naisekiNoWa );
			FreeFukusosuuMtx( naiseki );
			FreeDoubleVec( soukan );
			FreeCharMtx( tmpseq1 );
			FreeCharMtx( tmpseq2 );
#endif
#if RND
			FreeCharMtx( rndseq1 );
			FreeCharMtx( rndseq2 );
#endif
		}


		tmpseq1 = AllocateCharMtx( njob, nlen );
		tmpseq2 = AllocateCharMtx( njob, nlen );
		naisekiNoWa = AllocateFukusosuuVec( nlen );
		naiseki = AllocateFukusosuuMtx( n20or4or2, nlen );
		seqVector1 = AllocateFukusosuuMtx( n20or4or2+1, nlen+1 );
		seqVector2 = AllocateFukusosuuMtx( n20or4or2+1, nlen+1 );
		soukan = AllocateDoubleVec( nlen+1 );

#if RND
		rndseq1 = AllocateCharMtx( njob, nlen );
		rndseq2 = AllocateCharMtx( njob, nlen );
		for( i=0; i<njob; i++ )
		{
			generateRndSeq( rndseq1[i], nlen );
			generateRndSeq( rndseq2[i], nlen );
		}
#endif
		localalloclen = nlen;
	}
	
	for( j=0; j<clus1; j++ ) strcpy( tmpseq1[j], seq1[j] );
	for( j=0; j<clus2; j++ ) strcpy( tmpseq2[j], seq2[j] );

#if 0
fftfp = fopen( "input_of_Falign", "w" );
fprintf( fftfp, "nlen = %d\n", nlen );
fprintf( fftfp, "seq1: ( %d sequences ) \n", clus1 );
for( i=0; i<clus1; i++ )
	fprintf( fftfp, "%s\n", seq1[i] );
fprintf( fftfp, "seq2: ( %d sequences ) \n", clus2 );
for( i=0; i<clus2; i++ )
	fprintf( fftfp, "%s\n", seq2[i] );
fclose( fftfp );
system( "less input_of_Falign < /dev/tty > /dev/tty" );
#endif

	if( fftkeika ) fprintf( stderr,  " FFT ... " );

	for( j=0; j<n20or4or2; j++ ) vec_init( seqVector1[j], nlen );
	if( fftscore && scoremtx != -1 )
	{
		for( i=0; i<clus1; i++ )
		{
			seq_vec_2( seqVector1[0], polarity, eff1[i], tmpseq1[i] );
			seq_vec_2( seqVector1[1], volume,   eff1[i], tmpseq1[i] );
		}
	}
	else
	{
#if 0
		for( i=0; i<clus1; i++ ) for( j=0; j<n20or4or2; j++ ) 
			seq_vec( seqVector1[j], amino[j], eff1[i], tmpseq1[i] );
#else
		for( i=0; i<clus1; i++ )
			seq_vec_3( seqVector1, eff1[i], tmpseq1[i] );
#endif
	}
#if RND
	for( i=0; i<clus1; i++ )
	{
		vec_init2( seqVector1, rndseq1[i], eff1[i], len1, nlen );
	}
#endif
#if 0
fftfp = fopen( "seqVec", "w" );
fprintf( fftfp, "before transform\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "nlen=%d\n", nlen );
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector1[k][l].R, seqVector1[k][l].I );
}
fclose( fftfp );
system( "less seqVec < /dev/tty > /dev/tty" );
#endif

	for( j=0; j<n20or4or2; j++ ) vec_init( seqVector2[j], nlen );
	if( fftscore && scoremtx != -1 )
	{
		for( i=0; i<clus2; i++ )
		{
			seq_vec_2( seqVector2[0], polarity, eff2[i], tmpseq2[i] );
			seq_vec_2( seqVector2[1], volume,   eff2[i], tmpseq2[i] );
		}
	}
	else
	{
#if 0
		for( i=0; i<clus2; i++ ) for( j=0; j<n20or4or2; j++ ) 
			seq_vec( seqVector2[j], amino[j], eff2[i], tmpseq2[i] );
#else
		for( i=0; i<clus2; i++ )
			seq_vec_3( seqVector2, eff2[i], tmpseq2[i] );
#endif
	}
#if RND
	for( i=0; i<clus2; i++ )
	{
		vec_init2( seqVector2, rndseq2[i], eff2[i], len2, nlen );
	}
#endif

#if 0
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "before fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

	for( j=0; j<n20or4or2; j++ )
	{
		fft( nlen, seqVector2[j], 0 );
		fft( nlen, seqVector1[j], 0 );
	}
#if 0
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "#after fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "#%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

	for( k=0; k<n20or4or2; k++ ) 
	{
		for( l=0; l<nlen; l++ ) 
			calcNaiseki( naiseki[k]+l, seqVector1[k]+l, seqVector2[k]+l );
	}
	for( l=0; l<nlen; l++ ) 
	{
		naisekiNoWa[l].R = 0.0;
		naisekiNoWa[l].I = 0.0;
		for( k=0; k<n20or4or2; k++ ) 
		{
			naisekiNoWa[l].R += naiseki[k][l].R;
			naisekiNoWa[l].I += naiseki[k][l].I;
		}
	}

#if 0
fftfp = fopen( "naisekiNoWa", "w" );
fprintf( fftfp, "#Before fft\n" );
for( l=0; l<nlen; l++ )
	fprintf( fftfp, "%d  %f %f\n", l, naisekiNoWa[l].R, naisekiNoWa[l].I ); 
fclose( fftfp );
system( "less naisekiNoWa < /dev/tty > /dev/tty " );
#endif

	fft( -nlen, naisekiNoWa, 0 );

	for( m=0; m<=nlen2; m++ ) 
		soukan[m] = naisekiNoWa[nlen2-m].R;
	for( m=nlen2+1; m<nlen; m++ ) 
		soukan[m] = naisekiNoWa[nlen+nlen2-m].R;

#if 0
fftfp = fopen( "naisekiNoWa", "w" );
fprintf( fftfp, "#After fft\n" );
for( l=0; l<nlen; l++ )
	fprintf( fftfp, "%d  %f\n", l, naisekiNoWa[l].R ); 
fclose( fftfp );
fftfp = fopen( "list.plot", "w"  );
fprintf( fftfp, "plot 'naisekiNoWa'\npause -1" );
fclose( fftfp );
system( "/usr/bin/gnuplot list.plot &" );
#endif
#if 0
fprintf( stderr, "frt write start\n" );
fftfp = fopen( "frt", "w" );
for( l=0; l<nlen; l++ )
	fprintf( fftfp, "%d  %f\n", l-nlen2, soukan[l] ); 
fclose( fftfp );
system( "less frt < /dev/tty > /dev/tty" );
#if 0
fftfp = fopen( "list.plot", "w"  );
fprintf( fftfp, "plot 'frt'\n pause +1" );
fclose( fftfp );
system( "/usr/bin/gnuplot list.plot" );
#endif
#endif


	getKouho( kouho, NKOUHO, soukan, nlen );

#if 0
	for( i=0; i<NKOUHO; i++ )
	{
		fprintf( stdout, "kouho[%d] = %d\n", i, kouho[i] );
	}
#endif

#if KEIKA
	fprintf( stderr, "Searching anchors ... " );
#endif
	count = 0;



#define CAND 0
#if CAND
	fftfp = fopen( "cand", "w" );
	fclose( fftfp );
#endif

	for( k=0; k<NKOUHO; k++ ) 
	{

		lag = kouho[k];
		zurasu2( lag, clus1, clus2, seq1, seq2, tmpptr1, tmpptr2 );
#if CAND
		fftfp = fopen( "cand", "a" );
		fprintf( fftfp, ">Candidate No.%d lag = %d\n", k+1, lag );
		fprintf( fftfp, "%s\n", tmpptr1[0] );
		fprintf( fftfp, ">Candidate No.%d lag = %d\n", k+1, lag );
		fprintf( fftfp, "%s\n", tmpptr2[0] );
		fprintf( fftfp, ">\n", k+1, lag );
		fclose( fftfp );
#endif
		tmpint = alignableReagion( clus1, clus2, tmpptr1, tmpptr2, eff1, eff2, segment+count );
		
		if( count+tmpint > MAXSEG -3 ) ErrorExit( "TOO MANY SEGMENTS.\n" );


		if( tmpint == 0 ) break; // 060430 iinoka ?
		while( tmpint-- > 0 )
		{
			if( lag > 0 )
			{
				segment1[count].start  = segment[count].start ;
				segment1[count].end    = segment[count].end   ;
				segment1[count].center = segment[count].center;
				segment1[count].score  = segment[count].score;

				segment2[count].start  = segment[count].start  + lag;
				segment2[count].end    = segment[count].end    + lag;
				segment2[count].center = segment[count].center + lag;
				segment2[count].score  = segment[count].score       ;
			}
			else
			{
				segment1[count].start  = segment[count].start  - lag;
				segment1[count].end    = segment[count].end    - lag;
				segment1[count].center = segment[count].center - lag;
				segment1[count].score  = segment[count].score       ;

				segment2[count].start  = segment[count].start ;
				segment2[count].end    = segment[count].end   ;
				segment2[count].center = segment[count].center;
				segment2[count].score  = segment[count].score ;
			}
#if 0
			fprintf( stderr, "Goukaku=%dko\n", tmpint ); 
			fprintf( stderr, "in 1 %d\n", segment1[count].center );
			fprintf( stderr, "in 2 %d\n", segment2[count].center );
#endif
			segment1[count].pair = &segment2[count];
			segment2[count].pair = &segment1[count];
			count++;
#if 0
			fprintf( stderr, "count=%d\n", count );
#endif
		}
	}

#if 1
	fprintf( stderr, "done. (%d anchors)\r", count );
#endif
	if( !count && fftNoAnchStop )
		ErrorExit( "Cannot detect anchor!" );
#if 0
	fprintf( stdout, "RESULT before sort:\n" );
	for( l=0; l<count+1; l++ )
	{
		fprintf( stdout, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( stdout, "%d score = %f\n", segment2[l].center, segment1[l].score );
	}
	exit( 1 );
#endif

#if KEIKA
	fprintf( stderr, "Aligning anchors ... " );
#endif
	for( i=0; i<count; i++ )
	{
		sortedseg1[i] = &segment1[i];
		sortedseg2[i] = &segment2[i];
	}

	{
		mymergesort( 0, count-1, sortedseg1 ); 
		mymergesort( 0, count-1, sortedseg2 ); 
		for( i=0; i<count; i++ ) sortedseg1[i]->number = i;
		for( i=0; i<count; i++ ) sortedseg2[i]->number = i;

		if( crossscoresize < count+2 )
		{
			crossscoresize = count+2;
			fprintf( stderr, "####################################################################################################################################allocating crossscore, size = %d\n", crossscoresize );
			if( crossscore ) FreeDoubleMtx( crossscore );
			crossscore = AllocateDoubleMtx( crossscoresize, crossscoresize );
		}

		for( i=0; i<count+2; i++ ) for( j=0; j<count+2; j++ )
			crossscore[i][j] = 0.0;
		for( i=0; i<count; i++ )
		{
			crossscore[segment1[i].number+1][segment1[i].pair->number+1] = segment1[i].score;
			cut1[i+1] = sortedseg1[i]->center;
			cut2[i+1] = sortedseg2[i]->center;
		}

#if DEBUG
		fprintf( stderr, "AFTER SORT\n" );
		for( i=0; i<count; i++ ) fprintf( stderr, "%d, %d\n", segment1[i].start, segment2[i].start );
#endif

		crossscore[0][0] = 10000000.0;
		cut1[0] = 0; 
		cut2[0] = 0;
		crossscore[count+1][count+1] = 10000000.0;
		cut1[count+1] = len1;
		cut2[count+1] = len2;
		count += 2;
		count0 = count;

		blockAlign2( cut1, cut2, sortedseg1, sortedseg2, crossscore, &count );
	}
	if( fftkeika )
	{
		if( count0 > count )
		{
			fprintf( stderr, "REPEAT!? \n" ); 
			if( fftRepeatStop ) exit( 1 );
		}
#if KEIKA
		else 
			fprintf( stderr, "done\n" );
			fprintf( stderr, "done. (%d anchors)\n", count );
#endif
	}

#if 0
	fftfp = fopen( "fft", "a" );
	fprintf( fftfp, "RESULT after sort:\n" );
	for( l=0; l<count; l++ )
	{
		fprintf( fftfp, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( fftfp, "%d\n", segment2[l].center );
	}
	fclose( fftfp );
#endif

#if 0
	fftfp = fopen( "fft", "a" );
	fprintf( fftfp, "RESULT after sort:\n" );
	for( l=0; l<count; l++ )
	{
		fprintf( fftfp, "cut : %d %d\n", cut1[l], cut2[l] );
	}
	fclose( fftfp );
#endif

#if KEIKA
	fprintf( trap_g, "Devided to %d segments\n", count-1 );
	fprintf( trap_g, "%d  %d forg\n", MIN( clus1, clus2 ), count-1 );
#endif

	totallen = 0;
	for( j=0; j<clus1; j++ ) result1[j][0] = 0;
	for( j=0; j<clus2; j++ ) result2[j][0] = 0;
	for( i=0; i<count-1; i++ )
	{
		if( i == 0 ) headgp = outgap; else headgp = 1;
		if( i == count-2 ) tailgp = outgap; else tailgp = 1;
		
#if DEBUG
		fprintf( stderr, "DP %03d / %03d %4d to ", i+1, count-1, totallen );
#else
#if KEIKA
		fprintf( stderr, "DP %03d / %03d\r", i+1, count-1 );
#endif
#endif
		for( j=0; j<clus1; j++ )
		{
			strncpy( tmpres1[j], seq1[j]+cut1[i], cut1[i+1]-cut1[i] );
			tmpres1[j][cut1[i+1]-cut1[i]] = 0;
		}
		for( j=0; j<clus2; j++ )
		{
			strncpy( tmpres2[j], seq2[j]+cut2[i], cut2[i+1]-cut2[i] );
			tmpres2[j][cut2[i+1]-cut2[i]] = 0;
		}
		switch( alg )
		{
			case( 'a' ):
				Aalign( tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen );
				break;
			case( 'M' ):
					MSalignmm( n_dynamicmtx, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, headgp, tailgp, NULL, NULL, NULL, 0.0, 0.0 );
				break;
			case( 'A' ):
				if( clus1 == 1 && clus2 == 1 )
					G__align11( n_dynamicmtx, tmpres1, tmpres2, alloclen, headgp, tailgp );
				else
					A__align( n_dynamicmtx, penalty, penalty_ex, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, headgp, tailgp, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
				break;
			default:
				fprintf( stderr, "alg = %c\n", alg );
				ErrorExit( "ERROR IN SOURCE FILE Falign.c" );
				break;
		}

		nlen = strlen( tmpres1[0] );
		if( totallen + nlen > alloclen ) ErrorExit( "LENGTH OVER in Falign\n " );
		for( j=0; j<clus1; j++ ) strcat( result1[j], tmpres1[j] );
		for( j=0; j<clus2; j++ ) strcat( result2[j], tmpres2[j] );
		totallen += nlen;
#if 0
		fprintf( stderr, "%4d\r", totallen );
		fprintf( stderr, "\n\n" );
		for( j=0; j<clus1; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres1[j] );
		}
		fprintf( stderr, "-------\n" );
		for( j=0; j<clus2; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres2[j] );
		}
#endif
	}
#if KEIKA
	fprintf( stderr, "DP ... done   \n" );
#endif

	for( j=0; j<clus1; j++ ) strcpy( seq1[j], result1[j] );
	for( j=0; j<clus2; j++ ) strcpy( seq2[j], result2[j] );
#if 0
	for( j=0; j<clus1; j++ ) 
	{
		fprintf( stderr, "%s\n", result1[j] );
	}
	fprintf( stderr, "- - - - - - - - - - -\n" );
	for( j=0; j<clus2; j++ ) 
	{
		fprintf( stderr, "%s\n", result2[j] );
	}
#endif
	return( 0.0 );
}



/* work/v8s: MAFFT_V8CHECK=1 (timing builds) -- every parallel pass added by work/v8s (FFT input
   vectors, anchor profiles, row copies) is redone by the original serial code and compared bitwise. */
#ifdef PROGTIMING
int v8check( void )
{
	static int v = -1;
	if( v < 0 ) v = ( getenv( "MAFFT_V8CHECK" ) != NULL );
	return( v );
}
static long v8check_n[4];
static void v8check_report( void ) __attribute__(( destructor ));
static void v8check_report( void )
{
	if( v8check() > 0 ) reporterr( "\nV8CHECK: %ld FFT-vector builds and transform pairs, %ld anchor-profile builds, %ld row-copy passes, %ld in-place Falign segments verified, no mismatch\n", v8check_n[0], v8check_n[1], v8check_n[2], v8check_n[3] );
}
void v8check_count( int k ) { __atomic_add_fetch( v8check_n + k, 1, __ATOMIC_RELAXED ); }
#endif

static int fa_inplace_on( void ) // work/v8s: MAFFT_FAINPLACE=0 (timing builds) turns the in-place single segment off (A/B)
{
#ifdef PROGTIMING
	static int v = -1;
	if( v < 0 ) { char *e = getenv( "MAFFT_FAINPLACE" ); v = !( e && atoi( e ) == 0 ); }
	return( v );
#else
	return( 1 );
#endif
}

/* work/v8s: the two forward transforms of a Falign call are independent; with an idle helper they run
   at the same time.  fft() keeps its tables per thread (computed from n alone, the same on every
   thread), so each transform is the same computation whichever thread runs it. */
typedef struct { int n; Fukusosuu *a, *b; } fft2ctx_t;
static void fft2_body( void *vctx, int lo, int hi )
{
	fft2ctx_t *c = (fft2ctx_t *)vctx;
	int k;
	for( k=lo; k<hi; k++ ) fft( c->n, k ? c->b : c->a, 0 );
}
static int Falign_fft2_par( int n, Fukusosuu *a, Fukusosuu *b )
{
#ifdef enablemultithread
	fft2ctx_t c;
	if( n < 4096 || dphelp_idle() <= 0 ) return( 0 );
	c.n = n; c.a = a; c.b = b;
#ifdef PROGTIMING
	if( v8check() ) // the two transforms again, serially, on copies of the inputs, compared bitwise
	{
		Fukusosuu *ca = malloc( sizeof( Fukusosuu ) * n ), *cb = malloc( sizeof( Fukusosuu ) * n );
		memcpy( ca, a, sizeof( Fukusosuu ) * n ); memcpy( cb, b, sizeof( Fukusosuu ) * n );
		dphelp_for( 2, 1, (double)n, 0.0, fft2_body, &c );
		fft( n, ca, 0 ); fft( n, cb, 0 );
		if( memcmp( ca, a, sizeof( Fukusosuu ) * n ) || memcmp( cb, b, sizeof( Fukusosuu ) * n ) ) { reporterr( "\nV8CHECK MISMATCH in the parallel FFTs (n=%d)\n", n ); exit( 1 ); }
		free( ca ); free( cb );
		v8check_count( 0 );
		return( 1 );
	}
#endif
	dphelp_for( 2, 1, (double)n, 0.0, fft2_body, &c );
	return( 1 );
#else
	return( 0 );
#endif
}

/* work/v8s: the FFT input vectors of a large group, built by column ranges with idle helper threads.
   seq_vec_5 / seq_vec_3 add to entry j, for each member i in ascending order, a term that depends only
   on member i's letter in column j (the same statement here, so the same contraction).  Split by
   columns, every entry receives the same additions in the same order, whichever thread makes them.
   The rows of a group are aligned, all len long, so "up to the terminator" is "up to len".  On the
   critical chain of real rRNA sets (a group of ~1,800 rows x 5,000-8,000 columns plus a few
   sequences, step after step) this pass was ~20% of each step. */
#define FVEC_MINWORK 1.0e6
typedef struct { Fukusosuu **vec; char **seq; double *eff; int clus, nch, five; } fvecctx_t;
static void fvec_cols( void *vctx, int lo, int hi )
{
	fvecctx_t *c = (fvecctx_t *)vctx;
	int i, j, n;
	for( i=0; i<c->clus; i++ )
	{
		const char *s = c->seq[i];
		const double incr = c->eff[i];
		if( c->five ) // = seq_vec_5( vec[0], polarity, volume, incr, s ), columns lo..hi-1
		{
			Fukusosuu *result = c->vec[0];
			for( j=lo; j<hi; j++ )
			{
				n = amino_n[(int)s[j]];
				if( n > 20 ) continue;
				result[j].R += incr * polarity[n];
				result[j].I += incr * volume[n];
			}
		}
		else // = seq_vec_3( vec, incr, s ) with n20or4or2 = nch, columns lo..hi-1
		{
			Fukusosuu **result = c->vec;
			const int nch = c->nch;
			for( j=lo; j<hi; j++ )
			{
				n = amino_n[(int)s[j]];
				if( n < nch && n >= 0 ) result[n][j].R += incr;
			}
		}
	}
}
static void seq_vec_3( Fukusosuu **result, double incr, char *seq );
static void seq_vec_5( Fukusosuu *result, double *score1, double *score2, double incr, char *seq );
/* returns 0 (nothing done) when the caller should run the original loop */
static int Falign_vectors_par( Fukusosuu **vec, char **seq, double *eff, int clus, int len, int nch, int five )
{
#ifdef enablemultithread
	fvecctx_t c;
	if( (double)clus * len < FVEC_MINWORK || dphelp_idle() <= 0 ) return( 0 );
	c.vec = vec; c.seq = seq; c.eff = eff; c.clus = clus; c.nch = nch; c.five = five;
	dphelp_for( len, MAX( 256, len / 64 ), (double)clus * len, FVEC_MINWORK, fvec_cols, &c );
#ifdef PROGTIMING
	if( v8check() ) // the original loop into fresh vectors, compared bitwise (entries 0..len-1; the rest stay 0.0 in both)
	{
		int k, i, nv = five ? 1 : nch;
		Fukusosuu **ref = malloc( sizeof( Fukusosuu * ) * nv );
		for( k=0; k<nv; k++ ) ref[k] = calloc( len + 1, sizeof( Fukusosuu ) );
		for( i=0; i<clus; i++ )
		{
			if( five ) seq_vec_5( ref[0], polarity, volume, eff[i], seq[i] );
			else seq_vec_3( ref, eff[i], seq[i] );
		}
		for( k=0; k<nv; k++ )
		{
			if( memcmp( ref[k], vec[k], sizeof( Fukusosuu ) * len ) ) { reporterr( "\nV8CHECK MISMATCH in the FFT vectors (channel %d, clus=%d, len=%d)\n", k, clus, len ); exit( 1 ); }
			free( ref[k] );
		}
		free( ref );
		v8check_count( 0 );
	}
#endif
	return( 1 );
#else
	return( 0 );
#endif
}

/* work/v8s: Falign's row copies (segment extraction, appending, the final copy back) for large groups,
   row by row with idle helper threads.  Each row gets exactly the bytes the serial loop gives it (the
   same strncpy / memcpy / strcpy call), only possibly on another thread; rows are independent. */
#define RCPY_MINWORK 2.0e6
typedef struct { char **dst, **src; int off, n, mode; } rcpyctx_t;
static void rcpy_rows( void *vctx, int lo, int hi )
{
	rcpyctx_t *c = (rcpyctx_t *)vctx;
	int j;
	if( c->mode == 0 ) for( j=lo; j<hi; j++ ) { strncpy( c->dst[j], c->src[j]+c->off, c->n ); c->dst[j][c->n] = 0; }
	else if( c->mode == 1 ) for( j=lo; j<hi; j++ ) memcpy( c->dst[j]+c->off, c->src[j], c->n+1 );
	else for( j=lo; j<hi; j++ ) strcpy( c->dst[j], c->src[j] );
}
/* mode 0: dst[j] = src[j][off..off+n-1] + terminator; 1: src[j] (n chars + terminator) appended at dst[j]+off;
   2: strcpy( dst[j], src[j] ) with rows of n chars.  Returns 0 (nothing done) when the caller should
   run its own loop. */
static int rows_copy_par( char **dst, char **src, int clus, int off, int n, int mode )
{
#ifdef enablemultithread
	rcpyctx_t c;
	if( (double)clus * n < RCPY_MINWORK || dphelp_idle() <= 0 ) return( 0 );
	c.dst = dst; c.src = src; c.off = off; c.n = n; c.mode = mode;
	dphelp_for( clus, MAX( 8, clus / 64 ), (double)clus * n, RCPY_MINWORK, rcpy_rows, &c );
#ifdef PROGTIMING
	if( v8check() )
	{
		int j;
		for( j=0; j<clus; j++ )
		{
			int bad = 0;
			if( mode == 0 ) bad = memcmp( dst[j], src[j]+off, n ) || dst[j][n] != 0;
			else if( mode == 1 ) bad = memcmp( dst[j]+off, src[j], n+1 );
			else bad = strcmp( dst[j], src[j] );
			if( bad ) { reporterr( "\nV8CHECK MISMATCH in a row copy (mode %d, row %d)\n", mode, j ); exit( 1 ); }
		}
		v8check_count( 2 );
	}
#endif
	return( 1 );
#else
	return( 0 );
#endif
}

/* work/prog2: per-thread arena for Falign's per-call row buffers (result1/2 and tmpres1/2, clus x
   alloclen each, used to be calloc'ed row by row and freed on every call).  Every row is written
   with a terminator before it is read (tmpres by strncpy + explicit 0, result by result[j][0] = 0
   and terminated appends); the first byte of each row is still set to 0. */
static TLS char *fa_blk = NULL;
static TLS size_t fa_cap = 0;
static TLS char **fa_rows = NULL;
static TLS int fa_rowcap = 0;

static void Falign_rows( int clus1, int clus2, int alloclen, char ***r1, char ***r2, char ***t1, char ***t2 )
{
	int n = 2 * ( clus1 + clus2 ), i;
	size_t need = (size_t)n * alloclen;
	char **p;
	if( need > fa_cap )
	{
#ifdef PROGTIMING
		pt_acc[43] += 1.0; // work/v8s: arena reallocations
#endif
		/* work/v8s: realloc instead of free + malloc (the content does not matter): for a large
		   (mmapped) block glibc moves the pages with mremap, so the pages already touched need not be
		   faulted in again along the critical chain, where the arena grows every few steps */
		fa_cap = need + need / 2;
		{
			char *nb = realloc( fa_blk, fa_cap );
			if( nb == NULL ) ErrorExit( "Cannot allocate the work area in Falign.\n" );
			fa_blk = nb;
		}
	}
	if( n + 4 > fa_rowcap )
	{
		free( fa_rows );
		fa_rowcap = n + 4 + ( n + 4 ) / 2;
		fa_rows = malloc( fa_rowcap * sizeof( char * ) );
		if( fa_rows == NULL ) ErrorExit( "Cannot allocate the work area in Falign.\n" );
	}
	p = fa_rows;
	*r1 = p; for( i=0; i<clus1; i++ ) *p++ = fa_blk + (size_t)( i ) * alloclen; *p++ = NULL;
	*r2 = p; for( i=0; i<clus2; i++ ) *p++ = fa_blk + (size_t)( clus1 + i ) * alloclen; *p++ = NULL;
	*t1 = p; for( i=0; i<clus1; i++ ) *p++ = fa_blk + (size_t)( clus1 + clus2 + i ) * alloclen; *p++ = NULL;
	*t2 = p; for( i=0; i<clus2; i++ ) *p++ = fa_blk + (size_t)( 2*clus1 + clus2 + i ) * alloclen; *p++ = NULL;
	for( i=0; i<n; i++ ) fa_blk[(size_t)i*alloclen] = 0;
}

static void Falign_rows_free( void )
{
	free( fa_blk ); fa_blk = NULL; fa_cap = 0;
	free( fa_rows ); fa_rows = NULL; fa_rowcap = 0;
}

void Falign_trim( size_t keep ) // work/largen: give a large row arena back (see A__align_trim)
{
	if( fa_cap > keep ) { free( fa_blk ); fa_blk = NULL; fa_cap = 0; }
}

/* work/prog3: the segments of a Falign call in parallel.  After the anchors are fixed, every
   segment between two cut points is aligned on its own: its inputs are the original columns
   cut1[i]..cut1[i+1]-1 (and cut2[...]) of seq1/seq2 plus the gap states at the cut points
   (getkyokaigap reads the original rows), and its output is appended to the result in segment order.
   So segments can be aligned by any thread in any order, each into its own buffers, and then
   appended in order: every row of the result gets the same bytes as in the sequential loop.  The
   aligners keep their working state per thread (TLS), the totalscore sum is formed in segment order.
   Used for alg 'A' without constraint/codon options, when idle helper threads exist. */
typedef struct
{
	char **seq1, **seq2;
	int clus1, clus2, nseg, alloclen;
	const int *cut1, *cut2;
	double **n_dynamicmtx;
	double *eff1, *eff2;
	char ***r1, ***r2;      // per segment: the rows it aligns in place
	char **sg;              // per segment: sgap1, sgap2, egap1, egap2 (4 consecutive blocks)
	double *score;
} psegctx_t;

static void pseg_body( void *vctx, int lo, int hi )
{
	psegctx_t *c = (psegctx_t *)vctx;
	int i, j;
	for( i=lo; i<hi; i++ )
	{
		const int *cut1 = c->cut1, *cut2 = c->cut2;
		int clus1 = c->clus1, clus2 = c->clus2, headgp, tailgp;
		char **tmpres1 = c->r1[i], **tmpres2 = c->r2[i];
		char *sgap1 = c->sg[i], *sgap2 = sgap1 + clus1, *egap1 = sgap2 + clus2, *egap2 = egap1 + clus1;
		int len1 = cut1[c->nseg], len2 = cut2[c->nseg];
		double dumdb = 0.0;
		if( i == 0 ) headgp = outgap; else headgp = 1;
		if( i == c->nseg-1 ) tailgp = outgap; else tailgp = 1;
		if( cut1[i] )
			getkyokaigap( sgap1, c->seq1, cut1[i]-1, clus1 );
		else
			for( j=0; j<clus1; j++ ) sgap1[j] = 'o';
		if( cut2[i] )
			getkyokaigap( sgap2, c->seq2, cut2[i]-1, clus2 );
		else
			for( j=0; j<clus2; j++ ) sgap2[j] = 'o';
		if( cut1[i+1] != len1 )
			getkyokaigap( egap1, c->seq1, cut1[i+1], clus1 );
		else
			for( j=0; j<clus1; j++ ) egap1[j] = 'o';
		if( cut2[i+1] != len2 )
			getkyokaigap( egap2, c->seq2, cut2[i+1], clus2 );
		else
			for( j=0; j<clus2; j++ ) egap2[j] = 'o';
		for( j=0; j<clus1; j++ )
		{
			strncpy( tmpres1[j], c->seq1[j]+cut1[i], cut1[i+1]-cut1[i] );
			tmpres1[j][cut1[i+1]-cut1[i]] = 0;
		}
		for( j=0; j<clus2; j++ )
		{
			strncpy( tmpres2[j], c->seq2[j]+cut2[i], cut2[i+1]-cut2[i] );
			tmpres2[j][cut2[i+1]-cut2[i]] = 0;
		}
		if( clus1 == 1 && clus2 == 1 )
			c->score[i] = G__align11( c->n_dynamicmtx, tmpres1, tmpres2, c->alloclen, headgp, tailgp );
		else
			c->score[i] = A__align( c->n_dynamicmtx, penalty, penalty_ex, tmpres1, tmpres2, c->eff1, c->eff2, clus1, clus2, c->alloclen, 0, &dumdb, sgap1, sgap2, egap1, egap2, NULL, 0, NULL, headgp, tailgp, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
	}
}

#ifdef PROGTIMING
static TLS int pseg_check_pending = 0;
static TLS char **pseg_check_save1, **pseg_check_save2;
static TLS double pseg_check_ts;
static long pseg_checked = 0;
static void pseg_report( void ) __attribute__(( destructor ));
static void pseg_report( void ) { if( getenv( "MAFFT_PSEGCHECK" ) ) reporterr( "\nPSEGCHECK: %ld parallel Falign calls verified, no mismatch\n", pseg_checked ); }
#endif
static int pseg_force( void ) // MAFFT_PSEGFORCE=1 (timing builds): every multi-segment call takes the parallel code (testing)
{
#ifdef PROGTIMING
	static int v = -1;
	if( v < 0 ) v = ( getenv( "MAFFT_PSEGFORCE" ) != NULL );
	return( v );
#else
	return( 0 );
#endif
}
static int pseg_v = -1;
static int pseg_on( void )
{
	if( pseg_v < 0 )
	{
		/* work/prog3: off by default.  The segments' private buffers hold a second copy of both groups
		   while the segments run, on several threads at once: +0.7-0.8 GB peak RSS on 16S.T at 10
		   threads (above the baseline's), for a gain that did not show clearly.  MAFFT_PSEG=1 turns it
		   on in timing builds (A/B, PSEGCHECK). */
#ifdef PROGTIMING
		pseg_v = ( getenv( "MAFFT_PSEG" ) != NULL );
#else
		pseg_v = 0;
#endif
	}
	return( pseg_v );
}

/* Aligns the nseg segments (cut points cut1[0..nseg], cut2[0..nseg]) and appends them, in order, to
   result1/result2 (which start empty); returns totalscore, as the sequential loop does. */
static double Falign_psegs( char **seq1, char **seq2, int clus1, int clus2, int nseg, const int *cut1, const int *cut2,
	double **n_dynamicmtx, double *eff1, double *eff2, int alloclen, char **result1, char **result2 )
{
	psegctx_t c;
	size_t rowbytes = 0, k;
	char *blk, *p;
	int i, j, totallen = 0, off1 = 0, off2 = 0;
	double totalscore = 0.0;
	for( i=0; i<nseg; i++ ) rowbytes += ( cut1[i+1]-cut1[i] ) + ( cut2[i+1]-cut2[i] ) + 2;
	blk = malloc( rowbytes * ( clus1 + clus2 ) + (size_t)nseg * 2 * ( clus1 + clus2 ) );
	c.r1 = malloc( sizeof( char ** ) * nseg ); c.r2 = malloc( sizeof( char ** ) * nseg );
	c.sg = malloc( sizeof( char * ) * nseg ); c.score = malloc( sizeof( double ) * nseg );
	if( !blk || !c.r1 || !c.r2 || !c.sg || !c.score ) ErrorExit( "Cannot allocate the segment buffers in Falign.\n" );
	p = blk;
	for( i=0; i<nseg; i++ )
	{
		k = ( cut1[i+1]-cut1[i] ) + ( cut2[i+1]-cut2[i] ) + 2;
		c.r1[i] = malloc( sizeof( char * ) * clus1 ); c.r2[i] = malloc( sizeof( char * ) * clus2 );
		for( j=0; j<clus1; j++ ) { c.r1[i][j] = p; p += k; }
		for( j=0; j<clus2; j++ ) { c.r2[i][j] = p; p += k; }
		c.sg[i] = p; p += 2 * ( clus1 + clus2 );
	}
	c.seq1 = seq1; c.seq2 = seq2; c.clus1 = clus1; c.clus2 = clus2; c.nseg = nseg; c.alloclen = alloclen;
	c.cut1 = cut1; c.cut2 = cut2; c.n_dynamicmtx = n_dynamicmtx; c.eff1 = eff1; c.eff2 = eff2;
	dphelp_for( nseg, 1, 1.0, 0.0, pseg_body, &c ); // one segment per claim, any thread, any order
	for( i=0; i<nseg; i++ ) // append in order, exactly as the sequential loop
	{
		int nlen = strlen( c.r1[i][0] ), add2 = strlen( c.r2[i][0] );
		totalscore += c.score[i];
		if( totallen + nlen > alloclen )
		{
			fprintf( stderr, "totallen=%d +  nlen=%d > alloclen = %d\n", totallen, nlen, alloclen );
			ErrorExit( "LENGTH OVER in Falign\n " );
		}
		for( j=0; j<clus1; j++ ) memcpy( result1[j]+off1, c.r1[i][j], nlen+1 );
		for( j=0; j<clus2; j++ ) memcpy( result2[j]+off2, c.r2[i][j], add2+1 );
		off1 += nlen;
		off2 += add2;
		totallen += nlen;
	}
	for( i=0; i<nseg; i++ ) { free( c.r1[i] ); free( c.r2[i] ); }
	free( c.r1 ); free( c.r2 ); free( c.sg ); free( c.score ); free( blk );
	return( totalscore );
}

double Falign( int **whichmtx, double ***scoringmatrices, double **n_dynamicmtx,
			  char  **seq1, char  **seq2, 
			  double *eff1, double *eff2, 
			  double **eff1s, double **eff2s,
			  int    clus1, int    clus2,
			  int alloclen, int *fftlog,
			  int *chudanpt, int chudanref, int *chudanres )
{
	int i, j, k, l, m, maxk;
	int nlen, nlen2, nlen4;
	static TLS int crossscoresize = 0;
	char **tmpseq1 = NULL;
	char **tmpseq2 = NULL;
	char **tmpptr1 = NULL;
	char **tmpptr2 = NULL;
	char **tmpres1 = NULL;
	char **tmpres2 = NULL;
	char **result1 = NULL;
	char **result2 = NULL;
#if RND
	char **rndseq1 = NULL;
	char **rndseq2 = NULL;
#endif
	static TLS Fukusosuu **seqVector1 = NULL;
	static TLS Fukusosuu **seqVector2 = NULL;
	static TLS Fukusosuu **naiseki = NULL;   
	static TLS Fukusosuu *naisekiNoWa = NULL; 
	static TLS double *soukan = NULL;
	static TLS double **crossscore = NULL;
	int nlentmp;
	static TLS int *kouho = NULL;
	static TLS Segment *segment = NULL;
	static TLS Segment *segment1 = NULL;
	static TLS Segment *segment2 = NULL;
	static TLS Segment **sortedseg1 = NULL;
	static TLS Segment **sortedseg2 = NULL;
	static TLS int *cut1 = NULL;
	static TLS int *cut2 = NULL;
	char *sgap1, *egap1, *sgap2, *egap2;
	static TLS int localalloclen = 0;
	int lag;
	int tmpint;
	int count, count0;
	int len1, len2;
	int totallen;
	int off1, off2; // work/prog: current length of the rows of result1 / result2
	int inpl; // work/v8s: one segment, aligned in place
	double totalscore;
	double dumdb = 0.0;
	int headgp, tailgp;
	static TLS double *gstart = NULL;
	static TLS double *gend = NULL;
	static TLS double **codonscoremtx;
	PT_DECL;


	if( seq1 == NULL )
	{
		if( kouho ) 
		{
//			fprintf( stderr, "Freeing localarrays in Falign\n" );
			localalloclen = 0;
			crossscoresize = 0;
			mymergesort( 0, 0, NULL );
			alignableReagion( 0, 0, NULL, NULL, NULL, NULL, NULL );
			alignableReagion_prep_free();
			fft( 0, NULL, 1 );
			A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
			D__align( NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
			A__align_variousdist( NULL, NULL, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
			D__align_variousdist( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
			G__align11( NULL, NULL, NULL, 0, 0, 0 );
			G__align11psg( NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL );
			blockAlign2( NULL, NULL, NULL, NULL, NULL, NULL );
			if( crossscore ) FreeDoubleMtx( crossscore );
			crossscore = NULL;
			free( kouho );
			kouho = NULL;
			free( cut1 );
			free( cut2 );
			free( segment );
			free( segment1 );
			free( segment2 );
			free( sortedseg1 );
			free( sortedseg2 );
			if( gstart ) free( gstart );
			if( gend ) free( gend );
			if( codonscoremtx ) FreeDoubleMtx( codonscoremtx );
			Falign_rows_free();
			if( !kobetsubunkatsu )
			{
				FreeFukusosuuMtx ( seqVector1 );
				FreeFukusosuuMtx ( seqVector2 );
				FreeFukusosuuVec( naisekiNoWa );
				FreeFukusosuuMtx( naiseki );
				FreeDoubleVec( soukan );
			}
		}
		else
		{
//			fprintf( stderr, "Did not allocate localarrays in Falign\n" );
		}

		return( 0.0 );
	}


	PT_START();
	len1 = strlen( seq1[0] );
	len2 = strlen( seq2[0] );
	nlentmp = MAX( len1, len2 );

	nlen = 1;
	while( nlentmp >= nlen ) nlen <<= 1;
#if 0
	fprintf( stderr, "###   nlen    = %d\n", nlen );
#endif

	nlen2 = nlen/2; nlen4 = nlen2 / 2;

#if DEBUG
	fprintf( stderr, "len1 = %d, len2 = %d\n", len1, len2 );
	fprintf( stderr, "nlentmp = %d, nlen = %d\n", nlentmp, nlen );
#endif


	/* work/largen: rows of len1 + len2 + 100 bytes (at most alloclen) instead of alloclen.  Each row
	   here holds a merged row of this step or an aligned pair of segments of it, at most len1 + len2
	   characters plus the terminator; alloclen is the longest alignment so far in the whole pass, so
	   below the root the rows were mostly unused (and, with threads, several steps' arenas exist at
	   once). */
	/* work/v8s: the arena is set up only when the segment loop needs it (below; not for the in-place
	   single segment), after the anchor search */
	sgap1 = AllocateCharVec( clus1 );
	egap1 = AllocateCharVec( clus1 );
	sgap2 = AllocateCharVec( clus2 );
	egap2 = AllocateCharVec( clus2 );
	tmpptr1 = calloc( clus1, sizeof( char * ) );
	tmpptr2 = calloc( clus2, sizeof( char * ) );
#if RND
	tmpseq1 = AllocateCharMtx( clus1, nlen );
	tmpseq2 = AllocateCharMtx( clus2, nlen );
#endif
#if RND
	rndseq1 = AllocateCharMtx( clus1, nlen );
	rndseq2 = AllocateCharMtx( clus2, nlen );
	for( i=0; i<clus1; i++ )
		generateRndSeq( rndseq1[i], nlen );
	for( i=0; i<clus2; i++ )
		generateRndSeq( rndseq2[i], nlen );
#endif


	if( !localalloclen )
	{
		kouho = AllocateIntVec( NKOUHO );
		cut1 = AllocateIntVec( MAXSEG );
		cut2 = AllocateIntVec( MAXSEG );
//		crossscore = AllocateDoubleMtx( MAXSEG, MAXSEG );
		segment = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		segment1 = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		segment2 = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		sortedseg1 = (Segment **)calloc( MAXSEG, sizeof( Segment * ) );
		sortedseg2 = (Segment **)calloc( MAXSEG, sizeof( Segment * ) );
		if( !( segment && segment1 && segment2 && sortedseg1 && sortedseg2 ) )
			ErrorExit( "Allocation error\n" );

		if     ( scoremtx == -1 ) n20or4or2 = 1;
		else if( fftscore )       n20or4or2 = 1;
		else                      n20or4or2 = 20;


		gstart = NULL;
		gend = NULL;
	 	if( codonpos ) // gstart to gend ha sudeni allocate sareteiru kamo
		{
			FILE *cfp;
			int p;
			char *buf = calloc( sizeof( char ), 1000 );

			if( dorp != 'd' )
			{
				reporterr( "\n\nThe --codonpos and --codonscore options are only for DNA data.\n\n" );
				exit( 1 );
			}

//			reporterr( "\nIn Falign, loading position-specific gap costs\n" );
			gstart = calloc( sizeof(double), len1+1 );
			for( i=0; i<len1+1; i++ ) gstart[i] = 1.0 * 0.5; // init
			gend = calloc( sizeof(double), len1+1 );
			for( i=0; i<len1+1; i++ ) gend[i] = 1.0 * 0.5; // init
#define STRONG 3.0

			i = 0;
			cfp = fopen( "_codonpos", "r" );	
			if( cfp == NULL )
			{
				reporterr( "Cannot open _codonpos file\n" );
				exit( 1 );
			}
			while( fgets( buf, 1000, cfp ) )
			{
				if( i == len1 )
				{
					reporterr( "\n\nNumber of lines in the codonposition file must be the same as the length of the first sequences (%d).\n\n", len1 );
					exit( 1 );
				}

				if( buf[0] == '#' )
				{
					continue;
				}
				else if( buf[0] == '0' )
				{
					gstart[i] = 1.0 * 0.5; // noncoding
					gend[i] = 1.0 * 0.5; // noncoding
				}
				else if( buf[0] == '1' )
				{
					gstart[i] = STRONG * 0.5;
					gend[i] = 1.0 * 0.5;
				}
				else if( buf[0] == '2' )
				{
					gstart[i] = STRONG * 0.5;
					gend[i] = STRONG * 0.5;
				}
				else if( buf[0] == '3' )
				{
					gstart[i] = 1.0 * 0.5;
					gend[i] = STRONG * 0.5;
				}
				else
				{
					reporterr( "In the codonposition file, 1st letter in a line must be either of:\n" );
					reporterr( "	0 (noncoding)\n" );
					reporterr( "	1 (1st position in codon)\n" );
					reporterr( "	2 (2nd position in codon)\n" );
					reporterr( "	3 (3rd position in codon)\n" );
					reporterr( "When mutliple difference frames are used, set 2.\n" );
					exit( 1 );
				}
				i++;
			}
			fclose( cfp );
			free( buf );
			if( i < len1 )
			{
				reporterr( "\n\nNumber of lines in the codonposition file (%d) is less than the length of the first sequences (%d).\n\n", i, len1 );
				exit( 1 );
			}

		}

		if( codonscore )
		{
			int i, j;
			FILE *cfp;
			double *codonfreq = calloc( sizeof( double ), 64 );
			double totalcount, codonscore0av, codonscore1av;
			if( !codonpos )
			{
				reporterr( "\n\n --codonpos is necessary for --codonscore\n\n" );
				exit( 1 );
			}
			codonscoremtx = AllocateDoubleMtx( 64, 64 );
			cfp = fopen( "_codonscore", "r" );	
			loadcodonscore( cfp, codonscoremtx );
			fclose( cfp );
			for( i=0; i<64; i++ ) codonfreq[i] = 0.0;
			totalcount = 0.0;
			for( i=3; i<len1; i++ )
			{
//				reporterr( "i=%d\n", i );
				if( gstart[i-2] == STRONG * 0.5 && gend[i-2] == 1.0 * 0.5    && 
				    gstart[i-1] == STRONG * 0.5 && gend[i-1] == STRONG * 0.5 &&
				    gstart[i-0] == 1.0 * 0.5    && gend[i-0] == STRONG * 0.5 )
				{
//					reporterr( "codon=%.3s, id=%d\n", seq1[0]+i-2, codon2id(seq1[0]+i-2) );
					codonfreq[codon2id(seq1[0]+i-2)] += 1.0;
					totalcount += 1.0;
				}
			}
			for( i=0; i<64; i++ ) codonfreq[i] /= totalcount;
//			for( i=0; i<64; i++ ) reporterr( "%d, %f\n", i, codonfreq[i] );

			codonscore0av = 0.0;
			for( i=0; i<64; i++ ) for( j=0; j<64; j++ ) codonscore0av += codonfreq[i]*codonfreq[j] * codonscoremtx[i][j];
//			reporterr( "0av=%f\n", codonscore0av );
			for( i=0; i<64; i++ ) for( j=0; j<64; j++ ) codonscoremtx[i][j] -= codonscore0av;

			codonscore1av = 0.0;
			for( i=0; i<64; i++ ) codonscore1av += codonfreq[i] * codonscoremtx[i][i];
//			reporterr( "1av=%f\n", codonscore1av );
			for( i=0; i<64; i++ ) for( j=0; j<64; j++ ) codonscoremtx[i][j] /= codonscore1av;

#if 0
			codonscore1av = 0.0;
			for( i=0; i<64; i++ ) codonscore1av += codonfreq[i] * codonscoremtx[i][i];
			reporterr( "1av=%f\n", codonscore1av );

			codonscore0av = 0.0;
			for( i=0; i<64; i++ ) for( j=0; j<64; j++ ) codonscore0av += codonfreq[i]*codonfreq[j] * codonscoremtx[i][j];
			reporterr( "0av=%f\n", codonscore0av );
#endif

			free( codonfreq );
		}
		else
			codonscoremtx = NULL;
	}

	if( localalloclen < nlen )
	{
		if( localalloclen )
		{
#if 1
			if( !kobetsubunkatsu )
			{
				FreeFukusosuuMtx ( seqVector1 );
				FreeFukusosuuMtx ( seqVector2 );
				FreeFukusosuuVec( naisekiNoWa );
				FreeFukusosuuMtx( naiseki );
				FreeDoubleVec( soukan );
			}
#endif
		}

		if( !kobetsubunkatsu )
		{
			naisekiNoWa = AllocateFukusosuuVec( nlen );
			naiseki = AllocateFukusosuuMtx( n20or4or2, nlen );
			seqVector1 = AllocateFukusosuuMtx( n20or4or2+1, nlen+1 );
			seqVector2 = AllocateFukusosuuMtx( n20or4or2+1, nlen+1 );
			soukan = AllocateDoubleVec( nlen+1 );
		}
		localalloclen = nlen;
	}

#if RND
	for( j=0; j<clus1; j++ ) strcpy( tmpseq1[j], seq1[j] );
	for( j=0; j<clus2; j++ ) strcpy( tmpseq2[j], seq2[j] );
#else
	/* work/prog2: the FFT input vectors are built from seq1/seq2 directly; tmpseq1/2 were
	   byte-for-byte copies that seq_vec_3/seq_vec_5 only read up to the terminator. */
	tmpseq1 = seq1;
	tmpseq2 = seq2;
#endif
	PT_ADD(0);

#if 0
fftfp = fopen( "input_of_Falign", "w" );
fprintf( fftfp, "nlen = %d\n", nlen );
fprintf( fftfp, "seq1: ( %d sequences ) \n", clus1 );
for( i=0; i<clus1; i++ )
	fprintf( fftfp, "%s\n", seq1[i] );
fprintf( fftfp, "seq2: ( %d sequences ) \n", clus2 );
for( i=0; i<clus2; i++ )
	fprintf( fftfp, "%s\n", seq2[i] );
fclose( fftfp );
system( "less input_of_Falign < /dev/tty > /dev/tty" );
#endif
	if( !kobetsubunkatsu )
	{
		if( fftkeika ) fprintf( stderr,  " FFT ... " );
#ifdef PROGTIMING
		double _v0 = pt_now();
#endif

		for( j=0; j<n20or4or2; j++ ) vec_init( seqVector1[j], nlen );
		if( Falign_vectors_par( seqVector1, tmpseq1, eff1, clus1, len1, n20or4or2, fftscore && scoremtx != -1 ) ) // work/v8s
			;
		else if( fftscore && scoremtx != -1 )
		{
			for( i=0; i<clus1; i++ )
			{
#if 1
				seq_vec_5( seqVector1[0], polarity, volume, eff1[i], tmpseq1[i] );
#else
				seq_vec_2( seqVector1[0], polarity, eff1[i], tmpseq1[i] );
				seq_vec_2( seqVector1[1], volume,   eff1[i], tmpseq1[i] );
#endif
			}
		}
		else
		{
#if 0
			for( i=0; i<clus1; i++ ) for( j=0; j<n20or4or2; j++ ) 
				seq_vec( seqVector1[j], amino[j], eff1[i], tmpseq1[i] );
#else
			for( i=0; i<clus1; i++ )
				seq_vec_3( seqVector1, eff1[i], tmpseq1[i] );
#endif
		}
#if RND
		for( i=0; i<clus1; i++ )
		{
			vec_init2( seqVector1, rndseq1[i], eff1[i], len1, nlen );
		}
#endif
#if 0
fftfp = fopen( "seqVec", "w" );
fprintf( fftfp, "before transform\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "nlen=%d\n", nlen );
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector1[k][l].R, seqVector1[k][l].I );
}
fclose( fftfp );
system( "less seqVec < /dev/tty > /dev/tty" );
#endif

		for( j=0; j<n20or4or2; j++ ) vec_init( seqVector2[j], nlen );
		if( Falign_vectors_par( seqVector2, tmpseq2, eff2, clus2, len2, n20or4or2, fftscore && scoremtx != -1 ) ) // work/v8s
			;
		else if( fftscore && scoremtx != -1 )
		{
			for( i=0; i<clus2; i++ )
			{
#if 1
				seq_vec_5( seqVector2[0], polarity, volume, eff2[i], tmpseq2[i] );
#else
				seq_vec_2( seqVector2[0], polarity, eff2[i], tmpseq2[i] );
				seq_vec_2( seqVector2[1], volume,   eff2[i], tmpseq2[i] );
#endif
			}
		}
		else
		{
#if 0
			for( i=0; i<clus2; i++ ) for( j=0; j<n20or4or2; j++ ) 
				seq_vec( seqVector2[j], amino[j], eff2[i], tmpseq2[i] );
#else
			for( i=0; i<clus2; i++ )
				seq_vec_3( seqVector2, eff2[i], tmpseq2[i] );
#endif
		}
#if RND
		for( i=0; i<clus2; i++ )
		{
			vec_init2( seqVector2, rndseq2[i], eff2[i], len2, nlen );
		}
#endif

#if 0
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "before fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

#ifdef PROGTIMING
		pt_acc[39] += pt_now() - _v0; // work/v8s: FFT input vectors
#endif
		for( j=0; j<n20or4or2; j++ )
		{
			if( Falign_fft2_par( nlen, seqVector2[j], seqVector1[j] ) ) continue; // work/v8s: the two transforms at once
			fft( nlen, seqVector2[j], 0 );
			fft( nlen, seqVector1[j], 0 );
		}
#if 0
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "#after fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "#%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
	   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

		for( k=0; k<n20or4or2; k++ ) 
		{
			for( l=0; l<nlen; l++ ) 
				calcNaiseki( naiseki[k]+l, seqVector1[k]+l, seqVector2[k]+l );
		}
		for( l=0; l<nlen; l++ ) 
		{
			naisekiNoWa[l].R = 0.0;
			naisekiNoWa[l].I = 0.0;
			for( k=0; k<n20or4or2; k++ ) 
			{
				naisekiNoWa[l].R += naiseki[k][l].R;
				naisekiNoWa[l].I += naiseki[k][l].I;
			}
		}
	
#if 0
	fftfp = fopen( "naisekiNoWa", "w" );
	fprintf( fftfp, "#Before fft\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f %f\n", l, naisekiNoWa[l].R, naisekiNoWa[l].I ); 
	fclose( fftfp );
	system( "less naisekiNoWa < /dev/tty > /dev/tty " );
#endif

		fft( -nlen, naisekiNoWa, 0 );
	
		for( m=0; m<=nlen2; m++ ) 
			soukan[m] = naisekiNoWa[nlen2-m].R;
		for( m=nlen2+1; m<nlen; m++ ) 
			soukan[m] = naisekiNoWa[nlen+nlen2-m].R;

#if 0
	fftfp = fopen( "naisekiNoWa", "w" );
	fprintf( fftfp, "#After fft\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f\n", l, naisekiNoWa[l].R ); 
	fclose( fftfp );
	fftfp = fopen( "list.plot", "w"  );
	fprintf( fftfp, "plot 'naisekiNoWa'\npause -1" );
	fclose( fftfp );
	system( "/usr/bin/gnuplot list.plot &" );
#endif
#if 0
	fprintf( stderr, "soukan\n" );
	for( l=0; l<nlen; l++ )
		fprintf( stderr, "%d  %f\n", l-nlen2, soukan[l] ); 
#if 0
	fftfp = fopen( "list.plot", "w"  );
	fprintf( fftfp, "plot 'frt'\n pause +1" );
	fclose( fftfp );
	system( "/usr/bin/gnuplot list.plot" );
#endif
#endif


		getKouho( kouho, NKOUHO, soukan, nlen );
		PT_ADD(1);

#if 0
		for( i=0; i<NKOUHO; i++ )
		{
			fprintf( stderr, "kouho[%d] = %d\n", i, kouho[i] );
		}
#endif
	}

#if KEIKA
	fprintf( stderr, "Searching anchors ... " );
#endif
	count = 0;




#define CAND 0
#if CAND
	fftfp = fopen( "cand", "w" );
	fclose( fftfp );
#endif
	if( kobetsubunkatsu )
	{
		maxk = 1;
		kouho[0] = 0;
	}
	else
	{
		maxk = NKOUHO;
	}

	/* work/prog2: the column profiles for alignableReagion are built once for all candidate lags
	   (alignableReagion_prep / alignableReagion_lag in fftFunctions.c) */
#ifdef PROGTIMING
	{ double _a0 = pt_now();
#endif
#ifdef enablemultithread
#ifdef PROGTIMING
	{ int _nc = alignableReagion_prep( clus1, clus2, seq1, seq2, eff1, eff2, dphelp_for, dphelp_idle, v8check() ); // work/v8s: parallel for large groups
	  while( _nc-- > 0 ) v8check_count( 1 ); }
#else
	alignableReagion_prep( clus1, clus2, seq1, seq2, eff1, eff2, dphelp_for, dphelp_idle, 0 ); // work/v8s: parallel for large groups
#endif
#else
	alignableReagion_prep( clus1, clus2, seq1, seq2, eff1, eff2, NULL, NULL, 0 );
#endif
#ifdef PROGTIMING
	pt_acc[36] += pt_now() - _a0; } // work/v8s: anchors, profiles
#endif
	for( k=0; k<maxk; k++ ) 
	{
		lag = kouho[k];
		if( lag <= -len1 || len2 <= lag ) continue;
		zurasu2( lag, clus1, clus2, seq1, seq2, tmpptr1, tmpptr2 );
#if CAND
		fftfp = fopen( "cand", "a" );
		fprintf( fftfp, ">Candidate No.%d lag = %d\n", k+1, lag );
		fprintf( fftfp, "%s\n", tmpptr1[0] );
		fprintf( fftfp, ">Candidate No.%d lag = %d\n", k+1, lag );
		fprintf( fftfp, "%s\n", tmpptr2[0] );
		fprintf( fftfp, ">\n", k+1, lag );
		fclose( fftfp );
#endif

//		fprintf( stderr, "lag = %d\n", lag );
#ifdef PROGTIMING
		{ double _a0 = pt_now();
#endif
		tmpint = alignableReagion_lag( lag > 0 ? 0 : -lag, lag > 0 ? lag : 0, segment+count ); // = alignableReagion( clus1, clus2, tmpptr1, tmpptr2, eff1, eff2, segment+count )
#ifdef PROGTIMING
		pt_acc[37] += pt_now() - _a0; pt_acc[38] += 1.0; } // work/v8s: anchors, per-lag scoring and scan
#endif
#ifdef PROGTIMING
		if( getenv( "MAFFT_ARCHECK" ) ) // work/prog2: compare with the original, segment by segment, bitwise
		{
			static TLS Segment *chkseg = NULL;
			int t, ref;
			if( chkseg == NULL ) chkseg = calloc( MAXSEG, sizeof( Segment ) );
			for( t=0; t<tmpint+2 && count+t<MAXSEG; t++ ) chkseg[t] = segment[count+t];
			ref = alignableReagion( clus1, clus2, tmpptr1, tmpptr2, eff1, eff2, segment+count );
			if( ref != tmpint ) { reporterr( "\nARCHECK MISMATCH: %d vs %d segments\n", tmpint, ref ); exit( 1 ); }
			for( t=0; t<tmpint+2 && count+t<MAXSEG; t++ )
				if( chkseg[t].start != segment[count+t].start || chkseg[t].end != segment[count+t].end || chkseg[t].center != segment[count+t].center
				 || memcmp( &chkseg[t].score, &segment[count+t].score, sizeof( double ) ) || chkseg[t].skipForeward != segment[count+t].skipForeward || chkseg[t].skipBackward != segment[count+t].skipBackward )
				{
					reporterr( "\nARCHECK MISMATCH in segment %d of %d\n", t, tmpint ); exit( 1 );
				}
			pt_acc[27] += 0.0;
		}
#endif
#ifdef PROGTIMING
		pt_acc[27] += 1.0;
#endif

//		if( lag == -50 ) exit( 1 );
		
		if( count+tmpint > MAXSEG -3 ) ErrorExit( "TOO MANY SEGMENTS.\n" );


		if( tmpint == 0 ) break; // 060430 iinoka ?
		while( tmpint-- > 0 )
		{
#if 0
			if( segment[count].end - segment[count].start < fftWinSize )
			{
				count++;
				continue;
			}
#endif
			if( lag > 0 )
			{
				segment1[count].start  = segment[count].start ;
				segment1[count].end    = segment[count].end   ;
				segment1[count].center = segment[count].center;
				segment1[count].score  = segment[count].score;

				segment2[count].start  = segment[count].start  + lag;
				segment2[count].end    = segment[count].end    + lag;
				segment2[count].center = segment[count].center + lag;
				segment2[count].score  = segment[count].score       ;
			}
			else
			{
				segment1[count].start  = segment[count].start  - lag;
				segment1[count].end    = segment[count].end    - lag;
				segment1[count].center = segment[count].center - lag;
				segment1[count].score  = segment[count].score       ;

				segment2[count].start  = segment[count].start ;
				segment2[count].end    = segment[count].end   ;
				segment2[count].center = segment[count].center;
				segment2[count].score  = segment[count].score ;
			}
#if 0
			fprintf( stderr, "in 1 %d\n", segment1[count].center );
			fprintf( stderr, "in 2 %d\n", segment2[count].center );
#endif
			segment1[count].pair = &segment2[count];
			segment2[count].pair = &segment1[count];
			count++;
		}
	}
#if 0
	if( !kobetsubunkatsu && fftkeika )
		fprintf( stderr, "%d anchors found\r", count );
#endif
	if( !count && fftNoAnchStop )
		ErrorExit( "Cannot detect anchor!" );
#if 0
	fprintf( stderr, "RESULT before sort:\n" );
	for( l=0; l<count+1; l++ )
	{
		fprintf( stderr, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( stderr, "%d score = %f\n", segment2[l].center, segment1[l].score );
	}
#endif

#if KEIKA
	fprintf( stderr, "done. (%d anchors)\n", count );
	fprintf( stderr, "Aligning anchors ... " );
#endif
	for( i=0; i<count; i++ )
	{
		sortedseg1[i] = &segment1[i];
		sortedseg2[i] = &segment2[i];
	}
#if 0
	tmpsort( count, sortedseg1 ); 
	tmpsort( count, sortedseg2 ); 
	qsort( sortedseg1, count, sizeof( Segment * ), segcmp );
	qsort( sortedseg2, count, sizeof( Segment * ), segcmp );
#else
	mymergesort( 0, count-1, sortedseg1 ); 
	mymergesort( 0, count-1, sortedseg2 ); 
#endif
	for( i=0; i<count; i++ ) sortedseg1[i]->number = i;
	for( i=0; i<count; i++ ) sortedseg2[i]->number = i;


	if( kobetsubunkatsu )
	{
		for( i=0; i<count; i++ )
	    {
			cut1[i+1] = sortedseg1[i]->center;
			cut2[i+1] = sortedseg2[i]->center;
		}
		cut1[0] = 0;
		cut2[0] = 0;
		cut1[count+1] = len1;
		cut2[count+1] = len2;
		count += 2;
	}
	else
	{
		if( crossscoresize < count+2 )
		{
			crossscoresize = count+2;
#if 1
			if( fftkeika ) fprintf( stderr, "######allocating crossscore, size = %d\n", crossscoresize );
#endif
			if( crossscore ) FreeDoubleMtx( crossscore );
			crossscore = AllocateDoubleMtx( crossscoresize, crossscoresize );
		}
		for( i=0; i<count+2; i++ ) for( j=0; j<count+2; j++ )
			crossscore[i][j] = 0.0;
		for( i=0; i<count; i++ )
		{
			crossscore[segment1[i].number+1][segment1[i].pair->number+1] = segment1[i].score;
			cut1[i+1] = sortedseg1[i]->center;
			cut2[i+1] = sortedseg2[i]->center;
		}

#if 0
		fprintf( stderr, "AFTER SORT\n" );
		for( i=0; i<count+1; i++ ) fprintf( stderr, "%d, %d\n", cut1[i], cut2[i] );
		fprintf( stderr, "crossscore = \n" );
		for( i=0; i<count+1; i++ )
		{
			for( j=0; j<count+1; j++ )
				fprintf( stderr, "%.0f ", crossscore[i][j] );
			fprintf( stderr, "\n" );
		}
#endif

		crossscore[0][0] = 10000000.0;
		cut1[0] = 0; 
		cut2[0] = 0;
		crossscore[count+1][count+1] = 10000000.0;
		cut1[count+1] = len1;
		cut2[count+1] = len2;
		count += 2;
		count0 = count;
	
		blockAlign2( cut1, cut2, sortedseg1, sortedseg2, crossscore, &count );

//		if( count-count0 )
//			fprintf( stderr, "%d unused anchors\n", count0-count );

		if( !kobetsubunkatsu && fftkeika )
			fprintf( stderr, "%d anchors found\n", count );
		if( fftkeika )
		{
			if( count0 > count )
			{
#if 0
				fprintf( stderr, "\7 REPEAT!? \n" ); 
#else
				fprintf( stderr, "REPEAT!? \n" ); 
#endif
				if( fftRepeatStop ) exit( 1 );
			}
#if KEIKA
			else fprintf( stderr, "done\n" );
#endif
		}
	}

#if 0
	fftfp = fopen( "fft", "a" );
	fprintf( fftfp, "RESULT after sort:\n" );
	for( l=0; l<count; l++ )
	{
		fprintf( fftfp, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( fftfp, "%d\n", segment2[l].center );
	}
	fclose( fftfp );
#endif

#if 0
	fprintf( stderr, "RESULT after blckalign:\n" );
	for( l=0; l<count+1; l++ )
	{
		fprintf( stderr, "cut : %d %d\n", cut1[l], cut2[l] );
	}
#endif

#if 0
	fprintf( trap_g, "Devided to %d segments\n", count-1 );
	fprintf( trap_g, "%d  %d forg\n", MIN( clus1, clus2 ), count-1 );
#endif

#ifdef PROGTIMING
	pt_acc[28] += 1.0;
#endif
	PT_ADD(2);
	inpl = ( count - 1 == 1 && cut1[0] == 0 && cut2[0] == 0 && cut1[1] == len1 && cut2[1] == len2
	    && alg == 'A' && !scoringmatrices && !kobetsubunkatsu && !codonpos && !codonscore && !constraint && chudanpt == NULL
	    && !( clus1 == 1 && clus2 == 1 ) && fa_inplace_on() ); // work/v8s: see below
#ifdef PROGTIMING
	if( !inpl || v8check() ) // (the in-place check aligns the loop's copies in tmpres)
#else
	if( !inpl )
#endif
	{
		int rowlen = len1 + len2 + 100;
		if( rowlen > alloclen ) rowlen = alloclen;
#ifdef PROGTIMING
		double _fr0 = pt_now();
#endif
		Falign_rows( clus1, clus2, rowlen, &result1, &result2, &tmpres1, &tmpres2 ); // work/prog2: per-thread arena
#ifdef PROGTIMING
		pt_acc[42] += pt_now() - _fr0; // work/v8s
#endif
		for( j=0; j<clus1; j++ ) result1[j][0] = 0;
		for( j=0; j<clus2; j++ ) result2[j][0] = 0;
	}
	totallen = 0;
	off1 = off2 = 0;
	totalscore = 0.0;
	*fftlog = -1;
//	reporterr( "\nin Falign(), *fftlog = %d\n", *fftlog );
	/* work/v8s: one segment covering both groups (no anchor kept: most Falign calls of the large
	   groups on the critical chain of real rRNA sets) is aligned in place, in seq1/seq2.  The loop
	   below would copy every row into tmpres (strncpy of the whole row), align tmpres in place,
	   append it to the empty result rows and copy those back into seq1/seq2: the same bytes end up
	   in seq1/seq2 either way, and A__align gets byte-identical rows and the same arguments (headgp
	   = tailgp = outgap, boundary gap patterns all 'o' because the cut points are 0 and len).  The
	   rows hold what the final strcpy would have written (<= alloclen + 1 bytes, the LENGTH OVER
	   check is kept).  Not used with chudan (an aborted A__align must leave seq1/seq2 untouched).
	   MAFFT_V8CHECK (timing builds) also aligns the copies the loop would make and compares;
	   MAFFT_FAINPLACE=0 (timing builds) turns it off (A/B). */
	if( inpl )
	{
#ifdef PROGTIMING
		int chk = v8check();
		if( chk ) // the loop's copies, aligned below for comparison
		{
			for( j=0; j<clus1; j++ ) { strncpy( tmpres1[j], seq1[j], len1 ); tmpres1[j][len1] = 0; }
			for( j=0; j<clus2; j++ ) { strncpy( tmpres2[j], seq2[j], len2 ); tmpres2[j][len2] = 0; }
		}
#endif
		*fftlog += 1;
		for( j=0; j<clus1; j++ ) sgap1[j] = egap1[j] = 'o';
		for( j=0; j<clus2; j++ ) sgap2[j] = egap2[j] = 'o';
		PT_ADD(3);
		totalscore += A__align( n_dynamicmtx, penalty, penalty_ex, seq1, seq2, eff1, eff2, clus1, clus2, alloclen, 0, &dumdb, sgap1, sgap2, egap1, egap2, chudanpt, chudanref, chudanres, outgap, outgap, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
		PT_ADD(4);
		nlen = strlen( seq1[0] );
		if( nlen > alloclen )
		{
			fprintf( stderr, "totallen=%d +  nlen=%d > alloclen = %d\n", 0, nlen, alloclen );
			ErrorExit( "LENGTH OVER in Falign\n " );
		}
#ifdef PROGTIMING
		if( chk )
		{
			double ts = 0.0, dd = 0.0;
			ts += A__align( n_dynamicmtx, penalty, penalty_ex, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, 0, &dd, sgap1, sgap2, egap1, egap2, chudanpt, chudanref, chudanres, outgap, outgap, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
			for( j=0; j<clus1; j++ ) if( strcmp( tmpres1[j], seq1[j] ) ) { reporterr( "\nV8CHECK MISMATCH in the in-place Falign segment, group 1 row %d\n", j ); exit( 1 ); }
			for( j=0; j<clus2; j++ ) if( strcmp( tmpres2[j], seq2[j] ) ) { reporterr( "\nV8CHECK MISMATCH in the in-place Falign segment, group 2 row %d\n", j ); exit( 1 ); }
			if( memcmp( &ts, &totalscore, sizeof( double ) ) ) { reporterr( "\nV8CHECK MISMATCH in the in-place Falign score\n" ); exit( 1 ); }
			v8check_count( 3 );
		}
#endif
		goto inplace_done;
	}
#ifdef enablemultithread
	/* work/prog3: segments in parallel (Falign_psegs) when helpers are idle and there is enough work */
	if( count - 1 >= 2 && pseg_on() && alg == 'A' && !scoringmatrices && !kobetsubunkatsu && !codonpos && !codonscore && chudanpt == NULL
	    && ( pseg_force() || ( dphelp_idle() > 0 && (double)( clus1 + clus2 ) * ( len1 + len2 ) >= 2.0e5 ) ) )
	{
		totalscore = Falign_psegs( seq1, seq2, clus1, clus2, count-1, cut1, cut2, n_dynamicmtx, eff1, eff2, alloclen, result1, result2 );
		*fftlog = count - 2; // as after count-1 increments from -1
#ifdef PROGTIMING
		if( getenv( "MAFFT_PSEGCHECK" ) ) // the sequential loop again, into copies, compared bitwise
		{
			char **save1 = malloc( sizeof( char * ) * clus1 ), **save2 = malloc( sizeof( char * ) * clus2 );
			double ts = totalscore;
			for( j=0; j<clus1; j++ ) { save1[j] = strdup( result1[j] ); result1[j][0] = 0; }
			for( j=0; j<clus2; j++ ) { save2[j] = strdup( result2[j] ); result2[j][0] = 0; }
			pseg_check_pending = 1;
			pseg_check_save1 = save1; pseg_check_save2 = save2; pseg_check_ts = ts;
			totalscore = 0.0;
			*fftlog = -1;
			goto pseg_sequential;
		}
#endif
		goto pseg_done;
	}
#endif
#ifdef PROGTIMING
	pseg_sequential:
#endif
	for( i=0; i<count-1; i++ )
	{
		*fftlog += 1;
		if( i == 0 ) headgp = outgap; else headgp = 1;
		if( i == count-2 ) tailgp = outgap; else tailgp = 1;


#if 0
		if( cut1[i] ) // chuui
		{
//			getkyokaigap( sgap1, tmpres1, nlen-1, clus1 );
//			getkyokaigap( sgap2, tmpres2, nlen-1, clus2 );
			getkyokaigap( sgap1, tmpres1, nlen-1, clus1 );
			getkyokaigap( sgap2, tmpres2, nlen-1, clus2 );
		}
		else
		{
			for( j=0; j<clus1; j++ ) sgap1[j] = 'o';
			for( j=0; j<clus2; j++ ) sgap2[j] = 'o';
		}
		if( cut1[i+1] != len1 ) // chuui
		{       
			getkyokaigap( egap1, seq1, cut1[i+1], clus1 );
			getkyokaigap( egap2, seq2, cut2[i+1], clus2 );
		}       
		else    
		{       
			for( j=0; j<clus1; j++ ) egap1[j] = 'o';
			for( j=0; j<clus2; j++ ) egap2[j] = 'o';
		}
#else
		if( cut1[i] )
			getkyokaigap( sgap1, seq1, cut1[i]-1, clus1 );
		else
			for( j=0; j<clus1; j++ ) sgap1[j] = 'o';
		if( cut2[i] )
			getkyokaigap( sgap2, seq2, cut2[i]-1, clus2 );
		else
			for( j=0; j<clus2; j++ ) sgap2[j] = 'o';

		if( cut1[i+1] != len1 )
			getkyokaigap( egap1, seq1, cut1[i+1], clus1 );
		else    
			for( j=0; j<clus1; j++ ) egap1[j] = 'o';
		if( cut2[i+1] != len2 )
			getkyokaigap( egap2, seq2, cut2[i+1], clus2 );
		else    
			for( j=0; j<clus2; j++ ) egap2[j] = 'o';
#endif
#if 0
		{
			fprintf( stderr, "kyokkaigap1(%d)=", cut1[i]-1 );
			for( j=0; j<clus1; j++ )
				fprintf( stderr, "%c", sgap1[j] );
			fprintf( stderr, "=kyokkaigap1-start\n" );
		}
		{
			fprintf( stderr, "kyokkaigap2(%d)=", cut2[i]-1 );
			for( j=0; j<clus2; j++ )
				fprintf( stderr, "%c", sgap2[j] );
			fprintf( stderr, "=kyokkaigap2-start\n" );
		}
		{
			fprintf( stderr, "kyokkaigap1(%d)=", cut1[i]-1 );
			for( j=0; j<clus1; j++ )
				fprintf( stderr, "%c", egap1[j] );
			fprintf( stderr, "=kyokkaigap1-end\n" );
		}
		{
			fprintf( stderr, "kyokkaigap2(%d)=", cut2[i]-1 );
			for( j=0; j<clus2; j++ )
				fprintf( stderr, "%c", egap2[j] );
			fprintf( stderr, "=kyokkaigap2-end\n" );
		}
#endif

#if DEBUG
		fprintf( stderr, "DP %03d / %03d %4d to ", i+1, count-1, totallen );
#else
#if KEIKA
		fprintf( stderr, "DP %03d / %03d\r", i+1, count-1 );
#endif
#endif

		//reporterr( "cut1[] = %d\n", cut1[i] );
		//reporterr( "cut2[] = %d\n", cut2[i] );

		if( !rows_copy_par( tmpres1, seq1, clus1, cut1[i], cut1[i+1]-cut1[i], 0 ) ) // work/v8s
		for( j=0; j<clus1; j++ )
		{
			strncpy( tmpres1[j], seq1[j]+cut1[i], cut1[i+1]-cut1[i] );
			tmpres1[j][cut1[i+1]-cut1[i]] = 0;
		}
		if( kobetsubunkatsu && fftkeika ) commongappick( clus1, tmpres1 ); //dvtditr に呼ばれたとき fftkeika=1
//		if( kobetsubunkatsu ) commongappick( clus1, tmpres1 );
		if( !rows_copy_par( tmpres2, seq2, clus2, cut2[i], cut2[i+1]-cut2[i], 0 ) ) // work/v8s
		for( j=0; j<clus2; j++ )
		{
			strncpy( tmpres2[j], seq2[j]+cut2[i], cut2[i+1]-cut2[i] );
			tmpres2[j][cut2[i+1]-cut2[i]] = 0;
		}
		if( kobetsubunkatsu && fftkeika ) commongappick( clus2, tmpres2 ); //dvtditr に呼ばれたとき fftkeika=1
//		if( kobetsubunkatsu ) commongappick( clus2, tmpres2 );

		PT_ADD(3);
		if( constraint )
		{
			fprintf( stderr, "Not supported\n" );
			exit( 1 );
		}
#if 0
		fprintf( stderr, "i=%d, before alignment", i );
		fprintf( stderr, "%4d\n", totallen );
		fprintf( stderr, "\n\n" );
		for( j=0; j<clus1; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres1[j] );
		}
		fprintf( stderr, "-------\n" );
		for( j=0; j<clus2; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres2[j] );
		}
#endif

#if 0
		fprintf( stdout, "writing input\n" );
		for( j=0; j<clus1; j++ )
		{
			fprintf( stdout, ">%d of GROUP1\n", j );
			fprintf( stdout, "%s\n", tmpres1[j] );
		}
		for( j=0; j<clus2; j++ )
		{
			fprintf( stdout, ">%d of GROUP2\n", j );
			fprintf( stdout, "%s\n", tmpres2[j] );
		}
		fflush( stdout );
#endif
		switch( alg )
		{
			case( 'a' ):
				totalscore += Aalign( tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen );
				break;
			case( 'M' ):
					if( scoringmatrices ) // called by tditeration.c
						totalscore += MSalignmm_variousdist( NULL, scoringmatrices, NULL, tmpres1, tmpres2, eff1, eff2, eff1s, eff2s, clus1, clus2, alloclen, sgap1, sgap2, egap1, egap2, chudanpt, chudanref, chudanres, headgp, tailgp );
					else
						totalscore += MSalignmm( n_dynamicmtx, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, sgap1, sgap2, egap1, egap2, chudanpt, chudanref, chudanres, headgp, tailgp, NULL, NULL, NULL, 0.0, 0.0 );
//						totalscore += MSalignmm( n_dis_consweight_multi, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, sgap1, sgap2, egap1, egap2, chudanpt, chudanref, chudanres, headgp, tailgp );
				break;
			case( 'd' ):
				if( clus1 == 1 && clus2 == 1 )
				{
					totalscore += G__align11( n_dynamicmtx, tmpres1, tmpres2, alloclen, headgp, tailgp );
				}
				else
				{
					if( scoringmatrices ) // called by tditeration.c
					{
						totalscore += D__align_variousdist( whichmtx, scoringmatrices, NULL, tmpres1, tmpres2, eff1, eff2, eff1s, eff2s, clus1, clus2, alloclen, 0, &dumdb, sgap1, sgap2, egap1, egap2, chudanpt, chudanref, chudanres, headgp, tailgp );
					}
					else
					totalscore += D__align( n_dynamicmtx, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, 0, &dumdb, sgap1, sgap2, egap1, egap2, chudanpt, chudanref, chudanres, headgp, tailgp );
				}
				break;
			case( 'A' ):
				if( clus1 == 1 && clus2 == 1 )
				{
					if( codonpos || codonscore )
					{
//						reporterr( "calling G__align11psg\n" );
						totalscore += G__align11psg( codonscoremtx, n_dynamicmtx, tmpres1, tmpres2, alloclen, headgp, tailgp, gstart+cut1[i], gend+cut1[i] );
					}
					else
						totalscore += G__align11( n_dynamicmtx, tmpres1, tmpres2, alloclen, headgp, tailgp );
				}
				else
				{
					if( codonpos )
					{
						reporterr( "\n\ncodonpos will be soon supported for a reference MSA. For now, use a single sequence as reference.\n\n\n" );
						exit( 1 );
					}
					if( scoringmatrices ) // called by tditeration.c
					{
						totalscore += A__align_variousdist( whichmtx, scoringmatrices, NULL, penalty, penalty_ex, tmpres1, tmpres2, eff1, eff2, eff1s, eff2s, clus1, clus2, alloclen, 0, &dumdb, sgap1, sgap2, egap1, egap2, chudanpt, chudanref, chudanres, headgp, tailgp );
					}
					else
						totalscore += A__align( n_dynamicmtx, penalty, penalty_ex, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, 0, &dumdb, sgap1, sgap2, egap1, egap2, chudanpt, chudanref, chudanres, headgp, tailgp, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
				}
				break;
			default:
				fprintf( stderr, "alg = %c\n", alg );
				ErrorExit( "ERROR IN SOURCE FILE Falign.c" );
				break;
		}

#ifdef enablemultithread
		if( chudanres && *chudanres )
		{
//			fprintf( stderr, "\n\n## CHUUDAN!!! at Falign_localhom\n" );
//			Added 2021/Jul/25.
			// work/prog2: result1/2 and tmpres1/2 live in a per-thread arena; tmpseq1/2 alias seq1/2
#if RND
			FreeCharMtx( tmpseq1 );
			FreeCharMtx( tmpseq2 );
#endif
			free( sgap1 );
			free( egap1 );
			free( sgap2 );
			free( egap2 );
			free( tmpptr1 );
			free( tmpptr2 );
#if RND
			FreeCharMtx( rndseq1 );
			FreeCharMtx( rndseq2 );
#endif
			return( -1.0 );
		}
#endif

		PT_ADD(4);
		nlen = strlen( tmpres1[0] );
		if( totallen + nlen > alloclen )
		{
			fprintf( stderr, "totallen=%d +  nlen=%d > alloclen = %d\n", totallen, nlen, alloclen );
			ErrorExit( "LENGTH OVER in Falign\n " );
		}
		/* work/prog: append at the known end of each row instead of strcat(), which rescanned
		   the whole partial result of every row for every segment (quadratic in the number of
		   segments).  Within a group all rows are aligned, hence equally long, so the bytes are
		   the same as with strcat(). */
		{
			int add1 = strlen( tmpres1[0] ), add2 = strlen( tmpres2[0] );
			if( !rows_copy_par( result1, tmpres1, clus1, off1, add1, 1 ) ) // work/v8s
			for( j=0; j<clus1; j++ ) memcpy( result1[j]+off1, tmpres1[j], add1+1 );
			if( !rows_copy_par( result2, tmpres2, clus2, off2, add2, 1 ) )
			for( j=0; j<clus2; j++ ) memcpy( result2[j]+off2, tmpres2[j], add2+1 );
			off1 += add1;
			off2 += add2;
		}
		PT_ADD(5);
		totallen += nlen;
#if 0
		fprintf( stderr, "$#####$$$$ i=%d", i );
		fprintf( stderr, "%4d\n", totallen );
		fprintf( stderr, "\n\n" );
		for( j=0; j<clus1; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres1[j] );
		}
		fprintf( stderr, "-------\n" );
		for( j=0; j<clus2; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres2[j] );
		}
#endif
	}
//	reporterr( "\nafter Falign(), *fftlog = %d\n", *fftlog );

#if KEIKA
	fprintf( stderr, "DP ... done   \n" );
#endif
#ifdef PROGTIMING
	if( pseg_check_pending )
	{
		pseg_check_pending = 0;
		for( j=0; j<clus1; j++ ) if( strcmp( pseg_check_save1[j], result1[j] ) ) { reporterr( "\nPSEGCHECK MISMATCH, group 1 row %d\n", j ); exit( 1 ); }
		for( j=0; j<clus2; j++ ) if( strcmp( pseg_check_save2[j], result2[j] ) ) { reporterr( "\nPSEGCHECK MISMATCH, group 2 row %d\n", j ); exit( 1 ); }
		if( memcmp( &pseg_check_ts, &totalscore, sizeof( double ) ) ) { reporterr( "\nPSEGCHECK MISMATCH in totalscore\n" ); exit( 1 ); }
		for( j=0; j<clus1; j++ ) free( pseg_check_save1[j] );
		for( j=0; j<clus2; j++ ) free( pseg_check_save2[j] );
		free( pseg_check_save1 ); free( pseg_check_save2 );
		__atomic_add_fetch( &pseg_checked, 1, __ATOMIC_RELAXED );
	}
#endif
#ifdef enablemultithread
	pseg_done:
#endif

	if( !rows_copy_par( seq1, result1, clus1, 0, off1, 2 ) ) // work/v8s
	for( j=0; j<clus1; j++ ) strcpy( seq1[j], result1[j] );
	if( !rows_copy_par( seq2, result2, clus2, 0, off2, 2 ) )
	for( j=0; j<clus2; j++ ) strcpy( seq2[j], result2[j] );
	inplace_done: // work/v8s
#if 0
	for( j=0; j<clus1; j++ ) 
	{
		fprintf( stderr, "in Falign, %s\n", result1[j] );
	}
	fprintf( stderr, "- - - - - - - - - - -\n" );
	for( j=0; j<clus2; j++ ) 
	{
		fprintf( stderr, "in Falign, %s\n", result2[j] );
	}
#endif

	// work/prog2: result1/2 and tmpres1/2 live in a per-thread arena; tmpseq1/2 alias seq1/2
#if RND
	FreeCharMtx( tmpseq1 );
	FreeCharMtx( tmpseq2 );
#endif
	free( sgap1 );
	free( egap1 );
	free( sgap2 );
	free( egap2 );
	free( tmpptr1 );
	free( tmpptr2 );
#if RND
	FreeCharMtx( rndseq1 );
	FreeCharMtx( rndseq2 );
#endif

	PT_ADD(6);
	return( totalscore );
}



// 2019/Jun
static double estimategapfreq( int n, char **s )
{
	int i, j, len;
	char *seq, *pt1, *pt2;
	double f, fv;

	len = strlen( s[0] );
	seq = calloc( len+1, sizeof( char ) );

	fv = 0.0;
	for( i=0; i<n; i++  )
	{
		pt1 = s[i];
		while( *pt1 == '-' ) pt1++;
		pt2 = seq;
		while( *pt1 != 0 ) *pt2++ = *pt1++;
		*pt2 = *pt1; // 0
		pt1 = pt2-1;
		while( *pt1 == '-' ) pt1--;
		*(pt1+1) = 0;	
//		reporterr( "seq[i]=%s\n", s[i] );
//		reporterr( "seq=%s\n", seq );
		len = pt1-seq+1;
		f = 0.0;
		for( j=0; j<len; j++ )
			if( seq[j] == '-' ) f+=1.0;
		fv += f/(double)len;
	}
	free( seq );
	return( fv/(double)n );
}

static int getterminalmargin( int lshorter, double groupsizefac )
{
	return ( lshorter * 2.0 + terminalmargin ) * groupsizefac;
//	return ( lshorter * 1.1 + 10 ) * groupsizefac;
}

#if 0
static int estimatenogaplen( int n, char **s, int start, int end )
{
	int i, j, l, d, o;
	int minl;
	char *c;

	if( start < end ) d = 1; else d = -1;

	reporterr( "\nin estimatenogaplen, d=%d\n", d );

	minl = (end-start)*d;
	c = calloc( minl, sizeof( char ) );

	for( j=start; j!=end; j+=d )
	{
		o = 0;
		for( i=0; i<n; i++ )
		{
			if( s[i][j] != '-' ) o++;
			if( o >= n / 2 ) break;
		}
		if( o >= n/2 ) c[(j-start)*d] = 'o';
		else c[(j-start)*d] = '-';
	}
	c[(j-start)*d] = 0;

	reporterr( "c=%s\n", c );

	l = 0;
	for( j=start; j!=end; j+=d )
		if( c[j] == 'o' ) l++;
	reporterr( "l=%d\n", l );
	free( c );
	return( l );
}

static int nogapmargin( int n, char **s, int start, int end, int m )
{
	int i, j, l, d;
	int minl;

	if( start < end ) d = 1; else d = -1;

//	reporterr( "\nin nogapmargin, d=%d\n", d );

	minl = (end-start)*d;
	for( i=0; i<n; i++ )
	{
		l = 0;
		for( j=start; j!=end; j+=d )
		{
			if( s[i][j] != '-' ) l++;
			if( l>m ) break;
		}
//		reporterr( "i=%d, l=%d, j=%d\n", i, l, j );
		if( (j-start)*d < minl ) minl = (j-start)*d;
	}
	minl += 1;
//	reporterr( "minl=%d, so returning %d\n", minl, start+minl*d );
	return( start + minl*d );
}
#endif

double Falign_givenanchors( ExtAnch *pairanch, 
			  int **whichmtx, double ***scoringmatrices, 
			  double **n_dynamicmtx,
			  char  **seq1, char  **seq2, 
			  double *eff1, double *eff2, 
			  double **eff1s, double **eff2s,
			  int    clus1, int    clus2,
			  int alloclen, int *fftlog )
{

	int i, j;
	int nlen, nlen2, nlen4;
	static TLS int prevalloclen = 0;
	//static TLS int crossscoresize = 0;
	//static TLS char **tmpseq1 = NULL;
	//static TLS char **tmpseq2 = NULL;
	//static TLS char **tmpptr1 = NULL;
	//static TLS char **tmpptr2 = NULL;
	static TLS char **tmpres1 = NULL;
	static TLS char **tmpres2 = NULL;
	static TLS char **result1 = NULL;
	static TLS char **result2 = NULL;
#if RND
	//static TLS char **rndseq1 = NULL;
	//static TLS char **rndseq2 = NULL;
#endif
	//static TLS Fukusosuu **seqVector1 = NULL;
	//static TLS Fukusosuu **seqVector2 = NULL;
	//static TLS Fukusosuu **naiseki = NULL;   
	//static TLS Fukusosuu *naisekiNoWa = NULL; 
	//static TLS double *soukan = NULL;
	//static TLS double **crossscore = NULL;
	int nlentmp;
	//static TLS int *kouho = NULL;
	//static TLS Segment *segment = NULL;
	//static TLS Segment *segment1 = NULL;
	//static TLS Segment *segment2 = NULL;
	//static TLS Segment **sortedseg1 = NULL;
	//static TLS Segment **sortedseg2 = NULL;
	static TLS int *alignorcopy = NULL;
	static TLS int *cut1 = NULL;
	static TLS int *cut2 = NULL;
	static TLS char *sgap1, *egap1, *sgap2, *egap2;
	static TLS int localalloclen = 0;
//	int lag;
//	int tmpint;
	int count, count0;
	int len1, len2;
	int totallen;
	double totalscore;
//	int nkouho = 0;
	int headgp, tailgp;
//	double dumfl = 0.0;
	int orilen1, orilen2;
	int cutadd;
	int starttermcut1, starttermcut2, endtermcut1, endtermcut2;
	double marginfac1, marginfac2;

	if( seq1 == NULL )
	{
		if( result1 ) 
		{
//			fprintf( stderr, "### Freeing localarrays in Falign\n" );
			localalloclen = 0;
			prevalloclen = 0;
			//crossscoresize = 0;
			mymergesort( 0, 0, NULL );
			//alignableReagion( 0, 0, NULL, NULL, NULL, NULL, NULL );
			//fft( 0, NULL, 1 );
			A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
			A__align_variousdist( NULL, NULL, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
			D__align_variousdist( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
			G__align11( NULL, NULL, NULL, 0, 0, 0 );
			//blockAlign2( NULL, NULL, NULL, NULL, NULL, NULL );
			//if( crossscore ) FreeDoubleMtx( crossscore );
			//crossscore = NULL; // reallocate sareru kanousei ga arunode.
			FreeCharMtx( result1 ); result1 = NULL;
			FreeCharMtx( result2 );
			FreeCharMtx( tmpres1 );
			FreeCharMtx( tmpres2 );
			//FreeCharMtx( tmpseq1 );
			//FreeCharMtx( tmpseq2 );
			free( sgap1 );
			free( egap1 );
			free( sgap2 );
			free( egap2 );
			//free( kouho );
			free( alignorcopy );
			free( cut1 );
			free( cut2 );
			//free( tmpptr1 );
			//free( tmpptr2 );
			//free( segment );
			//free( segment1 );
			//free( segment2 );
			//free( sortedseg1 );
			//free( sortedseg2 );
			//if( !kobetsubunkatsu )
			//{
			//	FreeFukusosuuMtx ( seqVector1 );
			//	FreeFukusosuuMtx ( seqVector2 );
			//	FreeFukusosuuVec( naisekiNoWa );
			//	FreeFukusosuuMtx( naiseki );
			//	FreeDoubleVec( soukan );
			//}
		}
		else
		{
//			fprintf( stderr, "Did not allocate localarrays in Falign\n" );
		}

		return( 0.0 );
	}



	len1 = strlen( seq1[0] );
	len2 = strlen( seq2[0] );
	nlentmp = MAX( len1, len2 );

	nlen = 1;
	while( nlentmp >= nlen ) nlen <<= 1;
#if 0
	fprintf( stderr, "###   nlen    = %d\n", nlen );
#endif

	nlen2 = nlen/2; nlen4 = nlen2 / 2;

#if 0
	fprintf( stderr, "len1 = %d, len2 = %d\n", len1, len2 );
	fprintf( stderr, "nlentmp = %d, nlen = %d\n", nlentmp, nlen );
#endif

	if( prevalloclen != alloclen ) // Falign_noudp mo kaeru
	{
		if( prevalloclen )
		{
			FreeCharMtx( result1 );
			FreeCharMtx( result2 );
			FreeCharMtx( tmpres1 );
			FreeCharMtx( tmpres2 );
		}
//		fprintf( stderr, "\n\n\nreallocating ...\n" ); 
		result1 = AllocateCharMtx( njob, alloclen ); // ato de loca nseq ni kakihaosu
		result2 = AllocateCharMtx( njob, alloclen );
		tmpres1 = AllocateCharMtx( njob, alloclen );
		tmpres2 = AllocateCharMtx( njob, alloclen );
		prevalloclen = alloclen;
	}

	if( !localalloclen )
	{
		sgap1 = AllocateCharVec( njob );
		egap1 = AllocateCharVec( njob );
		sgap2 = AllocateCharVec( njob );
		egap2 = AllocateCharVec( njob );
		//kouho = AllocateIntVec( NKOUHO_LONG );
		alignorcopy = AllocateIntVec( MAXSEG_GIVENANCHORS );
		cut1 = AllocateIntVec( MAXSEG_GIVENANCHORS );
		cut2 = AllocateIntVec( MAXSEG_GIVENANCHORS );
		//tmpptr1 = AllocateCharMtx( njob, 0 );
		//tmpptr2 = AllocateCharMtx( njob, 0 );
		//segment = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		//segment1 = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		//segment2 = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		//sortedseg1 = (Segment **)calloc( MAXSEG, sizeof( Segment * ) );
		//sortedseg2 = (Segment **)calloc( MAXSEG, sizeof( Segment * ) );
		//if( !( segment && segment1 && segment2 && sortedseg1 && sortedseg2 ) )
		//	ErrorExit( "Allocation error\n" );

		//if     ( scoremtx == -1 ) n20or4or2 = 1;
		//else if( fftscore )       n20or4or2 = 1;
		//else                      n20or4or2 = 20;
	}
	if( localalloclen < nlen )
	{
		if( localalloclen )
		{
#if 1
			//if( !kobetsubunkatsu )
			//{
			//	FreeFukusosuuMtx ( seqVector1 );
			//	FreeFukusosuuMtx ( seqVector2 );
			//	FreeFukusosuuVec( naisekiNoWa );
			//	FreeFukusosuuMtx( naiseki );
			//	FreeDoubleVec( soukan );
			//}
			//FreeCharMtx( tmpseq1 );
			//FreeCharMtx( tmpseq2 );
#endif
#if RND
			//FreeCharMtx( rndseq1 );
			//FreeCharMtx( rndseq2 );
#endif
		}


		//tmpseq1 = AllocateCharMtx( njob, nlen );
		//tmpseq2 = AllocateCharMtx( njob, nlen );
		//if( !kobetsubunkatsu )
		//{
		//	naisekiNoWa = AllocateFukusosuuVec( nlen );
		//	naiseki = AllocateFukusosuuMtx( n20or4or2, nlen );
		//	seqVector1 = AllocateFukusosuuMtx( n20or4or2, nlen+1 );
		//	seqVector2 = AllocateFukusosuuMtx( n20or4or2, nlen+1 );
		//	soukan = AllocateDoubleVec( nlen+1 );
		//}
#if RND
		//rndseq1 = AllocateCharMtx( njob, nlen );
		//rndseq2 = AllocateCharMtx( njob, nlen );
		//for( i=0; i<njob; i++ )
		//{
		//	generateRndSeq( rndseq1[i], nlen );
		//	generateRndSeq( rndseq2[i], nlen );
		//}
#endif
		localalloclen = nlen;
	}
	

	marginfac1 = 1.0 + estimategapfreq( clus1, seq1 );
	marginfac2 = 1.0 + estimategapfreq( clus2, seq2 );
	starttermcut1 = starttermcut2 = 0;
	endtermcut1 = endtermcut2 = 0;
//	reporterr( "marginfac1=%f\n", marginfac1 );
//	reporterr( "marginfac2=%f\n", marginfac2 );

//	reporterr( "length1,length2=%d,%d\n", len1, len2  );
//	reporterr( "pairanch when uwagaki: %d:%d\n", pairanch[0].starti, pairanch[0].startj );
//	reporterr( "pairanch when uwagaki: i=%d, j=%d\n", pairanch[0].i, pairanch[0].j );
	count = count0 = 0;
	cut1[0] = 0;
	cut2[0] = 0;
	alignorcopy[0] = 'a';
//	while( pairanch[count].i == 0 &&  pairanch[count].j == 0 ) // ato de kentou
	while( pairanch[count0].i > -1 )
	{
		if( pairanch[count0].starti == -1 ) 
		{
			count0++;
			continue;
		}
		if( count+2 > MAXSEG_GIVENANCHORS -3 ) ErrorExit( "TOO MANY SEGMENTS.\n" );
#if 1 // mattan no tansaku hann'i wo seigen
		if( count == 0 ) 	
		{
//			if( pairanch[count0].starti - pairanch[count0].startj > TERMINALSEGMENTLENGTH ) // you kentou
//			nogaplen1 = estimatenogaplen( clus1, seq1, pairanch[count0].starti, 0 );
//			nogaplen2 = estimatenogaplen( clus2, seq2, pairanch[count0].startj, 0 );
			if( pairanch[count0].starti > getterminalmargin(pairanch[count0].startj,marginfac1) )
			{
//				alignorcopy[1] = 'A';
				reporterr( "check 1, because starti=%d > startj=%d -> %d (clus1=%d)\n", pairanch[count0].starti, pairanch[count0].startj, getterminalmargin(pairanch[count0].startj,marginfac1), clus1 );
				cutadd = pairanch[count0].starti - getterminalmargin(pairanch[count0].startj,marginfac1);
				reporterr( "cutadd(1)=%d\n", cutadd );
//				if( 0 || cutadd > TERMINALMARGIN(0) ) // iranai
				{
					cut1[1] = cutadd;
					cut2[1] = 0;
					count += 1;
					alignorcopy[1] = 'A';
					starttermcut1 = 1;
				}
			}
			else if( pairanch[count0].startj > getterminalmargin(pairanch[count0].starti, marginfac2) )
			{
//				alignorcopy[1] = 'A';
				reporterr( "check 2, because startj=%d > starti=%d -> %d (clus2=%d)\n", pairanch[count0].startj, pairanch[count0].starti, getterminalmargin(pairanch[count0].starti,marginfac2), clus2 );
				cutadd = pairanch[count0].startj - getterminalmargin( pairanch[count0].starti, marginfac2 );
				reporterr( "cutadd(2)=%d\n", cutadd );
				{
					cut1[1] = 0;
					cut2[1] = cutadd;
					count += 1;
					alignorcopy[1] = 'A';
					starttermcut2 = 1;
				}
			}
		}
#endif
#if 1
//		reporterr( "pairanch when uwagaki: %d:%d\n", pairanch[count0].starti, pairanch[count0].startj );
		cut1[count+1] = pairanch[count0].starti;
		cut2[count+1] = pairanch[count0].startj;
		alignorcopy[count+1] = 'c';
		count += 1;

#if 1
#define DIVIDEANCHOR 100
		if( pairanch[count0].endi - cut1[count] == pairanch[count0].endj - cut2[count] )
		while( pairanch[count0].endi+1 - cut1[count] > DIVIDEANCHOR && pairanch[count0].endj+1 - cut2[count] > DIVIDEANCHOR )
		{
			reporterr( "added an anchor, because the length is %d,%d > DIVIDEANCHOR.    \r", pairanch[count0].endi+1 - cut1[count], pairanch[count0].endj+1 - cut2[count] );
			cut1[count+1] = cut1[count] + DIVIDEANCHOR;
			cut2[count+1] = cut2[count] + DIVIDEANCHOR;
			alignorcopy[count+1] = 'c';
			count += 1;
		}
#endif

		cut1[count+1] = pairanch[count0].endi+1;
		cut2[count+1] = pairanch[count0].endj+1;
		alignorcopy[count+1] = 'a';
//		reporterr( "\n###cut1 at %d / %d\n", cut1[count+1], len1 );
//		reporterr( "###cut2 at %d / %d\n", cut2[count+1], len2 );
//		reporterr( "sa1=%d, sa2=%d\n", cut1[count+1]-cut1[count], cut2[count+1]-cut2[count] );
		count += 1;
		count0++;
	}
	reporterr( "\n" );

#if 1 // mattan no tansaku hanni wo seigen
	alignorcopy[count] = 'a';
//	if( count > 1 && (len1-cut1[count]) > (len2-cut2[count]) + 2*TERMINALSEGMENTLENGTH ) // 2 ha tekitou
	if( count > 1 && (len1-cut1[count]) > getterminalmargin(len2-cut2[count],marginfac1) )
	{
//		reporterr( "last\n" );
//		alignorcopy[count] = 'A'; // mae no wo uwagaki
		//reporterr( "insert one anchor to restrict terminal gap length, 1, cut1[count]=%d, cut2[count]=%d\n", cut1[count], cut2[count] );
		//alignorcopy[count] = 'A'; // mae no wo uwagaki
//		cut1[count+1] = cut1[count] + TERMINALSEGMENTLENGTH;
//		cut1[count+1] = cut1[count] + (len2-cut2[count]) + TERMINALSEGMENTLENGTH;

		cutadd = len1 - 1 - ( (len1-cut1[count]) - getterminalmargin(len2-cut2[count], marginfac1) ); // wakarinikuikedo

//		if( 1 || len1-1 - cutadd > TERMINALMARGIN(0) )
		{
			alignorcopy[count] = 'A'; // mae no wo uwagaki
			cut1[count+1] = cutadd;
			cut2[count+1] = len2;
			alignorcopy[count+1] = 'a';
			cut1[count+2] = len1;
			cut2[count+2] = len2;
			alignorcopy[count+2] = 'c'; // tsukawanai
			count += 1;
			endtermcut1 = 1;
		}
	}
//	else if( count > 1 && (len2-cut2[count]) > (len1-cut1[count])  + 2*TERMINALSEGMENTLENGTH ) // 2 ha tekitou
	else if( count > 1 && (len2-cut2[count]) > getterminalmargin(len1-cut1[count],marginfac2) )
	{
//		reporterr( "last\n" );
//		alignorcopy[count] = 'A'; // mae no wo uwagaki
		//reporterr( "insert one anchor to restrict terminal gap length, 2, cut1[count]=%d, cut2[count]=%d\n", cut1[count], cut2[count] );
		//alignorcopy[count] = 'A'; // mae no wo uwagaki
		cutadd = len2 - 1 - ( (len2-cut2[count]) - getterminalmargin(len1-cut1[count], marginfac2) );

//		if( 1 || len2-1 - cutadd > TERMINALMARGIN(0) ) // iranai
		{
			alignorcopy[count] = 'A'; // mae no wo uwagaki
			cut1[count+1] = len1;
			cut2[count+1] = cutadd;
			alignorcopy[count+1] = 'a';
			cut1[count+2] = len1;
			cut2[count+2] = len2;
			alignorcopy[count+2] = 'c'; // tsukawanai
			count += 1;
			endtermcut2 = 1;
		}
	}
#endif

	if( cut1[count] != len1 || cut2[count] != len2 )
	{
		cut1[count+1] = len1;
		cut2[count+1] = len2;
		alignorcopy[count+1] = 'c'; // tsukawanai kedo
		count += 1;
	}
	count += 1;



#if 0
	for( i=0; i<count; i++ )
	{
		reporterr( "cut1[%d],cut2[%d]=%d,%d, alignorcopy=%c\n", i, i, cut1[i], cut2[i], alignorcopy[i] );
	}
#endif


#else
		reporterr( "pairanch when uwagaki: %d:%d\n", pairanch[count0].starti, pairanch[count0].startj );
		cut1[count+1] = pairanch[count0].starti;
		cut2[count+1] = pairanch[count0].startj;
		reporterr( "###cut at %d and %d\n", cut1[count+1], cut1[count+2] );
		reporterr( "###cut at %d and %d\n", cut2[count+1], cut2[count+2] );
		count += 1;
		count0++;
	}
	cut1[count+1] = len1;
	cut2[count+1] = len2;
	count += 2;
#endif


	totallen = 0;
	for( j=0; j<clus1; j++ ) result1[j][0] = 0;
	for( j=0; j<clus2; j++ ) result2[j][0] = 0;
	totalscore = 0.0;
	*fftlog = -1;
//exit( 1 );
	for( i=0; i<count-1; i++ )
	{
//		reporterr( "\ni=%d / %d \n\n", i, count );
		*fftlog += 1;
		if( i == 0 || ( i == 1 && alignorcopy[1] == 'A' ) ) headgp = outgap; else headgp = 1;
		if( i == count-2 || ( i == count-3 && alignorcopy[count-3] == 'A' ) ) tailgp = outgap; else tailgp = 1;
		//reporterr( "i=%d, headgp=%d\n", i, headgp );
		//reporterr( "i=%d, tailgp=%d\n", i, tailgp );

//		for( j=0; j<clus1; j++ ) reporterr( "Cut seq 0-%d at %d\n", j, cut1[i] );
//		for( j=0; j<clus2; j++ ) reporterr( "Cut seq 1-%d at %d\n", j, cut2[i] );

#if 0 // hoka mo atode henkou suru
		if( cut1[i] )
		{
//			getkyokaigap( sgap1, seq1, cut1[i]-1, clus1 );
//			getkyokaigap( sgap2, seq2, cut2[i]-1, clus2 );
			getkyokaigap( sgap1, tmpres1, nlen-1, clus1 );
			getkyokaigap( sgap2, tmpres2, nlen-1, clus2 );
		}
		else
		{
			for( j=0; j<clus1; j++ ) sgap1[j] = 'o';
			for( j=0; j<clus2; j++ ) sgap2[j] = 'o';
		}
		if( cut1[i+1] != len1 )
		{       
			getkyokaigap( egap1, seq1, cut1[i+1], clus1 );
			getkyokaigap( egap2, seq2, cut2[i+1], clus2 );
		}       
		else    
		{       
			for( j=0; j<clus1; j++ ) egap1[j] = 'o';
			for( j=0; j<clus2; j++ ) egap2[j] = 'o';
		}
#else
		if( cut1[i] )
			getkyokaigap( sgap1, seq1, cut1[i]-1, clus1 );
		else
			for( j=0; j<clus1; j++ ) sgap1[j] = 'o';
		if( cut2[i] )
			getkyokaigap( sgap2, seq2, cut2[i]-1, clus2 );
		else
			for( j=0; j<clus2; j++ ) sgap2[j] = 'o';

		if( cut1[i+1] != len1 )
			getkyokaigap( egap1, seq1, cut1[i+1], clus1 );
		else    
			for( j=0; j<clus1; j++ ) egap1[j] = 'o';
		if( cut2[i+1] != len2 )
			getkyokaigap( egap2, seq2, cut2[i+1], clus2 );
		else    
			for( j=0; j<clus2; j++ ) egap2[j] = 'o';
#endif
#if DEBUG
		fprintf( stderr, "DP %03d / %03d %4d to ", i+1, count-1, totallen );
#else
#if 1
		if( 1 || fftkeika ) fprintf( stderr, "DP %05d / %05d \r", i+1, count-1 );
#endif
#endif
		for( j=0; j<clus1; j++ )
		{
//			if( cut1[i+1]-cut1[i] <= 0 )
//				fprintf( stderr, "### cut1[i+1]=%d, cut1[i]=%d\n", cut1[i+1], cut1[i] );
			strncpy( tmpres1[j], seq1[j]+cut1[i], cut1[i+1]-cut1[i] );
			tmpres1[j][cut1[i+1]-cut1[i]] = 0;
		}
		orilen1 = strlen( tmpres1[0] );
//		if( kobetsubunkatsu && fftkeika && i%2==0 ) commongappick( clus1, tmpres1 ); // for last-restricted alignment
		if( kobetsubunkatsu && fftkeika ) commongappick( clus1, tmpres1 ); //dvtditr に呼ばれたとき fftkeika=1
//		if( kobetsubunkatsu ) commongappick( clus1, tmpres1 );
		for( j=0; j<clus2; j++ )
		{
//			fprintf( stderr, "### cut2[i+1]-cut2[i] = %d\n", cut2[i+1]-cut2[i] );
//			if( cut2[i+1]-cut2[i] <= 0 )
//				fprintf( stderr, "### cut2[i+1]=%d, cut2[i]=%d\n", cut2[i+1], cut2[i] );
			strncpy( tmpres2[j], seq2[j]+cut2[i], cut2[i+1]-cut2[i] );
			tmpres2[j][cut2[i+1]-cut2[i]] = 0;
		}
		orilen2 = strlen( tmpres2[0] );
//		if( kobetsubunkatsu && fftkeika && i%2==0 ) commongappick( clus2, tmpres2 ); // for last-restricted alignment
		if( kobetsubunkatsu && fftkeika ) commongappick( clus2, tmpres2 ); //dvtditr に呼ばれたとき fftkeika=1
//		if( kobetsubunkatsu ) commongappick( clus2, tmpres2 );

		if( constraint )
		{
			fprintf( stderr, "Not supported\n" );
			exit( 1 );
		}
#if 0
		fprintf( stderr, "i=%d, before alignment", i );
		fprintf( stderr, "%4d\n", totallen );
		fprintf( stderr, "\n\n" );
		for( j=0; j<clus1; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres1[j] );
		}
		fprintf( stderr, "-------\n" );
		for( j=0; j<clus2; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres2[j] );
		}
#endif

#if 0
		fprintf( stdout, "writing input\n" );
		for( j=0; j<clus1; j++ )
		{
			fprintf( stdout, ">%d of GROUP1\n", j );
			fprintf( stdout, "%s\n", tmpres1[j] );
		}
		for( j=0; j<clus2; j++ )
		{
			fprintf( stdout, ">%d of GROUP2\n", j );
			fprintf( stdout, "%s\n", tmpres2[j] );
		}
		fflush( stdout );
#endif
//		reporterr( "i=%d, orilen1=%d, len1=%d, strlen(tmpseq1[0])=%d\n", i, orilen1, len1, strlen(tmpres1[0]) );
//		if( i%2 == 1 && orilen1==len1 && orilen1==orilen2 && orilen1==strlen( tmpres1[0] ) ) // zenchou itchi no toki nomi
//		if( 0 && i%2 == 1 && orilen1==orilen2 && orilen1==strlen( tmpres1[0] ) && !strcmp( tmpres1[0], tmpres2[0] ) ) // ato de fukkatsu saseru
		if( alignorcopy[i] == 'c'  && orilen1==orilen2 && orilen1==strlen( tmpres1[0] ) && !strcmp( tmpres1[0], tmpres2[0] ) ) // ato de fukkatsu saseru
		{
//			checklength = 1;
#if 0
			reporterr( "\ncopying\n" );
			for( j=0; j<clus1; j++ )
				reporterr( "tmpres1[j] = %s\n", tmpres1[j] );
			for( j=0; j<clus2; j++ )
				reporterr( "tmpres2[j] = %s\n", tmpres2[j] );
#endif
//			if( strlen( tmpres1[0] ) != strlen( tmpres2[0] ) )
//			{
//				reporterr( "Length differs? %d != %d\n", strlen( tmpres1[0] ), strlen( tmpres2[0] ) );
//				reporterr( "i=%d / %d\n", i, count-1 );
//				exit( 1 );
//			}
		}
		else
		{
//			reporterr( "\naligning %d x %d\n", strlen(tmpres1[0]), strlen(tmpres2[0]) );
			{
				//reporterr( "alignorcopy=%c, i=%d, seq1=%s\n", alignorcopy[i], i, tmpres1[0] );
				//reporterr( "alignorcopy=%c, i=%d, seq2=%s\n", alignorcopy[i], i, tmpres2[0] );
				if( alignorcopy[i] == 'A' ) // penalty_exx ..
//				if( 0 && alignorcopy[i] == 'A' ) // CHUUI! 'A' mukou
				{
					int penalty_exx;
					if( penalty_ex == 0 ) penalty_exx = -179; // == --exp 0.1 toriaezu
					else penalty_exx = penalty_ex;
					//reporterr( "i=%d, nomemsave, TERMINAL, penalty_exx = %d\n", i, penalty_exx );
					if( scoringmatrices ) // called by tditeration.c
						totalscore += A__align_variousdist( whichmtx, scoringmatrices, NULL, penalty, penalty_exx, tmpres1, tmpres2, eff1, eff2, eff1s, eff2s, clus1, clus2, alloclen, 0, NULL, sgap1, sgap2, egap1, egap2, NULL, 0, NULL, headgp, tailgp );
					else
						totalscore += A__align( n_dynamicmtx, penalty, penalty_exx, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, 0, NULL, sgap1, sgap2, egap1, egap2, NULL, 0, NULL, headgp, tailgp, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
					//use_getrusage();
				}
				else if( alg == 'A' )
//				else if( alg == 'A' || alignorcopy[i] == 'A' ) // hitomazu
				{
					if( scoringmatrices ) // called by tditeration.c
						totalscore += A__align_variousdist( whichmtx, scoringmatrices, NULL, penalty, penalty_ex, tmpres1, tmpres2, eff1, eff2, eff1s, eff2s, clus1, clus2, alloclen, 0, NULL, sgap1, sgap2, egap1, egap2, NULL, 0, NULL, headgp, tailgp );
					else
						totalscore += A__align( n_dynamicmtx, penalty, penalty_ex, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, 0, NULL, sgap1, sgap2, egap1, egap2, NULL, 0, NULL, headgp, tailgp, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
				}
				else if( alg == 'M' )
				{
					if( scoringmatrices ) // called by tditeration.c
						totalscore += MSalignmm_variousdist( NULL, scoringmatrices, NULL, tmpres1, tmpres2, eff1, eff2, eff1s, eff2s, clus1, clus2, alloclen, sgap1, sgap2, egap1, egap2, NULL, 0, NULL, headgp, tailgp );
					else
						totalscore += MSalignmm( n_dynamicmtx, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, sgap1, sgap2, egap1, egap2, NULL, 0, NULL, headgp, tailgp, NULL, NULL, NULL, 0.0, 0.0 );
				}
				else
				{
					fprintf( stderr, "alg = %c\n", alg );
					ErrorExit( "ERROR IN SOURCE FILE Falign.c" );
				}
			}
		}

		if( i == 1 && alignorcopy[i] == 'A' )
		{
			if( starttermcut1 )
			{
				for( j=0; j<clus2; j++ ) if( tmpres2[j][0] != '-' ) break;
				if( j<clus2 ) reporterr( "There may be a problem at the 5' end (1). Please contact katoh@ifrec.osaka-u.ac.jp\n", tmpres2[j] );
			}
			else if( starttermcut2 )
			{
				for( j=0; j<clus1; j++ ) if( tmpres1[j][0] != '-' ) break;
				if( j<clus1 ) reporterr( "There may be a problem at the 5' end (2).  Please contact katoh@ifrec.osaka-u.ac.jp\n" );
			}
			else
				break;
		}
		else if( alignorcopy[i] == 'A' )
		{	
			int tmplen;
			if( endtermcut1 )
			{
				tmplen = strlen( tmpres2[0] );
				for( j=0; j<clus2; j++ ) if( tmpres2[j][tmplen-1] != '-' ) break;
				if( j<clus2 ) reporterr( "There may be a problem at the 3' end (3).  Please contact katoh@ifrec.osaka-u.ac.jp\n" );
			}
			else if( endtermcut2 )
			{
				tmplen = strlen( tmpres1[0] );
				for( j=0; j<clus1; j++ ) if( tmpres1[j][tmplen-1] != '-' ) break;
				if( j<clus1 ) reporterr( "There may be a problem at the 3' end (4).  Please contact katoh@ifrec.osaka-u.ac.jp\n" );
			}
			else
				break;
		}

		nlen = strlen( tmpres1[0] );

		if( totallen + nlen > alloclen )
		{
			fprintf( stderr, "totallen=%d +  nlen=%d > alloclen = %d\n", totallen, nlen, alloclen );
			ErrorExit( "LENGTH OVER in Falign\n " );
		}
#if 0
		for( j=0; j<clus1; j++ ) strcat( result1[j], tmpres1[j] );
		for( j=0; j<clus2; j++ ) strcat( result2[j], tmpres2[j] );
#else
		for( j=0; j<clus1; j++ ) strcat( result1[j]+totallen, tmpres1[j] );
		for( j=0; j<clus2; j++ ) strcat( result2[j]+totallen, tmpres2[j] );
#endif
		totallen += nlen;
#if 0
		fprintf( stderr, "i=%d", i );
		fprintf( stderr, "%4d\n", totallen );
		fprintf( stderr, "\n\n" );
		for( j=0; j<clus1; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres1[j] );
		}
		fprintf( stderr, "-------\n" );
		for( j=0; j<clus2; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres2[j] );
		}
#endif
	}
#if KEIKA
	fprintf( stderr, "DP ... done   \n" );
#endif

	for( j=0; j<clus1; j++ ) strcpy( seq1[j], result1[j] );
	for( j=0; j<clus2; j++ ) strcpy( seq2[j], result2[j] );
#if 0
	fprintf( stderr, "keika \n\n" );
	for( j=0; j<clus1; j++ ) 
	{
		fprintf( stderr, ">group1-%d\n%100.100s\n", j, result1[j] );
	}
	fprintf( stderr, "- - - - - - - - - - -\n" );
	for( j=0; j<clus2; j++ ) 
	{
		fprintf( stderr, ">group2-%d\n%100.100s\n", j, result2[j] );
	}
//	if( clus1 == 1 && clus2 == 5 ) exit( 1 );
#endif
	return( totalscore );
}





/*
sakujo wo kentou (2010/10/05)
*/
double Falign_udpari_long(
			  int **whichmtx, double ***scoringmatrices, 
			  double **n_dynamicmtx,
			  char  **seq1, char  **seq2, 
			  double *eff1, double *eff2, 
			  double **eff1s, double **eff2s,
			  int    clus1, int    clus2,
			  int alloclen, int *fftlog )
{
	int i, j, k, l, m, maxk;
	int nlen, nlen2, nlen4;
	static TLS int crossscoresize = 0;
	char **tmpseq1 = NULL;
	char **tmpseq2 = NULL;
	char **tmpptr1 = NULL;
	char **tmpptr2 = NULL;
	char **tmpres1 = NULL;
	char **tmpres2 = NULL;
	char **result1 = NULL;
	char **result2 = NULL;
#if RND
	char **rndseq1 = NULL;
	char **rndseq2 = NULL;
#endif
	static TLS Fukusosuu **seqVector1 = NULL;
	static TLS Fukusosuu **seqVector2 = NULL;
	static TLS Fukusosuu **naiseki = NULL;   
	static TLS Fukusosuu *naisekiNoWa = NULL; 
	static TLS double *soukan = NULL;
	static TLS double **crossscore = NULL;
	int nlentmp;
	static TLS int *kouho = NULL;
	static TLS Segment *segment = NULL;
	static TLS Segment *segment1 = NULL;
	static TLS Segment *segment2 = NULL;
	static TLS Segment **sortedseg1 = NULL;
	static TLS Segment **sortedseg2 = NULL;
	static TLS int *cut1 = NULL;
	static TLS int *cut2 = NULL;
	char *sgap1, *egap1, *sgap2, *egap2;
	static TLS int localalloclen = 0;
	int lag;
	int tmpint;
	int count, count0;
	int len1, len2;
	int totallen;
	double totalscore;
	int nkouho = 0;
	int headgp, tailgp;
//	double dumfl = 0.0;

	if( seq1 == NULL )
	{
		if( kouho ) 
		{
//			fprintf( stderr, "### Freeing localarrays in Falign\n" );
			localalloclen = 0;
			crossscoresize = 0;
			mymergesort( 0, 0, NULL );
			alignableReagion( 0, 0, NULL, NULL, NULL, NULL, NULL );
			fft( 0, NULL, 1 );
			A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
			A__align_variousdist( NULL, NULL, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
			D__align_variousdist( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
			G__align11( NULL, NULL, NULL, 0, 0, 0 );
			blockAlign2( NULL, NULL, NULL, NULL, NULL, NULL );
			if( crossscore ) FreeDoubleMtx( crossscore );
			crossscore = NULL; // reallocate sareru kanousei ga arunode.
			free( kouho );
			kouho = NULL;
			free( cut1 );
			free( cut2 );
			free( segment );
			free( segment1 );
			free( segment2 );
			free( sortedseg1 );
			free( sortedseg2 );
			if( !kobetsubunkatsu )
			{
				FreeFukusosuuMtx ( seqVector1 );
				FreeFukusosuuMtx ( seqVector2 );
				FreeFukusosuuVec( naisekiNoWa );
				FreeFukusosuuMtx( naiseki );
				FreeDoubleVec( soukan );
			}
		}
		else
		{
//			fprintf( stderr, "Did not allocate localarrays in Falign\n" );
		}

		return( 0.0 );
	}

	len1 = strlen( seq1[0] );
	len2 = strlen( seq2[0] );
	nlentmp = MAX( len1, len2 );

	nlen = 1;
	while( nlentmp >= nlen ) nlen <<= 1;
#if 0
	fprintf( stderr, "###   nlen    = %d\n", nlen );
#endif

	nlen2 = nlen/2; nlen4 = nlen2 / 2;

#if 0
	fprintf( stderr, "len1 = %d, len2 = %d\n", len1, len2 );
	fprintf( stderr, "nlentmp = %d, nlen = %d\n", nlentmp, nlen );
#endif

	result1 = AllocateCharMtx( clus1, alloclen );
	result2 = AllocateCharMtx( clus2, alloclen );
	tmpres1 = AllocateCharMtx( clus1, alloclen );
	tmpres2 = AllocateCharMtx( clus2, alloclen );
	sgap1 = AllocateCharVec( clus1 );
	egap1 = AllocateCharVec( clus1 );
	sgap2 = AllocateCharVec( clus2 );
	egap2 = AllocateCharVec( clus2 );

	tmpseq1 = AllocateCharMtx( clus1, nlen );
	tmpseq2 = AllocateCharMtx( clus2, nlen );
	tmpptr1 = calloc( clus1, sizeof(char*) );
	tmpptr2 = calloc( clus2, sizeof(char*) );

#if RND
	rndseq1 = AllocateCharMtx( clus1, nlen );
	rndseq2 = AllocateCharMtx( clus2, nlen );
	for( i=0; i<clus1; i++ )
		generateRndSeq( rndseq1[i], nlen );
	for( i=0; i<clus2; i++ )
		generateRndSeq( rndseq2[i], nlen );
#endif

	if( !localalloclen )
	{
		kouho = AllocateIntVec( NKOUHO_LONG );
		cut1 = AllocateIntVec( MAXSEG );
		cut2 = AllocateIntVec( MAXSEG );
		segment = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		segment1 = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		segment2 = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
		sortedseg1 = (Segment **)calloc( MAXSEG, sizeof( Segment * ) );
		sortedseg2 = (Segment **)calloc( MAXSEG, sizeof( Segment * ) );
		if( !( segment && segment1 && segment2 && sortedseg1 && sortedseg2 ) )
			ErrorExit( "Allocation error\n" );

		if     ( scoremtx == -1 ) n20or4or2 = 1;
		else if( fftscore )       n20or4or2 = 1;
		else                      n20or4or2 = 20;
	}

	if( localalloclen < nlen )
	{
		if( localalloclen )
		{
#if 1
			if( !kobetsubunkatsu )
			{
				FreeFukusosuuMtx ( seqVector1 );
				FreeFukusosuuMtx ( seqVector2 );
				FreeFukusosuuVec( naisekiNoWa );
				FreeFukusosuuMtx( naiseki );
				FreeDoubleVec( soukan );
			}
#endif
		}


		if( !kobetsubunkatsu )
		{
			naisekiNoWa = AllocateFukusosuuVec( nlen );
			naiseki = AllocateFukusosuuMtx( n20or4or2, nlen );
			seqVector1 = AllocateFukusosuuMtx( n20or4or2, nlen+1 );
			seqVector2 = AllocateFukusosuuMtx( n20or4or2, nlen+1 );
			soukan = AllocateDoubleVec( nlen+1 );
		}
		localalloclen = nlen;
	}
	
	for( j=0; j<clus1; j++ ) strcpy( tmpseq1[j], seq1[j] );
	for( j=0; j<clus2; j++ ) strcpy( tmpseq2[j], seq2[j] );

#if 0
fftfp = fopen( "input_of_Falign", "w" );
fprintf( fftfp, "nlen = %d\n", nlen );
fprintf( fftfp, "seq1: ( %d sequences ) \n", clus1 );
for( i=0; i<clus1; i++ )
	fprintf( fftfp, "%s\n", seq1[i] );
fprintf( fftfp, "seq2: ( %d sequences ) \n", clus2 );
for( i=0; i<clus2; i++ )
	fprintf( fftfp, "%s\n", seq2[i] );
fclose( fftfp );
system( "less input_of_Falign < /dev/tty > /dev/tty" );
#endif
	if( !kobetsubunkatsu )
	{
		if( fftkeika ) fprintf( stderr,  " FFT ... " );

		for( j=0; j<n20or4or2; j++ ) vec_init( seqVector1[j], nlen );
		if( scoremtx == -1 )
		{
			for( i=0; i<clus1; i++ )
				seq_vec_4( seqVector1[0], eff1[i], tmpseq1[i] );
		}
		else if( fftscore )
		{
			for( i=0; i<clus1; i++ )
			{
#if 0
				seq_vec_2( seqVector1[0], polarity, eff1[i], tmpseq1[i] );
				seq_vec_2( seqVector1[1], volume,   eff1[i], tmpseq1[i] );
#else
				seq_vec_5( seqVector1[0], polarity, volume, eff1[i], tmpseq1[i] );
#endif
			}
		}
		else
		{
			for( i=0; i<clus1; i++ )
				seq_vec_3( seqVector1, eff1[i], tmpseq1[i] );
		}
#if RND
		for( i=0; i<clus1; i++ )
		{
			vec_init2( seqVector1, rndseq1[i], eff1[i], len1, nlen );
		}
#endif
#if 0
fftfp = fopen( "seqVec", "w" );
fprintf( fftfp, "before transform\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "nlen=%d\n", nlen );
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector1[k][l].R, seqVector1[k][l].I );
}
fclose( fftfp );
system( "less seqVec < /dev/tty > /dev/tty" );
#endif

		for( j=0; j<n20or4or2; j++ ) vec_init( seqVector2[j], nlen );
		if( scoremtx == -1 )
		{
			for( i=0; i<clus2; i++ )
				seq_vec_4( seqVector2[0], eff2[i], tmpseq2[i] );
		}
		else if( fftscore )
		{
			for( i=0; i<clus2; i++ )
			{
#if 0
				seq_vec_2( seqVector2[0], polarity, eff2[i], tmpseq2[i] );
				seq_vec_2( seqVector2[1], volume,   eff2[i], tmpseq2[i] );
#else
				seq_vec_5( seqVector2[0], polarity, volume, eff2[i], tmpseq2[i] );
#endif
			}
		}
		else
		{
			for( i=0; i<clus2; i++ )
				seq_vec_3( seqVector2, eff2[i], tmpseq2[i] );
		}
#if RND
		for( i=0; i<clus2; i++ )
		{
			vec_init2( seqVector2, rndseq2[i], eff2[i], len2, nlen );
		}
#endif

#if 0
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "before fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

		for( j=0; j<n20or4or2; j++ )
		{
			fft( nlen, seqVector2[j], 0 );
			fft( nlen, seqVector1[j], 0 );
		}
#if 0
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "#after fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "#%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
	   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

		for( k=0; k<n20or4or2; k++ ) 
		{
			for( l=0; l<nlen; l++ ) 
				calcNaiseki( naiseki[k]+l, seqVector1[k]+l, seqVector2[k]+l );
		}
		for( l=0; l<nlen; l++ ) 
		{
			naisekiNoWa[l].R = 0.0;
			naisekiNoWa[l].I = 0.0;
			for( k=0; k<n20or4or2; k++ ) 
			{
				naisekiNoWa[l].R += naiseki[k][l].R;
				naisekiNoWa[l].I += naiseki[k][l].I;
			}
		}
	
#if 0
	fftfp = fopen( "naisekiNoWa", "w" );
	fprintf( fftfp, "#Before fft\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f %f\n", l, naisekiNoWa[l].R, naisekiNoWa[l].I ); 
	fclose( fftfp );
	system( "less naisekiNoWa < /dev/tty > /dev/tty " );
#endif

		fft( -nlen, naisekiNoWa, 0 );
	
		for( m=0; m<=nlen2; m++ ) 
			soukan[m] = naisekiNoWa[nlen2-m].R;
		for( m=nlen2+1; m<nlen; m++ ) 
			soukan[m] = naisekiNoWa[nlen+nlen2-m].R;

#if 0
	fftfp = fopen( "naisekiNoWa", "w" );
	fprintf( fftfp, "#After fft\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f\n", l, naisekiNoWa[l].R ); 
	fclose( fftfp );
	fftfp = fopen( "list.plot", "w"  );
	fprintf( fftfp, "plot 'naisekiNoWa'\npause -1" );
	fclose( fftfp );
	system( "/usr/bin/gnuplot list.plot &" );
#endif
#if 0
	fprintf( stderr, "soukan\n" );
	for( l=0; l<nlen; l++ )
		fprintf( stderr, "%d  %f\n", l-nlen2, soukan[l] ); 
#if 0
	fftfp = fopen( "list.plot", "w"  );
	fprintf( fftfp, "plot 'frt'\n pause +1" );
	fclose( fftfp );
	system( "/usr/bin/gnuplot list.plot" );
#endif
#endif


		nkouho = getKouho( kouho, NKOUHO_LONG, soukan, nlen );

#if 0
		for( i=0; i<nkouho; i++ )
		{
			fprintf( stderr, "kouho[%d] = %d\n", i, kouho[i] );
		}
#endif
	}

#if KEIKA
	fprintf( stderr, "Searching anchors ... " );
#endif
	count = 0;



#define CAND 0
#if CAND
	fftfp = fopen( "cand", "w" );
	fclose( fftfp );
#endif
	if( kobetsubunkatsu )
	{
		maxk = 1;
		kouho[0] = 0;
	}
	else
	{
		maxk = nkouho;
	}

	for( k=0; k<maxk; k++ ) 
	{
		lag = kouho[k];
		if( lag <= -len1 || len2 <= lag ) continue;
//		fprintf( stderr, "k=%d, lag=%d\n", k, lag );
		zurasu2( lag, clus1, clus2, seq1, seq2, tmpptr1, tmpptr2 );
#if CAND
		fftfp = fopen( "cand", "a" );
		fprintf( fftfp, ">Candidate No.%d lag = %d\n", k+1, lag );
		fprintf( fftfp, "%s\n", tmpptr1[0] );
		fprintf( fftfp, ">Candidate No.%d lag = %d\n", k+1, lag );
		fprintf( fftfp, "%s\n", tmpptr2[0] );
		fprintf( fftfp, ">\n", k+1, lag );
		fclose( fftfp );
#endif

//		fprintf( stderr, "lag = %d\n", lag );
		tmpint = alignableReagion( clus1, clus2, tmpptr1, tmpptr2, eff1, eff2, segment+count );
//		fprintf( stderr, "lag = %d, %d found\n", lag, tmpint );

//		if( lag == -50 ) exit( 1 );
		
		if( count+tmpint > MAXSEG -3 ) ErrorExit( "TOO MANY SEGMENTS.\n" );

//		fprintf( stderr, "##### k=%d / %d\n", k, maxk );
//		if( tmpint == 0 ) break; // 060430 iinoka ? // 090530 yameta
		while( tmpint-- > 0 )
		{
#if 0
			if( segment[count].end - segment[count].start < fftWinSize )
			{
				count++;
				continue;
			}
#endif
			if( lag > 0 )
			{
				segment1[count].start  = segment[count].start ;
				segment1[count].end    = segment[count].end   ;
				segment1[count].center = segment[count].center;
				segment1[count].score  = segment[count].score;

				segment2[count].start  = segment[count].start  + lag;
				segment2[count].end    = segment[count].end    + lag;
				segment2[count].center = segment[count].center + lag;
				segment2[count].score  = segment[count].score       ;
			}
			else
			{
				segment1[count].start  = segment[count].start  - lag;
				segment1[count].end    = segment[count].end    - lag;
				segment1[count].center = segment[count].center - lag;
				segment1[count].score  = segment[count].score       ;

				segment2[count].start  = segment[count].start ;
				segment2[count].end    = segment[count].end   ;
				segment2[count].center = segment[count].center;
				segment2[count].score  = segment[count].score ;
			}
#if 0
			fprintf( stderr, "##### k=%d / %d\n", k, maxk );
			fprintf( stderr, "anchor %d, score = %f\n", count, segment1[count].score );
			fprintf( stderr, "in 1 %d\n", segment1[count].center );
			fprintf( stderr, "in 2 %d\n", segment2[count].center );
#endif
			segment1[count].pair = &segment2[count];
			segment2[count].pair = &segment1[count];
			count++;
#if 0
			fprintf( stderr, "count=%d\n", count );
#endif
		}
	}
#if 1
	if( !kobetsubunkatsu )
		if( fftkeika ) fprintf( stderr, "done. (%d anchors) ", count );
#endif
	if( !count && fftNoAnchStop )
		ErrorExit( "Cannot detect anchor!" );
#if 0
	fprintf( stderr, "RESULT before sort:\n" );
	for( l=0; l<count+1; l++ )
	{
		fprintf( stderr, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( stderr, "%d score = %f\n", segment2[l].center, segment1[l].score );
	}
#endif

	for( i=0; i<count; i++ )
	{
		sortedseg1[i] = &segment1[i];
		sortedseg2[i] = &segment2[i];
	}
#if 0
	tmpsort( count, sortedseg1 ); 
	tmpsort( count, sortedseg2 ); 
	qsort( sortedseg1, count, sizeof( Segment * ), segcmp );
	qsort( sortedseg2, count, sizeof( Segment * ), segcmp );
#else
	mymergesort( 0, count-1, sortedseg1 ); 
	mymergesort( 0, count-1, sortedseg2 ); 
#endif
	for( i=0; i<count; i++ ) sortedseg1[i]->number = i;
	for( i=0; i<count; i++ ) sortedseg2[i]->number = i;



	if( kobetsubunkatsu )
	{
		for( i=0; i<count; i++ )
	    {
			cut1[i+1] = sortedseg1[i]->center;
			cut2[i+1] = sortedseg2[i]->center;
		}
		cut1[0] = 0;
		cut2[0] = 0;
		cut1[count+1] = len1;
		cut2[count+1] = len2;
		count += 2;
	}

	else
	{
		if( count < 5000 )
		{
			if( crossscoresize < count+2 )
			{
				crossscoresize = count+2;
#if 1
				if( fftkeika ) fprintf( stderr, "######allocating crossscore, size = %d\n", crossscoresize );
#endif
				if( crossscore ) FreeDoubleMtx( crossscore );
				crossscore = AllocateDoubleMtx( crossscoresize, crossscoresize );
			}
			for( i=0; i<count+2; i++ ) for( j=0; j<count+2; j++ )
				crossscore[i][j] = 0.0;
			for( i=0; i<count; i++ )
			{
				crossscore[segment1[i].number+1][segment1[i].pair->number+1] = segment1[i].score;
				cut1[i+1] = sortedseg1[i]->center;
				cut2[i+1] = sortedseg2[i]->center;
			}
	
#if 0
			fprintf( stderr, "AFTER SORT\n" );
			for( i=0; i<count+1; i++ ) fprintf( stderr, "%d, %d\n", cut1[i], cut2[i] );
			fprintf( stderr, "crossscore = \n" );
			for( i=0; i<count+1; i++ )
			{
				for( j=0; j<count+1; j++ )
					fprintf( stderr, "%.0f ", crossscore[i][j] );
				fprintf( stderr, "\n" );
			}
#endif

			crossscore[0][0] = 10000000.0;
			cut1[0] = 0; 
			cut2[0] = 0;
			crossscore[count+1][count+1] = 10000000.0;
			cut1[count+1] = len1;
			cut2[count+1] = len2;
			count += 2;
			count0 = count;
		
//			fprintf( stderr, "\n\n\ncalling blockAlign2\n\n\n\n" );
			blockAlign2( cut1, cut2, sortedseg1, sortedseg2, crossscore, &count );
	
//			if( count-count0 )
//				fprintf( stderr, "%d unused anchors\n", count0-count );
	
			if( !kobetsubunkatsu && fftkeika )
				fprintf( stderr, "%d anchors found\n", count );
			if( fftkeika )
			{
				if( count0 > count )
				{
#if 0
					fprintf( stderr, "\7 REPEAT!? \n" ); 
#else
					fprintf( stderr, "REPEAT!? \n" ); 
#endif
					if( fftRepeatStop ) exit( 1 );
				}
#if KEIKA
				else fprintf( stderr, "done\n" );
#endif
			}
		}


		else
		{
			fprintf( stderr, "\nMany anchors were found. The upper-level DP is skipped.\n\n" );

			cut1[0] = 0; 
			cut2[0] = 0;
			count0 = 0;
			for( i=0; i<count; i++ )
			{
//				fprintf( stderr, "i=%d, %d-%d ?\n", i, sortedseg1[i]->center, sortedseg1[i]->pair->center );
				if( sortedseg1[i]->center > cut1[count0]
				 && sortedseg1[i]->pair->center > cut2[count0] )
				{
					count0++;
					cut1[count0] = sortedseg1[i]->center;
					cut2[count0] = sortedseg1[i]->pair->center;
				}
				else
				{
					if( i && sortedseg1[i]->score > sortedseg1[i-1]->score )
					{
						if( sortedseg1[i]->center > cut1[count0-1]
						 && sortedseg1[i]->pair->center > cut2[count0-1] )
						{
							cut1[count0] = sortedseg1[i]->center;
							cut2[count0] = sortedseg1[i]->pair->center;
						}
						else
						{
//							count0--;
						}
					}
				}
			}
//			if( count-count0 )
//				fprintf( stderr, "%d anchors unused\n", count-count0 );
			cut1[count0+1] = len1;
			cut2[count0+1] = len2;
			count = count0 + 2;
			count0 = count;
	
		}
	}

//	exit( 0 );

#if 0
	fftfp = fopen( "fft", "a" );
	fprintf( fftfp, "RESULT after sort:\n" );
	for( l=0; l<count; l++ )
	{
		fprintf( fftfp, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( fftfp, "%d\n", segment2[l].center );
	}
	fclose( fftfp );
#endif

#if 0
	fprintf( stderr, "RESULT after blckalign:\n" );
	for( l=0; l<count+1; l++ )
	{
		fprintf( stderr, "cut : %d %d\n", cut1[l], cut2[l] );
	}
#endif

#if 0
	fprintf( trap_g, "Devided to %d segments\n", count-1 );
	fprintf( trap_g, "%d  %d forg\n", MIN( clus1, clus2 ), count-1 );
#endif

	totallen = 0;
	for( j=0; j<clus1; j++ ) result1[j][0] = 0;
	for( j=0; j<clus2; j++ ) result2[j][0] = 0;
	totalscore = 0.0;
	*fftlog = -1;
//	reporterr( "\nin Falign_udpari(), *fftlog = %d\n", *fftlog );
	for( i=0; i<count-1; i++ )
	{
		*fftlog += 1;
		if( i == 0 ) headgp = outgap; else headgp = 1;
		if( i == count-2 ) tailgp = outgap; else tailgp = 1;

#if 0
		if( cut1[i] )
		{
//			getkyokaigap( sgap1, seq1, cut1[i]-1, clus1 );
//			getkyokaigap( sgap2, seq2, cut2[i]-1, clus2 );
			getkyokaigap( sgap1, tmpres1, nlen-1, clus1 );
			getkyokaigap( sgap2, tmpres2, nlen-1, clus2 );
		}
		else
		{
			for( j=0; j<clus1; j++ ) sgap1[j] = 'o';
			for( j=0; j<clus2; j++ ) sgap2[j] = 'o';
		}
		if( cut1[i+1] != len1 )
		{       
			getkyokaigap( egap1, seq1, cut1[i+1], clus1 );
			getkyokaigap( egap2, seq2, cut2[i+1], clus2 );
		}       
		else    
		{       
			for( j=0; j<clus1; j++ ) egap1[j] = 'o';
			for( j=0; j<clus2; j++ ) egap2[j] = 'o';
		}
#else
		if( cut1[i] )
			getkyokaigap( sgap1, seq1, cut1[i]-1, clus1 );
		else
			for( j=0; j<clus1; j++ ) sgap1[j] = 'o';
		if( cut2[i] )
			getkyokaigap( sgap2, seq2, cut2[i]-1, clus2 );
		else
			for( j=0; j<clus2; j++ ) sgap2[j] = 'o';

		if( cut1[i+1] != len1 )
			getkyokaigap( egap1, seq1, cut1[i+1], clus1 );
		else    
			for( j=0; j<clus1; j++ ) egap1[j] = 'o';
		if( cut2[i+1] != len2 )
			getkyokaigap( egap2, seq2, cut2[i+1], clus2 );
		else    
			for( j=0; j<clus2; j++ ) egap2[j] = 'o';
#endif

#if DEBUG
		fprintf( stderr, "DP %03d / %03d %4d to ", i+1, count-1, totallen );
#else
#if 1
		if( 1 || fftkeika ) fprintf( stderr, "DP %05d / %05d \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", i+1, count-1 );
#endif
#endif
		for( j=0; j<clus1; j++ )
		{
			strncpy( tmpres1[j], seq1[j]+cut1[i], cut1[i+1]-cut1[i] );
			tmpres1[j][cut1[i+1]-cut1[i]] = 0;
		}
		if( kobetsubunkatsu && fftkeika ) commongappick( clus1, tmpres1 ); //dvtditr に呼ばれたとき fftkeika=1
//		if( kobetsubunkatsu ) commongappick( clus1, tmpres1 );
		for( j=0; j<clus2; j++ )
		{
//			fprintf( stderr, "### cut2[i+1]-cut2[i] = %d\n", cut2[i+1]-cut2[i] );
			if( cut2[i+1]-cut2[i] <= 0 )
				fprintf( stderr, "### cut2[i+1]=%d, cut2[i]=%d\n", cut2[i+1], cut2[i] );
			strncpy( tmpres2[j], seq2[j]+cut2[i], cut2[i+1]-cut2[i] );
			tmpres2[j][cut2[i+1]-cut2[i]] = 0;
		}
		if( kobetsubunkatsu && fftkeika ) commongappick( clus2, tmpres2 ); //dvtditr に呼ばれたとき fftkeika=1
//		if( kobetsubunkatsu ) commongappick( clus2, tmpres2 );

		if( constraint )
		{
			fprintf( stderr, "Not supported\n" );
			exit( 1 );
		}
#if 0
		fprintf( stderr, "i=%d, before alignment", i );
		fprintf( stderr, "%4d\n", totallen );
		fprintf( stderr, "\n\n" );
		for( j=0; j<clus1; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres1[j] );
		}
		fprintf( stderr, "-------\n" );
		for( j=0; j<clus2; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres2[j] );
		}
#endif

#if 0
		fprintf( stdout, "writing input\n" );
		for( j=0; j<clus1; j++ )
		{
			fprintf( stdout, ">%d of GROUP1\n", j );
			fprintf( stdout, "%s\n", tmpres1[j] );
		}
		for( j=0; j<clus2; j++ )
		{
			fprintf( stdout, ">%d of GROUP2\n", j );
			fprintf( stdout, "%s\n", tmpres2[j] );
		}
		fflush( stdout );
#endif
		switch( alg )
		{
			case( 'M' ):
					if( scoringmatrices ) // called by tditeration.c
						totalscore += MSalignmm_variousdist( NULL, scoringmatrices, NULL, tmpres1, tmpres2, eff1, eff2, eff1s, eff2s, clus1, clus2, alloclen, sgap1, sgap2, egap1, egap2, NULL, 0, NULL, headgp, tailgp );
					else
						totalscore += MSalignmm( n_dynamicmtx, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, sgap1, sgap2, egap1, egap2, NULL, 0, NULL, headgp, tailgp, NULL, NULL, NULL, 0.0, 0.0 );
//						totalscore += G__align11( n_dynamicmtx, tmpres1, tmpres2, alloclen, headgp, tailgp ); // CHUUI!!!
				break;
			default:
				fprintf( stderr, "alg = %c\n", alg );
				ErrorExit( "ERROR IN SOURCE FILE Falign.c" );
				break;
		}

		nlen = strlen( tmpres1[0] );
		if( totallen + nlen > alloclen )
		{
			fprintf( stderr, "totallen=%d +  nlen=%d > alloclen = %d\n", totallen, nlen, alloclen );
			ErrorExit( "LENGTH OVER in Falign\n " );
		}
		for( j=0; j<clus1; j++ ) strcat( result1[j], tmpres1[j] );
		for( j=0; j<clus2; j++ ) strcat( result2[j], tmpres2[j] );
		totallen += nlen;
#if 0
		fprintf( stderr, "i=%d", i );
		fprintf( stderr, "%4d\n", totallen );
		fprintf( stderr, "\n\n" );
		for( j=0; j<clus1; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres1[j] );
		}
		fprintf( stderr, "-------\n" );
		for( j=0; j<clus2; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres2[j] );
		}
#endif
	}


//	reporterr( "\nafter Falign_udpari(), *fftlog = %d\n", *fftlog );

#if KEIKA
	fprintf( stderr, "DP ... done   \n" );
#endif

	for( j=0; j<clus1; j++ ) strcpy( seq1[j], result1[j] );
	for( j=0; j<clus2; j++ ) strcpy( seq2[j], result2[j] );
#if 0
	for( j=0; j<clus1; j++ ) 
	{
		fprintf( stderr, "%s\n", result1[j] );
	}
	fprintf( stderr, "- - - - - - - - - - -\n" );
	for( j=0; j<clus2; j++ ) 
	{
		fprintf( stderr, "%s\n", result2[j] );
	}
#endif

	FreeCharMtx( result1 );
	FreeCharMtx( result2 );
	FreeCharMtx( tmpres1 );
	FreeCharMtx( tmpres2 );
	free( sgap1 );
	free( egap1 );
	free( sgap2 );
	free( egap2 );
	FreeCharMtx( tmpseq1 );
	FreeCharMtx( tmpseq2 );
	free( tmpptr1 );
	free( tmpptr2 );
#if RND
	FreeCharMtx( rndseq1 );
	FreeCharMtx( rndseq2 );
#endif
	return( totalscore );
}
