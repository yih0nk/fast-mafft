#include "mltaln.h"
#include "progtiming.h"
#ifdef PROGTIMING
extern int pt_nofastdp;
#define FASTDP_ON ( pt_nofastdp < 0 ? ( pt_nofastdp = ( getenv( "MAFFT_NOFASTDP" ) != NULL ) ) == 0 : !pt_nofastdp )
#else
#define FASTDP_ON 1
#endif
#include "dp.h"

#define DEBUG 0
#define XXXXXXX    0
#define USE_PENALTY_EX  1

#define TERMGAPFAC 0.0
#define TERMGAPFAC_EX 0.0

#if 1

#if 0
static void match_calc_mtx_codon( double **codonmtx, double **mtx, double *match, char **s1, char **s2, int **codonseq1, int **codonseq2, int i1, int lgth2 )  
// gstart, gend wo tsukatte coding ka kakuninn suru beki
{
	char *seq2 = s2[0];
	double *doubleptr = mtx[(unsigned char)s1[0][i1]];
	int codonid1, codonid2;
//	int codonid1p, codonid2p;
//	int codonid1pp, codonid2pp;
	double codonmatch;
//	int *pt0, *pt1, *pt2;
	int *pt0;

	codonid1 = codonseq1[0][i1];
//	codonid1p = codonseq1[1][i1];
//	codonid1pp = codonseq1[2][i1];
	pt0 = codonseq2[0];
//	pt1 = codonseq2[1];
//	pt2 = codonseq2[2];
//	posin2 = 0;
	while( lgth2-- )
	{
//		reporterr( "%c-%c -> %f\n", s1[0][i1], *seq2, doubleptr[(unsigned char)*seq2] );
		
// track ha at de kokoromiru
		codonid2 = *pt0++;
//		codonid2p = *pt1++;
//		codonid2pp = *pt2++;
//		posin2++;
//
//		reporterr( "%c%c%c(%d)-%c%c%c(%d)\n", s1[0][i1], s1[0][i1+1], s1[0][i1+2], codonid1, *seq2, *(seq2+1), *(seq2+2), codonid2 );
		codonmatch = 0.0;
		if( codonid1 > -1 && codonid2 > -1 ) codonmatch = codonmtx[codonid1][codonid2] * 600.0;
//		if( codonid1p > -1 && codonid2p > -1 ) codonmatch -= codonmtx[codonid1p][codonid2p] * 600;
//		if( codonid1pp > -1 && codonid2pp > -1 ) codonmatch -= codonmtx[codonid1pp][codonid2pp] * 600;
		*match++ = doubleptr[(unsigned char)*seq2++] + codonmatch;
	}
}
#else
static void match_calc_mtx_codon( double **codonmtx, double **mtx, double *match, char **s1, char **s2, int **codonseq1, int **codonseq2, int i1, int lgth2 )  
// gstart, gend wo tsukatte coding ka kakuninn suru beki
{
	char *seq2 = s2[0];
	double *doubleptr = mtx[(unsigned char)s1[0][i1]];
	int codonid1, codonid2;
	double codonmatch;
	int *pt0;

	codonid1 = codonseq1[0][i1];
	pt0 = codonseq2[0];
	while( lgth2-- )
	{
//		reporterr( "%c-%c -> %f\n", s1[0][i1], *seq2, doubleptr[(unsigned char)*seq2] );
		
		codonid2 = *pt0++;
//		reporterr( "%c%c%c(%d)-%c%c%c(%d)\n", s1[0][i1], s1[0][i1+1], s1[0][i1+2], codonid1, *seq2, *(seq2+1), *(seq2+2), codonid2 );
		codonmatch = 0.0;
		if( codonid1 > -1 && codonid2 > -1 ) codonmatch = codonmtx[codonid1][codonid2] * 600.0; // codonid1>-1 ha iranai.
		*match++ = doubleptr[(unsigned char)*seq2++] + codonmatch;
	}
}
#endif

static void match_calc_mtx( double **mtx, double *match, char **s1, char **s2, int i1, int lgth2 ) 
{
	char *seq2 = s2[0];
	double *doubleptr = mtx[(unsigned char)s1[0][i1]];

	while( lgth2-- )
		*match++ = doubleptr[(unsigned char)*seq2++];
}
#else
static void match_calc( double *match, char **s1, char **s2, int i1, int lgth2 )
{
	int j;

	for( j=0; j<lgth2; j++ )
		match[j] = amino_dis[(*s1)[i1]][(*s2)[j]];
}
#endif

static double Atracking( double *lasthorizontalw, double *lastverticalw, 
						char **seq1, char **seq2, 
                        char **mseq1, char **mseq2, 
                        int **ijp,
						int tailgp,
						int *warpis, int *warpjs, int warpbase )
{
	int i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k, limk;
//	char gap[] = "-";
	char *gap;
	gap = newgapstr;
	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );
	double wm, g;
	double fpenalty = (double)penalty;
	double fpenalty_ex = (double)penalty_ex;


#if 0
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif
 
    for( i=0; i<lgth1+1; i++ ) 
    {
        ijp[i][0] = i + 1;
    }
    for( j=0; j<lgth2+1; j++ ) 
    {
        ijp[0][j] = -( j + 1 );
    }


//	if( tailgp == 1 || ijp[lgth1][lgth2] >= warpbase )
	if( tailgp == 1 )
		;
	else
	{
#if 1
//		reporterr( "lastverticalw[lgth1-1] = %f\n", lastverticalw[lgth1-1] );
//		reporterr( "lasthorizontalw[lgth2-1] = %f\n", lasthorizontalw[lgth2-1] );
		wm = lasthorizontalw[lgth2-1] - 1.0; // lasthorizontalw[lgth2-1] yori kanarazu chiisai.
		for( j=lgth2-2; j>=0; j-- )
		{
			if( (g=lasthorizontalw[j]+ ( fpenalty * TERMGAPFAC + fpenalty_ex * (lgth2-1-j) * TERMGAPFAC_EX ) ) > wm )
			{
				wm = g;
				iin = lgth1-1; jin = j;
				ijp[lgth1][lgth2] = -( lgth2 - j );
			}
		}
		for( i=lgth1-2; i>=0; i-- )
		{
			if( ( g=lastverticalw[i]+ ( fpenalty * TERMGAPFAC + fpenalty_ex * (lgth1-1-i) * TERMGAPFAC_EX ) ) > wm )
			{
				wm = g;
				iin = i; jin = lgth2-1;
				ijp[lgth1][lgth2] = +( lgth1 - i );
			}
		}
		if( lasthorizontalw[lgth2-1] > wm )  // score ga onaji baai erabarenai
		{
			wm = lasthorizontalw[lgth2-1];
			iin = lgth1-1; jin = lgth2-1;
			ijp[lgth1][lgth2] = 0;
		}
#else
		wm = lastverticalw[0];
		for( i=0; i<lgth1; i++ )
		{
			if( lastverticalw[i] >= wm )
			{
				wm = lastverticalw[i];
				iin = i; jin = lgth2-1;
				ijp[lgth1][lgth2] = +( lgth1 - i );
			}
		}
		for( j=0; j<lgth2; j++ )
		{
			if( lasthorizontalw[j] >= wm )
			{
				wm = lasthorizontalw[j];
				iin = lgth1-1; jin = j;
				ijp[lgth1][lgth2] = -( lgth2 - j );
			}
		}
#endif
	}



	mseq1[0] += lgth1+lgth2;
	*mseq1[0] = 0;
	mseq2[0] += lgth1+lgth2;
	*mseq2[0] = 0;



	iin = lgth1; jin = lgth2;
	limk = lgth1+lgth2 + 1;
	for( k=0; k<limk; k++ ) 
	{
		if( ijp[iin][jin] >= warpbase )
		{
//			fprintf( stderr, "WARP!\n" );
			ifi = warpis[ijp[iin][jin]-warpbase]; 
			jfi = warpjs[ijp[iin][jin]-warpbase];
		}
		else if( ijp[iin][jin] < 0 ) 
		{
			ifi = iin-1; jfi = jin+ijp[iin][jin];
		}
		else if( ijp[iin][jin] > 0 )
		{
			ifi = iin-ijp[iin][jin]; jfi = jin-1;
		}
		else
		{
			ifi = iin-1; jfi = jin-1;
		}

		if( ifi == -warpbase && jfi == -warpbase )
		{
			l = iin;
			while( --l >= 0 ) 
			{
				*--mseq1[0] = seq1[0][l];
				*--mseq2[0] = *gap;
				k++;
			}
			l= jin;
			while( --l >= 0 )
			{
				*--mseq1[0] = *gap;
				*--mseq2[0] = seq2[0][l];
				k++;
			}
			break;
		}
		else
		{
			l = iin - ifi;
			while( --l > 0 ) 
			{
				*--mseq1[0] = seq1[0][ifi+l];
				*--mseq2[0] = *gap;
				k++;
			}
			l= jin - jfi;
			while( --l > 0 )
			{
				*--mseq1[0] = *gap;
				*--mseq2[0] = seq2[0][jfi+l];
				k++;
			}
		}
		if( iin <= 0 || jin <= 0 ) break;
		*--mseq1[0] = seq1[0][ifi];
		*--mseq2[0] = seq2[0][jfi];
		k++;
		iin = ifi; jin = jfi;
	}

//	fprintf( stderr, "%s\n", mseq1[0] );
//	fprintf( stderr, "%s\n", mseq2[0] );
	return( wm );
}


double G__align11psg( double **codonmtx, double **n_dynamicmtx, char **seq1, char **seq2, int alloclen, int headgp, int tailgp, double *gstart, double *gend )
{



//	int k;
	register int i, j;
	int lasti;                      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lastj;
	int lgth1, lgth2;
	int resultlen;
	double wm, wmo;   /* int ?????? */
	double g;
	double *currentw, *previousw;
	double fpenalty = (double)penalty;
	double fpenalty_shift = (double)penalty_shift;
	double fpenalty_tmp;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penalty_ex;
	double fpenalty_ex_i;
#endif
#if 1
	double *wtmp;
	int *ijppt;
	double *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	static TLS double mi = 0.0;
	static TLS double *m = NULL;
	static TLS int **ijp = NULL;
	static TLS int mpi = 0;
	static TLS int *mp = NULL;
	static TLS double *w1 = NULL;
	static TLS double *w2 = NULL;
	static TLS double *match = NULL;
	static TLS double *initverticalw = NULL;    /* kufuu sureba iranai */
	static TLS double *lastverticalw = NULL;    /* kufuu sureba iranai */
	static TLS char **mseq1 = NULL;
	static TLS char **mseq2 = NULL;
	static TLS char **mseq = NULL;
	static TLS int **intwork = NULL;
	static TLS double **doublework = NULL;
	static TLS int orlgth1 = 0, orlgth2 = 0;
	static TLS double **amino_dynamicmtx = NULL; // ??
	static TLS int **codonseq1 = NULL;
	static TLS int **codonseq2 = NULL;

	int *warpis = NULL;
	int *warpjs = NULL;
	int *warpi = NULL;
	int *warpj = NULL;
	int *prevwarpi = NULL;
	int *prevwarpj = NULL;
	double *wmrecords = NULL;
	double *prevwmrecords = NULL;
	int warpn = 0;
	int warpbase;
	double curm = 0.0;
	double *wmrecordspt, *wmrecords1pt, *prevwmrecordspt;
	int *warpipt, *warpjpt;

	if( seq1 == NULL )
	{
		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			orlgth1 = 0;
			orlgth2 = 0;
			if( mseq1 ) free( mseq1 ); mseq1 = NULL;
			if( mseq2 ) free( mseq2 ); mseq2 = NULL;
			if( w1 ) FreeFloatVec( w1 ); w1 = NULL;
			if( w2 ) FreeFloatVec( w2 ); w2 = NULL;
			if( match ) FreeFloatVec( match ); match = NULL;
			if( initverticalw ) FreeFloatVec( initverticalw ); initverticalw = NULL;
			if( lastverticalw ) FreeFloatVec( lastverticalw ); lastverticalw = NULL;

			if( m ) FreeFloatVec( m ); m = NULL;
			if( mp ) FreeIntVec( mp ); mp = NULL;

			if( mseq ) FreeCharMtx( mseq ); mseq = NULL;



			if( doublework ) FreeFloatMtx( doublework ); doublework = NULL;
			if( intwork ) FreeIntMtx( intwork ); intwork = NULL;

			if( amino_dynamicmtx ) FreeDoubleMtx( amino_dynamicmtx ); amino_dynamicmtx = NULL;

			if( codonseq1 ) FreeIntMtx( codonseq1 ); codonseq1 = NULL;
			if( codonseq2 ) FreeIntMtx( codonseq2 ); codonseq2 = NULL;
		}
		orlgth1 = 0;
		orlgth2 = 0;
		return( 0.0 );
	}



	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );



	warpbase = lgth1 + lgth2;
	warpis = NULL;
	warpjs = NULL;
	warpn = 0;
	if( trywarp )
	{
//		fprintf( stderr, "IN G__align11\n" );
		if( headgp == 0 || tailgp == 0 )
		{
			fprintf( stderr, "At present, headgp and tailgp must be 1.\n" );
			exit( 1 );
		}

		wmrecords = AllocateFloatVec( lgth2+1 );
		warpi = AllocateIntVec( lgth2+1 );
		warpj = AllocateIntVec( lgth2+1 );
		prevwmrecords = AllocateFloatVec( lgth2+1 );
		prevwarpi = AllocateIntVec( lgth2+1 );
		prevwarpj = AllocateIntVec( lgth2+1 );
		for( i=0; i<lgth2+1; i++ ) prevwmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) wmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwarpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) prevwarpj[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpj[i] = -warpbase;
	}

#if 0
	if( lgth1 <= 0 || lgth2 <= 0 )
	{
		fprintf( stderr, "WARNING (g11): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif

#if 1
	if( lgth1 == 0 && lgth2 == 0 )
		return( 0.0 );

	if( lgth1 == 0 )
	{
		seq1[0][lgth2] = 0;
		while( lgth2 ) seq1[0][--lgth2] = *newgapstr;
//		fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
		return( 0.0 );
	}

	if( lgth2 == 0 )
	{
		seq2[0][lgth1] = 0;
		while( lgth1 ) seq2[0][--lgth1] = *newgapstr;
//		fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
		return( 0.0 );
	}
#endif


	wm = 0.0;

	if( orlgth1 == 0 )
	{
		mseq1 = AllocateCharMtx( 2, 0 ); // 2020/Apr
		mseq2 = AllocateCharMtx( 2, 0 ); // 2020/Apr
	}



	if( lgth1 > orlgth1 || lgth2 > orlgth2 )
	{
		int ll1, ll2;

		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );

			FreeCharMtx( mseq );



			FreeFloatMtx( doublework );
			FreeIntMtx( intwork );
			FreeDoubleMtx( amino_dynamicmtx );

			if( codonseq1 ) FreeIntMtx( codonseq1 );
			if( codonseq2 ) FreeIntMtx( codonseq2 );
		}

		ll1 = MAX( (int)(1.3*lgth1), orlgth1 ) + 100;
		ll2 = MAX( (int)(1.3*lgth2), orlgth2 ) + 100;

#if DEBUG
		fprintf( stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2 );
#endif

		w1 = AllocateFloatVec( ll2+2 );
		w2 = AllocateFloatVec( ll2+2 );
		match = AllocateFloatVec( ll2+2 );

		initverticalw = AllocateFloatVec( ll1+2 );
		lastverticalw = AllocateFloatVec( ll1+2 );

		m = AllocateFloatVec( ll2+2 );
		mp = AllocateIntVec( ll2+2 );

		mseq = AllocateCharMtx( 2, ll1+ll2 ); // 2020/Apr


		doublework = AllocateFloatMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		intwork = AllocateIntMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		amino_dynamicmtx = AllocateDoubleMtx( 0x100, 0x100 );

		if( codonscore )
		{
			codonseq1 = AllocateIntMtx( 3, ll1 ); // Only codonseq1[0] is used at this point
			codonseq2 = AllocateIntMtx( 3, ll2 ); // Only codonseq2[0] is used at this point
		}

#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;
	}
    for( i=0; i<nalphabets; i++) for( j=0; j<nalphabets; j++ )
		amino_dynamicmtx[(unsigned char)amino[i]][(unsigned char)amino[j]] = (double)n_dynamicmtx[i][j];


	mseq1[0] = mseq[0];
	mseq2[0] = mseq[1];


	if( orlgth1 > commonAlloc1 || orlgth2 > commonAlloc2 )
	{
		int ll1, ll2;

		if( commonAlloc1 && commonAlloc2 )
		{
			FreeIntMtx( commonIP );
		}

		ll1 = MAX( orlgth1, commonAlloc1 );
		ll2 = MAX( orlgth2, commonAlloc2 );

#if DEBUG
		fprintf( stderr, "\n\ntrying to allocate %dx%d matrices ... ", ll1+1, ll2+1 );
#endif

		commonIP = AllocateIntMtx( ll1+10, ll2+10 );

#if DEBUG
		fprintf( stderr, "succeeded\n\n" );
#endif

		commonAlloc1 = ll1;
		commonAlloc2 = ll2;
	}
	ijp = commonIP;

	if( codonscore )
	{
// codonseq[1..2] are not used yet.
		for( i=0; i<3&&i<lgth1; i++ )
		{
				codonseq1[0][i] = -1;
				codonseq1[1][i] = -1;
				codonseq1[2][i] = -1;
		}
		for( i=0; i<3&&i<lgth2; i++ )
		{
				codonseq2[0][i] = -1;
				codonseq2[1][i] = -1;
				codonseq2[2][i] = -1;
		}
		for( i=3; i<lgth1; i++ )
		{
			codonseq1[2][i] = codon2id( seq1[0]+i );
			codonseq1[1][i] = codon2id( seq1[0]+i-1 );
			codonseq1[0][i] = codon2id( seq1[0]+i-2 );
		}
		for( i=3; i<lgth2; i++ )
		{
			codonseq2[2][i] = codon2id( seq2[0]+i );
			codonseq2[1][i] = codon2id( seq2[0]+i-1 );
			codonseq2[0][i] = codon2id( seq2[0]+i-2 );
		}
	}
	else
	{
		codonseq1 = NULL;
		codonseq2 = NULL;
	}

#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

	currentw = w1;
	previousw = w2;


	if( codonscore )
	{

//		for( i=0; i<5; i++ ) for( j=0; j<5; j++ ) reporterr( "score[%c][%c]=%f\n", amino[i], amino[j], amino_dynamicmtx[amino[i]][amino[j]] );
//		exit( 1 );
		match_calc_mtx( amino_dynamicmtx, initverticalw, seq2, seq1, 0, lgth1 );
//		if( !codonpos || ! (gstart[0] == 0.5 && gend[0] == 0.5) ) // consider codon when unknown or coding
//			match_calc_mtx_codon( codonmtx, amino_dynamicmtx, currentw, seq1, seq2, codonseq1, codonseq2, 0, lgth2 );
#if 1
		if( gstart[0] == 0.5 && gend[0] > 0.5 ) // only at 3rd position
			match_calc_mtx_codon( codonmtx, amino_dynamicmtx, currentw, seq1, seq2, codonseq1, codonseq2, 0, lgth2 );
		else
			match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, 0, lgth2 );
#else
		match_calc_mtx_codon( codonmtx, amino_dynamicmtx, currentw, seq1, seq2, codonseq1, codonseq2, 0, lgth2, gstart[0] == 0.5 && gend[0] > 0.5 );
#endif
	}
	else
	{
		match_calc_mtx( amino_dynamicmtx, initverticalw, seq2, seq1, 0, lgth1 );
		match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, 0, lgth2 );
	}

	if( headgp == 1 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += fpenalty;
#if USE_PENALTY_EX
//			initverticalw[i] += fpenalty_ex * i; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += fpenalty;
#if USE_PENALTY_EX
//			currentw[j] += fpenalty_ex * j; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
		}
	}
	else
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += fpenalty * TERMGAPFAC;
#if USE_PENALTY_EX // 2018/Apr/22
			initverticalw[i] += fpenalty_ex * i * TERMGAPFAC_EX;
//			reporterr( "added _ex\n" );
#endif
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += fpenalty * TERMGAPFAC;
#if USE_PENALTY_EX // 2018/Apr/22
			currentw[j] += fpenalty_ex * j * TERMGAPFAC_EX;
//			reporterr( "added _ex\n" );
#endif
		}
	}


	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1]; mp[j] = 0;
	}

	if( lgth2 == 0 )
		lastverticalw[0] = 0.0;               // lgth2==0 no toki error
	else
		lastverticalw[0] = currentw[lgth2-1]; // lgth2==0 no toki error

	if( tailgp ) lasti = lgth1+1; else lasti = lgth1;
	lastj = lgth2+1;

#if XXXXXXX
fprintf( stderr, "currentw = \n" );
for( i=0; i<lgth1+1; i++ )
{
	fprintf( stderr, "%5.2f ", currentw[i] );
}
fprintf( stderr, "\n" );
fprintf( stderr, "initverticalw = \n" );
for( i=0; i<lgth2+1; i++ )
{
	fprintf( stderr, "%5.2f ", initverticalw[i] );
}
fprintf( stderr, "\n" );
#endif


	for( i=1; i<lasti; i++ )
	{
//		reporterr( "seq1[0][%d]=%c, gstart[i]=%f, gend[i]=%f\n", i, seq1[0][i], gstart[i], gend[i] );
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];


//		if( codonscore && ( !codonpos || !( gstart[i] == 0.5 && gend[i] == 0.5 ) ) ) // unknown or coding
//			match_calc_mtx_codon( codonmtx, amino_dynamicmtx, currentw, seq1, seq2, codonseq1, codonseq2, i, lgth2 );
#if 1
		if( codonscore && ( gstart[i] == 0.5 && gend[i] > 0.5 ) ) // only at 3rd position in seq1
			match_calc_mtx_codon( codonmtx, amino_dynamicmtx, currentw, seq1, seq2, codonseq1, codonseq2, i, lgth2 );
		else
			match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, i, lgth2 );
#else
		if( codonscore ) // only at 3rd position in seq1
			match_calc_mtx_codon( codonmtx, amino_dynamicmtx, currentw, seq1, seq2, codonseq1, codonseq2, i, lgth2, gstart[i] == 0.5 && gend[i] > 0.5 );
		else
			match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, i, lgth2 );
#endif
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
		currentw[0] = initverticalw[i];

		mi = previousw[0]; mpi = 0;

		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;
		if( trywarp )
		{
			prevwmrecordspt = prevwmrecords;
			wmrecordspt = wmrecords+1;
			wmrecords1pt = wmrecords;
			warpipt = warpi + 1;
			warpjpt = warpj + 1;
		}

		if( i < lgth1 ) fpenalty_ex_i = fpenalty_ex;
		else fpenalty_ex_i = 0.0; // 2018/May/11

		for( j=1; j<lastj; j++ )
		{

			wm = *prept;
			*ijppt = 0;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=mi+fpenalty*gend[i]) > wm )
			{
				wm = g;
				*ijppt = -( j - mpi );
//				reporterr( "hit! jump from %d to %d: %c->%c\n", j, mpi, seq2[0][j], seq2[0][mpi] );
			}
			if( (g=*prept+fpenalty*gstart[i-1]) >= mi )
//			if( (g=*prept) > mi )
			{
				mi = g;
				mpi = j-1;
//				reporterr( "hit! jump to %d: ->%c:%c\n", mpi, seq1[0][i-1], seq2[0][mpi] );
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex_i;
//			mi += fpenalty_ex;
#endif

#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=*mjpt + fpenalty*gend[i]) > wm )
			{
				wm = g;
				*ijppt = +( i - *mpjpt );
			}
			if( (g=*prept+fpenalty*gstart[i-1]) >= *mjpt )
//			if( (g=*prept) > *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}
#if USE_PENALTY_EX
			if( j < lgth2 ) // 2018/May/11
				m[j] += fpenalty_ex;
#endif
#if 1
			if( trywarp )
			{
				fpenalty_tmp = fpenalty_shift + fpenalty_ex * ( i - prevwarpi[j-1] + j - prevwarpj[j-1] );
//				fprintf( stderr, "fpenalty_shift = %f\n", fpenalty_tmp );

//				fprintf( stderr, "\n\n\nwarp to %c-%c (%d-%d) from %c-%c (%d-%d) ? prevwmrecords[%d] = %f + %f <- wm = %f\n", seq1[0][prevwarpi[j-1]], seq2[0][prevwarpj[j-1]], prevwarpi[j-1], prevwarpj[j-1], seq1[0][i], seq2[0][j], i, j, j, prevwmrecords[j-1], fpenalty_tmp, wm );
//				if( (g=prevwmrecords[j-1] + fpenalty_shift )> wm )
				if( ( g=*prevwmrecordspt++ + fpenalty_tmp )> wm ) // naka ha osokute kamawanai
				{
//					fprintf( stderr, "Yes! Warp!! from %d-%d (%c-%c) to  %d-%d (%c-%c) fpenalty_tmp = %f! warpn = %d\n", i, j, seq1[0][i], seq2[0][j-1], prevwarpi[j-1], prevwarpj[j-1],seq1[0][prevwarpi[j-1]], seq2[0][prevwarpj[j-1]], fpenalty_tmp, warpn );
					if( warpn && prevwarpi[j-1] == warpis[warpn-1] && prevwarpj[j-1] == warpjs[warpn-1] )
					{
						*ijppt = warpbase + warpn - 1;
					}
					else
					{
						*ijppt = warpbase + warpn;
						warpis = realloc( warpis, sizeof(int) * ( warpn+1 ) );
						warpjs = realloc( warpjs, sizeof(int) * ( warpn+1 ) );
						warpis[warpn] = prevwarpi[j-1];
						warpjs[warpn] = prevwarpj[j-1];
						warpn++;
					}
					wm = g;
				}
				else
				{
				}

				curm = *curpt + wm;
	
//				fprintf( stderr, "###### curm = %f at %c-%c, i=%d, j=%d\n", curm, seq1[0][i], seq2[0][j], i, j ); 
	
//				fprintf( stderr, "copy from i, j-1? %f > %f?\n", wmrecords[j-1], curm );
//				if( wmrecords[j-1] > wmrecords[j] )
				if( *wmrecords1pt > *wmrecordspt )
				{
//					fprintf( stderr, "yes\n" );
//					wmrecords[j] = wmrecords[j-1];
					*wmrecordspt = *wmrecords1pt;
//					warpi[j] = warpi[j-1];
//					warpj[j] = warpj[j-1];
					*warpipt  = *(warpipt-1);
					*warpjpt  = *(warpjpt-1);
//					fprintf( stderr, "warpi[j]=%d, warpj[j]=%d wmrecords[j] = %f\n", warpi[j], warpj[j], wmrecords[j] );
				}
//				else
//				{
//					fprintf( stderr, "no\n" );
//				}
	
//				fprintf( stderr, " curm = %f at %c-%c\n", curm, seq1[0][i], seq2[0][j] ); 
//				fprintf( stderr, " wmrecords[%d] = %f\n", j, wmrecords[j] ); 
//				fprintf( stderr, "replace?\n" );
	
//				if( curm > wmrecords[j] )
				if( curm > *wmrecordspt )
				{
//					fprintf( stderr, "yes at %d-%d (%c-%c), replaced warp: warpi[j]=%d, warpj[j]=%d warpn=%d, wmrecords[j] = %f -> %f\n", i, j, seq1[0][i], seq2[0][j], i, j, warpn, wmrecords[j], curm );
//					wmrecords[j] = curm;
					*wmrecordspt = curm;
//					warpi[j] = i;
//					warpj[j] = j;
					*warpipt = i;
					*warpjpt = j;
				}
//				else
//				{
//					fprintf( stderr, "No! warpi[j]=%d, warpj[j]=%d wmrecords[j] = %f\n", warpi[j], warpj[j], wmrecords[j] );
//				}
//				fprintf( stderr, "%d-%d (%c-%c) curm = %5.0f, wmrecords[j]=%f\n", i, j, seq1[0][i], seq2[0][j], curm, wmrecords[j] );
				wmrecordspt++;
				wmrecords1pt++;
				warpipt++;
				warpjpt++;
			}
#endif

			*curpt++ += wm;
			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
		}
		lastverticalw[i] = currentw[lgth2-1]; // lgth2==0 no toki error

		if( trywarp )
		{
			fltncpy( prevwmrecords, wmrecords, lastj );
			intncpy( prevwarpi, warpi, lastj );
			intncpy( prevwarpj, warpj, lastj );
		}
	}


	if( trywarp )
	{
//		fprintf( stderr, "\nwm = %f\n", wm );
//		fprintf( stderr, "warpn = %d\n", warpn );
		free( wmrecords );
		free( prevwmrecords );
		free( warpi );
		free( warpj );
		free( prevwarpi );
		free( prevwarpj );
	}

	wmo = Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, tailgp, warpis, warpjs, warpbase );
	if( !tailgp ) wm = wmo;

//	reporterr( "wm (after tracking) = %f\n", wm );
	if( warpis ) free( warpis );
	if( warpjs ) free( warpjs );


	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}


	strcpy( seq1[0], mseq1[0] );
	strcpy( seq2[0], mseq2[0] );
#if 0
	fprintf( stderr, "\n" );
	fprintf( stderr, ">\n%s\n", mseq1[0] );
	fprintf( stderr, ">\n%s\n", mseq2[0] );
	fprintf( stderr, "wm = %f\n", wm );
#endif
	return( wm );
}
/* work/prog: fast evaluation of one DP row of G__align11 (trywarp == 0).  Same cells,
   expressions and comparison order as the original loop; the column recurrence runs in its
   own vectorizable pass, the row recurrence stays in registers, selects are branch-free. */
/* work/gccpar: C ? A : B for an int C in {0,1}, written as bit arithmetic so that gcc folds it
   into a straight-line select (csel/cmov) instead of keeping a data-dependent branch; see
   ADP_CELL in Salignmm.c.  The cells below compute the same values in the same order as v3's. */
#define ISEL( C, A, B ) ( (B) ^ ( ( (B) ^ (A) ) & -(C) ) )
/* C >= X ? C : X for finite doubles; on x86 spelled as an IEEE max (one maxsd), see Salignmm.c */
#if defined( __x86_64__ ) || defined( __i386__ )
#define DSEL_GE( C, X ) ( (X) > (C) ? (X) : (C) )
#else
#define DSEL_GE( C, X ) ( (C) >= (X) ? (C) : (X) )
#endif

static void Gdp_row_fast( int i, int lastj, int lgth2, double * restrict colw, int * restrict colij, double * restrict colg2,
	double * restrict mloc, int * restrict mploc, double * restrict curw, int * restrict ijpi,
	const double * restrict prew, double fpenalty, double fpenalty_ex, double fpenalty_ex_i,
	double *mip, int *mpip, double *wmp )
{
	int j, ijv;
	double mil = *mip;
	int mpil = *mpip;
	double wmv = *wmp;

	for( j=1; j<lastj; j++ )
	{
		double mj = mloc[j];
		int mpj = mploc[j];
		double pw = prew[j-1];
		int take = ( pw >= mj );
		colw[j] = mj + fpenalty;
		colij[j] = i - mpj;
		mloc[j] = DSEL_GE( pw, mj );
		mploc[j] = ISEL( take, i-1, mpj );
	}
#if USE_PENALTY_EX
	for( j=1; j<lastj && j<lgth2; j++ ) // 2018/May/11: not for j == lgth2
		mloc[j] += fpenalty_ex;
	if( fpenalty_ex_i == 0.0 && lastj > 1 )
	{
		/* See Adp_row_fast (Salignmm.c): with a +0.0 extension penalty, mi + 0.0 == mi bitwise once
		   mi is not -0.0 (true after column 1), so the add moves off the serial chain into a
		   vectorized pass over the candidates.  Column 1 is done exactly as before. */
		double * restrict pwc = colg2;
		double gg, pw;
		for( j=1; j<lastj; j++ )
			pwc[j] = prew[j-1] + fpenalty_ex_i;
		j = 1;
		pw = prew[j-1];
		wmv = pw;
		ijv = 0;
		gg = mil + fpenalty;
		if( gg > wmv ) { wmv = gg; ijv = -( j - mpil ); }
		if( pw >= mil ) { mil = pw; mpil = j-1; }
		mil += fpenalty_ex_i;
		{
			double cw = colw[j];
			int cij = colij[j];
			int take = ( cw > wmv );
			ijv = take ? cij : ijv;
			wmv = take ? cw : wmv;
		}
		curw[j] += wmv;
		ijpi[j] = ijv;
		for( j=2; j<lastj; j++ )
		{
			double pc = pwc[j], pwj = prew[j-1];
			int take2, th;
			gg = mil + fpenalty;
			th = ( gg > pwj );
			wmv = th ? gg : pwj;
			ijv = -( j - mpil ) & -th;
			take2 = ( pc >= mil );
			mil = DSEL_GE( pc, mil );
			mpil = ISEL( take2, j-1, mpil );
			{
				double cw = colw[j];
				int cij = colij[j];
				int take = ( cw > wmv );
				ijv = ISEL( take, cij, ijv );
				wmv = take ? cw : wmv;
			}
			curw[j] += wmv;
			ijpi[j] = ijv;
		}
		*mip = mil;
		*mpip = mpil;
		*wmp = wmv;
		return;
	}
#endif
	for( j=1; j<lastj; j++ )
	{
		double gg, pw = prew[j-1];
		int th;
		gg = mil + fpenalty;
		th = ( gg > pw );
		wmv = th ? gg : pw;
		ijv = -( j - mpil ) & -th;
		th = ( pw >= mil );
		mil = DSEL_GE( pw, mil );
		mpil = ISEL( th, j-1, mpil );
#if USE_PENALTY_EX
		mil += fpenalty_ex_i;
#endif
		{
			double cw = colw[j];
			int cij = colij[j];
			int take = ( cw > wmv );
			ijv = ISEL( take, cij, ijv );
			wmv = take ? cw : wmv;
		}
		curw[j] += wmv;
		ijpi[j] = ijv;
	}
	*mip = mil;
	*mpip = mpil;
	*wmp = wmv;
}

/* work/prog2: one DP row of G__align11 in a single pass, with the match score of the cell
   (match_calc_mtx) read inline instead of in a separate pass, for fpenalty_ex == 0.0.  Every cell
   performs the operations of Gdp_row_fast's fast path (which are those of the original loop, see
   there) in the same order; the column update, previously a separate loop, is independent across
   columns, and the only other change is that currentw[j] = match + wm is formed from the match
   value directly instead of from the value match_calc_mtx stored in currentw[j].  Column lgth2
   has no match value (match_calc_mtx fills 0..lgth2-1) and keeps the buffer's old value, as
   before. */
#define GDP_CELL( FIRST, HASMATCH ) \
	{ \
		double mj = mloc[j], pw = prew[j-1], base, cw, mnew, gg, wh; \
		int mpj = mploc[j], take, th, ijv; \
		if( HASMATCH ) base = mrow[s2[j]]; else base = curw[j]; \
		cw = mj + fpenalty; \
		take = ( pw >= mj ); \
		mnew = DSEL_GE( pw, mj ); \
		if( HASMATCH ) mnew += fpenalty_ex;  /* the original adds it for j < lgth2 only */ \
		mloc[j] = mnew; \
		mploc[j] = ISEL( take, i-1, mpj ); \
		gg = mil + fpenalty; \
		th = ( gg > pw ); \
		wh = th ? gg : pw; \
		ijv = -( j - mpil ) & -th; \
		if( FIRST ) \
		{ \
			if( pw >= mil ) { mil = pw; mpil = j-1; } \
			mil += fpenalty_ex_i; \
		} \
		else \
		{ \
			double pc = pw + fpenalty_ex_i; \
			int take2 = ( pc >= mil ); \
			mil = DSEL_GE( pc, mil ); \
			mpil = ISEL( take2, j-1, mpil ); \
		} \
		take = ( cw > wh ); \
		ijv = ISEL( take, i - mpj, ijv ); \
		wmv = take ? cw : wh; \
		curw[j] = base + wmv; \
		ijpi[j] = ijv; \
	}

/* work/prog3: GDP_CELL( 0, 1 ) with the match score given (prow[j] == mrow[s2[j]]) */
/* work/v7: with gccpar's branch-free selects (ISEL/DSEL_GE, the integer mask for the horizontal-gap
   candidate), i.e. exactly GDP_CELL( 0, 1 ) above with base = PEXPR.  Same operations on the same
   operands in the same order as prog3's version (cw, gg, pc, mnew + fpenalty_ex, base + wmv; the
   same >= / > comparisons and tie rules; -( j - mpil ) taken before mpil is updated), so gcc no
   longer has data-dependent branches in it (the tails of the vector rows, and the scalar reference
   of MAFFT_VECCHECK). */
#define GDP_CELL_P( PEXPR ) \
	{ \
		double mj = mloc[j], pw = prew[j-1], base, cw, mnew, gg, wh, pc; \
		int mpj = mploc[j], take, th, ijv, take2; \
		base = (PEXPR); \
		cw = mj + fpenalty; \
		take = ( pw >= mj ); \
		mnew = DSEL_GE( pw, mj ); \
		mnew += fpenalty_ex; \
		mloc[j] = mnew; \
		mploc[j] = ISEL( take, i-1, mpj ); \
		gg = mil + fpenalty; \
		th = ( gg > pw ); \
		wh = th ? gg : pw; \
		ijv = -( j - mpil ) & -th; \
		pc = pw + fpenalty_ex_i; \
		take2 = ( pc >= mil ); \
		mil = DSEL_GE( pc, mil ); \
		mpil = ISEL( take2, j-1, mpil ); \
		take = ( cw > wh ); \
		ijv = ISEL( take, i - mpj, ijv ); \
		wmv = take ? cw : wh; \
		curw[j] = base + wmv; \
		ijpi[j] = ijv; \
	}

static void Gdp_cells_scalar( int i, int ja, int jb, const double * restrict prow,
	double * restrict mloc, int * restrict mploc, double * restrict curw, int * restrict ijpi,
	const double * restrict prew, double fpenalty, double fpenalty_ex, double fpenalty_ex_i,
	double *mip, int *mpip, double *wmp )
{
	double mil = *mip, wmv = *wmp;
	int mpil = *mpip, j;
	for( j=ja; j<jb; j++ )
		GDP_CELL_P( prow[j] )
	*mip = mil;
	*mpip = mpil;
	*wmp = wmv;
}

#if defined( __aarch64__ ) || defined( __ARM_NEON ) || defined( __x86_64__ ) || defined( __i386__ )
#define GDPV_HAVE 1
#if ( defined( __x86_64__ ) || defined( __i386__ ) ) && !defined( __AVX2__ ) && ( defined( __GNUC__ ) || defined( __clang__ ) )
#define GDPV_X86DISPATCH 1
#define GDPV_ATTR4 __attribute__(( target( "avx2" ) ))
#else
#define GDPV_ATTR4
#endif
#define GDPV_NAME Gdp_cells_v2
#define GDPV_VB 16
#define GDPV_ATTR
#include "gdpvec_body.h"
#undef GDPV_NAME
#undef GDPV_VB
#undef GDPV_ATTR
#define GDPV_NAME Gdp_cells_v4
#define GDPV_VB 32
#define GDPV_ATTR GDPV_ATTR4
#include "gdpvec_body.h"
#undef GDPV_NAME
#undef GDPV_VB
#undef GDPV_ATTR
#endif

typedef void (*gdp_cells_fn)( int, int, int, const double *, double *, int *, double *, int *, const double *, double, double, double, double *, int *, double * );
static gdp_cells_fn gdp_cells = NULL;
static gdp_cells_fn gdp_cells_get( void ) // MAFFT_ADPVEC as in Salignmm.c (0 scalar, 2, 4; default: 4 on x86 with AVX2, 2 on other x86, scalar on ARM)
{
	gdp_cells_fn f = __atomic_load_n( &gdp_cells, __ATOMIC_ACQUIRE );
	if( f == NULL )
	{
		char *e = getenv( "MAFFT_ADPVEC" );
		int w = e ? atoi( e ) : -1;
		f = Gdp_cells_scalar;
#ifdef GDPV_HAVE
		if( w < 0 )
		{
#if defined( __x86_64__ ) || defined( __i386__ )
#ifdef GDPV_X86DISPATCH
			__builtin_cpu_init();
			w = __builtin_cpu_supports( "avx2" ) ? 4 : 2;
#else
			w = 4;
#endif
#else
			w = 0;
#endif
		}
#ifdef GDPV_X86DISPATCH
		if( w == 4 ) { __builtin_cpu_init(); if( !__builtin_cpu_supports( "avx2" ) ) w = 2; }
#endif
		if( w == 2 ) f = Gdp_cells_v2;
		else if( w == 4 ) f = Gdp_cells_v4;
#endif
		__atomic_store_n( &gdp_cells, f, __ATOMIC_RELEASE );
	}
	return( f );
}

#ifdef PROGTIMING
static int g_veccheck_v = -1;
static TLS double *gvc_m = NULL, *gvc_cur = NULL; static TLS int *gvc_mp = NULL, *gvc_ijp = NULL; static TLS int gvc_cap = 0;
static long gvc_ranges = 0, gvc_cells = 0;
static void gvc_report( void ) __attribute__(( destructor ));
static void gvc_report( void ) { if( g_veccheck_v > 0 ) reporterr( "\nVECCHECK(G): %ld cell ranges, %ld cells verified, no mismatch\n", gvc_ranges, gvc_cells ); }
#endif

static inline void Gdp_cells( int i, int ja, int jb, const double * restrict prow,
	double * restrict mloc, int * restrict mploc, double * restrict curw, int * restrict ijpi,
	const double * restrict prew, double fpenalty, double fpenalty_ex, double fpenalty_ex_i,
	double *mip, int *mpip, double *wmp )
{
	gdp_cells_fn f = gdp_cells_get();
#ifdef PROGTIMING
	if( g_veccheck_v < 0 ) g_veccheck_v = ( getenv( "MAFFT_VECCHECK" ) != NULL );
	if( g_veccheck_v && jb > ja ) // MAFFT_VECCHECK: the scalar kernel on shadow copies, compared bitwise
	{
		double mi2 = *mip, wm2 = *wmp; int mpi2 = *mpip;
		if( jb + 1 > gvc_cap )
		{
			gvc_cap = jb + 1 + jb / 2;
			gvc_m = realloc( gvc_m, sizeof( double ) * gvc_cap ); gvc_cur = realloc( gvc_cur, sizeof( double ) * gvc_cap );
			gvc_mp = realloc( gvc_mp, sizeof( int ) * gvc_cap ); gvc_ijp = realloc( gvc_ijp, sizeof( int ) * gvc_cap );
		}
		memcpy( gvc_m + ja, mloc + ja, sizeof( double ) * ( jb - ja ) );
		memcpy( gvc_mp + ja, mploc + ja, sizeof( int ) * ( jb - ja ) );
		Gdp_cells_scalar( i, ja, jb, prow, gvc_m, gvc_mp, gvc_cur, gvc_ijp, prew, fpenalty, fpenalty_ex, fpenalty_ex_i, &mi2, &mpi2, &wm2 );
		f( i, ja, jb, prow, mloc, mploc, curw, ijpi, prew, fpenalty, fpenalty_ex, fpenalty_ex_i, mip, mpip, wmp );
		if( memcmp( gvc_m + ja, mloc + ja, sizeof( double ) * ( jb - ja ) ) || memcmp( gvc_mp + ja, mploc + ja, sizeof( int ) * ( jb - ja ) )
		 || memcmp( gvc_cur + ja, curw + ja, sizeof( double ) * ( jb - ja ) ) || memcmp( gvc_ijp + ja, ijpi + ja, sizeof( int ) * ( jb - ja ) )
		 || memcmp( &mi2, mip, sizeof( double ) ) || mpi2 != *mpip || memcmp( &wm2, wmp, sizeof( double ) ) )
		{
			reporterr( "\nVECCHECK MISMATCH in G__align11 row %d, cells %d..%d\n", i, ja, jb-1 );
			exit( 1 );
		}
		__atomic_add_fetch( &gvc_ranges, 1, __ATOMIC_RELAXED );
		__atomic_add_fetch( &gvc_cells, (long)( jb - ja ), __ATOMIC_RELAXED );
		return;
	}
#endif
	f( i, ja, jb, prow, mloc, mploc, curw, ijpi, prew, fpenalty, fpenalty_ex, fpenalty_ex_i, mip, mpip, wmp );
}

static void Gdp_row_fused( int i, int lgth2, const double * restrict mrow, const unsigned char * restrict s2, const double * restrict prow,
	double * restrict mloc, int * restrict mploc, double * restrict curw, int * restrict ijpi,
	const double * restrict prew, double fpenalty, double fpenalty_ex, double fpenalty_ex_i,
	double *mip, int *mpip, double *wmp )
{
	// requires lgth2 >= 2 and fpenalty_ex_i == 0.0 (see Gdp_row_fast)
	double mil = *mip, wmv = *wmp;
	int mpil = *mpip, j;
	j = 1;
	GDP_CELL( 1, 1 )
	if( prow ) // work/prog3: the inner cells from the letter's precomputed match row, vectorized
	{
		if( lgth2 > 2 )
			Gdp_cells( i, 2, lgth2, prow, mloc, mploc, curw, ijpi, prew, fpenalty, fpenalty_ex, fpenalty_ex_i, &mil, &mpil, &wmv );
	}
	else
	for( j=2; j<lgth2; j++ )
		GDP_CELL( 0, 1 )
	j = lgth2;
	GDP_CELL( 0, 0 )
	*mip = mil;
	*mpip = mpil;
	*wmp = wmv;
}

#ifdef PROGTIMING
static TLS double *g_chk_m = NULL, *g_chk_cur = NULL;
static TLS int *g_chk_mp = NULL, *g_chk_ijp = NULL;
static TLS int g_chk_pending = 0, g_chk_mpi;
static TLS double g_chk_mi, g_chk_wm, g_chk_t;
static int g_dpcheck_v = -1;
static int g_dpcheck( void ) { if( g_dpcheck_v < 0 ) g_dpcheck_v = ( getenv( "MAFFT_DPCHECK" ) != NULL ); return( g_dpcheck_v ); }
#endif

double G__align11( double **n_dynamicmtx, char **seq1, char **seq2, int alloclen, int headgp, int tailgp )
{

//	int k;
	register int i, j;
	int lasti;                      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lastj;
	int lgth1, lgth2;
	int resultlen;
	double wm, wmo;   /* int ?????? */
	double g;
	double *currentw, *previousw;
	double fpenalty = (double)penalty;
	double fpenalty_shift = (double)penalty_shift;
	int gfused; // work/prog2
	static TLS double *gpbuf = NULL; static TLS size_t gpcap = 0; // work/prog3
	int gpidx[0x100], gpstride = 0;
	double *gprow = NULL;
	double fpenalty_tmp;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penalty_ex;
	double fpenalty_ex_i;
#endif
#if 1
	double *wtmp;
	int *ijppt;
	double *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	static TLS double mi = 0.0;
	static TLS double *m = NULL;
	static TLS int **ijp = NULL;
	static TLS int mpi = 0;
	static TLS int *mp = NULL;
	static TLS double *w1 = NULL;
	static TLS double *w2 = NULL;
	static TLS double *match = NULL;
	static TLS double *initverticalw = NULL;    /* kufuu sureba iranai */
	static TLS double *lastverticalw = NULL;    /* kufuu sureba iranai */
	static TLS char **mseq1 = NULL;
	static TLS char **mseq2 = NULL;
	static TLS char **mseq = NULL;
	static TLS int **intwork = NULL;
	static TLS double **doublework = NULL;
	static TLS int orlgth1 = 0, orlgth2 = 0;
	static TLS double **amino_dynamicmtx = NULL; // ??
	static TLS double *dpcolw = NULL; // fast DP path (work/prog)
	static TLS int *dpcolij = NULL;
	static TLS double *dpcolg2 = NULL;
	static TLS int dpcolalloc = 0;

	int *warpis = NULL;
	int *warpjs = NULL;
	int *warpi = NULL;
	int *warpj = NULL;
	int *prevwarpi = NULL;
	int *prevwarpj = NULL;
	double *wmrecords = NULL;
	double *prevwmrecords = NULL;
	int warpn = 0;
	int warpbase;
	double curm = 0.0;
	double *wmrecordspt, *wmrecords1pt, *prevwmrecordspt;
	int *warpipt, *warpjpt;


	if( seq1 == NULL )
	{
		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			orlgth1 = 0;
			orlgth2 = 0;
			if( mseq1 ) free( mseq1 ); mseq1 = NULL;
			if( mseq2 ) free( mseq2 ); mseq2 = NULL;
			if( w1 ) FreeFloatVec( w1 ); w1 = NULL;
			if( w2 ) FreeFloatVec( w2 ); w2 = NULL;
			if( match ) FreeFloatVec( match ); match = NULL;
			if( initverticalw ) FreeFloatVec( initverticalw ); initverticalw = NULL;
			if( lastverticalw ) FreeFloatVec( lastverticalw ); lastverticalw = NULL;

			if( m ) FreeFloatVec( m ); m = NULL;
			if( mp ) FreeIntVec( mp ); mp = NULL;

			if( mseq ) FreeCharMtx( mseq ); mseq = NULL;



			if( doublework ) FreeFloatMtx( doublework ); doublework = NULL;
			if( intwork ) FreeIntMtx( intwork ); intwork = NULL;

			if( amino_dynamicmtx ) FreeDoubleMtx( amino_dynamicmtx ); amino_dynamicmtx = NULL;
		}
		orlgth1 = 0;
		orlgth2 = 0;
		if( dpcolalloc )
		{
			free( dpcolw ); dpcolw = NULL;
			free( dpcolij ); dpcolij = NULL;
			free( dpcolg2 ); dpcolg2 = NULL;
			dpcolalloc = 0;
		}
		return( 0.0 );
	}



	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );

	warpbase = lgth1 + lgth2;
	warpis = NULL;
	warpjs = NULL;
	warpn = 0;
	if( trywarp )
	{
//		fprintf( stderr, "IN G__align11\n" );
		if( headgp == 0 || tailgp == 0 )
		{
			fprintf( stderr, "At present, headgp and tailgp must be 1.\n" );
			exit( 1 );
		}

		wmrecords = AllocateFloatVec( lgth2+1 );
		warpi = AllocateIntVec( lgth2+1 );
		warpj = AllocateIntVec( lgth2+1 );
		prevwmrecords = AllocateFloatVec( lgth2+1 );
		prevwarpi = AllocateIntVec( lgth2+1 );
		prevwarpj = AllocateIntVec( lgth2+1 );
		for( i=0; i<lgth2+1; i++ ) prevwmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) wmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwarpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) prevwarpj[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpj[i] = -warpbase;
	}

#if 0
	if( lgth1 <= 0 || lgth2 <= 0 )
	{
		fprintf( stderr, "WARNING (g11): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif

#if 1
	if( lgth1 == 0 && lgth2 == 0 )
		return( 0.0 );

	if( lgth1 == 0 )
	{
		seq1[0][lgth2] = 0;
		while( lgth2 ) seq1[0][--lgth2] = *newgapstr;
//		fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
		return( 0.0 );
	}

	if( lgth2 == 0 )
	{
		seq2[0][lgth1] = 0;
		while( lgth1 ) seq2[0][--lgth1] = *newgapstr;
//		fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
		return( 0.0 );
	}
#endif


	wm = 0.0;

	if( orlgth1 == 0 )
	{
		mseq1 = AllocateCharMtx( 2, 0 ); // 2020/Apr
		mseq2 = AllocateCharMtx( 2, 0 ); // 2020/Apr
	}



	if( lgth1 > orlgth1 || lgth2 > orlgth2 )
	{
		int ll1, ll2;

		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );

			FreeCharMtx( mseq );



			FreeFloatMtx( doublework );
			FreeIntMtx( intwork );
			FreeDoubleMtx( amino_dynamicmtx );
		}

		ll1 = MAX( (int)(1.3*lgth1), orlgth1 ) + 100;
		ll2 = MAX( (int)(1.3*lgth2), orlgth2 ) + 100;

#if DEBUG
		fprintf( stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2 );
#endif

		w1 = AllocateFloatVec( ll2+2 );
		w2 = AllocateFloatVec( ll2+2 );
		match = AllocateFloatVec( ll2+2 );

		initverticalw = AllocateFloatVec( ll1+2 );
		lastverticalw = AllocateFloatVec( ll1+2 );

		m = AllocateFloatVec( ll2+2 );
		mp = AllocateIntVec( ll2+2 );

		mseq = AllocateCharMtx( 2, ll1+ll2 ); // 2020/Apr


		doublework = AllocateFloatMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		intwork = AllocateIntMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		amino_dynamicmtx = AllocateDoubleMtx( 0x100, 0x100 );

#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;
	}
    for( i=0; i<nalphabets; i++) for( j=0; j<nalphabets; j++ )
		amino_dynamicmtx[(unsigned char)amino[i]][(unsigned char)amino[j]] = (double)n_dynamicmtx[i][j];


	mseq1[0] = mseq[0];
	mseq2[0] = mseq[1];


	if( orlgth1 > commonAlloc1 || orlgth2 > commonAlloc2 )
	{
		int ll1, ll2;

		if( commonAlloc1 && commonAlloc2 )
		{
			FreeIntMtx( commonIP );
		}

		ll1 = MAX( orlgth1, commonAlloc1 );
		ll2 = MAX( orlgth2, commonAlloc2 );

#if DEBUG
		fprintf( stderr, "\n\ntrying to allocate %dx%d matrices ... ", ll1+1, ll2+1 );
#endif

		commonIP = AllocateIntMtx( ll1+10, ll2+10 );

#if DEBUG
		fprintf( stderr, "succeeded\n\n" );
#endif

		commonAlloc1 = ll1;
		commonAlloc2 = ll2;
	}
	ijp = commonIP;


#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

	currentw = w1;
	previousw = w2;


	match_calc_mtx( amino_dynamicmtx, initverticalw, seq2, seq1, 0, lgth1 );
	match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, 0, lgth2 );

	if( headgp == 1 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += fpenalty;
#if USE_PENALTY_EX
//			initverticalw[i] += fpenalty_ex * i; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += fpenalty;
#if USE_PENALTY_EX
//			currentw[j] += fpenalty_ex * j; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
		}
	}
	else
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += fpenalty * TERMGAPFAC;
#if USE_PENALTY_EX // 2018/Apr/22
			initverticalw[i] += fpenalty_ex * i * TERMGAPFAC_EX;
//			reporterr( "added _ex\n" );
#endif
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += fpenalty * TERMGAPFAC;
#if USE_PENALTY_EX // 2018/Apr/22
			currentw[j] += fpenalty_ex * j * TERMGAPFAC_EX;
//			reporterr( "added _ex\n" );
#endif
		}
	}


	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1]; mp[j] = 0;
	}

	if( lgth2 == 0 )
		lastverticalw[0] = 0.0;               // lgth2==0 no toki error
	else
		lastverticalw[0] = currentw[lgth2-1]; // lgth2==0 no toki error

	if( tailgp ) lasti = lgth1+1; else lasti = lgth1;
	lastj = lgth2+1;

#if XXXXXXX
fprintf( stderr, "currentw = \n" );
for( i=0; i<lgth1+1; i++ )
{
	fprintf( stderr, "%5.2f ", currentw[i] );
}
fprintf( stderr, "\n" );
fprintf( stderr, "initverticalw = \n" );
for( i=0; i<lgth2+1; i++ )
{
	fprintf( stderr, "%5.2f ", initverticalw[i] );
}
fprintf( stderr, "\n" );
#endif

	if( !trywarp && dpcolalloc < lgth2+2 )
	{
		if( dpcolalloc ) { free( dpcolw ); free( dpcolij ); free( dpcolg2 ); }
		dpcolalloc = lgth2 + 2 + 100;
		dpcolw = AllocateFloatVec( dpcolalloc );
		dpcolij = AllocateIntVec( dpcolalloc );
		dpcolg2 = AllocateFloatVec( dpcolalloc );
	}
#ifdef PROGTIMING
	if( g_dpcheck() )
	{
		g_chk_m = realloc( g_chk_m, sizeof( double ) * ( lgth2+2 ) );
		g_chk_mp = realloc( g_chk_mp, sizeof( int ) * ( lgth2+2 ) );
		g_chk_cur = realloc( g_chk_cur, sizeof( double ) * ( lgth2+2 ) );
		g_chk_ijp = realloc( g_chk_ijp, sizeof( int ) * ( lgth2+2 ) );
	}
#endif

	gfused = ( !trywarp && FASTDP_ON && fpenalty_ex == 0.0 && lgth2 >= 2 ); // work/prog2: Gdp_row_fused
	/* work/prog3: for the vector kernels, the match row of every letter of seq1 against seq2
	   (gprow + gpidx[letter] * gpstride; the same values as amino_dynamicmtx[letter][seq2[j]]) */
	gprow = NULL;
	if( gfused && gdp_cells_get() != Gdp_cells_scalar && lgth2 > 2 )
	{
		int nl = 0, c;
		for( c=0; c<0x100; c++ ) gpidx[c] = -1;
		for( i=1; i<=lgth1; i++ ) { c = (unsigned char)seq1[0][i]; if( gpidx[c] < 0 ) gpidx[c] = nl++; } // rows 1..lasti-1 (lasti <= lgth1+1; seq1[0][lgth1] is the terminator, as in match_calc_mtx)
		gpstride = lgth2 + 4;
		if( (size_t)nl * gpstride > gpcap ) { free( gpbuf ); gpcap = (size_t)nl * gpstride + gpstride * 4; gpbuf = malloc( sizeof( double ) * gpcap ); }
		for( c=0; c<0x100; c++ )
			if( gpidx[c] >= 0 )
			{
				const double *mr = amino_dynamicmtx[c];
				double *pr = gpbuf + (size_t)gpidx[c] * gpstride;
				const unsigned char *s2 = (const unsigned char *)seq2[0];
				for( j=0; j<lgth2; j++ ) pr[j] = mr[s2[j]];
			}
		gprow = gpbuf;
	}
	for( i=1; i<lasti; i++ )
	{
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

#ifdef PROGTIMING
		if( !gfused || g_dpcheck() ) // the check mode also runs the original loop, which needs the match row
#else
		if( !gfused ) // the fused row kernel reads the match scores itself
#endif
		match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, i, lgth2 );
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
		currentw[0] = initverticalw[i];

		mi = previousw[0]; mpi = 0;

		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;
		if( trywarp )
		{
			prevwmrecordspt = prevwmrecords;
			wmrecordspt = wmrecords+1;
			wmrecords1pt = wmrecords;
			warpipt = warpi + 1;
			warpjpt = warpj + 1;
		}

		if( i < lgth1 ) fpenalty_ex_i = fpenalty_ex;
		else fpenalty_ex_i = 0.0; // 2018/May/11

		if( !trywarp && FASTDP_ON )
		{
#ifdef PROGTIMING
			if( g_dpcheck() )
			{
				double t0, t1;
				g_chk_mi = mi; g_chk_wm = wm; g_chk_mpi = mpi;
				memcpy( g_chk_m, m, sizeof( double ) * lastj );
				memcpy( g_chk_mp, mp, sizeof( int ) * lastj );
				memcpy( g_chk_cur, currentw, sizeof( double ) * lastj );
				t0 = pt_now();
				if( gfused )
					Gdp_row_fused( i, lgth2, amino_dynamicmtx[(unsigned char)seq1[0][i]], (const unsigned char *)seq2[0], gprow ? gprow + (size_t)gpidx[(unsigned char)seq1[0][i]] * gpstride : NULL, g_chk_m, g_chk_mp, g_chk_cur, g_chk_ijp, previousw, fpenalty, fpenalty_ex, fpenalty_ex_i, &g_chk_mi, &g_chk_mpi, &g_chk_wm );
				else
					Gdp_row_fast( i, lastj, lgth2, dpcolw, dpcolij, dpcolg2, g_chk_m, g_chk_mp, g_chk_cur, g_chk_ijp, previousw, fpenalty, fpenalty_ex, fpenalty_ex_i, &g_chk_mi, &g_chk_mpi, &g_chk_wm );
				t1 = pt_now();
				pt_acc[18] += t1 - t0;
				g_chk_pending = 1; g_chk_t = t1;
				goto originalloop;
			}
#endif
			if( gfused )
				Gdp_row_fused( i, lgth2, amino_dynamicmtx[(unsigned char)seq1[0][i]], (const unsigned char *)seq2[0], gprow ? gprow + (size_t)gpidx[(unsigned char)seq1[0][i]] * gpstride : NULL, m, mp, currentw, ijp[i], previousw, fpenalty, fpenalty_ex, fpenalty_ex_i, &mi, &mpi, &wm );
			else
				Gdp_row_fast( i, lastj, lgth2, dpcolw, dpcolij, dpcolg2, m, mp, currentw, ijp[i], previousw, fpenalty, fpenalty_ex, fpenalty_ex_i, &mi, &mpi, &wm );
			lastverticalw[i] = currentw[lgth2-1]; // lgth2==0 no toki error
			continue;
		}
#ifdef PROGTIMING
		originalloop:
#endif

		for( j=1; j<lastj; j++ )
		{

			wm = *prept;
			*ijppt = 0;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=mi+fpenalty) > wm )
			{
				wm = g;
				*ijppt = -( j - mpi );
			}
			if( (g=*prept) >= mi )
//			if( (g=*prept) > mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex_i;
//			mi += fpenalty_ex;
#endif

#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=*mjpt + fpenalty) > wm )
			{
				wm = g;
				*ijppt = +( i - *mpjpt );
			}
			if( (g=*prept) >= *mjpt )
//			if( (g=*prept) > *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}
#if USE_PENALTY_EX
			if( j < lgth2 ) // 2018/May/11
				m[j] += fpenalty_ex;
#endif
#if 1
			if( trywarp )
			{
				fpenalty_tmp = fpenalty_shift + fpenalty_ex * ( i - prevwarpi[j-1] + j - prevwarpj[j-1] );
//				fprintf( stderr, "fpenalty_shift = %f\n", fpenalty_tmp );

//				fprintf( stderr, "\n\n\nwarp to %c-%c (%d-%d) from %c-%c (%d-%d) ? prevwmrecords[%d] = %f + %f <- wm = %f\n", seq1[0][prevwarpi[j-1]], seq2[0][prevwarpj[j-1]], prevwarpi[j-1], prevwarpj[j-1], seq1[0][i], seq2[0][j], i, j, j, prevwmrecords[j-1], fpenalty_tmp, wm );
//				if( (g=prevwmrecords[j-1] + fpenalty_shift )> wm )
				if( ( g=*prevwmrecordspt++ + fpenalty_tmp )> wm ) // naka ha osokute kamawanai
				{
//					fprintf( stderr, "Yes! Warp!! from %d-%d (%c-%c) to  %d-%d (%c-%c) fpenalty_tmp = %f! warpn = %d\n", i, j, seq1[0][i], seq2[0][j-1], prevwarpi[j-1], prevwarpj[j-1],seq1[0][prevwarpi[j-1]], seq2[0][prevwarpj[j-1]], fpenalty_tmp, warpn );
					if( warpn && prevwarpi[j-1] == warpis[warpn-1] && prevwarpj[j-1] == warpjs[warpn-1] )
					{
						*ijppt = warpbase + warpn - 1;
					}
					else
					{
						*ijppt = warpbase + warpn;
						warpis = realloc( warpis, sizeof(int) * ( warpn+1 ) );
						warpjs = realloc( warpjs, sizeof(int) * ( warpn+1 ) );
						warpis[warpn] = prevwarpi[j-1];
						warpjs[warpn] = prevwarpj[j-1];
						warpn++;
					}
					wm = g;
				}
				else
				{
				}

				curm = *curpt + wm;
	
//				fprintf( stderr, "###### curm = %f at %c-%c, i=%d, j=%d\n", curm, seq1[0][i], seq2[0][j], i, j ); 
	
//				fprintf( stderr, "copy from i, j-1? %f > %f?\n", wmrecords[j-1], curm );
//				if( wmrecords[j-1] > wmrecords[j] )
				if( *wmrecords1pt > *wmrecordspt )
				{
//					fprintf( stderr, "yes\n" );
//					wmrecords[j] = wmrecords[j-1];
					*wmrecordspt = *wmrecords1pt;
//					warpi[j] = warpi[j-1];
//					warpj[j] = warpj[j-1];
					*warpipt  = *(warpipt-1);
					*warpjpt  = *(warpjpt-1);
//					fprintf( stderr, "warpi[j]=%d, warpj[j]=%d wmrecords[j] = %f\n", warpi[j], warpj[j], wmrecords[j] );
				}
//				else
//				{
//					fprintf( stderr, "no\n" );
//				}
	
//				fprintf( stderr, " curm = %f at %c-%c\n", curm, seq1[0][i], seq2[0][j] ); 
//				fprintf( stderr, " wmrecords[%d] = %f\n", j, wmrecords[j] ); 
//				fprintf( stderr, "replace?\n" );
	
//				if( curm > wmrecords[j] )
				if( curm > *wmrecordspt )
				{
//					fprintf( stderr, "yes at %d-%d (%c-%c), replaced warp: warpi[j]=%d, warpj[j]=%d warpn=%d, wmrecords[j] = %f -> %f\n", i, j, seq1[0][i], seq2[0][j], i, j, warpn, wmrecords[j], curm );
//					wmrecords[j] = curm;
					*wmrecordspt = curm;
//					warpi[j] = i;
//					warpj[j] = j;
					*warpipt = i;
					*warpjpt = j;
				}
//				else
//				{
//					fprintf( stderr, "No! warpi[j]=%d, warpj[j]=%d wmrecords[j] = %f\n", warpi[j], warpj[j], wmrecords[j] );
//				}
//				fprintf( stderr, "%d-%d (%c-%c) curm = %5.0f, wmrecords[j]=%f\n", i, j, seq1[0][i], seq2[0][j], curm, wmrecords[j] );
				wmrecordspt++;
				wmrecords1pt++;
				warpipt++;
				warpjpt++;
			}
#endif

			*curpt++ += wm;
			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
		}
#ifdef PROGTIMING
		if( g_chk_pending )
		{
			pt_acc[19] += pt_now() - g_chk_t;
			g_chk_pending = 0;
			if( memcmp( g_chk_m, m, sizeof( double ) * lastj ) || memcmp( g_chk_mp, mp, sizeof( int ) * lastj ) || memcmp( g_chk_cur, currentw, sizeof( double ) * lastj )
			 || memcmp( g_chk_ijp+1, ijp[i]+1, sizeof( int ) * ( lastj-1 ) ) || memcmp( &g_chk_mi, &mi, sizeof( double ) ) || g_chk_mpi != mpi || memcmp( &g_chk_wm, &wm, sizeof( double ) ) )
			{
				reporterr( "\nDPCHECK MISMATCH in G__align11 row %d (lgth1=%d, lgth2=%d)\n", i, lgth1, lgth2 );
				exit( 1 );
			}
			pt_acc[20] += 1.0;
		}
#endif
		lastverticalw[i] = currentw[lgth2-1]; // lgth2==0 no toki error

		if( trywarp )
		{
			fltncpy( prevwmrecords, wmrecords, lastj );
			intncpy( prevwarpi, warpi, lastj );
			intncpy( prevwarpj, warpj, lastj );
		}
	}


	if( trywarp )
	{
//		fprintf( stderr, "\nwm = %f\n", wm );
//		fprintf( stderr, "warpn = %d\n", warpn );
		free( wmrecords );
		free( prevwmrecords );
		free( warpi );
		free( warpj );
		free( prevwarpi );
		free( prevwarpj );
	}

	wmo = Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, tailgp, warpis, warpjs, warpbase );
	if( !tailgp ) wm = wmo;

//	reporterr( "wm (after tracking) = %f\n", wm );
	if( warpis ) free( warpis );
	if( warpjs ) free( warpjs );


	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}


	strcpy( seq1[0], mseq1[0] );
	strcpy( seq2[0], mseq2[0] );
#if 0
	fprintf( stderr, "\n" );
	fprintf( stderr, ">\n%s\n", mseq1[0] );
	fprintf( stderr, ">\n%s\n", mseq2[0] );
	fprintf( stderr, "wm = %f\n", wm );
#endif

	return( wm );
}

double G__align11_noalign( double **n_dynamicmtx, int penal, int penal_ex, char **seq1, char **seq2, int alloclen )
/* warp mitaiou */
{
//	int k;
	register int i, j;
	int lasti;                      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lgth1, lgth2;
//	int resultlen;
	double wm;   /* int ?????? */
	double g;
	double *currentw, *previousw;
	double fpenalty = (double)penal;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penal_ex;
	double fpenalty_ex_i;
#endif
#if 1
	double *wtmp;
	double *mjpt, *prept, *curpt;
//	int *mpjpt;
#endif
	static TLS double mi, *m;
	static TLS double *w1, *w2;
	static TLS double *match;
	static TLS double *initverticalw;    /* kufuu sureba iranai */
	static TLS double *lastverticalw;    /* kufuu sureba iranai */
	static TLS int **intwork;
	static TLS double **doublework;
	static TLS int orlgth1 = 0, orlgth2 = 0;
	static TLS double **amino_dynamicmtx;

	if( seq1 == NULL )
	{
		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			orlgth1 = 0;
			orlgth2 = 0;
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );
			free( m );


			FreeFloatMtx( doublework );
			FreeIntMtx( intwork );
			FreeDoubleMtx( amino_dynamicmtx );
		}
		return( 0.0 );
	}


	wm = 0.0;



	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );



#if 0
	if( lgth1 <= 0 || lgth2 <= 0 )
	{
		fprintf( stderr, "WARNING (g11): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif

	if( lgth1 > orlgth1 || lgth2 > orlgth2 )
	{
		int ll1, ll2;

		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );




			FreeFloatMtx( doublework );
			FreeIntMtx( intwork );

			FreeDoubleMtx( amino_dynamicmtx );
		}

		ll1 = MAX( (int)(1.3*lgth1), orlgth1 ) + 100;
		ll2 = MAX( (int)(1.3*lgth2), orlgth2 ) + 100;

#if DEBUG
		fprintf( stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2 );
#endif

		w1 = AllocateFloatVec( ll2+2 );
		w2 = AllocateFloatVec( ll2+2 );
		match = AllocateFloatVec( ll2+2 );

		initverticalw = AllocateFloatVec( ll1+2 );
		lastverticalw = AllocateFloatVec( ll1+2 );

		m = AllocateFloatVec( ll2+2 );



		doublework = AllocateFloatMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		intwork = AllocateIntMtx( nalphabets, MAX( ll1, ll2 )+2 ); 


//		amino_dynamicmtx = AllocateDoubleMtx( 0x80, 0x80 ); 
		amino_dynamicmtx = AllocateDoubleMtx( 0x100, 0x100 );  // 2017/Nov.  constants.c no 'charsize' wo global hensuu nishita houga yoi?
#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;
	}


    for( i=0; i<nalphabets; i++) for( j=0; j<nalphabets; j++ )
		amino_dynamicmtx[(int)amino[i]][(int)amino[j]] = (double)n_dynamicmtx[i][j];


#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

	currentw = w1;
	previousw = w2;


	match_calc_mtx( amino_dynamicmtx, initverticalw, seq2, seq1, 0, lgth1 );


	match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, 0, lgth2 );

	if( 1 ) // tsuneni outgap-1
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += fpenalty;
#if USE_PENALTY_EX // 2018/Apr/23
//			initverticalw[i] += fpenalty_ex * i; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += fpenalty;
#if USE_PENALTY_EX // 2018/Apr/23
//			currentw[j] += fpenalty_ex * j; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
		}
	}

	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1];
	}

	if( lgth2 == 0 )
		lastverticalw[0] = 0.0;               // lgth2==0 no toki error
	else
		lastverticalw[0] = currentw[lgth2-1]; // lgth2==0 no toki error

	if( 1 ) lasti = lgth1+1; else lasti = lgth1; // tsuneni outgap-1

#if XXXXXXX
fprintf( stderr, "currentw = \n" );
for( i=0; i<lgth1+1; i++ )
{
	fprintf( stderr, "%5.2f ", currentw[i] );
}
fprintf( stderr, "\n" );
fprintf( stderr, "initverticalw = \n" );
for( i=0; i<lgth2+1; i++ )
{
	fprintf( stderr, "%5.2f ", initverticalw[i] );
}
fprintf( stderr, "\n" );
#endif

	for( i=1; i<lasti; i++ )
	{
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

		match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, i, lgth2 );
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
		currentw[0] = initverticalw[i];

		mi = previousw[0];

		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;

		if( i < lgth1 ) fpenalty_ex_i = fpenalty_ex;
		else fpenalty_ex_i = 0.0; // 2018/May/11

		for( j=1; j<lgth2+1; j++ )
		{
			wm = *prept;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=mi+fpenalty) > wm )
			{
				wm = g;
			}
//			if( (g=*prept) >= mi )
			if( (g=*prept) > mi ) // onaji hazu
			{
				mi = g;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex_i;
#endif

#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=*mjpt + fpenalty) > wm )
			{
				wm = g;
			}
//			if( (g=*prept) >= *mjpt )
			if( (g=*prept) > *mjpt ) // onaji hazu
			{
				*mjpt = g;
			}
#if USE_PENALTY_EX
			if( j < lgth2 ) // 2018/May/11
				m[j] += fpenalty_ex;
#endif

#if 0
			fprintf( stderr, "%5.0f ", wm );
#endif
			*curpt++ += wm;
			mjpt++;
			prept++;
		}
		lastverticalw[i] = currentw[lgth2-1]; // lgth2==0 no toki error
	}

#if 0
	fprintf( stderr, "\n" );
	fprintf( stderr, ">\n%s\n", mseq1[0] );
	fprintf( stderr, ">\n%s\n", mseq2[0] );
	fprintf( stderr, "wm (noalign) = %f\n", wm );
#endif

	return( wm );
}
