#ifndef _GNU_SOURCE
#define _GNU_SOURCE // work/v8s: glibc declares posix_memalign (used since work/prog2), nanosleep only with a feature macro under -std=c99; gcc >= 14 rejects implicit declarations
#endif
#include "mltaln.h"
#include "progtiming.h"
#include "dphelp.h"
#include "lazygap.h"
#include <time.h>
#ifdef PROGTIMING
int pt_nofastdp = -1; /* A/B switch for timing builds: MAFFT_NOFASTDP=1 disables the fast DP path */
static int pt_fastdp_enabled( void ) { if( pt_nofastdp < 0 ) pt_nofastdp = ( getenv( "MAFFT_NOFASTDP" ) != NULL ); return( !pt_nofastdp ); }
#define FASTDP_ON ( pt_fastdp_enabled() )
#else
#define FASTDP_ON 1
#endif
#include "dp.h"

#define MACHIGAI 0
#define OUTGAP0TRY 0
#define DEBUG 0
#define XXXXXXX    0
#define USE_PENALTY_EX  1
#define FASTMATCHCALC 1
#define SLOW 0

#define TERMGAPFAC 0.0
#define TERMGAPFAC_EX 0.0


static TLS double **impmtx = NULL;
static TLS int impalloclen = 0;
double imp_match_out_sc( int i1, int j1 )
{
//	fprintf( stderr, "imp+match = %f\n", impmtx[i1][j1] * fastathreshold );
//	fprintf( stderr, "val = %f\n", impmtx[i1][j1] );
	return( impmtx[i1][j1] );
}

#if 0
static void imp_match_out_vead_gapmap( double *imp, int i1, int lgth2, int *gapmap2 )
{
#if FASTMATCHCALC
	double *pt = impmtx[i1];
	int *gapmappt = gapmap2;
	while( lgth2-- )
		*imp++ += pt[*gapmappt++];
#else
	int j;
	double *pt = impmtx[i1];
	for( j=0; j<lgth2; j++ )
		*imp++ += pt[gapmap2[j]];
#endif
}


static void imp_match_out_vead_tate_gapmap( double *imp, int j1, int lgth1, int *gapmap1 )
{
#if FASTMATCHCALC
	int *gapmappt = gapmap1;
	while( lgth1-- )
		*imp++ += impmtx[*gapmappt++][j1];
#else
	int i;
	for( i=0; i<lgth1; i++ )
		*imp++ += impmtx[gapmap1[i]][j1];
#endif
}
#endif

static void imp_match_out_vead( double *imp, int i1, int lgth2 )
{
#if FASTMATCHCALC 
	double *pt = impmtx[i1];
	while( lgth2-- )
		*imp++ += *pt++;
#else
	int j;
	double *pt = impmtx[i1];
	for( j=0; j<lgth2; j++ )
		*imp++ += pt[j];
#endif
}
static void imp_match_out_vead_tate( double *imp, int j1, int lgth1 )
{
	int i;
	for( i=0; i<lgth1; i++ )
		*imp++ += impmtx[i][j1];
}

void imp_rna( int nseq1, int nseq2, char **seq1, char **seq2, double *eff1, double *eff2, RNApair ***grouprna1, RNApair ***grouprna2, int *gapmap1, int *gapmap2, RNApair *pair )
{
	foldrna( nseq1, nseq2, seq1, seq2, eff1, eff2, grouprna1, grouprna2, impmtx, gapmap1, gapmap2, pair );
}

void imp_match_init_strict( double *imp, int clus1, int clus2, int lgth1, int lgth2, char **seq1, char **seq2, double *eff1, double *eff2, double *eff1_kozo, double *eff2_kozo, LocalHom ***localhom, char *swaplist, int forscore, int *orinum1, int *orinum2, int *uselh, int *seedinlh1, int *seedinlh2, int nodeid, int nfiles )
{
//	int i, j, k1, k2, tmpint, start1, start2, end1, end2;
//	double effij;
//	double effij_kozo;
//	double effijx;
//	char *pt, *pt1, *pt2;
//	static TLS char *nocount1 = NULL;
//	static TLS char *nocount2 = NULL;
//	LocalHom *tmpptr;

	if( seq1 == NULL )
	{
		if( impmtx ) FreeFloatMtx( impmtx );
		impmtx = NULL;
//		if( nocount1 ) free( nocount1 );
//		nocount1 = NULL;
//		if( nocount2 ) free( nocount2 );
//		nocount2 = NULL;
		
		return;
	}

	if( impalloclen < lgth1 + 2 || impalloclen < lgth2 + 2 )
	{
		if( impmtx ) FreeFloatMtx( impmtx );
//		if( nocount1 ) free( nocount1 );
//		if( nocount2 ) free( nocount2 );
		impalloclen = MAX( lgth1, lgth2 ) + 2;
		impmtx = AllocateFloatMtx( impalloclen, impalloclen );
//		nocount1 = AllocateCharVec( impalloclen );
//		nocount2 = AllocateCharVec( impalloclen );
	}

	if( nodeid == -1 )
		fillimp( impmtx, imp, clus1, clus2, lgth1, lgth2, seq1, seq2, eff1, eff2, eff1_kozo, eff2_kozo, localhom, swaplist, forscore, orinum1, orinum2 ); // uselh -> target -> localhomtable. seedinlh12 -> localhom ni haitteiru.
	else
		fillimp_file( impmtx, imp, clus1, clus2, lgth1, lgth2, seq1, seq2, eff1, eff2, eff1_kozo, eff2_kozo, localhom, swaplist, forscore, orinum1, orinum2, uselh, seedinlh1, seedinlh2, nodeid, nfiles );
}


static void match_calc_del( int **which, double ***matrices, double *match, int n1, char **seq1, double *eff1, int n2, char **seq2, double *eff2, int i1, int lgth2, int mid, int nmask, int *mask1, int *mask2 ) 
{
// osoi!
	int i, j, k, m;
	int c1, c2;
//	fprintf( stderr, "\nmatch_calc_dynamicmtx... %d", i1 );
//	fprintf( stderr, "\nseq1[0]=%s\n", seq1[0] );
//	fprintf( stderr, "\nseq2[0]=%s\n", seq2[0] );
//	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
//	{
//		if( flip ) reporterr( "in match_calc_slow, which[%d][%d] = %d\n", j, i, which[j][i] );
//		else       reporterr( "in match_calc_slow, which[%d][%d] = %d\n", i, j, which[i][j] );
//	}
	for( k=0; k<lgth2; k++ )
	{
		for( m=0; m<nmask; m++ )
		{
			i = mask1[m];
			j = mask2[m];
//			reporterr( "Deleting %d-%d (c=%d)\n", i, j, mid );
//			if( k==0 ) fprintf( stderr, "pairoffset[%d][%d] = %f\n", i, j, po );
			c1 = amino_n[(unsigned int)seq1[i][i1]];
			c2 = amino_n[(unsigned int)seq2[j][k]];
//			reporterr( "k=%d, c1=%d, c2=%d, seq1[i][i1]=%c, seq2[%d][%d]=%c\n", k, c1, c2, seq1[i][i1], j, k, seq2[j][k] );
			if( seq1[i][i1] == '-' || seq2[j][k] == '-' ) continue;
			if( c1 < 0 || c2 < 0 ) continue;
//			fprintf( stderr, "c1=%d, c2=%d\n", c1, c2 );
//			fprintf( stderr, "match[k] = %f -> ", match[k], mid );
			match[k] -= matrices[mid][c1][c2] * eff1[i] * eff2[j];
//			fprintf( stderr, "match[k] = %f (mid=%d)\n", match[k], mid );
		}
	}
//	fprintf( stderr, "done\n" );
	return;
}


#if SLOW
static void match_calc_slow( int **which, double ***matrices, double *match, int n1, char **seq1, double *eff1, int n2, char **seq2, double *eff2, int i1, int lgth2, double **doublework, int **intwork, int initialize, int flip ) 
{
// osoi!
	int i, j, k;
	int c1, c2;
	int mid;
//	fprintf( stderr, "\nmatch_calc_dynamicmtx... %d", i1 );
//	fprintf( stderr, "\nseq1[0]=%s\n", seq1[0] );
//	fprintf( stderr, "\nseq2[0]=%s\n", seq2[0] );
//	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
//	{
//		if( flip ) reporterr( "in match_calc_slow, which[%d][%d] = %d\n", j, i, which[j][i] );
//		else       reporterr( "in match_calc_slow, which[%d][%d] = %d\n", i, j, which[i][j] );
//	}
	for( k=0; k<lgth2; k++ )
	{
		match[k] = 0.0;
		for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
		{
			if( flip ) mid = which[j][i];
			else       mid = which[i][j];
//			if( k==0 ) fprintf( stderr, "pairoffset[%d][%d] = %f\n", i, j, po );
			c1 = amino_n[(unsigned int)seq1[i][i1]];
			c2 = amino_n[(unsigned int)seq2[j][k]];
			if( seq1[i][i1] == '-' || seq2[j][k] == '-' ) continue;
			if( c1 < 0 || c2 < 0 ) continue;
//			fprintf( stderr, "c1=%d, c2=%d\n", c1, c2 );
			if( flip ) 
				match[k] += matrices[mid][c1][c2] * eff1[i] * eff2[j];
			else
				match[k] += matrices[mid][c1][c2] * eff1[i] * eff2[j];
//			fprintf( stderr, "match[k] = %f (which=%d)\n", match[k], mid );
		}
	}
//	fprintf( stderr, "done\n" );
	return;
}
#endif

static void fillzero( double *s, int l )
{
	while( l-- ) *s++ = 0.0;
}

static void match_calc_add( double **scoreingmtx, double *match, double **cpmx1, double **cpmx2, int i1, int lgth2, double **doublework, int **intwork, int initialize )
{
#if FASTMATCHCALC
//	fprintf( stderr, "\nmatch_calc... %d", i1 );
	int j, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	double *matchpt, *cpmxpdpt, **cpmxpdptpt;
	int *cpmxpdnpt, **cpmxpdnptpt;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
	if( initialize )
	{
		int count = 0;
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][j] )
				{
					cpmxpd[j][count] = cpmx2[l][j];
					cpmxpdn[j][count] = l;
					count++;
				}
			}
			cpmxpdn[j][count] = -1;
		}
	}

	{
		for( l=0; l<nalphabets; l++ )
		{
			scarr[l] = 0.0;
			for( j=0; j<nalphabets; j++ )
//				scarr[l] += n_dis[j][l] * cpmx1[j][i1];
//				scarr[l] += n_dis_consweight_multi[j][l] * cpmx1[j][i1];
				scarr[l] += scoreingmtx[j][l] * cpmx1[j][i1];
		}
		matchpt = match;
		cpmxpdnptpt = cpmxpdn;
		cpmxpdptpt = cpmxpd;
		while( lgth2-- )
		{
//			*matchpt = 0.0;
			cpmxpdnpt = *cpmxpdnptpt++;
			cpmxpdpt = *cpmxpdptpt++;
			while( *cpmxpdnpt>-1 )
				*matchpt += scarr[*cpmxpdnpt++] * *cpmxpdpt++;
			matchpt++;
		} 
	}
	free( scarr );
//	fprintf( stderr, "done\n" );
#else
	int j, k, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
// simple
	if( initialize )
	{
		int count = 0;
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][j] )
				{
					cpmxpd[count][j] = cpmx2[l][j];
					cpmxpdn[count][j] = l;
					count++;
				}
			}
			cpmxpdn[count][j] = -1;
		}
	}
	for( l=0; l<nalphabets; l++ )
	{
		scarr[l] = 0.0;
		for( k=0; k<nalphabets; k++ )
//			scarr[l] += n_dis[k][l] * cpmx1[k][i1];
//			scarr[l] += n_dis_consweight_multi[k][l] * cpmx1[k][i1];
			scarr[l] += scoreingmtx[k][l] * cpmx1[k][i1];
	}
	for( j=0; j<lgth2; j++ )
	{
		match[j] = 0.0;
		for( k=0; cpmxpdn[k][j]>-1; k++ )
			match[j] += scarr[cpmxpdn[k][j]] * cpmxpd[k][j];
	} 
	free( scarr );
#endif
}
static void match_calc( double **n_dynamicmtx, double *match, double **cpmx1, double **cpmx2, int i1, int lgth2, double **doublework, int **intwork, int initialize )
{
#if FASTMATCHCALC
//	fprintf( stderr, "\nmatch_calc... %d", i1 );
	int j, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	double *matchpt, *cpmxpdpt, **cpmxpdptpt;
	int *cpmxpdnpt, **cpmxpdnptpt;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
	if( initialize )
	{
		int count = 0;
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][j] )
				{
					cpmxpd[j][count] = cpmx2[l][j];
					cpmxpdn[j][count] = l;
					count++;
				}
			}
			cpmxpdn[j][count] = -1;
		}
	}

	{
		for( l=0; l<nalphabets; l++ )
		{
			scarr[l] = 0.0;
			for( j=0; j<nalphabets; j++ )
//				scarr[l] += n_dis[j][l] * cpmx1[j][i1];
//				scarr[l] += n_dis_consweight_multi[j][l] * cpmx1[j][i1];
				scarr[l] += n_dynamicmtx[j][l] * cpmx1[j][i1];
		}
		matchpt = match;
		cpmxpdnptpt = cpmxpdn;
		cpmxpdptpt = cpmxpd;
		while( lgth2-- )
		{
			*matchpt = 0.0;
			cpmxpdnpt = *cpmxpdnptpt++;
			cpmxpdpt = *cpmxpdptpt++;
			while( *cpmxpdnpt>-1 )
				*matchpt += scarr[*cpmxpdnpt++] * *cpmxpdpt++;
			matchpt++;
		} 
	}
	free( scarr );
//	fprintf( stderr, "done\n" );
#else
	int j, k, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
// simple
	if( initialize )
	{
		int count = 0;
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][j] )
				{
					cpmxpd[count][j] = cpmx2[l][j];
					cpmxpdn[count][j] = l;
					count++;
				}
			}
			cpmxpdn[count][j] = -1;
		}
	}
	for( l=0; l<nalphabets; l++ )
	{
		scarr[l] = 0.0;
		for( k=0; k<nalphabets; k++ )
//			scarr[l] += n_dis[k][l] * cpmx1[k][i1];
//			scarr[l] += n_dis_consweight_multi[k][l] * cpmx1[k][i1];
			scarr[l] += n_dynamicmtx[k][l] * cpmx1[k][i1];
	}
	for( j=0; j<lgth2; j++ )
	{
		match[j] = 0.0;
		for( k=0; cpmxpdn[k][j]>-1; k++ )
			match[j] += scarr[cpmxpdn[k][j]] * cpmxpd[k][j];
	} 
	free( scarr );
#endif
}

/* work/prog: dense variant of match_calc() for the per-row calls of the DP.
   match_calc() sums, for each column j, scarr[l]*cpmx2[l][j] over the letters l with a nonzero
   weight in that column, in ascending l, as a chain of fused multiply-adds starting from 0.0.
   Here the same chain runs over every letter that is nonzero in *any* column (ascending l).
   The extra terms have cpmx2[l][j] == 0.0, and fma(s, 0.0, acc) == acc exactly for finite s
   whenever acc is not -0.0; acc starts at +0.0 and can never become -0.0 (weights are >= 0 and
   exact cancellation rounds to +0.0), so every match[j] is bit-identical.  The loop over j is
   independent per column and vectorizes. */
#define DENSE_MAXU 32 // work/prog3: was 8; every letter of the alphabet (26) fits now
static TLS int dense_nu = -1;
static TLS int dense_u[DENSE_MAXU];

static void match_calc_dense_prep( double **cpmx2, int lgth2 )
{
	int l, j, n = 0;
	for( l=0; l<nalphabets; l++ )
	{
		double *c = cpmx2[l];
		for( j=0; j<lgth2; j++ ) if( c[j] ) break;
		if( j < lgth2 )
		{
			if( n == DENSE_MAXU ) { dense_nu = -1; return; }
			dense_u[n++] = l;
		}
	}
	dense_nu = n;
}

static void match_calc_dense( double **n_dynamicmtx, double * restrict match, double **cpmx1, double **cpmx2, int i1, int lgth2 )
{
	double scu[DENSE_MAXU];
	int u, j, k;
	for( u=0; u<dense_nu; u++ )
	{
		int l = dense_u[u];
		double sc = 0.0;
		for( k=0; k<nalphabets; k++ )
			sc += n_dynamicmtx[k][l] * cpmx1[k][i1];
		scu[u] = sc;
	}
	if( dense_nu == 4 ) // the usual case for nucleotides
	{
		const double s0 = scu[0], s1 = scu[1], s2 = scu[2], s3 = scu[3];
		const double * restrict c0 = cpmx2[dense_u[0]];
		const double * restrict c1 = cpmx2[dense_u[1]];
		const double * restrict c2 = cpmx2[dense_u[2]];
		const double * restrict c3 = cpmx2[dense_u[3]];
		for( j=0; j<lgth2; j++ )
		{
			double a = 0.0;
			a += s0 * c0[j];
			a += s1 * c1[j];
			a += s2 * c2[j];
			a += s3 * c3[j];
			match[j] = a;
		}
		return;
	}
	for( j=0; j<lgth2; j++ ) match[j] = 0.0;
	for( u=0; u<dense_nu; u++ )
	{
		const double sc = scu[u];
		const double * restrict c = cpmx2[dense_u[u]];
		for( j=0; j<lgth2; j++ )
			match[j] += sc * c[j];
	}
}

static void Atracking_localhom( double *impwmpt, double *lasthorizontalw, double *lastverticalw, 
						char **seq1, char **seq2, 
                        char **mseq1, char **mseq2, 
                        int **ijp, int icyc, int jcyc,
						int *warpis, int *warpjs, int warpbase,
						int *ngap1, int *ngap2, int reuseprofiles,
						char **gt1, char **gt2 )
{
	int i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k, limk;
	double wm;
	char *gaptable1, *gt1bk;
	char *gaptable2, *gt2bk;
	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );

	if( gt1 == NULL )
	{
		gt1bk = AllocateCharVec( lgth1+lgth2+1 );
		gt2bk = AllocateCharVec( lgth1+lgth2+1 );
	}
	else
	{
		gt1bk = *gt1;
		gt2bk = *gt2;
	}

#if 0
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif
 
	if( outgap == 1 )
		;
	else
	{
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
	}

    for( i=0; i<lgth1+1; i++ ) 
    {
        ijp[i][0] = i + 1;
    }
    for( j=0; j<lgth2+1; j++ ) 
    {
        ijp[0][j] = -( j + 1 );
    }

	gaptable1 = gt1bk + lgth1+lgth2;
	*gaptable1 = 0;
	gaptable2 = gt2bk + lgth1+lgth2;
	*gaptable2 = 0;
	//*ngap1 = *ngap2 = 0; // shita de keisan

	iin = lgth1; jin = lgth2;
	limk = lgth1+lgth2 + 1;
	*impwmpt = 0.0;
	for( k=0; k<limk; k++ ) 
	{
		if( ijp[iin][jin] >= warpbase )
		{
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
				*--gaptable1 = 'o';
				*--gaptable2 = '-';
				k++;
				//*ngap2 = 1; // shita de keisan
			}
			l= jin;
			while( --l >= 0 )
			{
				*--gaptable1 = '-';
				*--gaptable2 = 'o';
				k++;
				//*ngap1 = 1;
			}
			break;
		}
		else
		{
			l = iin - ifi;
			while( --l ) 
			{
				*--gaptable1 = 'o';
				*--gaptable2 = '-';
				k++;
				//*ngap2 = 1;
			}
			l= jin - jfi;
			while( --l )
			{
				*--gaptable1 = '-';
				*--gaptable2 = 'o';
				k++;
				//*ngap1 = 1;
			}
		}
		if( iin == lgth1 || jin == lgth2 )
			;
		else
		{
			*impwmpt += (double)imp_match_out_sc( iin, jin );

//		fprintf( stderr, "impwm = %f (iin=%d, jin=%d) seq1=%c, seq2=%c\n", *impwmpt, iin, jin, seq1[0][iin], seq2[0][jin] );
		}
		if( iin <= 0 || jin <= 0 ) break;
		*--gaptable1 = 'o';
		*--gaptable2 = 'o';
		k++;
		iin = ifi; jin = jfi;
	}


#if 1 // Atracking() to onaji.  Wazukani hayai hazu.  Test mada.
	if( strchr( gaptable1, '-' ) ) *ngap1 = 1; else *ngap1 = 0;
	if( strchr( gaptable2, '-' ) ) *ngap2 = 1; else *ngap2 = 0;

	if( *ngap1 == 0 && reuseprofiles )
		;
	else if( *ngap1 == 0 )
	{
		limk = gt1bk + lgth1+lgth2 - gaptable1;
		for( i=0; i<icyc; i++ ) 
		{
			strncpy( mseq1[i], seq1[i], limk );
			mseq1[i][limk] = 0;
		}
	}
	else
	{
		for( i=0; i<icyc; i++ ) gapireru( mseq1[i], seq1[i], gaptable1 );
	}

	if( *ngap2 == 0 && reuseprofiles )
		;
	else if( *ngap2 == 0 )
	{
		limk = gt2bk + lgth1+lgth2 - gaptable2; // nakutemo
		for( j=0; j<jcyc; j++ ) 
		{
			strncpy( mseq2[j], seq2[j], limk );
			mseq2[j][limk] = 0;
		}
	}
	else
	{
		for( j=0; j<jcyc; j++ ) gapireru( mseq2[j], seq2[j], gaptable2 );
	}
#else
	if( *ngap1 || !reuseprofiles )
		for( i=0; i<icyc; i++ ) gapireru( mseq1[i], seq1[i], gaptable1 );
//		gapmatometeireru( icyc, mseq1, seq1, gaptable1 );

	if( *ngap2 || !reuseprofiles )
		for( j=0; j<jcyc; j++ ) gapireru( mseq2[j], seq2[j], gaptable2 );
//		gapmatometeireru( jcyc, mseq2, seq2, gaptable2 );
#endif

	if( gt1 == NULL )
	{
		free( gt1bk );
		free( gt2bk );
	}
	else
	{
		*gt1 = gaptable1;
		*gt2 = gaptable2;
	}
}

typedef struct { double **cpmxresult; int limk, alen; double eff1, eff2; double ***cpmx1, ***cpmx2; char *gaptable1, *gaptable2; int usehist1, usehist2; } createcpmxctx_t;
static void createcpmx_letters( void *vctx, int lo, int hi ) // the per-letter body of createcpmxresult (unchanged)
{
	createcpmxctx_t *c = (createcpmxctx_t *)vctx;
	double **cpmxresult = c->cpmxresult, ***cpmx1 = c->cpmx1, ***cpmx2 = c->cpmx2;
	double eff1 = c->eff1, eff2 = c->eff2;
	char *gaptable1 = c->gaptable1, *gaptable2 = c->gaptable2;
	int limk = c->limk, alen = c->alen, usehist1 = c->usehist1, usehist2 = c->usehist2;
	int i, j, p;
	for( i=lo; i<hi; i++ )
	{
		cpmxresult[i] = calloc( limk+1, sizeof( double ) );
		for( j=0; j<alen; j++ ) cpmxresult[i][j] = 0.0;
		for( j=0,p=0; j<alen; j++ ) 
		{
			if( gaptable1[j] == '-' )
				;
				//if( amino_n['-'] == i ) cpmxresult[amino_n['-']][j] = eff1; // tsukawanai
			else
				cpmxresult[i][j] += (*cpmx1)[i][p++] * eff1;
		}

		for( j=0,p=0; j<alen; j++ ) 
		{
			if( gaptable2[j] == '-' )
				;
				//if( amino_n['-'] == i ) cpmxresult[amino_n['-']][j] = eff2; // tsukawanai
			else
				cpmxresult[i][j] += (*cpmx2)[i][p++] * eff2;
		}

		if( usehist1 ) 
		{
			free( (*cpmx1)[i] );
			(*cpmx1)[i] = NULL;
		}
		if( usehist2 ) 
		{
			free( (*cpmx2)[i] );
			(*cpmx2)[i] = NULL;
		}
	}
}

static void createcpmxresult( double **cpmxresult, int limk, double eff1, double eff2, double ***cpmx1, double ***cpmx2, char *gaptable1, char *gaptable2, int usehist1, int usehist2 ) // allocate, free ryouhou suru.
{
	int i, j, p;
	int alen = strlen( gaptable1 );

//	reporterr( "eff1 = %f, eff2=%f\n", eff1, eff2 );

#if 1 // sukoshi osoi
	/* work/prog3: the letters are independent (each writes and frees only its own rows); large
	   profiles are shared with idle helper threads (createcpmx_letters, same operations per entry) */
	createcpmxctx_t c;
	c.cpmxresult = cpmxresult; c.limk = limk; c.alen = alen; c.eff1 = eff1; c.eff2 = eff2;
	c.cpmx1 = cpmx1; c.cpmx2 = cpmx2; c.gaptable1 = gaptable1; c.gaptable2 = gaptable2; c.usehist1 = usehist1; c.usehist2 = usehist2;
#ifdef enablemultithread
	dphelp_for( nalphabets, 1, (double)nalphabets * alen, 2.0e5, createcpmx_letters, &c );
#else
	createcpmx_letters( &c, 0, nalphabets );
#endif
	(void)i; (void)j; (void)p;
#if 1
	if( usehist1 ) 
	{
		free( *cpmx1 );
		*cpmx1 = NULL;
	}
	if( usehist2 )
	{
		free( *cpmx2 );
		*cpmx2 = NULL;
	}
#endif
#endif
}

static void creategapfreqresult( double **gapfresult, int limk, double eff1, double eff2, double *gapf1, double *gapf2, char *gaptable1, char *gaptable2 ) // allocate nomi
{
	int j, p;
	int alen = strlen( gaptable1 );
	(*gapfresult) = calloc( limk+1, sizeof( double ) );

#if 1 // sukoshi osoi
	for( j=0; j<alen+1; j++ ) (*gapfresult)[j] = 0.0;
	for( j=0,p=0; j<alen+1; j++ ) 
	{
		if( gaptable1[j] == '-' )
			;
		else
			(*gapfresult)[j] += gapf1[p++] * eff1;
	}

	for( j=0,p=0; j<alen; j++ ) 
	{
		if( gaptable2[j] == '-' )
			;
		else
			(*gapfresult)[j] += gapf2[p++] * eff2;
	}
	(*gapfresult)[j] = 1.0;

//	free( gapf1 ); -> yobidashimoto
//	gapf1 = NULL;
//	free( gapf2 );
//	gapf1 = NULL;

#if 0
	reporterr( "gaptable1=%s\n", gaptable1 );
	reporterr( "gaptable2=%s\n", gaptable2 );
	reporterr( "result (gapfreq) = " );
	for( j=0; j<alen+1; j++ ) reporterr( "%4.2f ", (*gapfresult)[j] );
	reporterr( "\n" );
#endif
#endif
}

static void createogresult( double **gapfresult, int limk, double eff1, double eff2, double *ori1, double *ori2, double *gf1, double *gf2, char *gaptable1, char *gaptable2 ) // allocate nomi
{
	int j, p;
	int alen = strlen( gaptable1 );

	*gapfresult = calloc( limk+1, sizeof( double ) );
#if 1 // sukoshi osoi
	for( j=0; j<alen; j++ ) (*gapfresult)[j] = 0.0;
	for( j=0,p=0; j<alen; j++ ) 
	{
		if( gaptable1[j] == '-' )
		{
			if( j==0 )
			{
//				reporterr( "hit, j=%d, result += %f\n",  eff1 );
				(*gapfresult)[j] += 1.0 * eff1;
			}
			else if ( j && gaptable1[j-1] != '-' )
			{
//				reporterr( "hit, j=%d, p=%d, gf1[p]=%f, result += %f (1)\n", j, p, gf1[p], (gf1[p]) * eff1 );
				(*gapfresult)[j] += (gf1[p-1]) * eff1; // p-1 daijoubu?
			}
		}
		else
		{
			if( j==0 || ( j && gaptable1[j-1] != '-' ) )
			{
//				reporterr( "copying. j=%d, ori1[]=%f\n", j, ori1[p] );
				(*gapfresult)[j] += ori1[p] * eff1;
			}
			p++;
		}
	}

	for( j=0,p=0; j<alen; j++ ) 
	{
		if( gaptable2[j] == '-' )
		{
			if( j==0 )
			{
//				reporterr( "hit, j=%d, result += %f\n",  eff2 );
				(*gapfresult)[j] += 1.0 * eff2;
			}
			else if ( j && gaptable2[j-1] != '-' )
			{
//				reporterr( "hit, j=%d, p=%d, gf2[p]=%f, result += %f (2)\n", j, p, gf2[p], (gf2[p]) * eff2 );
				(*gapfresult)[j] += (gf2[p-1]) * eff2; // p-1 daijoubu?
			}
		}
		else
		{
			if( j==0 || ( j && gaptable2[j-1] != '-' ) )
			{
//				reporterr( "copying. j=%d, ori2[]=%f\n", j, ori2[p] );
				(*gapfresult)[j] += ori2[p] * eff2;
			}
			p++;
		}
	}
//	free( ori1 );
//	ori1 = NULL;
//	free( ori2 );
//	ori2 = NULL;
#if 0
	reporterr( "createogresult (o) = " );
	for( j=0; j<alen+1; j++ ) reporterr( "%4.2f ", (*gapfresult)[j] );
	reporterr( "\n" );
#endif
#endif
}

static void createfgresult( double **gapfresult, int limk, double eff1, double eff2, double *ori1, double *ori2, double *gf1, double *gf2, char *gaptable1, char *gaptable2 ) // allocate nomi
{
	int j, p;
	int alen = strlen( gaptable1 );
	(*gapfresult) = calloc( limk+1, sizeof( double ) );

	for( j=0; j<alen; j++ ) (*gapfresult)[j] = 0.0;
	for( j=0,p=0; j<alen; j++ ) 
	{
		if( gaptable1[j] == '-' )
		{
			if( j==alen-1 )
			{
				(*gapfresult)[j] += eff1;
			}
			else if( gaptable1[j+1] != '-' )
			{
//				reporterr( "hit, j=%d, p=%d, gf1[p]=%f, result += %f (f1)\n", j, p, gf1[p], (gf1[p]) * eff1 );
				(*gapfresult)[j] += gf1[p] * eff1;
			}
		}
		else
		{
			if( gaptable1[j+1] != '-' )
				(*gapfresult)[j] += ori1[p] * eff1;
			p++;
		}
	}

	for( j=0,p=0; j<alen; j++ ) 
	{
		if( gaptable2[j] == '-' )
		{
			if( j==alen-1 )
			{
				(*gapfresult)[j] += eff2;
			}
			else if( gaptable2[j+1] != '-' )
			{
//				reporterr( "hit, j=%d, p=%d, gf2[p]=%f, result += %f (f2)\n", j, p, gf1[p], (gf1[p]) * eff1 );
				(*gapfresult)[j] += gf2[p] * eff2;
			}
		}
		else
		{
			if( gaptable2[j+1] != '-' )
				(*gapfresult)[j] += ori2[p] * eff2;
			p++;
		}
	}

//	free( ori1 );
//	ori1 = NULL;
//	free( ori2 );
//	ori2 = NULL;
#if 0
	reporterr( "createogresult (f) = " );
	for( j=0; j<alen+1; j++ ) reporterr( "%4.2f ", (*gapfresult)[j] );
	reporterr( "\n" );
#endif
}

TLS a_lazy_t A_lz; // work/prog3: see lazygap.h
static TLS int atr_len1 = -1, atr_len2 = -1; // work/prog3: lengths for Atracking when rows are not current
/* work/v8s: band-major traceback of a tiled DP.  Band b (columns J[b] .. J[b+1]-1) keeps its codes in a
   region of its own, rows of w[b] ints one after the other (row i, column j at base[b][i*w[b] + j-J[b]]),
   instead of a slice of every row of the row-major matrix: a band then writes one contiguous stream,
   not 48-200 ints into a different page every row.  On CARC (32 threads over two sockets) the strided
   writes were ~2/3 of the tiled DPs' time (MAFFT_TILEIJPTEST, see notes/v8s-lever.md).  The same codes
   at the same (i, j); Atracking reads them through ATR_IJP.  The regions are carved from the same
   traceback block (rows 0 .. nrows+1 of every band: Atracking's corner may be row lgth1 = nrows+1). */
#define IJT_MAXB 32
typedef struct { int nb; int J[IJT_MAXB+1]; int *base[IJT_MAXB]; int w[IJT_MAXB]; } ijtile_t;
static TLS ijtile_t atr_tile;
static TLS int atr_tiled = 0; // set by Adp_tiled for the Atracking of the same A__align call
static inline int *ijt_cell( const ijtile_t *t, int i, int j )
{
	int lo = 0, hi = t->nb - 1;
	while( lo < hi ) { int mid = ( lo + hi + 1 ) / 2; if( t->J[mid] <= j ) lo = mid; else hi = mid - 1; }
	return( t->base[lo] + (size_t)i * t->w[lo] + ( j - t->J[lo] ) );
}

static double Atracking( double *lasthorizontalw, double *lastverticalw, 
						double fpenalty, double fpenalty_ex,
						char **seq1, char **seq2, 
                        char **mseq1, char **mseq2, 
                        int **ijp, int icyc, int jcyc,
						int tailgp,
						int *warpis, int *warpjs, int warpbase,
						int *ngap1, int *ngap2, int reuseprofiles,
#if 0
						double *eff1, double *eff2,
						double ***cpmxresult, 
						double **cpmx1, double **cpmx2, 
						double *gapfreq1, double *gapfreq2, 
						double *ogcp1o, double *ogcp2o, 
						double *fgcp1o, double *fgcp2o, 
						double orieff1, double orieff2,
						int usehist1, int usehist2 
#else
						char **gt1, char **gt2
#endif
						)
{
	int i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k, limk;
	double wm;
	char *gaptable1, *gt1bk;
	char *gaptable2, *gt2bk;
	lgth1 = ( atr_len1 >= 0 ) ? atr_len1 : strlen( seq1[0] ); // work/prog3: the rows may be stale (lazygap.h)
	lgth2 = ( atr_len2 >= 0 ) ? atr_len2 : strlen( seq2[0] );


	if( gt1 == NULL )
	{
		gt1bk = AllocateCharVec( lgth1+lgth2+1 );
		gt2bk = AllocateCharVec( lgth1+lgth2+1 );
//		gaptable1 = gt1bk + lgth1+lgth2;
//		gaptable2 = gt2bk + lgth1+lgth2;
	}
	else
	{
//		gaptable1 = *gt1 + lgth1+lgth2;
//		gaptable2 = *gt2 + lgth1+lgth2;
		gt1bk = *gt1;
		gt2bk = *gt2;
	}
//	*gaptable1 = 0;
//	*gaptable2 = 0;

#if 0
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif
 
	if( tailgp == 1 )
		;
	else
	{
#if 1
//		double fpenalty = (double)penalty;
//		double fpenalty_ex = (double)penalty_ex;
		double g;
//		reporterr( "in S, lastverticalw[lgth1-1] = %f\n", lastverticalw[lgth1-1] );
//		reporterr( "in S, lasthorizontalw[lgth2-1] = %f\n", lasthorizontalw[lgth2-1] );
		wm = lasthorizontalw[lgth2-1] - 1.0; // lasthorizontalw[lgth2-1] yori kanarazu chiisai.
		for( j=lgth2-2; j>=0; j-- )
		{
			if( (g=lasthorizontalw[j]+ ( fpenalty * TERMGAPFAC + fpenalty_ex * (lgth2-1-j) * TERMGAPFAC_EX ) ) > wm )
			{
				wm = g;
				iin = lgth1-1; jin = j;
				*( atr_tiled ? ijt_cell( &atr_tile, lgth1, lgth2 ) : &ijp[lgth1][lgth2] ) = -( lgth2 - j );
			}
		}
		for( i=lgth1-2; i>=0; i-- )
		{
			if( ( g=lastverticalw[i]+ ( fpenalty * TERMGAPFAC + fpenalty_ex * (lgth1-1-i) * TERMGAPFAC_EX ) ) > wm )
			{
				wm = g;
				iin = i; jin = lgth2-1;
				*( atr_tiled ? ijt_cell( &atr_tile, lgth1, lgth2 ) : &ijp[lgth1][lgth2] ) = +( lgth1 - i );
			}
		}
		if( lasthorizontalw[lgth2-1] > wm )  // score ga onaji baai erabarenai
		{
			wm = lasthorizontalw[lgth2-1];
			iin = lgth1-1; jin = lgth2-1;
			*( atr_tiled ? ijt_cell( &atr_tile, lgth1, lgth2 ) : &ijp[lgth1][lgth2] ) = 0;
		}
#else
		wm = lastverticalw[0];
		for( i=0; i<lgth1; i++ )
		{
			if( lastverticalw[i] >= wm )
			{
				wm = lastverticalw[i];
				iin = i; jin = lgth2-1;
				*( atr_tiled ? ijt_cell( &atr_tile, lgth1, lgth2 ) : &ijp[lgth1][lgth2] ) = +( lgth1 - i );
			}
		}
		for( j=0; j<lgth2; j++ )
		{
			if( lasthorizontalw[j] >= wm )
			{
				wm = lasthorizontalw[j];
				iin = lgth1-1; jin = j;
				*( atr_tiled ? ijt_cell( &atr_tile, lgth1, lgth2 ) : &ijp[lgth1][lgth2] ) = -( lgth2 - j );
			}
		}
#endif
	}

	/* work/prog3: row 0 and column 0 of ijp are no longer filled in (lgth1 + lgth2 scattered writes,
	   one per row of a matrix that can have 10^4 rows); the walk below uses the values they would
	   hold: ijp[0][j] = -( j + 1 ) (written last, so also ijp[0][0] = -1) and ijp[i][0] = i + 1. */
#define ATR_IJP( ii, jj ) ( (ii) == 0 ? -( (jj) + 1 ) : (jj) == 0 ? (ii) + 1 : atr_tiled ? *ijt_cell( &atr_tile, ii, jj ) : ijp[ii][jj] ) // work/v8s: atr_tiled

	gaptable1 = gt1bk + lgth1+lgth2;
	*gaptable1 = 0;
	gaptable2 = gt2bk + lgth1+lgth2;
	*gaptable2 = 0;

	iin = lgth1; jin = lgth2;
	limk = lgth1+lgth2 + 1;
	for( k=0; k<limk; k++ ) 
	{
		int ijv = ATR_IJP( iin, jin );
		if( ijv >= warpbase )
		{
			ifi = warpis[ijv-warpbase];
			jfi = warpjs[ijv-warpbase];
		}
		else if( ijv < 0 ) 
		{
			ifi = iin-1; jfi = jin+ijv;
		}
		else if( ijv > 0 )
		{
			ifi = iin-ijv; jfi = jin-1;
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
				*--gaptable1 = 'o';
				*--gaptable2 = '-';
				k++;
			}
			l= jin;
			while( --l >= 0 )
			{
				*--gaptable1 = '-';
				*--gaptable2 = 'o';
				k++;
			}
			break;
		}
		else
		{
			l = iin - ifi;
			while( --l ) 
			{
				*--gaptable1 = 'o';
				*--gaptable2 = '-';
				k++;
			}
			l= jin - jfi;
			while( --l )
			{
				*--gaptable1 = '-';
				*--gaptable2 = 'o';
				k++;
			}
		}
		if( iin <= 0 || jin <= 0 ) break;
		*--gaptable1 = 'o';
		*--gaptable2 = 'o';
		k++;
		iin = ifi; jin = jfi;
	}

#undef ATR_IJP
//	reporterr( "gaptable1=%s\n", gaptable1 );
//	reporterr( "gaptable2=%s\n", gaptable2 );
	if( strchr( gaptable1, '-' ) ) *ngap1 = 1; else *ngap1 = 0;
	if( strchr( gaptable2, '-' ) ) *ngap2 = 1; else *ngap2 = 0;
#if 1
	if( mseq1 == NULL ) // work/prog2: the caller inserts the gaps in place (A__align, gapireru_inplace)
		;
	else if( *ngap1 == 0 && reuseprofiles )
		;
	else if( *ngap1 == 0 )
	{
		limk = gt1bk + lgth1+lgth2 - gaptable1;
		for( i=0; i<icyc; i++ ) 
		{
			strncpy( mseq1[i], seq1[i], limk );
			mseq1[i][limk] = 0;
		}
	}
	else
	{
		for( i=0; i<icyc; i++ ) gapireru( mseq1[i], seq1[i], gaptable1 );
	}

	if( mseq2 == NULL ) // work/prog2: see above
		;
	else if( *ngap2 == 0 && reuseprofiles )
		;
	else if( *ngap2 == 0 )
	{
		limk = gt2bk + lgth1+lgth2 - gaptable2; // nakutemo
		for( j=0; j<jcyc; j++ ) 
		{
			strncpy( mseq2[j], seq2[j], limk );
			mseq2[j][limk] = 0;
		}
	}
	else
	{
		for( j=0; j<jcyc; j++ ) gapireru( mseq2[j], seq2[j], gaptable2 );
	}
#else
	if( *ngap1 || !reuseprofiles )
		for( i=0; i<icyc; i++ ) gapireru( mseq1[i], seq1[i], gaptable1 );

	if( *ngap2 || !reuseprofiles )
		for( j=0; j<jcyc; j++ ) gapireru( mseq2[j], seq2[j], gaptable2 );
#endif


	if( gt1 == NULL )
	{
		free( gt1bk );
		free( gt2bk );
	}
	else
	{
		*gt1 = gaptable1;
		*gt2 = gaptable2;
	}

	return( wm );
}

/* work/prog: fast evaluation of one DP row of A__align (trywarp == 0).  Same cells, same
   floating-point expressions (each a*b+c statement is contracted to one fmadd, exactly as in the
   original loop) and the same comparison order; only the evaluation is reorganized:
   (1) the column recurrence (m[], mp[]) and the candidate it offers to each cell are independent
       across j, so they are computed first in a separate, vectorizable pass;
   (2) the row recurrence (mi, mpi) is kept in registers instead of TLS statics;
   (3) all selects are branch-free and ijp is written once per cell. */
/* work/gccpar: C ? A : B for an int C in {0,1}, as bit arithmetic (see ADP_CELL below) */
#define ISEL( C, A, B ) ( (B) ^ ( ( (B) ^ (A) ) & -(C) ) )
/* work/gccpar: C >= X ? C : X for finite (non-NaN) doubles, bit-identical in both spellings (see
   ADP_CELL below).  x86: the "X > C ? X : C" spelling is an IEEE max, one maxsd/vmaxsd; arm64:
   the ">=" spelling lets clang share the compare with the matching integer select. */
#if defined( __x86_64__ ) || defined( __i386__ )
#define DSEL_GE( C, X ) ( (X) > (C) ? (X) : (C) )
#else
#define DSEL_GE( C, X ) ( (C) >= (X) ? (C) : (X) )
#endif

static void Adp_row_fast( int i, int lastj, double * restrict colw, int * restrict colij, double * restrict colg2,
	double * restrict mloc, int * restrict mploc, double * restrict curw, int * restrict ijpi,
	const double * restrict prew, const double * restrict gf2, const double * restrict fg2, const double * restrict og2,
	double fgcp1va, double ogcp1va, double gf1va, double gf1vapre, double fpenalty_ex,
	double *mip, int *mpip, double *wmp )
{
	int j, ijv;
#ifdef PROGTIMING
	double _t0 = pt_now(), _t1;
#endif
	double mil = *mip;
	int mpil = *mpip;
	double wmv = *wmp;

	for( j=1; j<lastj; j++ )
	{
		double mj = mloc[j];
		int mpj = mploc[j];
		double g3 = mj + fgcp1va * gf2[j];         // = *mjpt + fgcp1va * *gf2pt
		double g4 = prew[j-1] + ogcp1va * gf2[j-1]; // = *prept + ogcp1va * *gf2ptpre
		int take = ( g4 >= mj );
		colw[j] = g3;
		colij[j] = i - mpj;
		mloc[j] = DSEL_GE( g4, mj ) + fpenalty_ex;
		mploc[j] = ISEL( take, i-1, mpj );
	}
#ifdef PROGTIMING
	_t1 = pt_now(); pt_acc[21] += _t1 - _t0; _t0 = _t1;
#endif
	if( fpenalty_ex == 0.0 && lastj > 1 )
	{
		/* The row recurrence is  mi = ( g2 >= mi ? g2 : mi ) + fpenalty_ex.  Adding +0.0 only turns
		   -0.0 into +0.0, so after the first column mi is already "canonical" and mi + 0.0 == mi
		   bitwise.  Precomputing g2 + 0.0 (vectorized, below) and selecting between it and mi
		   therefore reproduces every mi bit-for-bit while taking the add off the serial chain.
		   Comparisons are unaffected by the sign of zero.  Column 1 is done exactly as before. */
		double * restrict g2c = colg2;
		for( j=1; j<lastj; j++ )
		{
			double g2 = prew[j-1] + og2[j] * gf1vapre; // = *prept + *ogcp2pt * gf1vapre
			g2c[j] = g2 + fpenalty_ex;
		}
#ifdef PROGTIMING
		_t1 = pt_now(); pt_acc[22] += _t1 - _t0; _t0 = _t1;
#endif
		j = 1;
		{
			double gg;
			wmv = prew[j-1];
			ijv = 0;
			gg = mil + fg2[j-1] * gf1va;
			if( gg > wmv ) { wmv = gg; ijv = -( j - mpil ); }
			gg = prew[j-1] + og2[j] * gf1vapre;
			if( gg >= mil ) { mil = gg; mpil = j-1; }
			mil += fpenalty_ex;
			{
				double cw = colw[j];
				int cij = colij[j];
				int take = ( cw > wmv );
				ijv = take ? cij : ijv;
				wmv = take ? cw : wmv;
			}
			curw[j] += wmv;
			ijpi[j] = ijv;
		}
		for( j=2; j<lastj; j++ )
		{
			double gg, g2 = g2c[j];
			int take2;
			double pw = prew[j-1];
			int th;
			gg = mil + fg2[j-1] * gf1va;                // = mi + *fgcp2pt * gf1va
			th = ( gg > pw );
			wmv = th ? gg : pw;
			ijv = -( j - mpil ) & -th;
			take2 = ( g2 >= mil );
			mil = DSEL_GE( g2, mil );
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
#ifdef PROGTIMING
		_t1 = pt_now(); pt_acc[23] += _t1 - _t0;
#endif
		return;
	}
	for( j=1; j<lastj; j++ )
	{
		double gg, pw = prew[j-1];
		int th;
		gg = mil + fg2[j-1] * gf1va;                // = mi + *fgcp2pt * gf1va
		th = ( gg > pw );
		wmv = th ? gg : pw;
		ijv = -( j - mpil ) & -th;
		gg = prew[j-1] + og2[j] * gf1vapre;         // = *prept + *ogcp2pt * gf1vapre
		th = ( gg >= mil );
		mil = DSEL_GE( gg, mil );
		mpil = ISEL( th, j-1, mpil );
		mil += fpenalty_ex;
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

/* work/prog: one DP row of A__align with the match scores fused in, for the common nucleotide
   case (exactly 4 letters in use in seq2's profile, see match_calc_dense) with a zero gap-extension
   penalty.  Per cell it performs the same operations as match_calc() + the original loop:
     match = (((0 + s0*c0) + s1*c1) + s2*c2) + s3*c3            (fmadd chain, see match_calc_dense)
     column candidate / m[] update, row candidate / mi update  (as in Adp_row_fast)
     currentw[j] = match + wm
   Column lgth2 has no match value in the original either (match_calc fills 0..lgth2-1 only), so
   its currentw entry keeps the value already in the buffer, exactly as before.  One scalar pass
   lets the out-of-order core overlap the independent per-column work with the serial mi chain. */
/* work/gccpar: the selects of the cell are written so that gcc compiles them without branches
   and x86 gets one instruction per double select.
   (1) Integer selects: gcc gimplifies "c ? a : b" into an if/else whenever an arm is not a plain
       value (i-1, -(j-mpil), i-mpj), sinks the arm into the conditional block, and its RTL
       if-conversion then left three of the four as data-dependent branches per cell (aarch64,
       gcc 15; v3's fused row took 1.9-2.1x the time per cell under gcc as under clang).  ISEL
       writes the same select as bit arithmetic, which gcc folds into a straight-line COND_EXPR
       that expands to csel/cmov on every target.  "if( gg > wm ) { wm = gg; ijv = ...; }" is
       now a double select plus an integer mask.
   (2) Double selects "c >= x ? c : x" go through DSEL_GE, which on x86 spells them
       "x > c ? x : c": for non-NaN operands both pick c exactly when c >= x, so the bits are the
       same (signed zeros included, since equal values select c either way); the second form is
       an IEEE max that x86 compiles to one maxsd (gcc: ix86_expand_sse_fp_minmax; clang
       likewise) instead of cmplesd/andpd/andnpd/orpd or blendvpd (clang x86-64 SSE2: 78 -> 59
       instructions per cell together with (1)).  DP values are finite, never NaN.
   Same comparisons for the traceback codes, same tie rules, same FP operations in the same
   order: every value is bit-identical to v3's (and to the original loop's). */
#define ADP_CELL( FIRST, HASMATCH ) ADP_CELL_PW( FIRST, HASMATCH, prew[j-1] )
#define ADP_CELL_PW( FIRST, HASMATCH, PWVAL ) \
	{ \
		double mj = mloc[j], pw = (PWVAL), base, g3, g4, gg, g2, wh; \
		int mpj = mploc[j], takec, take, th; \
		if( HASMATCH ) \
		{ \
			base = 0.0; \
			base += s0 * c0[j]; \
			base += s1 * c1[j]; \
			base += s2 * c2[j]; \
			base += s3 * c3[j]; \
		} \
		else \
			base = curw[j]; \
		g3 = mj + fgcp1va * gf2[j]; \
		g4 = pw + ogcp1va * gf2[j-1]; \
		takec = ( g4 >= mj ); \
		mloc[j] = DSEL_GE( g4, mj ) + fpenalty_ex; \
		mploc[j] = ISEL( takec, i-1, mpj ); \
		gg = mil + fg2[j-1] * gf1va; \
		th = ( gg > pw ); \
		wh = th ? gg : pw; \
		ijv = -( j - mpil ) & -th; \
		g2 = pw + og2[j] * gf1vapre; \
		if( FIRST ) \
		{ \
			if( g2 >= mil ) { mil = g2; mpil = j-1; } \
			mil += fpenalty_ex; \
		} \
		else \
		{ \
			double g2c = g2 + fpenalty_ex; \
			int take2 = ( g2c >= mil ); \
			mil = DSEL_GE( g2c, mil ); \
			mpil = ISEL( take2, j-1, mpil ); \
		} \
		take = ( g3 > wh ); \
		ijv = ISEL( take, i - mpj, ijv ); \
		wmv = take ? g3 : wh; \
		curw[j] = base + wmv; \
		ijpi[j] = ijv; \
	}

/* work/prog3: the inner cells of a row, j in [ja, jb), scalar (the reference for the vector kernels) */
static void Adp_cells_scalar( int i, int ja, int jb,
	double * restrict mloc, int * restrict mploc, double * restrict curw, int * restrict ijpi,
	const double * restrict prew, const double * restrict gf2, const double * restrict fg2, const double * restrict og2,
	const double * restrict c0, const double * restrict c1, const double * restrict c2, const double * restrict c3,
	const double *scu,
	double fgcp1va, double ogcp1va, double gf1va, double gf1vapre, double fpenalty_ex,
	double *mip, int *mpip, double *wmp )
{
	const double s0 = scu[0], s1 = scu[1], s2 = scu[2], s3 = scu[3];
	double mil = *mip, wmv = *wmp;
	int mpil = *mpip, ijv, j;
	for( j=ja; j<jb; j++ )
		ADP_CELL( 0, 1 )
	*mip = mil;
	*mpip = mpil;
	*wmp = wmv;
}

/* work/prog3: vector kernels for the inner cells (adpvec_body.h explains why they are exact).
   Adp_cells_v{2,4}: 4-letter inline match (as Adp_row_fused4); Adp_cellsm_v{2,4}: match in curw. */
#if defined( __aarch64__ ) || defined( __ARM_NEON ) || defined( __x86_64__ ) || defined( __i386__ )
#define ADPV_HAVE 1
#if ( defined( __x86_64__ ) || defined( __i386__ ) ) && !defined( __AVX2__ ) && ( defined( __GNUC__ ) || defined( __clang__ ) )
#define ADPV_X86DISPATCH 1
#define ADPV_ATTR4 __attribute__(( target( "avx2" ) )) // runtime dispatch below; FMA only if the build flags enable it
#else
#define ADPV_ATTR4
#endif
#define ADPV_ATTR
#define ADPV_VB 16
#define ADPV_MROW 0
#define ADPV_NAME Adp_cells_v2
#include "adpvec_body.h"
#undef ADPV_NAME
#undef ADPV_MROW
#define ADPV_MROW 1
#define ADPV_NAME Adp_cellsm_v2
#include "adpvec_body.h"
#undef ADPV_NAME
#undef ADPV_MROW
#undef ADPV_VB
#undef ADPV_ATTR
#define ADPV_ATTR ADPV_ATTR4
#define ADPV_VB 32
#define ADPV_MROW 0
#define ADPV_NAME Adp_cells_v4
#include "adpvec_body.h"
#undef ADPV_NAME
#undef ADPV_MROW
#define ADPV_MROW 1
#define ADPV_NAME Adp_cellsm_v4
#include "adpvec_body.h"
#undef ADPV_NAME
#undef ADPV_MROW
#undef ADPV_VB
#undef ADPV_ATTR
#endif

typedef void (*adp_cells_fn)( int, int, int, double *, int *, double *, int *,
	const double *, const double *, const double *, const double *,
	const double *, const double *, const double *, const double *, const double *,
	double, double, double, double, double, double *, int *, double * );

/* work/prog3: the inner cells with the match score already in curw[j] (ADP_CELL( 0, 0 )), scalar */
static void Adp_cellsm_scalar( int i, int ja, int jb,
	double * restrict mloc, int * restrict mploc, double * restrict curw, int * restrict ijpi,
	const double * restrict prew, const double * restrict gf2, const double * restrict fg2, const double * restrict og2,
	const double * restrict c0, const double * restrict c1, const double * restrict c2, const double * restrict c3,
	const double *scu,
	double fgcp1va, double ogcp1va, double gf1va, double gf1vapre, double fpenalty_ex,
	double *mip, int *mpip, double *wmp )
{
	const double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0; // (unused: ADP_CELL( 0, 0 ) reads curw)
	double mil = *mip, wmv = *wmp;
	int mpil = *mpip, ijv, j;
	for( j=ja; j<jb; j++ )
		ADP_CELL( 0, 0 )
	*mip = mil;
	*mpip = mpil;
	*wmp = wmv;
}

typedef void (*adp_mrow_fn)( int, int, unsigned int, const double *, const double * const *, const unsigned int *, double * );
static adp_mrow_fn adp_mrow = NULL; // set by adp_cells_init
static void Adp_mrow_blk_g( int ja, int jb, unsigned int rowmask, const double *su, const double * const *cu, const unsigned int *bm, double * restrict curw );
#ifdef ADPV_X86DISPATCH
static void Adp_mrow_blk_avx2( int ja, int jb, unsigned int rowmask, const double *su, const double * const *cu, const unsigned int *bm, double * restrict curw );
#endif

/* Which kernels: MAFFT_ADPVEC=0 scalar, 2 two lanes, 4 four lanes.  Default: four lanes on x86 with
   AVX2 (checked at run time unless the build targets AVX2), two lanes on other x86, scalar on ARM
   (on Apple M-series cores the scalar fused loop is as fast as two NEON lanes: measured, see
   notes/progressive-lever-3.md). */
static adp_cells_fn adp_cells4 = NULL, adp_cellsm = NULL;
static void adp_cells_init( void )
{
	adp_cells_fn f4 = Adp_cells_scalar, fm = Adp_cellsm_scalar;
#ifdef ADPV_HAVE
	char *e = getenv( "MAFFT_ADPVEC" );
	int w = e ? atoi( e ) : -1;
	if( w < 0 )
	{
#if defined( __x86_64__ ) || defined( __i386__ )
#ifdef ADPV_X86DISPATCH
		__builtin_cpu_init();
		w = __builtin_cpu_supports( "avx2" ) ? 4 : 2;
#else
		w = 4;
#endif
#else
		w = 0;
#endif
	}
#ifdef ADPV_X86DISPATCH
	if( w == 4 ) { __builtin_cpu_init(); if( !__builtin_cpu_supports( "avx2" ) ) w = 2; }
#endif
	if( w == 2 ) { f4 = Adp_cells_v2; fm = Adp_cellsm_v2; }
	else if( w == 4 ) { f4 = Adp_cells_v4; fm = Adp_cellsm_v4; }
#endif
	{
		adp_mrow_fn fr = Adp_mrow_blk_g;
#ifdef ADPV_X86DISPATCH
		__builtin_cpu_init();
		if( __builtin_cpu_supports( "avx2" ) ) fr = Adp_mrow_blk_avx2;
#endif
		__atomic_store_n( &adp_mrow, fr, __ATOMIC_RELAXED );
	}
	__atomic_store_n( &adp_cellsm, fm, __ATOMIC_RELAXED );
	__atomic_store_n( &adp_cells4, f4, __ATOMIC_RELEASE ); // every thread computes the same values
}
static inline adp_cells_fn adp_cells_get( int mrow )
{
	if( __atomic_load_n( &adp_cells4, __ATOMIC_ACQUIRE ) == NULL ) adp_cells_init();
	return( mrow ? adp_cellsm : adp_cells4 );
}

#ifdef PROGTIMING
/* MAFFT_VECCHECK=1: every vectorized cell range is also computed by the scalar kernel on shadow
   copies of the row state, and everything both produce is compared bitwise. */
static int pt_veccheck_v = -1;
static int pt_veccheck( void ) { if( pt_veccheck_v < 0 ) pt_veccheck_v = ( getenv( "MAFFT_VECCHECK" ) != NULL ); return( pt_veccheck_v ); }
static TLS double *vc_m = NULL, *vc_cur = NULL; static TLS int *vc_mp = NULL, *vc_ijp = NULL; static TLS int vc_cap = 0;
static long vc_ranges = 0, vc_cells = 0; // atomic counters, reported at exit
static void vc_report( void ) __attribute__(( destructor ));
static void vc_report( void ) { if( pt_veccheck_v > 0 ) reporterr( "\nVECCHECK: %ld cell ranges, %ld cells verified, no mismatch\n", vc_ranges, vc_cells ); }
#endif

#ifdef PROGTIMING
#define ADPV_CHECKING pt_veccheck()
#else
#define ADPV_CHECKING 0
#endif
/* the inner cells [ja, jb) of row i; mrow: the match scores are in curw (else 4 inline letters) */
static inline void Adp_cells( int mrow, int i, int ja, int jb,
	double * restrict mloc, int * restrict mploc, double * restrict curw, int * restrict ijpi,
	const double * restrict prew, const double * restrict gf2, const double * restrict fg2, const double * restrict og2,
	const double * restrict c0, const double * restrict c1, const double * restrict c2, const double * restrict c3,
	const double *scu,
	double fgcp1va, double ogcp1va, double gf1va, double gf1vapre, double fpenalty_ex,
	double *mip, int *mpip, double *wmp )
{
	adp_cells_fn f = adp_cells_get( mrow );
#ifdef PROGTIMING
	if( pt_veccheck() && jb > ja )
	{
		double mi2 = *mip, wm2 = *wmp; int mpi2 = *mpip;
		if( jb + 1 > vc_cap )
		{
			vc_cap = jb + 1 + jb / 2;
			vc_m = realloc( vc_m, sizeof( double ) * vc_cap ); vc_cur = realloc( vc_cur, sizeof( double ) * vc_cap );
			vc_mp = realloc( vc_mp, sizeof( int ) * vc_cap ); vc_ijp = realloc( vc_ijp, sizeof( int ) * vc_cap );
		}
		memcpy( vc_m + ja, mloc + ja, sizeof( double ) * ( jb - ja ) );
		memcpy( vc_mp + ja, mploc + ja, sizeof( int ) * ( jb - ja ) );
		memcpy( vc_cur + ja, curw + ja, sizeof( double ) * ( jb - ja ) ); // (the match scores, if mrow)
		( mrow ? Adp_cellsm_scalar : Adp_cells_scalar )( i, ja, jb, vc_m, vc_mp, vc_cur, vc_ijp, prew, gf2, fg2, og2, c0, c1, c2, c3, scu,
			fgcp1va, ogcp1va, gf1va, gf1vapre, fpenalty_ex, &mi2, &mpi2, &wm2 );
		f( i, ja, jb, mloc, mploc, curw, ijpi, prew, gf2, fg2, og2, c0, c1, c2, c3, scu,
			fgcp1va, ogcp1va, gf1va, gf1vapre, fpenalty_ex, mip, mpip, wmp );
		if( memcmp( vc_m + ja, mloc + ja, sizeof( double ) * ( jb - ja ) ) || memcmp( vc_mp + ja, mploc + ja, sizeof( int ) * ( jb - ja ) )
		 || memcmp( vc_cur + ja, curw + ja, sizeof( double ) * ( jb - ja ) ) || memcmp( vc_ijp + ja, ijpi + ja, sizeof( int ) * ( jb - ja ) )
		 || memcmp( &mi2, mip, sizeof( double ) ) || mpi2 != *mpip || memcmp( &wm2, wmp, sizeof( double ) ) )
		{
			reporterr( "\nVECCHECK MISMATCH in A__align row %d, cells %d..%d (mrow=%d)\n", i, ja, jb-1, mrow );
			exit( 1 );
		}
		__atomic_add_fetch( &vc_ranges, 1, __ATOMIC_RELAXED );
		__atomic_add_fetch( &vc_cells, (long)( jb - ja ), __ATOMIC_RELAXED );
		return;
	}
#endif
	f( i, ja, jb, mloc, mploc, curw, ijpi, prew, gf2, fg2, og2, c0, c1, c2, c3, scu,
		fgcp1va, ogcp1va, gf1va, gf1vapre, fpenalty_ex, mip, mpip, wmp );
}

static void Adp_row_fused4( int i, int lgth2,
	double * restrict mloc, int * restrict mploc, double * restrict curw, int * restrict ijpi,
	const double * restrict prew, const double * restrict gf2, const double * restrict fg2, const double * restrict og2,
	const double * restrict c0, const double * restrict c1, const double * restrict c2, const double * restrict c3,
	const double *scu,
	double fgcp1va, double ogcp1va, double gf1va, double gf1vapre, double fpenalty_ex,
	double *mip, int *mpip, double *wmp )
{
	// requires lgth2 >= 1 and fpenalty_ex == 0.0 (see Adp_row_fast for why mi + 0.0 can leave the chain)
	const double s0 = scu[0], s1 = scu[1], s2 = scu[2], s3 = scu[3];
	double mil = *mip, wmv = *wmp;
	int mpil = *mpip, ijv, j;

	j = 1;
	if( 1 < lgth2 ) ADP_CELL( 1, 1 ) else ADP_CELL( 1, 0 )
	if( lgth2 > 2 ) // work/prog3: the inner cells 2..lgth2-1, vectorized
	{
		if( adp_cells_get( 0 ) == Adp_cells_scalar && !ADPV_CHECKING ) // the loop of v3, inline (as compiled before)
		{
			for( j=2; j<lgth2; j++ )
				ADP_CELL( 0, 1 )
		}
		else
			Adp_cells( 0, i, 2, lgth2, mloc, mploc, curw, ijpi, prew, gf2, fg2, og2, c0, c1, c2, c3, scu,
				fgcp1va, ogcp1va, gf1va, gf1vapre, fpenalty_ex, &mil, &mpil, &wmv );
	}
	if( lgth2 >= 2 )
	{
		j = lgth2;
		ADP_CELL( 0, 0 )
	}
	*mip = mil;
	*mpip = mpil;
	*wmp = wmv;
}

/* work/prog3: MAFFT_NOGENERAL=1 (timing builds) turns the general-letter path off (A/B) */
static int adp_general_v = -1;
static int adp_general_on( void )
{
#ifdef PROGTIMING
	if( adp_general_v < 0 ) adp_general_v = ( getenv( "MAFFT_NOGENERAL" ) == NULL );
	return( adp_general_v );
#else
	return( 1 );
#endif
}

/* work/prog3: one DP row for any number of letters (the rows Adp_row_fused4 cannot take: more than 4
   letters with a nonzero row score, e.g. IUPAC codes in nucleotide profiles, or protein).  For each
   block of columns, the match scores are formed first, in curw (exactly the chain of
   match_calc_dense / Adp_row_fused4: from 0.0, one a += s*c statement per letter, in ascending
   letter order), then the cells of the block run with the match read back (ADP_CELL( ., 0 ), the
   operations of the original loop, see Adp_cells / adpvec_body.h).  Column lgth2 has no match score
   and keeps the value in the buffer, as in the original. */
#define ADP_MBLOCK 512
#define ADPM_NAME Adp_mrow_blk_g
#define ADPM_ATTR
#include "mrowblk_body.h"
#undef ADPM_NAME
#undef ADPM_ATTR
#ifdef ADPV_X86DISPATCH
#define ADPM_NAME Adp_mrow_blk_avx2
#define ADPM_ATTR __attribute__(( target( "avx2" ) ))
#include "mrowblk_body.h"
#undef ADPM_NAME
#undef ADPM_ATTR
#endif
static void Adp_mrow( int ja, int jb, unsigned int rowmask, const double *su, const double * const *cu, const unsigned int *bm, double * restrict curw )
{
	if( __atomic_load_n( &adp_cells4, __ATOMIC_ACQUIRE ) == NULL ) adp_cells_init();
	if( jb > ja ) adp_mrow( ja, jb, rowmask, su, cu, bm, curw );
}

/* the block masks of the profile (see mrowblk_body.h): bm[b] = the letters u with a nonzero weight in
   some column 4b..4b+3 (columns 0..lgth2-1) */
static TLS unsigned int *gen_bm = NULL; // work/prog3: block masks of the current A__align's profile 2
static TLS int gen_bmcap = 0;
static void Adp_blockmasks( unsigned int *bm, const double * const *cu, int nu, int lgth2 )
{
	int u, j, nb = lgth2 / 4 + 2;
	for( j=0; j<nb; j++ ) bm[j] = 0;
	for( u=0; u<nu; u++ )
	{
		const double *c = cu[u];
		for( j=0; j<lgth2; j++ ) if( c[j] ) bm[j>>2] |= 1u << u;
	}
}

#ifdef PROGTIMING
/* work/v8s: MAFFT_V8STAT=1 (timing builds) -- which row kernel the A__align DP rows take, how many
   letters they carry, and why a DP was not tiled; printed at exit.  Counters only, no effect on
   the computation. */
static int v8s_stat_v = -1;
static int v8s_stat_on( void ) { if( v8s_stat_v < 0 ) v8s_stat_v = ( getenv( "MAFFT_V8STAT" ) != NULL ); return( v8s_stat_v ); }
static double v8s_cells[2][3][DENSE_MAXU+1]; // [serial/tiled][rowfused 0/1/2][letters]: cells
static double v8s_eff4[2], v8s_eff16[2];       // gen rows: letter*cells of the 4- and 16-column block chains
static double v8s_notile[8], v8s_notilen[8];   // cells / DPs not tiled, by reason
static double v8s_tiled[2];                    // tiled DPs: count, cells
static double v8s_gentime[2], v8s_f4time;      // serial gen rows: time, cells; serial fused4 rows: time
static double v8s_tile[16];                    // tiled DPs: elapsed, participation, row time, claims, threads, band changes, bands, elapsed*bands; [8..13] the same for DPs >= 2e7 cells
static void v8s_add( double *p, double v ) // atomic add on a double (CAS loop)
{
	long long *q = (long long *)p, oi, ni;
	double o, n;
	oi = __atomic_load_n( q, __ATOMIC_RELAXED );
	do { memcpy( &o, &oi, sizeof o ); n = o + v; memcpy( &ni, &n, sizeof ni ); }
	while( !__atomic_compare_exchange_n( q, &oi, ni, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED ) );
}
static void v8s_stat_row( int tiled, int kind, int nl, int lgth2, unsigned int rowmask, const unsigned int *bm )
{
	if( nl > DENSE_MAXU ) nl = DENSE_MAXU;
	v8s_add( &v8s_cells[tiled][kind][nl], (double)lgth2 );
	if( kind == 2 && bm )
	{
		int b, nb = ( lgth2 + 3 ) / 4;
		double e4 = 0.0, e16 = 0.0;
		for( b=0; b<nb; b++ ) e4 += 4.0 * __builtin_popcount( rowmask & bm[b] );
		for( b=0; b<nb; b+=4 ) e16 += 16.0 * __builtin_popcount( rowmask & ( bm[b] | ( b+1 < nb ? bm[b+1] : 0 ) | ( b+2 < nb ? bm[b+2] : 0 ) | ( b+3 < nb ? bm[b+3] : 0 ) ) );
		v8s_add( &v8s_eff4[tiled], e4 );
		v8s_add( &v8s_eff16[tiled], e16 );
	}
}
static void v8s_stat_report( void ) __attribute__(( destructor ));
static void v8s_stat_report( void )
{
	int s, k, n;
	static const char *why[8] = { "nrows<2|cells<min|dense_nu<0", "no idle helper", "lgth2<2*minw", "nb<2", "NOGENERAL", "alloc", "", "" };
	if( v8s_stat_v <= 0 ) return;
	for( s=0; s<2; s++ )
	{
		double tot = 0.0;
		for( k=0; k<3; k++ ) for( n=0; n<=DENSE_MAXU; n++ ) tot += v8s_cells[s][k][n];
		reporterr( "\nV8STAT %s rows: %.3g cells\n", s ? "tiled" : "serial", tot );
		for( k=0; k<3; k++ )
		{
			double t = 0.0;
			for( n=0; n<=DENSE_MAXU; n++ ) t += v8s_cells[s][k][n];
			if( t == 0.0 ) continue;
			reporterr( "V8STAT   kernel %s: %.3g cells (%.1f%%); by letters:", k == 0 ? "two-pass" : k == 1 ? "fused4" : "general", t, 100.0 * t / tot );
			for( n=0; n<=DENSE_MAXU; n++ ) if( v8s_cells[s][k][n] > 0.0 ) reporterr( " %d:%.1f%%", n, 100.0 * v8s_cells[s][k][n] / t );
			reporterr( "\n" );
			if( k == 2 ) reporterr( "V8STAT   general rows: letters per cell, 4-col chains %.2f, 16-col chains %.2f\n", v8s_eff4[s] / t, v8s_eff16[s] / t );
		}
	}
	reporterr( "V8STAT tiled DPs: %.0f, %.3g cells\n", v8s_tiled[0], v8s_tiled[1] );
	for( k=0; k<8; k++ ) if( v8s_notilen[k] > 0.0 ) reporterr( "V8STAT not tiled (%s): %.0f DPs, %.3g cells\n", why[k], v8s_notilen[k], v8s_notile[k] );
	if( v8s_tile[0] > 0.0 )
	{
		reporterr( "V8STAT tiled DPs: elapsed %.3f s, thread participation %.3f s, row work %.3f s (%.0f%% of participation), %.2f threads and %.1f bands per DP, %.0f claims (%.0f band changes)\n",
			v8s_tile[0], v8s_tile[1], v8s_tile[2], 100.0 * v8s_tile[2] / v8s_tile[1], v8s_tile[4] / v8s_tiled[0], v8s_tile[6] / v8s_tiled[0], v8s_tile[3], v8s_tile[5] );
		if( v8s_tile[14] > 0.0 ) reporterr( "V8STAT critical tilings: %.0f DPs, %.2f participants per DP\n", v8s_tile[14], v8s_tile[15] / v8s_tile[14] );
		if( v8s_tile[11] > 0.0 ) reporterr( "V8STAT tiled DPs >= 2e7 cells: %.0f DPs, elapsed %.3f s, participation %.3f s, row work %.3f s (%.0f%%), %.2f threads, %.1f bands per DP\n",
			v8s_tile[11], v8s_tile[8], v8s_tile[10], v8s_tile[9], 100.0 * v8s_tile[9] / v8s_tile[10], v8s_tile[12] / v8s_tile[11], v8s_tile[13] / v8s_tile[11] );
	}
	reporterr( "V8STAT serial general rows: %.3f s for %.3g cells (%.2f ns/cell); serial fused4 rows %.3f s\n", v8s_gentime[0], v8s_gentime[1], v8s_gentime[1] > 0.0 ? 1e9 * v8s_gentime[0] / v8s_gentime[1] : 0.0, v8s_f4time );
}
#define V8S_NOTILE( r, c ) do { if( v8s_stat_on() ) { v8s_add( &v8s_notile[r], (c) ); v8s_add( &v8s_notilen[r], 1.0 ); } } while( 0 )
#else
#define V8S_NOTILE( r, c )
#endif

static void Adp_row_denseK( int i, int lgth2, unsigned int rowmask, const double *su, const double * const *cu, const unsigned int *bm,
	double * restrict mloc, int * restrict mploc, double * restrict curw, int * restrict ijpi,
	const double * restrict prew, const double * restrict gf2, const double * restrict fg2, const double * restrict og2,
	double fgcp1va, double ogcp1va, double gf1va, double gf1vapre, double fpenalty_ex,
	double *mip, int *mpip, double *wmp )
{
	// requires lgth2 >= 1 and fpenalty_ex == 0.0 (as Adp_row_fused4)
	const double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0; // (unused: ADP_CELL( ., 0 ) reads curw)
	const double *c0 = NULL, *c1 = NULL, *c2 = NULL, *c3 = NULL;
	double mil = *mip, wmv = *wmp;
	int mpil = *mpip, ijv, j, b0, b1;

	for( b0=1; b0<=lgth2; b0=b1 )
	{
		b1 = MIN( b0 + ADP_MBLOCK, lgth2 + 1 );
		Adp_mrow( b0, MIN( b1, lgth2 ), rowmask, su, cu, bm, curw ); // columns < lgth2 only
		j = b0;
		if( j == 1 )
		{
			ADP_CELL( 1, 0 )
			j = 2;
		}
		if( MIN( b1, lgth2 ) > j )
		{
			Adp_cells( 1, i, j, MIN( b1, lgth2 ), mloc, mploc, curw, ijpi, prew, gf2, fg2, og2, c0, c1, c2, c3, NULL,
				fgcp1va, ogcp1va, gf1va, gf1vapre, fpenalty_ex, &mil, &mpil, &wmv );
			j = MIN( b1, lgth2 );
		}
		if( j == lgth2 && j < b1 && lgth2 >= 2 )
			ADP_CELL( 0, 0 )
	}
	*mip = mil;
	*mpip = mpil;
	*wmp = wmv;
}

static void dense_scores( double **n_dynamicmtx, double **cpmx1, int i1, double *scu )
{
	int u, k;
	for( u=0; u<dense_nu; u++ ) // same chain as match_calc's scarr[l]
	{
		int l = dense_u[u];
		double sc = 0.0;
		for( k=0; k<nalphabets; k++ )
			sc += n_dynamicmtx[k][l] * cpmx1[k][i1];
		scu[u] = sc;
	}
}

/* work/prog2: the same sums over only the letters k that have a nonzero weight in some row
   1..lasti-1 of cpmx1 (dense1_k, ascending; dense1_prep).  Every skipped term is
   n_dynamicmtx[k][l] * 0.0 added to the running sum, and fma( x, 0.0, sc ) == sc exactly because
   sc starts at +0.0 and never becomes -0.0 (see match_calc_dense), so each sum is bit-identical
   while the loop drops from nalphabets to a handful of letters. */
static TLS int dense1_n = -1;
static TLS int dense1_k[0x100];

static void dense1_prep( double **cpmx1, int ncol ) // letters with a nonzero weight in columns 0..ncol-1
{
	int k, i;
	dense1_n = 0;
	for( k=0; k<nalphabets; k++ )
	{
		const double *c = cpmx1[k];
		for( i=0; i<ncol; i++ ) if( c[i] != 0.0 ) break;
		if( i < ncol ) dense1_k[dense1_n++] = k;
	}
}

/* work/prog2: initverticalw without the sparse tables.  The original call
     match_calc( n_dynamicmtx, initverticalw, cpmx2, cpmx1, 0, lgth1, doublework, intwork, 1 )
   forms scarr[l] = sum over k < nalphabets of n_dynamicmtx[k][l] * cpmx2[k][0] and then, for each
   column i of cpmx1, the chain 0.0 + scarr[l] * cpmx1[l][i] over the letters l with a nonzero
   weight in that column, ascending.  Running the chain over every letter of dense1_k (a superset,
   ascending) adds only exact zeros (see match_calc_dense). */
static void dense_initvertical( double **n_dynamicmtx, double *match, double **cpmx2, double **cpmx1, int lgth1 )
{
	double sc[0x100];
	int t, k, i;
	for( t=0; t<dense1_n; t++ )
	{
		int l = dense1_k[t];
		double a = 0.0;
		for( k=0; k<nalphabets; k++ )
			a += n_dynamicmtx[k][l] * cpmx2[k][0];
		sc[t] = a;
	}
	for( i=0; i<lgth1; i++ )
	{
		double a = 0.0;
		for( t=0; t<dense1_n; t++ )
			a += sc[t] * cpmx1[dense1_k[t]][i];
		match[i] = a;
	}
}

/* the initialize part of match_calc(): sparse lists of the nonzero letters of every column */
static void match_calc_sparse_init( double **cpmx2, int lgth2, double **doublework, int **intwork )
{
	int j, l, count;
	for( j=0; j<lgth2; j++ )
	{
		count = 0;
		for( l=0; l<nalphabets; l++ )
		{
			if( cpmx2[l][j] )
			{
				doublework[j][count] = cpmx2[l][j];
				intwork[j][count] = l;
				count++;
			}
		}
		intwork[j][count] = -1;
	}
}

static void dense_scores1( double **n_dynamicmtx, double **cpmx1, int i1, double *scu )
{
	int u, t;
	for( u=0; u<dense_nu; u++ )
	{
		int l = dense_u[u];
		double sc = 0.0;
		for( t=0; t<dense1_n; t++ )
		{
			int k = dense1_k[t];
			sc += n_dynamicmtx[k][l] * cpmx1[k][i1];
		}
		scu[u] = sc;
	}
}

#ifdef enablemultithread
/* ---------------------------------------------------------------------------------------------
   work/prog2: helper pool (see dphelp.h) */
dphelppool_t *dphelppool = NULL;

int dphelp_idle( void )
{
	dphelppool_t *pool = dphelppool;
	return( pool ? pool->nidle : 0 );
}

TLS double dphelp_myprio = 0.0; // work/v8s (dphelp.h)
static int dphelp_crit_v = -1;
int dphelp_crit_on( void )
{
	if( dphelp_crit_v < 0 )
	{
#ifdef PROGTIMING
		char *e = getenv( "MAFFT_CRITHELP" ); // off by default (no gain measured, see notes/v8s-lever.md); =1 turns it on
		dphelp_crit_v = ( e && atoi( e ) != 0 );
#else
		dphelp_crit_v = 0;
#endif
	}
	return( dphelp_crit_v );
}
int dphelp_iscritical( void )
{
	dphelppool_t *pool = dphelppool;
	return( pool && dphelp_crit_on() && dphelp_myprio > pool->readyprio );
}
void dphelp_post_crit( dphelpjob_t *job, int maxslots )
{
	dphelppool_t *pool = dphelppool;
	job->slots = 0;
	job->refs = 0;
	job->next = NULL;
	job->posted = 0;
	job->prio = dphelp_myprio;
	if( pool == NULL || maxslots <= 0 ) return;
	pthread_mutex_lock( pool->mutex );
	job->slots = maxslots;
	job->next = pool->jobs;
	pool->jobs = job;
	job->posted = 1;
	if( pool->nidle > 0 ) pthread_cond_broadcast( pool->cond );
	pthread_mutex_unlock( pool->mutex );
}
int dphelp_serve_crit_locked( double readyprio )
{
	dphelppool_t *pool = dphelppool;
	dphelpjob_t *job, **pp;
	if( pool == NULL ) return( 0 );
	for( pp=&pool->jobs; ( job = *pp ) != NULL; pp=&job->next )
		if( job->prio > readyprio && job->slots > 0 ) break;
	if( job == NULL ) return( 0 );
	if( --job->slots <= 0 ) *pp = job->next;
	__atomic_add_fetch( &job->refs, 1, __ATOMIC_ACQ_REL );
	pthread_mutex_unlock( pool->mutex );
	job->run( job->arg );
	__atomic_sub_fetch( &job->refs, 1, __ATOMIC_RELEASE ); // the poster may free the job after this
	pthread_mutex_lock( pool->mutex );
	return( 1 );
}

void dphelp_post( dphelpjob_t *job, int maxslots )
{
	dphelppool_t *pool = dphelppool;
	job->slots = 0;
	job->refs = 0;
	job->next = NULL;
	job->posted = 0;
	job->prio = dphelp_myprio;
	if( pool == NULL || maxslots <= 0 || pool->nidle <= 0 ) return;
	pthread_mutex_lock( pool->mutex );
	job->slots = MIN( maxslots, pool->nidle );
	if( job->slots > 0 )
	{
		job->next = pool->jobs;
		pool->jobs = job;
		job->posted = 1;
		pthread_cond_broadcast( pool->cond );
	}
	pthread_mutex_unlock( pool->mutex );
}

void dphelp_finish( dphelpjob_t *job )
{
	dphelppool_t *pool = dphelppool;
	int spins = 0;
	if( pool == NULL ) return;
	/* work/v8s: RACE FIX.  This used to take the mutex only when it saw job->slots > 0 or refs > 0,
	   reading slots without the mutex.  A helper in dphelp_serve_locked decrements slots (to 0, taking
	   the job off the list) and only then increments refs, both under the mutex; a poster reading in
	   between saw slots == 0 and refs == 0, returned, and its stack frame (the job, and dphelp_for's
	   loop state) died while the helper went on to call job->run from it.  Observed as a SIGBUS (a call
	   into a heap address from dphelp_serve_locked) once in ~100 local rnasim1k runs after work/v8s
	   added more dphelp_for calls; latent in every earlier version with helpers.  Now a posted job is
	   always withdrawn under the mutex: whatever a helper did with it happened before its unlock, so
	   refs is seen here. */
	if( job->posted )
	{
		dphelpjob_t **pp;
		pthread_mutex_lock( pool->mutex );
		for( pp=&pool->jobs; *pp; pp=&(*pp)->next ) if( *pp == job ) { *pp = job->next; break; }
		job->slots = 0;
		pthread_mutex_unlock( pool->mutex );
	}
	while( __atomic_load_n( &job->refs, __ATOMIC_ACQUIRE ) > 0 )
	{
		if( ++spins <= 4096 ) sched_yield();
		else { struct timespec ts; ts.tv_sec = 0; ts.tv_nsec = 100000; nanosleep( &ts, NULL ); } // work/v8s: never spin for long
	}
}

int dphelp_serve_locked( void )
{
	dphelppool_t *pool = dphelppool;
	dphelpjob_t *job;
	if( pool == NULL || ( job = pool->jobs ) == NULL ) return( 0 );
	if( --job->slots <= 0 ) pool->jobs = job->next;
	__atomic_add_fetch( &job->refs, 1, __ATOMIC_ACQ_REL );
	pthread_mutex_unlock( pool->mutex );
	job->run( job->arg );
	__atomic_sub_fetch( &job->refs, 1, __ATOMIC_RELEASE ); // the poster may free the job after this
	pthread_mutex_lock( pool->mutex );
	return( 1 );
}

typedef struct
{
	void (*body)( void *ctx, int lo, int hi );
	void *ctx;
	int n, chunk;
	int next; // atomic
} dphelpfor_t;

static void dphelp_for_run( void *arg )
{
	dphelpfor_t *f = (dphelpfor_t *)arg;
	int lo;
	while( ( lo = __atomic_fetch_add( &f->next, f->chunk, __ATOMIC_RELAXED ) ) < f->n )
		f->body( f->ctx, lo, MIN( lo + f->chunk, f->n ) );
}

static int dphelp_for_min = -1;
void dphelp_for( int n, int chunk, double work, double minwork, void (*body)( void *ctx, int lo, int hi ), void *ctx )
{
	dphelpfor_t f;
	dphelpjob_t job;
	int nchunk;
	if( dphelp_for_min < 0 )
	{
#ifdef PROGTIMING
		dphelp_for_min = getenv( "MAFFT_NOHELPFOR" ) ? 1 : 0;
#else
		dphelp_for_min = 0;
#endif
	}
	if( chunk < 1 ) chunk = 1;
	nchunk = ( n + chunk - 1 ) / chunk;
	if( dphelp_for_min || nchunk < 2 || work < minwork || dphelp_idle() <= 0 )
	{
		if( n > 0 ) body( ctx, 0, n );
		return;
	}
	f.body = body; f.ctx = ctx; f.n = n; f.chunk = chunk; f.next = 0;
	job.run = dphelp_for_run;
	job.arg = &f;
	dphelp_post( &job, nchunk - 1 );
	dphelp_for_run( &f );
	dphelp_finish( &job );
}

/* ---------------------------------------------------------------------------------------------
   work/prog2: the fused DP of A__align, split into column bands and row blocks ("tiles") that
   the calling thread and idle helper threads process in any order the dependencies allow.

   Row i of band b (columns J[b] .. J[b+1]-1) needs
     - row i-1 of the same band: m[], mp[], and the previous row's values prew[j-1] for
       j > J[b] (all in the band's own columns of the shared arrays), and
     - from band b-1 in row i: the row carry (mi, mpi) after column J[b]-1, and
       prew[J[b]-1] = the value of row i-1 in column J[b]-1.
   Band b-1 publishes those three values per row in bmi/bmpi/bw before it advances its row
   counter (release); band b reads them after seeing the counter (acquire).  Columns never
   cross bands, so each cell is computed by exactly the same ADP_CELL code, with the same
   operands, as in Adp_row_fused4; only the thread and the moment differ.  Rows live in the
   same two alternating buffers as in the serial loop (row i in buffer i&1), so even the unused
   stale value in column lgth2 (which the kernels never give a match score) evolves exactly as
   before.  A band is owned by at most one thread at a time (busy flag, taken with an acquiring
   compare-and-swap and released with a release store), so its column state moves between
   threads safely. */
#define ADPB_MAXB 32
#ifdef PROGTIMING
static int adpb_ijptest = 0; // work/v8s: diagnosis switch MAFFT_TILEIJPTEST, see adpb_row
#endif
#if defined( __x86_64__ ) || defined( __i386__ )
#define ADPB_RELAX() __builtin_ia32_pause()
#elif defined( __aarch64__ )
#define ADPB_RELAX() __asm__ __volatile__( "yield" )
#else
#define ADPB_RELAX() do { } while( 0 )
#endif
/* work/v8s: every wait of the tiled DP (a band for its predecessor, band 0 for its rows' selections, a
   thread for a band) spins briefly, then yields, then sleeps, so that waiting threads never keep a
   descheduled thread that they wait for off the cores.  Pure spinning (with a yield every few
   thousand polls) once made a 16S.T run take 658 s instead of ~10 s on this Mac with other load. */
static void adpb_wait( int *spins )
{
	int s = ++*spins;
	if( s <= 256 ) { int k; for( k=0; k<4; k++ ) ADPB_RELAX(); }
	else if( s <= 1024 ) sched_yield();
	else
	{
		struct timespec ts;
		ts.tv_sec = 0;
		ts.tv_nsec = ( s <= 1200 ) ? 50000 : 200000;
		nanosleep( &ts, NULL );
	}
}
typedef struct { int nu, n1; int u[DENSE_MAXU]; int k1[0x100]; } adpb_letters_t; // the poster's dense_u / dense1_k (TLS)
typedef struct
{
	int done;   // rows 1..done of this band are finished (atomic)
	int busy;   // a thread is working on this band (atomic)
	double wmlast; // work/v8s: the last band's wm after each row (see adpb_t)
	char pad[48];
} __attribute__(( aligned( 64 ) )) adpb_band_t; // work/v8s: one cache line per band (the array was not line-aligned)

typedef struct
{
	int nrows, lgth2, nb, R;
	int J[ADPB_MAXB+1];
	double *wb[2];                  // row i lives in wb[i&1]
	double *m; int *mp; int **ijp;
	const double *gf1, *gf2, *fgcp1, *ogcp1, *fgcp2, *ogcp2;
	const double *initverticalw; double *lastverticalw;
	const double *rs;               // per row: the row scores of the (first 4) selected letters (stride 4)
	const double **rc;              // per row: their profile columns (stride 4)
	const int *rn;                  // work/prog3: per row: number of letters; <= 4: fused kernel (padded to 4)
	const double *rsu;              // work/prog3: per row, the scores of every letter u of the profile (stride ks)
	const unsigned int *rmask;      // work/prog3: per row, the letters with a nonzero score
	const double * const *cu;       // work/prog3: the profile column of each letter u
	const unsigned int *bm;         // work/prog3: block masks of the profile (mrowblk_body.h)
	int ks;
	int lazyrows;                   // work/prog3: band 0 selects each row's letters (adpb_rowsetup)
	double **dynmtx, **cpmx1, **cpmx2;
	adpb_letters_t letters;
	double fpenalty_ex;
	double *bw, *bmi; int *bmpi;    // [b*(nrows+1)+i]: values handed from band b-1 to band b
	/* work/v8s: wmlast used to be here, written after every row by the last band, in the cache line
	   holding bw, bmi, bmpi and fpenalty_ex, which every band reads for every row: on CARC that line
	   moved between cores (often sockets) once per row and band.  It now lives in the last band's
	   own line (adpb_band_t), and setup_next (written during the selections) in a line of its own. */
	int setup_next __attribute__(( aligned( 64 ) )); // work/v8s: next row offset whose letter selection is unclaimed (atomic)
	char pad_setup[60];
	int sch;                        // work/v8s: rows per selection chunk
	int *ready;                     // work/v8s: per selection chunk, 1 when its rows are set up (atomic)
	int sticky;                     // work/v8s: band ownership policy (see adpb_run)
	int batch;                      // work/v8s: a band starts when its predecessor is this many rows ahead (or done)
	int localsetup;                 // work/v8s: every band makes its rows' letter selections itself (adpb_row)
	int claimcap;                   // work/v8s: rows a thread processes on a band before it picks again (adpb_run)
	int bandmajor;                  // work/v8s: traceback codes in band-major regions (ijt_cell)
	int *ijb[ADPB_MAXB]; int ijw[ADPB_MAXB];
	adpb_band_t band[ADPB_MAXB];
#ifdef PROGTIMING
	double st_busy, st_part, st_claims, st_nthr, st_migr; // work/v8s (MAFFT_V8STAT): row time, participation time, claims, threads, band changes
#endif
} adpb_t;

/* work/prog3: the per-row letter selection of the tiled DP (see Adp_tiled); returns the count */
/* work/v8s: the selection of row i into rs4[0..3], rc4[0..3], *rnp, scu[0..nu-1], *maskp (adpb_rowsetup's
   body, with the outputs of one row passed directly) */
static int adpb_rowsetup_to( double **n_dynamicmtx, double **cpmx1, double **cpmx2, int i,
	double *rs4, const double **rc4, int *rnp, double *scu, unsigned int *maskp, const adpb_letters_t *L )
{
	unsigned int mask = 0;
	int u, t, nsel = 0, n;
	for( u=0; u<L->nu; u++ ) // = dense_scores1(), with the letter lists passed explicitly (this may run on a helper thread)
	{
		int l = L->u[u];
		double sc = 0.0;
		for( t=0; t<L->n1; t++ )
		{
			int k = L->k1[t];
			sc += n_dynamicmtx[k][l] * cpmx1[k][i];
		}
		scu[u] = sc;
	}
	for( u=0; u<L->nu; u++ )
	{
		if( scu[u] == 0.0 ) continue;
		if( nsel < 4 )
		{
			rs4[nsel] = scu[u];
			rc4[nsel] = cpmx2[L->u[u]];
		}
		mask |= 1u << u;
		nsel++;
	}
	*rnp = n = nsel;
	*maskp = mask;
	for( ; nsel<4; nsel++ )
	{
		rs4[nsel] = 0.0;
		rc4[nsel] = cpmx2[0];
	}
	return( n );
}

static int adpb_rowsetup( double **n_dynamicmtx, double **cpmx1, double **cpmx2, int i, int ks,
	double *rs, const double **rc, int *rn, double *rsu, unsigned int *rmask, const adpb_letters_t *L )
{
	return( adpb_rowsetup_to( n_dynamicmtx, cpmx1, cpmx2, i, rs + 4*i, rc + 4*i, rn + i, rsu + ks*i, rmask + i, L ) );
}

static void adpb_row( adpb_t *p, int b, int i )
{
	/* work/v8s: with localsetup every band makes row i's letter selection itself, into thread-local
	   scratch (the same function on the same inputs, so the same values): no band reads another
	   core's selection arrays (rs, rc, rn, rsu, rmask: 4-5 cache lines per row) */
	static TLS double l_rs[4], l_rsu[DENSE_MAXU];
	static TLS const double *l_rc[4];
	static TLS int l_rn;
	static TLS unsigned int l_mask;
	const double *rsp, *rsup;
	const double * const *rcp;
	unsigned int rmaskv;
	int nselv;
	if( p->localsetup )
	{
		adpb_rowsetup_to( p->dynmtx, p->cpmx1, p->cpmx2, i, l_rs, l_rc, &l_rn, l_rsu, &l_mask, &p->letters );
		rsp = l_rs; rcp = l_rc; rsup = l_rsu; rmaskv = l_mask; nselv = l_rn;
	}
	else
	{
		if( b == 0 && p->lazyrows ) // work/prog3: band 0 sees every row first (the others wait for it)
		{
		/* work/v8s: the selections are made by every participant in chunks before it takes a band
		   (adpb_run); band 0 waits until row i's chunk is published.  It used to make them itself, one
		   row at a time, which made band 0 the slowest stage of the pipeline (~300 ns of letter
		   scores per row against ~400 ns for a 128-column band row). */
		const int k = ( i - 1 ) / p->sch;
		int sp = 0;
		while( !__atomic_load_n( p->ready + k, __ATOMIC_ACQUIRE ) )
			adpb_wait( &sp );
		}
		rsp = p->rs + 4*i; rcp = p->rc + 4*i; rsup = p->rsu + p->ks*i; rmaskv = p->rmask[i]; nselv = p->rn[i];
	}
#ifdef PROGTIMING
	if( b == 0 && v8s_stat_on() ) v8s_stat_row( 1, nselv > 4 ? 2 : 1, nselv, p->lgth2, rmaskv, p->bm );
#endif
	{
	const int lgth2 = p->lgth2, j0 = p->J[b], j1 = p->J[b+1], jm = MIN( j1, lgth2 );
	double * restrict prew = p->wb[(i-1)&1];
	double * restrict curw = p->wb[i&1];
	double * restrict mloc = p->m;
	int * restrict mploc = p->mp;
#ifdef PROGTIMING
	/* work/v8s: MAFFT_TILEIJPTEST=1 (timing builds, DIAGNOSIS ONLY, the output is then wrong): every band
	   writes its traceback codes into a private scratch row instead of the shared matrix, to measure
	   what the shared matrix's cache lines and pages cost the tiled DP */
	static TLS int *ijptest_row = NULL; static TLS int ijptest_cap = 0;
	int * restrict ijpi;
	if( adpb_ijptest )
	{
		if( ijptest_cap < p->lgth2 + 2 ) { free( ijptest_row ); ijptest_cap = p->lgth2 + 2 + p->lgth2 / 2; ijptest_row = malloc( sizeof( int ) * ijptest_cap ); }
		ijpi = ijptest_row;
	}
	else
		ijpi = p->bandmajor ? p->ijb[b] + (size_t)i * p->ijw[b] - j0 : p->ijp[i];
#else
	int * restrict ijpi = p->bandmajor ? p->ijb[b] + (size_t)i * p->ijw[b] - j0 : p->ijp[i]; // work/v8s: band-major (ijt_cell)
#endif
	const double * restrict gf2 = p->gf2;
	const double * restrict fg2 = p->fgcp2;
	const double * restrict og2 = p->ogcp2;
	const int nsel = nselv, gen = ( nsel > 4 ); // work/prog3: gen: match scores formed in curw first
	const double * restrict c0 = rcp[0], * restrict c1 = rcp[1], * restrict c2 = rcp[2], * restrict c3 = rcp[3];
	const double s0 = rsp[0], s1 = rsp[1], s2 = rsp[2], s3 = rsp[3];
	const double fgcp1va = p->fgcp1[i-1], ogcp1va = p->ogcp1[i], gf1va = p->gf1[i], gf1vapre = p->gf1[i-1];
	const double fpenalty_ex = p->fpenalty_ex;
	double mil, wmv = 0.0;
	int mpil, ijv, j;

	if( gen && jm > j0 ) // work/prog3: the match scores of the band's columns (< lgth2) first
		Adp_mrow( j0, jm, rmaskv, rsup, p->cu, p->bm, curw );
	if( b == 0 )
	{
		prew[0] = p->initverticalw[i-1];                       // previousw[0] = initverticalw[i-1];
		curw[0] = p->initverticalw[i];                         // currentw[0] = initverticalw[i];
		mil = prew[0] + p->ogcp2[1] * p->gf1[i-1]; mpil = 0;   // mi = previousw[0] + ogcp2[1] * gapfreq1pt[i-1]; mpi=0;
		j = 1;
		if( 1 < lgth2 && !gen ) ADP_CELL( 1, 1 ) else ADP_CELL( 1, 0 )
		j = 2;
	}
	else
	{
		const int nr1 = p->nrows + 1;
		mil = p->bmi[b*nr1+i];
		mpil = p->bmpi[b*nr1+i];
		j = j0;
		if( j < lgth2 && !gen ) ADP_CELL_PW( 0, 1, p->bw[b*nr1+i-1] ) else ADP_CELL_PW( 0, 0, p->bw[b*nr1+i-1] )
		j = j0 + 1;
	}
	if( jm > j ) // work/prog3: the band's inner cells, vectorized
	{
		double scu4[4]; scu4[0] = s0; scu4[1] = s1; scu4[2] = s2; scu4[3] = s3;
		Adp_cells( gen, i, j, jm, mloc, mploc, curw, ijpi, prew, gf2, fg2, og2, c0, c1, c2, c3, scu4,
			fgcp1va, ogcp1va, gf1va, gf1vapre, fpenalty_ex, &mil, &mpil, &wmv );
		j = jm;
	}
	if( j == lgth2 && j < j1 )
	{
		ADP_CELL( 0, 0 )
		j++;
	}
	if( b < p->nb - 1 )
	{
		const int nr1 = p->nrows + 1;
		p->bw[(b+1)*nr1+i] = curw[j1-1];
		p->bmi[(b+1)*nr1+i] = mil;
		p->bmpi[(b+1)*nr1+i] = mpil;
	}
	else
		p->band[b].wmlast = wmv; // work/v8s (was p->wmlast)
	if( j0 <= lgth2-1 && lgth2-1 < j1 ) p->lastverticalw[i] = curw[lgth2-1]; // lastverticalw[i] = currentw[lgth2-1];
}
}

/* work/v8s: the participants first make the per-row letter selections of the tiled DP in chunks
   (lazy rows), then work on bands.  Band policy (sticky, default): a thread keeps the band it took
   and processes every row its predecessor has finished, re-reading the predecessor's progress only
   when it has caught up; only then does it look for another band that has rows ready and no owner
   (earliest anti-diagonal first), and if there is none it waits on its own band.  The earlier policy
   gave the band up after at most R rows or whatever the predecessor had done (often a few rows) and
   re-picked the earliest anti-diagonal: on CARC (EPYC 7513, 32 threads, 16S.T) 53% of ~3.5 million
   claims changed band, a claim covered ~900 cells, and every change pulled the band's columns
   (m, mp, the two row buffers) from another core's cache: the tiled DPs ran at 18.5 ns per cell and
   thread against 2.8 ns serial.  Which thread computes a row of a band, and when, never changes a
   value: each cell is the same ADP_CELL code on the same operands, ordered by the done counters as
   before.  MAFFT_TILEOLD=1 (timing builds) restores the old policy (A/B). */
static int adpb_pick( adpb_t *p, int R )
{
	const int nb = p->nb, nrows = p->nrows;
	int b, best = -1, bestkey = 0, alldone = 1;
	for( b=0; b<nb; b++ )
	{
		int d = __atomic_load_n( &p->band[b].done, __ATOMIC_ACQUIRE ), avail, key;
		if( d >= nrows ) continue;
		alldone = 0;
		if( __atomic_load_n( &p->band[b].busy, __ATOMIC_RELAXED ) ) continue;
		avail = b ? __atomic_load_n( &p->band[b-1].done, __ATOMIC_ACQUIRE ) : nrows;
		if( avail <= d || ( avail - d < p->batch && avail < nrows ) ) continue; // work/v8s: a batch of rows, or the rest
		key = d / R + b; // wavefront order: earliest anti-diagonal first
		if( best < 0 || key < bestkey ) { best = b; bestkey = key; }
	}
	return( alldone ? -2 : best );
}

/* every unfinished band has an owner (racy; a helper that leaves on a wrong answer only costs time) */
static int adpb_allowned( adpb_t *p )
{
	int b;
	for( b=0; b<p->nb; b++ )
		if( __atomic_load_n( &p->band[b].done, __ATOMIC_ACQUIRE ) < p->nrows && !__atomic_load_n( &p->band[b].busy, __ATOMIC_RELAXED ) ) return( 0 );
	return( 1 );
}
static void adpb_run_core( adpb_t *p, int helper );
static void adpb_run( void *arg ) { adpb_run_core( (adpb_t *)arg, 1 ); } // a helper (dphelp job)
static void adpb_run_core( adpb_t *p, int helper )
{
	const int nrows = p->nrows, R = p->R;
	int spins = 0, cur = -1;
#ifdef PROGTIMING
	const int stat = v8s_stat_on();
	double st_t0 = stat ? pt_now() : 0.0, st_busy = 0.0, st_claims = 0.0, st_migr = 0.0;
	int st_last = -1;
#endif
	if( p->lazyrows && !p->localsetup ) // the letter selections, chunk by chunk (any thread, any order; published per chunk)
	{
		int k;
		while( ( k = __atomic_fetch_add( &p->setup_next, p->sch, __ATOMIC_RELAXED ) ) < nrows )
		{
			int i, hi = MIN( k + p->sch, nrows );
			for( i=k+1; i<=hi; i++ )
				adpb_rowsetup( p->dynmtx, p->cpmx1, p->cpmx2, i, p->ks, (double *)p->rs, (const double **)p->rc, (int *)p->rn, (double *)p->rsu, (unsigned int *)p->rmask, &p->letters );
			__atomic_store_n( p->ready + k / p->sch, 1, __ATOMIC_RELEASE );
		}
	}
	while( 1 )
	{
		int d, avail, end, i;
		if( cur < 0 )
		{
			int best = adpb_pick( p, R ), zero = 0;
			if( best == -2 ) break;
			if( best < 0 )
			{
				/* work/v8s: a helper with nothing to take while every unfinished band has an owner goes
				   back to the pool (the poster stays until all bands are done) */
				if( helper && adpb_allowned( p ) ) break;
				adpb_wait( &spins );
				continue;
			}
			if( !__atomic_compare_exchange_n( &p->band[best].busy, &zero, 1, 0, __ATOMIC_ACQUIRE, __ATOMIC_RELAXED ) ) continue;
			spins = 0;
			cur = best;
		}
		d = __atomic_load_n( &p->band[cur].done, __ATOMIC_ACQUIRE );
		if( d >= nrows ) // finished: give the band up
		{
			__atomic_store_n( &p->band[cur].busy, 0, __ATOMIC_RELEASE );
			cur = -1;
			continue;
		}
		avail = cur ? __atomic_load_n( &p->band[cur-1].done, __ATOMIC_ACQUIRE ) : nrows;
		if( avail > d && ( avail - d >= p->batch || avail >= nrows ) ) // work/v8s: wait for a batch of rows (fewer hand-offs across cores)
		{
			end = p->sticky ? MIN( avail, d + p->claimcap ) : MIN( avail, d + R ); // work/v8s: at most claimcap rows per claim
#ifdef PROGTIMING
			{ double _r0 = stat ? pt_now() : 0.0;
#endif
			for( i=d+1; i<=end; i++ )
			{
				adpb_row( p, cur, i );
				__atomic_store_n( &p->band[cur].done, i, __ATOMIC_RELEASE );
			}
#ifdef PROGTIMING
			if( stat ) { st_busy += pt_now() - _r0; st_claims += 1.0; if( cur != st_last ) st_migr += 1.0; st_last = cur; } }
#endif
			spins = 0;
			/* work/v8s: after claimcap rows the band is put back and the thread picks again (earliest
			   anti-diagonal first): a thread that runs slowly (another process on its core, an
			   efficiency core) holds a band, and every band behind it, for at most claimcap rows */
			if( p->sticky && end < d + p->claimcap ) continue;
			__atomic_store_n( &p->band[cur].busy, 0, __ATOMIC_RELEASE ); // old policy: one claim, then re-pick
			cur = -1;
			continue;
		}
		/* caught up with the predecessor: another band with rows ready and no owner, if any (work/v8s:
		   looked for on every 16th check only; scanning reads every band's line) */
		if( ( spins & 15 ) == 0 )
		{
			int best = adpb_pick( p, R ), zero = 0;
			if( best >= 0 && __atomic_compare_exchange_n( &p->band[best].busy, &zero, 1, 0, __ATOMIC_ACQUIRE, __ATOMIC_RELAXED ) )
			{
				__atomic_store_n( &p->band[cur].busy, 0, __ATOMIC_RELEASE );
				cur = best;
				spins = 0;
				continue;
			}
		}
		adpb_wait( &spins ); // work/v8s
	}
#ifdef PROGTIMING
	if( stat )
	{
		v8s_add( &p->st_busy, st_busy ); v8s_add( &p->st_part, pt_now() - st_t0 ); v8s_add( &p->st_claims, st_claims );
		v8s_add( &p->st_nthr, 1.0 ); v8s_add( &p->st_migr, st_migr );
	}
#endif
}

/* Returns 1 if the DP of rows 1..lasti-1 was done here (tiled, possibly with helpers), 0 if the
   caller must run its own loop (row not fusable, too small, or no helper idle).  On success the
   final rows are in w0/w1 exactly as the serial loop leaves them and *wmp is its final wm. */
static TLS double *adpb_rs = NULL;
static TLS const double **adpb_rc = NULL;
static TLS int *adpb_rn = NULL;
static TLS double *adpb_rsu = NULL;
static TLS unsigned int *adpb_rmask = NULL;
static TLS double *adpb_cu[DENSE_MAXU];
static TLS double *adpb_bw = NULL, *adpb_bmi = NULL;
static TLS int *adpb_bmpi = NULL;
static TLS int adpb_rowcap = 0, adpb_bcap = 0;
static int adpb_minwidth = -1, adpb_mincells = -1, adpb_force = 0, adpb_widthset = 0;
static int adpb_sticky = 1, adpb_nbfac = 1; // work/v8s: band policy (adpb_run) and bands per expected participant
static int adpb_align = 1, adpb_batch = 8;  // work/v8s: 16-column band boundaries; rows a follower waits for (adpb_run)
static int adpb_claimcap = 256;             // work/v8s: rows per claim of a band (adpb_run)
static int adpb_critnb = 8, adpb_critmincells = 4000000; // work/v8s: critical tiling (bands, minimum cells)
static int adpb_localsetup = 0; // work/v8s: see adpb_row
static int adpb_bandmajor = 1;  // work/v8s: band-major traceback (ijt_cell); MAFFT_TILEBANDMAJOR=0 (timing builds) turns it off
static size_t A_ijp_capacity( void ); // (after A_ijp_get)
static TLS int *adpb_ready = NULL;          // work/v8s: per selection chunk (adpb_t.ready)
static TLS int adpb_readycap = 0;

static void adpb_free( void )
{
	free( adpb_rs ); adpb_rs = NULL;
	free( adpb_rc ); adpb_rc = NULL;
	free( adpb_rn ); adpb_rn = NULL;
	free( adpb_rsu ); adpb_rsu = NULL;
	free( adpb_rmask ); adpb_rmask = NULL;
	free( gen_bm ); gen_bm = NULL; gen_bmcap = 0;
	free( adpb_bw ); adpb_bw = NULL;
	free( adpb_bmi ); adpb_bmi = NULL;
	free( adpb_bmpi ); adpb_bmpi = NULL;
	free( adpb_ready ); adpb_ready = NULL; adpb_readycap = 0;
	adpb_rowcap = adpb_bcap = 0;
}

static int Adp_tiled( int lgth1, int lgth2, int lasti, double **n_dynamicmtx, double **cpmx1, double **cpmx2,
	double *w0, double *w1, double *m, int *mp, int **ijp,
	const double *gf1, const double *gf2, const double *fgcp1, const double *ogcp1, const double *fgcp2, const double *ogcp2,
	const double *initverticalw, double *lastverticalw, double fpenalty_ex, double *wmp )
{
	adpb_t *p;
	dphelpjob_t job;
	int nrows = lasti - 1, nb, b, i, idle, maxb, ks, minw, crit;
	adpb_letters_t letters;
	double cells = (double)nrows * lgth2;

	if( adpb_minwidth < 0 )
	{
#ifdef PROGTIMING
		char *e; // tuning and test switches of the timing builds
		adpb_widthset = ( getenv( "MAFFT_TILEWIDTH" ) != NULL );
		adpb_minwidth = ( e = getenv( "MAFFT_TILEWIDTH" ) ) ? atoi( e ) : 384;
		adpb_mincells = ( e = getenv( "MAFFT_TILEMINCELLS" ) ) ? atoi( e ) : 2000000;
		adpb_force = getenv( "MAFFT_TILEFORCE" ) != NULL; // testing: tile even without helpers
		if( getenv( "MAFFT_NOTILE" ) ) adpb_minwidth = 1 << 30;
		if( getenv( "MAFFT_TILEOLD" ) ) { adpb_sticky = 0; adpb_nbfac = 3; } // work/v8s: the old band policy (A/B)
		if( ( e = getenv( "MAFFT_TILENBFAC" ) ) ) adpb_nbfac = atoi( e );
		if( ( e = getenv( "MAFFT_TILEALIGN" ) ) ) adpb_align = atoi( e );
		if( ( e = getenv( "MAFFT_TILEBATCH" ) ) ) adpb_batch = atoi( e );
		if( ( e = getenv( "MAFFT_TILECLAIM" ) ) ) adpb_claimcap = atoi( e );
		if( ( e = getenv( "MAFFT_CRITNB" ) ) ) adpb_critnb = atoi( e );
		adpb_ijptest = ( getenv( "MAFFT_TILEIJPTEST" ) != NULL );
		if( ( e = getenv( "MAFFT_TILELOCALSETUP" ) ) ) adpb_localsetup = atoi( e );
		if( ( e = getenv( "MAFFT_TILEBANDMAJOR" ) ) ) adpb_bandmajor = atoi( e );
		if( ( e = getenv( "MAFFT_CRITMINCELLS" ) ) ) adpb_critmincells = atoi( e );
#else
		adpb_minwidth = 384;
		adpb_mincells = 2000000;
		adpb_force = 0;
#endif
	}
	if( nrows < 2 || cells < adpb_mincells || dense_nu < 0 ) { V8S_NOTILE( 0, cells ); return( 0 ); }
	idle = dphelp_idle();
	/* work/v8s: critical tiling.  A large DP of a step on the critical path is tiled even when (almost)
	   no thread is idle; its bands are posted for threads that finish their steps (they serve it
	   before ready steps of lower priority, dphelp.h).  Without joiners the poster runs the bands one
	   after the other: the same cells, a little more per-row overhead. */
	crit = ( idle < adpb_critnb - 1 && cells >= adpb_critmincells && dphelp_iscritical() );
	if( idle <= 0 && !adpb_force && !crit ) { V8S_NOTILE( 1, cells ); return( 0 ); }
	/* work/prog3: narrower bands when many threads are idle (the merges at the top of the tree, which
	   often have only ~1500 columns on the band side: 384-column bands allowed only 4 participants) */
	minw = adpb_minwidth;
	/* work/v8s: narrower bands (48 / 96 / 192 columns; prog3 had 128 / 192 / 384).  On CARC (32 threads,
	   EPYC 7513, the job's cores spread over both sockets) a band row costs ~3x the serial cell rate
	   whatever the band width, so more bands, i.e. more threads, shorten the tiled DPs: the tiled DPs
	   of 16S.T took 2.3-2.5 s in all with 48-column bands, 3.1-3.3 s with 64-96, 4.0-4.5 s with 128 */
	if( !adpb_widthset ) { if( idle >= 12 ) minw = 48; else if( idle >= 4 || crit ) minw = 96; else minw = 192; }
	if( lgth2 < 2 * minw ) { V8S_NOTILE( 2, cells ); return( 0 ); }

	maxb = MIN( ADPB_MAXB, lgth2 / minw );
	nb = MIN( maxb, adpb_nbfac * ( 1 + MAX( idle, 1 ) ) ); // work/v8s: one band per expected participant (was 3, see adpb_run)
	if( crit ) nb = MAX( nb, MIN( maxb, adpb_critnb ) );
	if( nb < 2 ) { V8S_NOTILE( 3, cells ); return( 0 ); }

	ks = MAX( 1, dense_nu );
	if( nrows + 1 > adpb_rowcap )
	{
		free( adpb_rs ); free( adpb_rc ); free( adpb_rn ); free( adpb_rsu ); free( adpb_rmask );
		adpb_rowcap = nrows + 1 + ( nrows + 1 ) / 2;
		adpb_rs = malloc( sizeof( double ) * 4 * adpb_rowcap );
		adpb_rc = malloc( sizeof( double * ) * 4 * adpb_rowcap );
		adpb_rn = malloc( sizeof( int ) * adpb_rowcap );
		adpb_rsu = malloc( sizeof( double ) * DENSE_MAXU * adpb_rowcap );
		adpb_rmask = malloc( sizeof( unsigned int ) * adpb_rowcap );
	}
	letters.nu = dense_nu; letters.n1 = dense1_n;
	memcpy( letters.u, dense_u, sizeof( int ) * dense_nu );
	memcpy( letters.k1, dense1_k, sizeof( int ) * dense1_n );
	/* the row scores and profile columns of each row, selected exactly as in the serial loop
	   (work/prog3: rows with more than 4 letters take the general kernel, as rowfused == 2 there).
	   work/prog3: they are formed by band 0 just before it processes the row (adpb_rowsetup), off the
	   setup path; only the timing builds' A/B switch MAFFT_NOGENERAL needs them all in advance. */
	{
		int u;
		if( !adp_general_on() )
			for( i=1; i<=nrows; i++ )
				if( adpb_rowsetup( n_dynamicmtx, cpmx1, cpmx2, i, ks, adpb_rs, adpb_rc, adpb_rn, adpb_rsu, adpb_rmask, &letters ) > 4 ) { V8S_NOTILE( 4, cells ); return( 0 ); }
		if( dense_nu > 4 ) // some row may take the general kernel
		{
			for( u=0; u<dense_nu; u++ ) adpb_cu[u] = cpmx2[dense_u[u]];
			if( lgth2 / 4 + 2 > gen_bmcap ) { free( gen_bm ); gen_bmcap = lgth2 / 4 + 2 + lgth2 / 8; gen_bm = malloc( sizeof( unsigned int ) * gen_bmcap ); }
			Adp_blockmasks( gen_bm, (const double * const *)adpb_cu, dense_nu, lgth2 );
		}
	}
	if( (size_t)nb * ( nrows + 1 ) > (size_t)adpb_bcap )
	{
		free( adpb_bw ); free( adpb_bmi ); free( adpb_bmpi );
		adpb_bcap = nb * ( nrows + 1 ) + nb * ( nrows + 1 ) / 2;
		adpb_bw = malloc( sizeof( double ) * adpb_bcap );
		adpb_bmi = malloc( sizeof( double ) * adpb_bcap );
		adpb_bmpi = malloc( sizeof( int ) * adpb_bcap );
	}

	if( posix_memalign( (void **)&p, 64, sizeof( adpb_t ) ) ) { V8S_NOTILE( 5, cells ); return( 0 ); }
#ifdef PROGTIMING
	if( v8s_stat_on() ) { v8s_add( &v8s_tiled[0], 1.0 ); v8s_add( &v8s_tiled[1], cells ); }
#endif
	memset( p, 0, sizeof( adpb_t ) );
	p->nrows = nrows; p->lgth2 = lgth2; p->nb = nb;
	for( b=0; b<=nb; b++ ) p->J[b] = 1 + (int)( (long long)lgth2 * b / nb );
	if( adpb_align ) // work/v8s: inner boundaries on multiples of 16 columns (no shared cache lines, see adp_avec)
	{
		for( b=1; b<nb; b++ )
		{
			int jb = ( p->J[b] + 8 ) & ~15;
			if( jb <= p->J[b-1] ) jb = p->J[b-1] + 16;
			p->J[b] = jb;
		}
		for( b=1; b<nb; b++ ) if( p->J[b] >= p->J[b+1] ) { for( b=0; b<=nb; b++ ) p->J[b] = 1 + (int)( (long long)lgth2 * b / nb ); break; } // (cannot happen with bands >= 64 columns)
	}
	p->R = MAX( 8, 24576 / ( lgth2 / nb ) );
	p->wb[0] = w0; p->wb[1] = w1;
	p->m = m; p->mp = mp; p->ijp = ijp;
	p->gf1 = gf1; p->gf2 = gf2; p->fgcp1 = fgcp1; p->ogcp1 = ogcp1; p->fgcp2 = fgcp2; p->ogcp2 = ogcp2;
	p->initverticalw = initverticalw; p->lastverticalw = lastverticalw;
	p->rs = adpb_rs; p->rc = adpb_rc; p->rn = adpb_rn; p->ks = ks;
	p->rsu = adpb_rsu; p->rmask = adpb_rmask; p->cu = (const double * const *)adpb_cu; p->bm = gen_bm;
	p->lazyrows = adp_general_on(); p->dynmtx = n_dynamicmtx; p->cpmx1 = cpmx1; p->cpmx2 = cpmx2;
	p->letters = letters;
	p->fpenalty_ex = fpenalty_ex;
	p->bw = adpb_bw; p->bmi = adpb_bmi; p->bmpi = adpb_bmpi;
	p->sticky = adpb_sticky; // work/v8s
	p->batch = adpb_sticky ? MAX( 1, adpb_batch ) : 1;
	p->claimcap = MAX( p->batch, adpb_claimcap );
	p->localsetup = ( adpb_localsetup && p->lazyrows );
	/* work/v8s: band-major traceback regions (see ijt_cell), carved from the traceback block that
	   A_ijp_get made for this DP: rows 0 .. nrows+1 of every band, each region starting on a 64-byte
	   line; falls back to the row-major matrix if the block were too small (cannot happen: the block
	   holds ( lgth1 + 2 ) rows of >= lgth2 + 2 ints) */
	p->bandmajor = 0;
	if( adpb_bandmajor )
	{
		size_t off = 0, cap = A_ijp_capacity();
		int *blk = ijp[0];
		for( b=0; b<nb; b++ )
		{
			p->ijw[b] = p->J[b+1] - p->J[b];
			p->ijb[b] = blk + off;
			off += (size_t)( nrows + 2 ) * p->ijw[b];
			off = ( off + 15 ) & ~(size_t)15;
		}
		p->bandmajor = ( off <= cap && ( (size_t)blk & 63 ) == 0 );
	}
	p->sch = 64; // rows per selection chunk
	p->setup_next = 0;
	if( nrows / p->sch + 2 > adpb_readycap ) { free( adpb_ready ); adpb_readycap = nrows / p->sch + 2 + nrows / p->sch / 2; adpb_ready = malloc( sizeof( int ) * adpb_readycap ); }
	memset( adpb_ready, 0, sizeof( int ) * ( nrows / p->sch + 2 ) );
	p->ready = adpb_ready;
	for( b=1; b<nb; b++ ) p->bw[b*(nrows+1)+0] = w0[p->J[b]-1]; // row 0, last column of band b-1
	p->band[nb-1].wmlast = *wmp;

	job.run = adpb_run;
	job.arg = p;
#ifdef PROGTIMING
	{ double _e0 = pt_now();
#endif
	if( crit ) dphelp_post_crit( &job, nb - 1 ); else dphelp_post( &job, nb - 1 );
	adpb_run_core( p, 0 ); // the poster stays until every band is done
	dphelp_finish( &job );
	if( p->bandmajor ) // work/v8s: Atracking of this A__align call reads the band-major regions
	{
		atr_tile.nb = nb;
		for( b=0; b<=nb; b++ ) atr_tile.J[b] = p->J[b];
		for( b=0; b<nb; b++ ) { atr_tile.base[b] = p->ijb[b]; atr_tile.w[b] = p->ijw[b]; }
		atr_tiled = 1;
	}
#ifdef PROGTIMING
	if( adpb_ijptest ) { for( i=1; i<=nrows; i++ ) memset( ijp[i], 0, sizeof( int ) * ( lgth2 + 1 ) ); atr_tiled = 0; } // (diagnosis: a harmless traceback)
	if( v8s_stat_on() && crit ) { v8s_add( &v8s_tile[14], 1.0 ); v8s_add( &v8s_tile[15], p->st_nthr ); }
	if( v8s_stat_on() )
	{
		double el = pt_now() - _e0;
		v8s_add( &v8s_tile[0], el ); v8s_add( &v8s_tile[1], p->st_part ); v8s_add( &v8s_tile[2], p->st_busy );
		v8s_add( &v8s_tile[3], p->st_claims ); v8s_add( &v8s_tile[4], p->st_nthr ); v8s_add( &v8s_tile[5], p->st_migr );
		v8s_add( &v8s_tile[6], (double)nb ); v8s_add( &v8s_tile[7], el * nb );
		if( cells >= 2.0e7 ) { v8s_add( &v8s_tile[8], el ); v8s_add( &v8s_tile[9], p->st_busy ); v8s_add( &v8s_tile[10], p->st_part ); v8s_add( &v8s_tile[11], 1.0 ); v8s_add( &v8s_tile[12], p->st_nthr ); v8s_add( &v8s_tile[13], (double)nb ); }
	} }
#endif
	*wmp = p->band[nb-1].wmlast;
	free( p );
	return( 1 );
}
#endif

#ifdef PROGTIMING
static TLS double *chk_m = NULL, *chk_cur = NULL, *chk_match = NULL;
static TLS int *chk_mp = NULL, *chk_ijp = NULL;
static TLS int pt_chkrow_pending = 0;
static TLS double pt_chk_mi, pt_chk_wm, pt_chk_t;
static TLS int pt_chk_mpi;
static int pt_dpcheck_v = -1;
static int pt_dpcheck( void ) { if( pt_dpcheck_v < 0 ) pt_dpcheck_v = ( getenv( "MAFFT_DPCHECK" ) != NULL ); return( pt_dpcheck_v ); }
#endif

#if defined( PROGTIMING ) && defined( enablemultithread )
/* work/prog2: MAFFT_TILECHECK=1 -- for every tiled DP, save the initial state, keep the tiled
   result, rerun the original row loop from the initial state and compare every value it
   produces (both row buffers including the unused column lgth2, m, mp, the traceback cells,
   lastverticalw, wm) bitwise.  Exits on the first difference. */
static int pt_tilecheck_v = -1;
static int pt_tilecheck( void ) { if( pt_tilecheck_v < 0 ) pt_tilecheck_v = ( getenv( "MAFFT_TILECHECK" ) != NULL ); return( pt_tilecheck_v ); }
static TLS int tilechk_pending = 0;
static TLS double *tilechk_cur0, *tilechk_prev0;          // the two physical row buffers
static TLS double *tc_i0, *tc_i1, *tc_im, *tc_r0, *tc_r1, *tc_rm, *tc_rlv;
static TLS int *tc_imp, *tc_rmp, *tc_rijp;
static TLS double tc_iwm, tc_rwm;
static void tilechk_save_init( int lgth2, double *cur, double *prev, double *m, int *mp, double wm )
{
	size_t n = lgth2 + 1;
	tilechk_cur0 = cur; tilechk_prev0 = prev;
	tc_i0 = realloc( tc_i0, n * sizeof( double ) ); memcpy( tc_i0, cur, n * sizeof( double ) );
	tc_i1 = realloc( tc_i1, n * sizeof( double ) ); memcpy( tc_i1, prev, n * sizeof( double ) );
	tc_im = realloc( tc_im, n * sizeof( double ) ); memcpy( tc_im, m, n * sizeof( double ) );
	tc_imp = realloc( tc_imp, n * sizeof( int ) ); memcpy( tc_imp, mp, n * sizeof( int ) );
	tc_iwm = wm;
	tilechk_pending = 0;
}
static void tilechk_save_result( int lgth2, int lasti, double *m, int *mp, int **ijp, double *lastverticalw, double wm )
{
	size_t n = lgth2 + 1;
	int i;
	tc_r0 = realloc( tc_r0, n * sizeof( double ) ); memcpy( tc_r0, tilechk_cur0, n * sizeof( double ) );
	tc_r1 = realloc( tc_r1, n * sizeof( double ) ); memcpy( tc_r1, tilechk_prev0, n * sizeof( double ) );
	tc_rm = realloc( tc_rm, n * sizeof( double ) ); memcpy( tc_rm, m, n * sizeof( double ) );
	tc_rmp = realloc( tc_rmp, n * sizeof( int ) ); memcpy( tc_rmp, mp, n * sizeof( int ) );
	tc_rlv = realloc( tc_rlv, (size_t)lasti * sizeof( double ) ); memcpy( tc_rlv, lastverticalw, (size_t)lasti * sizeof( double ) );
	tc_rijp = realloc( tc_rijp, (size_t)lasti * lgth2 * sizeof( int ) );
	if( atr_tiled ) { int j; for( i=1; i<lasti; i++ ) for( j=1; j<=lgth2; j++ ) tc_rijp[(size_t)i * lgth2 + j-1] = *ijt_cell( &atr_tile, i, j ); } // work/v8s: band-major
	else
	for( i=1; i<lasti; i++ ) memcpy( tc_rijp + (size_t)i * lgth2, ijp[i] + 1, lgth2 * sizeof( int ) );
	tc_rwm = wm;
	tilechk_pending = 1;
}
static void tilechk_restore_init( int lgth2, double *m, int *mp, double *wm )
{
	size_t n = lgth2 + 1;
	memcpy( tilechk_cur0, tc_i0, n * sizeof( double ) );
	memcpy( tilechk_prev0, tc_i1, n * sizeof( double ) );
	memcpy( m, tc_im, n * sizeof( double ) );
	memcpy( mp, tc_imp, n * sizeof( int ) );
	*wm = tc_iwm;
}
static void tilechk_compare( int lgth1, int lgth2, int lasti, double *m, int *mp, int **ijp, double *lastverticalw, double wm )
{
	size_t n = lgth2 + 1;
	int i, bad = 0;
	if( memcmp( tc_r0, tilechk_cur0, n * sizeof( double ) ) ) bad |= 1;
	if( memcmp( tc_r1, tilechk_prev0, n * sizeof( double ) ) ) bad |= 2;
	if( memcmp( tc_rm, m, n * sizeof( double ) ) ) bad |= 4;
	if( memcmp( tc_rmp, mp, n * sizeof( int ) ) ) bad |= 8;
	if( memcmp( tc_rlv + 1, lastverticalw + 1, (size_t)( lasti - 1 ) * sizeof( double ) ) ) bad |= 16;
	if( memcmp( &tc_rwm, &wm, sizeof( double ) ) ) bad |= 32;
	for( i=1; i<lasti; i++ ) if( memcmp( tc_rijp + (size_t)i * lgth2, ijp[i] + 1, lgth2 * sizeof( int ) ) ) { bad |= 64; break; }
	if( bad )
	{
		reporterr( "\nTILECHECK MISMATCH (code %d) in A__align, lgth1=%d lgth2=%d lasti=%d row %d\n", bad, lgth1, lgth2, lasti, i );
		exit( 1 );
	}
	pt_acc[31] += 1.0;
	tilechk_pending = 0;
}
#endif

/* work/prog2: the profile of one group in a single pass over its rows.  Replaces
     cpmx_calc_new( seq, cpmx, eff, len, clus )
     st_OpeningGapCount( ogcp, ... ) or new_OpeningGapCount( ogcp, ..., sgappat )  (st: sgappat == NULL)
     st_FinalGapCount( fgcp, ... ) or new_FinalGapCount( fgcp, ..., egappat )
     gapcountf( gapfreq, seq, clus, eff, len )
   Each of those zeroes its array and then adds eff[j] to entries member by member (j ascending);
   here the four arrays receive exactly the same additions, member by member in the same order,
   while each row is read once.  Opening gaps: column k gains eff[j] when it is a gap and the
   previous column (or, for k = 0, sgappat[j], or no gap for the st_ version) is not.  Final
   gaps: column k gains eff[j] when it is a gap and column k+1 is not; column len is the
   terminator, never a gap, so the egappat argument of new_FinalGapCount never matters (its
   last test compares against the terminator's state).  st_OpeningGapCount also sets
   ogcp[len] = 0.0.  Rows must be terminated at len, as A__align ensures.  Large profiles are split
   by columns among idle helper threads (dphelp_for). */
typedef struct
{
	char **seq; int clus; double *eff; int len;
	double **cpmx; double *ogcp, *fgcp, *gapfreq;
	const char *sgappat;
} profctx_t;

/* columns lo..hi-1 of all four arrays (every addition to an entry comes from its own column range,
   member by member in ascending order, so ranges can be done by different threads) */
static void profile_fused_cols( void *vctx, int lo, int hi )
{
	profctx_t *c = (profctx_t *)vctx;
	char **seq = c->seq; double *eff = c->eff; double **cpmx = c->cpmx;
	double *ogcp = c->ogcp, *fgcp = c->fgcp, *gapfreq = c->gapfreq;
	const char *sgappat = c->sgappat;
	int j, k;
	for( k=0; k<nalphabets; k++ ) { double *cp = cpmx[k]; int i; for( i=lo; i<hi; i++ ) cp[i] = 0.0; }
	for( k=lo; k<hi; k++ ) { ogcp[k] = 0.0; fgcp[k] = 0.0; gapfreq[k] = 0.0; }
	for( j=0; j<c->clus; j++ )
	{
		const unsigned char *sp = (const unsigned char *)seq[j];
		const double e = eff[j];
		int prev = lo ? ( sp[lo-1] == '-' ) : ( sgappat && sgappat[j] == '-' ); // opening: previous column, or the boundary pattern
		int g = ( sp[lo] == '-' ), gn;
		for( k=lo; k<hi; k++ )
		{
			gn = ( sp[k+1] == '-' ); // column len is the terminator: never a gap
			cpmx[(unsigned char)amino_n[sp[k]]][k] += e;
			if( g )
			{
				gapfreq[k] += e;
				if( !prev ) ogcp[k] += e;
				if( !gn ) fgcp[k] += e;
			}
			prev = g;
			g = gn;
		}
	}
}

static void profile_fused( char **seq, int clus, double *eff, int len, double **cpmx, double *ogcp, double *fgcp, double *gapfreq, const char *sgappat )
{
	profctx_t c;
	c.seq = seq; c.clus = clus; c.eff = eff; c.len = len; c.cpmx = cpmx; c.ogcp = ogcp; c.fgcp = fgcp; c.gapfreq = gapfreq; c.sgappat = sgappat;
	if( sgappat == NULL ) ogcp[len] = 0.0; // st_OpeningGapCount
	if( len <= 0 ) return;
#ifdef enablemultithread
	dphelp_for( len, MAX( 256, len / 64 ), (double)clus * len, 4.0e6, profile_fused_cols, &c ); // work/v8s: was MAX( 512, len / 16 ): up to 64 claims (32 threads)
#else
	profile_fused_cols( &c, 0, len );
#endif
}

#ifdef PROGTIMING
/* MAFFT_PROFCHECK=1: recompute with the original functions and compare bitwise */
static void profile_check( char **seq, int clus, double *eff, int len, double **cpmx, double *ogcp, double *fgcp, double *gapfreq, char *sgappat, char *egappat )
{
	static TLS double **c2 = NULL, *o2 = NULL, *f2 = NULL, *g2 = NULL;
	static TLS int cap = 0;
	int k;
	if( !getenv( "MAFFT_PROFCHECK" ) ) return;
	if( len + 2 > cap )
	{
		if( c2 ) FreeDoubleMtx( c2 );
		free( o2 ); free( f2 ); free( g2 );
		cap = len + 2;
		c2 = AllocateDoubleMtx( nalphabets, cap );
		o2 = calloc( cap, sizeof( double ) ); f2 = calloc( cap, sizeof( double ) ); g2 = calloc( cap, sizeof( double ) );
	}
	cpmx_calc_new( seq, c2, eff, len, clus );
	if( sgappat )
	{
		new_OpeningGapCount( o2, clus, seq, eff, len, sgappat );
		new_FinalGapCount( f2, clus, seq, eff, len, egappat );
	}
	else
	{
		st_OpeningGapCount( o2, clus, seq, eff, len );
		st_FinalGapCount( f2, clus, seq, eff, len );
	}
	gapcountf( g2, seq, clus, eff, len );
	for( k=0; k<nalphabets; k++ ) if( memcmp( c2[k], cpmx[k], sizeof( double ) * len ) ) { reporterr( "\nPROFCHECK MISMATCH cpmx %d\n", k ); exit( 1 ); }
	if( memcmp( o2, ogcp, sizeof( double ) * ( len + ( sgappat == NULL ) ) ) ) { reporterr( "\nPROFCHECK MISMATCH ogcp\n" ); exit( 1 ); }
	if( memcmp( f2, fgcp, sizeof( double ) * len ) ) { reporterr( "\nPROFCHECK MISMATCH fgcp\n" ); exit( 1 ); }
	if( memcmp( g2, gapfreq, sizeof( double ) * len ) ) { reporterr( "\nPROFCHECK MISMATCH gapfreq\n" ); exit( 1 ); }
	pt_acc[28] += 0.0;
}
#else
#define profile_check( a, b, c, d, e, f, g, h, i, j )
#endif

/* work/prog2: gapireru( res, ori, gt ) with res == ori.  gapireru writes, for p = 0.. , a gap where
   gt[p] == '-' and otherwise the next character of ori.  Going backwards from the end, the
   character that lands at p is ori[q] with q <= p (q counts the non-gap entries before p), and
   every ori[q'] with q' > q has already been moved, so the row can be expanded in place.  The
   result is the same string gapireru would have produced in a separate buffer. */
static void gapireru_inplace( char *s, int lori, const char *gt, int lnew )
{
	const char gapchar = *newgapstr;
	int p, q = lori - 1;
	s[lnew] = 0;
	for( p=lnew-1; p>=0; p-- )
		s[p] = ( gt[p] == '-' ) ? gapchar : s[q--];
}

typedef struct { char **seq; const char *gt; int lori, lnew; } gapexpctx_t;

static void gapexp_rows( void *vctx, int lo, int hi )
{
	gapexpctx_t *c = (gapexpctx_t *)vctx;
	int i;
	for( i=lo; i<hi; i++ ) gapireru_inplace( c->seq[i], c->lori, c->gt, c->lnew );
}

static void gapexp_all( char **seq, int n, const char *gt, int lori, int lnew )
{
	gapexpctx_t c;
	c.seq = seq; c.gt = gt; c.lori = lori; c.lnew = lnew;
#ifdef enablemultithread
	dphelp_for( n, MAX( 16, n / 32 ), (double)n * lnew, 4.0e6, gapexp_rows, &c );
#else
	gapexp_rows( &c, 0, n );
#endif
}

/* work/prog2: per-thread arenas for A__align's per-call buffers.
   (1) The traceback matrix.  The original used the shared thread-local matrix commonIP, sized to
       (1.3*lgth+110) x (1.3*lgth+110) ints with every row calloc'ed separately, and freed and
       reallocated it (and zero-filled every row again) whenever a larger alignment came along.
       Here the rows are packed with a stride of lgth2+2 ints into one block that is only replaced
       when a call needs more than it holds; the block is sized for the largest alignment that fits
       in alloclen (lgth1+lgth2 < alloclen), so that happens only when alloclen itself grows.
       A__align writes every cell that it or Atracking reads (the DP fills rows 1..lasti-1,
       columns 1..lastj-1; Atracking sets row 0, column 0 and [lgth1][lgth2]), and the original
       reused commonIP across calls without clearing it, so stale contents never matter.
   (2) mseq, the rows that receive the gapped result before it is copied back: rows of
       lgth1+lgth2+100 chars in one block.  Every row that is read was written with a terminator
       first, except in the reuseprofiles case, where mseq1[0] may be read by strlen() unwritten;
       the first byte of every row is set to 0 so such a read still sees "" as with calloc(). */
static TLS int *aijp_blk = NULL;
static TLS size_t aijp_cap = 0;
static TLS int **aijp_rows = NULL;
static TLS int aijp_rowcap = 0;
static TLS char *amseq_blk = NULL;
static TLS size_t amseq_cap = 0;
static TLS char **amseq_rows = NULL;
static TLS int amseq_rowcap = 0;

/* Large traceback blocks are shared between threads through a small pool: a thread returns its
   block after a large alignment (A_ijp_release), and the next large alignment, on whichever
   thread, takes a pooled block that is big enough instead of allocating (and page-faulting in)
   a fresh one.  The contents of a block never matter (see above). */
#define AIJP_POOLMAX 16
#define AIJP_POOLMIN ( (size_t)16 << 20 ) // ints; smaller blocks stay with their thread
#ifdef enablemultithread
static pthread_mutex_t aijp_poolmutex = PTHREAD_MUTEX_INITIALIZER;
#endif
static int *aijp_poolblk[AIJP_POOLMAX];
static size_t aijp_poolcap[AIJP_POOLMAX];
static int aijp_pooln = 0;

/* work/largen: with more than 20000 sequences the original aligns one step at a time.  To keep the
   threaded pass's memory close to that, disttbfast() turns the pool off for those passes
   (A__align_nopool) and each thread gives its large blocks back when a step ends (A__align_trim). */
int A__align_nopool = 0;
void A__align_trim( size_t keep )
{
	if( aijp_cap * sizeof( int ) > keep ) { free( aijp_blk ); aijp_blk = NULL; aijp_cap = 0; }
	if( amseq_cap > keep ) { free( amseq_blk ); amseq_blk = NULL; amseq_cap = 0; }
}

static void A_ijp_release( void ) // give a large block back to the pool
{
#ifdef enablemultithread
	if( A__align_nopool ) return; // work/largen
	if( aijp_cap < AIJP_POOLMIN ) return;
	pthread_mutex_lock( &aijp_poolmutex );
	if( aijp_pooln < AIJP_POOLMAX )
	{
		aijp_poolblk[aijp_pooln] = aijp_blk;
		aijp_poolcap[aijp_pooln] = aijp_cap;
		aijp_pooln++;
		aijp_blk = NULL;
		aijp_cap = 0;
	}
	pthread_mutex_unlock( &aijp_poolmutex );
#endif
}

static void A_ijp_poolfree( void )
{
#ifdef enablemultithread
	int k;
	pthread_mutex_lock( &aijp_poolmutex );
	for( k=0; k<aijp_pooln; k++ ) free( aijp_poolblk[k] );
	aijp_pooln = 0;
	pthread_mutex_unlock( &aijp_poolmutex );
#endif
}

/* work/v8s: the arrays the tiled DP's bands write row after row (the two row buffers, m, mp and the
   traceback rows) start on a 64-byte boundary and the band boundaries are multiples of 16 columns
   (Adp_tiled), so two bands never write the same cache line: on CARC every band boundary otherwise
   moved 3-4 lines between cores (often across sockets) per row.  Layout only: the same values at
   the same indices. */
static void *adp_avec( size_t n, size_t size )
{
	void *p = NULL;
	if( posix_memalign( &p, 64, n * size + 64 ) ) ErrorExit( "Cannot allocate a DP vector in A__align.\n" );
	memset( p, 0, n * size + 64 );
	return( p );
}

static int **A_ijp_get( int lgth1, int lgth2, int alloclen )
{
	size_t stride = ( (size_t)lgth2 + 2 + 15 ) & ~(size_t)15, need = ( (size_t)lgth1 + 2 ) * stride; // work/v8s: rows of 16-int multiples (64 bytes, see adp_avec)
	int i;
#ifdef enablemultithread
	if( need > aijp_cap && need >= AIJP_POOLMIN / 4 && aijp_pooln > 0 )
	{
		int k, best = -1;
		pthread_mutex_lock( &aijp_poolmutex );
		for( k=0; k<aijp_pooln; k++ ) // the smallest pooled block that is big enough
			if( aijp_poolcap[k] >= need && ( best < 0 || aijp_poolcap[k] < aijp_poolcap[best] ) ) best = k;
		if( best >= 0 )
		{
			int *b = aijp_poolblk[best];
			size_t c = aijp_poolcap[best];
			if( aijp_blk && aijp_cap >= AIJP_POOLMIN ) { aijp_poolblk[best] = aijp_blk; aijp_poolcap[best] = aijp_cap; }
			else
			{
				free( aijp_blk );
				aijp_poolblk[best] = aijp_poolblk[aijp_pooln-1];
				aijp_poolcap[best] = aijp_poolcap[aijp_pooln-1];
				aijp_pooln--;
			}
			aijp_blk = b;
			aijp_cap = c;
		}
		pthread_mutex_unlock( &aijp_poolmutex );
	}
#endif
	if( need > aijp_cap )
	{
		size_t half = (size_t)( alloclen > 0 ? alloclen : 0 ) / 2 + 4;
		size_t cap = half * half;
		if( cap < need + need / 4 ) cap = need + need / 4;
		free( aijp_blk );
		aijp_blk = NULL;
		if( posix_memalign( (void **)&aijp_blk, 64, cap * sizeof( int ) ) ) aijp_blk = NULL; // work/v8s: 64-byte aligned rows (adp_avec)
		if( aijp_blk == NULL ) ErrorExit( "Cannot allocate the traceback matrix in A__align.\n" );
		aijp_cap = cap;
	}
	if( lgth1 + 2 > aijp_rowcap )
	{
		free( aijp_rows );
		aijp_rowcap = lgth1 + 2 + ( lgth1 + 2 ) / 2;
		aijp_rows = malloc( aijp_rowcap * sizeof( int * ) );
		if( aijp_rows == NULL ) ErrorExit( "Cannot allocate the traceback matrix in A__align.\n" );
	}
	for( i=0; i<lgth1+2; i++ ) aijp_rows[i] = aijp_blk + (size_t)i * stride;
	return( aijp_rows );
}

static size_t A_ijp_capacity( void ) { return( aijp_cap ); } // work/v8s: ints in this thread's traceback block

static char **A_mseq_get( int nrow, int rowlen )
{
	size_t need = (size_t)nrow * rowlen;
	int i;
	if( need > amseq_cap )
	{
		free( amseq_blk );
		amseq_cap = need + need / 2;
		amseq_blk = malloc( amseq_cap );
		if( amseq_blk == NULL ) ErrorExit( "Cannot allocate mseq in A__align.\n" );
	}
	if( nrow + 1 > amseq_rowcap )
	{
		free( amseq_rows );
		amseq_rowcap = nrow + 1 + ( nrow + 1 ) / 2;
		amseq_rows = malloc( amseq_rowcap * sizeof( char * ) );
		if( amseq_rows == NULL ) ErrorExit( "Cannot allocate mseq in A__align.\n" );
	}
	for( i=0; i<nrow; i++ )
	{
		amseq_rows[i] = amseq_blk + (size_t)i * rowlen;
		amseq_rows[i][0] = 0;
	}
	amseq_rows[nrow] = NULL;
	return( amseq_rows );
}

static void A_arenas_free( void )
{
	A_ijp_poolfree();
	free( aijp_blk ); aijp_blk = NULL; aijp_cap = 0;
	free( aijp_rows ); aijp_rows = NULL; aijp_rowcap = 0;
	free( amseq_blk ); amseq_blk = NULL; amseq_cap = 0;
	free( amseq_rows ); amseq_rows = NULL; amseq_rowcap = 0;
}

#if FASTMATCHCALC // work/v8s: match_calc's sparse tables, allocated at their first use in A__align
#define AW_ENSURE() do { if( doublework == NULL ) { doublework = AllocateFloatMtx( aw_rows, nalphabets ); intwork = AllocateIntMtx( aw_rows, nalphabets+1 ); } } while( 0 )
#else
#define AW_ENSURE() do { if( doublework == NULL ) { doublework = AllocateFloatMtx( nalphabets, aw_rows ); intwork = AllocateIntMtx( nalphabets, aw_rows ); } } while( 0 )
#endif
double A__align( double **n_dynamicmtx, int penalty_l, int penalty_ex_l, char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen, int constraint, double *impmatch, char *sgap1, char *sgap2, char *egap1, char *egap2, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp, int firstmem, int calledbyfulltreebase, double ***cpmxchild0, double ***cpmxchild1, double ***cpmxresult, double orieff1, double orieff2 )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{

	int reuseprofiles;
	static TLS int previousfirstlen; // 2016/Feb/1 // MEMBER NO CHECK GA HITSUYOU!!!!
	static TLS int previousicyc; // 2016/Feb/1 // MEMBER NO CHECK GA HITSUYOU!!!!
	static TLS int previousfirstmem;
	static TLS int previouscall;
	int ngap1, ngap2;
//	int k;
	register int i, j;
	int lasti, lastj;      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lgth1, lgth2;
	int resultlen;
	double wm = 0.0;   /* int ?????? */
	double wmo = 0.0;
	double g;
	double *currentw, *previousw;
//	double fpenalty = (double)penalty;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penalty_ex_l;
#endif
#if 1
	double *wtmp;
	int *ijppt;
	double *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	static TLS double mi, *m;
	static TLS int **ijp;
	static TLS int mpi, *mp;
	static TLS double *w1, *w2;
	static TLS double *match;
	static TLS double *initverticalw;    /* kufuu sureba iranai */
	static TLS double *lastverticalw;    /* kufuu sureba iranai */
	char **mseq1;
	char **mseq2;
	char **mseq;
	static TLS double *ogcp1, *ogcp1o;
	static TLS double *ogcp2, *ogcp2o;
	static TLS double *fgcp1, *fgcp1o;
	static TLS double *fgcp2, *fgcp2o;
	double *ogcp1opt, *ogcp2opt, *fgcp1opt, *fgcp2opt;
	static TLS double **cpmx1;
	double ***cpmx1pt = NULL;
	static TLS double **cpmx2;
	double ***cpmx2pt = NULL;
	static TLS int **intwork;
	static TLS double **doublework;
	static TLS int aw_rows = 0; // work/v8s: rows of doublework / intwork when they are allocated (AW_ENSURE)
	static TLS int orlgth1 = 0, orlgth2 = 0;
	static TLS double *gapfreq1;
	double *gapfreq1pt;
	static TLS double *gapfreq2;
	double *gapfreq2pt;
	static TLS double *dpcolw = NULL; // fast DP path (work/prog): column-candidate scores of the current row
	static TLS int *dpcolij = NULL;
	static TLS double *dpcolg2 = NULL;
	static TLS int dpcolalloc = 0;
	double fpenalty = (double)penalty_l;
	double fpenalty_shift = (double)penalty_shift;
	double *fgcp2pt;
	double *ogcp2pt;
	double fgcp1va;
	double ogcp1va;
	double *gf2pt;
	double *gf2ptpre;
	double gf1va;
	double gf1vapre;
	double headgapfreq1;
	double headgapfreq2;

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
	char *gt1, *gt2, *gt1bk, *gt2bk;
	int usefused, rowfused;
	a_lazy_t lz; // work/prog3
	int sparse2_ready = 1; // work/prog2: match_calc's sparse tables of cpmx2 are built
	int inplace = 0;       // work/prog2: the result rows were written in place (no mseq)
	double fsel_s[4];
	double *fsel_c[4];
	double fsel_su[DENSE_MAXU]; // work/prog3: the row scores of every letter of the profile
	unsigned int fsel_mask = 0;
	int fsel_n = 0;
	double *gen_cu[DENSE_MAXU];
	int gen_ready = 0;
	PT_DECL;

//	for( i=0; i<icyc; i++ ) fprintf( stderr, "%s, %f\n", seq1[i], eff1[i] );
//	for( i=0; i<jcyc; i++ ) fprintf( stderr, "%s, %f\n", seq2[i], eff2[i] );


	if( seq1 == NULL )
	{
		if( orlgth1 )
		{
//			fprintf( stderr, "## Freeing local arrays in A__align\n" );
			orlgth1 = 0;
			orlgth2 = 0;

			imp_match_init_strict( NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL, 0, 0 );

			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );


			FreeFloatVec( ogcp1 );
			FreeFloatVec( ogcp1o );
			FreeFloatVec( ogcp2 );
			FreeFloatVec( ogcp2o );
			FreeFloatVec( fgcp1 );
			FreeFloatVec( fgcp1o );
			FreeFloatVec( fgcp2 );
			FreeFloatVec( fgcp2o );


			FreeFloatMtx( cpmx1 );
			FreeFloatMtx( cpmx2 );

			FreeFloatVec( gapfreq1 );
			FreeFloatVec( gapfreq2 );

			if( doublework ) FreeFloatMtx( doublework ); // work/v8s (AW_ENSURE)
			if( intwork ) FreeIntMtx( intwork );
			doublework = NULL; intwork = NULL;
		}
		else
		{
//			fprintf( stderr, "## Not allocated\n" );
		}
		if( dpcolalloc )
		{
			free( dpcolw ); dpcolw = NULL;
			free( dpcolij ); dpcolij = NULL;
			free( dpcolg2 ); dpcolg2 = NULL;
			dpcolalloc = 0;
		}
		A_arenas_free();
		return( 0.0 );
	}


	PT_START();
	/* work/prog3: lazy gap insertion (lazygap.h); a stale side's rows are neither read nor written */
	lz = A_lz;
	A_lz.on = 0;
	if( !lz.on ) { lz.stale[0] = lz.stale[1] = lz.defer[0] = lz.defer[1] = 0; }
	lgth1 = lz.stale[0] ? lz.len[0] : strlen( seq1[0] );
	lgth2 = lz.stale[1] ? lz.len[1] : strlen( seq2[0] );
	if( lz.on && ( lgth1 == 0 || lgth2 == 0 || trywarp || constraint || !cpmxresult || specificityconsideration != 0.0
	           || ( lz.stale[0] && !( cpmxchild0 && *cpmxchild0 ) ) || ( lz.stale[1] && !( cpmxchild1 && *cpmxchild1 ) ) ) )
	{
		reporterr( "\nA__align: unsupported lazy-gap request (lgth %d %d)\n", lgth1, lgth2 ); // the caller checks all of these
		exit( 1 );
	}
	if( lz.on ) { lz.defer[0] |= lz.stale[0]; lz.defer[1] |= lz.stale[1]; }
#if 0
	if( lgth1 == 0 || lgth2 == 0 )
	{
		fprintf( stderr, "WARNING (Aalignmm): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif
	if( lgth1 == 0 && lgth2 == 0 )
		return( 0.0 );

	if( lgth1 == 0 )
	{
		for( i=0; i<icyc; i++ )
		{
			j = lgth2;
			seq1[i][j] = 0;
			while( j ) seq1[i][--j] = *newgapstr;
//			fprintf( stderr, "seq1[i] = %s\n", seq1[i] );
		}
		return( 0.0 );
	}

	if( lgth2 == 0 )
	{
		for( i=0; i<jcyc; i++ )
		{
			j = lgth1;
			seq2[i][j] = 0;
			while( j ) seq2[i][--j] = *newgapstr;
//			fprintf( stderr, "seq2[i] = %s\n", seq2[i] );
		}
		return( 0.0 );
	}

	warpbase = lgth1 + lgth2;
	warpis = NULL;
	warpjs = NULL;
	warpn = 0;



	if( trywarp )
	{
//		fprintf( stderr, "IN A__align, penalty_shift = %d\n", penalty_shift );
		if( headgp == 0 || tailgp == 0 )
		{
			fprintf( stderr, "At present, headgp and tailgp must be 1 to allow shift.\n" );
			exit( 1 );
		}
		wmrecords = AllocateFloatVec( lgth2+1 );
		warpi = AllocateIntVec( lgth2+1 );
		warpj = AllocateIntVec( lgth2+1 );
		prevwmrecords = AllocateFloatVec( lgth2+1 );
		prevwarpi = AllocateIntVec( lgth2+1 );
		prevwarpj = AllocateIntVec( lgth2+1 );
		for( i=0; i<lgth2+1; i++ ) wmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwarpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) prevwarpj[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpj[i] = -warpbase;
	}


#if 0
	fprintf( stderr, "####  eff in SA+++align\n" );
	fprintf( stderr, "####  seq1[0] = %s\n", seq1[0] );
	fprintf( stderr, "####  strlen( seq1[0] ) = %d\n", strlen( seq1[0] ) );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "eff1[%d] = %f\n", i, eff1[i] );
	fprintf( stderr, "####  seq2[0] = %s\n", seq2[0] );
	fprintf( stderr, "####  strlen( seq2[0] ) = %d\n", strlen( seq2[0] ) );
	for( i=0; i<jcyc; i++ ) fprintf( stderr, "eff2[%d] = %f\n", i, eff2[i] );
#endif

	if( constraint ) // work/prog2: per-thread arena (was AllocateCharMtx per call); without a constraint the gaps are inserted in place
	{
		mseq = A_mseq_get( icyc+jcyc, lgth1+lgth2+100 );
		mseq1 = mseq;
		mseq2 = mseq + icyc;
	}
	else
		mseq = mseq1 = mseq2 = NULL;
#ifdef PROGTIMING
	{ double _n = pt_now(); pt_acc[24] += _n - _pt_t; }
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
			FreeIntVec( mp );

			FreeFloatVec( ogcp1 );
			FreeFloatVec( ogcp1o );
			FreeFloatVec( ogcp2 );
			FreeFloatVec( ogcp2o );
			FreeFloatVec( fgcp1 );
			FreeFloatVec( fgcp1o );
			FreeFloatVec( fgcp2 );
			FreeFloatVec( fgcp2o );


			FreeFloatMtx( cpmx1 );
			FreeFloatMtx( cpmx2 );

			FreeFloatVec( gapfreq1 );
			FreeFloatVec( gapfreq2 );

			if( doublework ) FreeFloatMtx( doublework ); // work/v8s: may never have been allocated (AW_ENSURE)
			if( intwork ) FreeIntMtx( intwork );
		}

		ll1 = MAX( (int)(1.3*lgth1), orlgth1 ) + 100;
		ll2 = MAX( (int)(1.3*lgth2), orlgth2 ) + 100;

#if DEBUG
		fprintf( stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2 );
#endif

		w1 = adp_avec( ll2+2, sizeof( double ) ); // work/v8s: 64-byte aligned (zeroed, as calloc; freed with free)
		w2 = adp_avec( ll2+2, sizeof( double ) );
		match = AllocateFloatVec( ll2+2 );

		initverticalw = AllocateFloatVec( ll1+2 );
		lastverticalw = AllocateFloatVec( ll1+2 );

		m = adp_avec( ll2+2, sizeof( double ) ); // work/v8s: see w1
		mp = adp_avec( ll2+2, sizeof( int ) );


		ogcp1 = AllocateFloatVec( ll1+2 );
		ogcp1o = AllocateFloatVec( ll1+2 );
		ogcp2 = AllocateFloatVec( ll2+2 );
		ogcp2o = AllocateFloatVec( ll2+2 );
		fgcp1 = AllocateFloatVec( ll1+2 );
		fgcp1o = AllocateFloatVec( ll1+2 );
		fgcp2 = AllocateFloatVec( ll2+2 );
		fgcp2o = AllocateFloatVec( ll2+2 );

		cpmx1 = AllocateFloatMtx( nalphabets, ll1+2 );
		cpmx2 = AllocateFloatMtx( nalphabets, ll2+2 );
		previousfirstlen = -1;
		previousicyc = -1;

		gapfreq1 = AllocateFloatVec( ll1+2 );
		gapfreq2 = AllocateFloatVec( ll2+2 );

		/* work/v8s: match_calc's sparse tables (~2 x MAX( ll1, ll2 ) small rows, ~50,000 callocs for the
		   steps near the root, 5-8 ms whenever a large step lands on a thread whose arrays are smaller)
		   are only needed off the dense path; they are allocated at their first use (AW_ENSURE), with
		   the same size and the same (zeroed) content as here before */
		aw_rows = MAX( ll1, ll2 )+2;
		doublework = NULL;
		intwork = NULL;

#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;
	}


	for( i=0; i<icyc; i++ )
	{
		if( mseq ) mseq1[i] = mseq[i];
		if( !lz.stale[0] ) seq1[i][lgth1] = 0;
	}
	for( j=0; j<jcyc; j++ )
	{
		if( mseq ) mseq2[j] = mseq[icyc+j];
		if( !lz.stale[1] ) seq2[j][lgth2] = 0;
	}


	ijp = A_ijp_get( lgth1, lgth2, alloclen ); // work/prog2: per-thread arena (was the shared commonIP matrix)
#ifdef PROGTIMING
	double _pt_s1 = pt_now();
#endif

#if 0
	{
		double t = 0.0;
		for( i=0; i<icyc; i++ )
			t += eff1[i];
	fprintf( stderr, "## totaleff = %f\n", t );
	}
#endif

#if 0
	reporterr( "calledbyfulltreebase=%d\n", calledbyfulltreebase );
	reporterr( "previouscall=%d\n", previouscall );
	reporterr( "firstmem=%d\n", firstmem );
	reporterr( "previousfirstmem=%d\n", previousfirstmem );
	reporterr( "lgth1=%d\n", lgth1 );
	reporterr( "previousfirstlen=%d\n", previousfirstlen );
	reporterr( "icyc=%d\n", icyc );
	reporterr( "previousicyc=%d\n", previousicyc );
#endif
//	if( calledbyfulltreebase && previouscall && firstmem >= 0 && firstmem == previousfirstmem && lgth1 == previousfirstlen && nthread == 0 && icyc == previousicyc+1 ) // mouhitotsu jouken ga iru. firstmem >= 0 && firstmem == previousfirstmem && calledbytreebase && prevcalledbytreebase
	if( calledbyfulltreebase == 1 && previouscall && firstmem >= 0 && firstmem == previousfirstmem && lgth1 == previousfirstlen && nthread == 0 && icyc == previousicyc+1 ) // mouhitotsu jouken ga iru. firstmem >= 0 && firstmem == previousfirstmem && calledbytreebase && prevcalledbytreebase
		reuseprofiles = 1;
	else
		reuseprofiles = 0;
// Atracking de gapireru() wo skip surutameni reuseprofile flag wo tsukau.

#if 0
	reporterr( "nthread=%d\n", nthread );
	reporterr( "reuseprofiles=%d\n", reuseprofiles );
	reporterr( "calledbyfulltreebase=%d\n", calledbyfulltreebase );
	reporterr( "previouscall=%d\n", previouscall );
	reporterr( "firstmem=%d, prefiousfirstmem=%d\n", firstmem, previousfirstmem );
	reporterr( "lgth1=%d, previousfirstlen=%d\n", lgth1, previousfirstlen );
	if( reuseprofiles ) exit(1 );
#endif

	if( n_dis[0][amino_n['-']] != 0 )
	{
		reporterr( "Bug probably in versions >7.36.  Please report this issue to katoh@ifrec.osaka-u.ac.jp\n" );
		exit( 1 );
	}


	if( cpmxresult && specificityconsideration == 0.0 ) // n_dynamicmtx ga henka suru toki profile ha sairiyou dekinai.
	{
		if( sgap1 )
		{
			reporterr( "The combination of sgap1 and cpmxhit is not supported. See Salignmm.c\n" );
			exit( 1 );
		}


		if( cpmxchild0 && *cpmxchild0 )
		{
//			reporterr( "\nUse cpmxhist for child 0!\n" );
			cpmx1pt = (cpmxchild0);
			gapfreq1pt = (*cpmxchild0)[nalphabets];
			ogcp1opt = (*cpmxchild0)[nalphabets+1];
			fgcp1opt = (*cpmxchild0)[nalphabets+2];
		}
		else
		{
//			reporterr( "\nDo not use cpmxhist for child 0!\n" );
			cpmx1pt = &cpmx1;
			gapfreq1pt = gapfreq1;
			ogcp1opt = ogcp1o;
			fgcp1opt = fgcp1o;
			/* work/prog2: cpmx_calc_new + gapcountf + st_OpeningGapCount + st_FinalGapCount in one pass */
			profile_fused( seq1, icyc, eff1, lgth1, *cpmx1pt, ogcp1opt, fgcp1opt, gapfreq1pt, NULL );
			profile_check( seq1, icyc, eff1, lgth1, *cpmx1pt, ogcp1opt, fgcp1opt, gapfreq1pt, NULL, NULL );
			for( i=0; i<lgth1; i++ ) gapfreq1pt[i] = 1.0 - gapfreq1pt[i];
		}
	
		if( cpmxchild1 && *cpmxchild1 )
		{
//			reporterr( "\nUse cpmxhist for child 1!\n" );
			cpmx2pt = (cpmxchild1);
			gapfreq2pt = (*cpmxchild1)[nalphabets];
			ogcp2opt = (*cpmxchild1)[nalphabets+1];
			fgcp2opt = (*cpmxchild1)[nalphabets+2];
		}
		else
		{
//			reporterr( "\nDo not use cpmxhist for child 1!\n" );
			cpmx2pt = &cpmx2;
			gapfreq2pt = gapfreq2;
			ogcp2opt = ogcp2o;
			fgcp2opt = fgcp2o;
			profile_fused( seq2, jcyc, eff2, lgth2, *cpmx2pt, ogcp2opt, fgcp2opt, gapfreq2pt, NULL ); // work/prog2
			profile_check( seq2, jcyc, eff2, lgth2, *cpmx2pt, ogcp2opt, fgcp2opt, gapfreq2pt, NULL, NULL );
			for( i=0; i<lgth2; i++ ) gapfreq2pt[i] = 1.0 - gapfreq2pt[i];
		}

		//sgap1 == 1 ni taiou shite iani node, legacygap niyorazu, headgap to tailgap ha nai to minasu.

		headgapfreq1 = 1.0;
		headgapfreq2 = 1.0;
		gapfreq1pt[lgth1] = 1.0;
		gapfreq2pt[lgth2] = 1.0;
	}
	else
	{
		cpmx1pt = &cpmx1;
		cpmx2pt = &cpmx2;
		gapfreq1pt = gapfreq1;
		gapfreq2pt = gapfreq2;
		ogcp1opt = ogcp1o;
		ogcp2opt = ogcp2o;
		fgcp1opt = fgcp1o;
		fgcp2opt = fgcp2o;

		if( !reuseprofiles && legacygapcost == 0 )
		{
			/* work/prog2: cpmx_calc_new, the opening/final gap counts and gapcountf of each group in
			   one pass (profile_fused); the remaining steps exactly as below */
			profile_fused( seq1, icyc, eff1, lgth1, *cpmx1pt, ogcp1opt, fgcp1opt, gapfreq1pt, sgap1 ? sgap1 : NULL );
			profile_check( seq1, icyc, eff1, lgth1, *cpmx1pt, ogcp1opt, fgcp1opt, gapfreq1pt, sgap1, egap1 );
			profile_fused( seq2, jcyc, eff2, lgth2, *cpmx2pt, ogcp2opt, fgcp2opt, gapfreq2pt, sgap1 ? sgap2 : NULL );
			profile_check( seq2, jcyc, eff2, lgth2, *cpmx2pt, ogcp2opt, fgcp2opt, gapfreq2pt, sgap1 ? sgap2 : NULL, egap2 );
			if( sgap1 )
			{
				outgapcount( &headgapfreq1, icyc, sgap1, eff1 );
				outgapcount( &headgapfreq2, jcyc, sgap2, eff2 );
				outgapcount( gapfreq1pt+lgth1, icyc, egap1, eff1 );
				outgapcount( gapfreq2pt+lgth2, jcyc, egap2, eff2 );
			}
			else
			{
				headgapfreq1 = 0.0;
				headgapfreq2 = 0.0;
				gapfreq1pt[lgth1] = 0.0;
				gapfreq2pt[lgth2] = 0.0;
			}
			for( i=0; i<lgth1+1; i++ ) gapfreq1pt[i] = 1.0 - gapfreq1pt[i];
			for( i=0; i<lgth2+1; i++ ) gapfreq2pt[i] = 1.0 - gapfreq2pt[i];
			headgapfreq1 = 1.0 - headgapfreq1;
			headgapfreq2 = 1.0 - headgapfreq2;
		}
		else
		{
		if( reuseprofiles )
		{
//			reporterr( "reuse" );
			cpmx_calc_add( seq1, *cpmx1pt, eff1, lgth1, icyc );
		}
		else
		{
//			reporterr( "new profile" );
			cpmx_calc_new( seq1, *cpmx1pt, eff1, lgth1, icyc );
		}
		cpmx_calc_new( seq2, *cpmx2pt, eff2, lgth2, jcyc );

		if( sgap1 )
		{
			new_OpeningGapCount( ogcp1opt, icyc, seq1, eff1, lgth1, sgap1 );
			new_FinalGapCount( fgcp1opt, icyc, seq1, eff1, lgth1, egap1 );
	
			new_OpeningGapCount( ogcp2opt, jcyc, seq2, eff2, lgth2, sgap2 );
			new_FinalGapCount( fgcp2opt, jcyc, seq2, eff2, lgth2, egap2 );
	
			outgapcount( &headgapfreq1, icyc, sgap1, eff1 );
			outgapcount( &headgapfreq2, jcyc, sgap2, eff2 );
			outgapcount( gapfreq1pt+lgth1, icyc, egap1, eff1 );
			outgapcount( gapfreq2pt+lgth2, jcyc, egap2, eff2 );
		}
		else
		{
			if( reuseprofiles )
			{
				st_OpeningGapAdd( ogcp1opt, icyc, seq1, eff1, lgth1 );
				st_FinalGapAdd( fgcp1opt, icyc, seq1, eff1, lgth1 );
			}
			else
			{
				st_OpeningGapCount( ogcp1opt, icyc, seq1, eff1, lgth1 );
				st_FinalGapCount( fgcp1opt, icyc, seq1, eff1, lgth1 );
			}

			st_OpeningGapCount( ogcp2opt, jcyc, seq2, eff2, lgth2 );
			st_FinalGapCount( fgcp2opt, jcyc, seq2, eff2, lgth2 );

			headgapfreq1 = 0.0;
			headgapfreq2 = 0.0;
			gapfreq1pt[lgth1] = 0.0;
			gapfreq2pt[lgth2] = 0.0;
		}

		if( legacygapcost == 0 ) // sgap1 no toki headgapfreq, gapfreq1pt[lgth] ha 1 toha kagiranai.
		{
			if( reuseprofiles )
				gapcountadd( gapfreq1pt, seq1, icyc, eff1, lgth1 );
			else
				gapcountf( gapfreq1pt, seq1, icyc, eff1, lgth1 );
			gapcountf( gapfreq2pt, seq2, jcyc, eff2, lgth2 );
			for( i=0; i<lgth1+1; i++ ) gapfreq1pt[i] = 1.0 - gapfreq1pt[i];
			for( i=0; i<lgth2+1; i++ ) gapfreq2pt[i] = 1.0 - gapfreq2pt[i];
			headgapfreq1 = 1.0 - headgapfreq1;
			headgapfreq2 = 1.0 - headgapfreq2;
		}
		else
		{
			for( i=0; i<lgth1+1; i++ ) gapfreq1pt[i] = 1.0;
			for( i=0; i<lgth2+1; i++ ) gapfreq2pt[i] = 1.0;
			headgapfreq1 = 1.0;
			headgapfreq2 = 1.0;
		}
		}
	}


#if 0
	reporterr( "\n" );
	reporterr( "seq1 = \n" );
	for( i=0; i<icyc; i++ )
		reporterr( "%s\n", seq1[i] );

	reporterr( "seq2 = \n" );
	for( i=0; i<jcyc; i++ )
		reporterr( "%s\n", seq2[i] );

	reporterr( "\n    " );
	for( j=0; j<nalphabets; j++ ) 
		reporterr( "%4c ", amino[j] );
	reporterr( "\ncpmx1pt = \n" );
	for( i=0; i<lgth1; i++ ) 
	{
		reporterr( "%3d ", i );
		for( j=0; j<nalphabets; j++ ) 
			reporterr( "%4.2f ", (*cpmx1pt)[j][i] );
		reporterr( "\n" );
	}
#endif

#if 0
	reporterr( "\n    " );
	for( j=0; j<nalphabets; j++ ) 
		reporterr( "%4c ", amino[j] );
	reporterr( "\ncpmx2pt = \n" );
	for( i=0; i<lgth2; i++ ) 
	{
		reporterr( "%3d ", i );
		for( j=0; j<nalphabets; j++ ) 
			reporterr( "%4.2f ", (*cpmx2pt)[j][i] );
		reporterr( "\n" );
	}
	reporterr( "\ngapfreq1 = " );
	for( i=0; i<lgth1+1; i++ ) reporterr( "%4.2f ", gapfreq1pt[i] );
	reporterr( "\ngapfreq2 = " );
	for( i=0; i<lgth2+1; i++ ) reporterr( "%4.2f ", gapfreq2pt[i] );
	reporterr( "\n" );

	reporterr( "\nogcp1opt = " );
	for( i=0; i<lgth1+0; i++ ) reporterr( "%4.2f ", ogcp1opt[i] );
	reporterr( "\nfgcp1opt = " );
	for( i=0; i<lgth1+0; i++ ) reporterr( "%4.2f ", fgcp1opt[i] );
	reporterr( "\n" );

	reporterr( "\nogcp2opt = " );
	for( i=0; i<lgth2+0; i++ ) reporterr( "%4.2f ", ogcp2opt[i] );
	reporterr( "\nfgcp2opt = " );
	for( i=0; i<lgth2+0; i++ ) reporterr( "%4.2f ", fgcp2opt[i] );
	reporterr( "\n" );
#endif


#if 0
	fprintf( stderr, "\ngapfreq1[] =" );
	for( i=0; i<lgth1; i++ ) fprintf( stderr, "%5.2f ", gapfreq1[i] );
	fprintf( stderr, "\n" );

	fprintf( stderr, "\ngapfreq2[] =" );
	for( i=0; i<lgth2; i++ ) fprintf( stderr, "%5.2f ", gapfreq2[i] );
	fprintf( stderr, "\n" );
#endif
	

#ifdef PROGTIMING
	double _pt_s2 = pt_now(); pt_acc[25] += _pt_s1 - _pt_t; pt_acc[26] += _pt_s2 - _pt_s1; // (25 includes 24)
#endif
	for( i=0; i<lgth1; i++ ) 
	{
		ogcp1[i] = 0.5 * ( 1.0 - ogcp1opt[i] ) * fpenalty * ( gapfreq1pt[i] );
		fgcp1[i] = 0.5 * ( 1.0 - fgcp1opt[i] ) * fpenalty * ( gapfreq1pt[i] );
	}

	for( i=0; i<lgth2; i++ ) 
	{
		ogcp2[i] = 0.5 * ( 1.0 - ogcp2opt[i] ) * fpenalty * ( gapfreq2pt[i] );
		fgcp2[i] = 0.5 * ( 1.0 - fgcp2opt[i] ) * fpenalty * ( gapfreq2pt[i] );
	}
#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

	currentw = w1;
	previousw = w2;

	/* work/prog2: with the dense letter lists (at most 8 letters on either side), the first column
	   and the first row are formed densely (bit-identical, see dense_initvertical and
	   match_calc_dense), and the sparse tables of match_calc are built only if a row later needs
	   them (sparse2_ready). */
	if( FASTDP_ON && !constraint ) match_calc_dense_prep( *cpmx2pt, lgth2 );
	else dense_nu = -1;
	if( dense_nu >= 0 ) dense1_prep( *cpmx1pt, lgth1+1 ); // columns 0..lgth1: the first column and every DP row
	if( dense_nu >= 0 && dense1_n <= DENSE_MAXU )
	{
		dense_initvertical( n_dynamicmtx, initverticalw, *cpmx2pt, *cpmx1pt, lgth1 );
		match_calc_dense( n_dynamicmtx, currentw, *cpmx1pt, *cpmx2pt, 0, lgth2 );
		sparse2_ready = 0;
#ifdef PROGTIMING
		if( getenv( "MAFFT_PROFCHECK" ) ) // compare with the original sparse computation, bitwise
		{
			double *a1 = malloc( sizeof( double ) * ( lgth1 + 1 ) ), *a2 = malloc( sizeof( double ) * ( lgth2 + 1 ) );
			AW_ENSURE();
			match_calc( n_dynamicmtx, a1, *cpmx2pt, *cpmx1pt, 0, lgth1, doublework, intwork, 1 );
			match_calc( n_dynamicmtx, a2, *cpmx1pt, *cpmx2pt, 0, lgth2, doublework, intwork, 1 );
			sparse2_ready = 1;
			if( memcmp( a1, initverticalw, sizeof( double ) * lgth1 ) || memcmp( a2, currentw, sizeof( double ) * lgth2 ) )
			{
				reporterr( "\nPROFCHECK MISMATCH in the first column/row (lgth1=%d lgth2=%d)\n", lgth1, lgth2 );
				exit( 1 );
			}
			free( a1 ); free( a2 );
		}
#endif
	}
	else
	{
	AW_ENSURE();
	match_calc( n_dynamicmtx, initverticalw, *cpmx2pt, *cpmx1pt, 0, lgth1, doublework, intwork, 1 );
	if( constraint )
		imp_match_out_vead_tate( initverticalw, 0, lgth1 ); // 060306

	match_calc( n_dynamicmtx, currentw, *cpmx1pt, *cpmx2pt, 0, lgth2, doublework, intwork, 1 );
	if( constraint )
		imp_match_out_vead( currentw, 0, lgth2 ); // 060306
	sparse2_ready = 1;
	if( dense_nu >= 0 ) dense1_prep( *cpmx1pt, lgth1+1 );
	}
#if 0 // -> tbfast.c // impossible
	if( localhom )
		imp_match_calc( n_dynamicmtx, currentw, icyc, jcyc, lgth1, lgth2, seq1, seq2, eff1, eff2, localhom, 1, 0 );

#endif
	if( headgp == 1 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
//			initverticalw[i] += ( ogcp1[0] + fgcp1[i-1] ) ;
			initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2pt[0] ) ;
#if USE_PENALTY_EX // 2018/Apr/23
			initverticalw[i] += fpenalty_ex * i;
//			reporterr( "added _ex\n" );
#endif
		}
		for( j=1; j<lgth2+1; j++ )
		{
//			currentw[j] += ( ogcp2[0] + fgcp2[j-1] ) ;
			currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1pt[0] ) ;
#if USE_PENALTY_EX // 2018/Apr/23
			currentw[j] += fpenalty_ex * j;
//			reporterr( "added _ex\n" );
#endif
		}
	}
	else
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += fpenalty * TERMGAPFAC; // motto tanjun de yoi?
#if USE_PENALTY_EX // 2018/Apr/22
			initverticalw[i] += fpenalty_ex * i * TERMGAPFAC_EX;
//			reporterr( "added _ex\n" );
#endif
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += fpenalty * TERMGAPFAC; // motto tanjun de yoi?
#if USE_PENALTY_EX // 2018/Apr/22
			currentw[j] += fpenalty_ex * j * TERMGAPFAC_EX;
//			reporterr( "added _ex\n" );
#endif
		}
	}


#if OUTGAP0TRY
	else
	{
		fprintf( stderr, "offset = %d\n", offset );
		for( j=1; j<lgth2+1; j++ )
			currentw[j] -= offset * j / 2.0;
		for( i=1; i<lgth1+1; i++ )
			initverticalw[i] -= offset * i / 2.0;
	}
#endif


#if 0
	fprintf( stderr, "\n   " );
	for( j=0; j<lgth2+1; j++ ) fprintf( stderr, "    %c ", seq2[0][j] );
	fprintf( stderr, "\n%c ", seq1[0][0]  );
	for( j=0; j<lgth2+1; j++ )
	{
		fprintf( stderr, "%5.0f ", currentw[j] );
	}
	fprintf( stderr, "\n"  );
#endif


	for( j=1; j<lgth2+1; ++j ) 
	{
//		m[j] = currentw[j-1] + ogcp1[1]; mp[j] = 0;
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2pt[j-1]; mp[j] = 0;;
	}
	if( lgth2 == 0 )
		lastverticalw[0] = 0.0; // Falign kara yobaretatoki kounarukanousei ari
	else
		lastverticalw[0] = currentw[lgth2-1];

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
fprintf( stderr, "fcgp\n" );
for( i=0; i<lgth1; i++ ) 
	fprintf( stderr, "fgcp1[%d]=%f\n", i, ogcp1[i] );
for( i=0; i<lgth2; i++ ) 
	fprintf( stderr, "fgcp2[%d]=%f\n", i, ogcp2[i] );
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
	if( pt_dpcheck() )
	{
		chk_m = realloc( chk_m, sizeof( double ) * ( lgth2+2 ) );
		chk_mp = realloc( chk_mp, sizeof( int ) * ( lgth2+2 ) );
		chk_cur = realloc( chk_cur, sizeof( double ) * ( lgth2+2 ) );
		chk_ijp = realloc( chk_ijp, sizeof( int ) * ( lgth2+2 ) );
		chk_match = realloc( chk_match, sizeof( double ) * ( lgth2+2 ) );
	}
#endif

	usefused = ( FASTDP_ON && !trywarp && !constraint && dense_nu >= 0 && fpenalty_ex == 0.0 && lgth2 >= 1 );

	PT_ADD(7);
#ifdef PROGTIMING
	pt_acc[34] += 1.0; // work/v8s: A__align DPs of this step; cells: 32 fused not tiled, 33 tiled, 35 not fused
	if( usefused && chudanpt == NULL && lasti > 2 ) pt_acc[32] += (double)( lasti - 1 ) * lgth2; else pt_acc[35] += (double)( lasti - 1 ) * lgth2;
#endif
#ifdef enablemultithread
	/* work/prog2: large fused DPs are split into tiles that idle helper threads share (Adp_tiled);
	   the loop below is skipped when that path ran.  MAFFT_TILECHECK=1 (timing builds) reruns the
	   serial loop afterwards on the saved initial state and compares everything bitwise. */
	if( usefused && chudanpt == NULL && lasti > 2 )
	{
#ifdef PROGTIMING
		if( pt_tilecheck() ) tilechk_save_init( lgth2, currentw, previousw, m, mp, wm );
		else tilechk_pending = 0;
#endif
		if( Adp_tiled( lgth1, lgth2, lasti, n_dynamicmtx, *cpmx1pt, *cpmx2pt, currentw, previousw, m, mp, ijp,
			gapfreq1pt, gapfreq2pt, fgcp1, ogcp1, fgcp2, ogcp2, initverticalw, lastverticalw, fpenalty_ex, &wm ) )
		{
			if( ( lasti - 1 ) & 1 ) { wtmp = currentw; currentw = previousw; previousw = wtmp; }
#ifdef PROGTIMING
			pt_acc[29] += 1.0; pt_acc[30] += (double)( lasti - 1 ) * lgth2;
			pt_acc[32] -= (double)( lasti - 1 ) * lgth2; pt_acc[33] += (double)( lasti - 1 ) * lgth2;
			if( pt_tilecheck() )
			{
				tilechk_save_result( lgth2, lasti, m, mp, ijp, lastverticalw, wm );
				atr_tiled = 0; // work/v8s: the original loop writes the row-major matrix
				tilechk_restore_init( lgth2, m, mp, &wm ); // the initial state again, then the original loop
				currentw = tilechk_cur0; previousw = tilechk_prev0;
				goto adp_rowloop;
			}
#endif
			goto adp_rowsdone;
		}
	}
#endif
#if defined( PROGTIMING ) && defined( enablemultithread )
	adp_rowloop:
#endif
	for( i=1; i<lasti; i++ )
	{

#ifdef enablemultithread
//		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
		if( chudanpt && *chudanpt != chudanref ) 
		{
//			fprintf( stderr, "\n\n## CHUUDAN!!! S\n" );
			*chudanres = 1;
			// work/prog2: mseq lives in a per-thread arena
			return( -1.0 );
		}
#endif
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

		rowfused = 0;
		if( usefused )
		{
			/* Letters whose row score is exactly 0.0 (e.g. the gap letter, which scores 0 against
			   everything) only add exact zeros to the match sums, so they are dropped; if at most 4
			   letters remain, the fused kernel is used, padded with 0.0 * (weights >= 0) terms,
			   which are exact no-ops as well (see match_calc_dense).  work/prog3: otherwise the
			   general kernel (Adp_row_denseK) takes the row with all its letters (rowfused = 2). */
			int u;
			fsel_n = 0;
			fsel_mask = 0;
			dense_scores1( n_dynamicmtx, *cpmx1pt, i, fsel_su );
			for( u=0; u<dense_nu; u++ )
			{
				if( fsel_su[u] == 0.0 ) continue;
				if( fsel_n < 4 )
				{
					fsel_s[fsel_n] = fsel_su[u];
					fsel_c[fsel_n] = (*cpmx2pt)[dense_u[u]];
				}
				fsel_mask |= 1u << u;
				fsel_n++;
			}
			if( fsel_n <= 4 )
			{
				for( u=fsel_n; u<4; u++ )
				{
					fsel_s[u] = 0.0;
					fsel_c[u] = (*cpmx2pt)[0];
				}
				rowfused = 1;
			}
			else if( adp_general_on() )
			{
				rowfused = 2;
				if( !gen_ready ) // the profile's letter lists (once per call)
				{
					for( u=0; u<dense_nu; u++ ) gen_cu[u] = (*cpmx2pt)[dense_u[u]];
					if( lgth2 / 4 + 2 > gen_bmcap ) { free( gen_bm ); gen_bmcap = lgth2 / 4 + 2 + lgth2 / 8; gen_bm = malloc( sizeof( unsigned int ) * gen_bmcap ); }
					Adp_blockmasks( gen_bm, (const double * const *)gen_cu, dense_nu, lgth2 );
					gen_ready = 1;
				}
			}
			else
				rowfused = 0;
		}

#ifdef PROGTIMING
		if( v8s_stat_on() ) v8s_stat_row( 0, rowfused, rowfused ? fsel_n : 0, lgth2, fsel_mask, rowfused == 2 ? gen_bm : NULL );
		if( !rowfused || pt_dpcheck() ) { double _mc0 = pt_now();
#else
		if( !rowfused ) {
#endif
#ifdef PROGTIMING
		if( dense_nu >= 0 && pt_dpcheck() )
		{
			double t0 = pt_now(), t1, t2;
			match_calc_dense( n_dynamicmtx, chk_match, *cpmx1pt, *cpmx2pt, i, lgth2 );
			t1 = pt_now();
			if( !sparse2_ready ) { AW_ENSURE(); match_calc_sparse_init( *cpmx2pt, lgth2, doublework, intwork ); sparse2_ready = 1; }
			match_calc( n_dynamicmtx, currentw, *cpmx1pt, *cpmx2pt, i, lgth2, doublework, intwork, 0 );
			t2 = pt_now();
			pt_acc[16] += t1 - t0; pt_acc[17] += t2 - t1;
			if( memcmp( chk_match, currentw, sizeof( double ) * lgth2 ) )
			{
				reporterr( "\nDPCHECK MISMATCH in match_calc_dense row %d (lgth2=%d, nu=%d)\n", i, lgth2, dense_nu );
				exit( 1 );
			}
		}
		else
#endif
		if( dense_nu >= 0 )
			match_calc_dense( n_dynamicmtx, currentw, *cpmx1pt, *cpmx2pt, i, lgth2 );
		else
		{
			if( !sparse2_ready ) { AW_ENSURE(); match_calc_sparse_init( *cpmx2pt, lgth2, doublework, intwork ); sparse2_ready = 1; }
			match_calc( n_dynamicmtx, currentw, *cpmx1pt, *cpmx2pt, i, lgth2, doublework, intwork, 0 );
		}
#ifdef PROGTIMING
		pt_acc[11] += pt_now() - _mc0;
#endif
		}
#ifdef PROGTIMING
		pt_acc[12] += (double)lgth2;
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
		if( constraint )
		{
//			fprintf( stderr, "Calling imp_match_calc (o) lgth = %d, i = %d\n", lgth1, i );
#if  0
			imp_match_out_vead( currentw, i, lgth2 );
#else
			imp_match_out_vead( currentw, i, lgth2 );
#endif
		}
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

#if 0
		fprintf( stderr, "%c ", seq1[0][i] );
		for( j=0; j<lgth2+1; j++ )
		{
			fprintf( stderr, "%5.0f ", currentw[j] );
		}
		fprintf( stderr, "\n"  );
#endif
	
//		mi = previousw[0] + ogcp2[1]; mpi = 0;
		mi = previousw[0] + ogcp2[1] * gapfreq1pt[i-1]; mpi=0;
		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;
		fgcp2pt = fgcp2;
		ogcp2pt = ogcp2 + 1;
		fgcp1va = fgcp1[i-1];
		ogcp1va = ogcp1[i];
		gf1va = gapfreq1pt[i];
		gf1vapre = gapfreq1pt[i-1];
		gf2pt = gapfreq2pt+1;
		gf2ptpre = gapfreq2pt;

		if( trywarp )
		{
			prevwmrecordspt = prevwmrecords;
			wmrecordspt = wmrecords+1;
			wmrecords1pt = wmrecords;
			warpipt = warpi + 1;
			warpjpt = warpj + 1;
		}


		if( rowfused )
		{
#ifdef PROGTIMING
			if( pt_dpcheck() ) // verification mode: fused row on shadow copies, then the original match_calc (above) + loop
			{
				double t0, t1;
				double chk_mi = mi, chk_wm = wm; int chk_mpi = mpi;
				memcpy( chk_m, m, sizeof( double ) * ( lastj ) );
				memcpy( chk_mp, mp, sizeof( int ) * ( lastj ) );
				memcpy( chk_cur, currentw, sizeof( double ) * ( lastj ) );
				t0 = pt_now();
				if( rowfused == 2 )
					Adp_row_denseK( i, lgth2, fsel_mask, fsel_su, (const double * const *)gen_cu, gen_bm, chk_m, chk_mp, chk_cur, chk_ijp, previousw, gapfreq2pt, fgcp2, ogcp2,
						fgcp1va, ogcp1va, gf1va, gf1vapre, fpenalty_ex, &chk_mi, &chk_mpi, &chk_wm );
				else
				Adp_row_fused4( i, lgth2, chk_m, chk_mp, chk_cur, chk_ijp, previousw, gapfreq2pt, fgcp2, ogcp2,
					fsel_c[0], fsel_c[1], fsel_c[2], fsel_c[3], fsel_s,
					fgcp1va, ogcp1va, gf1va, gf1vapre, fpenalty_ex, &chk_mi, &chk_mpi, &chk_wm );
				t1 = pt_now();
				chk_cur[0] = currentw[0]; // column 0 is set outside the row kernels
				pt_acc[13] += t1 - t0;
				pt_chkrow_pending = 1;
				pt_chk_mi = chk_mi; pt_chk_mpi = chk_mpi; pt_chk_wm = chk_wm; pt_chk_t = t1;
				goto originalloop;
			}
			{ double _f0 = pt_now();
#endif
			if( rowfused == 2 ) // work/prog3
				Adp_row_denseK( i, lgth2, fsel_mask, fsel_su, (const double * const *)gen_cu, gen_bm, m, mp, currentw, ijp[i], previousw, gapfreq2pt, fgcp2, ogcp2,
					fgcp1va, ogcp1va, gf1va, gf1vapre, fpenalty_ex, &mi, &mpi, &wm );
			else
			Adp_row_fused4( i, lgth2, m, mp, currentw, ijp[i], previousw, gapfreq2pt, fgcp2, ogcp2,
				fsel_c[0], fsel_c[1], fsel_c[2], fsel_c[3], fsel_s,
				fgcp1va, ogcp1va, gf1va, gf1vapre, fpenalty_ex, &mi, &mpi, &wm );
#ifdef PROGTIMING
			{ double _f1 = pt_now(); pt_acc[23] += _f1 - _f0;
			if( v8s_stat_on() ) { if( rowfused == 2 ) { v8s_add( &v8s_gentime[0], _f1 - _f0 ); v8s_add( &v8s_gentime[1], (double)lgth2 ); } else v8s_add( &v8s_f4time, _f1 - _f0 ); } } }
#endif
			lastverticalw[i] = currentw[lgth2-1];
			continue;
		}

		if( !trywarp && FASTDP_ON )
		{
#ifdef PROGTIMING
			if( pt_dpcheck() ) // verification mode: run the fast row on shadow copies, then the original loop, and compare
			{
				double t0 = pt_now(), t1;
				double chk_mi = mi, chk_wm = wm; int chk_mpi = mpi;
				memcpy( chk_m, m, sizeof( double ) * ( lastj ) );
				memcpy( chk_mp, mp, sizeof( int ) * ( lastj ) );
				memcpy( chk_cur, currentw, sizeof( double ) * ( lastj ) );
				t0 = pt_now();
				Adp_row_fast( i, lastj, dpcolw, dpcolij, dpcolg2, chk_m, chk_mp, chk_cur, chk_ijp, previousw, gapfreq2pt, fgcp2, ogcp2, fgcp1va, ogcp1va, gf1va, gf1vapre, fpenalty_ex, &chk_mi, &chk_mpi, &chk_wm );
				t1 = pt_now();
				pt_acc[13] += t1 - t0;
				pt_chkrow_pending = 1;
				pt_chk_mi = chk_mi; pt_chk_mpi = chk_mpi; pt_chk_wm = chk_wm; pt_chk_t = t1;
				goto originalloop;
			}
#endif
			Adp_row_fast( i, lastj, dpcolw, dpcolij, dpcolg2, m, mp, currentw, ijp[i], previousw, gapfreq2pt, fgcp2, ogcp2, fgcp1va, ogcp1va, gf1va, gf1vapre, fpenalty_ex, &mi, &mpi, &wm );
			lastverticalw[i] = currentw[lgth2-1];
			continue;
		}
#ifdef PROGTIMING
		originalloop:
#endif

		for( j=1; j<lastj; j++ )
		{
#ifdef xxxenablemultithread
//			fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
			if( chudanpt && *chudanpt != chudanref ) 
			{
//				fprintf( stderr, "\n\n## CHUUDAN!!! S\n" );
				// work/prog2: mseq lives in a per-thread arena
				*chudanres = 1;
				return( -1.0 );
			}
#endif
			wm = *prept;
			*ijppt = 0;

#if 0
			fprintf( stderr, "\n i=%d, j=%d %c, %c", i, j, seq1[0][i], seq2[0][j] );
			fprintf( stderr, "%5.0f->", wm );
			fprintf( stderr, "%5.0f? (penal=%5.2f)", g=mi+*fgcp2pt*(1.0-gapfreq1pt[i]), *fgcp2pt*(1.0-gapfreq1pt[i]) );
#endif
			if( (g=mi+*fgcp2pt*gf1va) > wm )
			{
				wm = g;
				*ijppt = -( j - mpi );
//				fprintf( stderr, "Jump to %d (%c)!", mpi, seq2[0][mpi] );
			}
			if( (g=*prept+*ogcp2pt*gf1vapre) >= mi )
//			if( (g=*prept+*ogcp2pt*gf1vapre) > mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

#if 0 
			fprintf( stderr, "%5.0f->", wm );
			fprintf( stderr, "%5.0f? (penal=%5.2f)", g=*mjpt+fgcp1va*(1.0-gapfreq2pt[j]), fgcp1va*(1.0-gapfreq2pt[j]) );
#endif
			if( (g=*mjpt+ fgcp1va* *gf2pt) > wm )
			{
				wm = g;
				*ijppt = +( i - *mpjpt );
//				fprintf( stderr, "Jump to %d (%c)!", *mpjpt, seq1[0][*mpjpt] );
			}
			if( (g=*prept+ ogcp1va* *gf2ptpre) >= *mjpt )
//			if( (g=*prept+ ogcp1va* *gf2ptpre) > *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}
#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif


			if( trywarp )
			{
#if USE_PENALTY_EX
				if( ( g=*prevwmrecordspt++ + fpenalty_shift + fpenalty_ex * ( i - prevwarpi[j-1] + j - prevwarpj[j-1] ) ) > wm ) // naka ha osokute kamawanai
#else
				if( ( g=*prevwmrecordspt++ + fpenalty_shift ) > wm ) // naka ha osokute kamawanai
#endif
				{
//					fprintf( stderr, "WARP!!\n" );
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
	
#if 0
				fprintf( stderr, "%5.0f ", wm );
#endif
				curm = *curpt + wm;
	
				if( *wmrecords1pt > *wmrecordspt )
				{
					*wmrecordspt = *wmrecords1pt;
					*warpipt  = *(warpipt-1);
					*warpjpt  = *(warpjpt-1);
				}
				if( curm > *wmrecordspt )
				{
					*wmrecordspt = curm;
					*warpipt = i;
					*warpjpt = j;
				}
				wmrecordspt++;
				wmrecords1pt++;
				warpipt++;
				warpjpt++;
			}

			*curpt++ += wm;
			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
			fgcp2pt++;
			ogcp2pt++;
			gf2ptpre++;
			gf2pt++;
		}
#ifdef PROGTIMING
		if( pt_chkrow_pending )
		{
			pt_acc[14] += pt_now() - pt_chk_t;
			pt_chkrow_pending = 0;
			if( memcmp( chk_m, m, sizeof( double ) * lastj ) || memcmp( chk_mp, mp, sizeof( int ) * lastj ) || memcmp( chk_cur, currentw, sizeof( double ) * lastj )
			 || memcmp( chk_ijp+1, ijp[i]+1, sizeof( int ) * ( lastj-1 ) ) || memcmp( &pt_chk_mi, &mi, sizeof( double ) ) || pt_chk_mpi != mpi || memcmp( &pt_chk_wm, &wm, sizeof( double ) ) )
			{
				reporterr( "\nDPCHECK MISMATCH in A__align row %d (lgth1=%d, lgth2=%d)\n", i, lgth1, lgth2 );
				exit( 1 );
			}
			pt_acc[15] += 1.0; // rows checked
		}
#endif
		lastverticalw[i] = currentw[lgth2-1];

		if( trywarp )
		{
			fltncpy( prevwmrecords, wmrecords, lastj );
			intncpy( prevwarpi, warpi, lastj );
			intncpy( prevwarpj, warpj, lastj );
		}
	}
	adp_rowsdone:
#if defined( PROGTIMING ) && defined( enablemultithread )
	if( tilechk_pending ) tilechk_compare( lgth1, lgth2, lasti, m, mp, ijp, lastverticalw, wm );
#endif

	if( trywarp )
	{
//		fprintf( stderr, "wm = %f\n", wm );
//		fprintf( stderr, "warpn = %d\n", warpn );
		free( wmrecords );
		free( prevwmrecords );
		free( warpi );
		free( warpj );
		free( prevwarpi );
		free( prevwarpj );
	}

#if OUTGAP0TRY
	if( !outgap )
	{
		for( j=1; j<lgth2+1; j++ )
			currentw[j] -= offset * ( lgth2 - j ) / 2.0;
		for( i=1; i<lgth1+1; i++ )
			lastverticalw[i] -= offset * ( lgth1 - i  / 2.0);
	}
#endif
		
	/*
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr,"%s\n", seq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr,"%s\n", seq2[j] );
	fprintf( stderr, "====>" );
	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], seq1[i] );
	for( j=0; j<jcyc; j++ ) strcpy( mseq2[j], seq2[j] );
	*/

	PT_ADD(8);
	gt1 = gt1bk = AllocateCharVec( lgth1+lgth2+1 );
	gt2 = gt2bk = AllocateCharVec( lgth1+lgth2+1 );

	if( constraint )
	{
		if( atr_tiled ) ErrorExit( "work/v8s: band-major traceback with a constraint.\n" ); // (tiling excludes constraints)
		Atracking_localhom( impmatch, currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, warpis, warpjs, warpbase, &ngap1, &ngap2, reuseprofiles, &gt1, &gt2 );
	}
	else
	{
//		wmo = Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, tailgp, warpis, warpjs, warpbase, &ngap1, &ngap2, reuseprofiles, eff1, eff2, cpmxresult, cpmx1pt, cpmx2pt, gapfreq1pt, gapfreq2pt, ogcp1opt, ogcp2opt, fgcp1opt, fgcp2opt, orieff1, orieff2, (cpmx1pt!=cpmx1), (cpmx2pt!=cpmx2) );
		/* work/prog2: Atracking only builds the gap tables here; the gaps are inserted into seq1/seq2
		   in place below (gapireru_inplace), instead of into mseq1/mseq2 followed by a strcpy back.
		   Rows without inserted gaps (ngap == 0) are left as they are: the original copied them to
		   mseq and back unchanged. */
		if( lz.on ) { atr_len1 = lgth1; atr_len2 = lgth2; }
		wmo = Atracking( currentw, lastverticalw, fpenalty, fpenalty_ex, seq1, seq2, NULL, NULL, ijp, icyc, jcyc, tailgp, warpis, warpjs, warpbase, &ngap1, &ngap2, reuseprofiles, &gt1, &gt2 );
		atr_len1 = atr_len2 = -1;
		atr_tiled = 0; // work/v8s
		if( !tailgp ) wm = wmo;
		inplace = 1;
		{
			int lnew1 = gt1bk + lgth1+lgth2 - gt1, lnew2 = gt2bk + lgth1+lgth2 - gt2;
			resultlen = ( ngap1 == 0 && reuseprofiles ) ? 0 : lnew1; // = strlen( mseq1[0] ) in the original
			if( alloclen < resultlen || resultlen > N )
			{
				fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
				ErrorExit( "LENGTH OVER!\n" );
			}
			if( lz.on ) // work/prog3: hand deferred gap tables to the caller (lazygap.h)
			{
				A_lz.outlen = lnew1;
				A_lz.used = 1;
				A_lz.outgt[0] = A_lz.outgt[1] = NULL;
				if( lz.defer[0] && ngap1 ) { A_lz.outgt[0] = malloc( lnew1 + 1 ); memcpy( A_lz.outgt[0], gt1, lnew1 ); A_lz.outgt[0][lnew1] = 0; }
				if( lz.defer[1] && ngap2 ) { A_lz.outgt[1] = malloc( lnew2 + 1 ); memcpy( A_lz.outgt[1], gt2, lnew2 ); A_lz.outgt[1][lnew2] = 0; }
			}
			if( ngap1 && !lz.defer[0] ) gapexp_all( seq1, icyc, gt1, lgth1, lnew1 );
			if( ngap2 && !lz.defer[1] ) gapexp_all( seq2, jcyc, gt2, lgth2, lnew2 );
		}
	}

	PT_ADD(9);
//	reporterr( "wm = %f\n", wm );

#if 1
	if( cpmxresult )
	{
		if( icyc + jcyc > 20 )
//		if( 0 )
//		if( 1 )
		{
			int limk = gt1bk + lgth1+lgth2 - gt1;
			limk = (int)(( (limk+1000) / 1000 ) * 1000); // kiriage
			//reporterr( "limk=%d\n", limk );

#if 1 // marume gosa wo teigenn suru tame
			double totaleff1 = 0.0;
			double totaleff2 = 0.0;
			for( i=0; i<icyc; i++ ) totaleff1 += eff1[i];
			for( j=0; j<jcyc; j++ ) totaleff2 += eff2[j];

			if( fabs(totaleff1-1.0)>0.001 || fabs(totaleff2-1.0)>0.001 )
			{
				reporterr( "totaleff1 = %50.40f\n", totaleff1 );
				reporterr( "totaleff2 = %50.40f\n", totaleff2 );
				exit( 1 );
			}
			totaleff1 = totaleff1 * orieff1 / (orieff1 + orieff2);
			totaleff2 = totaleff2 * orieff2 / (orieff1 + orieff2);
#else

//			reporterr( "totaleff1 (rescaled    ) = %20.10f\n", totaleff1 );
//			reporterr( "totaleff2 (rescaled    ) = %20.10f\n", totaleff2 );
//			double bk1 = totaleff1;
//			double bk2 = totaleff2;

			double totaleff1 = orieff1 / (orieff1 + orieff2);
			double totaleff2 = orieff2 / (orieff1 + orieff2);

//			reporterr( "totaleff1 (not rescaled) = %20.10f\n", totaleff1 );
//			reporterr( "totaleff2 (not rescaled) = %20.10f\n", totaleff2 );

//			reporterr( "ratio1 = %30.25f\n", totaleff1/bk1 );
//			reporterr( "ratio2 = %30.25f\n", totaleff2/bk2 );
#endif
	
			*cpmxresult = AllocateDoubleMtx( nalphabets+3, 0 ); // gapcount, opg, fng no bun
			createcpmxresult( *cpmxresult, limk, totaleff1, totaleff2, cpmx1pt, cpmx2pt, gt1, gt2, (cpmx1!=*cpmx1pt), (cpmx2!=*cpmx2pt) ); // naka de free
			creategapfreqresult( *cpmxresult+nalphabets, limk, totaleff1, totaleff2, gapfreq1pt, gapfreq2pt, gt1, gt2 ); // naka deha free shinai
			// gapfreq1, gapfreq2 ha mada tsukau
			createogresult( *cpmxresult+nalphabets+1, limk, totaleff1, totaleff2, ogcp1opt, ogcp2opt, gapfreq1pt, gapfreq2pt, gt1, gt2 ); // naka deha free shinai
			if( cpmx1!=*cpmx1pt ) free( ogcp1opt );
			if( cpmx2!=*cpmx2pt ) free( ogcp2opt );
			createfgresult( *cpmxresult+nalphabets+2, limk, totaleff1, totaleff2, fgcp1opt, fgcp2opt, gapfreq1pt, gapfreq2pt, gt1, gt2 ); // naka deha free shinai
			if( cpmx1!=*cpmx1pt ) free( fgcp1opt );
			if( cpmx2!=*cpmx2pt ) free( fgcp2opt );

			if( cpmx1!=*cpmx1pt ) free( gapfreq1pt );
			if( cpmx2!=*cpmx2pt ) free( gapfreq2pt );
#if 0
			reporterr( "\n" );
			for( j=0; j<nalphabets; j++ ) 
				reporterr( "%4c ", amino[j] );
			reporterr( "\ncpmxresult = \n" );
			for( i=0; i<strlen( gt1 ); i++ ) 
			{
				for( j=0; j<nalphabets; j++ ) 
					reporterr( "%4.2f ", (*cpmxresult)[j][i] );
				reporterr( "\n" );
			}
#endif
		}
		else
			*cpmxresult = NULL;
	}
#endif

	free( gt1bk );
	free( gt2bk );
#if 0 // 2021/Jun/24
	if( *cpmx1pt != cpmx1 ) 
	{
		if( cpmxchild0 && *cpmxchild0 )
		{
			reporterr( "freeing *cpmxchild0\n" );
			FreeDoubleMtx( *cpmxchild0 );
			*cpmxchild0 = NULL;
		}
	}
	if( *cpmx2pt != cpmx2 ) 
	{
		if( cpmxchild1 && *cpmxchild1 )
		{
			reporterr( "freeing *cpmxchild1\n" );
			FreeDoubleMtx( *cpmxchild1 );
			*cpmxchild1 = NULL;
		}
	}
#endif

//	reporterr( "reuseprofiles after Atracking = %d\n", reuseprofiles );

	if( warpis ) free( warpis );
	if( warpjs ) free( warpjs );

//	fprintf( stderr, "### impmatch = %f\n", *impmatch );

	if( !inplace ) // work/prog2: otherwise done in place after Atracking
	{
	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}


	if( ngap1 || !reuseprofiles )
		for( i=0; i<icyc; i++ ) strcpy( seq1[i], mseq1[i] );
	if( ngap2 || !reuseprofiles )
		for( j=0; j<jcyc; j++ ) strcpy( seq2[j], mseq2[j] );
	}



#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "%s\n", mseq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, "%s\n", mseq2[j] );
//exit( 1 );
#endif

//	fprintf( stdout, "firstmem=%d, icyc=%d, jcyc=%d, wm = %f\n", firstmem, icyc, jcyc, wm );

//	fprintf( stderr, "lgth1 = %d\n", lgth1 );
//	fprintf( stderr, "->      %d\n", strlen( seq1[0] ) );
	previousfirstlen = lgth1;
	previousfirstmem = firstmem;
	previousicyc = icyc;
	previouscall = calledbyfulltreebase;

	PT_ADD(10);
	// work/prog2: mseq lives in a per-thread arena
	A_ijp_release(); // work/prog2: large traceback blocks go back to the shared pool

	return( wm );
}

double A__align_gapmap( char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen, int constraint, double *impmatch, int *gapmap1, int *gapmap2 )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
	fprintf( stderr, "Unexpected error.  Please contact katoh@ifrec.osaka-u.ac.jp\n" );
	exit( 1 );
}


double A__align_variousdist( int **which, double ***matrices, double **n_dynamicmtx, int penalty_l, int penalty_ex_l, char **seq1, char **seq2, double *eff1, double *eff2, double **eff1s, double **eff2s, int icyc, int jcyc, int alloclen, int constraint, double *impmatch, char *sgap1, char *sgap2, char *egap1, char *egap2, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{


//	int k;
	register int i, j, c;
	int ngap1, ngap2;
	int lasti, lastj;      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lgth1, lgth2;
	int resultlen;
	double wm = 0.0;   /* int ?????? */
	double wmo = 0.0;   /* int ?????? */
	double g;
	double *currentw, *previousw;
//	double fpenalty = (double)penalty;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penalty_ex_l;
#endif
#if 1
	double *wtmp;
	int *ijppt;
	double *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	static TLS double mi, *m;
	static TLS int **ijp;
	static TLS int mpi, *mp;
	static TLS double *w1, *w2;
	static TLS double *match;
	static TLS double *initverticalw;    /* kufuu sureba iranai */
	static TLS double *lastverticalw;    /* kufuu sureba iranai */
	static TLS char **mseq1;
	static TLS char **mseq2;
	static TLS char **mseq;
	static TLS double *ogcp1;
	static TLS double *ogcp2;
	static TLS double *fgcp1;
	static TLS double *fgcp2;
	static TLS double ***cpmx1s;
	static TLS double ***cpmx2s;
	static TLS int ***intwork;
	static TLS double ***doublework;
	static TLS int orlgth1 = 0, orlgth2 = 0;
	static TLS double *gapfreq1;
	static TLS double *gapfreq2;
	double fpenalty = (double)penalty_l;
	double fpenalty_shift = (double)penalty_shift;
	double *fgcp2pt;
	double *ogcp2pt;
	double fgcp1va;
	double ogcp1va;
	double *gf2pt;
	double *gf2ptpre;
	double gf1va;
	double gf1vapre;
	double headgapfreq1;
	double headgapfreq2;

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
	int *nmask, **masklist1, **masklist2;


	if( seq1 == NULL )
	{
		if( orlgth1 )
		{
//			fprintf( stderr, "## Freeing local arrays in A__align\n" );
			orlgth1 = 0;
			orlgth2 = 0;

			imp_match_init_strict( NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL, 0, 0 );

			free( mseq1 );
			free( mseq2 );
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );

			FreeCharMtx( mseq );

			FreeFloatVec( ogcp1 );
			FreeFloatVec( ogcp2 );
			FreeFloatVec( fgcp1 );
			FreeFloatVec( fgcp2 );


			FreeFloatCub( cpmx1s );
			FreeFloatCub( cpmx2s );

			FreeFloatVec( gapfreq1 );
			FreeFloatVec( gapfreq2 );

			FreeFloatCub( doublework );
			FreeIntCub( intwork );

		}
		else
		{
//			fprintf( stderr, "## Not allocated\n" );
		}
		return( 0.0 );
	}


#if SLOW
	nmask = calloc( maxdistclass, sizeof( int ) );
#else
	masklist1 = AllocateIntMtx( maxdistclass, 0 );
	masklist2 = AllocateIntMtx( maxdistclass, 0 );
	nmask = calloc( maxdistclass, sizeof( int ) );

	for( c=0; c<maxdistclass; c++ )
	{
		for( i=0; i<icyc; i++ ) for( j=0; j<jcyc; j++ )
		{
			if( eff1s[c][i] * eff2s[c][j] != 0.0 )
			{

				if( c != which[i][j] )
				{
					masklist1[c] = realloc( masklist1[c], sizeof( int ) * (nmask[c]+1) );
					masklist2[c] = realloc( masklist2[c], sizeof( int ) * (nmask[c]+1) );

					masklist1[c][nmask[c]] = i;
					masklist2[c][nmask[c]] = j;
					nmask[c]++;
				}
			}
		}
	}
	for( c=0; c<maxdistclass; c++ ) if( nmask[c] ) break;
	if( c<maxdistclass ) reporterr( "Found a complex grouping. This step may be a bit slow.\n" );
#endif

	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );
#if 0
	if( lgth1 == 0 || lgth2 == 0 )
	{
		fprintf( stderr, "WARNING (Aalignmm): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif
	if( lgth1 == 0 && lgth2 == 0 )
		return( 0.0 );

	if( lgth1 == 0 )
	{
		for( i=0; i<icyc; i++ )
		{
			j = lgth2;
			seq1[i][j] = 0;
			while( j ) seq1[i][--j] = *newgapstr;
//			fprintf( stderr, "seq1[i] = %s\n", seq1[i] );
		}
		return( 0.0 );
	}

	if( lgth2 == 0 )
	{
		for( i=0; i<jcyc; i++ )
		{
			j = lgth1;
			seq2[i][j] = 0;
			while( j ) seq2[i][--j] = *newgapstr;
//			fprintf( stderr, "seq2[i] = %s\n", seq2[i] );
		}
		return( 0.0 );
	}

	warpbase = lgth1 + lgth2;
	warpis = NULL;
	warpjs = NULL;
	warpn = 0;

	if( trywarp )
	{
//		fprintf( stderr, "In A__align_variousdist !!!!!\n" );
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
		for( i=0; i<lgth2+1; i++ ) wmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwarpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) prevwarpj[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpj[i] = -warpbase;
	}

#if 0
	fprintf( stderr, "####  eff in SA+++align\n" );
	fprintf( stderr, "####  seq1[0] = %s\n", seq1[0] );
	fprintf( stderr, "####  strlen( seq1[0] ) = %d\n", strlen( seq1[0] ) );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "eff1[%d] = %f\n", i, eff1[i] );
	fprintf( stderr, "####  seq2[0] = %s\n", seq2[0] );
	fprintf( stderr, "####  strlen( seq2[0] ) = %d\n", strlen( seq2[0] ) );
	for( i=0; i<jcyc; i++ ) fprintf( stderr, "eff2[%d] = %f\n", i, eff2[i] );
#endif
	if( orlgth1 == 0 )
	{
		mseq1 = AllocateCharMtx( njob, 0 );
		mseq2 = AllocateCharMtx( njob, 0 );
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

			FreeFloatVec( ogcp1 );
			FreeFloatVec( ogcp2 );
			FreeFloatVec( fgcp1 );
			FreeFloatVec( fgcp2 );


			FreeFloatCub( cpmx1s );
			FreeFloatCub( cpmx2s );

			FreeFloatVec( gapfreq1 );
			FreeFloatVec( gapfreq2 );

			FreeFloatCub( doublework );
			FreeIntCub( intwork );
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

		mseq = AllocateCharMtx( njob, ll1+ll2 );

		ogcp1 = AllocateFloatVec( ll1+2 );
		ogcp2 = AllocateFloatVec( ll2+2 );
		fgcp1 = AllocateFloatVec( ll1+2 );
		fgcp2 = AllocateFloatVec( ll2+2 );

		cpmx1s = AllocateFloatCub( maxdistclass, nalphabets, ll1+2 );
		cpmx2s = AllocateFloatCub( maxdistclass, nalphabets, ll2+2 );

		gapfreq1 = AllocateFloatVec( ll1+2 );
		gapfreq2 = AllocateFloatVec( ll2+2 );

		doublework = AllocateFloatCub( maxdistclass, MAX( ll1, ll2 )+2, nalphabets ); 
		intwork = AllocateIntCub( maxdistclass, MAX( ll1, ll2 )+2, nalphabets+1 ); 

#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;
	}


	for( i=0; i<icyc; i++ )
	{
		mseq1[i] = mseq[i];
		seq1[i][lgth1] = 0;
	}
	for( j=0; j<jcyc; j++ )
	{
		mseq2[j] = mseq[icyc+j];
		seq2[j][lgth2] = 0;
	}


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
	{
		double t = 0.0;
		for( i=0; i<icyc; i++ )
			t += eff1[i];
	fprintf( stderr, "## totaleff = %f\n", t );
	}
#endif

#if SLOW
#else
//	cpmx_calc_new( seq1, cpmx1, eff1, lgth1, icyc );
//	cpmx_calc_new( seq2, cpmx2, eff2, lgth2, jcyc );
	for( c=0; c<maxdistclass; c++ )
	{
		cpmx_calc_new( seq1, cpmx1s[c], eff1s[c], lgth1, icyc );
		cpmx_calc_new( seq2, cpmx2s[c], eff2s[c], lgth2, jcyc );
	}
#endif

	if( sgap1 )
	{
		new_OpeningGapCount( ogcp1, icyc, seq1, eff1, lgth1, sgap1 );
		new_OpeningGapCount( ogcp2, jcyc, seq2, eff2, lgth2, sgap2 );
		new_FinalGapCount( fgcp1, icyc, seq1, eff1, lgth1, egap1 );
		new_FinalGapCount( fgcp2, jcyc, seq2, eff2, lgth2, egap2 );
		outgapcount( &headgapfreq1, icyc, sgap1, eff1 );
		outgapcount( &headgapfreq2, jcyc, sgap2, eff2 );
		outgapcount( gapfreq1+lgth1, icyc, egap1, eff1 );
		outgapcount( gapfreq2+lgth2, jcyc, egap2, eff2 );
	}
	else
	{
		st_OpeningGapCount( ogcp1, icyc, seq1, eff1, lgth1 );
		st_OpeningGapCount( ogcp2, jcyc, seq2, eff2, lgth2 );
		st_FinalGapCount( fgcp1, icyc, seq1, eff1, lgth1 );
		st_FinalGapCount( fgcp2, jcyc, seq2, eff2, lgth2 );
		headgapfreq1 = 0.0;
		headgapfreq2 = 0.0;
		gapfreq1[lgth1] = 0.0;
		gapfreq2[lgth2] = 0.0;
	}

	if( legacygapcost == 0 )
	{
		gapcountf( gapfreq1, seq1, icyc, eff1, lgth1 );
		gapcountf( gapfreq2, seq2, jcyc, eff2, lgth2 );
		for( i=0; i<lgth1+1; i++ ) gapfreq1[i] = 1.0 - gapfreq1[i];
		for( i=0; i<lgth2+1; i++ ) gapfreq2[i] = 1.0 - gapfreq2[i];
		headgapfreq1 = 1.0 - headgapfreq1;
		headgapfreq2 = 1.0 - headgapfreq2;
	}
	else
	{
		for( i=0; i<lgth1+1; i++ ) gapfreq1[i] = 1.0;
		for( i=0; i<lgth2+1; i++ ) gapfreq2[i] = 1.0;
		headgapfreq1 = 1.0;
		headgapfreq2 = 1.0;
	}

#if 0
	fprintf( stderr, "\ngapfreq1[] =" );
	for( i=0; i<lgth1; i++ ) fprintf( stderr, "%5.2f ", gapfreq1[i] );
	fprintf( stderr, "\n" );

	fprintf( stderr, "\ngapfreq2[] =" );
	for( i=0; i<lgth2; i++ ) fprintf( stderr, "%5.2f ", gapfreq2[i] );
	fprintf( stderr, "\n" );
#endif
	

	for( i=0; i<lgth1; i++ ) 
	{
		ogcp1[i] = 0.5 * ( 1.0 - ogcp1[i] ) * fpenalty * ( gapfreq1[i] );
		fgcp1[i] = 0.5 * ( 1.0 - fgcp1[i] ) * fpenalty * ( gapfreq1[i] );
	}

	for( i=0; i<lgth2; i++ ) 
	{
		ogcp2[i] = 0.5 * ( 1.0 - ogcp2[i] ) * fpenalty * ( gapfreq2[i] );
		fgcp2[i] = 0.5 * ( 1.0 - fgcp2[i] ) * fpenalty * ( gapfreq2[i] );
	}
#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

	currentw = w1;
	previousw = w2;

//	for( i=0; i<icyc; i++ ) fprintf( stderr, "seq1[i] = %s\n", seq1[i] );
//	for( j=0; j<jcyc; j++ ) fprintf( stderr, "seq2[j] = %s\n", seq2[j] );

#if SLOW
	match_calc_slow( which, matrices, initverticalw, jcyc, seq2, eff2, icyc, seq1, eff1, 0, lgth1, *doublework, *intwork, 1, 1 );
#else
	fillzero( initverticalw, lgth1 );
	for( c=0; c<maxdistclass; c++ )
	{
//		fprintf( stderr, "c=%d matrices[c][W][W] = %f\n", c, matrices[c][amino_n['W']][amino_n['W']] );
//		for( i=0; i<lgth1; i++ ) fprintf( stderr, "seq1[i] = %c, cpmx1s[c][3][%d] = %f\n", seq1[0][i], i, cpmx1s[c][3][i] );
//		for( i=0; i<lgth2; i++ ) fprintf( stderr, "seq2[i] = %c, cpmx2s[c][3][%d] = %f\n", seq2[0][i], i, cpmx2s[c][3][i] );
		match_calc_add( matrices[c], initverticalw, cpmx2s[c], cpmx1s[c], 0, lgth1, doublework[c], intwork[c], 1 );
//		for( i=0; i<lgth1; i++ ) fprintf( stderr, "c=%d, %d - %f\n", c, i, initverticalw[i] );
		if( nmask[c] ) match_calc_del( which, matrices, initverticalw, jcyc, seq2, eff2, icyc, seq1, eff1, 0, lgth1, c, nmask[c], masklist2[c], masklist1[c] );
	}
//	for( i=0; i<lgth1; i++ ) fprintf( stderr, "%d - %f\n", i, initverticalw[i] );
#endif

//	exit( 1 );

	if( constraint )
		imp_match_out_vead_tate( initverticalw, 0, lgth1 ); // 060306

#if SLOW
	match_calc_slow( which, matrices, currentw, icyc, seq1, eff1, jcyc, seq2, eff2, 0, lgth2, *doublework, *intwork, 1, 0 );
//	for( i=0; i<lgth2; i++ ) fprintf( stderr, "%d - %f\n", i, currentw[i] );
//	exit( 1 );
#else
	fillzero( currentw, lgth2 );
	for( c=0; c<maxdistclass; c++ )
	{
		match_calc_add( matrices[c], currentw, cpmx1s[c], cpmx2s[c], 0, lgth2, doublework[c], intwork[c], 1 );
		if( nmask[c] ) match_calc_del( which, matrices, currentw, icyc, seq1, eff1, jcyc, seq2, eff2, 0, lgth2, c, nmask[c], masklist1[c], masklist2[c] );
	}
//	for( i=0; i<lgth2; i++ ) fprintf( stderr, "%d - %f\n", i, currentw[i] );
//	exit( 1 );
#endif

	if( constraint )
		imp_match_out_vead( currentw, 0, lgth2 ); // 060306
#if 0 // -> tbfast.c // impossible
	if( localhom )
		imp_match_calc( n_dynamicmtx, currentw, icyc, jcyc, lgth1, lgth2, seq1, seq2, eff1, eff2, localhom, 1, 0 );

#endif

	if( headgp == 1 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
//			initverticalw[i] += ( ogcp1[0] + fgcp1[i-1] ) ;
			initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2[0] ) ;
#if USE_PENALTY_EX // 2018/Apr/23
			initverticalw[i] += fpenalty_ex * i;
//			reporterr( "added _ex\n" );
#endif
		}
		for( j=1; j<lgth2+1; j++ )
		{
//			currentw[j] += ( ogcp2[0] + fgcp2[j-1] ) ;
			currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1[0] ) ;
#if USE_PENALTY_EX // 2018/Apr/23
			currentw[j] += fpenalty_ex * j;
//			reporterr( "added _ex\n" );
#endif
		}
	}
	else
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += fpenalty * TERMGAPFAC; // motto tanjun de yoi?
#if USE_PENALTY_EX // 2018/Apr/22
			initverticalw[i] += fpenalty_ex * i * TERMGAPFAC_EX;
//			reporterr( "added _ex\n" );
#endif
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += fpenalty * TERMGAPFAC; // motto tanjun de yoi?
#if USE_PENALTY_EX // 2018/Apr/22
			currentw[j] += fpenalty_ex * j * TERMGAPFAC_EX;
//			reporterr( "added _ex\n" );
#endif
		}
	}


#if OUTGAP0TRY
	else
	{
		fprintf( stderr, "offset = %d\n", offset );
		for( j=1; j<lgth2+1; j++ )
			currentw[j] -= offset * j / 2.0;
		for( i=1; i<lgth1+1; i++ )
			initverticalw[i] -= offset * i / 2.0;
	}
#endif
#if 0
	fprintf( stderr, "\n   " );
	for( j=0; j<lgth2+1; j++ ) fprintf( stderr, "    %c ", seq2[0][j] );
	fprintf( stderr, "\n%c ", seq1[0][0]  );
	for( j=0; j<lgth2+1; j++ )
	{
		fprintf( stderr, "%5.0f ", currentw[j] );
	}
	fprintf( stderr, "\n"  );
#endif


	for( j=1; j<lgth2+1; ++j ) 
	{
//		m[j] = currentw[j-1] + ogcp1[1]; mp[j] = 0;
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2[j-1]; mp[j] = 0;;
	}
	if( lgth2 == 0 )
		lastverticalw[0] = 0.0; // Falign kara yobaretatoki kounarukanousei ari
	else
		lastverticalw[0] = currentw[lgth2-1];

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
fprintf( stderr, "fcgp\n" );
for( i=0; i<lgth1; i++ ) 
	fprintf( stderr, "fgcp1[%d]=%f\n", i, ogcp1[i] );
for( i=0; i<lgth2; i++ ) 
	fprintf( stderr, "fgcp2[%d]=%f\n", i, ogcp2[i] );
#endif

	for( i=1; i<lasti; i++ )
	{

#ifdef enablemultithread
//		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
		if( chudanpt && *chudanpt != chudanref ) 
		{
//			fprintf( stderr, "\n\n## CHUUDAN!!! S\n" );
			if( masklist1 ) freeintmtx( masklist1, maxdistclass ); masklist1 = NULL;
			if( masklist2 ) freeintmtx( masklist2, maxdistclass ); masklist2 = NULL;
			if( nmask ) free( nmask ); nmask = NULL;
			*chudanres = 1;
			return( -1.0 );
		}
#endif
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

#if SLOW
		match_calc_slow( which, matrices, currentw, icyc, seq1, eff1, jcyc, seq2, eff2, i, lgth2, *doublework, *intwork, 0, 0 );
#else
		fillzero( currentw, lgth2 );
		for( c=0; c<maxdistclass; c++ )
		{
			match_calc_add( matrices[c], currentw, cpmx1s[c], cpmx2s[c], i, lgth2, doublework[c], intwork[c], 0 );
			if( nmask[c] ) match_calc_del( which, matrices, currentw, icyc, seq1, eff1, jcyc, seq2, eff2, i, lgth2, c, nmask[c], masklist1[c], masklist2[c] );
		}
#endif
#if 0
		if( i == 1 )
		{
			fprintf( stderr, "\n" );
			for( j=0; j<lgth2; j++ ) fprintf( stderr, "%d - %f\n", j, currentw[j] );
			exit( 1 );
		}
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
		if( constraint )
		{
//			fprintf( stderr, "Calling imp_match_calc (o) lgth = %d, i = %d\n", lgth1, i );
#if  0
			imp_match_out_vead( currentw, i, lgth2 );
#else
			imp_match_out_vead( currentw, i, lgth2 );
#endif
		}
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

#if 0
		fprintf( stderr, "%c ", seq1[0][i] );
		for( j=0; j<lgth2+1; j++ )
		{
			fprintf( stderr, "%5.0f ", currentw[j] );
		}
		fprintf( stderr, "\n"  );
#endif
	
//		mi = previousw[0] + ogcp2[1]; mpi = 0;
		mi = previousw[0] + ogcp2[1] * gapfreq1[i-1]; mpi=0;
		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;
		fgcp2pt = fgcp2;
		ogcp2pt = ogcp2 + 1;
		fgcp1va = fgcp1[i-1];
		ogcp1va = ogcp1[i];
		gf1va = gapfreq1[i];
		gf1vapre = gapfreq1[i-1];
		gf2pt = gapfreq2+1;
		gf2ptpre = gapfreq2;

		if( trywarp )
		{
			prevwmrecordspt = prevwmrecords;
			wmrecordspt = wmrecords+1;
			wmrecords1pt = wmrecords;
			warpipt = warpi + 1;
			warpjpt = warpj + 1;
		}

		for( j=1; j<lastj; j++ )
		{
#ifdef xxxenablemultithread
//			fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
			if( chudanpt && *chudanpt != chudanref ) 
			{
//				fprintf( stderr, "\n\n## CHUUDAN!!! S\n" );
				if( masklist1 ) freeintmtx( masklist1, maxdistclass ); masklist1 = NULL;
				if( masklist2 ) freeintmtx( masklist2, maxdistclass ); masklist2 = NULL;
				if( nmask ) free( nmask ); nmask = NULL;
				*chudanres = 1;
				return( -1.0 );
			}
#endif
			wm = *prept;
			*ijppt = 0;

#if 0
			fprintf( stderr, "\n i=%d, j=%d %c, %c", i, j, seq1[0][i], seq2[0][j] );
			fprintf( stderr, "%5.0f->", wm );
			fprintf( stderr, "%5.0f? (penal=%5.2f)", g=mi+*fgcp2pt*(1.0-gapfreq1[i]), *fgcp2pt*(1.0-gapfreq1[i]) );
#endif
			if( (g=mi+*fgcp2pt*gf1va) > wm )
			{
				wm = g;
				*ijppt = -( j - mpi );
//				fprintf( stderr, "Jump to %d (%c)!", mpi, seq2[0][mpi] );
			}
			if( (g=*prept+*ogcp2pt*gf1vapre) >= mi )
//			if( (g=*prept+*ogcp2pt*gf1vapre) > mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

#if 0 
			fprintf( stderr, "%5.0f->", wm );
			fprintf( stderr, "%5.0f? (penal=%5.2f)", g=*mjpt+fgcp1va*(1.0-gapfreq2[j]), fgcp1va*(1.0-gapfreq2[j]) );
#endif
			if( (g=*mjpt+ fgcp1va* *gf2pt) > wm )
			{
				wm = g;
				*ijppt = +( i - *mpjpt );
//				fprintf( stderr, "Jump to %d (%c)!", *mpjpt, seq1[0][*mpjpt] );
			}
			if( (g=*prept+ ogcp1va* *gf2ptpre) >= *mjpt )
//			if( (g=*prept+ ogcp1va* *gf2ptpre) > *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}
#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif
			if( trywarp )
			{
#if USE_PENALTY_EX
				if( ( g=*prevwmrecordspt++ + fpenalty_shift + fpenalty_ex * ( i - prevwarpi[j-1] + j - prevwarpj[j-1] ) ) > wm ) // naka ha osokute kamawanai
#else
				if( ( g=*prevwmrecordspt++ + fpenalty_shift ) > wm ) // naka ha osokute kamawanai
#endif
				{
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
				curm = *curpt + wm;

				if( *wmrecords1pt > *wmrecordspt )
				{
					*wmrecordspt = *wmrecords1pt;
					*warpipt  = *(warpipt-1);
					*warpjpt  = *(warpjpt-1);
				}
				if( curm > *wmrecordspt )
				{
					*wmrecordspt = curm;
					*warpipt = i;
					*warpjpt = j;
				}
				wmrecordspt++;
				wmrecords1pt++;
				warpipt++;
				warpjpt++;
			}

#if 0
			fprintf( stderr, "%5.0f ", wm );
#endif
			*curpt++ += wm;
			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
			fgcp2pt++;
			ogcp2pt++;
			gf2ptpre++;
			gf2pt++;

		}
		lastverticalw[i] = currentw[lgth2-1];

		if( trywarp )
		{
			fltncpy( prevwmrecords, wmrecords, lastj );
			intncpy( prevwarpi, warpi, lastj );
			intncpy( prevwarpj, warpj, lastj );
		}
	}
	if( trywarp )
	{
//		fprintf( stderr, "wm = %f\n", wm );
//		fprintf( stderr, "warpn = %d\n", warpn );
		free( wmrecords );
		free( prevwmrecords );
		free( warpi );
		free( warpj );
		free( prevwarpi );
		free( prevwarpj );
	}


#if OUTGAP0TRY
	if( !outgap )
	{
		for( j=1; j<lgth2+1; j++ )
			currentw[j] -= offset * ( lgth2 - j ) / 2.0;
		for( i=1; i<lgth1+1; i++ )
			lastverticalw[i] -= offset * ( lgth1 - i  / 2.0);
	}
#endif
		
	/*
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr,"%s\n", seq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr,"%s\n", seq2[j] );
	fprintf( stderr, "====>" );
	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], seq1[i] );
	for( j=0; j<jcyc; j++ ) strcpy( mseq2[j], seq2[j] );
	*/
	if( constraint )
	{
		Atracking_localhom( impmatch, currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, warpis, warpjs, warpbase, &ngap1, &ngap2, 0, NULL, NULL );
	}
	else
	{
//		Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, tailgp, warpis, warpjs, warpbase, &ngap1, &ngap2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0.0, 0.0, 1, 1 ); // NULL x 11 ha atode awaseru.
		wmo = Atracking( currentw, lastverticalw, fpenalty, fpenalty_ex, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, tailgp, warpis, warpjs, warpbase, &ngap1, &ngap2, 0, NULL, NULL);
		if( !tailgp ) wm = wmo;
	}

	if( warpis ) free( warpis );
	if( warpjs ) free( warpjs );

//	fprintf( stderr, "### impmatch = %f\n", *impmatch );

	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}


	for( i=0; i<icyc; i++ ) strcpy( seq1[i], mseq1[i] );
	for( j=0; j<jcyc; j++ ) strcpy( seq2[j], mseq2[j] );
#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "%s\n", mseq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, "%s\n", mseq2[j] );
#endif

//	fprintf( stderr, "wm = %f\n", wm );


	if( masklist1 ) freeintmtx( masklist1, maxdistclass ); masklist1 = NULL;
	if( masklist2 ) freeintmtx( masklist2, maxdistclass ); masklist2 = NULL;
	if( nmask ) free( nmask ); nmask = NULL;

	return( wm );
}
