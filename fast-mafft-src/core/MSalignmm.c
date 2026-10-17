#ifndef _GNU_SOURCE
#define _GNU_SOURCE // work/v10w: glibc declares nanosleep, struct timespec and sched_yield only with a feature macro under -std=c99 (gcc >= 14 rejects implicit declarations)
#endif
#include "mltaln.h"
#include "dp.h"
#include "dphelp.h" // work/v10w
#include "progtiming.h" // work/v10w (MAFFT_MSSTAT in timing builds)
#include <sched.h> // work/v10w
#include <time.h> // work/v10w

#define MEMSAVE 1

#define DEBUG 0
#define USE_PENALTY_EX  0
#define FASTMATCHCALC 1
#define STOREWM 0

#define DPTANNI 100

#define ATO 1

static TLS int reccycle = 0;


// [seq][alphabet]
static void match_calc_add( double **scoringmtx, double *match, double **cpmx1, double **cpmx2, int i1, int lgth2, double **doublework, int **intwork, int initialize )
{
	int j, k, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	int count = 0;
	double *matchpt;
	double **cpmxpdpt;
	int **cpmxpdnpt;
	int cpkd;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
	if( initialize )
	{
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[j][l] )
				{
					cpmxpd[j][count] = cpmx2[j][l];
					cpmxpdn[j][count] = l;
					count++;
				}
			}
			cpmxpdn[j][count] = -1;
		}
	}

	for( l=0; l<nalphabets; l++ )
	{
		scarr[l] = 0.0;
		for( k=0; k<nalphabets; k++ )
		{
//			scarr[l] += n_dis[k][l] * cpmx1[i1][k];
			scarr[l] += scoringmtx[k][l] * cpmx1[i1][k];
		}
	}
#if 0 /* これを使うときはdoubleworkのアロケートを逆にする */
	{
		double *fpt, **fptpt, *fpt2;
		int *ipt, **iptpt;
		fpt2 = match;
		iptpt = cpmxpdn;
		fptpt = cpmxpd;
		while( lgth2-- )
		{
			*fpt2 = 0.0;
			ipt=*iptpt,fpt=*fptpt;
			while( *ipt > -1 )
				*fpt2 += scarr[*ipt++] * *fpt++;
			fpt2++,iptpt++,fptpt++;
		} 
	}
	for( j=0; j<lgth2; j++ )
	{
		match[j] = 0.0;
		for( k=0; cpmxpdn[j][k]>-1; k++ )
			match[j] += scarr[cpmxpdn[j][k]] * cpmxpd[j][k];
	} 
#else
	matchpt = match;
	cpmxpdnpt = cpmxpdn;
	cpmxpdpt = cpmxpd;
	while( lgth2-- )
	{
//		*matchpt = 0.0; // add dakara
		for( k=0; (cpkd=(*cpmxpdnpt)[k])>-1; k++ )
			*matchpt += scarr[cpkd] * (*cpmxpdpt)[k];
		matchpt++;
		cpmxpdnpt++;
		cpmxpdpt++;
	}
#endif
	free( scarr );
}

#if 0
// [seq][alphabet]
static void match_calc( double **n_dynamicmtx, double *match, double **cpmx1, double **cpmx2, int i1, int lgth2, double **doublework, int **intwork, int initialize )
{
	int j, k, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	int count = 0;
	double *matchpt;
	double **cpmxpdpt;
	int **cpmxpdnpt;
	int cpkd;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
	if( initialize )
	{
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[j][l] )
				{
					cpmxpd[j][count] = cpmx2[j][l];
					cpmxpdn[j][count] = l;
					count++;
				}
			}
			cpmxpdn[j][count] = -1;
		}
	}

	for( l=0; l<nalphabets; l++ )
	{
		scarr[l] = 0.0;
		for( k=0; k<nalphabets; k++ )
		{
//			scarr[l] += n_dis[k][l] * cpmx1[i1][k];
			scarr[l] += n_dynamicmtx[k][l] * cpmx1[i1][k];
		}
	}
#if 0 /* これを使うときはdoubleworkのアロケートを逆にする */
	{
		double *fpt, **fptpt, *fpt2;
		int *ipt, **iptpt;
		fpt2 = match;
		iptpt = cpmxpdn;
		fptpt = cpmxpd;
		while( lgth2-- )
		{
			*fpt2 = 0.0;
			ipt=*iptpt,fpt=*fptpt;
			while( *ipt > -1 )
				*fpt2 += scarr[*ipt++] * *fpt++;
			fpt2++,iptpt++,fptpt++;
		} 
	}
	for( j=0; j<lgth2; j++ )
	{
		match[j] = 0.0;
		for( k=0; cpmxpdn[j][k]>-1; k++ )
			match[j] += scarr[cpmxpdn[j][k]] * cpmxpd[j][k];
	} 
#else
	matchpt = match;
	cpmxpdnpt = cpmxpdn;
	cpmxpdpt = cpmxpd;
	while( lgth2-- )
	{
		*matchpt = 0.0;
		for( k=0; (cpkd=(*cpmxpdnpt)[k])>-1; k++ )
			*matchpt += scarr[cpkd] * (*cpmxpdpt)[k];
		matchpt++;
		cpmxpdnpt++;
		cpmxpdpt++;
	}
#endif
	free( scarr );
}
#endif


// [alphabet][seq]
static void match_calc_alphabet_seq( double **n_dynamicmtx, double *match, double **cpmx1, double **cpmx2, int i1, int start2, int lgth2, double **doublework, int **intwork, int initialize )
{
#if FASTMATCHCALC
//	fprintf( stderr, "\nmatch_calc... %d", i1 );
	int j, l, p;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	double *matchpt, *cpmxpdpt, **cpmxpdptpt;
	int *cpmxpdnpt, **cpmxpdnptpt;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );

//	reporterr( "lgth2=%d.  j=%d-%d, p=%d-%d\n", lgth2, 0, lgth2, start2, start2+lgth2 );
	if( initialize )
	{
		int count = 0;
		for( j=0,p=start2; j<lgth2; j++,p++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][p] )
				{
					cpmxpd[j][count] = cpmx2[l][p];
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
	int j, k, l, p;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
// simple
	if( initialize )
	{
		int count = 0;
		for( j=0,p=start2; j<lgth2; j++,p++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][p] )
				{
					cpmxpd[count][j] = cpmx2[l][p];
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

static void createcpmxresult( double **cpmxresult, double eff1, double eff2, double **cpmx1, double **cpmx2, char *gaptable1, char *gaptable2 )
{
	int i, j, p;
	int alen = strlen( gaptable1 );

//	reporterr( "eff1 = %f, eff2=%f\n", eff1, eff2 );

#if 1 // sukoshi osoi
	for( i=0; i<nalphabets; i++ )
	{
		for( j=0; j<alen; j++ ) cpmxresult[i][j] = 0.0;
		for( j=0,p=0; j<alen; j++ ) 
		{
			if( gaptable1[j] == '-' )
				;
			else
				cpmxresult[i][j] += cpmx1[i][p++] * eff1;
		}

		for( j=0,p=0; j<alen; j++ ) 
		{
			if( gaptable2[j] == '-' )
				;
			else
				cpmxresult[i][j] += cpmx2[i][p++] * eff2;
		}
	}
#endif
}

static void creategapfreqresult( double *gapfresult, double eff1, double eff2, double *gapf1, double *gapf2, char *gaptable1, char *gaptable2 )
{
	int j, p;
	int alen = strlen( gaptable1 );

#if 1 // sukoshi osoi
	for( j=0; j<alen+1; j++ ) gapfresult[j] = 0.0;
	for( j=0,p=0; j<alen+1; j++ ) 
	{
		if( gaptable1[j] == '-' )
			;
		else
			gapfresult[j] += gapf1[p++] * eff1;
	}

	for( j=0,p=0; j<alen; j++ ) 
	{
		if( gaptable2[j] == '-' )
			;
		else
			gapfresult[j] += gapf2[p++] * eff2;
	}
	gapfresult[j] = 1.0;

#if 0
	reporterr( "gaptable1=%s\n", gaptable1 );
	reporterr( "gaptable2=%s\n", gaptable2 );
	reporterr( "result (gapfreq) = " );
	for( j=0; j<alen+1; j++ ) reporterr( "%4.2f ", gapfresult[j] );
	reporterr( "\n" );
#endif
#endif
}

static void createogresult( double *gapfresult, double eff1, double eff2, double *ori1, double *ori2, double *gf1, double *gf2, char *gaptable1, char *gaptable2 )
{
	int j, p;
	int alen = strlen( gaptable1 );

#if 1 // sukoshi osoi
	for( j=0; j<alen; j++ ) gapfresult[j] = 0.0;
	for( j=0,p=0; j<alen; j++ ) 
	{
		if( gaptable1[j] == '-' )
		{
			if( j==0 )
			{
//				reporterr( "hit, j=%d, result += %f\n",  eff1 );
				gapfresult[j] += 1.0 * eff1;
			}
			else if ( j && gaptable1[j-1] != '-' )
			{
//				reporterr( "hit, j=%d, p=%d, gf1[p]=%f, result += %f (1)\n", j, p, gf1[p], (gf1[p]) * eff1 );
				gapfresult[j] += (gf1[p-1]) * eff1; // p-1 daijoubu?
			}
		}
		else
		{
			if( j==0 || ( j && gaptable1[j-1] != '-' ) )
			{
//				reporterr( "copying. j=%d, ori1[]=%f\n", j, ori1[p] );
				gapfresult[j] += ori1[p] * eff1;
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
				gapfresult[j] += 1.0 * eff2;
			}
			else if ( j && gaptable2[j-1] != '-' )
			{
//				reporterr( "hit, j=%d, p=%d, gf2[p]=%f, result += %f (2)\n", j, p, gf2[p], (gf2[p]) * eff2 );
				gapfresult[j] += (gf2[p-1]) * eff2; // p-1 daijoubu?
			}
		}
		else
		{
			if( j==0 || ( j && gaptable2[j-1] != '-' ) )
			{
//				reporterr( "copying. j=%d, ori2[]=%f\n", j, ori2[p] );
				gapfresult[j] += ori2[p] * eff2;
			}
			p++;
		}
	}
#if 0
	reporterr( "createogresult (o) = " );
	for( j=0; j<alen+1; j++ ) reporterr( "%4.2f ", gapfresult[j] );
	reporterr( "\n" );
#endif
#endif
}

static void createfgresult( double *gapfresult, double eff1, double eff2, double *ori1, double *ori2, double *gf1, double *gf2, char *gaptable1, char *gaptable2 )
{
	int j, p;
	int alen = strlen( gaptable1 );

	for( j=0; j<alen; j++ ) gapfresult[j] = 0.0;
	for( j=0,p=0; j<alen; j++ ) 
	{
		if( gaptable1[j] == '-' )
		{
			if( j==alen-1 )
			{
				gapfresult[j] += eff1;
			}
			else if( gaptable1[j+1] != '-' )
			{
//				reporterr( "hit, j=%d, p=%d, gf1[p]=%f, result += %f (f1)\n", j, p, gf1[p], (gf1[p]) * eff1 );
				gapfresult[j] += gf1[p] * eff1;
			}
		}
		else
		{
			if( gaptable1[j+1] != '-' )
				gapfresult[j] += ori1[p] * eff1;
			p++;
		}
	}

	for( j=0,p=0; j<alen; j++ ) 
	{
		if( gaptable2[j] == '-' )
		{
			if( j==alen-1 )
			{
				gapfresult[j] += eff2;
			}
			else if( gaptable2[j+1] != '-' )
			{
//				reporterr( "hit, j=%d, p=%d, gf2[p]=%f, result += %f (f2)\n", j, p, gf1[p], (gf1[p]) * eff1 );
				gapfresult[j] += gf2[p] * eff2;
			}
		}
		else
		{
			if( gaptable2[j+1] != '-' )
				gapfresult[j] += ori2[p] * eff2;
			p++;
		}
	}
#if 0
	reporterr( "createogresult (f) = " );
	for( j=0; j<alen+1; j++ ) reporterr( "%4.2f ", gapfresult[j] );
	reporterr( "\n" );
#endif
}


static double Atracking( double *lasthorizontalw, double *lastverticalw,
						char **seq1, char **seq2, 
                        char **mseq1, char **mseq2, 
                        int **ijp, int icyc, int jcyc,
						int ist, int ien, int jst, int jen, 
						int fulllen1, int fulllen2, int tailgp,
						char **gt1, char **gt2 )
{
	int i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k, klim;
	char *gaptable1, *gt1bk;
	char *gaptable2, *gt2bk;
	double wm;
	lgth1 = ien-ist+1;
	lgth2 = jen-jst+1;

	if( gt1 == NULL )
	{
		gt1bk = AllocateCharVec( lgth1+lgth2+1 );
		gt2bk = AllocateCharVec( lgth1+lgth2+1 );
		gaptable1 = gt1bk + lgth1+lgth2;
		gaptable2 = gt2bk + lgth1+lgth2;
	}
	else
	{
		gaptable1 = *gt1 + lgth1+lgth2;
		gaptable2 = *gt2 + lgth1+lgth2;
		gt1bk = *gt1;
		gt2bk = *gt2;
	}
	*gaptable1 = 0;
	*gaptable2 = 0;

#if 0
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif


//	fprintf( stderr, "in Atracking, tailgp=%d, ien=%d, fulllen1=%d, jen=%d, fulllen2=%d\n", tailgp, ien, fulllen1, jen, fulllen2 );

	if( tailgp == 1 )
		;
	else if( ien == fulllen1-1 || jen == fulllen2-1 )
	{
//		fprintf( stderr, "searching lasthorizontalw and lasthorizontalw\n" );
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
#if 0
	else if( jen == fulllen2-1 )
	{
		fprintf( stderr, "searching lastverticalw\n" );
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
	}
	else if( ien == fulllen1-1 )
	{
		fprintf( stderr, "searching lasthorizontalw\n" );
		wm = lasthorizontalw[0];
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
#endif
 
    for( i=0; i<lgth1+1; i++ ) 
    {
        ijp[i][0] = i + 1;
    }
    for( j=0; j<lgth2+1; j++ ) 
    {
        ijp[0][j] = -( j + 1 );
    }



//	if( lgth2 == 1 ) fprintf( stderr, "in Atracking, mseq1 = %s, mseq2 = %s\n", mseq1[0], mseq2[0] );

	iin = lgth1; jin = lgth2;
	klim = lgth1+lgth2;
	for( k=0; k<=klim; k++ ) 
	{
		if( ijp[iin][jin] < 0 ) 
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
		if( iin <= 0 || jin <= 0 ) break;
		*--gaptable1 = 'o';
		*--gaptable2 = 'o';
		k++;
		iin = ifi; jin = jfi;

	}

	klim = gt1bk + lgth1+lgth2 - gaptable1;
//	reporterr( "klim = %d, strlen=%d\n", klim, strlen( gaptable1 ) );
//	klim = strlen( gaptable1 );
	if( strchr( gaptable1, '-' ) )
		for( i=0; i<icyc; i++ ) gapireru( mseq1[i], seq1[i]+ist, gaptable1 );
	else
		for( i=0; i<icyc; i++ ) 
		{
			strncpy( mseq1[i], seq1[i]+ist, klim );
			mseq1[i][klim] = 0;
		}

	if( strchr( gaptable2, '-' ) )
		for( j=0; j<jcyc; j++ ) gapireru( mseq2[j], seq2[j]+jst, gaptable2 );
	else
		for( j=0; j<jcyc; j++ )
		{
			strncpy( mseq2[j], seq2[j]+jst, klim );
			mseq2[j][klim] = 0;
		}



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


//	fprintf( stderr, "in Atracking (owari), mseq1 = %s\n", mseq1[0] );
//	fprintf( stderr, "in Atracking (owari), mseq2 = %s\n", mseq2[0] );
	return( 0.0 );
}

static double MSalignmm_tanni( double **n_dynamicmtx, int icyc, int jcyc, double *eff1, double *eff2, char **seq1, char **seq2, double **cpmx1pt, double **cpmx2pt, int ist, int ien, int jst, int jen, int alloclen, int fulllen1, int fulllen2, char **mseq1, char **mseq2, char *mgt1, char *mgt2, double **gapinfo, int headgp, int tailgp, double headgapfreq1_g, double headgapfreq2_g )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
//	int k;
	register int i, j;
	int ll1, ll2;
	int lasti, lastj;
	double wm = 0.0;   /* int ?????? */
	double g;
	double *currentw, *previousw;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penalty_ex;
#endif
#if 1
	double *wtmp;
	int *ijppt;
	double *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	double mi, *m;
	int **ijp;
	int mpi, *mp;
	double *w1, *w2;
	double *initverticalw;    /* kufuu sureba iranai */
	double *lastverticalw;    /* kufuu sureba iranai */
	int **intwork;
	double **doublework;
	int **intmtx;
	int lgth1, lgth2;
	double *ogcp1;
	double *fgcp1;
	double *ogcp2;
	double *fgcp2;
	double *gapfreq1f;
	double *gapfreq2f;
//	char **aseq1;
//	char **aseq2;
//	char **aseq1bk, **aseq2bk;
	double headgapfreq1;
	double headgapfreq2;
	char *gt1, *gt1bk, *gt2, *gt2bk;


	ogcp1 = gapinfo[0] + ist;
	fgcp1 = gapinfo[1] + ist;
	ogcp2 = gapinfo[2] + jst;
	fgcp2 = gapinfo[3] + jst;
	gapfreq1f = gapinfo[4] + ist;
	gapfreq2f = gapinfo[5] + jst;

	if( ist > 0 ) headgapfreq1 = gapfreq1f[-1];
	else headgapfreq1 = headgapfreq1_g;
	if( jst > 0 ) headgapfreq2 = gapfreq2f[-1];
	else headgapfreq2 = headgapfreq2_g;

#if STOREWM
	char ttt1[10000], ttt2[10000];
#endif


	lgth1 = ien-ist+1;
	lgth2 = jen-jst+1;


#if STOREWM
	strncpy( ttt1, seq1[0]+ist, lgth1 ); ttt1[lgth1] = 0;
	strncpy( ttt2, seq2[0]+jst, lgth2 ); ttt2[lgth2] = 0;

	fprintf( stderr, "in _tanni ist,ien = %d,%d, lgth1=%d\n", ist, ien, lgth1 );
	fprintf( stderr, "in _tanni jst,jen = %d,%d, lgth2=%d\n", jst, jen, lgth2 );
	fprintf( stderr, "ttt1 = %s\n", ttt1 );
	fprintf( stderr, "ttt2 = %s\n", ttt2 );
#endif

#if 0
	fprintf( stderr, "in _tanni ist,ien = %d,%d, fulllen1=%d\n", ist, ien, fulllen1 );
	fprintf( stderr, "in _tanni jst,jen = %d,%d, fulllen2=%d\n", jst, jen, fulllen2 );
	fprintf( stderr, "in _tanni seq1[0] = %-*.*s\n", ien-ist+1, ien-ist+1, seq1[0]+ist );
	fprintf( stderr, "in _tanni seq2[0] = %-*.*s\n", jen-jst+1, jen-jst+1, seq2[0]+jst );
#endif


	ll1 = ( (int)(lgth1) ) + 100;
	ll2 = ( (int)(lgth2) ) + 100;

//	aseq1 = AllocateCharMtx( icyc, 0 );
//	aseq2 = AllocateCharMtx( jcyc, 0 );
//	aseq1bk = AllocateCharMtx( icyc, lgth1+lgth2+100 );
//	aseq2bk = AllocateCharMtx( jcyc, lgth1+lgth2+100 );
//	for( i=0; i<icyc; i++ ) aseq1[i] = aseq1bk[i];
//	for( i=0; i<jcyc; i++ ) aseq2[i] = aseq2bk[i];

	w1 = AllocateFloatVec( ll2+2 );
	w2 = AllocateFloatVec( ll2+2 );

	initverticalw = AllocateFloatVec( ll1+2 );
	lastverticalw = AllocateFloatVec( ll1+2 );

	m = AllocateFloatVec( ll2+2 );
	mp = AllocateIntVec( ll2+2 );

#if FASTMATCHCALC
	doublework = AllocateFloatMtx( MAX( ll1, ll2 )+2, nalphabets+1 ); 
	intwork = AllocateIntMtx( MAX( ll1, ll2 )+2, nalphabets+1 ); 
#else
	doublework = AllocateFloatMtx( nalphabets+1, MAX( ll1, ll2 )+2 );
	intwork = AllocateIntMtx( nalphabets+1, MAX( ll1, ll2 )+2 );
#endif

	intmtx = AllocateIntMtx( ll1+1, ll2+1 );

	ijp = intmtx;

	currentw = w1;
	previousw = w2;

//	match_calc( n_dynamicmtx, initverticalw, cpmx2+jst, cpmx1+ist, 0, lgth1, doublework, intwork, 1 );
	match_calc_alphabet_seq( n_dynamicmtx, initverticalw, cpmx2pt, cpmx1pt, jst, ist, lgth1, doublework, intwork, 1 );

//	match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, 0, lgth2, doublework, intwork, 1 );
	match_calc_alphabet_seq( n_dynamicmtx, currentw, cpmx1pt, cpmx2pt, ist, jst, lgth2, doublework, intwork, 1 );

	if( headgp || ist != 0 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2f[0] ) ;
		}
	}
	if( headgp || jst != 0 )
	{
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1f[0] ) ;
		}
	}

	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2f[j-1]; mp[j] = 0;;
	}

	lastverticalw[0] = currentw[lgth2-1];



	if( tailgp || jen != fulllen2-1 ) lasti = lgth1+1; else lasti = lgth1;
//	if( 1 ) lasti = lgth1+1; else lasti = lgth1;
	for( i=1; i<lasti; i++ )
	{
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

//		match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, i, lgth2, doublework, intwork, 0 );
		match_calc_alphabet_seq( n_dynamicmtx, currentw, cpmx1pt, cpmx2pt, ist+i, jst, lgth2, doublework, intwork, 0 );
		currentw[0] = initverticalw[i];

		mi = previousw[0] + ogcp2[1] * gapfreq1f[i-1];
		mpi = 0;

		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;
//		if( tailgp && jen != fulllen2-1 ) lastj = lgth2+1; else lastj = lgth2;
		lastj = lgth2+1; 
		for( j=1; j<lastj; j++ )
		{
			wm = *prept;
			*ijppt = 0;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
			g = mi + fgcp2[j-1]  * gapfreq1f[i];
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
				*ijppt = -( j - mpi );
			}
			g = *prept + ogcp2[j] * gapfreq1f[i-1];
			if( g >= mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

			g = *mjpt + fgcp1[i-1] * gapfreq2f[j];
#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
				*ijppt = +( i - *mpjpt );
			}
			g = *prept + ogcp1[i] * gapfreq2f[j-1];
			if( g >= *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}
#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif

#if 0
			fprintf( stderr, "%5.0f ", wm );
#endif
			*curpt += wm;


			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
			curpt++;
		}
		lastverticalw[i] = currentw[lgth2-1];
	}

//	fprintf( stderr, "wm = %f\n", wm );

	gt1 = gt1bk = AllocateCharVec( ien-ist+jen-jst+3 );
	gt2 = gt2bk = AllocateCharVec( ien-ist+jen-jst+3 );

	Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, ist, ien, jst, jen, fulllen1, fulllen2, tailgp, &gt1, &gt2 );
	strcpy( mgt1, gt1 );
	strcpy( mgt2, gt2 );


#if 0
	fprintf( stderr, "res after _tanni = %s\n", mseq1[0] );
	fprintf( stderr, "res after _tanni = %s\n", mseq2[0] );
	fprintf( stderr, "gt1 after _tanni = %s\n", gt1 );
	fprintf( stderr, "gt1 after _tanni = %s\n", gt2 );
#endif
	free( gt1bk );
	free( gt2bk );

//	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], aseq1[i] );
//	for( i=0; i<jcyc; i++ ) strcpy( mseq2[i], aseq2[i] );

//	fprintf( stderr, "in _tanni, aseq1 = %s\n", aseq1[0] );
//	fprintf( stderr, "in _tanni, mseq1 = %s\n", mseq1[0] );

	FreeFloatVec( w1 );
	FreeFloatVec( w2 );
	FreeFloatVec( initverticalw );
	FreeFloatVec( lastverticalw );

	FreeFloatVec( m );
	FreeIntVec( mp );


	FreeFloatMtx( doublework );
	FreeIntMtx( intwork );

	FreeIntMtx( intmtx );


//	FreeCharMtx( aseq1bk );
//	FreeCharMtx( aseq2bk );

//	free( aseq1 );
//	free( aseq2 );

	return( wm );

}

static void freearrays_rec1(
	double *w1, double *w2, double *initverticalw, double *lastverticalw,
	double *midw, double *midm, double *midn,
	int *jumpbacki, int *jumpbackj, int *jumpforwi, int *jumpforwj, int *jumpdummi, int *jumpdummj,
	double *m, int *mp,
	double **doublework, int **intwork
#if STOREWM
	, double **WMMTX, double **WMMTX2
#endif
)
{

	FreeFloatVec( w1 );
	FreeFloatVec( w2 );
	FreeFloatVec( initverticalw );
	FreeFloatVec( lastverticalw );
	FreeFloatVec( midw );
	FreeFloatVec( midm );
	FreeFloatVec( midn );

	FreeIntVec( jumpbacki );
	FreeIntVec( jumpbackj );
	FreeIntVec( jumpforwi );
	FreeIntVec( jumpforwj );
	FreeIntVec( jumpdummi );
	FreeIntVec( jumpdummj );

	FreeFloatVec( m );
	FreeIntVec( mp );

	FreeFloatMtx( doublework );
	FreeIntMtx( intwork );

#if STOREWM
	FreeFloatMtx( WMMTX );
	FreeFloatMtx( WMMTX2 );
#endif
}

static void freearrays_rec2( char *gaps, char **aseq1, char **aseq2 )
{
	free( gaps );
#if MEMSAVE
	free( aseq1 );
	free( aseq2 );
#else
	FreeCharMtx( aseq1 );
	FreeCharMtx( aseq2 );
#endif
}

static double MSalignmm_rec( double **n_dynamicmtx, int icyc, int jcyc, double *eff1, double *eff2, char **seq1, char **seq2, double **cpmx1pt, double **cpmx2pt, int ist, int ien, int jst, int jen, int alloclen, int fulllen1, int fulllen2, char **mseq1, char **mseq2, char *mgt1, char *mgt2, int depth, double **gapinfo, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp, double headgapfreq1_g, double headgapfreq2_g )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
//	int k;
	int alnlen;
	double value = 0.0;
	register int i, j;
	char **aseq1, **aseq2, *agt1, *agt2;
	int ll1, ll2, l, len;
	int lasti, lastj, imid;
	int jmid = 0; // by D.Mathog, a guess
	double wm = 0.0;   /* int ?????? */
	double g;
	double *currentw, *previousw;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penalty_ex;
#endif
	double *wtmp;
//	short *ijppt;
	int *mpjpt;
//	short **ijp;
	int *mp;
	int mpi;
	double *mjpt, *prept, *curpt;
	double mi;
	double *m;
	double *w1, *w2;
//	double *match;
	double *initverticalw;    /* kufuu sureba iranai */
	double *lastverticalw;    /* kufuu sureba iranai */
	int **intwork;
	double **doublework;
//	short **shortmtx;
#if STOREWM
	double **WMMTX;
	double **WMMTX2;
#endif
	double *midw;
	double *midm;
	double *midn;
	int lgth1, lgth2;
	double maxwm;
	int *jumpforwi;
	int *jumpforwj;
	int *jumpbacki;
	int *jumpbackj;
	int *jumpdummi; //muda
	int *jumpdummj = NULL; // by D.Mathog, a guess
	int jumpi, jumpj = 0; // by D.Mathog, a guess
	char *gaps;
	int ijpi, ijpj;
	double *ogcp1;
	double *fgcp1;
	double *ogcp2;
	double *fgcp2;
	double firstm;
	int firstmp;
#if STOREWM
	static TLS char ttt1[50000];
	static TLS char ttt2[50000];
#endif
	double *gapfreq1f;
	double *gapfreq2f;
	double headgapfreq1;
	double headgapfreq2;

#if 0
	int nglen1, nglen2;
	nglen1 = seqlen( seq1[0] );
	nglen2 = seqlen( seq2[0] );
#endif

//	fprintf( stderr, "fulllen1 = %d, fulllen2 = %d, headgp = %d, tailgp = %d\n", fulllen1, fulllen2, headgp, tailgp );

	ogcp1 = gapinfo[0] + ist;
	fgcp1 = gapinfo[1] + ist;
	ogcp2 = gapinfo[2] + jst;
	fgcp2 = gapinfo[3] + jst;
	gapfreq1f = gapinfo[4] + ist;
	gapfreq2f = gapinfo[5] + jst;

	if( ist > 0 ) headgapfreq1 = gapfreq1f[-1];
	else headgapfreq1 = headgapfreq1_g;
	if( jst > 0 ) headgapfreq2 = gapfreq2f[-1];
	else headgapfreq2 = headgapfreq2_g;

	depth++;
	reccycle++;

	lgth1 = ien-ist+1;
	lgth2 = jen-jst+1;

//	if( lgth1 < 5 )
//		fprintf( stderr, "\nWARNING: lgth1 = %d\n", lgth1 );
//	if( lgth2 < 5 )
//		fprintf( stderr, "\nWARNING: lgth2 = %d\n", lgth2 );
//


#if STOREWM
	fprintf( stderr, "==== MSalign (depth=%d, reccycle=%d), ist=%d, ien=%d, jst=%d, jen=%d\n", depth, reccycle, ist, ien, jst, jen );
	strncpy( ttt1, seq1[0]+ist, lgth1 );
	strncpy( ttt2, seq2[0]+jst, lgth2 );
	ttt1[lgth1] = 0;
	ttt2[lgth2] = 0;
	fprintf( stderr, "seq1 = %s\n", ttt1 );
	fprintf( stderr, "seq2 = %s\n", ttt2 );
#endif
	if( lgth2 <= 0 ) // lgth1 <= 0 ha?
	{
//		fprintf( stderr, "\n\n==== jimei\n\n" );
//		exit( 1 );
		for( i=0; i<icyc; i++ ) 
		{
			strncpy( mseq1[i], seq1[i]+ist, lgth1 );
			mseq1[i][lgth1] = 0;
		}
		for( j=0; j<lgth1; j++ ) mseq2[0][j] = *newgapstr;
		mseq2[0][lgth1] = 0;
		for( i=1; i<jcyc; i++ ) strcpy( mseq2[i], mseq2[0] );

		for( j=0; j<lgth1; j++ ) mgt1[j] = 'o';
		mgt1[lgth1] = 0;
		for( j=0; j<lgth1; j++ ) mgt2[j] = '-'; // *newgapstr deha nai
		mgt2[lgth1] = 0;

#if 0
		fprintf( stderr, "==== mseq1[0] (%d) = %s\n", depth, mseq1[0] );
		fprintf( stderr, "==== mseq2[0] (%d) = %s\n", depth, mseq2[0] );
		fprintf( stderr, "==== mgt1     (%d) = %s\n", depth, mgt1     );
		fprintf( stderr, "==== mgt2     (%d) = %s\n", depth, mgt2     );
#endif

		return( 0.0 );
	}

#if MEMSAVE
	aseq1 = AllocateCharMtx( icyc, 0 );
	aseq2 = AllocateCharMtx( jcyc, 0 );
	for( i=0; i<icyc; i++ ) aseq1[i] = mseq1[i];
	for( i=0; i<jcyc; i++ ) aseq2[i] = mseq2[i];
	agt1 = mgt1;
	agt2 = mgt2;
#else
	aseq1 = AllocateCharMtx( icyc, lgth1+lgth2+100 );
	aseq2 = AllocateCharMtx( jcyc, lgth1+lgth2+100 );
#endif

//	fprintf( stderr, "####(s) seq1[0] (%d) = \n%-*.*s\n (a%d-a%d)\n", depth, ien-ist+1, ien-ist+1, seq1[0]+ist, ist, ien );
//	fprintf( stderr, "####(s) seq2[0] (%d) = \n%-*.*s\n (b%d-b%d)\n", depth, jen-jst+1, jen-jst+1, seq2[0]+jst, jst, jen );

//  if( lgth1 < DPTANNI && lgth2 < DPTANNI ) // & dato lgth ==1 no kanousei ga arunode yokunai 
//    if( lgth1 < DPTANNI ) // kore mo lgth2 ga mijikasugiru kanousei ari
    if( lgth1 < DPTANNI || lgth2 < DPTANNI ) // zettai ni anzen ka?
	{
//		fprintf( stderr, "==== Going to _tanni\n" );

		value = MSalignmm_tanni( n_dynamicmtx, icyc, jcyc, eff1, eff2, seq1, seq2, cpmx1pt, cpmx2pt, ist, ien, jst, jen, alloclen, fulllen1, fulllen2, aseq1, aseq2, agt1, agt2, gapinfo, headgp, tailgp, headgapfreq1_g, headgapfreq2_g );	


#if MEMSAVE
		free( aseq1 );
		free( aseq2 );
#else
		for( i=0; i<icyc; i++ ) strcpy( mseq1[i], aseq1[i] );
		for( i=0; i<jcyc; i++ ) strcpy( mseq2[i], aseq2[i] );

		FreeCharMtx( aseq1 );
		FreeCharMtx( aseq2 );
#endif

//		fprintf( stderr, "value = %f\n", value );

		return( value );
	}
//	fprintf( stderr, "Trying to divide the mtx\n" );

	ll1 = ( (int)(lgth1) ) + 100;
	ll2 = ( (int)(lgth2) ) + 100;

//	fprintf( stderr, "ll1,ll2=%d,%d\n", ll1, ll2 );

	w1 = AllocateFloatVec( ll2+2 );
	w2 = AllocateFloatVec( ll2+2 );
//	match = AllocateFloatVec( ll2+2 );
	midw = AllocateFloatVec( ll2+2 );
	midn = AllocateFloatVec( ll2+2 );
	midm = AllocateFloatVec( ll2+2 );
	jumpbacki = AllocateIntVec( ll2+2 );
	jumpbackj = AllocateIntVec( ll2+2 );
	jumpforwi = AllocateIntVec( ll2+2 );
	jumpforwj = AllocateIntVec( ll2+2 );
	jumpdummi = AllocateIntVec( ll2+2 );
	jumpdummj = AllocateIntVec( ll2+2 );

	initverticalw = AllocateFloatVec( ll1+2 );
	lastverticalw = AllocateFloatVec( ll1+2 );

	m = AllocateFloatVec( ll2+2 );
	mp = AllocateIntVec( ll2+2 );
	gaps = AllocateCharVec( MAX( ll1, ll2 ) + 2 );


#if FASTMATCHCALC
	doublework = AllocateFloatMtx( MAX( ll1, ll2 )+2, nalphabets ); 
	intwork = AllocateIntMtx( MAX( ll1, ll2 )+2, nalphabets ); 
#else
	doublework = AllocateFloatMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
	intwork = AllocateIntMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
#endif

#if DEBUG
	fprintf( stderr, "succeeded\n" );
#endif

#if STOREWM
	WMMTX = AllocateFloatMtx( ll1, ll2 );
	WMMTX2 = AllocateFloatMtx( ll1, ll2 );
#endif
#if 0
	shortmtx = AllocateShortMtx( ll1, ll2 );

#if DEBUG
	fprintf( stderr, "succeeded\n\n" );
#endif

	ijp = shortmtx;
#endif

	currentw = w1;
	previousw = w2;

//	match_calc( n_dynamicmtx, initverticalw, cpmx2+jst, cpmx1+ist, 0, lgth1, doublework, intwork, 1 );
	match_calc_alphabet_seq( n_dynamicmtx, initverticalw, cpmx2pt, cpmx1pt, jst, ist, lgth1, doublework, intwork, 1 );

//	match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, 0, lgth2, doublework, intwork, 1 );
	match_calc_alphabet_seq( n_dynamicmtx, currentw, cpmx1pt, cpmx2pt, ist, jst, lgth2, doublework, intwork, 1 );

	for( i=1; i<lgth1+1; i++ )
	{
//		initverticalw[i] += ( ogcp1[0] + fgcp1[i-1] ) ;
		initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2f[0] ) ;
	}
	for( j=1; j<lgth2+1; j++ )
	{
		currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1f[0]) ;
	}

#if STOREWM
	WMMTX[0][0] = initverticalw[0];
	for( i=1; i<lgth1+1; i++ )
	{
		WMMTX[i][0] = initverticalw[i];
	}
	for( j=1; j<lgth2+1; j++ )
	{
		WMMTX[0][j] = currentw[j];
	}
#endif


	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2f[j-1];
//		m[j] = currentw[j-1];
		mp[j] = 0;
	}

	lastverticalw[0] = currentw[lgth2-1];

	imid = lgth1 * 0.5;

	jumpi = 0; // atode kawaru.
	lasti = lgth1+1;
#if STOREWM
	for( i=1; i<lasti; i++ )
#else
	for( i=1; i<=imid; i++ )
#endif
	{
#ifdef enablemultithread
//		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
		if( chudanpt && *chudanpt != chudanref ) 
		{
//			fprintf( stderr, "\n\n## CHUUDAN!!! zenhan\n" );
			*chudanres = 1;
			freearrays_rec1
			(
				w1, w2, initverticalw, lastverticalw, midw, midm, midn,
				jumpbacki, jumpbackj, jumpforwi, jumpforwj, jumpdummi, jumpdummj,
				m, mp,
				doublework, intwork
#if STOREWM
				, WMMTX, WMMTX2
#endif
			);
			freearrays_rec2( gaps, aseq1, aseq2 );
			return( -1.0 );
		}
#endif
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

//		match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, i, lgth2, doublework, intwork, 0 );
		match_calc_alphabet_seq( n_dynamicmtx, currentw, cpmx1pt, cpmx2pt, ist+i, jst, lgth2, doublework, intwork, 0 );
		currentw[0] = initverticalw[i];

		m[0] = ogcp1[i];
#if STOREM
		WMMTX2[i][0] = m[0];
#endif
		if( i == imid ) midm[0] = m[0];

		mi = previousw[0] + ogcp2[1] * gapfreq1f[i-1]; 
//		mi = previousw[0];
		mpi = 0;


//		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;


		lastj = lgth2+1;
		for( j=1; j<lastj; j++ )
		{

			wm = *prept;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
			g = mi + fgcp2[j-1] * gapfreq1f[i];
//			g = mi + fpenalty;
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
//				*ijppt = -( j - mpi );
			}
			g = *prept + ogcp2[j] * gapfreq1f[i-1];
//			g = *prept;
			if( g >= mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

			g = *mjpt + fgcp1[i-1] * gapfreq2f[j];
//			g = *mjpt + fpenalty;
#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
//				*ijppt = +( i - *mpjpt );
			}


			g = *prept + ogcp1[i] * gapfreq2f[j-1];
//			g = *prept;
			if( g >= *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}
#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif

#if 0
			fprintf( stderr, "%5.0f ", wm );
#endif
			*curpt += wm;

#if STOREWM
			WMMTX[i][j] = *curpt;
			WMMTX2[i][j] = *mjpt;
#endif

			if( i == imid ) //muda
			{	
				jumpbackj[j] = *mpjpt; // muda atode matomeru
				jumpbacki[j] = mpi; // muda atode matomeru
//				fprintf( stderr, "jumpbackj[%d] in forward dp is %d\n", j, *mpjpt );
//				fprintf( stderr, "jumpbacki[%d] in forward dp is %d\n", j, mpi );
				midw[j] = *curpt;
				midm[j] = *mjpt;
				midn[j] = mi;
			}

//			fprintf( stderr, "m[%d] = %f\n", j, m[j] );
			mjpt++;
			prept++;
			mpjpt++;
			curpt++;

		}
		lastverticalw[i] = currentw[lgth2-1];

#if STOREWM
		WMMTX2[i][lgth2] = m[lgth2-1];
#endif

#if 0  // ue
		if( i == imid )
		{
			for( j=0; j<lgth2; j++ ) midw[j] = currentw[j];
			for( j=0; j<lgth2; j++ ) midm[j] = m[j];
		}
#endif
	}
//	for( j=0; j<lgth2; j++ ) midw[j] = WMMTX[imid][j];
//	for( j=0; j<lgth2; j++ ) midm[j] = WMMTX2[imid][j];

#if 0
    for( i=0; i<lgth1; i++ )
    {
        for( j=0; j<lgth2; j++ )
        {
            fprintf( stderr, "% 10.2f ", WMMTX[i][j] );
        }
        fprintf( stderr, "\n" );
    }
	fprintf( stderr, "\n" );
	fprintf( stderr, "WMMTX2 = \n" );
    for( i=0; i<lgth1; i++ )
    {
        for( j=0; j<lgth2; j++ )
        {
            fprintf( stderr, "% 10.2f ", WMMTX2[i][j] );
        }
        fprintf( stderr, "\n" );
    }
	fprintf( stderr, "\n" );
#endif

// gyakudp

//	match_calc( n_dynamicmtx, initverticalw, cpmx2+jst, cpmx1+ist, lgth2-1, lgth1, doublework, intwork, 1 );
	match_calc_alphabet_seq( n_dynamicmtx, initverticalw, cpmx2pt, cpmx1pt, jst+lgth2-1, ist, lgth1, doublework, intwork, 1 );

//	match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, lgth1-1, lgth2, doublework, intwork, 1 );
	match_calc_alphabet_seq( n_dynamicmtx, currentw, cpmx1pt, cpmx2pt, ist+lgth1-1, jst, lgth2, doublework, intwork, 1 );

	for( i=0; i<lgth1-1; i++ )
	{
//		initverticalw[i] += ( fgcp1[lgth1-1] + ogcp1[i+1] );
		initverticalw[i] += ( fgcp1[lgth1-1] * gapfreq2f[lgth2] + ogcp1[i+1] * gapfreq2f[lgth2-1] );
	}
	for( j=0; j<lgth2-1; j++ )
	{
//		currentw[j] += ( fgcp2[lgth2-1] + ogcp2[j+1] );
		currentw[j] += ( fgcp2[lgth2-1] * gapfreq1f[lgth1] + ogcp2[j+1] * gapfreq1f[lgth1-1] );
	}

#if STOREWM
	for( i=0; i<lgth1-1; i++ )
	{
		WMMTX[i][lgth2-1] += ( fgcp1[lgth1-1] + ogcp1[i+1] );
		fprintf( stderr, "fgcp1[lgth1-1] + ogcp1[i+1] = %f\n", fgcp1[lgth1-1] + ogcp1[i+1] );
	}
	for( j=0; j<lgth2-1; j++ )
	{
		WMMTX[lgth1-1][j] += ( fgcp2[lgth2-1] + ogcp2[j+1] );
		fprintf( stderr, "fgcp2[lgth2-1] + ogcp2[j+1] = %f\n", fgcp2[lgth2-1] + ogcp2[j+1] );
	}
#endif






#if 0
	for( j=lgth2-1; j>0; --j )
	{
		m[j-1] = currentw[j] + fgcp2[lgth2-2];
//		m[j-1] = currentw[j];
		mp[j] = lgth1-1;
	}
#else
	for( j=lgth2-1; j>-1; --j )
	{
		m[j] = currentw[j+1] + fgcp1[lgth1-2] * gapfreq2f[j+1];
//		m[j-1] = currentw[j];
		mp[j] = lgth1-1;
	}
#endif

//	for( j=0; j<lgth2; j++ ) m[j] = 0.0;
	// m[lgth2-1] ha irunoka?


//	for( i=lgth1-2; i>=imid; i-- )
	firstm = -9999999.9;
//	firstmp = lgth1-1;
	firstmp = lgth1;
	for( i=lgth1-2; i>-1; i-- )
	{
#ifdef enablemultithread
//		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
		if( chudanpt && *chudanpt != chudanref ) 
		{
//			fprintf( stderr, "\n\n## CHUUDAN!!! kouhan\n" );
			*chudanres = 1;
			freearrays_rec1
			(
				w1, w2, initverticalw, lastverticalw, midw, midm, midn,
				jumpbacki, jumpbackj, jumpforwi, jumpforwj, jumpdummi, jumpdummj,
				m, mp,
				doublework, intwork
#if STOREWM
				, WMMTX, WMMTX2
#endif
			);
			freearrays_rec2( gaps, aseq1, aseq2 );
			return( -1.0 );
		}
#endif
		wtmp = previousw;
		previousw = currentw;
		currentw = wtmp;
		previousw[lgth2-1] = initverticalw[i+1];
//		match_calc( currentw, seq1, seq2, i, lgth2 );
//		match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, i, lgth2, doublework, intwork, 0 );
		match_calc_alphabet_seq( n_dynamicmtx, currentw, cpmx1pt, cpmx2pt, ist+i, jst, lgth2, doublework, intwork, 0 );

		currentw[lgth2-1] = initverticalw[i];

//		m[lgth2] = fgcp1[i];
//		WMMTX2[i][lgth2] += m[lgth2];
//		fprintf( stderr, "m[] = %f\n", m[lgth2] );

		mi = previousw[lgth2-1] + fgcp2[lgth2-2] * gapfreq1f[i+1];
//		mi = previousw[lgth2-1];
		mpi = lgth2 - 1;

		mjpt = m + lgth2 - 2;
		prept = previousw + lgth2 - 1;
		curpt = currentw + lgth2 - 2;
		mpjpt = mp + lgth2 - 2;


		for( j=lgth2-2; j>-1; j-- )
		{
			wm = *prept;
			ijpi = i+1;
			ijpj = j+1;

			g = mi + ogcp2[j+1] * gapfreq1f[i];
//			g = mi + fpenalty;
			if( g > wm )
			{
				wm = g;
				ijpj = mpi;
				ijpi = i+1;
			}

			g = *prept + fgcp2[j] * gapfreq1f[i+1];
//			g = *prept;
			if( g >= mi )
			{
//				fprintf( stderr, "i,j=%d,%d - renewed! mpi = %d\n", i, j, j+1 );
				mi = g;
				mpi = j + 1;
			}

#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

//			fprintf( stderr, "i,j=%d,%d *mpjpt = %d\n", i, j, *mpjpt );
			g = *mjpt + ogcp1[i+1] * gapfreq2f[j];
//			g = *mjpt + fpenalty;
			if( g > wm )
			{
				wm = g;
				ijpi = *mpjpt;
				ijpj = j+1;
			}

//			if( i == imid )fprintf( stderr, "i,j=%d,%d \n", i, j );
			g = *prept + fgcp1[i] * gapfreq2f[j+1];
//			g = *prept;
			if( g >= *mjpt )
			{
				*mjpt = g;
				*mpjpt = i + 1;
			}

#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif

			if( i == jumpi || i == imid - 1 )
			{
				jumpforwi[j] = ijpi; //muda
				jumpforwj[j] = ijpj; //muda
//				fprintf( stderr, "jumpfori[%d] = %d\n", j, ijpi );
//				fprintf( stderr, "jumpforj[%d] = %d\n", j, ijpj );
			}
			if( i == imid ) // muda
			{
				midw[j] += wm;
//				midm[j+1] += *mjpt + fpenalty; //??????
				midm[j+1] += *mjpt; //??????
			}
			if( i == imid - 1 )
			{
//				midn[j] += mi + fpenalty;  //????
				midn[j] += mi;  //????
			}
#if STOREWM
			WMMTX[i][j] += wm;
//			WMMTX2[i][j+1] += *mjpt + fpenalty;
			WMMTX2[i][j+1] += *mjpt;
#endif
			*curpt += wm;

			mjpt--;
			prept--;
			mpjpt--;
			curpt--;
		}
//		fprintf( stderr, "adding *mjpt (=%f) to WMMTX2[%d][%d]\n", *mjpt, i, j+1 );
		g = *prept + fgcp1[i];
		if( firstm < g ) 
		{
			firstm = g;
			firstmp = i + 1;
		}
#if STOREWM
		WMMTX2[i][j+1] += firstm;
#endif
		if( i == imid ) midm[j+1] += firstm;


		if( i == imid - 1 )	
		{
			maxwm = midw[1];
			jmid = 0;
//			if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			for( j=2; j<lgth2-1; j++ )
			{
				wm = midw[j];
				if( wm > maxwm )
				{
					jmid = j;
					maxwm = wm;
				}
//				if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			}
			for( j=0; j<lgth2+1; j++ )
			{
				wm = midm[j];
				if( wm > maxwm )
				{
					jmid = j;
					maxwm = wm;
				}
//				if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			}

//			if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );


//			fprintf( stderr, "### imid=%d, jmid=%d\n", imid, jmid );
			wm = midw[jmid];
			jumpi = imid-1;
			jumpj = jmid-1;
			if( jmid > 0 && midn[jmid-1] > wm ) //060413
			{
				jumpi = imid-1;
				jumpj = jumpbacki[jmid];
				wm = midn[jmid-1];
//				fprintf( stderr, "rejump (n)\n" );
			}
			if( midm[jmid] > wm )
			{
				jumpi = jumpbackj[jmid];
				jumpj = jmid-1;
				wm = midm[jmid];
//				fprintf( stderr, "rejump (m) jumpi=%d\n", jumpi );
			}


//			fprintf( stderr, "--> imid=%d, jmid=%d\n", imid, jmid );
//			fprintf( stderr, "--> jumpi=%d, jumpj=%d\n", jumpi, jumpj );
#if STOREWM
			fprintf( stderr, "imid = %d\n", imid );
			fprintf( stderr, "midn = \n" );
			for( j=0; j<lgth2; j++ )
			{
				fprintf( stderr, "% 7.1f ", midn[j] );
			}
			fprintf( stderr, "\n" );
			fprintf( stderr, "midw = \n" );
			for( j=0; j<lgth2; j++ )
			{
				fprintf( stderr, "% 7.1f ", midw[j] );
			}
			fprintf( stderr, "\n" );
			fprintf( stderr, "midm = \n" );
			for( j=0; j<lgth2; j++ )
			{
				fprintf( stderr, "% 7.1f ", midm[j] );
			}
			fprintf( stderr, "\n" );
#endif
//			fprintf( stderr, "maxwm = %f\n", maxwm );
		}
		if( i == jumpi ) //saki?
		{
//			fprintf( stderr, "#### FIRST i=%d, jumpi<imid=%d<%d, ist=%d, ien=%d, firstmp=%d, lgth1=%d\n", i, jumpi, imid, ist, ien, firstmp, lgth1 );
//			fprintf( stderr, "#### mark 1\n" );
//			fprintf( stderr, "imid, jumpi = %d,%d\n", imid, jumpi );
//			fprintf( stderr, "jmid, jumpj = %d,%d\n", jmid, jumpj );

			if( jmid == 0 )
			{
//				fprintf( stderr, "#### CHUI2!\n" );
				jumpj = 0; jmid = 1;
#if 0 // v6.823 made
				jumpi = firstmp-1;
				imid = firstmp;
#endif
#if 0
				jumpi = 0;
				imid = 1;
#else
//				if( 1 || firstmp > 100 ) // naze 100
				if( imid < firstmp-1 ) // naze 100
				{
					jumpi = firstmp;
					imid = firstmp+1;
				}
#if 0
				else
				{
					jumpi = 0;
					imid = 1;
				}
#endif
#endif
			}

#if 0
			else if( jmid == lgth2 ) 
			{
				fprintf( stderr, "CHUI1!\n" );
				jumpi=0; jumpj=0;
				imid=jumpforwi[0]; jmid=lgth2-1;
			}
#else // 060414
			else if( jmid >= lgth2 ) 
			{
//				fprintf( stderr, "CHUI1!\n" );
				jumpi=imid-1; jmid=lgth2;
				jumpj = lgth2-1;
			}
#endif
			else
			{
//				fprintf( stderr, "#### CHUI3!\n" );
				imid = jumpforwi[jumpj];
				jmid = jumpforwj[jumpj];
				if( imid == jumpi ) jumpi = imid-1;
			}
#if 0
			fprintf( stderr, "jumpi -> %d\n", jumpi );
			fprintf( stderr, "jumpj -> %d\n", jumpj );
			fprintf( stderr, "imid -> %d\n", imid );
			fprintf( stderr, "jmid -> %d\n", jmid );
#endif
//			fprintf( stderr, "#### FINAL i=%d, jumpi<imid=%d<%d, ist=%d, ien=%d, firstmp=%d\n", i, jumpi, imid, ist, ien, firstmp );

#if STOREWM
			break;
#else
			break;
#endif
		}
	}

//	fprintf( stderr, "imid = %d, but jumpi = %d\n", imid, jumpi );
//	fprintf( stderr, "jmid = %d, but jumpj = %d\n", jmid, jumpj );

//	for( j=0; j<lgth2; j++ ) midw[j] += currentw[j];
//	for( j=0; j<lgth2; j++ ) midm[j] += m[j+1];
//	for( j=0; j<lgth2; j++ ) midw[j] += WMMTX[imid][j];
//	for( j=0; j<lgth2; j++ ) midw[j] += WMMTX[imid][j];


#if STOREWM
	fprintf( stderr, "WMMTX = \n" );
    for( i=0; i<lgth1; i++ )
    {
        fprintf( stderr, "%d ", i );
        for( j=0; j<lgth2; j++ )
        {
            fprintf( stderr, "% 7.2f ", WMMTX[i][j] );
        }
        fprintf( stderr, "\n" );
    }
//	fprintf( stderr, "WMMTX2 = (p = %f)\n", fpenalty );
    for( i=0; i<lgth1; i++ )
    {
        fprintf( stderr, "%d ", i );
        for( j=0; j<lgth2+1; j++ )
        {
            fprintf( stderr, "% 7.2f ", WMMTX2[i][j] );
        }
        fprintf( stderr, "\n" );
    }

	fprintf( stderr, "jumpbacki = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpbacki[j] );
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "jumpbackj = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpbackj[j] );
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "jumpforwi = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpforwi[j] );
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "jumpforwj = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpforwj[j] );
	}
	fprintf( stderr, "\n" );


#endif


//	Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, cpmx1, cpmx2, ijp, icyc, jcyc );

#if 0 // irukamo
	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}
#endif



#if 0
	fprintf( stderr, "jumpi = %d, imid = %d\n", jumpi, imid );
	fprintf( stderr, "jumpj = %d, jmid = %d\n", jumpj, jmid );

	fprintf( stderr, "imid = %d\n", imid );
	fprintf( stderr, "jmid = %d\n", jmid );
#endif

	freearrays_rec1
	(
		w1, w2, initverticalw, lastverticalw, midw, midm, midn,
		jumpbacki, jumpbackj, jumpforwi, jumpforwj, jumpdummi, jumpdummj,
		m, mp,
		doublework, intwork
#if STOREWM
		, WMMTX, WMMTX2
#endif
	);


//	fprintf( stderr, "==== calling myself (first), depth=%d\n", depth );
#if 0
		fprintf( stderr, "seq1[0] = %.*s\n", lgth1, seq1[0] );
		fprintf( stderr, "seq2[0] = %.*s\n", lgth2, seq2[0] );
#endif

	value = MSalignmm_rec( n_dynamicmtx, icyc, jcyc, eff1, eff2, seq1, seq2, cpmx1pt, cpmx2pt, ist, ist+jumpi, jst, jst+jumpj, alloclen, fulllen1, fulllen2, aseq1, aseq2, agt1, agt2, depth, gapinfo, NULL, 0, NULL, headgp, tailgp, headgapfreq1_g, headgapfreq2_g ); // chudan mada	
#if 0
	reporterr( "length1=%d -> %d? %d?\n", lgth1, strlen(seq1[0]), strlen(aseq1[0]) );
		reporterr( "after first _rec\n" );
		if( strlen( aseq1[0] ) != strlen( agt1 ) ) reporterr( "WARNING\n" );
		fprintf( stderr, "aseq1[0] = %s\n", aseq1[0] );
		fprintf( stderr, "aseq2[0] = %s\n", aseq2[0] );
		fprintf( stderr, "agt1     = %s\n", agt1     );
		fprintf( stderr, "agt2     = %s\n", agt2     );
#endif
#if MEMSAVE
#else
	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], aseq1[i] );
	for( i=0; i<jcyc; i++ ) strcpy( mseq2[i], aseq2[i] );
#endif

//	fprintf( stderr, "====(f) aseq1[0] (%d) = %s (%d-%d)\n", depth, aseq1[0], ist, ien );
//	fprintf( stderr, "====(f) aseq2[0] (%d) = %s (%d-%d)\n", depth, aseq2[0], jst, jen );

	len = strlen( mseq1[0] );
//	fprintf( stderr, "len = %d\n", len );
	l = jmid - jumpj - 1;
	if( l > 0 )
	{
//		for( i=0; i<l; i++ ) gaps[i] = '-'; gaps[i] = 0;
		for( i=0; i<l; i++ ) gaps[i] = *newgapstr; gaps[i] = 0;
		for( i=0; i<icyc; i++ ) 
		{
			strcat( mseq1[i], gaps );
			mseq1[i][len+l] = 0;
		}
		for( j=0; j<jcyc; j++ )
		{
			strncat( mseq2[j], seq2[j]+jst+jumpj+1, l );
			mseq2[j][len+l] = 0;
		}
		for( i=0,j=len; i<l; i++,j++ ) mgt1[j] = '-'; mgt1[j] = 0;
		for( i=0,j=len; i<l; i++,j++ ) mgt2[j] = 'o'; mgt2[j] = 0;
//		fprintf( stderr, "penalizing (2) .. %f(%d), %f(%d)\n", ogcp2[jumpj+1], jumpj+1, fgcp2[jmid-1], jmid-1 );
		value +=  ( ogcp2[jumpj+1] + fgcp2[jmid-1] );
//		value += fpenalty;
	}
	len = strlen( mseq1[0] );
	l = imid - jumpi - 1;
//	fprintf( stderr, "l=%d\n", l );
	if( l > 0 )
	{
//		for( i=0; i<l; i++ ) gaps[i] = '-'; gaps[i] = 0;
		for( i=0; i<l; i++ ) gaps[i] = *newgapstr; gaps[i] = 0;
		for( i=0; i<icyc; i++ )
		{
			strncat( mseq1[i], seq1[i]+ist+jumpi+1, l );
			mseq1[i][len+l] = 0;
		}
		for( j=0; j<jcyc; j++ ) 
		{
			strcat( mseq2[j], gaps );
			mseq2[j][len+l] = 0;
		}
		for( i=0,j=len; i<l; i++,j++ ) mgt1[j] = 'o'; mgt1[j] = 0;
		for( i=0,j=len; i<l; i++,j++ ) mgt2[j] = '-'; mgt2[j] = 0;

//		for( i=0; i<lgth1; i++ ) fprintf( stderr, "ogcp1[%d] = %f\n", i, ogcp1[i] );
//		for( i=0; i<lgth1; i++ ) fprintf( stderr, "fgcp1[%d] = %f\n", i, fgcp1[i] );


//		fprintf( stderr, "penalizing (1) .. ogcp1[%d] = %f, fgcp1[%d] = %f\n", jumpi+1, ogcp1[jumpi+1], imid-1, fgcp1[imid-1] );
		value += ( ogcp1[jumpi+1] + fgcp1[imid-1] );
//		value += fpenalty;
	}
#if 0
	for( i=0; i<icyc; i++ ) fprintf( stderr, "after gapfill mseq1[%d]=%s\n", i, mseq1[i] );
	for( i=0; i<jcyc; i++ ) fprintf( stderr, "after gapfill mseq2[%d]=%s\n", i, mseq2[i] );
#endif

//	fprintf( stderr, "==== calling myself (second), depth=%d\n", depth );

#if MEMSAVE
	alnlen = strlen( aseq1[0] );
	for( i=0; i<icyc; i++ ) aseq1[i] += alnlen;
	for( i=0; i<jcyc; i++ ) aseq2[i] += alnlen;
	agt1 += alnlen;
	agt2 += alnlen;
#endif

#if 0
		fprintf( stderr, "seq1[0] = %.*s\n", lgth1, seq1[0] );
		fprintf( stderr, "seq2[0] = %.*s\n", lgth2, seq2[0] );
#endif
	value += MSalignmm_rec( n_dynamicmtx, icyc, jcyc, eff1, eff2, seq1, seq2, cpmx1pt, cpmx2pt, ist+imid, ien, jst+jmid, jen, alloclen, fulllen1, fulllen2, aseq1, aseq2, agt1, agt2, depth, gapinfo, NULL, 0, NULL, headgp, tailgp, headgapfreq1_g, headgapfreq2_g ); // chudan mada
#if 0
		reporterr( "after second _rec\n" );
		if( strlen( aseq1[0] ) != strlen( agt1 ) ) reporterr( "WARNING\n" );
		fprintf( stderr, "aseq1[0] = %s\n", aseq1[0] );
		fprintf( stderr, "aseq2[0] = %s\n", aseq2[0] );
		fprintf( stderr, "agt1     = %s\n", agt1     );
		fprintf( stderr, "agt2     = %s\n", agt2     );
#endif



#if DEBUG
	if( value - maxwm > 1 || maxwm - value > 1 )
	{
		fprintf( stderr, "WARNING value  = %f, but maxwm = %f\n", value, maxwm );
		for( i=0; i<icyc; i++ )
		{
			fprintf( stderr, ">1-%d\n%s\n", i, mseq1[i] );
			fprintf( stderr, "%s\n", aseq1[i] );
		}
		for( i=0; i<jcyc; i++ )
		{
			fprintf( stderr, ">2-%d\n%s\n", i, mseq2[i] );
			fprintf( stderr, "%s\n", aseq2[i] );
		}

//		exit( 1 );
	}
	else
	{
		fprintf( stderr, "value = %.0f, maxwm = %.0f -> ok\n", value, maxwm );
	}
#endif

#if MEMSAVE
#else
	for( i=0; i<icyc; i++ ) strcat( mseq1[i], aseq1[i] );
	for( i=0; i<jcyc; i++ ) strcat( mseq2[i], aseq2[i] );
#endif


//	fprintf( stderr, "====(s) aseq1[0] (%d) = \n%s\n (a%d-a%d)\n", depth, mseq1[0], ist, ien );
//	fprintf( stderr, "====(s) aseq2[0] (%d) = \n%s\n (b%d-b%d)\n", depth, mseq2[0], jst, jen );
//	reporterr( "agt1     = %s\n", agt1     );

#if 0
	reporterr( "At the end of _rec,\n" );
	reporterr( "aseq1[0] = %s\n", aseq1[0] );
	reporterr( "agt1     = %s\n", agt1     );
	if( strlen( aseq1[0] ) != strlen (agt1) )
		reporterr( "Dame\n" );
#endif

	freearrays_rec2( gaps, aseq1, aseq2 );

#if 0
	if( seqlen( seq1[0] ) != nglen1 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen1) seqlen(seq1[0])=%d but nglen1=%d\n", seqlen( seq1[0] ), nglen1 );
		fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
		exit( 1 );
	}
	else
		fprintf( stderr, "nglen1 is ok in _rec\n" );
	if( seqlen( seq2[0] ) != nglen2 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen2) seqlen(seq2[0])=%d but nglen2=%d\n", seqlen( seq2[0] ), nglen2 );
		exit( 1 );
	}
	else
		fprintf( stderr, "nglen2 is ok in _rec\n" );
#endif



	return( value );
}

/* =============================================================================================
   work/v10w: the memsave aligner, fast and exact.

   Once an alignment is wider than 30000 columns MAFFT switches every later progressive step to
   MSalignmm (directly, or per segment from Falign_udpari_long): a Myers-Miller / Hirschberg
   divide and conquer.  At each node of the recursion a forward sweep (rows 1..imid) and a
   backward sweep (rows lgth1-2 down to imid-1, or further down to jumpi) find the split
   (imid, jmid, jumpi, jumpj); the two halves are solved recursively, with up to two filler runs
   of gaps between them; a node with fewer than 100 rows or columns is solved by a full DP with a
   traceback (MSalignmm_tanni).  The result is a gap table per group ('o' residue, '-' gap) and the
   rows, which are the input rows with gaps inserted where the table says '-'.

   What is new here, and why every value is the same:
   1. Match scores.  The original forms each row of match scores with match_calc_alphabet_seq: a
      calloc'ed 26-entry table scarr[l] = sum over all 26 letters k (ascending) of
      n[k][l] * cpmx1[k][row], then per column the sum over the column's nonzero letters
      (ascending) of scarr[l] * cpmx2[l][col], every term one "a += x * y" statement from 0.0.
      Here scarr is computed once per row of profile 1 for the whole call (table S, only for the
      letters of profile 2) instead of once per row per sweep per recursion level, and the column
      sums run over a fixed list of letters.  The terms left out or added are exact no-ops: a term
      x * 0.0 or 0.0 * y is +-0.0, and adding +-0.0 to a sum that is not -0.0 returns the sum
      bit for bit; these sums start at +0.0 and never become -0.0 (an exact cancellation rounds to
      +0.0, and +0.0 + -0.0 is +0.0), with or without contraction to fma.  Letters whose S entry
      is 0.0 in every row (the gap letter, which scores 0 against everything) are dropped; the
      others keep their ascending order.  The same holds for the first row and column of each
      sweep (table T, the transposed product).
   2. The DP cells.  Every cell computes the same four candidates with the same statements
      (each "a + b * c" one statement, so a compiler that contracts to fma does so in both
      versions), compares them in the same order with the same strictness (> or >=), and stores
      the same values; the selects are branch-free (MSF_DSEL_GE / MSF_ISEL, as gccpar's
      DSEL_GE / ISEL in Salignmm.c).  Rows other than imid (forward) and imid, imid-1, jumpi
      (backward) need neither the row carry's column (mpi) nor the traceback-like indices, so the
      common rows skip them; the special rows run the original statements.  The stale value in
      column lgth2 of the forward rows (match_calc fills columns 0..lgth2-1 only; the loop goes
      to lgth2) evolves exactly as before in the same two zero-initialized buffers, because
      midw[lgth2] can be read in the jmid search.
   3. Rows.  Each recursion node of the original writes its piece of every row (gapireru) as it
      goes; here the recursion produces only the two gap tables, and the rows are made once at
      the end with the same gap character (the rows are exactly gapireru( input row, table )).
   4. Parallelism (with idle helper threads): the backward sweep's rows below imid+1 do not
      depend on the forward sweep, so the two sweeps of a large node run at the same time (the
      backward one waits for the forward one before row imid); the two halves of a node are
      independent subproblems, solved by any thread, their gap tables joined in order; the rows
      are made in parallel.  Each value is still computed by the same statements from the same
      operands; only the thread differs.  The score that MSalignmm returns (unused by the
      progressive pass, but kept) is summed in the original order: left half, fillers, right half.

   MAFFT_MSCHECK=1 runs the original MSalignmm on copies of the inputs for every call and compares
   rows, gap tables, cached profile and score bitwise (exit on a mismatch); MAFFT_MSFAST=0 turns
   the new code off. */
#define MSF_MAXL 32
#if defined( __x86_64__ ) || defined( __i386__ )
#define MSF_DSEL_GE( C, X ) ( (X) > (C) ? (X) : (C) )
#else
#define MSF_DSEL_GE( C, X ) ( (C) >= (X) ? (C) : (X) )
#endif
#define MSF_ISEL( C, A, B ) ( (B) ^ ( ( (B) ^ (A) ) & -(C) ) )

typedef struct
{
	double **n_dynamicmtx;
	double **cpmx1, **cpmx2;   // [letter][column]
	double **gapinfo;          // ogcp1, fgcp1, ogcp2, fgcp2, gapfreq1f, gapfreq2f (whole profiles)
	int fulllen1, fulllen2, headgp, tailgp;
	double hgf1g, hgf2g;       // headgapfreq1_g, headgapfreq2_g
	int nl1, k1[MSF_MAXL];     // letters with a nonzero weight in profile 1 (columns 0..fulllen1), ascending
	int nl2, k2[MSF_MAXL];     // the same for profile 2 (columns 0..fulllen2)
	int ns, su[MSF_MAXL];      // slots u (into k2) whose S entry is nonzero in some row, ascending
	double *S;                 // S[r*nl2+u], r = 0..fulllen1: scarr of row r for letter k2[u]
	double *T;                 // T[p*nl1+t], p = 0..fulllen2: the transposed scarr
	const double *c4[4];       // ns <= 4: profile-2 rows of the slots su[], padded with a zero row
	int s4[4];                 // their slots (padding: -1)
	double *zero;              // fulllen2+2 zeros
	int general;               // ns > 4: match rows from per-column lists
	int *cst;                  // per column p: entries cst[p] .. cst[p+1]-1 of cu/cw
	int *cu;                   // slot u of each entry (nonzero weight, letter in su[]), ascending
	double *cw;                // its weight
	char gapchar;
} msf_t;

static int msf_on( void )
{
	static int v = -1;
	if( v < 0 ) { char *e = getenv( "MAFFT_MSFAST" ); v = !( e && atoi( e ) == 0 ); }
	return( v );
}
static int mscheck_v = -1;
static int mscheck( void )
{
	if( mscheck_v < 0 ) mscheck_v = ( getenv( "MAFFT_MSCHECK" ) != NULL );
	return( mscheck_v );
}
static long mscheck_n[2];
static void mscheck_report( void ) __attribute__(( destructor ));
static void mscheck_report( void )
{
	if( mscheck_v > 0 ) reporterr( "\nMSCHECK: %ld MSalignmm calls (%ld rows) verified against the original, no mismatch\n", mscheck_n[0], mscheck_n[1] );
}

static void *msf_alloc( size_t n )
{
	void *p = malloc( n ? n : 1 );
	if( p == NULL ) ErrorExit( "Cannot allocate in MSalignmm (work/v10w).\n" );
	return( p );
}
static void *msf_calloc( size_t n, size_t s )
{
	void *p = calloc( n ? n : 1, s ? s : 1 );
	if( p == NULL ) ErrorExit( "Cannot allocate in MSalignmm (work/v10w).\n" );
	return( p );
}

/* returns 0 if the call is not supported (then the original runs) */
static int msf_init( msf_t *c, double **n_dynamicmtx, double **cpmx1, double **cpmx2, double **gapinfo, int lgth1, int lgth2, int headgp, int tailgp, double hgf1, double hgf2 )
{
	int r, p, t, u, k, seen[MSF_MAXL];
	if( nalphabets > MSF_MAXL || nalphabets <= 0 ) return( 0 );
	memset( c, 0, sizeof( msf_t ) );
	c->n_dynamicmtx = n_dynamicmtx; c->cpmx1 = cpmx1; c->cpmx2 = cpmx2; c->gapinfo = gapinfo;
	c->fulllen1 = lgth1; c->fulllen2 = lgth2; c->headgp = headgp; c->tailgp = tailgp;
	c->hgf1g = hgf1; c->hgf2g = hgf2;
	c->gapchar = *newgapstr;

	for( k=0; k<nalphabets; k++ ) seen[k] = 0;
	for( k=0; k<nalphabets; k++ ) for( r=0; r<=lgth1; r++ ) if( cpmx1[k][r] ) { seen[k] = 1; break; }
	c->nl1 = 0; for( k=0; k<nalphabets; k++ ) if( seen[k] ) c->k1[c->nl1++] = k;
	for( k=0; k<nalphabets; k++ ) seen[k] = 0;
	for( k=0; k<nalphabets; k++ ) for( p=0; p<=lgth2; p++ ) if( cpmx2[k][p] ) { seen[k] = 1; break; }
	c->nl2 = 0; for( k=0; k<nalphabets; k++ ) if( seen[k] ) c->k2[c->nl2++] = k;

	/* S: scarr of every row of profile 1 (match_calc_alphabet_seq: scarr[l] = 0.0; scarr[l] += n[j][l] * cpmx1[j][i1], j ascending) */
	c->S = msf_alloc( sizeof( double ) * (size_t)( lgth1 + 1 ) * ( c->nl2 ? c->nl2 : 1 ) );
	for( r=0; r<=lgth1; r++ )
		for( u=0; u<c->nl2; u++ )
		{
			int l = c->k2[u];
			double sc = 0.0;
			for( t=0; t<c->nl1; t++ )
			{
				int j = c->k1[t];
				sc += n_dynamicmtx[j][l] * cpmx1[j][r];
			}
			c->S[(size_t)r*c->nl2+u] = sc;
		}
	/* T: the same with the profiles swapped (the first row and column of each sweep) */
	c->T = msf_alloc( sizeof( double ) * (size_t)( lgth2 + 1 ) * ( c->nl1 ? c->nl1 : 1 ) );
	for( p=0; p<=lgth2; p++ )
		for( t=0; t<c->nl1; t++ )
		{
			int l = c->k1[t];
			double sc = 0.0;
			for( u=0; u<c->nl2; u++ )
			{
				int j = c->k2[u];
				sc += n_dynamicmtx[j][l] * cpmx2[j][p];
			}
			c->T[(size_t)p*c->nl1+t] = sc;
		}
	/* the slots whose S entry is nonzero somewhere */
	c->ns = 0;
	for( u=0; u<c->nl2; u++ )
	{
		for( r=0; r<=lgth1; r++ ) if( c->S[(size_t)r*c->nl2+u] != 0.0 ) break;
		if( r <= lgth1 ) c->su[c->ns++] = u;
	}
	c->zero = msf_calloc( lgth2 + 2, sizeof( double ) );
	if( c->ns <= 4 )
	{
		c->general = 0;
		for( t=0; t<4; t++ )
		{
			if( t < c->ns ) { c->s4[t] = c->su[t]; c->c4[t] = cpmx2[c->k2[c->su[t]]]; }
			else { c->s4[t] = -1; c->c4[t] = c->zero; }
		}
	}
	else
	{
		int ne = 0;
		c->general = 1;
		c->cst = msf_alloc( sizeof( int ) * ( lgth2 + 2 ) );
		for( p=0; p<lgth2; p++ ) for( t=0; t<c->ns; t++ ) if( cpmx2[c->k2[c->su[t]]][p] ) ne++;
		c->cu = msf_alloc( sizeof( int ) * ( ne + 1 ) );
		c->cw = msf_alloc( sizeof( double ) * ( ne + 1 ) );
		ne = 0;
		for( p=0; p<lgth2; p++ )
		{
			c->cst[p] = ne;
			for( t=0; t<c->ns; t++ )
			{
				double w = cpmx2[c->k2[c->su[t]]][p];
				if( w ) { c->cu[ne] = c->su[t]; c->cw[ne] = w; ne++; }
			}
		}
		c->cst[lgth2] = ne;
		c->cst[lgth2+1] = ne;
	}
	return( 1 );
}
static void msf_free( msf_t *c )
{
	free( c->S ); free( c->T ); free( c->zero );
	free( c->cst ); free( c->cu ); free( c->cw );
}

/* match_calc_alphabet_seq( n, out, cpmx1, cpmx2, r, jst, n ): out[j] = match of profile-1 column r with profile-2 column jst+j */
static void msf_matchrow( const msf_t *c, int r, int jst, int n, double * restrict out )
{
	const double *sr = c->S + (size_t)r * c->nl2;
	int j;
	if( !c->general )
	{
		const double s0 = c->s4[0] >= 0 ? sr[c->s4[0]] : 0.0, s1 = c->s4[1] >= 0 ? sr[c->s4[1]] : 0.0;
		const double s2 = c->s4[2] >= 0 ? sr[c->s4[2]] : 0.0, s3 = c->s4[3] >= 0 ? sr[c->s4[3]] : 0.0;
		const double *c0 = c->c4[0] + jst, *c1 = c->c4[1] + jst, *c2 = c->c4[2] + jst, *c3 = c->c4[3] + jst;
		if( c->c4[0] == c->zero ) c0 = c->zero; // (the zero row is not offset: it has fulllen2+2 entries)
		if( c->c4[1] == c->zero ) c1 = c->zero;
		if( c->c4[2] == c->zero ) c2 = c->zero;
		if( c->c4[3] == c->zero ) c3 = c->zero;
		for( j=0; j<n; j++ )
		{
			double base = 0.0;
			base += s0 * c0[j];
			base += s1 * c1[j];
			base += s2 * c2[j];
			base += s3 * c3[j];
			out[j] = base;
		}
	}
	else
	{
		const int *cst = c->cst + jst;
		for( j=0; j<n; j++ )
		{
			double base = 0.0;
			int e;
			for( e=cst[j]; e<cst[j+1]; e++ )
				base += sr[c->cu[e]] * c->cw[e];
			out[j] = base;
		}
	}
}
/* match_calc_alphabet_seq( n, out, cpmx2, cpmx1, p, ist, n ): out[i] = match of profile-2 column p with profile-1 column ist+i */
static void msf_matchcol( const msf_t *c, int p, int ist, int n, double * restrict out )
{
	const double *tp = c->T + (size_t)p * c->nl1;
	int i, t;
	for( i=0; i<n; i++ )
	{
		double base = 0.0;
		for( t=0; t<c->nl1; t++ )
			base += tp[t] * c->cpmx1[c->k1[t]][ist+i];
		out[i] = base;
	}
}

/* ---- the DP rows ---- */

/* One cell of a forward row (the original's statements, branch-free; see 2. above):
     wm = previousw[j-1]; row-gap candidate (>); mi update (>=); column-gap candidate (>);
     m[j] update (>=), mp[j] = i-1; currentw[j] = match + wm */
#define MSF_FCELL( J, PW, BASE ) \
	{ \
		double mj = m[J], pw = (PW), base = (BASE), gg, wh, g2, g3, wv, g4; \
		int takec; \
		gg = mil + fg2[(J)-1] * gf1va; \
		wh = ( gg > pw ) ? gg : pw; \
		g2 = pw + og2[J] * gf1vapre; \
		mil = MSF_DSEL_GE( g2, mil ); \
		g3 = mj + fgcp1va * gf2[J]; \
		wv = ( g3 > wh ) ? g3 : wh; \
		g4 = pw + ogcp1va * gf2[(J)-1]; \
		takec = ( g4 >= mj ); \
		m[J] = MSF_DSEL_GE( g4, mj ); \
		mp[J] = MSF_ISEL( takec, im1, mp[J] ); \
		cur[J] = base + wv; \
	}
/* the match score of column j: the chain of match_calc over the slots s4[] (ns <= 4, padded) */
#define MSF_MATCH4( J ) ( { double _b = 0.0; _b += s0 * c0[J]; _b += s1 * c1[J]; _b += s2 * c2[J]; _b += s3 * c3[J]; _b; } )

/* One cell of a backward row: wm = previousw[j+1]; row-gap candidate (>); mi update (>=);
   column-gap candidate (>); m[j] update (>=), mp[j] = i+1; currentw[j] = match + wm */
#define MSF_BCELL( J, PW, BASE ) \
	{ \
		double mj = m[J], pw = (PW), base = (BASE), gg, wh, g2, g3, wv, g4; \
		int takec; \
		gg = mil + og2[(J)+1] * gf1va; \
		wh = ( gg > pw ) ? gg : pw; \
		g2 = pw + fg2[J] * gf1vanext; \
		mil = MSF_DSEL_GE( g2, mil ); \
		g3 = mj + ogcp1vanext * gf2[J]; \
		wv = ( g3 > wh ) ? g3 : wh; \
		g4 = pw + fgcp1va * gf2[(J)+1]; \
		takec = ( g4 >= mj ); \
		m[J] = MSF_DSEL_GE( g4, mj ); \
		mp[J] = MSF_ISEL( takec, ip1, mp[J] ); \
		cur[J] = base + wv; \
	}



/* ---- x86: the plain-row cells, 4 (AVX2) or 2 (SSE2) per step ----
   Every value of a cell is computed by the same operations as in MSF_FCELL / MSF_BCELL, one lane per
   cell, written as the same C statements (vector operators on __m256d / __m128d), so a build that
   contracts a*b+c to fma (FMA in the target) does so in both, and one that does not (gcc's ISO mode,
   or no FMA in the target: the AVX2 functions are compiled with target("avx2") only) in neither.
   The selects are maxpd: MAXPD( a, b ) = a > b ? a : b, lane by lane, which is exactly
     wh = gg > pw ? gg : pw         -> max( gg, pw )
     wv = g3 > wh ? g3 : wh         -> max( g3, wh )
     m  = g4 >= mj ? g4 : mj        -> max( mj, g4 )     (equal: g4, as the scalar code)
     mi = g2 >= mi ? g2 : mi        -> max( mi, g2 )     (the scalar x86 code is the same maxsd)
   The one dependence between the cells of a row is the carry mi, a running "later element wins
   ties" maximum: max( earlier, later ) is associative for non-NaN values (both groupings return the
   latest of the maximal elements, equal values including +0.0/-0.0 counting as ties), so a prefix scan
   inside the vector (shifts by 1 and 2; a lane combined with itself is itself), combined with the
   incoming carry, gives each lane exactly the mi the scalar loop has before its cell, bit for bit
   (values are only selected, never recomputed).  The same argument as work/prog3's adpvec_body.h,
   without the column index (the plain rows do not need it).  Backward rows run with j descending:
   the lanes are the same four columns, in reverse sequence order. */
#if defined( __x86_64__ ) && ( defined( __GNUC__ ) || defined( __clang__ ) )
#define MSF_X86 1
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

typedef struct
{
	const double *c0, *c1, *c2, *c3;
	double s0, s1, s2, s3;
	int mrow;                   // match scores already in cur[] (general rows)
	double fgcp1va, ogcp1va, gf1va, gf1vx; // forward: gf1vx = gapfreq1f[i-1]; backward: ogcp1va = ogcp1[i+1], fgcp1va = fgcp1[i], gf1vx = gapfreq1f[i+1]
} msfv_t;

/* forward cells j = ja..jb-1 (each with a match score and pw = prev[j-1]); returns the carry */
#define MSF_FV4_BODY( MROW ) \
		for( ; j + 4 <= jb; j += 4 ) \
		{ \
			__m256d pw = _mm256_loadu_pd( prev + j - 1 ), mj = _mm256_loadu_pd( m + j ); \
			__m256d base, g2, s, t, ev, gg, wh, g3, wv, g4, mk; \
			__m128i mpo, mk4; \
			if( MROW ) base = _mm256_loadu_pd( cur + j ); \
			else \
			{ \
				base = _mm256_setzero_pd(); \
				base += vs0 * _mm256_loadu_pd( c0 + j ); \
				base += vs1 * _mm256_loadu_pd( c1 + j ); \
				base += vs2 * _mm256_loadu_pd( c2 + j ); \
				base += vs3 * _mm256_loadu_pd( c3 + j ); \
			} \
			g2 = pw + _mm256_loadu_pd( og2 + j ) * vgf1p; \
			s = _mm256_max_pd( _mm256_permute4x64_pd( g2, _MM_SHUFFLE( 2, 1, 0, 0 ) ), g2 ); \
			s = _mm256_max_pd( _mm256_permute4x64_pd( s, _MM_SHUFFLE( 1, 0, 1, 0 ) ), s ); \
			t = _mm256_blend_pd( _mm256_permute4x64_pd( s, _MM_SHUFFLE( 2, 1, 0, 0 ) ), cv, 1 ); \
			ev = _mm256_max_pd( cv, t ); \
			cv = _mm256_max_pd( cv, _mm256_permute4x64_pd( s, _MM_SHUFFLE( 3, 3, 3, 3 ) ) ); \
			gg = ev + _mm256_loadu_pd( fg2 + j - 1 ) * vgf1; \
			wh = _mm256_max_pd( gg, pw ); \
			g3 = mj + vfg1 * _mm256_loadu_pd( gf2 + j ); \
			wv = _mm256_max_pd( g3, wh ); \
			g4 = pw + vog1 * _mm256_loadu_pd( gf2 + j - 1 ); \
			_mm256_storeu_pd( m + j, _mm256_max_pd( mj, g4 ) ); \
			mk = _mm256_cmp_pd( g4, mj, _CMP_GE_OQ ); \
			mk4 = _mm256_castsi256_si128( _mm256_permutevar8x32_epi32( _mm256_castpd_si256( mk ), pk ) ); \
			mpo = _mm_loadu_si128( (const __m128i *)( mp + j ) ); \
			_mm_storeu_si128( (__m128i *)( mp + j ), _mm_blendv_epi8( mpo, vim1, mk4 ) ); \
			_mm256_storeu_pd( cur + j, base + wv ); \
		}

/* forward cells j = ja..jb-1 (each with a match score and pw = prev[j-1]); returns the carry */
__attribute__(( target( "avx2" ) ))
static double msf_fwd_v4( const msfv_t *v, int ja, int jb, int im1,
	double * restrict m, int * restrict mp, double * restrict cur, const double * restrict prev,
	const double * restrict og2, const double * restrict fg2, const double * restrict gf2, double mil )
{
	int j = ja;
	const double * restrict c0 = v->c0, * restrict c1 = v->c1, * restrict c2 = v->c2, * restrict c3 = v->c3;
	const int mrow = v->mrow;
	if( jb - ja >= 4 )
	{
		const __m256d vs0 = _mm256_set1_pd( v->s0 ), vs1 = _mm256_set1_pd( v->s1 ), vs2 = _mm256_set1_pd( v->s2 ), vs3 = _mm256_set1_pd( v->s3 );
		const __m256d vfg1 = _mm256_set1_pd( v->fgcp1va ), vog1 = _mm256_set1_pd( v->ogcp1va );
		const __m256d vgf1 = _mm256_set1_pd( v->gf1va ), vgf1p = _mm256_set1_pd( v->gf1vx );
		const __m128i vim1 = _mm_set1_epi32( im1 );
		const __m256i pk = _mm256_setr_epi32( 0, 2, 4, 6, 0, 2, 4, 6 );
		__m256d cv = _mm256_set1_pd( mil );
		if( mrow ) { MSF_FV4_BODY( 1 ) }
		else { MSF_FV4_BODY( 0 ) }
		mil = _mm256_cvtsd_f64( cv );
	}
	{
		const double s0 = v->s0, s1 = v->s1, s2 = v->s2, s3 = v->s3;
		const double fgcp1va = v->fgcp1va, ogcp1va = v->ogcp1va, gf1va = v->gf1va, gf1vapre = v->gf1vx;
		for( ; j<jb; j++ )
		{
			if( mrow ) MSF_FCELL( j, prev[j-1], cur[j] )
			else MSF_FCELL( j, prev[j-1], MSF_MATCH4( j ) )
		}
	}
	return( mil );
}

#define MSF_BV4_BODY( MROW ) \
		for( ; j - 3 >= ja; j -= 4 ) \
		{ \
			const int q = j - 3; \
			__m256d pw = _mm256_loadu_pd( prev + q + 1 ), mj = _mm256_loadu_pd( m + q ); \
			__m256d base, g2, s, t, ev, gg, wh, g3, wv, g4, mk; \
			__m128i mpo, mk4; \
			if( MROW ) base = _mm256_loadu_pd( cur + q ); \
			else \
			{ \
				base = _mm256_setzero_pd(); \
				base += vs0 * _mm256_loadu_pd( c0 + q ); \
				base += vs1 * _mm256_loadu_pd( c1 + q ); \
				base += vs2 * _mm256_loadu_pd( c2 + q ); \
				base += vs3 * _mm256_loadu_pd( c3 + q ); \
			} \
			g2 = pw + _mm256_loadu_pd( fg2 + q ) * vgf1n; \
			s = _mm256_max_pd( _mm256_permute4x64_pd( g2, _MM_SHUFFLE( 3, 3, 2, 1 ) ), g2 ); \
			s = _mm256_max_pd( _mm256_permute4x64_pd( s, _MM_SHUFFLE( 3, 2, 3, 2 ) ), s ); \
			t = _mm256_blend_pd( _mm256_permute4x64_pd( s, _MM_SHUFFLE( 3, 3, 2, 1 ) ), cv, 8 ); \
			ev = _mm256_max_pd( cv, t ); \
			cv = _mm256_max_pd( cv, _mm256_permute4x64_pd( s, _MM_SHUFFLE( 0, 0, 0, 0 ) ) ); \
			gg = ev + _mm256_loadu_pd( og2 + q + 1 ) * vgf1; \
			wh = _mm256_max_pd( gg, pw ); \
			g3 = mj + vog1n * _mm256_loadu_pd( gf2 + q ); \
			wv = _mm256_max_pd( g3, wh ); \
			g4 = pw + vfg1 * _mm256_loadu_pd( gf2 + q + 1 ); \
			_mm256_storeu_pd( m + q, _mm256_max_pd( mj, g4 ) ); \
			mk = _mm256_cmp_pd( g4, mj, _CMP_GE_OQ ); \
			mk4 = _mm256_castsi256_si128( _mm256_permutevar8x32_epi32( _mm256_castpd_si256( mk ), pk ) ); \
			mpo = _mm_loadu_si128( (const __m128i *)( mp + q ) ); \
			_mm_storeu_si128( (__m128i *)( mp + q ), _mm_blendv_epi8( mpo, vip1, mk4 ) ); \
			_mm256_storeu_pd( cur + q, base + wv ); \
		}

/* backward cells j = jb-1 down to ja (pw = prev[j+1]); lanes 0..3 = columns q..q+3, in sequence
   order 3, 2, 1, 0; returns the carry */
__attribute__(( target( "avx2" ) ))
static double msf_bwd_v4( const msfv_t *v, int ja, int jb, int ip1,
	double * restrict m, int * restrict mp, double * restrict cur, const double * restrict prev,
	const double * restrict og2, const double * restrict fg2, const double * restrict gf2, double mil )
{
	int j = jb - 1; // the highest column left
	const double * restrict c0 = v->c0, * restrict c1 = v->c1, * restrict c2 = v->c2, * restrict c3 = v->c3;
	const int mrow = v->mrow;
	if( jb - ja >= 4 )
	{
		const __m256d vs0 = _mm256_set1_pd( v->s0 ), vs1 = _mm256_set1_pd( v->s1 ), vs2 = _mm256_set1_pd( v->s2 ), vs3 = _mm256_set1_pd( v->s3 );
		const __m256d vfg1 = _mm256_set1_pd( v->fgcp1va ), vog1n = _mm256_set1_pd( v->ogcp1va );
		const __m256d vgf1 = _mm256_set1_pd( v->gf1va ), vgf1n = _mm256_set1_pd( v->gf1vx );
		const __m128i vip1 = _mm_set1_epi32( ip1 );
		const __m256i pk = _mm256_setr_epi32( 0, 2, 4, 6, 0, 2, 4, 6 );
		__m256d cv = _mm256_set1_pd( mil );
		if( mrow ) { MSF_BV4_BODY( 1 ) }
		else { MSF_BV4_BODY( 0 ) }
		mil = _mm256_cvtsd_f64( cv );
	}
	{
		const double s0 = v->s0, s1 = v->s1, s2 = v->s2, s3 = v->s3;
		const double fgcp1va = v->fgcp1va, ogcp1vanext = v->ogcp1va, gf1va = v->gf1va, gf1vanext = v->gf1vx;
		for( ; j>=ja; j-- )
		{
			if( mrow ) MSF_BCELL( j, prev[j+1], cur[j] )
			else MSF_BCELL( j, prev[j+1], MSF_MATCH4( j ) )
		}
	}
	return( mil );
}

/* the same, two lanes (SSE2, every x86-64) */
static double msf_fwd_v2( const msfv_t *v, int ja, int jb, int im1,
	double * restrict m, int * restrict mp, double * restrict cur, const double * restrict prev,
	const double * restrict og2, const double * restrict fg2, const double * restrict gf2, double mil )
{
	int j = ja;
	const double * restrict c0 = v->c0, * restrict c1 = v->c1, * restrict c2 = v->c2, * restrict c3 = v->c3;
	const int mrow = v->mrow;
	if( jb - ja >= 2 )
	{
		const __m128d vs0 = _mm_set1_pd( v->s0 ), vs1 = _mm_set1_pd( v->s1 ), vs2 = _mm_set1_pd( v->s2 ), vs3 = _mm_set1_pd( v->s3 );
		const __m128d vfg1 = _mm_set1_pd( v->fgcp1va ), vog1 = _mm_set1_pd( v->ogcp1va );
		const __m128d vgf1 = _mm_set1_pd( v->gf1va ), vgf1p = _mm_set1_pd( v->gf1vx );
		const __m128i vim1 = _mm_set1_epi32( im1 );
		__m128d cv = _mm_set1_pd( mil );
		for( ; j + 2 <= jb; j += 2 )
		{
			__m128d pw = _mm_loadu_pd( prev + j - 1 ), mj = _mm_loadu_pd( m + j );
			__m128d base, g2, s, t, ev, gg, wh, g3, wv, g4, mk;
			__m128i mpo, mk2, mpn;
			if( mrow ) base = _mm_loadu_pd( cur + j );
			else
			{
				base = _mm_setzero_pd();
				base += vs0 * _mm_loadu_pd( c0 + j );
				base += vs1 * _mm_loadu_pd( c1 + j );
				base += vs2 * _mm_loadu_pd( c2 + j );
				base += vs3 * _mm_loadu_pd( c3 + j );
			}
			g2 = pw + _mm_loadu_pd( og2 + j ) * vgf1p;
			s = _mm_max_pd( _mm_unpacklo_pd( g2, g2 ), g2 );          // (g0, max(g0,g1))
			t = _mm_unpacklo_pd( cv, s );                              // (c, s0)
			ev = _mm_max_pd( cv, t );
			cv = _mm_max_pd( cv, _mm_unpackhi_pd( s, s ) );
			gg = ev + _mm_loadu_pd( fg2 + j - 1 ) * vgf1;
			wh = _mm_max_pd( gg, pw );
			g3 = mj + vfg1 * _mm_loadu_pd( gf2 + j );
			wv = _mm_max_pd( g3, wh );
			g4 = pw + vog1 * _mm_loadu_pd( gf2 + j - 1 );
			_mm_storeu_pd( m + j, _mm_max_pd( mj, g4 ) );
			mk = _mm_cmpge_pd( g4, mj );
			mk2 = _mm_shuffle_epi32( _mm_castpd_si128( mk ), _MM_SHUFFLE( 2, 0, 2, 0 ) );
			mpo = _mm_loadl_epi64( (const __m128i *)( mp + j ) );
			mpn = _mm_or_si128( _mm_and_si128( mk2, vim1 ), _mm_andnot_si128( mk2, mpo ) );
			_mm_storel_epi64( (__m128i *)( mp + j ), mpn );
			_mm_storeu_pd( cur + j, base + wv );
		}
		mil = _mm_cvtsd_f64( cv );
	}
	{
		const double s0 = v->s0, s1 = v->s1, s2 = v->s2, s3 = v->s3;
		const double fgcp1va = v->fgcp1va, ogcp1va = v->ogcp1va, gf1va = v->gf1va, gf1vapre = v->gf1vx;
		for( ; j<jb; j++ )
		{
			if( mrow ) MSF_FCELL( j, prev[j-1], cur[j] )
			else MSF_FCELL( j, prev[j-1], MSF_MATCH4( j ) )
		}
	}
	return( mil );
}

static double msf_bwd_v2( const msfv_t *v, int ja, int jb, int ip1,
	double * restrict m, int * restrict mp, double * restrict cur, const double * restrict prev,
	const double * restrict og2, const double * restrict fg2, const double * restrict gf2, double mil )
{
	int j = jb - 1;
	const double * restrict c0 = v->c0, * restrict c1 = v->c1, * restrict c2 = v->c2, * restrict c3 = v->c3;
	const int mrow = v->mrow;
	if( jb - ja >= 2 )
	{
		const __m128d vs0 = _mm_set1_pd( v->s0 ), vs1 = _mm_set1_pd( v->s1 ), vs2 = _mm_set1_pd( v->s2 ), vs3 = _mm_set1_pd( v->s3 );
		const __m128d vfg1 = _mm_set1_pd( v->fgcp1va ), vog1n = _mm_set1_pd( v->ogcp1va );
		const __m128d vgf1 = _mm_set1_pd( v->gf1va ), vgf1n = _mm_set1_pd( v->gf1vx );
		const __m128i vip1 = _mm_set1_epi32( ip1 );
		__m128d cv = _mm_set1_pd( mil );
		for( ; j - 1 >= ja; j -= 2 )
		{
			const int q = j - 1; // lanes 0, 1 = columns q, j; sequence order 1, 0
			__m128d pw = _mm_loadu_pd( prev + q + 1 ), mj = _mm_loadu_pd( m + q );
			__m128d base, g2, s, t, ev, gg, wh, g3, wv, g4, mk;
			__m128i mpo, mk2, mpn;
			if( mrow ) base = _mm_loadu_pd( cur + q );
			else
			{
				base = _mm_setzero_pd();
				base += vs0 * _mm_loadu_pd( c0 + q );
				base += vs1 * _mm_loadu_pd( c1 + q );
				base += vs2 * _mm_loadu_pd( c2 + q );
				base += vs3 * _mm_loadu_pd( c3 + q );
			}
			g2 = pw + _mm_loadu_pd( fg2 + q ) * vgf1n;
			s = _mm_max_pd( _mm_unpackhi_pd( g2, g2 ), g2 );          // (max(g1,g0), g1)
			t = _mm_unpackhi_pd( s, cv );                              // (s1, c)
			ev = _mm_max_pd( cv, t );
			cv = _mm_max_pd( cv, _mm_unpacklo_pd( s, s ) );
			gg = ev + _mm_loadu_pd( og2 + q + 1 ) * vgf1;
			wh = _mm_max_pd( gg, pw );
			g3 = mj + vog1n * _mm_loadu_pd( gf2 + q );
			wv = _mm_max_pd( g3, wh );
			g4 = pw + vfg1 * _mm_loadu_pd( gf2 + q + 1 );
			_mm_storeu_pd( m + q, _mm_max_pd( mj, g4 ) );
			mk = _mm_cmpge_pd( g4, mj );
			mk2 = _mm_shuffle_epi32( _mm_castpd_si128( mk ), _MM_SHUFFLE( 2, 0, 2, 0 ) );
			mpo = _mm_loadl_epi64( (const __m128i *)( mp + q ) );
			mpn = _mm_or_si128( _mm_and_si128( mk2, vip1 ), _mm_andnot_si128( mk2, mpo ) );
			_mm_storel_epi64( (__m128i *)( mp + q ), mpn );
			_mm_storeu_pd( cur + q, base + wv );
		}
		mil = _mm_cvtsd_f64( cv );
	}
	{
		const double s0 = v->s0, s1 = v->s1, s2 = v->s2, s3 = v->s3;
		const double fgcp1va = v->fgcp1va, ogcp1vanext = v->ogcp1va, gf1va = v->gf1va, gf1vanext = v->gf1vx;
		for( ; j>=ja; j-- )
		{
			if( mrow ) MSF_BCELL( j, prev[j+1], cur[j] )
			else MSF_BCELL( j, prev[j+1], MSF_MATCH4( j ) )
		}
	}
	return( mil );
}
#endif

/* kernel width: MAFFT_ADPVEC=0/2/4 as for A__align (default: 4 on x86 with AVX2, 2 on other x86,
   scalar elsewhere); MAFFT_SIMD=scalar|sse2 lowers it too */
static int msf_vw = -1;
static int msf_vwidth( void )
{
	int w = __atomic_load_n( &msf_vw, __ATOMIC_ACQUIRE );
	if( w >= 0 ) return( w );
#ifdef MSF_X86
	{
		char *e = getenv( "MAFFT_ADPVEC" ), *s = getenv( "MAFFT_SIMD" );
		__builtin_cpu_init();
		w = __builtin_cpu_supports( "avx2" ) ? 4 : 2;
		if( s && !strcmp( s, "scalar" ) ) w = 0;
		else if( s && ( !strcmp( s, "sse2" ) || !strcmp( s, "popcnt" ) ) && w > 2 ) w = 2;
		if( e ) { int x = atoi( e ); w = ( x >= 4 && __builtin_cpu_supports( "avx2" ) ) ? 4 : ( x >= 2 ? 2 : 0 ); }
	}
#else
	w = 0;
#endif
	__atomic_store_n( &msf_vw, w, __ATOMIC_RELEASE );
	return( w );
}

/* forward row i, cells j = lo..hi-1 (1 <= lo, hi <= lgth2+1), previous row's value in column lo-1 =
   pwfirst, row carry mil before column lo; returns the carry after column hi-1.  Match scores
   inline (ns <= 4) or already in cur[j] (general).  Column lgth2 has no match score (cur[lgth2]
   keeps what the buffer holds, as in the original). */
static double msf_fwd_cells( const msf_t *c, int general, int i, int lo, int hi, int lgth2, int r, int jst,
	double * restrict m, int * restrict mp, double * restrict cur, const double * restrict prev, double pwfirst,
	const double * restrict og2, const double * restrict fg2, const double * restrict gf2,
	double fgcp1va, double ogcp1va, double gf1va, double gf1vapre, double mil )
{
	int j, hm = ( hi > lgth2 ) ? lgth2 : hi; // cells with a match score: lo..hm-1
	const int im1 = i - 1;
	const int vw = msf_vwidth();
	if( lo >= hi ) return( mil );
	if( !general )
	{
		const double *sr = c->S + (size_t)r * c->nl2;
		const double s0 = c->s4[0] >= 0 ? sr[c->s4[0]] : 0.0, s1 = c->s4[1] >= 0 ? sr[c->s4[1]] : 0.0;
		const double s2 = c->s4[2] >= 0 ? sr[c->s4[2]] : 0.0, s3 = c->s4[3] >= 0 ? sr[c->s4[3]] : 0.0;
		const double * restrict c0 = c->c4[0] + jst, * restrict c1 = c->c4[1] + jst, * restrict c2 = c->c4[2] + jst, * restrict c3 = c->c4[3] + jst;
		j = lo;
		if( j < hm ) { MSF_FCELL( j, pwfirst, MSF_MATCH4( j ) ) j++; }
		else { MSF_FCELL( j, pwfirst, cur[j] ) j++; } // (lo == lgth2: only the stale column)
#ifdef MSF_X86
		if( vw && hm - j >= vw )
		{
			msfv_t v;
			v.c0 = c0; v.c1 = c1; v.c2 = c2; v.c3 = c3; v.s0 = s0; v.s1 = s1; v.s2 = s2; v.s3 = s3; v.mrow = 0;
			v.fgcp1va = fgcp1va; v.ogcp1va = ogcp1va; v.gf1va = gf1va; v.gf1vx = gf1vapre;
			mil = ( vw == 4 ? msf_fwd_v4 : msf_fwd_v2 )( &v, j, hm, im1, m, mp, cur, prev, og2, fg2, gf2, mil );
			j = hm;
		}
#endif
		for( ; j<hm; j++ )
			MSF_FCELL( j, prev[j-1], MSF_MATCH4( j ) )
	}
	else
	{
		j = lo;
		MSF_FCELL( j, pwfirst, cur[j] )
		j = lo+1;
#ifdef MSF_X86
		if( vw && hm - j >= vw )
		{
			msfv_t v;
			memset( &v, 0, sizeof( v ) ); v.mrow = 1;
			v.fgcp1va = fgcp1va; v.ogcp1va = ogcp1va; v.gf1va = gf1va; v.gf1vx = gf1vapre;
			mil = ( vw == 4 ? msf_fwd_v4 : msf_fwd_v2 )( &v, j, hm, im1, m, mp, cur, prev, og2, fg2, gf2, mil );
			j = hm;
		}
#endif
		for( ; j<hm; j++ )
			MSF_FCELL( j, prev[j-1], cur[j] )
	}
	if( hi > lgth2 && lgth2 > lo ) // the stale column
	{
		j = lgth2;
		MSF_FCELL( j, prev[j-1], cur[j] )
	}
	return( mil );
}

/* backward row i (not imid, imid-1 or jumpi), cells j = hi-1 down to lo (0 <= lo, hi <= lgth2-1),
   next row's value in column hi = pwfirst, row carry before column hi-1 = mil; returns the carry
   after column lo */
static double msf_bwd_cells( const msf_t *c, int general, int i, int lo, int hi, int r, int jst,
	double * restrict m, int * restrict mp, double * restrict cur, const double * restrict prev, double pwfirst,
	const double * restrict og2, const double * restrict fg2, const double * restrict gf2,
	double fgcp1va, double ogcp1vanext, double gf1va, double gf1vanext, double mil )
{
	int j;
	const int ip1 = i + 1;
	const int vw = msf_vwidth();
	if( lo >= hi ) return( mil );
	if( !general )
	{
		const double *sr = c->S + (size_t)r * c->nl2;
		const double s0 = c->s4[0] >= 0 ? sr[c->s4[0]] : 0.0, s1 = c->s4[1] >= 0 ? sr[c->s4[1]] : 0.0;
		const double s2 = c->s4[2] >= 0 ? sr[c->s4[2]] : 0.0, s3 = c->s4[3] >= 0 ? sr[c->s4[3]] : 0.0;
		const double * restrict c0 = c->c4[0] + jst, * restrict c1 = c->c4[1] + jst, * restrict c2 = c->c4[2] + jst, * restrict c3 = c->c4[3] + jst;
		j = hi-1;
		MSF_BCELL( j, pwfirst, MSF_MATCH4( j ) )
		j = hi-2;
#ifdef MSF_X86
		if( vw && j - lo + 1 >= vw )
		{
			msfv_t v;
			v.c0 = c0; v.c1 = c1; v.c2 = c2; v.c3 = c3; v.s0 = s0; v.s1 = s1; v.s2 = s2; v.s3 = s3; v.mrow = 0;
			v.fgcp1va = fgcp1va; v.ogcp1va = ogcp1vanext; v.gf1va = gf1va; v.gf1vx = gf1vanext;
			mil = ( vw == 4 ? msf_bwd_v4 : msf_bwd_v2 )( &v, lo, j+1, ip1, m, mp, cur, prev, og2, fg2, gf2, mil );
			j = lo - 1;
		}
#endif
		for( ; j>=lo; j-- )
			MSF_BCELL( j, prev[j+1], MSF_MATCH4( j ) )
	}
	else
	{
		j = hi-1;
		MSF_BCELL( j, pwfirst, cur[j] )
		j = hi-2;
#ifdef MSF_X86
		if( vw && j - lo + 1 >= vw )
		{
			msfv_t v;
			memset( &v, 0, sizeof( v ) ); v.mrow = 1;
			v.fgcp1va = fgcp1va; v.ogcp1va = ogcp1vanext; v.gf1va = gf1va; v.gf1vx = gf1vanext;
			mil = ( vw == 4 ? msf_bwd_v4 : msf_bwd_v2 )( &v, lo, j+1, ip1, m, mp, cur, prev, og2, fg2, gf2, mil );
			j = lo - 1;
		}
#endif
		for( ; j>=lo; j-- )
			MSF_BCELL( j, prev[j+1], cur[j] )
	}
	return( mil );
}

/* the match scores of row r for columns lo..hi-1 (general rows) */
static void msf_matchrange( const msf_t *c, int r, int jst, int lo, int hi, double * restrict out )
{
	if( hi > lo ) msf_matchrow( c, r, jst + lo, hi - lo, out + lo );
}

/* A band of columns of a sweep's plain rows, processed by one thread at a time (see msf_phase_*) */
#define MSF_MAXB 32
typedef struct
{
	int lo, hi;                 // columns [lo, hi)
	volatile int done;          // the last row finished (forward: ascending; backward: descending)
	int off;                    // row i is entry i - off of bw[] and bmi[]
	double *bw;                 // per row: forward cur[hi-1], backward cur[lo] (the neighbour's pw)
	double *bmi;                // per row: the row carry after the band's last cell
	char pad[64];
} msband_t;
typedef struct msf_phase_s
{
	struct msf_sweep_s *w;
	int dir;                    // 0 forward (rows r0..r1 ascending), 1 backward (rows r0 down to r1)
	int r0, r1, nb;
	msband_t band[MSF_MAXB];
	int claim;                  // next band to take (atomic)
	int finished;               // bands done (atomic)
	int refs;                   // helpers inside msf_phase_work (atomic)
} msf_phase_t;

/* the state a node's two sweeps share */
typedef struct msf_sweep_s
{
	const msf_t *c;
	int ist, ien, jst, jen, lgth1, lgth2, imid;
	double *midw, *midm, *midn;
	int *jumpbacki, *jumpbackj;
	volatile int fwddone;      // the forward sweep has written the mid arrays (release / acquire)
	/* the forward sweep */
	double *fw[2], *finitv, *fm;
	int *fmp;
	/* the backward sweep, between msf_backward_start and msf_backward_finish */
	double *bw1, *bw2, *binitv, *bm, *bcur, *bprev, firstm;
	int *bmp, *jumpforwi, *jumpforwj, firstmp, bjumpi, bjumpj, bjmid, bimid;
	/* results of the backward sweep */
	int jmid, jumpi, jumpj, imidout;
} msf_sweep_t;

static void msf_spinwait( int *spins )
{
	if( ++*spins <= 64 ) ;
	else if( *spins <= 4096 ) sched_yield();
	else { struct timespec ts; ts.tv_sec = 0; ts.tv_nsec = 20000; nanosleep( &ts, NULL ); }
}

/* forward rows r0..r1 of band b: row i in buffer i&1 (as the serial loop's alternation); the cell in
   column lo reads the previous row's value in column lo-1 from the left neighbour's bw[], and the
   carry from its bmi[] */
static void msf_fwd_band( msf_phase_t *P, int b )
{
	msf_sweep_t *w = P->w;
	const msf_t *c = w->c;
	msband_t *mb = P->band + b, *nbr = b ? P->band + b - 1 : NULL;
	const int ist = w->ist, jst = w->jst, lgth2 = w->lgth2, lo = mb->lo, hi = mb->hi;
	const double *ogcp1 = c->gapinfo[0] + ist, *fgcp1 = c->gapinfo[1] + ist;
	const double *ogcp2 = c->gapinfo[2] + jst, *fgcp2 = c->gapinfo[3] + jst;
	const double *gapfreq1f = c->gapinfo[4] + ist, *gapfreq2f = c->gapinfo[5] + jst;
	const double *initverticalw = w->finitv;
	int i, spins = 0;
	for( i=P->r0; i<=P->r1; i++ )
	{
		double *currentw = w->fw[i&1], *previousw = w->fw[(i-1)&1], pwf, mi;
		if( nbr ) while( __atomic_load_n( &nbr->done, __ATOMIC_ACQUIRE ) < i ) msf_spinwait( &spins );
		spins = 0;
		if( b == 0 )
		{
			double pw0 = initverticalw[i-1]; // previousw[0] = initverticalw[i-1]
			pwf = pw0;
			mi = pw0 + ogcp2[1] * gapfreq1f[i-1];
			currentw[0] = initverticalw[i];
		}
		else
		{
			pwf = nbr->bw[i-1-nbr->off];
			mi = nbr->bmi[i-nbr->off];
		}
		if( c->general ) msf_matchrange( c, ist+i, jst, lo, hi > lgth2 ? lgth2 : hi, currentw );
		mi = msf_fwd_cells( c, c->general, i, lo, hi, lgth2, ist+i, jst, w->fm, w->fmp, currentw, previousw, pwf,
			ogcp2, fgcp2, gapfreq2f, fgcp1[i-1], ogcp1[i], gapfreq1f[i], gapfreq1f[i-1], mi );
		mb->bw[i-mb->off] = currentw[hi-1];
		mb->bmi[i-mb->off] = mi;
		__atomic_store_n( &mb->done, i, __ATOMIC_RELEASE );
	}
}

/* backward rows r0 down to r1 of band b (columns descending; band 0 is the rightmost): row i in
   buffer bw1/bw2 by ( lgth1-1-i ) & 1 (the serial loop's alternation, row lgth1-1 in bw1) */
static void msf_bwd_band( msf_phase_t *P, int b )
{
	msf_sweep_t *w = P->w;
	const msf_t *c = w->c;
	msband_t *mb = P->band + b, *nbr = b ? P->band + b - 1 : NULL;
	const int ist = w->ist, jst = w->jst, lgth1 = w->lgth1, lgth2 = w->lgth2, lo = mb->lo, hi = mb->hi;
	const double *ogcp1 = c->gapinfo[0] + ist, *fgcp1 = c->gapinfo[1] + ist;
	const double *ogcp2 = c->gapinfo[2] + jst, *fgcp2 = c->gapinfo[3] + jst;
	const double *gapfreq1f = c->gapinfo[4] + ist, *gapfreq2f = c->gapinfo[5] + jst;
	const double *initverticalw = w->binitv;
	double *buf[2];
	int i, spins = 0;
	buf[0] = w->bw1; buf[1] = w->bw2;
	for( i=P->r0; i>=P->r1; i-- )
	{
		double *currentw = buf[( lgth1-1-i ) & 1], *previousw = buf[( lgth1-i ) & 1], pwf, mi;
		if( nbr ) while( __atomic_load_n( &nbr->done, __ATOMIC_ACQUIRE ) > i ) msf_spinwait( &spins );
		spins = 0;
		if( b == 0 )
		{
			double pwl = initverticalw[i+1]; // previousw[lgth2-1] = initverticalw[i+1]
			pwf = pwl;
			mi = pwl + fgcp2[lgth2-2] * gapfreq1f[i+1];
			currentw[lgth2-1] = initverticalw[i];
		}
		else
		{
			pwf = nbr->bw[i+1-nbr->off];
			mi = nbr->bmi[i-nbr->off];
		}
		if( c->general ) msf_matchrange( c, ist+i, jst, lo, hi, currentw );
		mi = msf_bwd_cells( c, c->general, i, lo, hi, ist+i, jst, w->bm, w->bmp, currentw, previousw, pwf,
			ogcp2, fgcp2, gapfreq2f, fgcp1[i], ogcp1[i+1], gapfreq1f[i], gapfreq1f[i+1], mi );
		if( lo == 0 ) // the band with column 0: the running maximum of the original's row end
		{
			double g = previousw[0] + fgcp1[i];
			if( w->firstm < g )
			{
				w->firstm = g;
				w->firstmp = i + 1;
			}
		}
		mb->bw[i-mb->off] = currentw[lo];
		mb->bmi[i-mb->off] = mi;
		__atomic_store_n( &mb->done, i, __ATOMIC_RELEASE );
	}
}

/* take bands of a phase until none is left; bands are taken in increasing order, so the band a
   taker waits for is always being processed by somebody */
static void msf_phase_work( msf_phase_t *P )
{
	int b;
	while( ( b = __atomic_fetch_add( &P->claim, 1, __ATOMIC_ACQ_REL ) ) < P->nb )
	{
		if( P->dir == 0 ) msf_fwd_band( P, b );
		else msf_bwd_band( P, b );
		__atomic_add_fetch( &P->finished, 1, __ATOMIC_ACQ_REL );
	}
}

/* set up the bands of a phase: nb bands of the columns lo0..hi0-1 (forward: 1..lgth2+1, backward:
   0..lgth2-1), boundaries on multiples of 16 columns (no shared cache lines in the row arrays) */
static void msf_phase_init( msf_phase_t *P, msf_sweep_t *w, int dir, int r0, int r1, int nb, int lo0, int hi0, int nrowsalloc )
{
	int b, ncol = hi0 - lo0;
	memset( P, 0, sizeof( msf_phase_t ) );
	P->w = w; P->dir = dir; P->r0 = r0; P->r1 = r1;
	if( nb > MSF_MAXB ) nb = MSF_MAXB;
	if( nb < 1 ) nb = 1;
	P->nb = nb;
	for( b=0; b<nb; b++ )
	{
		/* forward: band 0 leftmost; backward: band 0 rightmost */
		int a = lo0 + (int)( (double)ncol * b / nb ), e = lo0 + (int)( (double)ncol * ( b + 1 ) / nb );
		if( b > 0 ) a = ( a + 15 ) & ~15;
		if( b < nb-1 ) e = ( e + 15 ) & ~15;
		if( a < lo0 ) a = lo0;
		if( e > hi0 ) e = hi0;
		if( dir == 0 ) { P->band[b].lo = a; P->band[b].hi = e; }
		else { P->band[nb-1-b].lo = a; P->band[nb-1-b].hi = e; }
	}
	for( b=0; b<nb; b++ ) // rows r0-1..r1 (forward) or r1..r0+1 (backward)
	{
		int nr = ( dir == 0 ? r1 - r0 : r0 - r1 ) + 3;
		if( nr < 3 ) nr = 3;
		P->band[b].done = dir == 0 ? r0 - 1 : r0 + 1;
		P->band[b].off = dir == 0 ? r0 - 1 : r1;
		P->band[b].bw = msf_alloc( sizeof( double ) * nr );
		P->band[b].bmi = msf_alloc( sizeof( double ) * nr );
	}
}
static void msf_phase_free( msf_phase_t *P )
{
	int b;
	for( b=0; b<P->nb; b++ ) { free( P->band[b].bw ); free( P->band[b].bmi ); }
}

/* ---- the forward sweep: set-up, plain rows 1..imid-1 (serial or banded), row imid ---- */
static void msf_forward_setup( msf_sweep_t *w )
{
	const msf_t *c = w->c;
	const int ist = w->ist, jst = w->jst, lgth1 = w->lgth1, lgth2 = w->lgth2;
	const double *ogcp1 = c->gapinfo[0] + ist, *fgcp1 = c->gapinfo[1] + ist;
	const double *ogcp2 = c->gapinfo[2] + jst, *fgcp2 = c->gapinfo[3] + jst;
	const double *gapfreq1f = c->gapinfo[4] + ist, *gapfreq2f = c->gapinfo[5] + jst;
	double headgapfreq1, headgapfreq2, *initverticalw, *currentw, *m;
	int *mp, i, j;

	if( ist > 0 ) headgapfreq1 = gapfreq1f[-1];
	else headgapfreq1 = c->hgf1g;
	if( jst > 0 ) headgapfreq2 = gapfreq2f[-1];
	else headgapfreq2 = c->hgf2g;

	w->fw[0] = msf_calloc( lgth2+2, sizeof( double ) ); // w1 (row 0, and every even row)
	w->fw[1] = msf_calloc( lgth2+2, sizeof( double ) ); // w2
	w->finitv = initverticalw = msf_calloc( lgth1+2, sizeof( double ) );
	w->fm = m = msf_calloc( lgth2+2, sizeof( double ) );
	w->fmp = mp = msf_calloc( lgth2+2, sizeof( int ) );
	currentw = w->fw[0];

	msf_matchcol( c, jst, ist, lgth1, initverticalw );
	msf_matchrow( c, ist, jst, lgth2, currentw );
	for( i=1; i<lgth1+1; i++ )
	{
		initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2f[0] ) ;
	}
	for( j=1; j<lgth2+1; j++ )
	{
		currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1f[0]) ;
	}
	for( j=1; j<lgth2+1; ++j )
	{
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2f[j-1];
		mp[j] = 0;
	}
}

/* rows 1..imid-1, one thread, all columns */
static void msf_forward_rows( msf_sweep_t *w )
{
	const msf_t *c = w->c;
	const int ist = w->ist, jst = w->jst, lgth2 = w->lgth2, imid = w->imid;
	const double *ogcp1 = c->gapinfo[0] + ist, *fgcp1 = c->gapinfo[1] + ist;
	const double *ogcp2 = c->gapinfo[2] + jst, *fgcp2 = c->gapinfo[3] + jst;
	const double *gapfreq1f = c->gapinfo[4] + ist, *gapfreq2f = c->gapinfo[5] + jst;
	const double *initverticalw = w->finitv;
	int i;
	for( i=1; i<imid; i++ )
	{
		double *currentw = w->fw[i&1], *previousw = w->fw[(i-1)&1], mi;
		previousw[0] = initverticalw[i-1];
		if( c->general ) msf_matchrow( c, ist+i, jst, lgth2, currentw );
		currentw[0] = initverticalw[i];
		mi = previousw[0] + ogcp2[1] * gapfreq1f[i-1];
		msf_fwd_cells( c, c->general, i, 1, lgth2+1, lgth2, ist+i, jst, w->fm, w->fmp, currentw, previousw, previousw[0],
			ogcp2, fgcp2, gapfreq2f, fgcp1[i-1], ogcp1[i], gapfreq1f[i], gapfreq1f[i-1], mi );
	}
}

/* row imid (the original statements) and the end of the sweep */
static void msf_forward_mid( msf_sweep_t *w )
{
	const msf_t *c = w->c;
	const int ist = w->ist, jst = w->jst, lgth2 = w->lgth2, imid = w->imid;
	const double *ogcp1 = c->gapinfo[0] + ist, *fgcp1 = c->gapinfo[1] + ist;
	const double *ogcp2 = c->gapinfo[2] + jst, *fgcp2 = c->gapinfo[3] + jst;
	const double *gapfreq1f = c->gapinfo[4] + ist, *gapfreq2f = c->gapinfo[5] + jst;
	const double *initverticalw = w->finitv;
	double *m = w->fm, *midw = w->midw, *midm = w->midm, *midn = w->midn, mi, g, wm;
	int *mp = w->fmp, *jumpbacki = w->jumpbacki, *jumpbackj = w->jumpbackj, i = imid, j, mpi;
	double *currentw = w->fw[i&1], *previousw = w->fw[(i-1)&1];

	previousw[0] = initverticalw[i-1];
	msf_matchrow( c, ist+i, jst, lgth2, currentw );
	currentw[0] = initverticalw[i];
	m[0] = ogcp1[i];
	midm[0] = m[0];
	mi = previousw[0] + ogcp2[1] * gapfreq1f[i-1];
	mpi = 0;
	for( j=1; j<lgth2+1; j++ )
	{
		wm = previousw[j-1];
		g = mi + fgcp2[j-1] * gapfreq1f[i];
		if( g > wm )
		{
			wm = g;
		}
		g = previousw[j-1] + ogcp2[j] * gapfreq1f[i-1];
		if( g >= mi )
		{
			mi = g;
			mpi = j-1;
		}
		g = m[j] + fgcp1[i-1] * gapfreq2f[j];
		if( g > wm )
		{
			wm = g;
		}
		g = previousw[j-1] + ogcp1[i] * gapfreq2f[j-1];
		if( g >= m[j] )
		{
			m[j] = g;
			mp[j] = i-1;
		}
		currentw[j] += wm;
		jumpbackj[j] = mp[j];
		jumpbacki[j] = mpi;
		midw[j] = currentw[j];
		midm[j] = m[j];
		midn[j] = mi;
	}
	free( w->fw[0] ); free( w->fw[1] ); free( w->finitv ); free( w->fm ); free( w->fmp );
	__atomic_store_n( &w->fwddone, 1, __ATOMIC_RELEASE );
}

static void msf_forward( msf_sweep_t *w )
{
	msf_forward_setup( w );
	msf_forward_rows( w );
	msf_forward_mid( w );
}

static void msf_waitfwd( msf_sweep_t *w )
{
	int spins = 0;
	while( !__atomic_load_n( &w->fwddone, __ATOMIC_ACQUIRE ) ) msf_spinwait( &spins );
}

/* ---- the backward sweep ---- */

/* one row i of the backward sweep (the loop body of MSalignmm_rec's second loop); 1: the loop ends here */
static int msf_bwd_row( msf_sweep_t *w, int i )
{
	const msf_t *c = w->c;
	const int ist = w->ist, jst = w->jst, lgth2 = w->lgth2;
	int imid = w->bimid, jumpi = w->bjumpi, jumpj = w->bjumpj, jmid = w->bjmid;
	const double *ogcp1 = c->gapinfo[0] + ist, *fgcp1 = c->gapinfo[1] + ist;
	const double *ogcp2 = c->gapinfo[2] + jst, *fgcp2 = c->gapinfo[3] + jst;
	const double *gapfreq1f = c->gapinfo[4] + ist, *gapfreq2f = c->gapinfo[5] + jst;
	double *currentw, *previousw, *m = w->bm, *initverticalw = w->binitv, mi, g, wm, maxwm;
	int *mp = w->bmp, *jumpforwi = w->jumpforwi, *jumpforwj = w->jumpforwj, j, mpi, ijpi, ijpj, ret = 0;
	double *midw = w->midw, *midm = w->midm, *midn = w->midn;
	int *jumpbacki = w->jumpbacki, *jumpbackj = w->jumpbackj;
	int special = ( i == imid || i == imid - 1 || i == jumpi );

	previousw = w->bcur;
	currentw = w->bprev;
	previousw[lgth2-1] = initverticalw[i+1];
	if( c->general || special ) msf_matchrow( c, ist+i, jst, lgth2, currentw );
	currentw[lgth2-1] = initverticalw[i];
	mi = previousw[lgth2-1] + fgcp2[lgth2-2] * gapfreq1f[i+1];
	if( !special )
	{
		msf_bwd_cells( c, c->general, i, 0, lgth2-1, ist+i, jst, m, mp, currentw, previousw, previousw[lgth2-1], ogcp2, fgcp2, gapfreq2f,
			fgcp1[i], ogcp1[i+1], gapfreq1f[i], gapfreq1f[i+1], mi );
	}
	else /* the original statements */
	{
		mpi = lgth2 - 1;
		for( j=lgth2-2; j>-1; j-- )
		{
			wm = previousw[j+1];
			ijpi = i+1;
			ijpj = j+1;
			g = mi + ogcp2[j+1] * gapfreq1f[i];
			if( g > wm )
			{
				wm = g;
				ijpj = mpi;
				ijpi = i+1;
			}
			g = previousw[j+1] + fgcp2[j] * gapfreq1f[i+1];
			if( g >= mi )
			{
				mi = g;
				mpi = j + 1;
			}
			g = m[j] + ogcp1[i+1] * gapfreq2f[j];
			if( g > wm )
			{
				wm = g;
				ijpi = mp[j];
				ijpj = j+1;
			}
			g = previousw[j+1] + fgcp1[i] * gapfreq2f[j+1];
			if( g >= m[j] )
			{
				m[j] = g;
				mp[j] = i + 1;
			}
			if( i == jumpi || i == imid - 1 )
			{
				jumpforwi[j] = ijpi;
				jumpforwj[j] = ijpj;
			}
			if( i == imid )
			{
				midw[j] += wm;
				midm[j+1] += m[j];
			}
			if( i == imid - 1 )
			{
				midn[j] += mi;
			}
			currentw[j] += wm;
		}
	}
	w->bcur = currentw;
	w->bprev = previousw;
	g = previousw[0] + fgcp1[i];
	if( w->firstm < g )
	{
		w->firstm = g;
		w->firstmp = i + 1;
	}
	if( i == imid ) midm[0] += w->firstm;

	if( i == imid - 1 )
	{
		maxwm = midw[1];
		jmid = 0;
		for( j=2; j<lgth2-1; j++ )
		{
			wm = midw[j];
			if( wm > maxwm )
			{
				jmid = j;
				maxwm = wm;
			}
		}
		for( j=0; j<lgth2+1; j++ )
		{
			wm = midm[j];
			if( wm > maxwm )
			{
				jmid = j;
				maxwm = wm;
			}
		}
		wm = midw[jmid];
		jumpi = imid-1;
		jumpj = jmid-1;
		if( jmid > 0 && midn[jmid-1] > wm ) //060413
		{
			jumpi = imid-1;
			jumpj = jumpbacki[jmid];
			wm = midn[jmid-1];
		}
		if( midm[jmid] > wm )
		{
			jumpi = jumpbackj[jmid];
			jumpj = jmid-1;
			wm = midm[jmid];
		}
	}
	if( i == jumpi )
	{
		if( jmid == 0 )
		{
			jumpj = 0; jmid = 1;
			if( imid < w->firstmp-1 )
			{
				jumpi = w->firstmp;
				imid = w->firstmp+1;
			}
		}
		else if( jmid >= lgth2 )
		{
			jumpi=imid-1; jmid=lgth2;
			jumpj = lgth2-1;
		}
		else
		{
			imid = jumpforwi[jumpj];
			jmid = jumpforwj[jumpj];
			if( imid == jumpi ) jumpi = imid-1;
		}
		ret = 1;
	}
	w->bimid = imid; w->bjumpi = jumpi; w->bjumpj = jumpj; w->bjmid = jmid;
	return( ret );
}

/* the backward sweep's set-up */
static void msf_backward_setup( msf_sweep_t *w )
{
	const msf_t *c = w->c;
	const int ist = w->ist, jst = w->jst, lgth1 = w->lgth1, lgth2 = w->lgth2;
	const double *ogcp1 = c->gapinfo[0] + ist, *fgcp1 = c->gapinfo[1] + ist;
	const double *ogcp2 = c->gapinfo[2] + jst, *fgcp2 = c->gapinfo[3] + jst;
	const double *gapfreq1f = c->gapinfo[4] + ist, *gapfreq2f = c->gapinfo[5] + jst;
	double *initverticalw, *m, *currentw;
	int *mp, i, j;

	w->bw1 = msf_calloc( lgth2+2, sizeof( double ) );
	w->bw2 = msf_calloc( lgth2+2, sizeof( double ) );
	w->binitv = initverticalw = msf_calloc( lgth1+2, sizeof( double ) );
	w->bm = m = msf_calloc( lgth2+2, sizeof( double ) );
	w->bmp = mp = msf_calloc( lgth2+2, sizeof( int ) );
	w->jumpforwi = msf_calloc( lgth2+2, sizeof( int ) );
	w->jumpforwj = msf_calloc( lgth2+2, sizeof( int ) );
	currentw = w->bw1;

	msf_matchcol( c, jst+lgth2-1, ist, lgth1, initverticalw );
	msf_matchrow( c, ist+lgth1-1, jst, lgth2, currentw );
	for( i=0; i<lgth1-1; i++ )
	{
		initverticalw[i] += ( fgcp1[lgth1-1] * gapfreq2f[lgth2] + ogcp1[i+1] * gapfreq2f[lgth2-1] );
	}
	for( j=0; j<lgth2-1; j++ )
	{
		currentw[j] += ( fgcp2[lgth2-1] * gapfreq1f[lgth1] + ogcp2[j+1] * gapfreq1f[lgth1-1] );
	}
	for( j=lgth2-1; j>-1; --j )
	{
		m[j] = currentw[j+1] + fgcp1[lgth1-2] * gapfreq2f[j+1]; // (m[lgth2-1] is never read; the original's value here is the forward sweep's stale column)
		mp[j] = lgth1-1;
	}
	w->bcur = w->bw1;
	w->bprev = w->bw2;
	w->firstm = -9999999.9;
	w->firstmp = lgth1;
	w->bjumpi = 0;
	w->bjumpj = 0;
	w->bjmid = 0;
	w->bimid = w->imid;
}

/* rows lgth1-2 .. imid+1, one thread (independent of the forward sweep) */
static void msf_backward_rows( msf_sweep_t *w )
{
	int i;
	for( i=w->lgth1-2; i>w->imid; i-- ) msf_bwd_row( w, i ); // (i == jumpi == 0 cannot happen above imid >= 1)
}

/* after the banded rows: the buffers as the serial loop leaves them after row imid+1 */
static void msf_backward_afterbands( msf_sweep_t *w )
{
	int i = w->imid + 1; // the last row done
	double *buf[2];
	buf[0] = w->bw1; buf[1] = w->bw2;
	w->bcur = buf[( w->lgth1-1-i ) & 1];
	w->bprev = buf[( w->lgth1-i ) & 1];
}

/* the rest: rows imid .. jumpi, after the forward sweep */
static void msf_backward_finish( msf_sweep_t *w )
{
	int i;
	msf_waitfwd( w );
	for( i=w->imid; i>-1; i-- ) if( msf_bwd_row( w, i ) ) break;
	w->jmid = w->bjmid; w->jumpi = w->bjumpi; w->jumpj = w->bjumpj; w->imidout = w->bimid;
	free( w->bw1 ); free( w->bw2 ); free( w->binitv ); free( w->bm ); free( w->bmp ); free( w->jumpforwi ); free( w->jumpforwj );
}

/* ---- MSalignmm_tanni: full DP with a traceback, gap tables out (the rows are made at the end) ---- */
static TLS int *msf_ijpblk = NULL;
static TLS size_t msf_ijpcap = 0;
#define MSF_IJPKEEP ( (size_t)1 << 20 ) // ints (4 MB): a per-thread block up to this size is kept between calls

/* the tables go to t1/t2 (lgth1+lgth2+1 chars each, terminated); returns what MSalignmm_tanni returns */
static double msf_tanni( const msf_t *c, int ist, int ien, int jst, int jen, char *t1, char *t2 )
{
	int i, j, l, iin, jin, ifi, jfi, k, klim;
	int lgth1, lgth2, lasti, mpi, stride;
	double wm = 0.0, g, mi;
	double *w1, *w2, *initverticalw, *lastverticalw, *m, *currentw, *previousw, *wtmp;
	int *mp, *ijpb;
	char *gt1bk, *gt2bk, *gaptable1, *gaptable2;
	const double *ogcp1 = c->gapinfo[0] + ist, *fgcp1 = c->gapinfo[1] + ist;
	const double *ogcp2 = c->gapinfo[2] + jst, *fgcp2 = c->gapinfo[3] + jst;
	const double *gapfreq1f = c->gapinfo[4] + ist, *gapfreq2f = c->gapinfo[5] + jst;
	double headgapfreq1, headgapfreq2;
	const int headgp = c->headgp, tailgp = c->tailgp, fulllen1 = c->fulllen1, fulllen2 = c->fulllen2;
	size_t need;

	if( ist > 0 ) headgapfreq1 = gapfreq1f[-1];
	else headgapfreq1 = c->hgf1g;
	if( jst > 0 ) headgapfreq2 = gapfreq2f[-1];
	else headgapfreq2 = c->hgf2g;

	lgth1 = ien-ist+1;
	lgth2 = jen-jst+1;

	w1 = msf_calloc( lgth2+102, sizeof( double ) );
	w2 = msf_calloc( lgth2+102, sizeof( double ) );
	initverticalw = msf_calloc( lgth1+102, sizeof( double ) );
	lastverticalw = msf_calloc( lgth1+102, sizeof( double ) );
	m = msf_calloc( lgth2+102, sizeof( double ) );
	mp = msf_calloc( lgth2+102, sizeof( int ) );

	/* the traceback cells (0..lgth1) x (0..lgth2), in a reused per-thread block: every cell the walk
	   below reads is written first (rows 1..lasti-1 by the DP, row and column 0 and, when the last
	   row is not computed, cell (lgth1, lgth2) by the walk's setup), as in MSalignmm_tanni */
	stride = lgth2 + 1;
	need = (size_t)( lgth1 + 1 ) * stride;
	if( need > MSF_IJPKEEP ) ijpb = msf_alloc( sizeof( int ) * need ); // large: only for this call (memory: one per thread otherwise)
	else
	{
		if( need > msf_ijpcap )
		{
			free( msf_ijpblk );
			msf_ijpcap = need + need / 4;
			if( msf_ijpcap > MSF_IJPKEEP ) msf_ijpcap = MSF_IJPKEEP;
			msf_ijpblk = msf_alloc( sizeof( int ) * msf_ijpcap );
		}
		ijpb = msf_ijpblk;
	}
#define MSF_IJP( I, J ) ijpb[(size_t)(I) * stride + (J)]

	currentw = w1;
	previousw = w2;

	msf_matchcol( c, jst, ist, lgth1, initverticalw );
	msf_matchrow( c, ist, jst, lgth2, currentw );

	if( headgp || ist != 0 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2f[0] ) ;
		}
	}
	if( headgp || jst != 0 )
	{
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1f[0] ) ;
		}
	}

	for( j=1; j<lgth2+1; ++j )
	{
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2f[j-1]; mp[j] = 0;;
	}

	lastverticalw[0] = currentw[lgth2-1];

	if( tailgp || jen != fulllen2-1 ) lasti = lgth1+1; else lasti = lgth1;
	for( i=1; i<lasti; i++ )
	{
		const double fgcp1va = fgcp1[i-1], ogcp1va = ogcp1[i], gf1va = gapfreq1f[i], gf1vapre = gapfreq1f[i-1];
		int *ijpi = &MSF_IJP( i, 0 );
		const int im1 = i - 1;
		int mpil;
		double mil;
		wtmp = previousw;
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

		msf_matchrow( c, ist+i, jst, lgth2, currentw );
		currentw[0] = initverticalw[i];

		mi = previousw[0] + ogcp2[1] * gapfreq1f[i-1];
		mpi = 0;
		mil = mi; mpil = mpi;
		for( j=1; j<lgth2+1; j++ )
		{
			/* the original cell, branch-free: wm = previousw[j-1]; ijp = 0; the row-gap candidate
			   (> , ijp = -(j-mpi) with the old mpi); mi update (>=); the column-gap candidate
			   (>, ijp = +(i-mp[j]) with the old mp[j]); m[j] update (>=); currentw[j] += wm */
			double mj = m[j], pw = previousw[j-1], gg, wh, g2, g3, g4;
			int th, take, takec, ijv, mpj = mp[j];
			gg = mil + fgcp2[j-1] * gf1va;
			th = ( gg > pw );
			wh = th ? gg : pw;
			ijv = -( j - mpil ) & -th;
			g2 = pw + ogcp2[j] * gf1vapre;
			take = ( g2 >= mil );
			mil = MSF_DSEL_GE( g2, mil );
			mpil = MSF_ISEL( take, j-1, mpil );
			g3 = mj + fgcp1va * gapfreq2f[j];
			th = ( g3 > wh );
			wm = th ? g3 : wh;
			ijv = MSF_ISEL( th, i - mpj, ijv );
			g4 = pw + ogcp1va * gapfreq2f[j-1];
			takec = ( g4 >= mj );
			m[j] = MSF_DSEL_GE( g4, mj );
			mp[j] = MSF_ISEL( takec, im1, mpj );
			currentw[j] += wm;
			ijpi[j] = ijv;
		}
		lastverticalw[i] = currentw[lgth2-1];
	}

	/* Atracking( currentw, lastverticalw, ..., tailgp, &gt1, &gt2 ): the gap tables */
	gt1bk = msf_alloc( lgth1+lgth2+1 );
	gt2bk = msf_alloc( lgth1+lgth2+1 );
	gaptable1 = gt1bk + lgth1+lgth2;
	gaptable2 = gt2bk + lgth1+lgth2;
	*gaptable1 = 0;
	*gaptable2 = 0;
	if( tailgp == 1 )
		;
	else if( ien == fulllen1-1 || jen == fulllen2-1 )
	{
		double wmt = lastverticalw[0];
		for( i=0; i<lgth1; i++ )
		{
			if( lastverticalw[i] >= wmt )
			{
				wmt = lastverticalw[i];
				MSF_IJP( lgth1, lgth2 ) = +( lgth1 - i );
			}
		}
		for( j=0; j<lgth2; j++ )
		{
			if( currentw[j] >= wmt )
			{
				wmt = currentw[j];
				MSF_IJP( lgth1, lgth2 ) = -( lgth2 - j );
			}
		}
	}
	for( i=0; i<lgth1+1; i++ ) MSF_IJP( i, 0 ) = i + 1;
	for( j=0; j<lgth2+1; j++ ) MSF_IJP( 0, j ) = -( j + 1 );

	iin = lgth1; jin = lgth2;
	klim = lgth1+lgth2;
	for( k=0; k<=klim; k++ )
	{
		int v = MSF_IJP( iin, jin );
		if( v < 0 )
		{
			ifi = iin-1; jfi = jin+v;
		}
		else if( v > 0 )
		{
			ifi = iin-v; jfi = jin-1;
		}
		else
		{
			ifi = iin-1; jfi = jin-1;
		}
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
		if( iin <= 0 || jin <= 0 ) break;
		*--gaptable1 = 'o';
		*--gaptable2 = 'o';
		k++;
		iin = ifi; jin = jfi;
	}
#undef MSF_IJP
	strcpy( t1, gaptable1 );
	strcpy( t2, gaptable2 );
	if( ijpb != msf_ijpblk ) free( ijpb );
	free( gt1bk ); free( gt2bk );
	free( w1 ); free( w2 ); free( initverticalw ); free( lastverticalw ); free( m ); free( mp );
	return( wm );
}

/* ---- the recursion ---- */

/* A node of the divide and conquer.  Its gap tables are: [left child][fillers][right child] for an
   inner node, or its own tables for a leaf (tanni, or lgth2 <= 0). */
typedef struct msnode_s
{
	int ist, ien, jst, jen;
	int leaf;
	int left, right;         // children (node indices), inner nodes
	int fill1, fill2;        // filler lengths (jmid-jumpj-1 gap columns in group 1, imid-jumpi-1 in group 2), or <= 0
	double pen1, pen2;       // the score terms of the fillers
	char *t1, *t2;           // leaf tables
	int tlen;                // leaf table length
	double value;            // leaf: its score; inner: filled in at the end
} msnode_t;

typedef struct
{
	const msf_t *c;
	msnode_t *nodes;
	int nnodes, capnodes;
	pthread_mutex_t mu;
	int *stack; int nstack, capstack; // nodes waiting to be processed
	int active;                       // nodes being processed (guarded by mu)
	int par;                          // threads may take part
} mstree_t;

static int mst_newnode( mstree_t *t, int ist, int ien, int jst, int jen ) // with t->mu held if par
{
	msnode_t *n;
	if( t->nnodes == t->capnodes )
	{
		t->capnodes = t->capnodes ? 2 * t->capnodes : 64;
		t->nodes = realloc( t->nodes, sizeof( msnode_t ) * t->capnodes );
		if( t->nodes == NULL ) ErrorExit( "Cannot allocate in MSalignmm (work/v10w).\n" );
	}
	n = t->nodes + t->nnodes;
	memset( n, 0, sizeof( msnode_t ) );
	n->ist = ist; n->ien = ien; n->jst = jst; n->jen = jen;
	n->left = n->right = -1;
	return( t->nnodes++ );
}

static void msf_sweeps( mstree_t *t, msf_sweep_t *w, int nthreadshint );

static int mskcheck_v = -1;
static long mskcheck_n = 0;
static void mskcheck_report( void ) __attribute__(( destructor ));
static void mskcheck_report( void )
{
	if( mskcheck_v > 0 ) reporterr( "\nMSKCHECK: %ld nodes' sweeps verified against the original loops, no mismatch\n", mskcheck_n );
}
static int mskcheck( void )
{
	if( mskcheck_v < 0 ) mskcheck_v = ( getenv( "MAFFT_MSKCHECK" ) != NULL );
	return( mskcheck_v );
}

/* MAFFT_MSKCHECK=1: the sweeps of a node by the original code (MSalignmm_rec's two loops, extracted
   with unifdef for this configuration: STOREWM 0, USE_PENALTY_EX 0, MEMSAVE 1, FASTMATCHCALC 1, no
   chudan; the base cases and the recursion removed), to compare every node's mid arrays and split
   with the fast sweeps' */
typedef struct { double *midw, *midm, *midn; int *jumpbacki, *jumpbackj; int imid, jmid, jumpi, jumpj; } msf_refout_t;
static void msf_ref_sweeps( double **n_dynamicmtx, double **cpmx1pt, double **cpmx2pt, int ist, int ien, int jst, int jen, double **gapinfo, double headgapfreq1_g, double headgapfreq2_g, msf_refout_t *out )
{
	int depth = 0;
//	int k;
	int alnlen;
	double value = 0.0;
	register int i, j;
	char **aseq1, **aseq2, *agt1, *agt2;
	int ll1, ll2, l, len;
	int lasti, lastj, imid;
	int jmid = 0; // by D.Mathog, a guess
	double wm = 0.0;   /* int ?????? */
	double g;
	double *currentw, *previousw;
	double *wtmp;
//	short *ijppt;
	int *mpjpt;
//	short **ijp;
	int *mp;
	int mpi;
	double *mjpt, *prept, *curpt;
	double mi;
	double *m;
	double *w1, *w2;
//	double *match;
	double *initverticalw;    /* kufuu sureba iranai */
	double *lastverticalw;    /* kufuu sureba iranai */
	int **intwork;
	double **doublework;
//	short **shortmtx;
	double *midw;
	double *midm;
	double *midn;
	int lgth1, lgth2;
	double maxwm;
	int *jumpforwi;
	int *jumpforwj;
	int *jumpbacki;
	int *jumpbackj;
	int *jumpdummi; //muda
	int *jumpdummj = NULL; // by D.Mathog, a guess
	int jumpi, jumpj = 0; // by D.Mathog, a guess
	char *gaps;
	int ijpi, ijpj;
	double *ogcp1;
	double *fgcp1;
	double *ogcp2;
	double *fgcp2;
	double firstm;
	int firstmp;
	double *gapfreq1f;
	double *gapfreq2f;
	double headgapfreq1;
	double headgapfreq2;

#if 0
	int nglen1, nglen2;
	nglen1 = seqlen( seq1[0] );
	nglen2 = seqlen( seq2[0] );
#endif

//	fprintf( stderr, "fulllen1 = %d, fulllen2 = %d, headgp = %d, tailgp = %d\n", fulllen1, fulllen2, headgp, tailgp );

	ogcp1 = gapinfo[0] + ist;
	fgcp1 = gapinfo[1] + ist;
	ogcp2 = gapinfo[2] + jst;
	fgcp2 = gapinfo[3] + jst;
	gapfreq1f = gapinfo[4] + ist;
	gapfreq2f = gapinfo[5] + jst;

	if( ist > 0 ) headgapfreq1 = gapfreq1f[-1];
	else headgapfreq1 = headgapfreq1_g;
	if( jst > 0 ) headgapfreq2 = gapfreq2f[-1];
	else headgapfreq2 = headgapfreq2_g;

	depth++;
	reccycle++;

	lgth1 = ien-ist+1;
	lgth2 = jen-jst+1;

//	if( lgth1 < 5 )
//		fprintf( stderr, "\nWARNING: lgth1 = %d\n", lgth1 );
//	if( lgth2 < 5 )
//		fprintf( stderr, "\nWARNING: lgth2 = %d\n", lgth2 );
//


	ll1 = ( (int)(lgth1) ) + 100;
	ll2 = ( (int)(lgth2) ) + 100;

//	fprintf( stderr, "ll1,ll2=%d,%d\n", ll1, ll2 );

	w1 = AllocateFloatVec( ll2+2 );
	w2 = AllocateFloatVec( ll2+2 );
//	match = AllocateFloatVec( ll2+2 );
	midw = AllocateFloatVec( ll2+2 );
	midn = AllocateFloatVec( ll2+2 );
	midm = AllocateFloatVec( ll2+2 );
	jumpbacki = AllocateIntVec( ll2+2 );
	jumpbackj = AllocateIntVec( ll2+2 );
	jumpforwi = AllocateIntVec( ll2+2 );
	jumpforwj = AllocateIntVec( ll2+2 );
	jumpdummi = AllocateIntVec( ll2+2 );
	jumpdummj = AllocateIntVec( ll2+2 );

	initverticalw = AllocateFloatVec( ll1+2 );
	lastverticalw = AllocateFloatVec( ll1+2 );

	m = AllocateFloatVec( ll2+2 );
	mp = AllocateIntVec( ll2+2 );


	doublework = AllocateFloatMtx( MAX( ll1, ll2 )+2, nalphabets ); 
	intwork = AllocateIntMtx( MAX( ll1, ll2 )+2, nalphabets ); 


#if 0
	shortmtx = AllocateShortMtx( ll1, ll2 );


	ijp = shortmtx;
#endif

	currentw = w1;
	previousw = w2;

//	match_calc( n_dynamicmtx, initverticalw, cpmx2+jst, cpmx1+ist, 0, lgth1, doublework, intwork, 1 );
	match_calc_alphabet_seq( n_dynamicmtx, initverticalw, cpmx2pt, cpmx1pt, jst, ist, lgth1, doublework, intwork, 1 );

//	match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, 0, lgth2, doublework, intwork, 1 );
	match_calc_alphabet_seq( n_dynamicmtx, currentw, cpmx1pt, cpmx2pt, ist, jst, lgth2, doublework, intwork, 1 );

	for( i=1; i<lgth1+1; i++ )
	{
//		initverticalw[i] += ( ogcp1[0] + fgcp1[i-1] ) ;
		initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2f[0] ) ;
	}
	for( j=1; j<lgth2+1; j++ )
	{
		currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1f[0]) ;
	}



	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2f[j-1];
//		m[j] = currentw[j-1];
		mp[j] = 0;
	}

	lastverticalw[0] = currentw[lgth2-1];

	imid = lgth1 * 0.5;

	jumpi = 0; // atode kawaru.
	lasti = lgth1+1;
	for( i=1; i<=imid; i++ )
	{
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

//		match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, i, lgth2, doublework, intwork, 0 );
		match_calc_alphabet_seq( n_dynamicmtx, currentw, cpmx1pt, cpmx2pt, ist+i, jst, lgth2, doublework, intwork, 0 );
		currentw[0] = initverticalw[i];

		m[0] = ogcp1[i];
		if( i == imid ) midm[0] = m[0];

		mi = previousw[0] + ogcp2[1] * gapfreq1f[i-1]; 
//		mi = previousw[0];
		mpi = 0;


//		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;


		lastj = lgth2+1;
		for( j=1; j<lastj; j++ )
		{

			wm = *prept;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
			g = mi + fgcp2[j-1] * gapfreq1f[i];
//			g = mi + fpenalty;
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
//				*ijppt = -( j - mpi );
			}
			g = *prept + ogcp2[j] * gapfreq1f[i-1];
//			g = *prept;
			if( g >= mi )
			{
				mi = g;
				mpi = j-1;
			}

			g = *mjpt + fgcp1[i-1] * gapfreq2f[j];
//			g = *mjpt + fpenalty;
#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
//				*ijppt = +( i - *mpjpt );
			}


			g = *prept + ogcp1[i] * gapfreq2f[j-1];
//			g = *prept;
			if( g >= *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}

#if 0
			fprintf( stderr, "%5.0f ", wm );
#endif
			*curpt += wm;


			if( i == imid ) //muda
			{	
				jumpbackj[j] = *mpjpt; // muda atode matomeru
				jumpbacki[j] = mpi; // muda atode matomeru
//				fprintf( stderr, "jumpbackj[%d] in forward dp is %d\n", j, *mpjpt );
//				fprintf( stderr, "jumpbacki[%d] in forward dp is %d\n", j, mpi );
				midw[j] = *curpt;
				midm[j] = *mjpt;
				midn[j] = mi;
			}

//			fprintf( stderr, "m[%d] = %f\n", j, m[j] );
			mjpt++;
			prept++;
			mpjpt++;
			curpt++;

		}
		lastverticalw[i] = currentw[lgth2-1];


#if 0  // ue
		if( i == imid )
		{
			for( j=0; j<lgth2; j++ ) midw[j] = currentw[j];
			for( j=0; j<lgth2; j++ ) midm[j] = m[j];
		}
#endif
	}
//	for( j=0; j<lgth2; j++ ) midw[j] = WMMTX[imid][j];
//	for( j=0; j<lgth2; j++ ) midm[j] = WMMTX2[imid][j];

#if 0
    for( i=0; i<lgth1; i++ )
    {
        for( j=0; j<lgth2; j++ )
        {
            fprintf( stderr, "% 10.2f ", WMMTX[i][j] );
        }
        fprintf( stderr, "\n" );
    }
	fprintf( stderr, "\n" );
	fprintf( stderr, "WMMTX2 = \n" );
    for( i=0; i<lgth1; i++ )
    {
        for( j=0; j<lgth2; j++ )
        {
            fprintf( stderr, "% 10.2f ", WMMTX2[i][j] );
        }
        fprintf( stderr, "\n" );
    }
	fprintf( stderr, "\n" );
#endif

// gyakudp

//	match_calc( n_dynamicmtx, initverticalw, cpmx2+jst, cpmx1+ist, lgth2-1, lgth1, doublework, intwork, 1 );
	match_calc_alphabet_seq( n_dynamicmtx, initverticalw, cpmx2pt, cpmx1pt, jst+lgth2-1, ist, lgth1, doublework, intwork, 1 );

//	match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, lgth1-1, lgth2, doublework, intwork, 1 );
	match_calc_alphabet_seq( n_dynamicmtx, currentw, cpmx1pt, cpmx2pt, ist+lgth1-1, jst, lgth2, doublework, intwork, 1 );

	for( i=0; i<lgth1-1; i++ )
	{
//		initverticalw[i] += ( fgcp1[lgth1-1] + ogcp1[i+1] );
		initverticalw[i] += ( fgcp1[lgth1-1] * gapfreq2f[lgth2] + ogcp1[i+1] * gapfreq2f[lgth2-1] );
	}
	for( j=0; j<lgth2-1; j++ )
	{
//		currentw[j] += ( fgcp2[lgth2-1] + ogcp2[j+1] );
		currentw[j] += ( fgcp2[lgth2-1] * gapfreq1f[lgth1] + ogcp2[j+1] * gapfreq1f[lgth1-1] );
	}







#if 0
	for( j=lgth2-1; j>0; --j )
	{
		m[j-1] = currentw[j] + fgcp2[lgth2-2];
//		m[j-1] = currentw[j];
		mp[j] = lgth1-1;
	}
#else
	for( j=lgth2-1; j>-1; --j )
	{
		m[j] = currentw[j+1] + fgcp1[lgth1-2] * gapfreq2f[j+1];
//		m[j-1] = currentw[j];
		mp[j] = lgth1-1;
	}
#endif

//	for( j=0; j<lgth2; j++ ) m[j] = 0.0;
	// m[lgth2-1] ha irunoka?


//	for( i=lgth1-2; i>=imid; i-- )
	firstm = -9999999.9;
//	firstmp = lgth1-1;
	firstmp = lgth1;
	for( i=lgth1-2; i>-1; i-- )
	{
		wtmp = previousw;
		previousw = currentw;
		currentw = wtmp;
		previousw[lgth2-1] = initverticalw[i+1];
//		match_calc( currentw, seq1, seq2, i, lgth2 );
//		match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, i, lgth2, doublework, intwork, 0 );
		match_calc_alphabet_seq( n_dynamicmtx, currentw, cpmx1pt, cpmx2pt, ist+i, jst, lgth2, doublework, intwork, 0 );

		currentw[lgth2-1] = initverticalw[i];

//		m[lgth2] = fgcp1[i];
//		WMMTX2[i][lgth2] += m[lgth2];
//		fprintf( stderr, "m[] = %f\n", m[lgth2] );

		mi = previousw[lgth2-1] + fgcp2[lgth2-2] * gapfreq1f[i+1];
//		mi = previousw[lgth2-1];
		mpi = lgth2 - 1;

		mjpt = m + lgth2 - 2;
		prept = previousw + lgth2 - 1;
		curpt = currentw + lgth2 - 2;
		mpjpt = mp + lgth2 - 2;


		for( j=lgth2-2; j>-1; j-- )
		{
			wm = *prept;
			ijpi = i+1;
			ijpj = j+1;

			g = mi + ogcp2[j+1] * gapfreq1f[i];
//			g = mi + fpenalty;
			if( g > wm )
			{
				wm = g;
				ijpj = mpi;
				ijpi = i+1;
			}

			g = *prept + fgcp2[j] * gapfreq1f[i+1];
//			g = *prept;
			if( g >= mi )
			{
//				fprintf( stderr, "i,j=%d,%d - renewed! mpi = %d\n", i, j, j+1 );
				mi = g;
				mpi = j + 1;
			}


//			fprintf( stderr, "i,j=%d,%d *mpjpt = %d\n", i, j, *mpjpt );
			g = *mjpt + ogcp1[i+1] * gapfreq2f[j];
//			g = *mjpt + fpenalty;
			if( g > wm )
			{
				wm = g;
				ijpi = *mpjpt;
				ijpj = j+1;
			}

//			if( i == imid )fprintf( stderr, "i,j=%d,%d \n", i, j );
			g = *prept + fgcp1[i] * gapfreq2f[j+1];
//			g = *prept;
			if( g >= *mjpt )
			{
				*mjpt = g;
				*mpjpt = i + 1;
			}


			if( i == jumpi || i == imid - 1 )
			{
				jumpforwi[j] = ijpi; //muda
				jumpforwj[j] = ijpj; //muda
//				fprintf( stderr, "jumpfori[%d] = %d\n", j, ijpi );
//				fprintf( stderr, "jumpforj[%d] = %d\n", j, ijpj );
			}
			if( i == imid ) // muda
			{
				midw[j] += wm;
//				midm[j+1] += *mjpt + fpenalty; //??????
				midm[j+1] += *mjpt; //??????
			}
			if( i == imid - 1 )
			{
//				midn[j] += mi + fpenalty;  //????
				midn[j] += mi;  //????
			}
			*curpt += wm;

			mjpt--;
			prept--;
			mpjpt--;
			curpt--;
		}
//		fprintf( stderr, "adding *mjpt (=%f) to WMMTX2[%d][%d]\n", *mjpt, i, j+1 );
		g = *prept + fgcp1[i];
		if( firstm < g ) 
		{
			firstm = g;
			firstmp = i + 1;
		}
		if( i == imid ) midm[j+1] += firstm;


		if( i == imid - 1 )	
		{
			maxwm = midw[1];
			jmid = 0;
//			if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			for( j=2; j<lgth2-1; j++ )
			{
				wm = midw[j];
				if( wm > maxwm )
				{
					jmid = j;
					maxwm = wm;
				}
//				if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			}
			for( j=0; j<lgth2+1; j++ )
			{
				wm = midm[j];
				if( wm > maxwm )
				{
					jmid = j;
					maxwm = wm;
				}
//				if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			}

//			if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );


//			fprintf( stderr, "### imid=%d, jmid=%d\n", imid, jmid );
			wm = midw[jmid];
			jumpi = imid-1;
			jumpj = jmid-1;
			if( jmid > 0 && midn[jmid-1] > wm ) //060413
			{
				jumpi = imid-1;
				jumpj = jumpbacki[jmid];
				wm = midn[jmid-1];
//				fprintf( stderr, "rejump (n)\n" );
			}
			if( midm[jmid] > wm )
			{
				jumpi = jumpbackj[jmid];
				jumpj = jmid-1;
				wm = midm[jmid];
//				fprintf( stderr, "rejump (m) jumpi=%d\n", jumpi );
			}


//			fprintf( stderr, "--> imid=%d, jmid=%d\n", imid, jmid );
//			fprintf( stderr, "--> jumpi=%d, jumpj=%d\n", jumpi, jumpj );
//			fprintf( stderr, "maxwm = %f\n", maxwm );
		}
		if( i == jumpi ) //saki?
		{
//			fprintf( stderr, "#### FIRST i=%d, jumpi<imid=%d<%d, ist=%d, ien=%d, firstmp=%d, lgth1=%d\n", i, jumpi, imid, ist, ien, firstmp, lgth1 );
//			fprintf( stderr, "#### mark 1\n" );
//			fprintf( stderr, "imid, jumpi = %d,%d\n", imid, jumpi );
//			fprintf( stderr, "jmid, jumpj = %d,%d\n", jmid, jumpj );

			if( jmid == 0 )
			{
//				fprintf( stderr, "#### CHUI2!\n" );
				jumpj = 0; jmid = 1;
#if 0 // v6.823 made
				jumpi = firstmp-1;
				imid = firstmp;
#endif
#if 0
				jumpi = 0;
				imid = 1;
#else
//				if( 1 || firstmp > 100 ) // naze 100
				if( imid < firstmp-1 ) // naze 100
				{
					jumpi = firstmp;
					imid = firstmp+1;
				}
#if 0
				else
				{
					jumpi = 0;
					imid = 1;
				}
#endif
#endif
			}

#if 0
			else if( jmid == lgth2 ) 
			{
				fprintf( stderr, "CHUI1!\n" );
				jumpi=0; jumpj=0;
				imid=jumpforwi[0]; jmid=lgth2-1;
			}
#else // 060414
			else if( jmid >= lgth2 ) 
			{
//				fprintf( stderr, "CHUI1!\n" );
				jumpi=imid-1; jmid=lgth2;
				jumpj = lgth2-1;
			}
#endif
			else
			{
//				fprintf( stderr, "#### CHUI3!\n" );
				imid = jumpforwi[jumpj];
				jmid = jumpforwj[jumpj];
				if( imid == jumpi ) jumpi = imid-1;
			}
#if 0
			fprintf( stderr, "jumpi -> %d\n", jumpi );
			fprintf( stderr, "jumpj -> %d\n", jumpj );
			fprintf( stderr, "imid -> %d\n", imid );
			fprintf( stderr, "jmid -> %d\n", jmid );
#endif
//			fprintf( stderr, "#### FINAL i=%d, jumpi<imid=%d<%d, ist=%d, ien=%d, firstmp=%d\n", i, jumpi, imid, ist, ien, firstmp );

			break;
		}
	}
	out->midw = midw; out->midm = midm; out->midn = midn; out->jumpbacki = jumpbacki; out->jumpbackj = jumpbackj;
	out->imid = imid; out->jmid = jmid; out->jumpi = jumpi; out->jumpj = jumpj;
	FreeFloatVec( w1 ); FreeFloatVec( w2 ); FreeFloatVec( initverticalw ); FreeFloatVec( lastverticalw );
	FreeIntVec( jumpforwi ); FreeIntVec( jumpforwj ); FreeIntVec( jumpdummi ); FreeIntVec( jumpdummj );
	FreeFloatVec( m ); FreeIntVec( mp ); FreeFloatMtx( doublework ); FreeIntMtx( intwork );
}


/* process one node: a leaf is solved; an inner node is split, and its children are returned */
static void mst_process( mstree_t *t, int id, int *kids, int nthreadshint )
{
	const msf_t *c = t->c;
	int ist, ien, jst, jen, lgth1, lgth2;
	if( t->par ) pthread_mutex_lock( &t->mu );
	ist = t->nodes[id].ist; ien = t->nodes[id].ien; jst = t->nodes[id].jst; jen = t->nodes[id].jen;
	if( t->par ) pthread_mutex_unlock( &t->mu );
	lgth1 = ien-ist+1;
	lgth2 = jen-jst+1;
	kids[0] = kids[1] = -1;
	if( lgth2 <= 0 )
	{
		int n = lgth1 > 0 ? lgth1 : 0, q;
		char *t1 = msf_alloc( n + 1 ), *t2 = msf_alloc( n + 1 );
		for( q=0; q<n; q++ ) { t1[q] = 'o'; t2[q] = '-'; }
		t1[n] = t2[n] = 0;
		if( t->par ) pthread_mutex_lock( &t->mu );
		t->nodes[id].leaf = 1; t->nodes[id].t1 = t1; t->nodes[id].t2 = t2; t->nodes[id].tlen = n; t->nodes[id].value = 0.0;
		if( t->par ) pthread_mutex_unlock( &t->mu );
		return;
	}
	if( lgth1 < DPTANNI || lgth2 < DPTANNI )
	{
		char *t1 = msf_alloc( ( lgth1 > 0 ? lgth1 : 0 ) + lgth2 + 2 ), *t2 = msf_alloc( ( lgth1 > 0 ? lgth1 : 0 ) + lgth2 + 2 );
		double v = msf_tanni( c, ist, ien, jst, jen, t1, t2 );
		if( t->par ) pthread_mutex_lock( &t->mu );
		t->nodes[id].leaf = 1; t->nodes[id].t1 = t1; t->nodes[id].t2 = t2; t->nodes[id].tlen = strlen( t1 ); t->nodes[id].value = v;
		if( t->par ) pthread_mutex_unlock( &t->mu );
		return;
	}
	{
		msf_sweep_t w;
		int imid, jmid, jumpi, jumpj, l;
		const double *ogcp1 = c->gapinfo[0] + ist, *fgcp1 = c->gapinfo[1] + ist;
		const double *ogcp2 = c->gapinfo[2] + jst, *fgcp2 = c->gapinfo[3] + jst;
		memset( &w, 0, sizeof( w ) );
		w.c = c; w.ist = ist; w.ien = ien; w.jst = jst; w.jen = jen; w.lgth1 = lgth1; w.lgth2 = lgth2;
		w.imid = lgth1 * 0.5;
		w.midw = msf_calloc( lgth2+2, sizeof( double ) );
		w.midm = msf_calloc( lgth2+2, sizeof( double ) );
		w.midn = msf_calloc( lgth2+2, sizeof( double ) );
		w.jumpbacki = msf_calloc( lgth2+2, sizeof( int ) );
		w.jumpbackj = msf_calloc( lgth2+2, sizeof( int ) );
		w.fwddone = 0;
		msf_sweeps( t, &w, nthreadshint );
		imid = w.imidout; jmid = w.jmid; jumpi = w.jumpi; jumpj = w.jumpj;
		if( mskcheck() ) // the original loops for this node, compared bitwise
		{
			msf_refout_t ro;
			msf_ref_sweeps( c->n_dynamicmtx, c->cpmx1, c->cpmx2, ist, ien, jst, jen, c->gapinfo, c->hgf1g, c->hgf2g, &ro );
			if( ro.imid != imid || ro.jmid != jmid || ro.jumpi != jumpi || ro.jumpj != jumpj
			 || memcmp( ro.midw, w.midw, sizeof( double ) * ( lgth2+1 ) ) || memcmp( ro.midm, w.midm, sizeof( double ) * ( lgth2+1 ) )
			 || memcmp( ro.midn, w.midn, sizeof( double ) * ( lgth2+1 ) ) || memcmp( ro.jumpbacki, w.jumpbacki, sizeof( int ) * ( lgth2+1 ) )
			 || memcmp( ro.jumpbackj, w.jumpbackj, sizeof( int ) * ( lgth2+1 ) ) )
			{
				reporterr( "\nMSKCHECK MISMATCH at node %d-%d x %d-%d: split %d %d %d %d vs %d %d %d %d\n", ist, ien, jst, jen, imid, jmid, jumpi, jumpj, ro.imid, ro.jmid, ro.jumpi, ro.jumpj );
				exit( 1 );
			}
			FreeFloatVec( ro.midw ); FreeFloatVec( ro.midm ); FreeFloatVec( ro.midn ); FreeIntVec( ro.jumpbacki ); FreeIntVec( ro.jumpbackj );
			__atomic_add_fetch( &mskcheck_n, 1, __ATOMIC_RELAXED );
		}
		free( w.midw ); free( w.midm ); free( w.midn ); free( w.jumpbacki ); free( w.jumpbackj );

		if( t->par ) pthread_mutex_lock( &t->mu );
		{
			int a = mst_newnode( t, ist, ist+jumpi, jst, jst+jumpj );
			int b = mst_newnode( t, ist+imid, ien, jst+jmid, jen );
			msnode_t *n = t->nodes + id;
			n->left = a; n->right = b;
			l = jmid - jumpj - 1;
			n->fill1 = l;
			if( l > 0 ) n->pen1 = ( ogcp2[jumpj+1] + fgcp2[jmid-1] );
			l = imid - jumpi - 1;
			n->fill2 = l;
			if( l > 0 ) n->pen2 = ( ogcp1[jumpi+1] + fgcp1[imid-1] );
			kids[0] = a; kids[1] = b;
		}
		if( t->par ) pthread_mutex_unlock( &t->mu );
	}
}

/* the score and the length of a subtree, in the original order: value = left; value += pen1;
   value += pen2; value += right */
static double mst_value( mstree_t *t, int id, size_t *len )
{
	msnode_t *n = t->nodes + id;
	double value;
	size_t a, b;
	if( n->leaf ) { *len = n->tlen; return( n->value ); }
	value = mst_value( t, n->left, &a );
	if( n->fill1 > 0 ) { value += n->pen1; a += n->fill1; }
	if( n->fill2 > 0 ) { value += n->pen2; a += n->fill2; }
	value += mst_value( t, n->right, &b );
	*len = a + b;
	return( value );
}
static void mst_tables( mstree_t *t, int id, char *g1, char *g2, size_t *pos )
{
	msnode_t *n = t->nodes + id;
	int q;
	if( n->leaf ) { memcpy( g1 + *pos, n->t1, n->tlen ); memcpy( g2 + *pos, n->t2, n->tlen ); *pos += n->tlen; return; }
	mst_tables( t, n->left, g1, g2, pos );
	for( q=0; q<n->fill1; q++ ) { g1[*pos] = '-'; g2[*pos] = 'o'; (*pos)++; }
	for( q=0; q<n->fill2; q++ ) { g1[*pos] = 'o'; g2[*pos] = '-'; (*pos)++; }
	mst_tables( t, n->right, g1, g2, pos );
}

/* ---- running the tree: serially, or with idle helper threads ---- */
#define MSF_PARCELLS 1.0e6   // a node's two sweeps run on two threads from this size (lgth1 x lgth2)

typedef struct { int type; int id; msf_phase_t *P; } mstask_t; // type 0: process node id; 2: take bands of phase P
typedef struct msjob_s { dphelpjob_t job; struct msjob_s *next; } msjob_t;

typedef struct
{
	mstree_t *t;
	mstask_t *tasks; int ntasks, cap;
	msjob_t *jobs;
	int done;          // the root and every node below it have been processed
	int outstanding;   // node tasks pushed or being processed (guarded by t->mu)
	int nhelped;       // tasks run by helper threads (statistics)
} msrun_t;

static TLS msrun_t *msf_run_cur = NULL; // the run a thread works for (to offer bands)

static void msr_push_locked( msrun_t *r, int type, int id, msf_phase_t *P )
{
	if( r->ntasks == r->cap )
	{
		r->cap = r->cap ? 2 * r->cap : 64;
		r->tasks = realloc( r->tasks, sizeof( mstask_t ) * r->cap );
		if( r->tasks == NULL ) ErrorExit( "Cannot allocate in MSalignmm (work/v10w).\n" );
	}
	r->tasks[r->ntasks].type = type; r->tasks[r->ntasks].id = id; r->tasks[r->ntasks].P = P;
	r->ntasks++;
}

static void msr_worker( void *arg );
static void msr_offer( msrun_t *r, int n ) // post n helper slots, if anybody is idle
{
#ifdef enablemultithread
	msjob_t *j;
	if( n <= 0 || dphelp_idle() <= 0 ) return;
	j = msf_alloc( sizeof( msjob_t ) );
	j->job.run = msr_worker;
	j->job.arg = r;
	pthread_mutex_lock( &r->t->mu );
	j->next = r->jobs; r->jobs = j;
	pthread_mutex_unlock( &r->t->mu );
	dphelp_post( &j->job, n );
#endif
}

#define MSF_BANDW 1024    // columns per band at least
#define MSF_BANDROWS 64   // a phase with fewer rows is not split into bands
static int msf_bandw = -1, msf_forcepar = -1; // test switches: MAFFT_MSBANDW (columns per band), MAFFT_MSFORCEPAR=1 (bands and tasks even without helpers)
static void msf_switches( void )
{
	if( msf_bandw < 0 )
	{
		char *e = getenv( "MAFFT_MSBANDW" );
		msf_forcepar = ( getenv( "MAFFT_MSFORCEPAR" ) != NULL );
		__atomic_store_n( &msf_bandw, ( e && atoi( e ) > 0 ) ? atoi( e ) : MSF_BANDW, __ATOMIC_RELEASE );
	}
}
static int msf_nbands( int ncol, int nrows )
{
	int nb;
	msf_switches();
	nb = ncol / msf_bandw;
	if( nrows < MSF_BANDROWS && !msf_forcepar ) return( 1 );
	if( nrows < 1 ) return( 1 );
	if( nb > 16 ) nb = 16;
	if( nb < 1 ) nb = 1;
	return( nb );
}

/* a node's two sweeps.  With helpers and a large node: the forward sweep's rows 1..imid-1 and the
   backward sweep's rows lgth1-2..imid+1 do not depend on each other; each is split into bands of
   columns (msf_phase_*), which this thread and any helpers take in order; then this thread does
   the forward sweep's row imid and the backward sweep's rows imid down to jumpi.  Every cell is
   computed by the same statements from the same operands as in the serial loops (a band's first
   cell reads its left/right neighbour's values from the row before, published per row). */
static void msf_sweeps( mstree_t *t, msf_sweep_t *w, int nthreadshint )
{
	msrun_t *r = msf_run_cur;
	msf_switches();
	if( r && t->par && ( (double)w->lgth1 * w->lgth2 >= MSF_PARCELLS || msf_forcepar ) )
	{
		msf_phase_t *Pf = msf_alloc( sizeof( msf_phase_t ) ), *Pb = msf_alloc( sizeof( msf_phase_t ) );
		int b, k, spins = 0, nbf, nbb, ntask;
		msf_forward_setup( w );
		msf_backward_setup( w );
		nbf = msf_nbands( w->lgth2, w->imid - 1 );
		nbb = msf_nbands( w->lgth2, w->lgth1 - 2 - w->imid );
		msf_phase_init( Pf, w, 0, 1, w->imid - 1, nbf, 1, w->lgth2 + 1, w->lgth1 + 1 );
		msf_phase_init( Pb, w, 1, w->lgth1 - 2, w->imid + 1, nbb, 0, w->lgth2 - 1, w->lgth1 + 1 );
		for( b=0; b<Pf->nb; b++ ) Pf->band[b].bw[0-Pf->band[b].off] = w->fw[0][Pf->band[b].hi-1];         // row 0, in w1
		for( b=0; b<Pb->nb; b++ ) Pb->band[b].bw[w->lgth1-1-Pb->band[b].off] = w->bw1[Pb->band[b].lo];    // row lgth1-1, in bw1
		ntask = Pf->nb + Pb->nb;
		pthread_mutex_lock( &t->mu );
		for( k=0; k<Pb->nb; k++ ) msr_push_locked( r, 2, -1, Pb );
		for( k=0; k<Pf->nb; k++ ) msr_push_locked( r, 2, -1, Pf );
		pthread_mutex_unlock( &t->mu );
		msr_offer( r, ntask - 1 );
		msf_phase_work( Pf );
		msf_phase_work( Pb );
		while( __atomic_load_n( &Pf->finished, __ATOMIC_ACQUIRE ) < Pf->nb || __atomic_load_n( &Pb->finished, __ATOMIC_ACQUIRE ) < Pb->nb ) msf_spinwait( &spins );
		pthread_mutex_lock( &t->mu ); // the tasks nobody took
		for( k=r->ntasks-1; k>=0; k-- )
			if( r->tasks[k].type == 2 && ( r->tasks[k].P == Pf || r->tasks[k].P == Pb ) )
			{
				memmove( r->tasks + k, r->tasks + k + 1, sizeof( mstask_t ) * ( r->ntasks - k - 1 ) );
				r->ntasks--;
			}
		pthread_mutex_unlock( &t->mu );
		spins = 0;
		while( __atomic_load_n( &Pf->refs, __ATOMIC_ACQUIRE ) > 0 || __atomic_load_n( &Pb->refs, __ATOMIC_ACQUIRE ) > 0 ) msf_spinwait( &spins );
		msf_forward_mid( w );
		msf_backward_afterbands( w );
		msf_backward_finish( w );
		msf_phase_free( Pf ); msf_phase_free( Pb );
		free( Pf ); free( Pb );
		return;
	}
	msf_forward( w );
	msf_backward_setup( w );
	msf_backward_rows( w );
	msf_backward_finish( w );
}

static void msr_loop( msrun_t *r, int poster )
{
	mstree_t *t = r->t;
	int spins = 0;
	msrun_t *saved = msf_run_cur;
	msf_run_cur = r;
	while( 1 )
	{
		mstask_t task;
		int kids[2], n, big;
		pthread_mutex_lock( &t->mu );
		if( r->ntasks == 0 )
		{
			int fin = ( r->outstanding == 0 );
			pthread_mutex_unlock( &t->mu );
			if( fin || !poster ) break; // helpers leave when there is nothing to take
			if( ++spins <= 64 ) ;
			else if( spins <= 4096 ) sched_yield();
			else { struct timespec ts; ts.tv_sec = 0; ts.tv_nsec = 50000; nanosleep( &ts, NULL ); }
			continue;
		}
		spins = 0;
		task = r->tasks[--r->ntasks];
		if( task.type == 2 ) __atomic_add_fetch( &task.P->refs, 1, __ATOMIC_ACQ_REL ); // (the phase lives until refs is 0)
		pthread_mutex_unlock( &t->mu );
		if( !poster ) __atomic_add_fetch( &r->nhelped, 1, __ATOMIC_RELAXED );
		if( task.type == 2 ) { msf_phase_work( task.P ); __atomic_sub_fetch( &task.P->refs, 1, __ATOMIC_RELEASE ); continue; }
		mst_process( t, task.id, kids, 0 );
		pthread_mutex_lock( &t->mu );
		n = 0; big = 0;
		for( int q=1; q>=0; q-- ) if( kids[q] >= 0 ) // the left child on top
		{
			msnode_t *nd = t->nodes + kids[q];
			msr_push_locked( r, 0, kids[q], NULL );
			r->outstanding++;
			n++;
			if( (double)( nd->ien - nd->ist + 1 ) * ( nd->jen - nd->jst + 1 ) >= MSF_PARCELLS ) big++;
		}
		r->outstanding--; // this node
		pthread_mutex_unlock( &t->mu );
		if( big > 1 ) msr_offer( r, big - 1 ); // this thread takes one of them itself
	}
	msf_run_cur = saved;
}
static void msr_worker( void *arg )
{
	msr_loop( (msrun_t *)arg, 0 );
}

/* the gap tables of the whole problem into mgt1/mgt2; returns the score */
static double msf_solve( msf_t *c, int lgth1, int lgth2, char *mgt1, char *mgt2 )
{
	mstree_t t;
	msrun_t r;
	int root, k;
	size_t len, pos = 0;
	double value;
	memset( &t, 0, sizeof( t ) );
	memset( &r, 0, sizeof( r ) );
	t.c = c;
#ifdef PROGTIMING
	double st0 = pt_now();
	static int msstat = -1;
	if( msstat < 0 ) msstat = ( getenv( "MAFFT_MSSTAT" ) != NULL );
#endif
#ifdef enablemultithread
	msf_switches();
	t.par = ( dphelppool != NULL || msf_forcepar );
#else
	t.par = 0;
#endif
	if( t.par ) pthread_mutex_init( &t.mu, NULL );
	root = mst_newnode( &t, 0, lgth1-1, 0, lgth2-1 );
	r.t = &t;
	if( t.par )
	{
		msr_push_locked( &r, 0, root, NULL );
		r.outstanding = 1;
		msr_loop( &r, 1 );
		while( r.jobs ) // every helper has left
		{
			msjob_t *j = r.jobs;
			r.jobs = j->next;
#ifdef enablemultithread
			dphelp_finish( &j->job );
#endif
			free( j );
		}
	}
	else
	{
		/* serially: depth first, left before right (the order does not matter: nodes are independent) */
		int *st = msf_alloc( sizeof( int ) * 64 ), nst = 0, cst = 64, kids[2];
		msrun_t *saved = msf_run_cur;
		msf_run_cur = NULL;
		st[nst++] = root;
		while( nst > 0 )
		{
			int id = st[--nst];
			mst_process( &t, id, kids, 0 );
			if( nst + 2 > cst ) { cst *= 2; st = realloc( st, sizeof( int ) * cst ); if( st == NULL ) ErrorExit( "Cannot allocate in MSalignmm (work/v10w).\n" ); }
			if( kids[1] >= 0 ) st[nst++] = kids[1];
			if( kids[0] >= 0 ) st[nst++] = kids[0];
		}
		free( st );
		msf_run_cur = saved;
	}
#ifdef PROGTIMING
	if( msstat && pt_now() - st0 > 0.05 ) reporterr( "\nMSSTAT %d x %d: %.3f s, %d nodes, %d tasks by helpers, idle now %d\n", lgth1, lgth2, pt_now() - st0, t.nnodes, r.nhelped, dphelp_idle() );
#endif
	value = mst_value( &t, root, &len );
	mst_tables( &t, root, mgt1, mgt2, &pos );
	mgt1[pos] = 0;
	mgt2[pos] = 0;
	for( k=0; k<t.nnodes; k++ ) { free( t.nodes[k].t1 ); free( t.nodes[k].t2 ); }
	free( t.nodes );
	free( r.tasks );
	if( t.par ) pthread_mutex_destroy( &t.mu );
	return( value );
}

/* the rows: seq[i] becomes gapireru( seq[i], gt ), in place (backwards: the result is not shorter) */
typedef struct { char **seq; const char *gt; int len, src; char gapchar; } msrows_t;
static void msf_rows_body( void *vctx, int lo, int hi )
{
	msrows_t *x = (msrows_t *)vctx;
	int i, p, s;
	for( i=lo; i<hi; i++ )
	{
		char *row = x->seq[i];
		s = x->src - 1;
		row[x->len] = 0;
		for( p=x->len-1; p>=0; p-- )
		{
			if( x->gt[p] == '-' ) row[p] = x->gapchar;
			else row[p] = row[s--];
		}
	}
}
static void msf_rows( char **seq, int n, int srclen, const char *gt, char gapchar )
{
	msrows_t x;
	x.seq = seq; x.gt = gt; x.len = strlen( gt ); x.src = srclen; x.gapchar = gapchar;
#ifdef enablemultithread
	dphelp_for( n, MAX( 4, n / 64 ), (double)n * x.len, 2.0e6, msf_rows_body, &x );
#else
	msf_rows_body( &x, 0, n );
#endif
}

/* ---- the profiles of a call with boundary gap patterns (Falign_udpari_long's segments), fused and
   by column ranges with helpers: cpmx_calc_new, new_OpeningGapCount, new_FinalGapCount and gapcountf
   in one pass.  Every entry receives the same additions ("+= eff[j]") for the same members j, in
   ascending j, as in the four functions; only the columns are split among threads. */
typedef struct { char **seq; int clus, len; double *eff; char *sgap; double **cpmx; double *og, *fg, *gf; } msprof_t;
static void msf_profile_cols( void *vctx, int lo, int hi )
{
	msprof_t *x = (msprof_t *)vctx;
	int j, k, c;
	for( k=0; k<nalphabets; k++ ) for( c=lo; c<hi; c++ ) x->cpmx[k][c] = 0.0;
	for( c=lo; c<hi; c++ ) { x->og[c] = 0.0; x->fg[c] = 0.0; x->gf[c] = 0.0; }
	for( j=0; j<x->clus; j++ )
	{
		const char *s = x->seq[j];
		const double feff = (double)x->eff[j];
		int gb = ( lo == 0 ) ? ( x->sgap[j] == '-' ) : ( s[lo-1] == '-' );
		for( c=lo; c<hi; c++ )
		{
			const int gc = ( s[c] == '-' ), gn = ( s[c+1] == '-' );
			x->cpmx[(unsigned char)amino_n[(unsigned char)s[c]]][c] += feff;
			if( !gb *  gc ) x->og[c] += feff;
			if( gc * !gn ) x->fg[c] += feff;
			if( gc ) x->gf[c] += feff;
			gb = gc;
		}
	}
}
static void msf_profile_sgap( char **seq, int clus, double *eff, int len, char *sgap, double **cpmx, double *ogcp, double *fgcp, double *gapfreq )
{
	msprof_t x;
	x.seq = seq; x.clus = clus; x.len = len; x.eff = eff; x.sgap = sgap; x.cpmx = cpmx; x.og = ogcp; x.fg = fgcp; x.gf = gapfreq;
#ifdef enablemultithread
	dphelp_for( len, MAX( 512, len / 64 ), (double)clus * len, 1.0e6, msf_profile_cols, &x );
#else
	msf_profile_cols( &x, 0, len );
#endif
}

static void freearrays( 
	double *ogcp1, 
	double *ogcp2, 
	double *ogcp1o, 
	double *ogcp2o, 
	double *fgcp1, 
	double *fgcp2, 
	double *fgcp1o, 
	double *fgcp2o, 
	double **cpmx1,
	double **cpmx2, 
	double *gapfreq1f,
	double *gapfreq2f,
	double **gapinfo, 
	char **mseq1, 
	char **mseq2,
	char *mgt1, 
	char *mgt2 
)
{
	FreeFloatVec( ogcp1 );
	FreeFloatVec( ogcp2 );
	FreeFloatVec( ogcp1o );
	FreeFloatVec( ogcp2o );
	FreeFloatVec( fgcp1 );
	FreeFloatVec( fgcp2 );
	FreeFloatVec( fgcp1o );
	FreeFloatVec( fgcp2o );
	FreeFloatMtx( cpmx1 );
	FreeFloatMtx( cpmx2 );
	FreeFloatVec( gapfreq1f );
	FreeFloatVec( gapfreq2f );
	free( (void *)gapinfo );

	FreeCharMtx( mseq1 );
	FreeCharMtx( mseq2 );
	free( mgt1 );
	free( mgt2 );
}

double MSalignmm( double **n_dynamicmtx, char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen, char *sgap1, char *sgap2, char *egap1, char *egap2, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp, double ***cpmxchild0, double ***cpmxchild1, double ***cpmxresult, double orieff1, double orieff2 )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{

//	int k;
	int i, j;
	int ll1, ll2;
	int lgth1, lgth2;
	double wm = 0.0;   /* int ?????? */
	char **mseq1;
	char **mseq2;
	char *mgt1;
	char *mgt2;
	double *ogcp1, *ogcp1o, *ogcp1opt;
	double *ogcp2, *ogcp2o, *ogcp2opt;
	double *fgcp1, *fgcp1o, *fgcp1opt;
	double *fgcp2, *fgcp2o, *fgcp2opt;
	double **cpmx1, **cpmx1pt;
	double **cpmx2, **cpmx2pt;
	double **gapinfo;
	double fpenalty = (double)penalty;
	double *gapfreq1f, *gapfreq1pt;
	double *gapfreq2f, *gapfreq2pt;
	int nglen1, nglen2;
	double headgapfreq1;
	double headgapfreq2;
	int msfast = ( msf_on() && chudanpt == NULL && nalphabets <= MSF_MAXL ); // work/v10w: see above
	char **ckseq1 = NULL, **ckseq2 = NULL, *ckgt1 = NULL, *ckgt2 = NULL; // work/v10w: MAFFT_MSCHECK

#if 0
	fprintf( stderr, "eff in SA+++align\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "eff1[%d] = %f\n", i, eff1[i] );
#endif

	nglen1 = seqlen( seq1[0] );
	nglen2 = seqlen( seq2[0] );


#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "seq1[%d] at root = %s\n", i, seq1[i] );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, "seq2[%d] at root = %s\n", j, seq2[j] );
	fprintf( stderr, "\n" );
#endif

	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );


	ll1 = ( (int)(lgth1) ) + 100;
	ll2 = ( (int)(lgth2) ) + 100;

	mseq1 = AllocateCharMtx( icyc, msfast ? 0 : ll1+ll2 ); // work/v10w: the fast path makes the rows from the gap tables
	mseq2 = AllocateCharMtx( jcyc, msfast ? 0 : ll1+ll2 );
	mgt1 = AllocateCharVec( ll1+ll2 );
	mgt2 = AllocateCharVec( ll1+ll2 );

	gapinfo = AllocateFloatMtx( 6, 0 );
	ogcp1 = AllocateFloatVec( ll1+2 );
	ogcp2 = AllocateFloatVec( ll2+2 );
	fgcp1 = AllocateFloatVec( ll1+2 );
	fgcp2 = AllocateFloatVec( ll2+2 );
	ogcp1o = AllocateFloatVec( ll1+2 );
	ogcp2o = AllocateFloatVec( ll2+2 );
	fgcp1o = AllocateFloatVec( ll1+2 );
	fgcp2o = AllocateFloatVec( ll2+2 );

	cpmx1 = AllocateFloatMtx( nalphabets, ll1+2 );
	cpmx2 = AllocateFloatMtx( nalphabets, ll2+2 );

	gapfreq1f = AllocateFloatVec( ll1+2 ); // must be filled with 0.0
	gapfreq2f = AllocateFloatVec( ll2+2 ); // must be filled with 0.0

	for( i=0; i<icyc; i++ ) 
	{
		if( strlen( seq1[i] ) != lgth1 )
		{
			fprintf( stderr, "i = %d / %d\n", i, icyc );
			fprintf( stderr, "bug! hairetsu ga kowareta!\n" );
			exit( 1 );
		}
	}
	for( j=0; j<jcyc; j++ )
	{
		if( strlen( seq2[j] ) != lgth2 )
		{
			fprintf( stderr, "j = %d / %d\n", j, jcyc );
			fprintf( stderr, "bug! hairetsu ga kowareta!\n" );
			exit( 1 );
		}
	}

//	cpmx_calc_new( seq1, cpmx1, eff1, lgth1, icyc );
//	cpmx_calc_new( seq2, cpmx2, eff2, lgth2, jcyc );

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
			cpmx1pt = *cpmxchild0;
#if ATO
			gapfreq1pt = (*cpmxchild0)[nalphabets];
			ogcp1opt = (*cpmxchild0)[nalphabets+1];
			fgcp1opt = (*cpmxchild0)[nalphabets+2];
#endif
		}
		else
		{
//			reporterr( "\nDo not use cpmxhist for child 0!\n" );
			cpmx1pt = cpmx1;
			cpmx_calc_new( seq1, cpmx1pt, eff1, lgth1, icyc );
#if ATO
			gapfreq1pt = gapfreq1f;
			gapcountf( gapfreq1pt, seq1, icyc, eff1, lgth1 );
			for( i=0; i<lgth1; i++ ) gapfreq1pt[i] = 1.0 - gapfreq1pt[i];

			ogcp1opt = ogcp1o;
			fgcp1opt = fgcp1o;
			st_OpeningGapCount( ogcp1opt, icyc, seq1, eff1, lgth1 );
			st_FinalGapCount( fgcp1opt, icyc, seq1, eff1, lgth1 );
#endif
		}
	
		if( cpmxchild1 && *cpmxchild1 )
		{
//			reporterr( "\nUse cpmxhist for child 1!\n" );
			cpmx2pt = *cpmxchild1;
#if ATO
			gapfreq2pt = (*cpmxchild1)[nalphabets];
			ogcp2opt = (*cpmxchild1)[nalphabets+1];
			fgcp2opt = (*cpmxchild1)[nalphabets+2];
#endif
		}
		else
		{
//			reporterr( "\nDo not use cpmxhist for child 1!\n" );
			cpmx2pt = cpmx2;
			cpmx_calc_new( seq2, cpmx2pt, eff2, lgth2, jcyc );
#if ATO
			gapfreq2pt = gapfreq2f;
			gapcountf( gapfreq2pt, seq2, jcyc, eff2, lgth2 );
			for( i=0; i<lgth2; i++ ) gapfreq2pt[i] = 1.0 - gapfreq2pt[i];

			ogcp2opt = ogcp2o;
			fgcp2opt = fgcp2o;
			st_OpeningGapCount( ogcp2opt, jcyc, seq2, eff2, lgth2 );
			st_FinalGapCount( fgcp2opt, jcyc, seq2, eff2, lgth2 );
#endif
		}

		//sgap1 == 1 ni taiou shite iani node, legacygap niyorazu, headgap to tailgap ha nai to minasu.
#if ATO
		headgapfreq1 = 1.0;
		headgapfreq2 = 1.0;
		gapfreq1pt[lgth1] = 1.0;
		gapfreq2pt[lgth2] = 1.0;
#endif
	}
	else
	{
		cpmx1pt = cpmx1;
		cpmx2pt = cpmx2;
#if ATO
		gapfreq1pt = gapfreq1f;
		gapfreq2pt = gapfreq2f;
		ogcp1opt = ogcp1o;
		ogcp2opt = ogcp2o;
		fgcp1opt = fgcp1o;
		fgcp2opt = fgcp2o;
#endif

		int msprof = ( msfast && sgap1 && legacygapcost == 0 ); // work/v10w: the four profile passes fused, by columns (msf_profile_sgap)
		if( msprof )
		{
			msf_profile_sgap( seq1, icyc, eff1, lgth1, sgap1, cpmx1pt, ogcp1opt, fgcp1opt, gapfreq1pt );
			msf_profile_sgap( seq2, jcyc, eff2, lgth2, sgap2, cpmx2pt, ogcp2opt, fgcp2opt, gapfreq2pt );
			if( mscheck() ) // the four original functions into fresh arrays, compared bitwise
			{
				int q, k;
				for( q=0; q<2; q++ )
				{
					char **sq = q ? seq2 : seq1, *sg = q ? sgap2 : sgap1, *eg = q ? egap2 : egap1;
					int n = q ? jcyc : icyc, len = q ? lgth2 : lgth1;
					double *ef = q ? eff2 : eff1, **cp = q ? cpmx2pt : cpmx1pt;
					double **rc = AllocateFloatMtx( nalphabets, len+2 ), *ro = AllocateFloatVec( len+2 ), *rf = AllocateFloatVec( len+2 ), *rg = AllocateFloatVec( len+2 );
					cpmx_calc_new( sq, rc, ef, len, n );
					new_OpeningGapCount( ro, n, sq, ef, len, sg );
					new_FinalGapCount( rf, n, sq, ef, len, eg );
					gapcountf( rg, sq, n, ef, len );
					for( k=0; k<nalphabets; k++ ) if( memcmp( rc[k], cp[k], sizeof( double ) * len ) ) { reporterr( "\nMSCHECK MISMATCH: profile, letter %d\n", k ); exit( 1 ); }
					if( memcmp( ro, q ? ogcp2opt : ogcp1opt, sizeof( double ) * len ) || memcmp( rf, q ? fgcp2opt : fgcp1opt, sizeof( double ) * len ) || memcmp( rg, q ? gapfreq2pt : gapfreq1pt, sizeof( double ) * len ) )
						{ reporterr( "\nMSCHECK MISMATCH: gap counts\n" ); exit( 1 ); }
					FreeFloatMtx( rc ); FreeFloatVec( ro ); FreeFloatVec( rf ); FreeFloatVec( rg );
				}
			}
		}
		else
		{
		cpmx_calc_new( seq1, cpmx1pt, eff1, lgth1, icyc );
		cpmx_calc_new( seq2, cpmx2pt, eff2, lgth2, jcyc );
		}

#if ATO
		if( sgap1 )
		{
			if( !msprof ) // work/v10w
			{
			new_OpeningGapCount( ogcp1opt, icyc, seq1, eff1, lgth1, sgap1 );
			new_FinalGapCount( fgcp1opt, icyc, seq1, eff1, lgth1, egap1 );
			new_OpeningGapCount( ogcp2opt, jcyc, seq2, eff2, lgth2, sgap2 );
			new_FinalGapCount( fgcp2opt, jcyc, seq2, eff2, lgth2, egap2 );
			}
	
			outgapcount( &headgapfreq1, icyc, sgap1, eff1 );
			outgapcount( &headgapfreq2, jcyc, sgap2, eff2 );
			outgapcount( gapfreq1pt+lgth1, icyc, egap1, eff1 );
			outgapcount( gapfreq2pt+lgth2, jcyc, egap2, eff2 );
		}
		else
		{
			st_OpeningGapCount( ogcp1opt, icyc, seq1, eff1, lgth1 );
			st_FinalGapCount( fgcp1opt, icyc, seq1, eff1, lgth1 );
			st_OpeningGapCount( ogcp2opt, jcyc, seq2, eff2, lgth2 );
			st_FinalGapCount( fgcp2opt, jcyc, seq2, eff2, lgth2 );

			headgapfreq1 = 0.0;
			headgapfreq2 = 0.0;
			gapfreq1pt[lgth1] = 0.0;
			gapfreq2pt[lgth2] = 0.0;
		}

		if( legacygapcost == 0 ) // sgap1 no toki headgapfreq, gapfreq1pt[lgth] ha 1 toha kagiranai.
		{
			if( !msprof ) // work/v10w
			{
			gapcountf( gapfreq1pt, seq1, icyc, eff1, lgth1 );
			gapcountf( gapfreq2pt, seq2, jcyc, eff2, lgth2 );
			}
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
#endif
	}



#if ATO
	for( i=0; i<lgth1; i++ ) 
	{
		ogcp1[i] = 0.5 * ( 1.0 - ogcp1opt[i] ) * fpenalty * ( gapfreq1pt[i] );
		fgcp1[i] = 0.5 * ( 1.0 - fgcp1opt[i] ) * fpenalty * ( gapfreq1pt[i] );
//		fprintf( stderr, "fgcp1[%d] = %f\n", i, fgcp1[i] );
	}
	for( i=0; i<lgth2; i++ ) 
	{
		ogcp2[i] = 0.5 * ( 1.0 - ogcp2opt[i] ) * fpenalty * ( gapfreq2pt[i] );
		fgcp2[i] = 0.5 * ( 1.0 - fgcp2opt[i] ) * fpenalty * ( gapfreq2pt[i] );
//		fprintf( stderr, "fgcp2[%d] = %f\n", i, fgcp2[i] );
	}
	gapinfo[0] = ogcp1;
	gapinfo[1] = fgcp1;
	gapinfo[2] = ogcp2;
	gapinfo[3] = fgcp2;
	gapinfo[4] = gapfreq1pt;
	gapinfo[5] = gapfreq2pt;
#else
	if( sgap1 )
	{
		new_OpeningGapCount( ogcp1, icyc, seq1, eff1, lgth1, sgap1 );
		new_OpeningGapCount( ogcp2, jcyc, seq2, eff2, lgth2, sgap2 );
		new_FinalGapCount( fgcp1, icyc, seq1, eff1, lgth1, egap2 );
		new_FinalGapCount( fgcp2, jcyc, seq2, eff2, lgth2, egap2 );
		outgapcount( &headgapfreq1, icyc, sgap1, eff1 );
		outgapcount( &headgapfreq2, jcyc, sgap2, eff2 );
		outgapcount( gapfreq1f+lgth1, icyc, egap1, eff1 );
		outgapcount( gapfreq2f+lgth2, jcyc, egap2, eff2 );
	}
	else
	{
		st_OpeningGapCount( ogcp1, icyc, seq1, eff1, lgth1 );
		st_OpeningGapCount( ogcp2, jcyc, seq2, eff2, lgth2 );
		st_FinalGapCount( fgcp1, icyc, seq1, eff1, lgth1 );
		st_FinalGapCount( fgcp2, jcyc, seq2, eff2, lgth2 );
		headgapfreq1 = 0.0;
		headgapfreq2 = 0.0;
		gapfreq1f[lgth1] = 0.0;
		gapfreq2f[lgth2] = 0.0;
	}

	if( legacygapcost == 0 )
	{
		gapcountf( gapfreq1f, seq1, icyc, eff1, lgth1 );
		gapcountf( gapfreq2f, seq2, jcyc, eff2, lgth2 );
		for( i=0; i<lgth1+1; i++ ) gapfreq1f[i] = 1.0 - gapfreq1f[i];
		for( i=0; i<lgth2+1; i++ ) gapfreq2f[i] = 1.0 - gapfreq2f[i];
		headgapfreq1 = 1.0 - headgapfreq1;
		headgapfreq2 = 1.0 - headgapfreq2;
	}
	else
	{
		for( i=0; i<lgth1+1; i++ ) gapfreq1f[i] = 1.0;
		for( i=0; i<lgth2+1; i++ ) gapfreq2f[i] = 1.0;
		headgapfreq1 = 1.0;
		headgapfreq2 = 1.0;
	}
	for( i=0; i<lgth1; i++ ) 
	{
		ogcp1[i] = 0.5 * ( 1.0 - ogcp1[i] ) * fpenalty * ( gapfreq1f[i] );
		fgcp1[i] = 0.5 * ( 1.0 - fgcp1[i] ) * fpenalty * ( gapfreq1f[i] );
//		fprintf( stderr, "fgcp1[%d] = %f\n", i, fgcp1[i] );
	}
	for( i=0; i<lgth2; i++ ) 
	{
		ogcp2[i] = 0.5 * ( 1.0 - ogcp2[i] ) * fpenalty * ( gapfreq2f[i] );
		fgcp2[i] = 0.5 * ( 1.0 - fgcp2[i] ) * fpenalty * ( gapfreq2f[i] );
//		fprintf( stderr, "fgcp2[%d] = %f\n", i, fgcp2[i] );
	}

	gapinfo[0] = ogcp1;
	gapinfo[1] = fgcp1;
	gapinfo[2] = ogcp2;
	gapinfo[3] = fgcp2;
	gapinfo[4] = gapfreq1f;
	gapinfo[5] = gapfreq2f;
#endif

#if 0
	fprintf( stdout, "in MSalignmm.c\n" );
	for( i=0; i<icyc; i++ )
	{
		fprintf( stdout, ">%d of GROUP1\n", i );
		fprintf( stdout, "%s\n", seq1[i] );
	}
	for( i=0; i<jcyc; i++ )
	{
		fprintf( stdout, ">%d of GROUP2\n", i );
		fprintf( stdout, "%s\n", seq2[i] );
	}
	fflush( stdout );
#endif

	if( msfast ) // work/v10w
	{
		msf_t msc;
		if( !msf_init( &msc, n_dynamicmtx, cpmx1pt, cpmx2pt, gapinfo, lgth1, lgth2, headgp, tailgp, headgapfreq1, headgapfreq2 ) ) ErrorExit( "msf_init\n" );
		wm = msf_solve( &msc, lgth1, lgth2, mgt1, mgt2 );
		msf_free( &msc );
		if( mscheck() ) // the original recursion into its own buffers, compared bitwise (the rows below)
		{
			double wref;
			ckseq1 = AllocateCharMtx( icyc, ll1+ll2 );
			ckseq2 = AllocateCharMtx( jcyc, ll1+ll2 );
			ckgt1 = AllocateCharVec( ll1+ll2 );
			ckgt2 = AllocateCharVec( ll1+ll2 );
			wref = MSalignmm_rec( n_dynamicmtx, icyc, jcyc, eff1, eff2, seq1, seq2, cpmx1pt, cpmx2pt, 0, lgth1-1, 0, lgth2-1, alloclen, lgth1, lgth2, ckseq1, ckseq2, ckgt1, ckgt2, 0, gapinfo, NULL, 0, NULL, headgp, tailgp, headgapfreq1, headgapfreq2 );
			if( memcmp( &wref, &wm, sizeof( double ) ) ) { reporterr( "\nMSCHECK MISMATCH: score %.17g vs %.17g (%d x %d)\n", wm, wref, lgth1, lgth2 ); exit( 1 ); }
			if( strcmp( ckgt1, mgt1 ) || strcmp( ckgt2, mgt2 ) ) { reporterr( "\nMSCHECK MISMATCH: gap tables (%d x %d)\n", lgth1, lgth2 ); exit( 1 ); }
		}
	}
	else
	wm = MSalignmm_rec( n_dynamicmtx, icyc, jcyc, eff1, eff2, seq1, seq2, cpmx1pt, cpmx2pt, 0, lgth1-1, 0, lgth2-1, alloclen, lgth1, lgth2, mseq1, mseq2, mgt1, mgt2, 0, gapinfo, chudanpt, chudanref, chudanres, headgp, tailgp, headgapfreq1, headgapfreq2 );
#ifdef enablemultithread
	if( chudanres && *chudanres ) 
	{
//		fprintf( stderr, "\n\n## CHUUDAN!!! relay\n" );
		*chudanres = 1;
		freearrays( ogcp1, ogcp2, ogcp1o, ogcp2o, fgcp1, fgcp2, fgcp1o, fgcp2o, cpmx1, cpmx2, gapfreq1f, gapfreq2f, gapinfo, mseq1, mseq2, mgt1, mgt2 );
		return( -1.0 );
	}
#endif


#if 1
	if( cpmxresult )
	{
		if( icyc + jcyc > 20 )
//		if( 0 )
//		if( 1 )
		{
#if 1 // marume gosa wo teigen suru tame
			double totaleff1 = 0.0;
			double totaleff2 = 0.0;
			for( i=0; i<icyc; i++ ) totaleff1 += eff1[i];
			for( j=0; j<jcyc; j++ ) totaleff2 += eff2[j];

			if( fabs(totaleff1-1.0)>0.001 || fabs(totaleff2-1.0)>0.001 )
			{
				reporterr( "Warning: rounding error may be large.  totaleff1 = %50.40f\n", totaleff1 );
				reporterr( "Warning: rounding error may be large.  totaleff2 = %50.40f\n", totaleff2 );
				exit( 1 );
			}
			totaleff1 = totaleff1 * orieff1 / (orieff1 + orieff2);
			totaleff2 = totaleff2 * orieff2 / (orieff1 + orieff2);
#else
			double totaleff1 = orieff1 / ( orieff1 + orieff2 );
			double totaleff2 = orieff2 / ( orieff1 + orieff2 );
#endif

	
			*cpmxresult = AllocateDoubleMtx( nalphabets+3, strlen( mgt1 )+1 ); // gapcount, opg, fng no bun de +3
			createcpmxresult( *cpmxresult, totaleff1, totaleff2, cpmx1pt, cpmx2pt, mgt1, mgt2 );
#if ATO
			creategapfreqresult( (*cpmxresult)[nalphabets], totaleff1, totaleff2, gapfreq1pt, gapfreq2pt, mgt1, mgt2 );
			createogresult( (*cpmxresult)[nalphabets+1], totaleff1, totaleff2, ogcp1o, ogcp2o, gapfreq1pt, gapfreq2pt, mgt1, mgt2 );
			createfgresult( (*cpmxresult)[nalphabets+2], totaleff1, totaleff2, fgcp1o, fgcp2o, gapfreq1pt, gapfreq2pt, mgt1, mgt2 );
#endif

#if 0
			reporterr( "\n" );
			for( j=0; j<nalphabets; j++ ) 
				reporterr( "%4c ", amino[j] );
			reporterr( "\ncpmxresult = \n" );
			for( i=0; i<strlen( gaptable1 ); i++ ) 
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

// matomete free
	if( cpmx1pt != cpmx1 && cpmxchild0 && *cpmxchild0 ) 
	{
//		reporterr( "freeing cpmxchild0\n" );
		FreeDoubleMtx( *cpmxchild0 );
		*cpmxchild0 = NULL;
	}
	if( cpmx2pt != cpmx2 && cpmxchild1 && *cpmxchild1 ) 
	{
//		reporterr( "freeing cpmxchild1\n" );
		FreeDoubleMtx( *cpmxchild1 );
		*cpmxchild1 = NULL;
	}


#if 0
		fprintf( stderr, "\n" );
		fprintf( stderr, " seq1[0] = %s\n", seq1[0] );
		fprintf( stderr, " seq2[0] = %s\n", seq2[0] );
		fprintf( stderr, "mseq1[0] = %s\n", mseq1[0] );
		fprintf( stderr, "mseq2[0] = %s\n", mseq2[0] );
		fprintf( stderr, "\n" );
#endif

//	fprintf( stderr, "wm = %f\n", wm );


	if( msfast ) // work/v10w: the rows from the gap tables (= gapireru of each input row, as the recursion writes them)
	{
		msf_rows( seq1, icyc, lgth1, mgt1, *newgapstr );
		msf_rows( seq2, jcyc, lgth2, mgt2, *newgapstr );
		if( ckseq1 ) // MAFFT_MSCHECK
		{
			for( i=0; i<icyc; i++ ) if( strcmp( seq1[i], ckseq1[i] ) ) { reporterr( "\nMSCHECK MISMATCH: group 1 row %d\n", i ); exit( 1 ); }
			for( i=0; i<jcyc; i++ ) if( strcmp( seq2[i], ckseq2[i] ) ) { reporterr( "\nMSCHECK MISMATCH: group 2 row %d\n", i ); exit( 1 ); }
			FreeCharMtx( ckseq1 ); FreeCharMtx( ckseq2 ); free( ckgt1 ); free( ckgt2 );
			__atomic_add_fetch( mscheck_n + 0, 1, __ATOMIC_RELAXED );
			__atomic_add_fetch( mscheck_n + 1, icyc + jcyc, __ATOMIC_RELAXED );
		}
	}
	else
	{
	for( i=0; i<icyc; i++ ) strcpy( seq1[i], mseq1[i] );
	for( i=0; i<jcyc; i++ ) strcpy( seq2[i], mseq2[i] );
	}

	if( seqlen( seq1[0] ) != nglen1 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen1) seqlen(seq1[0])=%d but nglen1=%d\n", seqlen( seq1[0] ), nglen1 );
		fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
		exit( 1 );
	}
	if( seqlen( seq2[0] ) != nglen2 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen2) seqlen(seq2[0])=%d but nglen2=%d\n", seqlen( seq2[0] ), nglen2 );
		exit( 1 );
	}

	freearrays( ogcp1, ogcp2, ogcp1o, ogcp2o, fgcp1, fgcp2, fgcp1o, fgcp2o, cpmx1, cpmx2, gapfreq1f, gapfreq2f, gapinfo, mseq1, mseq2, mgt1, mgt2 );

	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );
	for( i=0; i<icyc; i++ ) 
	{
		if( strlen( seq1[i] ) != lgth1 )
		{
			fprintf( stderr, "i = %d / %d\n", i, icyc );
			fprintf( stderr, "hairetsu ga kowareta (end of MSalignmm) !\n" );
			exit( 1 );
		}
	}
	for( j=0; j<jcyc; j++ )
	{
		if( strlen( seq2[j] ) != lgth2 )
		{
			fprintf( stderr, "j = %d / %d\n", j, jcyc );
			fprintf( stderr, "hairetsu ga kowareta (end of MSalignmm) !\n" );
			exit( 1 );
		}
	}

#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, " seq1[i] = %s\n", seq1[i] );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, " seq2[j] = %s\n", seq2[j] );
	fprintf( stderr, "\n" );
#endif

	return( wm );
}

// -------------------------------
//    variousdist
// -------------------------------

static void fillzero( double *s, int l )
{
	while( l-- ) *s++ = 0.0; 
}

static double MSalignmm_tanni_variousdist( double ***matrices, int icyc, int jcyc, char **seq1, char **seq2, double ***cpmx1s, double ***cpmx2s, int ist, int ien, int jst, int jen, int alloclen, int fulllen1, int fulllen2, char **mseq1, char **mseq2, double **gapinfo, int headgp, int tailgp, double headgapfreq1_g, double headgapfreq2_g )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
//	int k;
	register int i, j, c;
	int ll1, ll2;
	int lasti, lastj;
	double wm = 0.0;   /* int ?????? */
	double g;
	double *currentw, *previousw;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penalty_ex;
#endif
#if 1
	double *wtmp;
	int *ijppt;
	double *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	double mi, *m;
	int **ijp;
	int mpi, *mp;
	double *w1, *w2;
	double *initverticalw;    /* kufuu sureba iranai */
	double *lastverticalw;    /* kufuu sureba iranai */
	int ***intwork;
	double ***doublework;
	int **intmtx;
	int lgth1, lgth2;
	double *ogcp1;
	double *fgcp1;
	double *ogcp2;
	double *fgcp2;
	double *gapfreq1f;
	double *gapfreq2f;
//	char **aseq1;
//	char **aseq2;
//	char **aseq1bk, **aseq2bk;
	double headgapfreq1;
	double headgapfreq2;


	ogcp1 = gapinfo[0] + ist;
	fgcp1 = gapinfo[1] + ist;
	ogcp2 = gapinfo[2] + jst;
	fgcp2 = gapinfo[3] + jst;
	gapfreq1f = gapinfo[4] + ist;
	gapfreq2f = gapinfo[5] + jst;

	if( ist > 0 ) headgapfreq1 = gapfreq1f[-1];
	else headgapfreq1 = headgapfreq1_g;
	if( jst > 0 ) headgapfreq2 = gapfreq2f[-1];
	else headgapfreq2 = headgapfreq2_g;

#if STOREWM
	char ttt1[10000], ttt2[10000];
#endif


	lgth1 = ien-ist+1;
	lgth2 = jen-jst+1;

#if STOREWM
	strncpy( ttt1, seq1[0]+ist, lgth1 ); ttt1[lgth1] = 0;
	strncpy( ttt2, seq2[0]+jst, lgth2 ); ttt2[lgth2] = 0;

	fprintf( stderr, "in _tanni ist,ien = %d,%d, lgth1=%d\n", ist, ien, lgth1 );
	fprintf( stderr, "in _tanni jst,jen = %d,%d, lgth2=%d\n", jst, jen, lgth2 );
	fprintf( stderr, "ttt1 = %s\n", ttt1 );
	fprintf( stderr, "ttt2 = %s\n", ttt2 );
#endif

#if 0
	fprintf( stderr, "in _tanni ist,ien = %d,%d, fulllen1=%d\n", ist, ien, fulllen1 );
	fprintf( stderr, "in _tanni jst,jen = %d,%d, fulllen2=%d\n", jst, jen, fulllen2 );
	fprintf( stderr, "in _tanni seq1[0] = %-*.*s\n", ien-ist+1, ien-ist+1, seq1[0]+ist );
	fprintf( stderr, "in _tanni seq2[0] = %-*.*s\n", jen-jst+1, jen-jst+1, seq2[0]+jst );
#endif


	ll1 = ( (int)(lgth1) ) + 100;
	ll2 = ( (int)(lgth2) ) + 100;

//	aseq1 = AllocateCharMtx( icyc, 0 );
//	aseq2 = AllocateCharMtx( jcyc, 0 );
//	aseq1bk = AllocateCharMtx( icyc, lgth1+lgth2+100 );
//	aseq2bk = AllocateCharMtx( jcyc, lgth1+lgth2+100 );
//	for( i=0; i<icyc; i++ ) aseq1[i] = aseq1bk[i];
//	for( i=0; i<jcyc; i++ ) aseq2[i] = aseq2bk[i];

	w1 = AllocateFloatVec( ll2+2 );
	w2 = AllocateFloatVec( ll2+2 );

	initverticalw = AllocateFloatVec( ll1+2 );
	lastverticalw = AllocateFloatVec( ll1+2 );

	m = AllocateFloatVec( ll2+2 );
	mp = AllocateIntVec( ll2+2 );

	doublework = AllocateFloatCub( maxdistclass, MAX( ll1, ll2 )+2, nalphabets+1 ); 
	intwork = AllocateIntCub( maxdistclass, MAX( ll1, ll2 )+2, nalphabets+1 ); 


	intmtx = AllocateIntMtx( ll1+1, ll2+1 );

	ijp = intmtx;

	currentw = w1;
	previousw = w2;

#if 0
	match_calc( n_dynamicmtx, initverticalw, cpmx2+jst, cpmx1+ist, 0, lgth1, doublework, intwork, 1 );
	match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, 0, lgth2, doublework, intwork, 1 );
#else
	fillzero( initverticalw, lgth1 );
	for( c=0; c<maxdistclass; c++ )
		match_calc_add( matrices[c], initverticalw, cpmx2s[c]+jst, cpmx1s[c]+ist, 0, lgth1, doublework[c], intwork[c], 1 );

	fillzero( currentw, lgth2 );
	for( c=0; c<maxdistclass; c++ )
		match_calc_add( matrices[c], currentw, cpmx1s[c]+ist, cpmx2s[c]+jst, 0, lgth2, doublework[c], intwork[c], 1 );
#endif

	if( headgp || ist != 0 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2f[0] ) ;
		}
	}
	if( headgp || jst != 0 )
	{
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1f[0] ) ;
		}
	}

	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2f[j-1]; mp[j] = 0;;
	}

	lastverticalw[0] = currentw[lgth2-1];



	if( tailgp || jen != fulllen2-1 ) lasti = lgth1+1; else lasti = lgth1;
//	if( 1 ) lasti = lgth1+1; else lasti = lgth1;
	for( i=1; i<lasti; i++ )
	{
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

#if 0
		match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, i, lgth2, doublework, intwork, 0 );
#else
		fillzero( currentw, lgth2 );
		for( c=0; c<maxdistclass; c++ )
			match_calc_add( matrices[c], currentw, cpmx1s[c]+ist, cpmx2s[c]+jst, i, lgth2, doublework[c], intwork[c], 0 );
#endif
		currentw[0] = initverticalw[i];

		mi = previousw[0] + ogcp2[1] * gapfreq1f[i-1];
		mpi = 0;

		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;
//		if( tailgp && jen != fulllen2-1 ) lastj = lgth2+1; else lastj = lgth2;
		lastj = lgth2+1; 
		for( j=1; j<lastj; j++ )
		{
			wm = *prept;
			*ijppt = 0;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
			g = mi + fgcp2[j-1]  * gapfreq1f[i];
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
				*ijppt = -( j - mpi );
			}
			g = *prept + ogcp2[j] * gapfreq1f[i-1];
			if( g >= mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

			g = *mjpt + fgcp1[i-1] * gapfreq2f[j];
#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
				*ijppt = +( i - *mpjpt );
			}
			g = *prept + ogcp1[i] * gapfreq2f[j-1];
			if( g >= *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}
#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif

#if 0
			fprintf( stderr, "%5.0f ", wm );
#endif
			*curpt += wm;


			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
			curpt++;
		}
		lastverticalw[i] = currentw[lgth2-1];
	}

//	fprintf( stderr, "wm = %f\n", wm );

	Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, ist, ien, jst, jen, fulllen1, fulllen2, tailgp, NULL, NULL );
#if 0
	fprintf( stderr, "res in _tanni mseq1[0] = %s\n", mseq1[0] );
	fprintf( stderr, "res in _tanni mseq2[0] = %s\n", mseq2[0] );
#endif

//	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], aseq1[i] );
//	for( i=0; i<jcyc; i++ ) strcpy( mseq2[i], aseq2[i] );

//	fprintf( stderr, "in _tanni, aseq1 = %s\n", aseq1[0] );
//	fprintf( stderr, "in _tanni, mseq1 = %s\n", mseq1[0] );

	FreeFloatVec( w1 );
	FreeFloatVec( w2 );
	FreeFloatVec( initverticalw );
	FreeFloatVec( lastverticalw );

	FreeFloatVec( m );
	FreeIntVec( mp );


	FreeFloatCub( doublework );
	FreeIntCub( intwork );

	FreeIntMtx( intmtx );


//	FreeCharMtx( aseq1bk );
//	FreeCharMtx( aseq2bk );

//	free( aseq1 );
//	free( aseq2 );

	return( wm );

}

static void freearrays_rec1_variousdist(
	double *w1, double *w2, double *initverticalw, double *lastverticalw,
	double *midw, double *midm, double *midn,
	int *jumpbacki, int *jumpbackj, int *jumpforwi, int *jumpforwj, int *jumpdummi, int *jumpdummj,
	double *m, int *mp,
	double ***doublework, int ***intwork
#if STOREWM
	, double **WMMTX, double **WMMTX2
#endif
)
{

	FreeFloatVec( w1 );
	FreeFloatVec( w2 );
	FreeFloatVec( initverticalw );
	FreeFloatVec( lastverticalw );
	FreeFloatVec( midw );
	FreeFloatVec( midm );
	FreeFloatVec( midn );

	FreeIntVec( jumpbacki );
	FreeIntVec( jumpbackj );
	FreeIntVec( jumpforwi );
	FreeIntVec( jumpforwj );
	FreeIntVec( jumpdummi );
	FreeIntVec( jumpdummj );

	FreeFloatVec( m );
	FreeIntVec( mp );

	FreeFloatCub( doublework );
	FreeIntCub( intwork );

#if STOREWM
	FreeFloatMtx( WMMTX );
	FreeFloatMtx( WMMTX2 );
#endif
}

static double MSalignmm_rec_variousdist( double ***matrices, int icyc, int jcyc, char **seq1, char **seq2, double ***cpmx1s, double ***cpmx2s, int ist, int ien, int jst, int jen, int alloclen, int fulllen1, int fulllen2, char **mseq1, char **mseq2, int depth, double **gapinfo, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp, double headgapfreq1_g, double headgapfreq2_g )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
//	int k;
	int alnlen;
	double value = 0.0;
	register int i, j, c;
	char **aseq1, **aseq2;
	int ll1, ll2, l, len;
	int lasti, lastj, imid;
	int jmid = 0; // by D.Mathog, a guess
	double wm = 0.0;   /* int ?????? */
	double g;
	double *currentw, *previousw;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penalty_ex;
#endif
	double *wtmp;
//	short *ijppt;
	int *mpjpt;
//	short **ijp;
	int *mp;
	int mpi;
	double *mjpt, *prept, *curpt;
	double mi;
	double *m;
	double *w1, *w2;
//	double *match;
	double *initverticalw;    /* kufuu sureba iranai */
	double *lastverticalw;    /* kufuu sureba iranai */
	int ***intwork;
	double ***doublework;
//	short **shortmtx;
#if STOREWM
	double **WMMTX;
	double **WMMTX2;
#endif
	double *midw;
	double *midm;
	double *midn;
	int lgth1, lgth2;
	double maxwm;
	int *jumpforwi;
	int *jumpforwj;
	int *jumpbacki;
	int *jumpbackj;
	int *jumpdummi; //muda
	int *jumpdummj = NULL; // by D.Mathog, a guess
	int jumpi, jumpj = 0; // by D.Mathog, a guess
	char *gaps;
	int ijpi, ijpj;
	double *ogcp1;
	double *fgcp1;
	double *ogcp2;
	double *fgcp2;
	double firstm;
	int firstmp;
#if STOREWM
	static TLS char ttt1[50000];
	static TLS char ttt2[50000];
#endif
	double *gapfreq1f;
	double *gapfreq2f;
	double headgapfreq1;
	double headgapfreq2;

#if 0
	int nglen1, nglen2;
	nglen1 = seqlen( seq1[0] );
	nglen2 = seqlen( seq2[0] );
#endif

//	fprintf( stderr, "fulllen1 = %d, fulllen2 = %d, headgp = %d, tailgp = %d\n", fulllen1, fulllen2, headgp, tailgp );

	ogcp1 = gapinfo[0] + ist;
	fgcp1 = gapinfo[1] + ist;
	ogcp2 = gapinfo[2] + jst;
	fgcp2 = gapinfo[3] + jst;
	gapfreq1f = gapinfo[4] + ist;
	gapfreq2f = gapinfo[5] + jst;

	if( ist > 0 ) headgapfreq1 = gapfreq1f[-1];
	else headgapfreq1 = headgapfreq1_g;
	if( jst > 0 ) headgapfreq2 = gapfreq2f[-1];
	else headgapfreq2 = headgapfreq2_g;

	depth++;
	reccycle++;

	lgth1 = ien-ist+1;
	lgth2 = jen-jst+1;

//	if( lgth1 < 5 )
//		fprintf( stderr, "\nWARNING: lgth1 = %d\n", lgth1 );
//	if( lgth2 < 5 )
//		fprintf( stderr, "\nWARNING: lgth2 = %d\n", lgth2 );
//


#if STOREWM
	fprintf( stderr, "==== MSalign (depth=%d, reccycle=%d), ist=%d, ien=%d, jst=%d, jen=%d\n", depth, reccycle, ist, ien, jst, jen );
	strncpy( ttt1, seq1[0]+ist, lgth1 );
	strncpy( ttt2, seq2[0]+jst, lgth2 );
	ttt1[lgth1] = 0;
	ttt2[lgth2] = 0;
	fprintf( stderr, "seq1 = %s\n", ttt1 );
	fprintf( stderr, "seq2 = %s\n", ttt2 );
#endif
	if( lgth2 <= 0 ) // lgth1 <= 0 ha?
	{
//		fprintf( stderr, "\n\n==== jimei\n\n" );
//		exit( 1 );
		for( i=0; i<icyc; i++ ) 
		{
			strncpy( mseq1[i], seq1[i]+ist, lgth1 );
			mseq1[i][lgth1] = 0;
		}
		for( i=0; i<jcyc; i++ ) 
		{
			mseq2[i][0] = 0;
			for( j=0; j<lgth1; j++ )
//				strcat( mseq2[i], "-" );
				strcat( mseq2[i], newgapstr );
		}

//		fprintf( stderr, "==== mseq1[0] (%d) = %s\n", depth, mseq1[0] );
//		fprintf( stderr, "==== mseq2[0] (%d) = %s\n", depth, mseq2[0] );

		return( 0.0 );
	}

#if MEMSAVE
	aseq1 = AllocateCharMtx( icyc, 0 );
	aseq2 = AllocateCharMtx( jcyc, 0 );
	for( i=0; i<icyc; i++ ) aseq1[i] = mseq1[i];
	for( i=0; i<jcyc; i++ ) aseq2[i] = mseq2[i];
#else
	aseq1 = AllocateCharMtx( icyc, lgth1+lgth2+100 );
	aseq2 = AllocateCharMtx( jcyc, lgth1+lgth2+100 );
#endif

//	fprintf( stderr, "####(s) seq1[0] (%d) = \n%-*.*s\n (a%d-a%d)\n", depth, ien-ist+1, ien-ist+1, seq1[0]+ist, ist, ien );
//	fprintf( stderr, "####(s) seq2[0] (%d) = \n%-*.*s\n (b%d-b%d)\n", depth, jen-jst+1, jen-jst+1, seq2[0]+jst, jst, jen );

//  if( lgth1 < DPTANNI && lgth2 < DPTANNI ) // & dato lgth ==1 no kanousei ga arunode yokunai 
//    if( lgth1 < DPTANNI ) // kore mo lgth2 ga mijikasugiru kanousei ari
    if( lgth1 < DPTANNI || lgth2 < DPTANNI ) // zettai ni anzen ka?
	{
//		fprintf( stderr, "==== Going to _tanni\n" );

		value = MSalignmm_tanni_variousdist( matrices, icyc, jcyc, seq1, seq2, cpmx1s, cpmx2s, ist, ien, jst, jen, alloclen, fulllen1, fulllen2, aseq1, aseq2, gapinfo, headgp, tailgp, headgapfreq1_g, headgapfreq2_g );	


#if MEMSAVE
		free( aseq1 );
		free( aseq2 );
#else
		for( i=0; i<icyc; i++ ) strcpy( mseq1[i], aseq1[i] );
		for( i=0; i<jcyc; i++ ) strcpy( mseq2[i], aseq2[i] );

		FreeCharMtx( aseq1 );
		FreeCharMtx( aseq2 );
#endif

//		fprintf( stderr, "value = %f\n", value );

		return( value );
	}
//	fprintf( stderr, "Trying to divide the mtx\n" );

	ll1 = ( (int)(lgth1) ) + 100;
	ll2 = ( (int)(lgth2) ) + 100;

//	fprintf( stderr, "ll1,ll2=%d,%d\n", ll1, ll2 );

	w1 = AllocateFloatVec( ll2+2 );
	w2 = AllocateFloatVec( ll2+2 );
//	match = AllocateFloatVec( ll2+2 );
	midw = AllocateFloatVec( ll2+2 );
	midn = AllocateFloatVec( ll2+2 );
	midm = AllocateFloatVec( ll2+2 );
	jumpbacki = AllocateIntVec( ll2+2 );
	jumpbackj = AllocateIntVec( ll2+2 );
	jumpforwi = AllocateIntVec( ll2+2 );
	jumpforwj = AllocateIntVec( ll2+2 );
	jumpdummi = AllocateIntVec( ll2+2 );
	jumpdummj = AllocateIntVec( ll2+2 );

	initverticalw = AllocateFloatVec( ll1+2 );
	lastverticalw = AllocateFloatVec( ll1+2 );

	m = AllocateFloatVec( ll2+2 );
	mp = AllocateIntVec( ll2+2 );
	gaps = AllocateCharVec( MAX( ll1, ll2 ) + 2 );

	doublework = AllocateFloatCub( maxdistclass, MAX( ll1, ll2 )+2, nalphabets ); 
	intwork = AllocateIntCub( maxdistclass, MAX( ll1, ll2 )+2, nalphabets ); 

#if DEBUG
	fprintf( stderr, "succeeded\n" );
#endif

#if STOREWM
	WMMTX = AllocateFloatMtx( ll1, ll2 );
	WMMTX2 = AllocateFloatMtx( ll1, ll2 );
#endif
#if 0
	shortmtx = AllocateShortMtx( ll1, ll2 );

#if DEBUG
	fprintf( stderr, "succeeded\n\n" );
#endif

	ijp = shortmtx;
#endif

	currentw = w1;
	previousw = w2;

#if 0
	match_calc( n_dynamicmtx, initverticalw, cpmx2+jst, cpmx1+ist, 0, lgth1, doublework, intwork, 1 );
	match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, 0, lgth2, doublework, intwork, 1 );
#else
	fillzero( initverticalw, lgth1 );
	for( c=0; c<maxdistclass; c++ )
		match_calc_add( matrices[c], initverticalw, cpmx2s[c]+jst, cpmx1s[c]+ist, 0, lgth1, doublework[c], intwork[c], 1 );

	fillzero( currentw, lgth2 );
	for( c=0; c<maxdistclass; c++ )
		match_calc_add( matrices[c], currentw, cpmx1s[c]+ist, cpmx2s[c]+jst, 0, lgth2, doublework[c], intwork[c], 1 );
#endif


	for( i=1; i<lgth1+1; i++ )
	{
//		initverticalw[i] += ( ogcp1[0] + fgcp1[i-1] ) ;
		initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2f[0] ) ;
	}
	for( j=1; j<lgth2+1; j++ )
	{
		currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1f[0]) ;
	}

#if STOREWM
	WMMTX[0][0] = initverticalw[0];
	for( i=1; i<lgth1+1; i++ )
	{
		WMMTX[i][0] = initverticalw[i];
	}
	for( j=1; j<lgth2+1; j++ )
	{
		WMMTX[0][j] = currentw[j];
	}
#endif


	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2f[j-1];
//		m[j] = currentw[j-1];
		mp[j] = 0;
	}

	lastverticalw[0] = currentw[lgth2-1];

	imid = lgth1 * 0.5;

	jumpi = 0; // atode kawaru.
	lasti = lgth1+1;
#if STOREWM
	for( i=1; i<lasti; i++ )
#else
	for( i=1; i<=imid; i++ )
#endif
	{
#ifdef enablemultithread
//		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
		if( chudanpt && *chudanpt != chudanref ) 
		{
//			fprintf( stderr, "\n\n## CHUUDAN!!! zenhan\n" );
			*chudanres = 1;
			freearrays_rec1_variousdist
			(
				w1, w2, initverticalw, lastverticalw, midw, midm, midn,
				jumpbacki, jumpbackj, jumpforwi, jumpforwj, jumpdummi, jumpdummj,
				m, mp,
				doublework, intwork
#if STOREWM
				, WMMTX, WMMTX2
#endif
			);
			freearrays_rec2( gaps, aseq1, aseq2 );
			return( -1.0 );
		}
#endif
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

#if 0
		match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, i, lgth2, doublework, intwork, 0 );
#else
		fillzero( currentw, lgth2 );
		for( c=0; c<maxdistclass; c++ )
			match_calc_add( matrices[c], currentw, cpmx1s[c]+ist, cpmx2s[c]+jst, i, lgth2, doublework[c], intwork[c], 0 );
#endif
		currentw[0] = initverticalw[i];

		m[0] = ogcp1[i];
#if STOREM
		WMMTX2[i][0] = m[0];
#endif
		if( i == imid ) midm[0] = m[0];

		mi = previousw[0] + ogcp2[1] * gapfreq1f[i-1]; 
//		mi = previousw[0];
		mpi = 0;


//		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;


		lastj = lgth2+1;
		for( j=1; j<lastj; j++ )
		{

			wm = *prept;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
			g = mi + fgcp2[j-1] * gapfreq1f[i];
//			g = mi + fpenalty;
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
//				*ijppt = -( j - mpi );
			}
			g = *prept + ogcp2[j] * gapfreq1f[i-1];
//			g = *prept;
			if( g >= mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

			g = *mjpt + fgcp1[i-1] * gapfreq2f[j];
//			g = *mjpt + fpenalty;
#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
//				*ijppt = +( i - *mpjpt );
			}


			g = *prept + ogcp1[i] * gapfreq2f[j-1];
//			g = *prept;
			if( g >= *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}
#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif

#if 0
			fprintf( stderr, "%5.0f ", wm );
#endif
			*curpt += wm;

#if STOREWM
			WMMTX[i][j] = *curpt;
			WMMTX2[i][j] = *mjpt;
#endif

			if( i == imid ) //muda
			{	
				jumpbackj[j] = *mpjpt; // muda atode matomeru
				jumpbacki[j] = mpi; // muda atode matomeru
//				fprintf( stderr, "jumpbackj[%d] in forward dp is %d\n", j, *mpjpt );
//				fprintf( stderr, "jumpbacki[%d] in forward dp is %d\n", j, mpi );
				midw[j] = *curpt;
				midm[j] = *mjpt;
				midn[j] = mi;
			}

//			fprintf( stderr, "m[%d] = %f\n", j, m[j] );
			mjpt++;
			prept++;
			mpjpt++;
			curpt++;

		}
		lastverticalw[i] = currentw[lgth2-1];

#if STOREWM
		WMMTX2[i][lgth2] = m[lgth2-1];
#endif

#if 0  // ue
		if( i == imid )
		{
			for( j=0; j<lgth2; j++ ) midw[j] = currentw[j];
			for( j=0; j<lgth2; j++ ) midm[j] = m[j];
		}
#endif
	}
//	for( j=0; j<lgth2; j++ ) midw[j] = WMMTX[imid][j];
//	for( j=0; j<lgth2; j++ ) midm[j] = WMMTX2[imid][j];

#if 0
    for( i=0; i<lgth1; i++ )
    {
        for( j=0; j<lgth2; j++ )
        {
            fprintf( stderr, "% 10.2f ", WMMTX[i][j] );
        }
        fprintf( stderr, "\n" );
    }
	fprintf( stderr, "\n" );
	fprintf( stderr, "WMMTX2 = \n" );
    for( i=0; i<lgth1; i++ )
    {
        for( j=0; j<lgth2; j++ )
        {
            fprintf( stderr, "% 10.2f ", WMMTX2[i][j] );
        }
        fprintf( stderr, "\n" );
    }
	fprintf( stderr, "\n" );
#endif

// gyakudp

#if 0
	match_calc( n_dynamicmtx, initverticalw, cpmx2+jst, cpmx1+ist, lgth2-1, lgth1, doublework, intwork, 1 );
	match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, lgth1-1, lgth2, doublework, intwork, 1 );
#else
	fillzero( initverticalw, lgth1 );
	for( c=0; c<maxdistclass; c++ )
		match_calc_add( matrices[c], initverticalw, cpmx2s[c]+jst, cpmx1s[c]+ist, lgth2-1, lgth1, doublework[c], intwork[c], 1 );

	fillzero( currentw, lgth2 );
	for( c=0; c<maxdistclass; c++ )
		match_calc_add( matrices[c], currentw, cpmx1s[c]+ist, cpmx2s[c]+jst, lgth1-1, lgth2, doublework[c], intwork[c], 1 );
#endif

	for( i=0; i<lgth1-1; i++ )
	{
//		initverticalw[i] += ( fgcp1[lgth1-1] + ogcp1[i+1] );
		initverticalw[i] += ( fgcp1[lgth1-1] * gapfreq2f[lgth2] + ogcp1[i+1] * gapfreq2f[lgth2-1] );
	}
	for( j=0; j<lgth2-1; j++ )
	{
//		currentw[j] += ( fgcp2[lgth2-1] + ogcp2[j+1] );
		currentw[j] += ( fgcp2[lgth2-1] * gapfreq1f[lgth1] + ogcp2[j+1] * gapfreq1f[lgth1-1] );
	}

#if STOREWM
	for( i=0; i<lgth1-1; i++ )
	{
		WMMTX[i][lgth2-1] += ( fgcp1[lgth1-1] + ogcp1[i+1] );
		fprintf( stderr, "fgcp1[lgth1-1] + ogcp1[i+1] = %f\n", fgcp1[lgth1-1] + ogcp1[i+1] );
	}
	for( j=0; j<lgth2-1; j++ )
	{
		WMMTX[lgth1-1][j] += ( fgcp2[lgth2-1] + ogcp2[j+1] );
		fprintf( stderr, "fgcp2[lgth2-1] + ogcp2[j+1] = %f\n", fgcp2[lgth2-1] + ogcp2[j+1] );
	}
#endif






#if 0
	for( j=lgth2-1; j>0; --j )
	{
		m[j-1] = currentw[j] + fgcp2[lgth2-2];
//		m[j-1] = currentw[j];
		mp[j] = lgth1-1;
	}
#else
	for( j=lgth2-1; j>-1; --j )
	{
		m[j] = currentw[j+1] + fgcp1[lgth1-2] * gapfreq2f[j+1];
//		m[j-1] = currentw[j];
		mp[j] = lgth1-1;
	}
#endif

//	for( j=0; j<lgth2; j++ ) m[j] = 0.0;
	// m[lgth2-1] ha irunoka?


//	for( i=lgth1-2; i>=imid; i-- )
	firstm = -9999999.9;
//	firstmp = lgth1-1;
	firstmp = lgth1;
	for( i=lgth1-2; i>-1; i-- )
	{
#ifdef enablemultithread
//		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
		if( chudanpt && *chudanpt != chudanref ) 
		{
//			fprintf( stderr, "\n\n## CHUUDAN!!! kouhan\n" );
			*chudanres = 1;
			freearrays_rec1_variousdist
			(
				w1, w2, initverticalw, lastverticalw, midw, midm, midn,
				jumpbacki, jumpbackj, jumpforwi, jumpforwj, jumpdummi, jumpdummj,
				m, mp,
				doublework, intwork
#if STOREWM
				, WMMTX, WMMTX2
#endif
			);
			freearrays_rec2( gaps, aseq1, aseq2 );
			return( -1.0 );
		}
#endif
		wtmp = previousw;
		previousw = currentw;
		currentw = wtmp;
		previousw[lgth2-1] = initverticalw[i+1];
#if 0
		match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, i, lgth2, doublework, intwork, 0 );
#else
		fillzero( currentw, lgth2 );
		for( c=0; c<maxdistclass; c++ )
			match_calc_add( matrices[c], currentw, cpmx1s[c]+ist, cpmx2s[c]+jst, i, lgth2, doublework[c], intwork[c], 0 );
#endif

		currentw[lgth2-1] = initverticalw[i];

//		m[lgth2] = fgcp1[i];
//		WMMTX2[i][lgth2] += m[lgth2];
//		fprintf( stderr, "m[] = %f\n", m[lgth2] );

		mi = previousw[lgth2-1] + fgcp2[lgth2-2] * gapfreq1f[i+1];
//		mi = previousw[lgth2-1];
		mpi = lgth2 - 1;

		mjpt = m + lgth2 - 2;
		prept = previousw + lgth2 - 1;
		curpt = currentw + lgth2 - 2;
		mpjpt = mp + lgth2 - 2;


		for( j=lgth2-2; j>-1; j-- )
		{
			wm = *prept;
			ijpi = i+1;
			ijpj = j+1;

			g = mi + ogcp2[j+1] * gapfreq1f[i];
//			g = mi + fpenalty;
			if( g > wm )
			{
				wm = g;
				ijpj = mpi;
				ijpi = i+1;
			}

			g = *prept + fgcp2[j] * gapfreq1f[i+1];
//			g = *prept;
			if( g >= mi )
			{
//				fprintf( stderr, "i,j=%d,%d - renewed! mpi = %d\n", i, j, j+1 );
				mi = g;
				mpi = j + 1;
			}

#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

//			fprintf( stderr, "i,j=%d,%d *mpjpt = %d\n", i, j, *mpjpt );
			g = *mjpt + ogcp1[i+1] * gapfreq2f[j];
//			g = *mjpt + fpenalty;
			if( g > wm )
			{
				wm = g;
				ijpi = *mpjpt;
				ijpj = j+1;
			}

//			if( i == imid )fprintf( stderr, "i,j=%d,%d \n", i, j );
			g = *prept + fgcp1[i] * gapfreq2f[j+1];
//			g = *prept;
			if( g >= *mjpt )
			{
				*mjpt = g;
				*mpjpt = i + 1;
			}

#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif

			if( i == jumpi || i == imid - 1 )
			{
				jumpforwi[j] = ijpi; //muda
				jumpforwj[j] = ijpj; //muda
//				fprintf( stderr, "jumpfori[%d] = %d\n", j, ijpi );
//				fprintf( stderr, "jumpforj[%d] = %d\n", j, ijpj );
			}
			if( i == imid ) // muda
			{
				midw[j] += wm;
//				midm[j+1] += *mjpt + fpenalty; //??????
				midm[j+1] += *mjpt; //??????
			}
			if( i == imid - 1 )
			{
//				midn[j] += mi + fpenalty;  //????
				midn[j] += mi;  //????
			}
#if STOREWM
			WMMTX[i][j] += wm;
//			WMMTX2[i][j+1] += *mjpt + fpenalty;
			WMMTX2[i][j+1] += *mjpt;
#endif
			*curpt += wm;

			mjpt--;
			prept--;
			mpjpt--;
			curpt--;
		}
//		fprintf( stderr, "adding *mjpt (=%f) to WMMTX2[%d][%d]\n", *mjpt, i, j+1 );
		g = *prept + fgcp1[i];
		if( firstm < g ) 
		{
			firstm = g;
			firstmp = i + 1;
		}
#if STOREWM
		WMMTX2[i][j+1] += firstm;
#endif
		if( i == imid ) midm[j+1] += firstm;


		if( i == imid - 1 )	
		{
			maxwm = midw[1];
			jmid = 0;
//			if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			for( j=2; j<lgth2-1; j++ )
			{
				wm = midw[j];
				if( wm > maxwm )
				{
					jmid = j;
					maxwm = wm;
				}
//				if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			}
			for( j=0; j<lgth2+1; j++ )
			{
				wm = midm[j];
				if( wm > maxwm )
				{
					jmid = j;
					maxwm = wm;
				}
//				if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			}

//			if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );


//			fprintf( stderr, "### imid=%d, jmid=%d\n", imid, jmid );
			wm = midw[jmid];
			jumpi = imid-1;
			jumpj = jmid-1;
			if( jmid > 0 && midn[jmid-1] > wm ) //060413
			{
				jumpi = imid-1;
				jumpj = jumpbacki[jmid];
				wm = midn[jmid-1];
//				fprintf( stderr, "rejump (n)\n" );
			}
			if( midm[jmid] > wm )
			{
				jumpi = jumpbackj[jmid];
				jumpj = jmid-1;
				wm = midm[jmid];
//				fprintf( stderr, "rejump (m) jumpi=%d\n", jumpi );
			}


//			fprintf( stderr, "--> imid=%d, jmid=%d\n", imid, jmid );
//			fprintf( stderr, "--> jumpi=%d, jumpj=%d\n", jumpi, jumpj );
#if STOREWM
			fprintf( stderr, "imid = %d\n", imid );
			fprintf( stderr, "midn = \n" );
			for( j=0; j<lgth2; j++ )
			{
				fprintf( stderr, "% 7.1f ", midn[j] );
			}
			fprintf( stderr, "\n" );
			fprintf( stderr, "midw = \n" );
			for( j=0; j<lgth2; j++ )
			{
				fprintf( stderr, "% 7.1f ", midw[j] );
			}
			fprintf( stderr, "\n" );
			fprintf( stderr, "midm = \n" );
			for( j=0; j<lgth2; j++ )
			{
				fprintf( stderr, "% 7.1f ", midm[j] );
			}
			fprintf( stderr, "\n" );
#endif
//			fprintf( stderr, "maxwm = %f\n", maxwm );
		}
		if( i == jumpi ) //saki?
		{
//			fprintf( stderr, "#### FIRST i=%d, jumpi<imid=%d<%d, ist=%d, ien=%d, firstmp=%d, lgth1=%d\n", i, jumpi, imid, ist, ien, firstmp, lgth1 );
//			fprintf( stderr, "#### mark 1\n" );
//			fprintf( stderr, "imid, jumpi = %d,%d\n", imid, jumpi );
//			fprintf( stderr, "jmid, jumpj = %d,%d\n", jmid, jumpj );

			if( jmid == 0 )
			{
//				fprintf( stderr, "#### CHUI2!\n" );
				jumpj = 0; jmid = 1;
#if 0 // v6.823 made
				jumpi = firstmp-1;
				imid = firstmp;
#endif
#if 0
				jumpi = 0;
				imid = 1;
#else
//				if( 1 || firstmp > 100 ) // naze 100
				if( imid < firstmp-1 ) // naze 100
				{
					jumpi = firstmp;
					imid = firstmp+1;
				}
#if 0
				else
				{
					jumpi = 0;
					imid = 1;
				}
#endif
#endif
			}

#if 0
			else if( jmid == lgth2 ) 
			{
				fprintf( stderr, "CHUI1!\n" );
				jumpi=0; jumpj=0;
				imid=jumpforwi[0]; jmid=lgth2-1;
			}
#else // 060414
			else if( jmid >= lgth2 ) 
			{
//				fprintf( stderr, "CHUI1!\n" );
				jumpi=imid-1; jmid=lgth2;
				jumpj = lgth2-1;
			}
#endif
			else
			{
//				fprintf( stderr, "#### CHUI3!\n" );
				imid = jumpforwi[jumpj];
				jmid = jumpforwj[jumpj];
				if( imid == jumpi ) jumpi = imid-1;
			}
#if 0
			fprintf( stderr, "jumpi -> %d\n", jumpi );
			fprintf( stderr, "jumpj -> %d\n", jumpj );
			fprintf( stderr, "imid -> %d\n", imid );
			fprintf( stderr, "jmid -> %d\n", jmid );
#endif
//			fprintf( stderr, "#### FINAL i=%d, jumpi<imid=%d<%d, ist=%d, ien=%d, firstmp=%d\n", i, jumpi, imid, ist, ien, firstmp );

#if STOREWM
			break;
#else
			break;
#endif
		}
	}

//	fprintf( stderr, "imid = %d, but jumpi = %d\n", imid, jumpi );
//	fprintf( stderr, "jmid = %d, but jumpj = %d\n", jmid, jumpj );

//	for( j=0; j<lgth2; j++ ) midw[j] += currentw[j];
//	for( j=0; j<lgth2; j++ ) midm[j] += m[j+1];
//	for( j=0; j<lgth2; j++ ) midw[j] += WMMTX[imid][j];
//	for( j=0; j<lgth2; j++ ) midw[j] += WMMTX[imid][j];


#if STOREWM
	fprintf( stderr, "WMMTX = \n" );
    for( i=0; i<lgth1; i++ )
    {
        fprintf( stderr, "%d ", i );
        for( j=0; j<lgth2; j++ )
        {
            fprintf( stderr, "% 7.2f ", WMMTX[i][j] );
        }
        fprintf( stderr, "\n" );
    }
//	fprintf( stderr, "WMMTX2 = (p = %f)\n", fpenalty );
    for( i=0; i<lgth1; i++ )
    {
        fprintf( stderr, "%d ", i );
        for( j=0; j<lgth2+1; j++ )
        {
            fprintf( stderr, "% 7.2f ", WMMTX2[i][j] );
        }
        fprintf( stderr, "\n" );
    }

	fprintf( stderr, "jumpbacki = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpbacki[j] );
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "jumpbackj = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpbackj[j] );
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "jumpforwi = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpforwi[j] );
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "jumpforwj = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpforwj[j] );
	}
	fprintf( stderr, "\n" );


#endif


//	Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, cpmx1, cpmx2, ijp, icyc, jcyc );

#if 0 // irukamo
	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}
#endif



#if 0
	fprintf( stderr, "jumpi = %d, imid = %d\n", jumpi, imid );
	fprintf( stderr, "jumpj = %d, jmid = %d\n", jumpj, jmid );

	fprintf( stderr, "imid = %d\n", imid );
	fprintf( stderr, "jmid = %d\n", jmid );
#endif

	freearrays_rec1_variousdist
	(
		w1, w2, initverticalw, lastverticalw, midw, midm, midn,
		jumpbacki, jumpbackj, jumpforwi, jumpforwj, jumpdummi, jumpdummj,
		m, mp,
		doublework, intwork
#if STOREWM
		, WMMTX, WMMTX2
#endif
	);


//	fprintf( stderr, "==== calling myself (first)\n" );

	value = MSalignmm_rec_variousdist( matrices, icyc, jcyc, seq1, seq2, cpmx1s, cpmx2s, ist, ist+jumpi, jst, jst+jumpj, alloclen, fulllen1, fulllen2, aseq1, aseq2, depth, gapinfo, NULL, 0, NULL, headgp, tailgp, headgapfreq1_g, headgapfreq2_g ); // chudan mada	
#if 0
		fprintf( stderr, "aseq1[0] = %s\n", aseq1[0] );
		fprintf( stderr, "aseq2[0] = %s\n", aseq2[0] );
#endif
#if MEMSAVE
#else
	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], aseq1[i] );
	for( i=0; i<jcyc; i++ ) strcpy( mseq2[i], aseq2[i] );
#endif

//	fprintf( stderr, "====(f) aseq1[0] (%d) = %s (%d-%d)\n", depth, aseq1[0], ist, ien );
//	fprintf( stderr, "====(f) aseq2[0] (%d) = %s (%d-%d)\n", depth, aseq2[0], jst, jen );

	len = strlen( mseq1[0] );
//	fprintf( stderr, "len = %d\n", len );
	l = jmid - jumpj - 1;
//	fprintf( stderr, "l=%d\n", l );
	if( l > 0 )
	{
//		for( i=0; i<l; i++ ) gaps[i] = '-'; gaps[i] = 0;
		for( i=0; i<l; i++ ) gaps[i] = *newgapstr; gaps[i] = 0;
		for( i=0; i<icyc; i++ ) 
		{
			strcat( mseq1[i], gaps );
			mseq1[i][len+l] = 0;
		}
		for( j=0; j<jcyc; j++ )
		{
			strncat( mseq2[j], seq2[j]+jst+jumpj+1, l );
			mseq2[j][len+l] = 0;
		}
//		fprintf( stderr, "penalizing (2) .. %f(%d), %f(%d)\n", ogcp2[jumpj+1], jumpj+1, fgcp2[jmid-1], jmid-1 );
		value +=  ( ogcp2[jumpj+1] + fgcp2[jmid-1] );
//		value += fpenalty;
	}
	len = strlen( mseq1[0] );
	l = imid - jumpi - 1;
//	fprintf( stderr, "l=%d\n", l );
	if( l > 0 )
	{
//		for( i=0; i<l; i++ ) gaps[i] = '-'; gaps[i] = 0;
		for( i=0; i<l; i++ ) gaps[i] = *newgapstr; gaps[i] = 0;
		for( i=0; i<icyc; i++ )
		{
			strncat( mseq1[i], seq1[i]+ist+jumpi+1, l );
			mseq1[i][len+l] = 0;
		}
		for( j=0; j<jcyc; j++ ) 
		{
			strcat( mseq2[j], gaps );
			mseq2[j][len+l] = 0;
		}

//		for( i=0; i<lgth1; i++ ) fprintf( stderr, "ogcp1[%d] = %f\n", i, ogcp1[i] );
//		for( i=0; i<lgth1; i++ ) fprintf( stderr, "fgcp1[%d] = %f\n", i, fgcp1[i] );


//		fprintf( stderr, "penalizing (1) .. ogcp1[%d] = %f, fgcp1[%d] = %f\n", jumpi+1, ogcp1[jumpi+1], imid-1, fgcp1[imid-1] );
		value += ( ogcp1[jumpi+1] + fgcp1[imid-1] );
//		value += fpenalty;
	}
#if 0
	for( i=0; i<icyc; i++ ) fprintf( stderr, "after gapfill mseq1[%d]=%s\n", i, mseq1[i] );
	for( i=0; i<jcyc; i++ ) fprintf( stderr, "after gapfill mseq2[%d]=%s\n", i, mseq2[i] );
#endif

//	fprintf( stderr, "==== calling myself (second)\n" );

#if MEMSAVE
	alnlen = strlen( aseq1[0] );
	for( i=0; i<icyc; i++ ) aseq1[i] += alnlen;
	for( i=0; i<jcyc; i++ ) aseq2[i] += alnlen;
#endif

	value += MSalignmm_rec_variousdist( matrices, icyc, jcyc, seq1, seq2, cpmx1s, cpmx2s, ist+imid, ien, jst+jmid, jen, alloclen, fulllen1, fulllen2, aseq1, aseq2, depth, gapinfo, NULL, 0, NULL, headgp, tailgp, headgapfreq1_g, headgapfreq2_g ); // chudan mada
#if 0
		fprintf( stderr, "aseq1[0] = %s\n", aseq1[0] );
		fprintf( stderr, "aseq2[0] = %s\n", aseq2[0] );
#endif



#if DEBUG
	if( value - maxwm > 1 || maxwm - value > 1 )
	{
		fprintf( stderr, "WARNING value  = %f, but maxwm = %f\n", value, maxwm );
		for( i=0; i<icyc; i++ )
		{
			fprintf( stderr, ">1-%d\n%s\n", i, mseq1[i] );
			fprintf( stderr, "%s\n", aseq1[i] );
		}
		for( i=0; i<jcyc; i++ )
		{
			fprintf( stderr, ">2-%d\n%s\n", i, mseq2[i] );
			fprintf( stderr, "%s\n", aseq2[i] );
		}

//		exit( 1 );
	}
	else
	{
		fprintf( stderr, "value = %.0f, maxwm = %.0f -> ok\n", value, maxwm );
	}
#endif

#if MEMSAVE
#else
	for( i=0; i<icyc; i++ ) strcat( mseq1[i], aseq1[i] );
	for( i=0; i<jcyc; i++ ) strcat( mseq2[i], aseq2[i] );
#endif

//	fprintf( stderr, "====(s) aseq1[0] (%d) = \n%s\n (a%d-a%d)\n", depth, mseq1[0], ist, ien );
//	fprintf( stderr, "====(s) aseq2[0] (%d) = \n%s\n (b%d-b%d)\n", depth, mseq2[0], jst, jen );

	freearrays_rec2( gaps, aseq1, aseq2 );

#if 0
	if( seqlen( seq1[0] ) != nglen1 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen1) seqlen(seq1[0])=%d but nglen1=%d\n", seqlen( seq1[0] ), nglen1 );
		fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
		exit( 1 );
	}
	else
		fprintf( stderr, "nglen1 is ok in _rec\n" );
	if( seqlen( seq2[0] ) != nglen2 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen2) seqlen(seq2[0])=%d but nglen2=%d\n", seqlen( seq2[0] ), nglen2 );
		exit( 1 );
	}
	else
		fprintf( stderr, "nglen2 is ok in _rec\n" );
#endif

	return( value );
}

static void freearrays_variousdist( 
	double *ogcp1, 
	double *ogcp2, 
	double *fgcp1, 
	double *fgcp2, 
	double ***cpmx1s,
	double ***cpmx2s, 
	double *gapfreq1f,
	double *gapfreq2f,
	double **gapinfo, 
	char **mseq1, 
	char **mseq2 
)
{
	FreeFloatVec( ogcp1 );
	FreeFloatVec( ogcp2 );
	FreeFloatVec( fgcp1 );
	FreeFloatVec( fgcp2 );
	FreeFloatCub( cpmx1s );
	FreeFloatCub( cpmx2s );
	FreeFloatVec( gapfreq1f );
	FreeFloatVec( gapfreq2f );
	free( (void *)gapinfo );

	FreeCharMtx( mseq1 );
	FreeCharMtx( mseq2 );
}


double MSalignmm_variousdist( double **pairoffset, double ***matrices, double **dummy_mtx, char **seq1, char **seq2, double *eff1, double *eff2, double **eff1s, double **eff2s, int icyc, int jcyc, int alloclen, char *sgap1, char *sgap2, char *egap1, char *egap2, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
//	fprintf( stderr, "IN MSalignmm_variousdist\n" );
//	int k;
	int i, j, c;
	int ll1, ll2;
	int lgth1, lgth2;
	double wm = 0.0;   /* int ?????? */
	char **mseq1;
	char **mseq2;
	double *ogcp1;
	double *ogcp2;
	double *fgcp1;
	double *fgcp2;
	double ***cpmx1s;
	double ***cpmx2s;
	double **gapinfo;
	double fpenalty = (double)penalty;
	double *gapfreq1f;
	double *gapfreq2f;
	int nglen1, nglen2;
	double headgapfreq1;
	double headgapfreq2;

#if 0
	fprintf( stderr, "eff in SA+++align\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "eff1[%d] = %f\n", i, eff1[i] );
#endif

	nglen1 = seqlen( seq1[0] );
	nglen2 = seqlen( seq2[0] );

#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "seq1[%d] at root = %s\n", i, seq1[i] );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, "seq2[%d] at root = %s\n", j, seq2[j] );
	fprintf( stderr, "\n" );
#endif

	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );

	ll1 = ( (int)(lgth1) ) + 100;
	ll2 = ( (int)(lgth2) ) + 100;

	mseq1 = AllocateCharMtx( icyc, ll1+ll2 );
	mseq2 = AllocateCharMtx( jcyc, ll1+ll2 );

	gapinfo = AllocateFloatMtx( 6, 0 );
	ogcp1 = AllocateFloatVec( ll1+2 );
	ogcp2 = AllocateFloatVec( ll2+2 );
	fgcp1 = AllocateFloatVec( ll1+2 );
	fgcp2 = AllocateFloatVec( ll2+2 );


	cpmx1s = AllocateFloatCub( maxdistclass, ll1+2, nalphabets+1 );
	cpmx2s = AllocateFloatCub( maxdistclass, ll2+2, nalphabets+1 );

	gapfreq1f = AllocateFloatVec( ll1+2 ); // must be filled with 0.0
	gapfreq2f = AllocateFloatVec( ll2+2 ); // must be filled with 0.0

	for( i=0; i<icyc; i++ ) 
	{
		if( strlen( seq1[i] ) != lgth1 )
		{
			fprintf( stderr, "i = %d / %d\n", i, icyc );
			fprintf( stderr, "bug! hairetsu ga kowareta!\n" );
			exit( 1 );
		}
	}
	for( j=0; j<jcyc; j++ )
	{
		if( strlen( seq2[j] ) != lgth2 )
		{
			fprintf( stderr, "j = %d / %d\n", j, jcyc );
			fprintf( stderr, "bug! hairetsu ga kowareta!\n" );
			exit( 1 );
		}
	}

	for( c=0; c<maxdistclass; c++ )
	{
		MScpmx_calc_new( seq1, cpmx1s[c], eff1s[c], lgth1, icyc );
		MScpmx_calc_new( seq2, cpmx2s[c], eff2s[c], lgth2, jcyc );
	}


#if 1

	if( sgap1 )
	{
		new_OpeningGapCount( ogcp1, icyc, seq1, eff1, lgth1, sgap1 );
		new_OpeningGapCount( ogcp2, jcyc, seq2, eff2, lgth2, sgap2 );
		new_FinalGapCount( fgcp1, icyc, seq1, eff1, lgth1, egap2 );
		new_FinalGapCount( fgcp2, jcyc, seq2, eff2, lgth2, egap2 );
		outgapcount( &headgapfreq1, icyc, sgap1, eff1 );
		outgapcount( &headgapfreq2, jcyc, sgap2, eff2 );
		outgapcount( gapfreq1f+lgth1, icyc, egap1, eff1 );
		outgapcount( gapfreq2f+lgth2, jcyc, egap2, eff2 );
	}
	else
	{
		st_OpeningGapCount( ogcp1, icyc, seq1, eff1, lgth1 );
		st_OpeningGapCount( ogcp2, jcyc, seq2, eff2, lgth2 );
		st_FinalGapCount( fgcp1, icyc, seq1, eff1, lgth1 );
		st_FinalGapCount( fgcp2, jcyc, seq2, eff2, lgth2 );
		headgapfreq1 = 0.0;
		headgapfreq2 = 0.0;
		gapfreq1f[lgth1] = 0.0;
		gapfreq2f[lgth2] = 0.0;
	}

	if( legacygapcost == 0 )
	{
		gapcountf( gapfreq1f, seq1, icyc, eff1, lgth1 );
		gapcountf( gapfreq2f, seq2, jcyc, eff2, lgth2 );
		for( i=0; i<lgth1+1; i++ ) gapfreq1f[i] = 1.0 - gapfreq1f[i];
		for( i=0; i<lgth2+1; i++ ) gapfreq2f[i] = 1.0 - gapfreq2f[i];
		headgapfreq1 = 1.0 - headgapfreq1;
		headgapfreq2 = 1.0 - headgapfreq2;
	}
	else
	{
		for( i=0; i<lgth1+1; i++ ) gapfreq1f[i] = 1.0;
		for( i=0; i<lgth2+1; i++ ) gapfreq2f[i] = 1.0;
		headgapfreq1 = 1.0;
		headgapfreq2 = 1.0;
	}

#if 1
	for( i=0; i<lgth1; i++ ) 
	{
		ogcp1[i] = 0.5 * ( 1.0 - ogcp1[i] ) * fpenalty * ( gapfreq1f[i] );
		fgcp1[i] = 0.5 * ( 1.0 - fgcp1[i] ) * fpenalty * ( gapfreq1f[i] );
//		fprintf( stderr, "fgcp1[%d] = %f\n", i, fgcp1[i] );
	}
	for( i=0; i<lgth2; i++ ) 
	{
		ogcp2[i] = 0.5 * ( 1.0 - ogcp2[i] ) * fpenalty * ( gapfreq2f[i] );
		fgcp2[i] = 0.5 * ( 1.0 - fgcp2[i] ) * fpenalty * ( gapfreq2f[i] );
//		fprintf( stderr, "fgcp2[%d] = %f\n", i, fgcp2[i] );
	}
#else
	for( i=0; i<lgth1; i++ ) 
	{
		ogcp1[i] = 0.5 * fpenalty;
		fgcp1[i] = 0.5 * fpenalty;
	}
	for( i=0; i<lgth2; i++ ) 
	{
		ogcp2[i] = 0.5 * fpenalty;
		fgcp2[i] = 0.5 * fpenalty;
	}
#endif

	gapinfo[0] = ogcp1;
	gapinfo[1] = fgcp1;
	gapinfo[2] = ogcp2;
	gapinfo[3] = fgcp2;
	gapinfo[4] = gapfreq1f;
	gapinfo[5] = gapfreq2f;
#endif

#if 0
	fprintf( stdout, "in MSalignmm.c\n" );
	for( i=0; i<icyc; i++ )
	{
		fprintf( stdout, ">%d of GROUP1\n", i );
		fprintf( stdout, "%s\n", seq1[i] );
	}
	for( i=0; i<jcyc; i++ )
	{
		fprintf( stdout, ">%d of GROUP2\n", i );
		fprintf( stdout, "%s\n", seq2[i] );
	}
	fflush( stdout );
#endif

	wm = MSalignmm_rec_variousdist( matrices, icyc, jcyc, seq1, seq2, cpmx1s, cpmx2s, 0, lgth1-1, 0, lgth2-1, alloclen, lgth1, lgth2, mseq1, mseq2, 0, gapinfo, chudanpt, chudanref, chudanres, headgp, tailgp, headgapfreq1, headgapfreq2 );
#ifdef enablemultithread
	if( chudanres && *chudanres ) 
	{
//		fprintf( stderr, "\n\n## CHUUDAN!!! relay\n" );
		*chudanres = 1;
		freearrays_variousdist( ogcp1, ogcp2, fgcp1, fgcp2, cpmx1s, cpmx2s, gapfreq1f, gapfreq2f, gapinfo, mseq1, mseq2 );
		return( -1.0 );
	}
#endif

#if 0
		fprintf( stderr, "\n" );
		fprintf( stderr, " seq1[0] = %s\n", seq1[0] );
		fprintf( stderr, " seq2[0] = %s\n", seq2[0] );
		fprintf( stderr, "mseq1[0] = %s\n", mseq1[0] );
		fprintf( stderr, "mseq2[0] = %s\n", mseq2[0] );
		fprintf( stderr, "\n" );
#endif

//	fprintf( stderr, "wm = %f\n", wm );


	for( i=0; i<icyc; i++ ) strcpy( seq1[i], mseq1[i] );
	for( i=0; i<jcyc; i++ ) strcpy( seq2[i], mseq2[i] );

	if( seqlen( seq1[0] ) != nglen1 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen1) seqlen(seq1[0])=%d but nglen1=%d\n", seqlen( seq1[0] ), nglen1 );
		fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
		exit( 1 );
	}
	if( seqlen( seq2[0] ) != nglen2 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen2) seqlen(seq2[0])=%d but nglen2=%d\n", seqlen( seq2[0] ), nglen2 );
		exit( 1 );
	}


	freearrays_variousdist( ogcp1, ogcp2, fgcp1, fgcp2, cpmx1s, cpmx2s, gapfreq1f, gapfreq2f, gapinfo, mseq1, mseq2 );

	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );
	for( i=0; i<icyc; i++ ) 
	{
		if( strlen( seq1[i] ) != lgth1 )
		{
			fprintf( stderr, "i = %d / %d\n", i, icyc );
			fprintf( stderr, "hairetsu ga kowareta (end of MSalignmm) !\n" );
			exit( 1 );
		}
	}
	for( j=0; j<jcyc; j++ )
	{
		if( strlen( seq2[j] ) != lgth2 )
		{
			fprintf( stderr, "j = %d / %d\n", j, jcyc );
			fprintf( stderr, "hairetsu ga kowareta (end of MSalignmm) !\n" );
			exit( 1 );
		}
	}

#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, " seq1[i] = %s\n", seq1[i] );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, " seq2[j] = %s\n", seq2[j] );
	fprintf( stderr, "\n" );
#endif

	return( wm );
}
