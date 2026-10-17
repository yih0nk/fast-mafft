#include "mltaln.h"
#include "mtxutl.h"

/*
  from "C gengo niyoru saishin algorithm jiten" ISBN4-87408-414-1 Haruhiko Okumura
*/
static void make_sintbl(int n, double sintbl[])
{
        int i, n2, n4, n8;
        double c, s, dc, ds, t;

        n2 = n / 2;  n4 = n / 4;  n8 = n / 8;
        t = sin(PI / n);
        dc = 2 * t * t;  ds = sqrt(dc * (2 - dc));
        t = 2 * dc;  c = sintbl[n4] = 1;  s = sintbl[0] = 0;
        for (i = 1; i < n8; i++) {
                c -= dc;  dc += t * c;
                s += ds;  ds -= t * s;
                sintbl[i] = s;  sintbl[n4 - i] = c;
        }
        if (n8 != 0) sintbl[n8] = sqrt(0.5);
        for (i = 0; i < n4; i++)
                sintbl[n2 - i] = sintbl[i];
        for (i = 0; i < n2 + n4; i++)
                sintbl[i + n2] = - sintbl[i];
}
/*
  {\tt fft()}.
*/
static void make_bitrev(int n, int bitrev[])
{
        int i, j, k, n2;

        n2 = n / 2;  i = j = 0;
        for ( ; ; ) {
                bitrev[i] = j;
                if (++i >= n) break;
                k = n2;
                while (k <= j) {  j -= k;  k /= 2;  }
                j += k;
        }
}
/*
*/
/* work/v10w: for long transforms (n >= FFT_REORDER_MIN) each stage runs its butterflies block by
   block (i0 = 0, k2, 2*k2, ...; j = 0..k-1 inside, with the stage's twiddles copied to two contiguous
   arrays) instead of twiddle by twiddle with stride k2 through the whole vector.  The butterflies of a
   stage are independent and touch disjoint elements; each one is the same statements with the same
   c and s (sintbl[j*d + n4], +-sintbl[j*d]: h counts j*d in both orders), so every element gets the
   same value.  2.5-3.5x faster for n = 32768..131072 (the original's stride-k2 sweeps miss the cache);
   short transforms keep the original order.  MAFFT_FFTCHECK=1 reruns every long transform the
   original way on a copy and compares bitwise. */
#define FFT_REORDER_MIN 16384
#define FFT_REORDER_KSMALL 4 // stages with k < this keep the original order (a few sequential sweeps)
static int fft_body(int n, Fukusosuu *x, int freeflag, int reorder);
static int fftcheck_v = -1;
static long fftcheck_n = 0;
static void fftcheck_report( void ) __attribute__(( destructor ));
static void fftcheck_report( void ) { if( fftcheck_v > 0 ) reporterr( "\nFFTCHECK: %ld long transforms verified, no mismatch\n", fftcheck_n ); }
int fft(int n, Fukusosuu *x, int freeflag)
{
	int an = n < 0 ? -n : n;
	if( fftcheck_v < 0 ) fftcheck_v = ( getenv( "MAFFT_FFTCHECK" ) != NULL );
	if( !freeflag && an >= FFT_REORDER_MIN && fftcheck_v > 0 )
	{
		Fukusosuu *ref = malloc( sizeof( Fukusosuu ) * an );
		int r;
		if( ref == NULL ) ErrorExit( "Cannot allocate in fft.\n" );
		memcpy( ref, x, sizeof( Fukusosuu ) * an );
		fft_body( n, ref, 0, 0 );
		r = fft_body( n, x, 0, 1 );
		if( memcmp( ref, x, sizeof( Fukusosuu ) * an ) ) { reporterr( "\nFFTCHECK MISMATCH (n=%d)\n", n ); exit( 1 ); }
		free( ref );
		__atomic_add_fetch( &fftcheck_n, 1, __ATOMIC_RELAXED );
		return( r );
	}
	return( fft_body( n, x, freeflag, an >= FFT_REORDER_MIN ) );
}
static int fft_body(int n, Fukusosuu *x, int freeflag, int reorder)
{
        static TLS int    last_n = 0;    /*  {\tt n} */
        static TLS int   *bitrev = NULL; /*  */
        static TLS double *sintbl = NULL; /*  */
        static TLS double *twc = NULL, *tws = NULL; /* work/v10w: the twiddles of a stage, contiguous */
        static TLS int twcap = 0;
        int i, j, k, ik, h, d, k2, n4, inverse, i0;
        double t, s, c, dR, dI;

		if (freeflag)
		{
			if (bitrev) free(bitrev); bitrev = NULL;
			if (sintbl) free(sintbl); sintbl = NULL;
			last_n = 0;
			free( twc ); free( tws ); twc = tws = NULL; twcap = 0; // work/v10w
			return( 0 );
		}

        /*  */
        if (n < 0) {
                n = -n;  inverse = 1;  /*  */
        } else inverse = 0;
        n4 = n / 4;
        if (n != last_n || n == 0) {
                last_n = n;
#if 0
                if (sintbl != NULL) {
					free(sintbl);
					sintbl = NULL;
				}
                if (bitrev != NULL) {
					free(bitrev);
					bitrev = NULL;
				}
                if (n == 0) return 0;  /*  */
                sintbl = (double *)malloc((n + n4) * sizeof(double));
                bitrev = (int *)malloc(n * sizeof(int));
#else /* by T. Nishiyama */
				sintbl = realloc(sintbl, (n + n4) * sizeof(double));
				bitrev = realloc(bitrev, n * sizeof(int));
#endif
                if (sintbl == NULL || bitrev == NULL) {
                        fprintf(stderr, "\n");  return 1;
                }
                make_sintbl(n, sintbl);
                make_bitrev(n, bitrev);
        }
        for (i = 0; i < n; i++) {    /*  */
                j = bitrev[i];
                if (i < j) {
                        t = x[i].R;  x[i].R = x[j].R;  x[j].R = t;
                        t = x[i].I;  x[i].I = x[j].I;  x[j].I = t;
                }
        }
        if (reorder && twcap < n / 2) { /* work/v10w */
                free( twc ); free( tws );
                twcap = n / 2;
                twc = malloc( sizeof( double ) * twcap ); tws = malloc( sizeof( double ) * twcap );
                if (twc == NULL || tws == NULL) { fprintf(stderr, "\n");  return 1; }
        }
        for (k = 1; k < n; k = k2) {    /*  */
#if 0
				fprintf( stderr, "%d / %d\n", k, n );
#endif
                if (reorder && k >= FFT_REORDER_KSMALL) { /* work/v10w: see above */
                        k2 = k + k;  d = n / k2;
                        for (j = 0, h = 0; j < k; j++, h += d) {
                                twc[j] = sintbl[h + n4];
                                if (inverse) tws[j] = - sintbl[h];
                                else         tws[j] =   sintbl[h];
                        }
                        for (i0 = 0; i0 < n; i0 += k2) {
                                Fukusosuu *xa = x + i0, *xb = x + i0 + k;
                                for (j = 0; j < k; j++) {
                                        c = twc[j];  s = tws[j];
                                        dR = s * xb[j].I + c * xb[j].R;
                                        dI = c * xb[j].I - s * xb[j].R;
                                        xb[j].R = xa[j].R - dR;  xa[j].R += dR;
                                        xb[j].I = xa[j].I - dI;  xa[j].I += dI;
                                }
                        }
                        continue;
                }
                h = 0;  k2 = k + k;  d = n / k2;
                for (j = 0; j < k; j++) {
#if 0
					if( j % 1 == 0 )
						fprintf( stderr, "%d / %d\r", j, k );
#endif
                        c = sintbl[h + n4];
                        if (inverse) s = - sintbl[h];
                        else         s =   sintbl[h];
                        for (i = j; i < n; i += k2) {
#if 0
								if( k>=4194000 ) fprintf( stderr, "in loop %d -  %d < %d, k2=%d\r", j, i, n, k2 );
#endif
                                ik = i + k;
                                dR = s * x[ik].I + c * x[ik].R;
                                dI = c * x[ik].I - s * x[ik].R;
                                x[ik].R = x[i].R - dR;  x[i].R += dR;
                                x[ik].I = x[i].I - dI;  x[i].I += dI;
                        }
                        h += d;
                }
        }
        if (! inverse)    /* n */
                for (i = 0; i < n; i++) {  x[i].R /= n;  x[i].I /= n;  }
        return 0;  /*  */
}
