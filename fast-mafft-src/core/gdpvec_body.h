/* work/prog3: body of the vectorized G__align11 row kernel (template by inclusion; define
   GDPV_NAME, GDPV_VB (16 or 32), GDPV_ATTR).  Cells j in [ja, jb) of row i, all inner cells
   (2 <= ja, jb <= lgth2-... every cell has a match score, none is the row's first), exactly as
   GDP_CELL( 0, 1 ) with the match score read from prow[j] (= mrow[s2[j]], copied beforehand).
   The argument for exactness is that of adpvec_body.h: per-lane operations identical to the scalar
   cell, the row carry (mi, mpi) by an associative "later wins ties" prefix selection, indices as
   exact doubles.  (G__align11's scores are integers anyway.) */
#if GDPV_VB == 32
#define GDPV_V 4
#elif GDPV_VB == 16
#define GDPV_V 2
#else
#error GDPV_VB
#endif

GDPV_ATTR
static void GDPV_NAME( int i, int ja, int jb, const double * restrict prow,
	double * restrict mloc, int * restrict mploc, double * restrict curw, int * restrict ijpi,
	const double * restrict prew, double fpenalty, double fpenalty_ex, double fpenalty_ex_i,
	double *mip, int *mpip, double *wmp )
{
	typedef double vd_t __attribute__(( vector_size( GDPV_VB ) ));
	typedef int vi_t __attribute__(( vector_size( GDPV_VB / 2 ) ));
	typedef __typeof__( (vd_t){0} >= (vd_t){0} ) vl_t;
#define VLD( p ) ( { vd_t _v; memcpy( &_v, (p), sizeof( vd_t ) ); _v; } )
#define VST( p, v ) do { vd_t _v = (v); memcpy( (p), &_v, sizeof( vd_t ) ); } while( 0 )
#define VLDI( p ) ( { vi_t _v; memcpy( &_v, (p), sizeof( vi_t ) ); __builtin_convertvector( _v, vd_t ); } )
#define VSTI( p, v ) do { vi_t _v = __builtin_convertvector( (v), vi_t ); memcpy( (p), &_v, sizeof( vi_t ) ); } while( 0 )
#define VSEL( m, a, b ) ( (vd_t)( ( (vl_t)(a) & (m) ) | ( (vl_t)(b) & ~(m) ) ) )
#define VSPLAT( x ) ( (x) - (vd_t){0} )
	double mil = *mip, wmv = *wmp;
	int mpil = *mpip, j = ja;

	if( jb - ja >= GDPV_V )
	{
		const vd_t vzero = (vd_t){0};
		const vd_t vfp = VSPLAT( fpenalty ), vfpe = VSPLAT( fpenalty_ex ), vfpei = VSPLAT( fpenalty_ex_i );
		const vd_t vi = VSPLAT( (double)i ), vim1 = VSPLAT( (double)( i - 1 ) );
		vd_t cv = VSPLAT( mil ), cx = VSPLAT( (double)mpil );
		vd_t jv, wlast = vzero;
		int k;
		for( k=0; k<GDPV_V; k++ ) jv[k] = (double)( ja + k );
		for( ; j + GDPV_V <= jb; j += GDPV_V )
		{
			vd_t mj = VLD( mloc + j ), pw = VLD( prew + j - 1 ), mpj = VLDI( mploc + j );
			vd_t base = VLD( prow + j ), cw, mnew, gg, pc, wv, ijvv, ix, av, ax, sv, sx, pv, px, ev, ex, bv, bx;
			vl_t t;
			cw = mj + vfp;
			t = ( pw >= mj );
			mnew = VSEL( t, pw, mj );
			mnew += vfpe;
			VST( mloc + j, mnew );
			VSTI( mploc + j, VSEL( t, vim1, mpj ) );
			pc = pw + vfpei;
			ix = jv - 1.0;
#if GDPV_V == 4
			av = __builtin_shufflevector( pc, pc, 0, 0, 1, 2 );
			ax = __builtin_shufflevector( ix, ix, 0, 0, 1, 2 );
			t = ( pc >= av );
			sv = VSEL( t, pc, av ); sx = VSEL( t, ix, ax );
			av = __builtin_shufflevector( sv, sv, 0, 1, 0, 1 );
			ax = __builtin_shufflevector( sx, sx, 0, 1, 0, 1 );
			t = ( sv >= av );
			sv = VSEL( t, sv, av ); sx = VSEL( t, sx, ax );
			pv = __builtin_shufflevector( cv, sv, 0, 4, 5, 6 );
			px = __builtin_shufflevector( cx, sx, 0, 4, 5, 6 );
			bv = __builtin_shufflevector( sv, sv, 3, 3, 3, 3 );
			bx = __builtin_shufflevector( sx, sx, 3, 3, 3, 3 );
#else
			av = __builtin_shufflevector( pc, pc, 0, 0 );
			ax = __builtin_shufflevector( ix, ix, 0, 0 );
			t = ( pc >= av );
			sv = VSEL( t, pc, av ); sx = VSEL( t, ix, ax );
			pv = __builtin_shufflevector( cv, sv, 0, 2 );
			px = __builtin_shufflevector( cx, sx, 0, 2 );
			bv = __builtin_shufflevector( sv, sv, 1, 1 );
			bx = __builtin_shufflevector( sx, sx, 1, 1 );
#endif
			t = ( pv >= cv );
			ev = VSEL( t, pv, cv ); ex = VSEL( t, px, cx );
			t = ( bv >= cv );
			cv = VSEL( t, bv, cv ); cx = VSEL( t, bx, cx );
			gg = ev + vfp;
			t = ( gg > pw );
			wv = VSEL( t, gg, pw );
			ijvv = VSEL( t, ex - jv, vzero );       // -( j - mpi )
			t = ( cw > wv );
			ijvv = VSEL( t, vi - mpj, ijvv );       // i - mp[j]
			wv = VSEL( t, cw, wv );
			VST( curw + j, base + wv );
			VSTI( ijpi + j, ijvv );
			wlast = wv;
			jv += (double)GDPV_V;
		}
		mil = cv[0];
		mpil = (int)cx[0];
		wmv = wlast[GDPV_V-1];
	}
	for( ; j<jb; j++ )
		GDP_CELL_P( prow[j] )
	*mip = mil;
	*mpip = mpil;
	*wmp = wmv;
#undef VLD
#undef VST
#undef VLDI
#undef VSTI
#undef VSEL
#undef VSPLAT
}

#undef GDPV_V
