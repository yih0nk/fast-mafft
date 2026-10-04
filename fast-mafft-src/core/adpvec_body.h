/* work/prog3: body of the vectorized A__align row kernel; included from Salignmm.c once per
   vector width (template by inclusion).  Define before including:
     ADPV_NAME   function name
     ADPV_VB     vector size in bytes (16: 2 doubles, 32: 4 doubles)
     ADPV_ATTR   function attributes (e.g. __attribute__((target("avx2"))) or empty)
     ADPV_MROW   0: the match score is formed inline from 4 letters (ADP_CELL( 0, 1 ));
                 1: it was stored in curw[j] beforehand (ADP_CELL( 0, 0 ), the general case)

   It computes the cells j in [ja, jb) of row i (all of them "inner" cells: 2 <= ja, jb <= lgth2,
   so every cell has a match score and none is the row's first cell) exactly as ADP_CELL( 0, 1 )
   (or ADP_CELL( 0, 0 ) with ADPV_MROW) does, V cells per step, and the remaining jb-ja mod V
   cells with that macro itself.

   Why this is bit-identical to the scalar loop:
   - Every floating-point value of a cell (the match sum, g3, g4, g2, g2 + fpenalty_ex, gg, the new
     m[j], currentw[j]) is produced by the same operations on the same operands in the same order,
     one lane per cell, written as the same C statements as in ADP_CELL (so a compiler that
     contracts a*b+c into an fma, or does not, treats both alike).
   - The only dependence between cells of a row is the row carry (mi, mpi), updated per cell by
         if( g2c >= mi ) { mi = g2c; mpi = j-1; }
     i.e. "keep the later element on ties" selection over the sequence (carry, g2c[ja], g2c[ja+1], ...).
     That selection is associative: combining (a, b) gives b if b.v >= a.v else a, and both
     groupings of (a, b, c) give the latest element whose value is the maximum (>= is a total
     preorder on non-NaN doubles; equal values, including +0.0 and -0.0, count as ties and the
     later one wins in both).  So a prefix scan inside the vector (Hillis-Steele: shifts by 1, 2),
     combined with the carry, yields for every lane exactly the (mi, mpi) the scalar loop has before
     and after that cell.  Values are only selected, never recomputed, so even their bits match.
   - Selections use comparisons with the same operand order and strictness (>, >=) as ADP_CELL.
   - Indices (mpi, j, i, mp[j]) are carried as doubles; they are integers far below 2^53, so the
     differences -(j - mpi) and i - mp[j] are exact and convert back to the same ints. */

#if ADPV_VB == 32
#define ADPV_V 4
#elif ADPV_VB == 16
#define ADPV_V 2
#else
#error ADPV_VB
#endif

ADPV_ATTR
static void ADPV_NAME( int i, int ja, int jb,
	double * restrict mloc, int * restrict mploc, double * restrict curw, int * restrict ijpi,
	const double * restrict prew, const double * restrict gf2, const double * restrict fg2, const double * restrict og2,
	const double * restrict c0, const double * restrict c1, const double * restrict c2, const double * restrict c3,
	const double *scu,
	double fgcp1va, double ogcp1va, double gf1va, double gf1vapre, double fpenalty_ex,
	double *mip, int *mpip, double *wmp )
{
	typedef double vd_t __attribute__(( vector_size( ADPV_VB ) ));
	typedef __typeof__( (vd_t){0} >= (vd_t){0} ) vl_t; // the comparison result type (64-bit lanes)
	typedef int vi_t __attribute__(( vector_size( ADPV_VB / 2 ) ));
#define VLD( p ) ( { vd_t _v; memcpy( &_v, (p), sizeof( vd_t ) ); _v; } )
#define VST( p, v ) do { vd_t _v = (v); memcpy( (p), &_v, sizeof( vd_t ) ); } while( 0 )
#define VLDI( p ) ( { vi_t _v; memcpy( &_v, (p), sizeof( vi_t ) ); __builtin_convertvector( _v, vd_t ); } )
#define VSTI( p, v ) do { vi_t _v = __builtin_convertvector( (v), vi_t ); memcpy( (p), &_v, sizeof( vi_t ) ); } while( 0 )
#define VSEL( m, a, b ) ( (vd_t)( ( (vl_t)(a) & (m) ) | ( (vl_t)(b) & ~(m) ) ) )
#define VSPLAT( x ) ( (x) - (vd_t){0} ) // x - (+0.0) == x for every x, -0.0 included
#if ADPV_MROW
	const double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0; // (unused: the match scores are in curw)
#else
	const double s0 = scu[0], s1 = scu[1], s2 = scu[2], s3 = scu[3];
#endif
	double mil = *mip, wmv = *wmp;
	int mpil = *mpip, ijv, j = ja;

	if( jb - ja >= ADPV_V )
	{
		const vd_t vs0 = VSPLAT( s0 ), vs1 = VSPLAT( s1 ), vs2 = VSPLAT( s2 ), vs3 = VSPLAT( s3 );
		const vd_t vzero = (vd_t){0}; // +0.0 in every lane
		const vd_t vfg1 = VSPLAT( fgcp1va ), vog1 = VSPLAT( ogcp1va ), vgf1 = VSPLAT( gf1va ), vgf1p = VSPLAT( gf1vapre );
		const vd_t vfpe = VSPLAT( fpenalty_ex );
		const vd_t vi = VSPLAT( (double)i ), vim1 = VSPLAT( (double)( i - 1 ) );
		vd_t cv = VSPLAT( mil ), cx = VSPLAT( (double)mpil ); // the row carry (mi, mpi) in every lane
		vd_t jv, wlast = vzero;
		int k;
		for( k=0; k<ADPV_V; k++ ) jv[k] = (double)( ja + k );
		for( ; j + ADPV_V <= jb; j += ADPV_V )
		{
			vd_t mj = VLD( mloc + j ), pw = VLD( prew + j - 1 ), mpj = VLDI( mploc + j );
			vd_t base, g3, g4, g2, g2c, gg, wv, ijvv, ix, av, ax, sv, sx, pv, px, ev, ex, bv, bx;
			vl_t t;
#if ADPV_MROW
			base = VLD( curw + j );
#else
			base = vzero;
			base += vs0 * VLD( c0 + j );
			base += vs1 * VLD( c1 + j );
			base += vs2 * VLD( c2 + j );
			base += vs3 * VLD( c3 + j );
#endif
			g3 = mj + vfg1 * VLD( gf2 + j );
			g4 = pw + vog1 * VLD( gf2 + j - 1 );
			t = ( g4 >= mj );
			VST( mloc + j, VSEL( t, g4, mj ) + vfpe );
			VSTI( mploc + j, VSEL( t, vim1, mpj ) );
			g2 = pw + VLD( og2 + j ) * vgf1p;
			g2c = g2 + vfpe;
			/* inclusive prefix selection of (g2c, j-1) inside the vector; a lane combined with itself
			   is itself, so the shifts duplicate the low lanes instead of inserting an identity */
			ix = jv - 1.0;
#if ADPV_V == 4
			av = __builtin_shufflevector( g2c, g2c, 0, 0, 1, 2 );
			ax = __builtin_shufflevector( ix, ix, 0, 0, 1, 2 );
			t = ( g2c >= av );
			sv = VSEL( t, g2c, av ); sx = VSEL( t, ix, ax );
			av = __builtin_shufflevector( sv, sv, 0, 1, 0, 1 );
			ax = __builtin_shufflevector( sx, sx, 0, 1, 0, 1 );
			t = ( sv >= av );
			sv = VSEL( t, sv, av ); sx = VSEL( t, sx, ax );
			/* exclusive: lane k gets carry (+) lanes 0..k-1; lane 0 combines the carry with itself */
			pv = __builtin_shufflevector( cv, sv, 0, 4, 5, 6 );
			px = __builtin_shufflevector( cx, sx, 0, 4, 5, 6 );
			bv = __builtin_shufflevector( sv, sv, 3, 3, 3, 3 );
			bx = __builtin_shufflevector( sx, sx, 3, 3, 3, 3 );
#else
			av = __builtin_shufflevector( g2c, g2c, 0, 0 );
			ax = __builtin_shufflevector( ix, ix, 0, 0 );
			t = ( g2c >= av );
			sv = VSEL( t, g2c, av ); sx = VSEL( t, ix, ax );
			pv = __builtin_shufflevector( cv, sv, 0, 2 );
			px = __builtin_shufflevector( cx, sx, 0, 2 );
			bv = __builtin_shufflevector( sv, sv, 1, 1 );
			bx = __builtin_shufflevector( sx, sx, 1, 1 );
#endif
			t = ( pv >= cv );
			ev = VSEL( t, pv, cv ); ex = VSEL( t, px, cx ); // (mi, mpi) before each lane's cell
			t = ( bv >= cv );
			cv = VSEL( t, bv, cv ); cx = VSEL( t, bx, cx ); // after the last lane
			/* the cell: wm = previous row's value, then the row-gap and column-gap candidates */
			gg = ev + VLD( fg2 + j - 1 ) * vgf1;
			t = ( gg > pw );
			wv = VSEL( t, gg, pw );
			ijvv = VSEL( t, ex - jv, vzero );       // -( j - mpi )
			t = ( g3 > wv );
			ijvv = VSEL( t, vi - mpj, ijvv );       // +( i - mp[j] )
			wv = VSEL( t, g3, wv );
			VST( curw + j, base + wv );
			VSTI( ijpi + j, ijvv );
			wlast = wv;
			jv += (double)ADPV_V;
		}
		mil = cv[0];
		mpil = (int)cx[0];
		wmv = wlast[ADPV_V-1];
	}
	for( ; j<jb; j++ )
#if ADPV_MROW
		ADP_CELL( 0, 0 )
#else
		ADP_CELL( 0, 1 )
#endif
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

#undef ADPV_V
