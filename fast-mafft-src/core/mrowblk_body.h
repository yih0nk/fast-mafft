/* work/prog3: body of Adp_mrow_blk (template by inclusion, see adpvec_body.h); define ADPM_NAME,
   ADPM_ATTR.  Match scores of row i for the columns [ja, jb), into curw:
     curw[j] = chain over the letters u in ( rowmask & bm[j/4] ), ascending, of  a += s[u] * c[u][j]
   starting from a = 0.0.  rowmask holds the letters with a nonzero row score, bm[b] the letters with
   a nonzero weight in some column 4b..4b+3 of the profile.  Every letter with a nonzero weight in
   column j and a nonzero row score is in the chain, in ascending order, as in match_calc(); the
   others add s * 0.0 or 0.0 * w, i.e. +-0.0, to a sum that is never -0.0: exact no-ops (see
   match_calc_dense).  Four columns of a block share one chain of vector operations. */
ADPM_ATTR
static void ADPM_NAME( int ja, int jb, unsigned int rowmask, const double *su, const double * const *cu,
	const unsigned int *bm, double * restrict curw )
{
	typedef double v4_t __attribute__(( vector_size( 32 ) ));
	int j = ja;
	for( ; j<jb && ( j & 3 ); j++ )
	{
		unsigned int m = rowmask & bm[j>>2];
		double a = 0.0;
		while( m ) { int u = __builtin_ctz( m ); m &= m - 1; a += su[u] * cu[u][j]; }
		curw[j] = a;
	}
	/* four or two blocks at a time: independent chains (the add latency, not the throughput, limits a
	   single chain); the letters of all the blocks run over each, which only adds more exact zeros */
	for( ; j+16<=jb; j+=16 )
	{
		unsigned int m = rowmask & ( bm[j>>2] | bm[(j>>2)+1] | bm[(j>>2)+2] | bm[(j>>2)+3] );
		v4_t a = (v4_t){0}, b = (v4_t){0}, e = (v4_t){0}, f = (v4_t){0}; // +0.0
		while( m )
		{
			int u = __builtin_ctz( m );
			const double *cp = cu[u] + j;
			v4_t c, d, g, h, s = su[u] - (v4_t){0}; // splat, exact (x - +0.0 == x)
			m &= m - 1;
			memcpy( &c, cp, sizeof( v4_t ) );
			memcpy( &d, cp + 4, sizeof( v4_t ) );
			memcpy( &g, cp + 8, sizeof( v4_t ) );
			memcpy( &h, cp + 12, sizeof( v4_t ) );
			a += s * c;
			b += s * d;
			e += s * g;
			f += s * h;
		}
		memcpy( curw + j, &a, sizeof( v4_t ) );
		memcpy( curw + j + 4, &b, sizeof( v4_t ) );
		memcpy( curw + j + 8, &e, sizeof( v4_t ) );
		memcpy( curw + j + 12, &f, sizeof( v4_t ) );
	}
	for( ; j+8<=jb; j+=8 )
	{
		unsigned int m = rowmask & ( bm[j>>2] | bm[(j>>2)+1] );
		v4_t a = (v4_t){0}, b = (v4_t){0}; // +0.0
		while( m )
		{
			int u = __builtin_ctz( m );
			const double *cp = cu[u] + j;
			v4_t c, d, s = su[u] - (v4_t){0}; // splat, exact (x - +0.0 == x)
			m &= m - 1;
			memcpy( &c, cp, sizeof( v4_t ) );
			memcpy( &d, cp + 4, sizeof( v4_t ) );
			a += s * c;
			b += s * d;
		}
		memcpy( curw + j, &a, sizeof( v4_t ) );
		memcpy( curw + j + 4, &b, sizeof( v4_t ) );
	}
	for( ; j+4<=jb; j+=4 )
	{
		unsigned int m = rowmask & bm[j>>2];
		v4_t a = (v4_t){0}; // +0.0
		while( m )
		{
			int u = __builtin_ctz( m );
			v4_t c, s = su[u] - (v4_t){0}; // splat, exact (x - +0.0 == x)
			m &= m - 1;
			memcpy( &c, cu[u] + j, sizeof( v4_t ) );
			a += s * c;
		}
		memcpy( curw + j, &a, sizeof( v4_t ) );
	}
	for( ; j<jb; j++ )
	{
		unsigned int m = rowmask & bm[j>>2];
		double a = 0.0;
		while( m ) { int u = __builtin_ctz( m ); m &= m - 1; a += su[u] * cu[u][j]; }
		curw[j] = a;
	}
}
