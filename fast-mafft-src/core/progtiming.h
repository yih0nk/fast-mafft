/* Optional instrumentation for the progressive stage (fast-mafft, work/prog).
 * Compiled only with -DPROGTIMING; the default build is unaffected. */
#ifndef PROGTIMING_H
#define PROGTIMING_H
#ifdef PROGTIMING
#include <time.h>
static inline double pt_now( void )
{
	struct timespec ts;
	clock_gettime( CLOCK_MONOTONIC, &ts );
	return( ts.tv_sec + ts.tv_nsec * 1e-9 );
}
#define PT_NACC 32
extern TLS double pt_acc[PT_NACC];
#define PT_DECL double _pt_t = 0.0
#define PT_START() ( _pt_t = pt_now() )
#define PT_ADD(k) do { double _pt_n = pt_now(); pt_acc[k] += _pt_n - _pt_t; _pt_t = _pt_n; } while( 0 )
#else
#define PT_DECL
#define PT_START()
#define PT_ADD(k)
#endif
#endif
