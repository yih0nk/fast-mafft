/* work/prog2: idle worker threads of the progressive pass help with the large alignments.
 *
 * The pass runs one pthread per --thread; a thread that finds no ready guide-tree step waits on
 * the scheduler's condition variable.  A thread that is about to run a large piece of work
 * (the DP of a big alignment, profile construction, row copies) can post a job; waiting
 * threads take a slot, call job->run( job->arg ) alongside the poster, and go back to waiting.
 * run() must be safe to call from any number of threads at once and must return when no work
 * is left that the caller could claim.  The results never depend on how many threads take
 * part: every value is computed exactly as in the one-thread code, only by another thread.
 *
 * Without a pool (dphelppool == NULL: other programs, in-order dispatch, no threads) posting
 * does nothing and the poster does all the work itself.
 */
#ifndef DPHELP_H
#define DPHELP_H
#ifdef enablemultithread

typedef struct _dphelpjob
{
	void (*run)( void *arg );
	void *arg;
	int slots;               // helpers that may still join (guarded by the pool mutex)
	int refs;                // helpers currently inside run() (atomic)
	struct _dphelpjob *next; // (guarded by the pool mutex)
	int posted;              // work/v8s: dphelp_post put the job on the list (dphelp_finish must then take the mutex)
	double prio;             // work/v8s: the posting step's priority (dphelp_myprio)
} dphelpjob_t;

typedef struct
{
	pthread_mutex_t *mutex;  // the progressive pass's global mutex
	pthread_cond_t *cond;    // its condition variable (idle threads wait on it)
	int nidle;               // threads waiting on cond (guarded by mutex)
	dphelpjob_t *jobs;       // posted jobs with free slots (guarded by mutex)
	volatile int *pending;   // optional: > 0 while a guide-tree step is ready to run (read racily)
	volatile double readyprio; // work/v8s: priority of the best ready guide-tree step, -1 if none (written under mutex, read racily)
} dphelppool_t;

extern dphelppool_t *dphelppool;

/* poster */
int dphelp_idle( void );                          // idle threads right now (racy hint, 0 without a pool)
void dphelp_post( dphelpjob_t *job, int maxslots ); // offer up to maxslots helper slots
void dphelp_finish( dphelpjob_t *job );           // withdraw the job and wait until no helper is inside run()
/* idle worker, called with pool->mutex held; returns 1 if it helped (mutex held again on return) */
int dphelp_serve_locked( void );

/* work/v8s: critical helping.  A step on the critical path (its priority above that of every ready
   step) may post a large job even when no thread is idle; a thread that finishes a step serves such a
   job before it takes a ready step of lower priority (dphelp_serve_crit_locked, from the scheduler
   loop with the mutex held).  dphelp_myprio is the priority of the step the calling thread runs. */
extern TLS double dphelp_myprio;
int dphelp_crit_on( void );
int dphelp_iscritical( void );                       // racy hint: this thread's step beats every ready step
void dphelp_post_crit( dphelpjob_t *job, int maxslots ); // as dphelp_post, slots not limited to idle threads
int dphelp_serve_crit_locked( double readyprio );   // serve a posted job with prio > readyprio; 1 if helped

/* body( ctx, lo, hi ) for [0,n) in chunks, by the caller and any idle helpers; runs body( ctx, 0, n )
   directly when nobody is idle or work < minwork (work = n * cost per item, any unit).  The
   chunks must be independent: each one writes only its own items. */
void dphelp_for( int n, int chunk, double work, double minwork, void (*body)( void *ctx, int lo, int hi ), void *ctx );

#endif
#endif
