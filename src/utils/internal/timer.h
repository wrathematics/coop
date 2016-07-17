#ifndef __COOP_UTILS_INTERNAL_TIMER_H__
#define __COOP_UTILS_INTERNAL_TIMER_H__

#include <time.h>
#include <stdio.h>

clock_t _timer_start, _timer_stop;

#define TIMER(expression,nreps) \
  _timer_start = clock(); \
  for (int i=0; i<nreps; i++) \
    expression; \
  _timer_stop = clock(); \
  printf(#expression ":\t%f\n", (double)(_timer_stop-_timer_start)/CLOCKS_PER_SEC);

#endif
