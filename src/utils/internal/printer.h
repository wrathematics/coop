#ifndef __COOP_NAIVE_PRINTER_H__
#define __COOP_NAIVE_PRINTER_H__

static inline void matprinter(const int m, const int n, const double *const restrict x)
{
  int i, j;
  
  for (i=0; i<m; i++)
  {
    for (j=0; j<n; j++)
      printf("%.0f ", x[i+m*j]);
    
    putchar('\n');
  }
  
  putchar('\n');
}

#endif
