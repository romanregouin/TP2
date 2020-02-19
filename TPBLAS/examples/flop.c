#include <stdio.h>
#include <x86intrin.h>

// La frequence du processeur est de 2.6 GHZ
static const float duree_cycle = (float) 1 / (float) 2.6 ;
// duree du cycle en nano seconde 10^-9

static unsigned long long int residu ;

void init_flop ()
{
  unsigned long long int start, end ;

  start = _rdtsc () ;

  end =_rdtsc () ;

  residu = end - start ;
  
}


void calcul_flop (char *message, int nb_operations_flottantes, unsigned long long int cycles)
{
  printf ("%s %d operations %5.3f GFLOP/s\n", message, nb_operations_flottantes, ((float) nb_operations_flottantes) / (((float) (cycles - residu)) * duree_cycle)) ;
  return ;
}
