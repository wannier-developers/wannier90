
#include <stdlib.h>
#include <stdio.h>

void f90wrap_abort_()
{
  fprintf(stderr, "The program wants to abort\n");
  exit(0);
}
