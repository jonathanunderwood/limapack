#include <stdlib.h>

typedef struct _foo
{
  int x;
  void (*dtor) (struct _foo *self);
} foo_t;

void
foo_dtor(foo_t * foo)
{
  free(foo);
}

int
main ()
{
  foo_t *a = malloc(sizeof(foo_t));
  a->x = 1.0;
  a->dtor = foo_dtor;
  a->dtor(a);

  return 0;
}
