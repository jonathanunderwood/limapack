#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

/* This simple program finds the length that a string needs to be in
   order to hold the largest integer value. We use INT_MIN to ensure
   there's space for the minus sign. snprintf returns the number of
   bits that would have been written in the event that the string
   length wouldn't fit in the buffer passed to it. See:

   http://stackoverflow.com/questions/10536207/ansi-c-maximum-number-of-characters-printing-a-decimal-int
*/

int
main()
{
  char sbuf[2];
  int ndigits;

  ndigits = snprintf(sbuf, (size_t) 1, "%lld", (long long) INT_MIN);

  printf("ndigits: %d\n", ndigits);

  return 0;
}
