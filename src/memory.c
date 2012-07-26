#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <errno.h>

#include "memory.h"

/* Return 1 if an array of N objects, each of size S, cannot exist due
   to size arithmetic overflow.  S must be positive and N must be
   nonnegative.  This is a macro, not an inline function, so that it
   works correctly even when SIZE_MAX < N.

   By gnulib convention, SIZE_MAX represents overflow in size
   calculations, so the conservative dividend to use here is
   SIZE_MAX - 1, since SIZE_MAX might represent an overflowed value.
   However, malloc (SIZE_MAX) fails on all known hosts where
   sizeof (ptrdiff_t) <= sizeof (size_t), so do not bother to test for
   exactly-SIZE_MAX allocations on such hosts; this avoids a test and
   branch when S is known to be 1.  */
#ifndef xalloc_oversized
# define xalloc_oversized(n, s) \
    ((size_t) (sizeof (ptrdiff_t) <= sizeof (size_t) ? -1 : -2) / (s) < (n))
#endif

/**
 * memory_alloc:
 * @ptrptr: pointer to pointer for address of allocated memory
 * @size: number of bytes to allocate
 *
 * _allocate  'size' bytes of memory. Return the address of the
 * allocated memory in 'ptrptr'. The newly allocated memory is
 * filled with zeros.
 *
 * Returns -1 on failure to allocate, zero on success
 */
int memory_alloc(void *ptrptr, size_t size)
{
    *(void **)ptrptr = calloc(1, size);
    if (*(void **)ptrptr == NULL)
        return -1;
    return 0;
}

/**
 * memory_allocN:
 * @ptrptr: pointer to pointer for address of allocated memory
 * @size: number of bytes to allocate
 * @count: number of elements to allocate
 *
 * _allocate an array of memory 'count' elements long,
 * each with 'size' bytes. Return the address of the
 * allocated memory in 'ptrptr'.  The newly allocated
 * memory is filled with zeros.
 *
 * Returns -1 on failure to allocate, zero on success
 */
int memory_allocN(void *ptrptr, size_t size, size_t count)
{
    *(void**)ptrptr = calloc(count, size);
    if (*(void**)ptrptr == NULL)
        return -1;
    return 0;
}

/**
 * memory_reallocN:
 * @ptrptr: pointer to pointer for address of allocated memory
 * @size: number of bytes to allocate
 * @count: number of elements in array
 *
 * Resize the block of memory in 'ptrptr' to be an array of
 * 'count' elements, each 'size' bytes in length. Update 'ptrptr'
 * with the address of the newly allocated memory. On failure,
 * 'ptrptr' is not changed and still points to the original memory
 * block. The newly allocated memory is filled with zeros.
 *
 * Returns -1 on failure to allocate, zero on success
 */
int memory_reallocN(void *ptrptr, size_t size, size_t count)
{
    void *tmp;

    if (xalloc_oversized(count, size)) {
        errno = ENOMEM;
        return -1;
    }
    tmp = realloc(*(void**)ptrptr, size * count);
    if (!tmp && (size * count))
        return -1;
    *(void**)ptrptr = tmp;
    return 0;
}

/**
 * memory_free:
 * @ptrptr: pointer to pointer for address of memory to be freed
 *
 * Release the chunk of memory in the pointer pointed to by
 * the 'ptrptr' variable. After release, 'ptrptr' will be
 * updated to point to NULL.
 */
void memory_free(void *ptrptr)
{
    free(*(void**)ptrptr);
    *(void**)ptrptr = NULL;
}

/**
 * memory_oomerr:
 * @file: pointer to a string containing the current source filename. 
 *        Usually, the macro __FILE__ will be used
 * @line: line number of file
 *
 * Prints an OOM error message to stderr.
 */
void
memory_oomerr (const char *file, const int line)
{
  fprintf(stderr, "Failed to allocate memory at line %d in file %s.\n",
	  line, file);
}
