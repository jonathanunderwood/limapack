#ifndef __MEMORY_H__
#define __MEMORY_H__


#ifndef ATTRIBUTE_RETURN_CHECK
#if __GNUC_PREREQ (3, 4)
#define ATTRIBUTE_RETURN_CHECK __attribute__((__warn_unused_result__))
#else
#define ATTRIBUTE_RETURN_CHECK
#endif
#endif

/* Don't call these directly - use the macros below */
int memory_alloc(void *ptrptr, size_t size) ATTRIBUTE_RETURN_CHECK;
int memory_allocN(void *ptrptr, size_t size, size_t count) ATTRIBUTE_RETURN_CHECK;
int memory_reallocN(void *ptrptr, size_t size, size_t count) ATTRIBUTE_RETURN_CHECK;
void memory_free(void *ptrptr);
void memory_oomerr (const char *file, const int line);

/**
 * MEMORY_ALLOC:
 * @ptr: pointer to hold address of allocated memory
 *
 * _allocate sizeof(*ptr) bytes of memory and store
 * the address of allocated memory in 'ptr'. Fill the
 * newly allocated memory with zeros.
 *
 * Returns -1 on failure, 0 on success
 */
#define MEMORY_ALLOC(ptr) memory_alloc(&(ptr), sizeof(*(ptr)))

/**
 * MEMORY_ALLOC_N:
 * @ptr: pointer to hold address of allocated memory
 * @count: number of elements to allocate
 *
 * _allocate an array of 'count' elements, each sizeof(*ptr)
 * bytes long and store the address of allocated memory in
 * 'ptr'. Fill the newly allocated memory with zeros.
 *
 * Returns -1 on failure, 0 on success
 */
#define MEMORY_ALLOC_N(ptr, count) memory_allocN(&(ptr), sizeof(*(ptr)), (count))

/**
 * MEMORY_REALLOC_N:
 * @ptr: pointer to hold address of allocated memory
 * @count: number of elements to allocate
 *
 * Re-allocate an array of 'count' elements, each sizeof(*ptr)
 * bytes long and store the address of allocated memory in
 * 'ptr'. Fill the newly allocated memory with zeros
 *
 * Returns -1 on failure, 0 on success
 */
#define MEMORY_REALLOC_N(ptr, count) memory_reallocN(&(ptr), sizeof(*(ptr)), (count))

/**
 * MEMORY_FREE:
 * @ptr: pointer holding address to be freed
 *
 * _free the memory stored in 'ptr' and update to point
 * to NULL.
 */
#define MEMORY_FREE(ptr) memory_free(&(ptr))

/**
 * memory_oomerr:
 *
 * Prints an OOM error message to stderr indicating at which line in 
 * which file the error occured,
 */
#define MEMORY_OOMERR memory_oomerr(__FILE__, __LINE__)

#endif /* __MEMORY_H__ */
