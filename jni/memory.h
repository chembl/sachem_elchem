#ifndef MEMORY_H_
#define MEMORY_H_

#include <string.h>


static inline size_t align_size(size_t size)
{
    return (size + 3) & ~((size_t) 3u);
}


static inline void *alloc_memory(void **memory, size_t size)
{
    size = align_size(size);
    void *address = *memory;
    *memory += size;
    return address;
}


static inline void *alloc_memory_zero(void **memory, size_t size)
{
    size = align_size(size);
    void *address = *memory;
    bzero(address, size);
    *memory += size;
    return address;
}


static inline void *alloc_memory_one(void **memory, size_t size)
{
    size = align_size(size);
    void *address = *memory;
    memset(address, -1, size);
    *memory += size;
    return address;
}

#endif /* MEMORY_H_ */
