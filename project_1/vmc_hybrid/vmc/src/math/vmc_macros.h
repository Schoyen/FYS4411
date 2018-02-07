#ifndef VMC_MACROS_H
#define VMC_MACROS_H

#define SQUARE(x) (x*x)

/* Create a macro for creating a uniform double in the interval [0, 1) */
#define RANDOM_UNIFORM_DOUBLE \
    ((((((unsigned long) arc4random()) << 32) | arc4random()) \
     / ((double) UINT64_MAX)))

#endif
