#ifndef COMMON_H
#define COMMON_H


// Common declarations
// Math, constants, auxilary functions declaration put here
// Actual functions put into common.cpp

#ifdef USE_DOUBLES
typedef double my_float;
#else
typedef float my_float;
#endif

#ifndef PI
#define PI 3.141592653589793238462643383279
#endif //PI

// development time thingies

//#define DEVTALK
//#define SYNC_PROF
//#define DETECTOR_PARALLEL


#endif

