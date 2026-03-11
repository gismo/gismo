// Detects whether the compiler enables AVX for this build configuration.
// Used by CMake try_run() to decide if Cling's JIT needs EIGEN_DONT_VECTORIZE.
#include <cstdio>
int main() {
#ifdef __AVX__
    printf("AVX");
#else
    printf("NO_AVX");
#endif
    return 0;
}
