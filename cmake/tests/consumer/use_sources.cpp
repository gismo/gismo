#include <cstdio>
#include <gsdeptest/dep_lib.h>

// Exercises the MODE SOURCES target gismo_add_dependency() produces: a real
// link+run against the compiled dependency, not just a compile check.
int main()
{
  const int v = gsdeptest_lib_answer();
  if (v == 4242)
  {
    std::printf("OK sources 4242\n");
    return 0;
  }
  std::printf("FAIL sources %d\n", v);
  return 1;
}
