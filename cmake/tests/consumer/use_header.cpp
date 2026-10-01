#include <cstdio>
#include <gsdeptest/dep_header.h>

// Exercises the HEADER_ONLY target gismo_add_dependency() produces: a real
// link+run, not just a compile check, since main()'s exit code depends on
// the value the header actually returns.
int main()
{
  const int v = gsdeptest_header_answer();
  if (v == 42)
  {
    std::printf("OK header 42\n");
    return 0;
  }
  std::printf("FAIL header %d\n", v);
  return 1;
}
