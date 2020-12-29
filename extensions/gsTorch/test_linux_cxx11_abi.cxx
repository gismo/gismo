#include <type_traits>

int main()
{ 
  static_assert(_GLIBCXX_USE_CXX11_ABI, "Pre-cxx11 ABI");
  return 0;
}
