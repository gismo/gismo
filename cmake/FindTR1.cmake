# Check availability of C++ TR1 contents.
#
# Sets the following variables:
#
# TR1_SHARED_PTR_FOUND          -- std::tr1::shared_ptr1<T> available
# TR1_SHARED_PTR_USE_TR1_MEMORY -- #include <tr1/memory>
# TR1_SHARED_PTR_USE_MEMORY     -- #include <memory>

# We need to have at least this version to support the VERSION_LESS argument to 'if' (2.6.2) and unset (2.6.3)
cmake_policy(PUSH)
cmake_minimum_required(VERSION 2.6.3)
cmake_policy(POP)

INCLUDE (CheckCXXSourceCompiles)

# ---------------------------------------------------------------------------
# std::tr1::shared_ptr<T>
# ---------------------------------------------------------------------------

check_cxx_source_compiles(
    "
        #include <memory>
        int main() {
            std::shared_ptr<int> ptr;
            return 0;
        }
    "
    STD_SHARED_PTR_FOUND)

if (STD_SHARED_PTR_FOUND)
   mark_as_advanced (STD_SHARED_PTR_FOUND)
   return()
endif()

check_cxx_source_compiles(
    "
        #include <memory>
        int main() {
            std::tr1::shared_ptr<int> ptr;
            return 0;
        }
    "
    TR1_SHARED_PTR_FOUND)

if (TR1_SHARED_PTR_FOUND)
   mark_as_advanced (TR1_SHARED_PTR_FOUND)
   return()
endif()

check_cxx_source_compiles(
    "
        #include <tr1/memory>
        int main() {
            std::tr1::shared_ptr<int> ptr;
            return 0;
        }
    "
    TR1_SHARED_PTR_USE_TR1_MEMORY)

if (TR1_SHARED_PTR_USE_TR1_MEMORY)
   set (TR1_SHARED_PTR_FOUND TRUE)
   mark_as_advanced (TR1_SHARED_PTR_FOUND) 
   mark_as_advanced (TR1_SHARED_PTR_USE_TR1_MEMORY)
   return()
endif()


check_cxx_source_compiles(
    "
        #include <boost/shared_ptr.hpp>
        int main() {
            boost::shared_ptr<int> ptr;
            return 0;
        }
    "
    BOOST_SHARED_PTR_FOUND)

if (BOOST_SHARED_PTR_FOUND)
   mark_as_advanced (BOOST_SHARED_PTR_FOUND)
   return()
endif()


# ---------------------------------------------------------------------------
# std::tr1::unordered_map<K, V>
# ---------------------------------------------------------------------------

#check_cxx_source_compiles(
#    "
#        #include <tr1/unordered_map>
#        int main() {
#            std::tr1::unordered_map<int, int> m;
#            return 0;
#        }
#    "
#    TR1_UNORDERED_MAP_USE_TR1_UNORDERED_MAP)
#check_cxx_source_compiles(
#    "
#        #include <unordered_map>
#        int main() {
#            std::tr1::unordered_map<int, int> m;
#            return 0;
#        }
#    "
#    TR1_UNORDERED_MAP_USE_UNORDERED_MAP)
#
#set (TR1_UNORDERED_MAP -NOTFOUND)
#if (TR1_UNORDERED_MAP_USE_TR1_UNORDERED_MAP)
#  set (TR1_UNORDERED_MAP_FOUND TRUE)
#endif (TR1_UNORDERED_MAP_USE_TR1_UNORDERED_MAP)
#if (TR1_UNORDERED_MAP_USE_UNORDERED_MAP)
#  set (TR1_UNORDERED_MAP_FOUND TRUE)
#endif (TR1_UNORDERED_MAP_USE_UNORDERED_MAP)
#
#mark_as_advanced (TR1_UNORDERED_MAP_FOUND)
#mark_as_advanced (TR1_UNORDERED_MAP_USE_TR1_UNORDERED_MAP)
#mark_as_advanced (TR1_UNORDERED_MAP_USE_UNORDERED_MAP)

