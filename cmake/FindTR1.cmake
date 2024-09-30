# Check availability of C++ TR1 contents.
#
# Sets the following variables:
#
# TR1_SHARED_PTR_FOUND          -- std::tr1::shared_ptr1<T> available
# TR1_SHARED_PTR_USE_TR1_MEMORY -- #include <tr1/memory>
# TR1_SHARED_PTR_USE_MEMORY     -- #include <memory>

include(CheckCXXSourceCompiles)

# ---------------------------------------------------------------------------
# std::tr1::shared_ptr<T>
# ---------------------------------------------------------------------------

#check_cxx_source_compiles(
#    "
#        #include <memory>
#        int main() {
#            std::shared_ptr<int> ptr;
#            return 0;
#        }
#    "
#    STD_SHARED_PTR_FOUND)
#if (STD_SHARED_PTR_FOUND)
#   mark_as_advanced (STD_SHARED_PTR_FOUND)
#else()

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
else()

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
else()


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
endif(BOOST_SHARED_PTR_FOUND)
endif(TR1_SHARED_PTR_USE_TR1_MEMORY)
endif(TR1_SHARED_PTR_FOUND)
#endif(STD_SHARED_PTR_FOUND)

# ---------------------------------------------------------------------------
# std::unique_ptr<T>
# ---------------------------------------------------------------------------

#check_cxx_source_compiles(
#    "
#        #include <memory>
#        int main() {
#            std::unique_ptr<int> ptr;
#            return 0;
#        }
#    "
#    STD_UNIQUE_PTR_FOUND)
#if (STD_UNIQUE_PTR_FOUND)
#   mark_as_advanced (STD_UNIQUE_PTR_FOUND)
#else()

# boost::movelib::unique_ptr requires boost libraries >= 1.57.0
# see https://dieboostcppbibliotheken.de/boost.smartpointers (German)
check_cxx_source_compiles(
    "
        #include <boost/move/unique_ptr.hpp>
        int main() {
            boost::movelib::unique_ptr<int> ptr;
            return 0;
        }
    "
    BOOST_UNIQUE_PTR_FOUND)

if (BOOST_UNIQUE_PTR_FOUND)
   mark_as_advanced (BOOST_UNIQUE_PTR_FOUND)
endif(BOOST_UNIQUE_PTR_FOUND)
#endif(TR1_UNIQUE_PTR_USE_TR1_MEMORY)
#endif(TR1_UNIQUE_PTR_FOUND)
#endif(STD_UNIQUE_PTR_FOUND)

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

