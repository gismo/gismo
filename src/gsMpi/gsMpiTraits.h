/**
 * @file gsMpiTraits.h
 *
 * @brief Traits classes for mapping types onto MPI_Datatype.
 * @author C.Hofer, taken from DUNE - author Markus Blatt
 */

#pragma once

#include <utility>
#include <cstddef>

namespace gismo
{

  /**
   * @brief A traits class describing the mapping of types onto MPI_Datatypes.
   *
   * Specializations exist for the default types.
   * Specializations should provide a static method
   * \code
   * static MPI_Datatype getType();
   * \endcode
   */
  template<typename T>
  struct MPITraits
  {
  private:
    MPITraits(){}
    MPITraits(const MPITraits&){}
    static MPI_Datatype datatype;
    static MPI_Datatype vectortype;
  public:
    static inline MPI_Datatype getType()
    {
      if(datatype==MPI_DATATYPE_NULL) {
        MPI_Type_contiguous(sizeof(T),MPI_BYTE,&datatype);
        MPI_Type_commit(&datatype);
      }
      return datatype;
    }

  };
  template<class T>
  MPI_Datatype MPITraits<T>::datatype = MPI_DATATYPE_NULL;

#ifndef DOXYGEN

  // A Macro for defining traits for the primitive data types
#define ComposeMPITraits(p,m) \
  template<> \
  struct MPITraits<p>{ \
    static inline MPI_Datatype getType(){ \
      return m; \
    } \
  }

  // Identify inbuild MPI types with the C types
  ComposeMPITraits(char, MPI_CHAR);
  ComposeMPITraits(unsigned char,MPI_UNSIGNED_CHAR);
  ComposeMPITraits(short,MPI_SHORT);
  ComposeMPITraits(unsigned short,MPI_UNSIGNED_SHORT);
  ComposeMPITraits(int,MPI_INT);
  ComposeMPITraits(unsigned int,MPI_UNSIGNED);
  ComposeMPITraits(long,MPI_LONG);
  ComposeMPITraits(unsigned long,MPI_UNSIGNED_LONG);
  ComposeMPITraits(float,MPI_FLOAT);
  ComposeMPITraits(double,MPI_DOUBLE);
  ComposeMPITraits(long double,MPI_LONG_DOUBLE);


#undef ComposeMPITraits

 // FiedVector and bigunsignedint<k> not known is gismo, additionally it needs #include <cstdint>
/*
  template<class T, int n> class gsVector<T,n>; // = gsVector<K,n>

  template<class T, int n>
  struct MPITraits<gsVector<T,n> >
  {
    static MPI_Datatype datatype;
    static MPI_Datatype vectortype;

    static inline MPI_Datatype getType()
    {
      if(datatype==MPI_DATATYPE_NULL) {
        MPI_Type_contiguous(n, MPITraits<T>::getType(), &vectortype);
        MPI_Type_commit(&vectortype);
        gsVector<T,n> fvector;
        MPI_Aint base;
        MPI_Aint displ;
        MPI_Get_address(&fvector, &base);
        MPI_Get_address(&(fvector[0]), &displ);
        displ -= base;
        int length[1]={1};

        MPI_Type_create_struct(1, length, &displ, &vectortype, &datatype);
        MPI_Type_commit(&datatype);
      }
      return datatype;
    }

  };

  template<class T, int n>
  MPI_Datatype MPITraits<gsVector<T,n> >::datatype = MPI_DATATYPE_NULL;
  template<class T, int n>
  MPI_Datatype MPITraits<gsVector<T,n> >::vectortype = {MPI_DATATYPE_NULL};

  template<class T> class gsSparseMatrix<T>;

  template<class T, int n, int m>
  struct MPITraits<gsSparseMatrix<T> >
  {
    static MPI_Datatype datatype;
    static MPI_Datatype vectortype;

    static inline MPI_Datatype getType()
    {
      if(datatype==MPI_DATATYPE_NULL) {

          gsSparseMatrix<T,n,m> mat;
          MPI_Aint base;
          MPI_Address(mat.data), &base);
          MPI_Aint* displ = new int[n];
          for (int i=0; i<n; ++i)
          {
              MPI_Address(mat[i], &displ[i]);
              displ[i] -= base;
          }


          MPI_Type_hindexed(n, m, displ, MPITraits<T>::getType(), &datatype);
          MPI_Type_commit(&datatype);
      }
      return datatype;
    }

  };

  template<class T, int n>
  MPI_Datatype MPITraits<gsVector<T,n> >::datatype = MPI_DATATYPE_NULL;
  template<class T, int n>
  MPI_Datatype MPITraits<gsVector<T,n> >::vectortype = {MPI_DATATYPE_NULL};

  template<int k>
  class bigunsignedint;

  template<int k>
  struct MPITraits<bigunsignedint<k> >
  {
    static MPI_Datatype datatype;
    static MPI_Datatype vectortype;

    static inline MPI_Datatype getType()
    {
      if(datatype==MPI_DATATYPE_NULL) {
        MPI_Type_contiguous(bigunsignedint<k>::n, MPITraits<std::uint16_t>::getType(),
                            &vectortype);
        //MPI_Type_commit(&vectortype);
        bigunsignedint<k> data;
        MPI_Aint base;
        MPI_Aint displ;
        MPI_Get_address(&data, &base);
        MPI_Get_address(&(data.digit), &displ);
        displ -= base;
        int length[1]={1};
        MPI_Type_create_struct(1, length, &displ, &vectortype, &datatype);
        MPI_Type_commit(&datatype);
      }
      return datatype;
    }
  };

  template<int k>
  MPI_Datatype MPITraits<bigunsignedint<k> >::datatype = MPI_DATATYPE_NULL;
  template<int k>
  MPI_Datatype MPITraits<bigunsignedint<k> >::vectortype = MPI_DATATYPE_NULL;
*/

  template<typename T1, typename T2>
  struct MPITraits<std::pair<T1,T2 > >
  {
  public:
    inline static MPI_Datatype getType();
  private:
    static MPI_Datatype type;
  };
  template<typename T1, typename T2>
  MPI_Datatype MPITraits<std::pair<T1,T2> >::getType()
  {
    if(type==MPI_DATATYPE_NULL) 
    {
        int length[2] = {1, 1};
        MPI_Aint disp[2];
        MPI_Datatype types[2] = {MPITraits<T1>::getType(),
                                 MPITraits<T2>::getType()};
        
        typedef std::pair<T1, T2> Pair;
        //static_assert(std::is_standard_layout<Pair>::value, "offsetof() is only defined for standard layout types");
        disp[0] = offsetof(Pair, first);
        disp[1] = offsetof(Pair, second);
        
        MPI_Datatype tmp;
        MPI_Type_create_struct(2, length, disp, types, &tmp);
        
        MPI_Type_create_resized(tmp, 0, sizeof(Pair), &type);
        MPI_Type_commit(&type);
        MPI_Type_free(&tmp);
    }
    return type;
  }

  template<typename T1, typename T2>
  MPI_Datatype MPITraits<std::pair<T1,T2> >::type=MPI_DATATYPE_NULL;
}

#endif
