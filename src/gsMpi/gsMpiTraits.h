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

  /// Forward declaration of gsVector class
  template<class T, int _Rows, int _Options> class gsVector;

  /// Specialization for fixed-size gsVector class
  template<class T, int _Rows, int _Options>
  struct MPITraits<gsVector<T, _Rows, _Options> >
  {
  public:
    static inline MPI_Datatype getType()
    {
      if(datatype==MPI_DATATYPE_NULL)
        {
          MPI_Type_contiguous(_Rows, MPITraits<T>::getType(), &vectortype);
          MPI_Type_commit(&vectortype);
          gsVector<T, _Rows, _Options> fvector;
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

  private:
    static MPI_Datatype datatype;
    static MPI_Datatype vectortype;
  };

  template<class T, int _Rows, int _Options>
  MPI_Datatype MPITraits<gsVector<T, _Rows, _Options> >::datatype = MPI_DATATYPE_NULL;
  template<class T, int _Rows, int _Options>
  MPI_Datatype MPITraits<gsVector<T, _Rows, _Options> >::vectortype = {MPI_DATATYPE_NULL};



  /// Forward declaration of gsMatrix class
  template<class T, int _Rows, int _Cols, int _Options> class gsMatrix;

  /// Specialization for fixed-size gsMatrix class
  template<class T, int _Rows, int _Cols, int _Options>
  struct MPITraits<gsMatrix<T, _Rows, _Cols, _Options> >
  {
  public:
    static inline MPI_Datatype getType()
    {
      if(datatype==MPI_DATATYPE_NULL)
        {
          MPI_Type_contiguous(_Rows*_Cols, MPITraits<T>::getType(), &matrixtype);
          MPI_Type_commit(&matrixtype);
          gsMatrix<T, _Rows, _Cols, _Options> fmatrix;
          MPI_Aint base;
          MPI_Aint displ;
          MPI_Get_address(&fmatrix, &base);
          MPI_Get_address(&(fmatrix(0,0)), &displ);
          displ -= base;
          int length[1]={1};

          MPI_Type_create_struct(1, length, &displ, &matrixtype, &datatype);
          MPI_Type_commit(&datatype);
        }
      return datatype;
    }

  private:
    static MPI_Datatype datatype;
    static MPI_Datatype matrixtype;
  };

  template<class T, int _Rows, int _Cols, int _Options>
  MPI_Datatype MPITraits<gsMatrix<T, _Rows, _Cols, _Options> >::datatype = MPI_DATATYPE_NULL;
  template<class T, int _Rows, int _Cols, int _Options>
  MPI_Datatype MPITraits<gsMatrix<T, _Rows, _Cols, _Options> >::matrixtype = {MPI_DATATYPE_NULL};


  /// Forward declaration gsSparseMatrix
  template<typename T, int _Options, typename _Index> class gsSparseMatrix;

  template<typename T, int _Options, typename _Index>
  struct MPITraits<gsSparseMatrix<T, _Options, _Index> >
  {
  public:
    static inline MPI_Datatype getType()
    {
      if (datatype==MPI_DATATYPE_NULL)
        {
          gsSparseMatrix<T, _Options, _Index> mat;
          MPI_Aint base;
          MPI_Get_address(mat.data, &base);
          //MPI_Aint* displ = new int[n];
          //for (int i=0; i<n; ++i)
          //  {
          //    MPI_Get_address(mat[i], &displ[i]);
          //    displ[i] -= base;
          //  }

          //MPI_Type_hindexed(n, m, displ, MPITraits<T>::getType(), &datatype);
          //MPI_Type_commit(&datatype);
        }
      return datatype;
    }

  private:
    static MPI_Datatype datatype;
    static MPI_Datatype matrixtype;
  };

  template<class T, int _Options, typename _Index>
  MPI_Datatype MPITraits<gsSparseMatrix<T, _Options, _Index> >::datatype = MPI_DATATYPE_NULL;
  template<class T, int _Options, typename _Index>
  MPI_Datatype MPITraits<gsSparseMatrix<T, _Options, _Index> >::matrixtype = {MPI_DATATYPE_NULL};

  /*
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
  struct MPITraits<std::pair<T1,T2> >
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
