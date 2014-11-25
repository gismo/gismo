 
#pragma once

#include <iostream>

namespace gismo
{

  /** 
      Functions to compute the indices in the quad-tree
  */

  
/// lvl: current level, i: index in level lvl, lmax: max level
/// return the global index to be used in the quadtree
inline unsigned Qlocal2global( unsigned const & i, unsigned const & lvl,unsigned const & lmax)
{
    //Returns  i * 2^(lmax-lvl);
    return i << (lmax-lvl);

}

/// lvl: current level, I: global index, lmax: max level
/// return the local index of level lvl to be used in the char-matrix
inline unsigned Qglobal2local( unsigned const & I, unsigned const & lvl,unsigned const & lmax)
{
    // Returns I / 2^(lmax-lvl)
    return I >> (lmax-lvl);
}

/// lvl: some level, lmax: max level. Returns the step (stride) needed
/// to be added to a global index in order to to find the global index
/// of the next unique index of level lvl
inline unsigned Qstride(unsigned const & lvl,unsigned const & lmax)
{
    // Returns 2^(lmax-lvl)
    return 1 << (lmax-lvl);
}

inline unsigned pow2(unsigned const & i)
{
    // Returns 2^i
    return 1 << i;
}


/// lvl1: current level, i: index at current level, lvl2 new level
/// Translate the (unique) index i from  lvl1 to lvl2
inline unsigned Qtranslate( unsigned const & i, unsigned const & lvl1,unsigned const & lvl2)
{
    if ( lvl1< lvl2 )
        return i << (lvl2-lvl1);
    else // lvl1 >= lvl2
      {
	if ( lvl1< lvl2 )
	  return i >> (lvl1-lvl2);
	else
	  return i;
      }
}

/// lvl1: some level, i: index at that level, lvl2 new level
/// number of New knots less than i between lvl1 and lvl2
inline unsigned QnewBefore( unsigned const & i, unsigned const & lvl1, unsigned const & lvl2)
{
    if ( lvl1< lvl2 )
        // i * 2^lvl1 (2^{lvl2-lvl1}-1) = i* (2^lvl2 - 2^lvl1) 
        return 0;
    else
    //  i* (2^lvl2 - 2^lvl1)  / 2^{lvl2-lvl12^{lvl2-lvl1}-1}
        return 0;
}


}; // namespace gismo
