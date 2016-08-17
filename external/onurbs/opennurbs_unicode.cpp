/* $NoKeywords: $ */
/*
//
// Copyright (c) 1993-2012 Robert McNeel & Associates. All rights reserved.
// OpenNURBS, Rhinoceros, and Rhino3D are registered trademarks of Robert
// McNeel & Associates.
//
// THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY.
// ALL IMPLIED WARRANTIES OF FITNESS FOR ANY PARTICULAR PURPOSE AND OF
// MERCHANTABILITY ARE HEREBY DISCLAIMED.
//				
// For complete openNURBS copyright information see <http://www.opennurbs.org>.
//
////////////////////////////////////////////////////////////////
*/

#include "opennurbs.h"

int ON_IsValidUnicodeCodePoint(ON__UINT32 u)
{
  return ( u < 0xD800 || (u >= 0xE000 && u <= 0x10FFFF) );
}

enum ON_UnicodeEncoding ON_UnicodeNativeCPU_UTF16()
{
  return (ON::little_endian== ON::Endian()) ? ON_UTF_16LE : ON_UTF_16BE;
}

enum ON_UnicodeEncoding ON_UnicodeNativeCPU_UTF32()
{
  return (ON::little_endian== ON::Endian()) ? ON_UTF_32LE : ON_UTF_32BE;
}

enum ON_UnicodeEncoding ON_IsUTFByteOrderMark(
  const void* buffer,
  size_t sizeof_buffer
  )
{  
  if ( 0 != buffer && sizeof_buffer >= 2 )
  {
    const unsigned char* b = static_cast<const unsigned char*>(buffer);

    if ( 0 == b[0] )
    {
      if ( sizeof_buffer >= 4 && 0 == b[1] && 0xFE == b[2] && 0xFF == b[3] )
        return ON_UTF_32BE;
    }
    else if ( 0xEF == b[0] )
    {
      if ( sizeof_buffer >= 3 && 0xBB == b[1] && 0xBF == b[2] )
        return ON_UTF_8;
    }
    else if ( 0xFE == b[0] )
    {
      if ( 0xFF == b[1] )
        return ON_UTF_16BE;
    }
    else if ( 0xFF == b[0] && 0xFE == b[1] )
    {
      return ( sizeof_buffer >= 4 && 0 == b[2] && 0 == b[3] )
        ? ON_UTF_32LE
        : ON_UTF_16LE;
    }

  }

  return ON_UTF_unset;
}

int ON_IsUTF8ByteOrderMark(
  const ON__UINT8* sUTF8,
  int sUTF8_count
  )
{
  if ( 0 == sUTF8 )
    return 0;
  if ( -1 != sUTF8_count || sUTF8_count < 3 )
    return 0;
  return (0xEF == sUTF8[0] && 0xBB == sUTF8[1] && 0xBF == sUTF8[2]);
}

int ON_EncodeUTF8( ON__UINT32 u, ON__UINT8 sUTF8[6] )
{
  ON__UINT32 c;

  if ( u <= 0x7F )
  {
    // 1 byte UTF8 encoding: 0xxxxxxx
    sUTF8[0] = (ON__UINT8)u;
    return 1;
  }

  if ( u <= 0x7FF )
  {
    // 2 byte UTF8 encoding: 110xxxxx, 10xxxxxx
    c = (u / 0x40);  // c  = 000xxxxx
    c |= 0xC0;                      //   |= 11000000
    sUTF8[0] = (ON__UINT8)c;
    c = (u & 0x3F);
    c |= 0x80;
    sUTF8[1] = (ON__UINT8)c;
    return 2;
  }

  if ( u <= 0xFFFF )
  {
    // 3 byte UTF8 encoding: 1110xxxx, 10xxxxxx, 10xxxxxx
    c = (u / 0x1000); // c  = 0000xxxx
    c |= 0xE0;                       //   |= 11100000
    sUTF8[0] = (ON__UINT8)c;
    c = ((u & 0xFFF) / 0x40);
    c |= 0x80;
    sUTF8[1] = (ON__UINT8)c;
    c = u & 0x3F;
    c |= 0x80;
    sUTF8[2] = (ON__UINT8)c;
    return 3;
  }

  if ( u <= 0x1FFFFF )
  {
    // 4 byte UTF8 encoding: 11110xxx, 10xxxxxx, 10xxxxxx, 10xxxxxx
    c = (u / 0x40000);  // c  = 00000xxx
    c |= 0xF0;                         //   |= 11110000
    sUTF8[0] = (ON__UINT8)c;
    c = ((u & 0x3FFFF)/0x1000);
    c |= 0x80;
    sUTF8[1] = (ON__UINT8)c;
    c = ((u & 0xFFF) / 0x40);
    c |= 0x80;
    sUTF8[2] = (ON__UINT8)c;
    c = u & 0x3F;
    c |= 0x80;
    sUTF8[3] = (ON__UINT8)c;
    return 4;
  }

  if ( u <= 0x3FFFFFF )
  {
    // 5 byte UTF8 encoding: 111110xx, 10xxxxxx, 10xxxxxx, 10xxxxxx, 10xxxxxx
    c = (u / 0xFFFFFF); // c  = 000000xx
    c |= 0xF8;                         //   |= 11111000
    sUTF8[0] = (ON__UINT8)c;
    c = ((u & 0xFFFFFF)/0x40000);
    c |= 0x80;
    sUTF8[1] = (ON__UINT8)c;
    c = ((u & 0x3FFFF)/0x1000);
    c |= 0x80;
    sUTF8[2] = (ON__UINT8)c;
    c = ((u & 0xFFF) / 0x40);
    c |= 0x80;
    sUTF8[3] = (ON__UINT8)c;
    c = u & 0x3F;
    c |= 0x80;
    sUTF8[4] = (ON__UINT8)c;
    return 5;
  }

  if ( u <= 0x7FFFFFFF )
  {
    // 6 byte UTF8 encoding: 111110xx, 10xxxxxx, 10xxxxxx, 10xxxxxx, 10xxxxxx, 10xxxxxx
    c = (u / 0x40000000); // c  = 0000000x
    c |= 0xFC;                           //   |= 11111100
    sUTF8[0] = (ON__UINT8)c;
    c = ((u & 0x3FFFFFFF)/0x1000000);
    c |= 0x80;
    sUTF8[1] = (ON__UINT8)c;
    c = ((u & 0xFFFFFF)/0x40000);
    c |= 0x80;
    sUTF8[2] = (ON__UINT8)c;
    c = ((u & 0x3FFFF)/0x1000);
    c |= 0x80;
    sUTF8[3] = (ON__UINT8)c;
    c = ((u & 0xFFF) / 0x40);
    c |= 0x80;
    sUTF8[4] = (ON__UINT8)c;
    c = u & 0x3F;
    c |= 0x80;
    sUTF8[5] = (ON__UINT8)c;
    return 6;
  }

  return 0;
}

static int ON_DecodeUTF8Helper(
    const ON__UINT8* sUTF8,
    int sUTF8_count,
    ON__UINT32* value,
    unsigned int* error_status
    )
{
#define INPUT_BUFFER_TOO_SHORT 16
#define INVALID_CONTINUATION_VALUE 16
#define OVERLONG_ENCODING 8

  ON__UINT32 u;
  ON__UINT8 c;

  c = sUTF8[0];
   
  if ( 0 == (0x80 & c) )
  {
    // 1 byte ASCII encoding: 0xxxxxxx
    *value = c;
    return 1;
  }

  if ( 0xC0 == ( 0xE0 & c) )
  {
    // 2 byte character encoding: 10xxxxxx, 10xxxxxx
    if ( sUTF8_count < 2 )
    {
      *error_status |= INPUT_BUFFER_TOO_SHORT; // input buffer too short
      return 0;
    }
    u = (0x1F & c);
    c = sUTF8[1];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    if ( u <= 0x7F )
    {
      *error_status |= OVERLONG_ENCODING; // overlong 2 byte character encoding
    }
    *value = u;
    return 2;
  }

  if ( 0xE0 == ( 0xF0 & c) )
  {
    // 3 byte character encoding: 110xxxxx, 10xxxxxx, 10xxxxxx
    if ( sUTF8_count < 3 )
    {
      *error_status |= INPUT_BUFFER_TOO_SHORT; // input buffer too short
      return 0;
    }
    u = (0x0F & c);
    c = sUTF8[1];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    c = sUTF8[2];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    if ( u <= 0x7FF )
    {
      *error_status |= OVERLONG_ENCODING; // overlong 3 byte character encoding
    }
    *value = u;
    return 3;
  }

  if ( 0xF0 == ( 0xF8 & c) )
  {
    // 4 byte character encoding: 11110xxx, 10xxxxxx, 10xxxxxx, 10xxxxxx
    if ( sUTF8_count < 4 )
    {
      *error_status |= INPUT_BUFFER_TOO_SHORT; // input buffer too short
      return 0;
    }

    u = (0x07 & c);
    c = sUTF8[1];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    c = sUTF8[2];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    c = sUTF8[3];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    if ( u <= 0xFFFF )
    {
      *error_status |= OVERLONG_ENCODING; // overlong 4 byte character encoding
    }
    *value = u;
    return 4;
  }
  
  if ( 0xF8 == ( 0xFC & c) )
  {
    // 5 byte character encoding: 111110xx, 10xxxxxx, 10xxxxxx, 10xxxxxx, 10xxxxxx
    if ( sUTF8_count < 5 )
    {
      *error_status |= INPUT_BUFFER_TOO_SHORT; // input buffer too short
      return 0;
    }

    u = (0x03 & c);
    c = sUTF8[1];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    c = sUTF8[2];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    c = sUTF8[3];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    c = sUTF8[4];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    if ( u <= 0x1FFFFF )
    {
      *error_status |= OVERLONG_ENCODING; // overlong 5 byte character encoding
    }
    *value = u;
    return 5;
  }

  if ( 0xFC == ( 0xFE & c) )
  {
    // 6 byte character encoding: 110xxxxx, 10xxxxxx, 10xxxxxx, 10xxxxxx, 10xxxxxx, 10xxxxxx
    if ( sUTF8_count < 6 )
    {
      *error_status |= INPUT_BUFFER_TOO_SHORT; // input buffer too short
      return 0;
    }

    u = (0x01 & c);
    c = sUTF8[1];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    c = sUTF8[2];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    c = sUTF8[3];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    c = sUTF8[4];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    c = sUTF8[5];
    if (  0x80 != ( 0xC0 & c) )
    {
      *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 continuation value
      return 0;
    }
    u *= 64;
    u |= (0x3F & c);
    if ( u <= 0x3FFFFFF )
    {
      *error_status |= OVERLONG_ENCODING; // overlong 6 byte character encoding
    }
    *value = u;
    return 6;
  }

  *error_status |= INVALID_CONTINUATION_VALUE; // invalid UTF=8 start value
  return 0;

#undef INPUT_BUFFER_TOO_SHORT
#undef INVALID_CONTINUATION_VALUE
#undef OVERLONG_ENCODING
}

int ON_DecodeUTF8(
    const ON__UINT8* sUTF8,
    int sUTF8_count,
    struct ON_UnicodeErrorParameters* e,
    ON__UINT32* unicode_code_point
    )
{
  ON__UINT32 u0, u1;
  int i0, i1;
  unsigned int error_status;
  ON__UINT16 sUTF16[2];
  ON__UINT8 c;

  if (  0 == sUTF8 || sUTF8_count <= 0 || 0 == unicode_code_point )
    return 0;

  // special cases for most common unicode values
  // If any error conditions exist, then ON_DecodeUTF8Helper()
  // is used.
  if ( 0 == (0x80 & sUTF8[0]) )
  {
    *unicode_code_point = sUTF8[0];
    return 1;
  }
  
  c = sUTF8[0];
  if ( 0xC0 == ( 0xE0 & c) && sUTF8_count >= 2 )
  {
    // 2 byte character encoding: 10xxxxxx, 10xxxxxx
    u0 = (0x1F & c);
    c = sUTF8[1];
    if (  0x80 == ( 0xC0 & c) )
    {
      u0 *= 64;
      u0 |= (0x3F & c);
      if ( u0 > 0x7F )
      {
        *unicode_code_point = u0;
        return 2;
      }
    }
  }
  else if ( 0xE0 == ( 0xF0 & c) && sUTF8_count >= 3 )
  {
    // 3 byte character encoding: 110xxxxx, 10xxxxxx, 10xxxxxx
    u0 = (0x0F & c);
    c = sUTF8[1];
    if (  0x80 == ( 0xC0 & c) )
    {
      u0 *= 64;
      u0 |= (0x3F & c);
      c = sUTF8[2];
      if (  0x80 == ( 0xC0 & c) )
      {
        u0 *= 64;
        u0 |= (0x3F & c);
        if ( u0 >= 0x0800 && (u0 <= 0xD800 || u0 >= 0xE000) )
        {
          *unicode_code_point = u0;
          return 3;
        }
      }
    }
  }
  else if ( 0xF0 == ( 0xF8 & c) && sUTF8_count >= 4 )
  {
    // 4 byte character encoding: 11110xxx, 10xxxxxx, 10xxxxxx, 10xxxxxx
    u0 = (0x07 & c);
    c = sUTF8[1];
    if (  0x80 == ( 0xC0 & c) )
    {
      u0 *= 64;
      u0 |= (0x3F & c);
      c = sUTF8[2];
      if (  0x80 == ( 0xC0 & c) )
      {
        u0 *= 64;
        u0 |= (0x3F & c);
        c = sUTF8[3];
        if (  0x80 == ( 0xC0 & c) )
        {
          u0 *= 64;
          u0 |= (0x3F & c);
          if ( u0 >= 0x010000 && u0 <= 0x10FFFF )
          {
            *unicode_code_point = u0;
            return 4;
          }
        }
      }
    }
  }


  error_status = 0;
  u0 = 0xFFFFFFFF;
  i0 = ON_DecodeUTF8Helper(sUTF8,sUTF8_count,&u0,&error_status);
  if ( i0 > 0 && 0 == error_status && (u0 < 0xD800 || (u0 >= 0xE000 && u0 <= 0x10FFFF) ) )
  {
    // valid UTF-8 multibyte encoding parsed
    *unicode_code_point = u0;
    return i0;
  }

  // handle errors
  if ( 0 == e )
  {
    // no errors are masked.
    return 0;
  }

  // report error condition
  e->m_error_status |= error_status;

  if ( error_status != (error_status & e->m_error_mask) )
  {
    // this error is not masked
    return 0;
  }

  if ( i0 <= 0 )
  {
    i0 = 1;
    if ( ON_IsValidUnicodeCodePoint(e->m_error_code_point) )
    {
      // skip to next UTF-8 start elemement
      for ( /*empty for initializer*/; i0 < sUTF8_count; i0++ )
      {
        // Search for the next element of sUTF8[] that is the 
        // start of a UTF-8 encoding sequence.
        c = sUTF8[i0];
        if (    0 == (0x80 & c)     // ASCII 0 - 127
              || 0xC0 == ( 0xE0 & c) // 2 byte encoding first character
              || 0xE0 == ( 0xF0 & c) // 3 byte encoding first character
              || 0xF0 == ( 0xF8 & c) // 4 byte encoding first character
              || 0xF8 == ( 0xFC & c) // 5 byte encoding first character
              || 0xFC == ( 0xFE & c) // 6 byte encoding first character
            )
        {
          // resume parsing at this character
          break;
        }
      }
      *unicode_code_point = e->m_error_code_point;
    }
    return i0;
  }

  if ( ON_IsValidUnicodeCodePoint(u0) && 8 == error_status )
  {
    // overlong UTF-8 multibyte encoding of valid unicode code point
    *unicode_code_point = u0;
    return i0;
  }
  
  if ( i0 < sUTF8_count 
       && u0 >= 0xD800 && u0 <= 0xDBFF 
       && (0 == error_status || 8 == error_status)
       && 0 != (4 & e->m_error_mask) 
     )
  {
    // See if a UFT-16 surrogate pair was incorrectly encoded 
    // as two consecutive UTF-8 sequences.
    u1 = 0xFFFFFFFF;
    i1 = ON_DecodeUTF8Helper(sUTF8+i0,sUTF8_count-i0,&u1,&error_status);
    if ( i1 > 0 && (0 == error_status || 8 == error_status) )
    {
      error_status = 0;
      sUTF16[0] = (ON__UINT16)u0;
      sUTF16[1] = (ON__UINT16)u1;
      u0 = 0xFFFFFFFF;
      if ( 2 == ON_ConvertUTF16ToUTF32(false,sUTF16,2,&u0,1,&error_status,0,0,0) 
           && 0 == error_status 
           && ON_IsValidUnicodeCodePoint(u0)
         )
      {
        *unicode_code_point = u0;
        e->m_error_status |= 4;
        return i0+i1;
      }
    }
  }

  if ( ON_IsValidUnicodeCodePoint(e->m_error_code_point) )
  {
    *unicode_code_point = e->m_error_code_point;
    return i0;
  }

  return 0;
}

int ON_EncodeUTF16( ON__UINT32 unicode_code_point, ON__UINT16 sUTF16[2] )
{
  // put the most comman case first
  if ( unicode_code_point < 0xD800 )
  {
    // code point values U+0000 ... U+D7FF
    // = UTF-16 values
    sUTF16[0] = (ON__UINT16)unicode_code_point;
    sUTF16[1] = 0;
    return 1;
  }

  if ( unicode_code_point < 0xE000 )
  {
    // 0xD800 ... 0xDFFF are invalid unicode code point values
    sUTF16[0] = sUTF16[1] = 0;
    return 0;
  }

  if ( unicode_code_point <= 0xFFFF )
  {
    // code point values U+E000 ... U+FFFF
    // = UTF-16 values
    sUTF16[0] = (ON__UINT16)unicode_code_point;
    sUTF16[1] = 0;
    return 1;
  }

  if ( unicode_code_point <= 0x10FFFF )
  {
    // code point values U+10000 ... U+10FFFF
    // = surrogate pair UTF-16 values
    unicode_code_point -= 0x10000;
    sUTF16[0] = (ON__UINT16)(0xD800 + (unicode_code_point / 0x400)); // high surrogate value (0xD800 ... 0xDBFF)
    sUTF16[1] = (ON__UINT16)(0xDC00 + (unicode_code_point & 0x3FF)); // low surrogate value (0xDC00 ... 0xDFFF)
    return 2;
  }

  // 0x110000 ... 0xFFFFFFFF are invalid unicode code point values
  sUTF16[0] = sUTF16[1] = 0;
  return 0;
}

int ON_DecodeUTF16(
    const ON__UINT16* sUTF16,
    int sUTF16_count,
    struct ON_UnicodeErrorParameters* e,
    ON__UINT32* unicode_code_point
    )
{
  ON__UINT32 uhi, ulo;

  if ( 0 == sUTF16 || sUTF16_count <= 0 || 0 == unicode_code_point )
  {
    if ( e )
      e->m_error_status |= 1;
    return 0;
  }

  // special case for most common UTF-16 single element values
  if ( ( sUTF16[0] < 0xD800 ) || ( sUTF16[0] >= 0xE000 ) )
  {
    *unicode_code_point = sUTF16[0];
    return 1;
  }

  if ( sUTF16_count >= 2 && sUTF16[0] < 0xDC00 && sUTF16[1] >=  0xDC00 && sUTF16[1] < 0xE000 )
  {
    // UTF-16 surrogate pair
    uhi = sUTF16[0];
    ulo = sUTF16[1];
    *unicode_code_point = (uhi-0xD800)*0x400 + (ulo-0xDC00) + 0x10000;
    return 2;
  }

  // handle errors
  if ( 0 == e )
  {
    // no errors are masked.
    return 0;
  }

  // report error condition
  e->m_error_status |= 16;

  if ( 16 != (16 & e->m_error_mask) || !ON_IsValidUnicodeCodePoint(e->m_error_code_point) )
  {
    // this error is not masked
    return 0;
  }

  // Search for the next element of sUTF16[] that is a 
  // valid UTF-16 encoding sequence.
  int i;
  for ( i = 1; i < sUTF16_count; i++ )
  {
    if ( ( sUTF16[i] < 0xD800 ) || ( sUTF16[i] >= 0xE000 ) )
    {
      // valid single UTF-16 code unit
      break;
    }
    if ( i+1 < sUTF16_count 
         && sUTF16[i] >= 0xD800 && sUTF16[i] < 0xDC00 
         && sUTF16[i+1] >= 0xDC00 && sUTF16[i+1] < 0xE000
       )
    {
      // valid UTF-16 surrogate pair
      break;
    }
  }

  *unicode_code_point = e->m_error_code_point;

  return i;
}

int ON_DecodeSwapByteUTF16(
    const ON__UINT16* sUTF16,
    int sUTF16_count,
    struct ON_UnicodeErrorParameters* e,
    ON__UINT32* unicode_code_point
    )
{
  int i;
  ON__UINT32 uhi, ulo;
  ON__UINT16 w0, w1;
  const ON__UINT8* p;
  ON__UINT8* p0;
  ON__UINT8* p1;


  if ( 0 == sUTF16 || sUTF16_count <= 0 || 0 == unicode_code_point )
  {
    if ( e )
      e->m_error_status |= 1;
    return 0;
  }

  // special case for most common UTF-16 single element values
  // w0 = byte swapped sUTF16[0]
  p = (const ON__UINT8*)sUTF16;
  p0 = (ON__UINT8*)&w0;
  p0[1] = p[0];
  p0[0] = p[1];
  if ( ( w0 < 0xD800 ) || (w0 >= 0xE000 ) )
  {
    *unicode_code_point = w0;
    return 1;
  }

  if ( sUTF16_count >= 2 && w0 < 0xDC00 )
  {
    // w1 = byte swapped sUTF16[1]
    p1 = (ON__UINT8*)&w1;
    p1[1] = p[2];
    p1[0] = p[3];
    if ( w1 >=  0xDC00 && w1 < 0xE000 )
    {
      // UTF-16 surrogate pair
      uhi = w0;
      ulo = w1;
      *unicode_code_point = (uhi-0xD800)*0x400 + (ulo-0xDC00) + 0x10000;
      return 2;
    }
  }

  // handle errors
  if ( 0 == e )
  {
    // no errors are masked.
    return 0;
  }

  // report error condition
  e->m_error_status |= 16;

  if ( 16 != (16 & e->m_error_mask) || !ON_IsValidUnicodeCodePoint(e->m_error_code_point) )
  {
    // this error is not masked
    return 0;
  }

  // Search for the next element of sUTF16[] that is a 
  // valid UTF-16 encoding sequence.
  p1 = (ON__UINT8*)&w1;
  p += sizeof(sUTF16[0]);
  for ( i = 1; i < sUTF16_count; i++, p += sizeof(sUTF16[0]) )
  {
    // w0 = byte swapped sUTF16[i]
    p0[1] = p[0];
    p0[0] = p[1];
    if ( ( w0 < 0xD800 ) || ( w0 >= 0xE000 ) )
    {
      // valid single UTF-16 code unit
      break;
    }
    if ( i+1 < sUTF16_count && w0 >= 0xD800 && w0 < 0xDC00 )
    {
      // w1 = byte swapped sUTF16[i+1]
      p1[1] = p[sizeof(sUTF16[0])];
      p1[0] = p[sizeof(sUTF16[0])+1];
      if ( w1 >= 0xDC00 && w1 < 0xE000 )
      {
        // valid UTF-16 surrogate pair
        break;
      }
    }
  }

  *unicode_code_point = e->m_error_code_point;

  return i;
}

int ON_ConvertUTF8ToUTF16(
    int bTestByteOrder,
     const ON__UINT8* sUTF8,
    int sUTF8_count,
    ON__UINT16* sUTF16,
    int sUTF16_count,
    unsigned int* error_status,
    unsigned int error_mask,
    ON__UINT32 error_code_point,
    const ON__UINT8** sNextUTF8
    )
{
  int i, j, k, output_count;
  ON__UINT32 u = 0;
  ON__UINT16 w[2];
  struct ON_UnicodeErrorParameters e;

  if ( 0 != error_status )
    *error_status = 0;

  if ( -1 == sUTF8_count && 0 != sUTF8 )
  {
    for ( sUTF8_count = 0; 0 != sUTF8[sUTF8_count]; sUTF8_count++)
    {
      // empty for body
    }
  }

  if ( 0 == sUTF8 || sUTF8_count < 0 )
  {
    if ( 0 != error_status )
      *error_status |= 1;
    if ( sNextUTF8 )
      *sNextUTF8 = sUTF8;
    return 0;
  }

  if ( bTestByteOrder && ON_IsUTF8ByteOrderMark(sUTF8,sUTF8_count) )
  {
    // skip UTF-8 byte order element
    sUTF8_count -= 3;
    sUTF8 += 3;
  }

  if ( 0 == sUTF16_count )
  {
    sUTF16 = 0;
    sUTF16_count = 2147483647; // maximum value of a 32-bit signed int
  }
  else if ( 0 == sUTF16 )
  {
    if ( 0 != error_status )
      *error_status |= 1;
    if ( sNextUTF8 )
      *sNextUTF8 = sUTF8;
    return 0;
  }

  e.m_error_status = 0;
  e.m_error_mask = error_mask;
  e.m_error_code_point = error_code_point;

  output_count = 0;

  for ( i = 0; i < sUTF8_count; i += j )
  {
    j = ON_DecodeUTF8(sUTF8+i,sUTF8_count-i,&e,&u);
    if ( j <= 0 )
      break;
    k = ON_EncodeUTF16(u,w);
    if ( 0 != sUTF16 )
    {
      if ( output_count + k > sUTF16_count )
      {
        e.m_error_status |= 2;
        break;
      }
      sUTF16[output_count] = w[0];
      if ( 2 == k )
        sUTF16[output_count+1] = w[1];
    }
    output_count += k;
  }

  if ( 0 != sUTF16 && output_count < sUTF16_count)
    sUTF16[output_count] = 0;
  if ( sNextUTF8 )
    *sNextUTF8 = sUTF8+i;
  if ( error_status )
    *error_status = e.m_error_status;
  
  return output_count;
}

int ON_ConvertUTF8ToUTF32(
    int bTestByteOrder,
    const ON__UINT8* sUTF8,
    int sUTF8_count,
    ON__UINT32* sUTF32,
    int sUTF32_count,
    unsigned int* error_status,
    unsigned int error_mask,
    ON__UINT32 error_code_point,
    const ON__UINT8** sNextUTF8
    )
{
  int i, j, output_count;
  ON__UINT32 u = 0;
  struct ON_UnicodeErrorParameters e;

  if ( 0 != error_status )
    *error_status = 0;

  if ( -1 == sUTF8_count && 0 != sUTF8 )
  {
    for ( sUTF8_count = 0; 0 != sUTF8[sUTF8_count]; sUTF8_count++)
    {
      // empty for body
    }
  }

  if ( 0 == sUTF8 || sUTF8_count < 0 )
  {
    if ( 0 != error_status )
      *error_status |= 1;
    if ( sNextUTF8 )
      *sNextUTF8 = sUTF8;
    return 0;
  }

  if ( bTestByteOrder && ON_IsUTF8ByteOrderMark(sUTF8,sUTF8_count) )
  {
    // skip UTF-8 byte order element
    sUTF8_count -= 3;
    sUTF8 += 3;
  }

  if ( 0 == sUTF32_count )
  {
    sUTF32 = 0;
    sUTF32_count = 2147483647; // maximum value of a 32-bit signed int
  }
  else if ( 0 == sUTF32 )
  {
    if ( 0 != error_status )
      *error_status |= 1;
    if ( sNextUTF8 )
      *sNextUTF8 = sUTF8;
    return 0;
  }

  e.m_error_status = 0;
  e.m_error_mask = error_mask;
  e.m_error_code_point = error_code_point;

  output_count = 0;

  for ( i = 0; i < sUTF8_count; i += j )
  {
    j = ON_DecodeUTF8(sUTF8+i,sUTF8_count-i,&e,&u);
    if ( j <= 0 )
      break;
    if ( 0 != sUTF32 )
    {
      if ( output_count >= sUTF32_count )
      {
        e.m_error_status |= 2;
        break;
      }
      sUTF32[output_count] = u;
    }
    output_count++;
  }

  if ( 0 != sUTF32 && output_count < sUTF32_count)
    sUTF32[output_count] = 0;
  if ( sNextUTF8 )
    *sNextUTF8 = sUTF8+i;
  if ( error_status )
    *error_status = e.m_error_status;
  
  return output_count;
}

int ON_ConvertUTF16ToUTF8(
    int bTestByteOrder,
    const ON__UINT16* sUTF16,
    int sUTF16_count,
    ON__UINT8* sUTF8,
    int sUTF8_count,
    unsigned int* error_status,
    unsigned int error_mask,
    ON__UINT32 error_code_point,
    const ON__UINT16** sNextUTF16
    )
{
  int i, j, k, output_count, bSwapBytes;
  ON__UINT32 u;
  ON__UINT8 s[6];
  struct ON_UnicodeErrorParameters e;

  if ( 0 != error_status )
    *error_status = 0;

  if ( -1 == sUTF16_count && 0 != sUTF16 )
  {
    for ( sUTF16_count = 0; 0 != sUTF16[sUTF16_count]; sUTF16_count++)
    {
      // empty for body
    }
  }

  if ( 0 == sUTF16 || sUTF16_count < 0 )
  {
    if ( 0 != error_status )
      *error_status |= 1;
    if ( sNextUTF16 )
      *sNextUTF16 = sUTF16;
    return 0;
  }

  if ( 0 == sUTF8_count )
  {
    sUTF8 = 0;
    sUTF8_count = 2147483647; // maximum value of a 32-bit signed int
  }
  else if ( 0 == sUTF8 )
  {
    if ( 0 != error_status )
      *error_status |= 1;
    if ( sNextUTF16 )
      *sNextUTF16 = sUTF16;
    return 0;
  }

  bSwapBytes = false;
  if ( bTestByteOrder && sUTF16_count > 0 )
  {
    if ( 0xFEFF == sUTF16[0] )
    {
      // skip BOM
      sUTF16_count--;
      sUTF16++;
    }
    else if ( 0xFFFE == sUTF16[0] )
    {
      // skip BOM and swap bytes in rest of sUTF16
      bSwapBytes = true;
      sUTF16_count--;
      sUTF16++;
    }
  }

  e.m_error_status = 0;
  e.m_error_mask = error_mask;
  e.m_error_code_point = error_code_point;

  output_count = 0;

  if ( bSwapBytes )
  {
    for ( i = 0; i < sUTF16_count; i += j )
    {
      j = ON_DecodeSwapByteUTF16(sUTF16+i,sUTF16_count-i,&e,&u);
      if ( j <= 0 )
        break;
      k = ON_EncodeUTF8(u,s);
      if ( 0 != sUTF8 )
      {
        if ( output_count + k > sUTF8_count )
        {
          e.m_error_status |= 2;
          break;
        }
        memcpy(sUTF8+output_count,s,k*sizeof(sUTF8[0]));
      }
      output_count += k;
    }
  }
  else
  {
    for ( i = 0; i < sUTF16_count; i += j )
    {
      j = ON_DecodeUTF16(sUTF16+i,sUTF16_count-i,&e,&u);
      if ( j <= 0 )
        break;
      k = ON_EncodeUTF8(u,s);
      if ( 0 != sUTF8 )
      {
        if ( output_count + k > sUTF8_count )
        {
          e.m_error_status |= 2;
          break;
        }
        memcpy(sUTF8+output_count,s,k*sizeof(sUTF8[0]));
      }
      output_count += k;
    }
  }
  if ( 0 != sUTF8 && output_count < sUTF8_count)
    sUTF8[output_count] = 0;
  if ( sNextUTF16 )
    *sNextUTF16 = sUTF16+i;
  if ( error_status )
    *error_status = e.m_error_status;
  
  return output_count;
}

int ON_ConvertUTF16ToUTF32(
    int bTestByteOrder,
    const ON__UINT16* sUTF16,
    int sUTF16_count,
    unsigned int* sUTF32,
    int sUTF32_count,
    unsigned int* error_status,
    unsigned int error_mask,
    ON__UINT32 error_code_point,
    const ON__UINT16** sNextUTF16
    )
{
  int i, j, output_count, bSwapBytes;
  ON__UINT32 u;
  struct ON_UnicodeErrorParameters e;

  if ( 0 != error_status )
    *error_status = 0;

  if ( -1 == sUTF16_count && 0 != sUTF16 )
  {
    for ( sUTF16_count = 0; 0 != sUTF16[sUTF16_count]; sUTF16_count++)
    {
      // empty for body
    }
  }

  if ( 0 == sUTF16 || sUTF16_count < 0 )
  {
    if ( 0 != error_status )
      *error_status |= 1;
    if ( sNextUTF16 )
      *sNextUTF16 = sUTF16;
    return 0;
  }

  if ( 0 == sUTF32_count )
  {
    sUTF32 = 0;
    sUTF32_count = 2147483647; // maximum value of a 32-bit signed int
  }
  else if ( 0 == sUTF32 )
  {
    if ( 0 != error_status )
      *error_status |= 1;
    if ( sNextUTF16 )
      *sNextUTF16 = sUTF16;
    return 0;
  }

  bSwapBytes = false;
  if ( bTestByteOrder && sUTF16_count > 0 )
  {
    if ( 0xFEFF == sUTF16[0] )
    {
      // skip BOM
      sUTF16_count--;
      sUTF16++;
    }
    else if ( 0xFFFE == sUTF16[0] )
    {
      // skip BOM and swap bytes in rest of sUTF16
      bSwapBytes = true;
      sUTF16_count--;
      sUTF16++;
    }
  }

  e.m_error_status = 0;
  e.m_error_mask = error_mask;
  e.m_error_code_point = error_code_point;

  output_count = 0;

  if ( bSwapBytes )
  {
    for ( i = 0; i < sUTF16_count; i += j )
    {
      j = ON_DecodeSwapByteUTF16(sUTF16+i,sUTF16_count-i,&e,&u);
      if ( j <= 0 )
        break;
      if ( 0 != sUTF32 )
      {
        if ( output_count >= sUTF32_count )
        {
          e.m_error_status |= 2;
          break;
        }
        sUTF32[output_count] = u;
      }
      output_count++;
    }
  }
  else
  {
    for ( i = 0; i < sUTF16_count; i += j )
    {
      j = ON_DecodeUTF16(sUTF16+i,sUTF16_count-i,&e,&u);
      if ( j <= 0 )
        break;
      if ( 0 != sUTF32 )
      {
        if ( output_count >= sUTF32_count )
        {
          e.m_error_status |= 2;
          break;
        }
        sUTF32[output_count] = u;
      }
      output_count++;
    }
  }

  if ( 0 != sUTF32 && output_count < sUTF32_count)
    sUTF32[output_count] = 0;
  if ( sNextUTF16 )
    *sNextUTF16 = sUTF16+i;
  if ( error_status )
    *error_status = e.m_error_status;
  
  return output_count;
}

static ON__UINT32 SwapBytes32(ON__UINT32 u)
{
  ON__UINT8 b;
  ON__UINT8* p = (ON__UINT8*)&u;
  b = p[0]; p[0] = p[3]; p[3] = b;
  b = p[1]; p[1] = p[2]; p[2] = b;
  return u;
}

int ON_DecodeUTF32(
    const ON__UINT32* sUTF32,
    int sUTF32_count,
    struct ON_UnicodeErrorParameters* e,
    ON__UINT32* unicode_code_point
    )
{
  ON__UINT32 uhi, ulo;

  if ( 0 == sUTF32 || sUTF32_count <= 0 || 0 == unicode_code_point )
  {
    if ( e )
      e->m_error_status |= 1;
    return 0;
  }

  // special case for most common UTF-16 single element values
  if ( ( sUTF32[0] < 0xD800 ) || ( sUTF32[0] >= 0xE000 && sUTF32[0] <= 0x10FFFF) )
  {
    // valid UTF-32 encoding.
    *unicode_code_point = sUTF32[0];
    return 1;
  }

  // handle errors
  if ( 0 == e )
    return 0;

  if ( sUTF32_count >= 2 && sUTF32[0] < 0xDC00 && sUTF32[1] >=  0xDC00 && sUTF32[1] < 0xE000 )
  {
    // UTF-16 surrogate pair appears in UTF-32 array
    e->m_error_status |= 4;
    if ( 0 == (4 & e->m_error_mask) )
      return 0; // this error is not masked.

    uhi = sUTF32[0];
    ulo = sUTF32[1];
    *unicode_code_point = (uhi-0xD800)*0x400 + (ulo-0xDC00) + 0x10000;

    return 2; // error masked and reasonable value returned.
  }

  // bogus value
  e->m_error_status |= 16;
  if ( 16 != (16 & e->m_error_mask) || !ON_IsValidUnicodeCodePoint(e->m_error_code_point) )
  {
    // this error is not masked
    return 0;
  }

  *unicode_code_point = e->m_error_code_point;
  return 1; // error masked and e->m_error_code_point returnred.
}

int ON_DecodeSwapByteUTF32(
    const ON__UINT32* sUTF32,
    int sUTF32_count,
    struct ON_UnicodeErrorParameters* e,
    ON__UINT32* unicode_code_point
    )
{
  ON__UINT32 sUTF32swap[2];

  if ( 0 != sUTF32 && sUTF32_count > 0 )
  {
    sUTF32swap[0] = SwapBytes32(sUTF32[0]);
    if ( e && sUTF32_count > 1 )
    {
      // Get up to 2 elements to pass to the unswapped
      // decoder so that masked errors are uniformly
      // handled.
      sUTF32swap[1] = SwapBytes32(sUTF32[1]);
      sUTF32_count = 2;
    }
    sUTF32 = sUTF32swap;
  }

  return ON_DecodeUTF32(sUTF32,sUTF32_count,e,unicode_code_point);
}

int ON_ConvertUTF32ToUTF8(
    int bTestByteOrder,
    const ON__UINT32* sUTF32,
    int sUTF32_count,
    ON__UINT8* sUTF8,
    int sUTF8_count,
    unsigned int* error_status,
    unsigned int error_mask,
    ON__UINT32 error_code_point,
    const ON__UINT32** sNextUTF32
    )
{
  int i, k, output_count, bSwapBytes;
  ON__UINT32 u;
  ON__UINT8 s[6];
  struct ON_UnicodeErrorParameters e;

  if ( 0 != error_status )
    *error_status = 0;

  if ( -1 == sUTF32_count && 0 != sUTF32 )
  {
    for ( sUTF32_count = 0; 0 != sUTF32[sUTF32_count]; sUTF32_count++)
    {
      // empty for body
    }
  }

  if ( 0 == sUTF32 || sUTF32_count < 0 )
  {
    if ( 0 != error_status )
      *error_status |= 1;
    if ( sNextUTF32 )
      *sNextUTF32 = sUTF32;
    return 0;
  }

  if ( 0 == sUTF8_count )
  {
    sUTF8 = 0;
    sUTF8_count = 2147483647; // maximum value of a 32-bit signed int
  }
  else if ( 0 == sUTF8 )
  {
    if ( 0 != error_status )
      *error_status |= 1;
    if ( sNextUTF32 )
      *sNextUTF32 = sUTF32;
    return 0;
  }

  bSwapBytes = false;
  if ( bTestByteOrder && sUTF32_count > 0 )
  {
    if ( 0x0000FEFF == sUTF32[0] )
    {
      // skip BOM
      sUTF32_count--;
      sUTF32++;
    }
    else if ( 0xFFFE0000 == sUTF32[0] )
    {
      // skip BOM and swap bytes in rest of sUTF32
      bSwapBytes = true;
      sUTF32_count--;
      sUTF32++;
    }
  }

  e.m_error_status = 0;
  e.m_error_mask = error_mask;
  e.m_error_code_point = error_code_point;

  output_count = 0;

  for ( i = 0; i < sUTF32_count; i++ )
  {
    u = bSwapBytes ? SwapBytes32(sUTF32[i]) : sUTF32[i];
    if ( !ON_IsValidUnicodeCodePoint(u) )
    {
      e.m_error_status |= 16;
      if ( 16 != (16 & e.m_error_mask) )
        break;
      if ( !ON_IsValidUnicodeCodePoint(e.m_error_code_point) )
        break;
      u = e.m_error_code_point;
    }
    k = ON_EncodeUTF8(u,s);
    if ( 0 != sUTF8 )
    {
      if ( output_count + k > sUTF8_count )
      {
        e.m_error_status |= 2;
        break;
      }
      memcpy(sUTF8+output_count,s,k*sizeof(sUTF8[0]));
    }
    output_count += k;
  }

  if ( 0 != sUTF8 && output_count < sUTF8_count)
    sUTF8[output_count] = 0;
  if ( sNextUTF32 )
    *sNextUTF32 = sUTF32+i;
  if ( error_status )
    *error_status = e.m_error_status;
  
  return output_count;
}

int ON_ConvertUTF32ToUTF16(
    int bTestByteOrder,
    const ON__UINT32* sUTF32,
    int sUTF32_count,
    ON__UINT16* sUTF16,
    int sUTF16_count,
    unsigned int* error_status,
    unsigned int error_mask,
    ON__UINT32 error_code_point,
    const ON__UINT32** sNextUTF32
    )
{
  int i, k, output_count, bSwapBytes;
  ON__UINT32 u;
  ON__UINT16 w[2];
  struct ON_UnicodeErrorParameters e;

  if ( 0 != error_status )
    *error_status = 0;

  if ( -1 == sUTF32_count && 0 != sUTF32 )
  {
    for ( sUTF32_count = 0; 0 != sUTF32[sUTF32_count]; sUTF32_count++)
    {
      // empty for body
    }
  }

  if ( 0 == sUTF32 || sUTF32_count < 0 )
  {
    if ( 0 != error_status )
      *error_status |= 1;
    if ( sNextUTF32 )
      *sNextUTF32 = sUTF32;
    return 0;
  }

  if ( 0 == sUTF16_count )
  {
    sUTF16 = 0;
    sUTF16_count = 2147483647; // maximum value of a 32-bit signed int
  }
  else if ( 0 == sUTF16 )
  {
    if ( 0 != error_status )
      *error_status |= 1;
    if ( sNextUTF32 )
      *sNextUTF32 = sUTF32;
    return 0;
  }

  bSwapBytes = false;
  if ( bTestByteOrder && sUTF32_count > 0 )
  {
    if ( 0x0000FEFF == sUTF32[0] )
    {
      // skip BOM
      sUTF32_count--;
      sUTF32++;
    }
    else if ( 0xFFFE0000 == sUTF32[0] )
    {
      // skip BOM and swap bytes in rest of sUTF32
      bSwapBytes = true;
      sUTF32_count--;
      sUTF32++;
    }
  }

  e.m_error_status = 0;
  e.m_error_mask = error_mask;
  e.m_error_code_point = error_code_point;

  output_count = 0;

  for ( i = 0; i < sUTF32_count; i++ )
  {
    u = bSwapBytes ? SwapBytes32(sUTF32[i]) : sUTF32[i];
    if ( !ON_IsValidUnicodeCodePoint(u) )
    {
      e.m_error_status |= 16;
      if ( 16 != (16 & e.m_error_mask) )
        break;
      if ( !ON_IsValidUnicodeCodePoint(e.m_error_code_point) )
        break;
      u = e.m_error_code_point;
    }
    k = ON_EncodeUTF16(u,w);
    if ( 0 != sUTF16 )
    {
      if ( output_count + k > sUTF16_count )
      {
        e.m_error_status |= 2;
        break;
      }
      sUTF16[output_count] = w[0];
      if ( 2 == k )
        sUTF16[output_count+1] = w[1];
    }
    output_count += k;
  }

  if ( 0 != sUTF16 && output_count < sUTF16_count)
    sUTF16[output_count] = 0;
  if ( sNextUTF32 )
    *sNextUTF32 = sUTF32+i;
  if ( error_status )
    *error_status = e.m_error_status;
  
  return output_count;
}


int ON_ConvertUTF32ToUTF32(
    int bTestByteOrder,
    const ON__UINT32* sUTF16,
    int sUTF16_count,
    unsigned int* sUTF32,
    int sUTF32_count,
    unsigned int* error_status,
    unsigned int error_mask,
    ON__UINT32 error_code_point,
    const ON__UINT32** sNextUTF16
    )
{
  int i, j, output_count, bSwapBytes;
  ON__UINT32 u;
  struct ON_UnicodeErrorParameters e;

  if ( 0 != error_status )
    *error_status = 0;

  if ( -1 == sUTF16_count && 0 != sUTF16 )
  {
    for ( sUTF16_count = 0; 0 != sUTF16[sUTF16_count]; sUTF16_count++)
    {
      // empty for body
    }
  }

  if ( 0 == sUTF16 || sUTF16_count < 0 )
  {
    if ( 0 != error_status )
      *error_status |= 1;
    if ( sNextUTF16 )
      *sNextUTF16 = sUTF16;
    return 0;
  }

  if ( 0 == sUTF32_count )
  {
    sUTF32 = 0;
    sUTF32_count = 2147483647; // maximum value of a 32-bit signed int
  }
  else if ( 0 == sUTF32 )
  {
    if ( 0 != error_status )
      *error_status |= 1;
    if ( sNextUTF16 )
      *sNextUTF16 = sUTF16;
    return 0;
  }

  bSwapBytes = false;
  if ( bTestByteOrder && sUTF16_count > 0 )
  {
    if ( 0x0000FEFF == sUTF16[0] )
    {
      // skip BOM
      sUTF16_count--;
      sUTF16++;
    }
    else if ( 0xFFFE0000 == sUTF16[0])
    {
      // skip BOM and swap bytes in rest of sUTF16
      bSwapBytes = true;
      sUTF16_count--;
      sUTF16++;
    }
  }

  e.m_error_status = 0;
  e.m_error_mask = error_mask;
  e.m_error_code_point = error_code_point;

  output_count = 0;

  if ( bSwapBytes )
  {
    for ( i = 0; i < sUTF16_count; i += j )
    {
      j = ON_DecodeSwapByteUTF32(sUTF16+i,sUTF16_count-i,&e,&u);
      if ( j <= 0 )
        break;
      if ( 0 != sUTF32 )
      {
        if ( output_count >= sUTF32_count )
        {
          e.m_error_status |= 2;
          break;
        }
        sUTF32[output_count] = u;
      }
      output_count++;
    }
  }
  else
  {
    for ( i = 0; i < sUTF16_count; i += j )
    {
      j = ON_DecodeUTF32(sUTF16+i,sUTF16_count-i,&e,&u);
      if ( j <= 0 )
        break;
      if ( 0 != sUTF32 )
      {
        if ( output_count >= sUTF32_count )
        {
          e.m_error_status |= 2;
          break;
        }
        sUTF32[output_count] = u;
      }
      output_count++;
    }
  }

  if ( 0 != sUTF32 && output_count < sUTF32_count)
    sUTF32[output_count] = 0;
  if ( sNextUTF16 )
    *sNextUTF16 = sUTF16+i;
  if ( error_status )
    *error_status = e.m_error_status;
  
  return output_count;
}

ON_DECL
int ON_ConvertWideCharToUTF8(
    int bTestByteOrder,
    const wchar_t* sWideChar,
    int sWideChar_count,
    char* sUTF8,
    int sUTF8_count,
    unsigned int* error_status,
    unsigned int error_mask,
    ON__UINT32 error_code_point,
    const wchar_t** sNextWideChar
    )
{
  int rc;

  switch(sizeof(sWideChar[0]))
  {
  case sizeof(ON__UINT16):
    // assume wchar_t strings are UTF-16 encoded
    rc = ON_ConvertUTF16ToUTF8(
            bTestByteOrder,
            (const ON__UINT16*)sWideChar,sWideChar_count,
            (ON__UINT8*)sUTF8,sUTF8_count,
            error_status,error_mask,error_code_point,
            (const ON__UINT16**)sNextWideChar
            );
    break;

  case sizeof(ON__UINT32):
    // assume wchar_t strings are UTF-32 encoded
    rc = ON_ConvertUTF32ToUTF8(
            bTestByteOrder,
            (const ON__UINT32*)sWideChar,sWideChar_count,
            (ON__UINT8*)sUTF8,sUTF8_count,
            error_status,error_mask,error_code_point,
            (const ON__UINT32**)sNextWideChar
            );
    break;

  default:
    rc = 0;
  }

  return rc;
}


ON_DECL
int ON_ConvertWideCharToUTF32(
    int bTestByteOrder,
    const wchar_t* sWideChar,
    int sWideChar_count,
    ON__UINT32* sUTF32,
    int sUTF32_count,
    unsigned int* error_status,
    unsigned int error_mask,
    ON__UINT32 error_code_point,
    const wchar_t** sNextWideChar
    )
{
  int rc;

  switch(sizeof(sWideChar[0]))
  {
  case sizeof(ON__UINT8):
    // assume wchar_t strings are UTF-8 encoded
    rc = ON_ConvertUTF8ToUTF32(
            bTestByteOrder,
            (const ON__UINT8*)sWideChar,sWideChar_count,
            sUTF32,sUTF32_count,
            error_status,error_mask,error_code_point,
            (const ON__UINT8**)sNextWideChar
            );
    break;

  case sizeof(ON__UINT16):
    // assume wchar_t strings are UTF-16 encoded
    rc = ON_ConvertUTF16ToUTF32(
            bTestByteOrder,
            (const ON__UINT16*)sWideChar,sWideChar_count,
            sUTF32,sUTF32_count,
            error_status,error_mask,error_code_point,
            (const ON__UINT16**)sNextWideChar
            );
    break;

  case sizeof(ON__UINT32):
    // assume wchar_t strings are UTF-32 encoded
    rc = ON_ConvertUTF32ToUTF32(
            bTestByteOrder,
            (const ON__UINT32*)sWideChar,sWideChar_count,
            sUTF32,sUTF32_count,
            error_status,error_mask,error_code_point,
            (const ON__UINT32**)sNextWideChar
            );
     break;

  default:
    rc = 0;
  }

  return rc;
}


ON_DECL
int ON_ConvertUTF8ToWideChar(
    int bTestByteOrder,
    const char* sUTF8,
    int sUTF8_count,
    wchar_t* sWideChar,
    int sWideChar_count,
    unsigned int* error_status,
    unsigned int error_mask,
    ON__UINT32 error_code_point,
    const char** sNextUTF8
    )
{
  int rc;

  switch(sizeof(sWideChar[0]))
  {
  case sizeof(ON__UINT16):
    // assume wchar_t strings are UTF-16 encoded
    rc = ON_ConvertUTF8ToUTF16(
            bTestByteOrder,
            (const ON__UINT8*)sUTF8,sUTF8_count,
            (ON__UINT16*)sWideChar,sWideChar_count,
            error_status,error_mask,error_code_point,
            (const ON__UINT8**)sNextUTF8
            );
    break;

  case sizeof(ON__UINT32):
    // assume wchar_t strings are UTF-32 encoded
    rc = ON_ConvertUTF8ToUTF32(
            bTestByteOrder,
            (const ON__UINT8*)sUTF8,sUTF8_count,
            (ON__UINT32*)sWideChar,sWideChar_count,
            error_status,error_mask,error_code_point,
            (const ON__UINT8**)sNextUTF8
            );
    break;

  default:
    if (error_status)
      *error_status = 1;
    if (sNextUTF8)
      *sNextUTF8 = sUTF8;
    rc = 0;
  }

  return rc;
}
