/** @file gsBoundedPriorityQueue.h

    @brief Provides declaration of bounded priority queue

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): Keith Schwarz (htiek@cs.stanford.edu),
               minor adjustments to G+Smo by M. Möller

    The original source code can be found at
    @link http://www.keithschwarz.com/interesting/code/?dir=bounded-pqueue
*/

#pragma once

#include <map>
#include <algorithm>
#include <limits>

namespace gismo
{

/**
   \brief An implementation of the bounded priority queue abstraction.
   
   A bounded priority queue is in many ways like a regular priority
   queue.  It stores a collection of elements tagged with a real-
   valued priority, and allows for access to the element whose
   priority is the smallest.  However, unlike a regular priority
   queue, the number of elements in a bounded priority queue has
   a hard limit that is specified in the constructor.  Whenever an
   element is added to the bounded priority queue such that the
   size exceeds the maximum, the element with the highest priority
   value will be ejected from the bounded priority queue.  In this
   sense, a bounded priority queue is like a high score table for
   a video game that stores a fixed number of elements and deletes
   the least-important entry whenever a new value is inserted.
   
   When creating a bounded priority queue, you must specify the
   maximum number of elements to store in the queue as an argument
   to the constructor.

   For example:
   
   gsBoundedPriorityQueue<int> bpq(15); // Holds up to fifteen values.
   
   The maximum size of the bounded priority queue can be obtained
   using the maxSize() function, as in
   
   size_t k = bpq.maxSize();
   
   Beyond these restrictions, the bounded priority queue behaves
   similarly to other containers.  You can query its size using
   size() and check whether it is empty using empty().  You
   can enqueue an element into the bounded priority queue by
   writing
   
   bpq.enqueue(elem, priority);
   
   Note that after enqueuing the element, there is no guarantee
   that the value will actually be in the queue.  If the queue
   is full and the new element's priority exceeds the largest
   priority in the container, it will not be added.
   
   You can dequeue elements from a bounded priority queue using
   the dequeueMin() function, as in
   
   int val = bpq.dequeueMin();
   
   The bounded priority queue also allows you to query the min
   and max priorities of the values in the queue.  These values
   can be queried using the best() and worst() functions, which
   return the smallest and largest priorities in the queue,
   respectively.
*/
  template <class ValueType, class T=real_t>
class gsBoundedPriorityQueue {
public:
  // Constructor: gsBoundedPriorityQueue(size_t maxSize);
  // Usage: gsBoundedPriorityQueue<int> bpq(15);
  // --------------------------------------------------
  // Constructs a new, empty gsBoundedPriorityQueue with
  // maximum size equal to the constructor argument.
  ///
  explicit gsBoundedPriorityQueue(std::size_t maxSize);

  // void enqueue(const ValueType& value, T priority);
  // Usage: bpq.enqueue("Hi!", 2.71828);
  // --------------------------------------------------
  // Enqueues a new element into the gsBoundedPriorityQueue with
  // the specified priority. If this overflows the maximum
  // size of the queue, the element with the highest
  // priority will be deleted from the queue. Note that
  // this might be the element that was just added.
  void enqueue(const ValueType& value, T priority);

  // ValueType dequeueMin();
  // Usage: int val = bpq.dequeueMin();
  // --------------------------------------------------
  // Returns the element from the gsBoundedPriorityQueue with the
  // smallest priority value, then removes that element
  // from the queue.
  ValueType dequeueMin();

  // size_t size() const;
  // bool empty() const;
  // Usage: while (!bpq.empty()) { ... }
  // --------------------------------------------------
  // Returns the number of elements in the queue and whether
  // the queue is empty, respectively.
  std::size_t size() const;
  bool empty() const;

  // size_t maxSize() const;
  // Usage: size_t queueSize = bpq.maxSize();
  // --------------------------------------------------
  // Returns the maximum number of elements that can be
  // stored in the queue.
  std::size_t maxSize() const;

  // T best() const;
  // T worst() const;
  // Usage: T highestPriority = bpq.worst();
  // --------------------------------------------------
  // best() returns the smallest priority of an element
  // stored in the container (i.e. the priority of the
  // element that will be dequeued first using dequeueMin).
  // worst() returns the largest priority of an element
  // stored in the container.  If an element is enqueued
  // with a priority above this value, it will automatically
  // be deleted from the queue.  Both functions return
  // numeric_limits<T>::infinity() if the queue is
  // empty.
  T best()  const;
  T worst() const;

  /// Prints the object as a string.
  void print(std::ostream &os) const;

  /// Print (as string) operator
  friend std::ostream &operator<<(std::ostream &os, const gsBoundedPriorityQueue &obj)
  {
    obj.print(os);
    return os;
  }
    
private:
  // This class is layered on top of a multimap mapping from priorities
  // to elements with those priorities.
  std::multimap<T, ValueType> elems;
  std::size_t maximumSize;
    
}; // class gsBoundedPriorityQueue

// Constructor
template <class ValueType, class T>
gsBoundedPriorityQueue<ValueType, T>::gsBoundedPriorityQueue(std::size_t maxSize)
{
maximumSize = maxSize;
}

// enqueue adds the element to the map, then deletes the last element of the
// map if there size exceeds the maximum size.
template <class ValueType, class T>
void gsBoundedPriorityQueue<ValueType, T>::enqueue(const ValueType& value, T priority)
{
// Add the element to the collection.
elems.insert(std::make_pair(priority, value));

// If there are too many elements in the queue, drop off the last one.
if (size() > maxSize()) {
typename std::multimap<T, ValueType>::iterator last = elems.end();
--last; // Now points to highest-priority element
elems.erase(last);
}
}
  
// dequeueMin copies the lowest element of the map (the one pointed at by
// begin()) and then removes it.
template <class ValueType, class T>
ValueType gsBoundedPriorityQueue<ValueType, T>::dequeueMin()
{
// Copy the best value.
ValueType result = elems.begin()->second;

// Remove it from the map.
elems.erase(elems.begin());

return result;
}

// size() and empty() call directly down to the underlying map.
template <class ValueType, class T>
std::size_t gsBoundedPriorityQueue<ValueType, T>::size() const
{
return elems.size();
}

template <class ValueType, class T>
bool gsBoundedPriorityQueue<ValueType, T>::empty() const
{
return elems.empty();
}

// maxSize just returns the appropriate data member.
template <class ValueType, class T>
std::size_t gsBoundedPriorityQueue<ValueType, T>::maxSize() const
{
return maximumSize;
}

// The best() and worst() functions check if the queue is empty,
// and if so return infinity.
template <class ValueType, class T>
T gsBoundedPriorityQueue<ValueType, T>::best() const
{
return empty()? std::numeric_limits<T>::infinity() : elems.begin()->first;
}

template <class ValueType, class T>
T gsBoundedPriorityQueue<ValueType, T>::worst() const
{
return empty()? std::numeric_limits<T>::infinity() : elems.rbegin()->first;
}

template <class ValueType, class T>
void gsBoundedPriorityQueue<ValueType, T>::print(std::ostream& os) const
{
  os << "BPQ(" << size() << ") [" << best() << ".." << worst() << "].\n";
}
  
} // namespace gismo
