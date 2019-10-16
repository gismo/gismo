/** @file gsKDTree.h

    @brief Provides declaration of kd-tree interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): Junjie Dong (junjiedong.umich@gmail.com),
               minor adjustments to G+Smo by M. Möller

    The original source code can be found at
    https://github.com/junjiedong/KDTree
*/

#pragma once

#include <stdexcept>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm>

#include <gsUtils/gsBoundedPriorityQueue.h>

namespace gismo
{

template <class T>
struct gsKDTreeTraits
{
  static inline std::size_t size() { return 1; }
  
  static inline bool islhalf(const T& lhs, const T& rhs, std::size_t axis) { return lhs[axis] < rhs[axis]; }
  
  static inline double fabs(const T& lhs, const T& rhs, std::size_t axis) { return std::abs(lhs[axis] - rhs[axis]); }
  
  static inline double distance(const T& lhs, const T& rhs)
  {
    double result = 0.0;
    for (std::size_t i = 0; i < size(); ++i) {
      result += fabs(lhs, rhs, i);
    }
    return result;
  }
};
  
/**
   \brief An interface representing a kd-tree in some number of dimensions

   The tree can be constructed from a set of data, (key,value)-pairs,
   and then queried for membership and k-nearest neighbors.
 */  
template <class KeyType, class ValueType>
class gsKDTree {
public:

  // Constructs an empty gsKDTree.
  gsKDTree();

  // Efficiently build a balanced kd-tree from a large set of data
  gsKDTree(std::vector<std::pair<KeyType, ValueType> >& data);
  
  // Frees up all the dynamically allocated resources
  ~gsKDTree();

  // Frees up all the dynamically allocated resources
  void clear();
  
  // Deep-copies the contents of another gsKDTree into this one.
  gsKDTree(const gsKDTree& other);

  // Deep-copies the contents of another gsKDTree into this one.
  gsKDTree& operator=(const gsKDTree& other);
  
  // Returns the dimension of the data stored in this gsKDTree.
  std::size_t dimension() const;

  // Returns the number of elements in the kd-tree.
  std::size_t size() const;

  // Returns true if this gsKDTree is empty and false otherwise.
  bool empty() const;

  // Returns true if the specified key is contained in the gsKDTree.
  bool contains(const KeyType& key) const;
  
  /*
   * Inserts the data with the given key into the gsKDTree,
   * associating it with the specified value. If another data element
   * with the same key already existed in the tree, the new value will
   * overwrite the existing one.
   */
  void insert(const KeyType& key, const ValueType& value=ValueType());
  
  /*
   * Returns a reference to the value associated with the data stored
   * under the given key in the gsKDTree. If the key does not exist,
   * then it is added to the gsKDTree using the default value of
   * ValueType as its value.
   */
  ValueType& operator[](const KeyType& key);
  
  /*
   * Returns a reference to the value associated with the given
   * key. If the key is not in the tree, this function throws an
   * out_of_range exception.
   */
  ValueType& at(const KeyType& key);
  const ValueType& at(const KeyType& key) const;
  
  /*
   * Given a key and an integer k, finds the k data elements in the
   * gsKDTree nearest to the data element associated with the given
   * key and returns the most common value associated with those data
   * elements. In the event of a tie, one of the most frequent value
   * will be chosen.
   */
  ValueType kNNValue(const KeyType& key, std::size_t k) const;

  /*
   * Given a key and an integer k, finds the k data elements in the
   * gsKDTree nearest to the data element associated with the given
   * key and returns a reference to the most common value associated
   * with those data elements. In the event of a tie, one of the most
   * frequent value will be chosen.
   */
  ValueType& kNNValue(const KeyType& key, std::size_t k);

  /// Prints the object as a string.
  void print(std::ostream &os) const;

  /// Print (as string) operator
  friend std::ostream &operator<<(std::ostream &os, const gsKDTree &obj)
  {
    obj.print(os);
    return os;
  }
  
private:
  struct Node {
    KeyType point;
    Node *left;
    Node *right;
    int level;  // level of the node in the tree, starts at 0 for the root
    ValueType value;
    Node(const KeyType& _key, int _level, const ValueType& _value=ValueType()):
      point(_key), left(NULL), right(NULL), level(_level), value(_value) {}
  };
  
  // Root node of the gsKDTree
  Node* root_;
  
  // Number of points in the gsKDTree
  std::size_t size_;
  
  /*
   * Recursively build a subtree that satisfies the kd-tree invariant using points in [start, end)
   * At each level, we split points into two halves using the median of the points as pivot
   * The root of the subtree is at level 'currLevel'
   * O(n) time partitioning algorithm is used to locate the median element
   */
  Node* buildTree(typename std::vector<std::pair<KeyType, ValueType> >::iterator start,
                  typename std::vector<std::pair<KeyType, ValueType> >::iterator end,
                  int currLevel);
  
  /*
   * Returns the Node that contains element with given key if it is present in subtree 'currNode'
   * Returns the Node below where key should be inserted if key is not in the subtree
   */
  Node* findNode(Node* currNode, const KeyType& key) const;
  
  // Recursive helper method for kNNValue(key, k)
  void nearestNeighborRecurse(const Node* currNode,
                              const KeyType& key,
                              gsBoundedPriorityQueue<ValueType>& bpq) const;

  // Recursive helper method for kNNValue(key, k)
  void nearestNeighborRecurse(const Node* currNode,
                              const KeyType& key,
                              gsBoundedPriorityQueue<ValueType*>& bpq) const;
  
  /*
   * Recursive helper method for copy constructor and assignment operator
   * Deep copies tree 'root' and returns the root of the copied tree
   */
  Node* deepcopyTree(Node* root);
  
  // Recursively free up all resources of subtree rooted at 'currNode'
  void freeResource(Node* currNode);

}; // class gsKDTree

template <class KeyType, class ValueType>
gsKDTree<KeyType, ValueType>::gsKDTree() :
  root_(NULL), size_(0) { }

template <class KeyType, class ValueType>
typename gsKDTree<KeyType, ValueType>::Node*
gsKDTree<KeyType, ValueType>::deepcopyTree(typename gsKDTree<KeyType, ValueType>::Node* root)
{
  if (root == NULL) return NULL;
  Node* newRoot = new Node(*root);
  newRoot->left = deepcopyTree(root->left);
  newRoot->right = deepcopyTree(root->right);
  return newRoot;
}

template <class KeyType, class ValueType>
typename gsKDTree<KeyType, ValueType>::Node*
gsKDTree<KeyType, ValueType>::buildTree(typename std::vector<std::pair<KeyType, ValueType> >::iterator start,
                                        typename std::vector<std::pair<KeyType, ValueType>>::iterator  end,
                                        int currLevel)
{
  if (start >= end) return NULL; // empty tree
  
  int axis = currLevel % gsKDTreeTraits<KeyType>::size(); // the axis to split on
  auto cmp = [axis](const std::pair<KeyType, ValueType>& p1,
                    const std::pair<KeyType, ValueType>& p2) {
    return p1.first[axis] < p2.first[axis];
  };
  std::size_t len = end - start;
  auto mid = start + len / 2;
  std::nth_element(start, mid, end, cmp); // linear time partition
  
  // move left (if needed) so that all the equal points are to the right
  // The tree will still be balanced as long as there aren't many points that are equal along each axis
  while (mid > start && (mid - 1)->first[axis] == mid->first[axis]) {
    --mid;
  }
  
  Node* newNode = new Node(mid->first, currLevel, mid->second);
  newNode->left = buildTree(start, mid, currLevel + 1);
  newNode->right = buildTree(mid + 1, end, currLevel + 1);
  return newNode;
}

template <class KeyType, class ValueType>
gsKDTree<KeyType, ValueType>::gsKDTree(std::vector<std::pair<KeyType, ValueType> >& data)
{
  root_ = buildTree(data.begin(), data.end(), 0);
  size_ = data.size();
}

template <class KeyType, class ValueType>
gsKDTree<KeyType, ValueType>::gsKDTree(const gsKDTree& rhs)
{
  root_ = deepcopyTree(rhs.root_);
  size_ = rhs.size_;
}

template <class KeyType, class ValueType>
gsKDTree<KeyType, ValueType>& gsKDTree<KeyType, ValueType>::operator=(const gsKDTree& rhs)
{
  if (this != &rhs) { // make sure we don't self-assign
    freeResource(root_);
    root_ = deepcopyTree(rhs.root_);
    size_ = rhs.size_;
  }
  return *this;
}

template <class KeyType, class ValueType>
void gsKDTree<KeyType, ValueType>::freeResource(typename gsKDTree<KeyType, ValueType>::Node* currNode)
{
  if (currNode == NULL) return;
  freeResource(currNode->left);
  freeResource(currNode->right);
  delete currNode;
}
  
template <class KeyType, class ValueType>
gsKDTree<KeyType, ValueType>::~gsKDTree()
{
  clear();
}

template <class KeyType, class ValueType>
void gsKDTree<KeyType, ValueType>::clear()
{
  freeResource(root_);
}
  
template <class KeyType, class ValueType>
std::size_t gsKDTree<KeyType, ValueType>::dimension() const
{
  return gsKDTreeTraits<KeyType>::size();
}

template <class KeyType, class ValueType>
std::size_t gsKDTree<KeyType, ValueType>::size() const
{
  return size_;
}

template <class KeyType, class ValueType>
bool gsKDTree<KeyType, ValueType>::empty() const
{
  return size_ == 0;
}
  
template <class KeyType, class ValueType>
typename gsKDTree<KeyType, ValueType>::Node*
gsKDTree<KeyType, ValueType>::findNode(typename gsKDTree<KeyType, ValueType>::Node* currNode,
                                       const KeyType& key) const
{
  if (currNode == NULL || currNode->point == key) return currNode;
  
  const KeyType& currPoint = currNode->point;
  int currLevel = currNode->level;
  if (gsKDTreeTraits<KeyType>::islhalf(key, currPoint, currLevel%gsKDTreeTraits<KeyType>::size()))
    {
      // recurse to the left side
      return currNode->left == NULL ? currNode : findNode(currNode->left, key);
    } else {
    // recurse to the right side
    return currNode->right == NULL ? currNode : findNode(currNode->right, key);
  }
}

template <class KeyType, class ValueType>
bool gsKDTree<KeyType, ValueType>::contains(const KeyType& key) const
{
  auto node = findNode(root_, key);
  return node != NULL && node->point == key;
}
  
template <class KeyType, class ValueType>
void gsKDTree<KeyType, ValueType>::insert(const KeyType& key, const ValueType& value)
{
  auto targetNode = findNode(root_, key);
  if (targetNode == NULL) { // this means the tree is empty
    root_ = new Node(key, 0, value);
    size_ = 1;
  } else {
    if (targetNode->point == key) { // key is already in the tree, simply update its value
      targetNode->value = value;
    } else { // construct a new node and insert it to the right place (child of targetNode)
      int currLevel = targetNode->level;
      Node* newNode = new Node(key, currLevel + 1, value);
      if (gsKDTreeTraits<KeyType>::islhalf(key, targetNode->point, currLevel%gsKDTreeTraits<KeyType>::size())) {
        targetNode->left = newNode;
      } else {
        targetNode->right = newNode;
      }
      ++size_;
    }
  }
}

template <class KeyType, class ValueType>
const ValueType& gsKDTree<KeyType, ValueType>::at(const KeyType& key) const
{
  auto node = findNode(root_, key);
  if (node == NULL || node->point != key) {
    throw std::out_of_range("Key not found in gsKDTree");
  } else {
    return node->value;
  }
}

template <class KeyType, class ValueType>
ValueType& gsKDTree<KeyType, ValueType>::at(const KeyType& key)
{
  const gsKDTree<KeyType, ValueType>& constThis = *this;
  return const_cast<ValueType&>(constThis.at(key));
}
  
template <class KeyType, class ValueType>
ValueType& gsKDTree<KeyType, ValueType>::operator[](const KeyType& key)
{
  auto node = findNode(root_, key);
  if (node != NULL && node->point == key) { // key is already in the tree
    return node->value;
  } else { // insert key with default ValueType value, and return reference to the new ValueType
    insert(key);
    if (node == NULL) return root_->value; // the new node is the root
    else return (node->left != NULL && node->left->point == key) ? node->left->value: node->right->value;
  }
}
  
template <class KeyType, class ValueType>
void gsKDTree<KeyType, ValueType>::nearestNeighborRecurse(const typename gsKDTree<KeyType, ValueType>::Node* currNode,
                                                          const KeyType& key,
                                                          gsBoundedPriorityQueue<ValueType>& bpq) const
{
  if (currNode == NULL) return;
  const KeyType& currPoint = currNode->point;

  // Add the current point to the BPQ if it is closer to 'key' that some point in the BPQ
  bpq.enqueue(currNode->value, gsKDTreeTraits<KeyType>::distance(key, currPoint));

  // Recursively search the half of the tree that contains Point 'key'
  int currLevel = currNode->level;
  bool isLeftTree;
  if (gsKDTreeTraits<KeyType>::islhalf(key, currPoint, currLevel%gsKDTreeTraits<KeyType>::size())) {
    nearestNeighborRecurse(currNode->left, key, bpq);
    isLeftTree = true;
  } else {
    nearestNeighborRecurse(currNode->right, key, bpq);
    isLeftTree = false;
  }
  
  if (bpq.size() < bpq.maxSize() ||
      gsKDTreeTraits<KeyType>::fabs(key, currPoint, currLevel%gsKDTreeTraits<KeyType>::size()) < bpq.worst()) {
    // Recursively search the other half of the tree if necessary
    if (isLeftTree) nearestNeighborRecurse(currNode->right, key, bpq);
    else nearestNeighborRecurse(currNode->left, key, bpq);
  }
}

template <class KeyType, class ValueType>
void gsKDTree<KeyType, ValueType>::nearestNeighborRecurse(const typename gsKDTree<KeyType, ValueType>::Node* currNode,
                                                          const KeyType& key,
                                                          gsBoundedPriorityQueue<ValueType*>& bpq) const
{
  if (currNode == NULL) return;
  const KeyType& currPoint = currNode->point;

  // Add the current point to the BPQ if it is closer to 'key' that some point in the BPQ
  bpq.enqueue(const_cast<ValueType*>(&(currNode->value)), gsKDTreeTraits<KeyType>::distance(key, currPoint));

  // Recursively search the half of the tree that contains Point 'key'
  int currLevel = currNode->level;
  bool isLeftTree;
  if (gsKDTreeTraits<KeyType>::islhalf(key, currPoint, currLevel%gsKDTreeTraits<KeyType>::size())) {
    nearestNeighborRecurse(currNode->left, key, bpq);
    isLeftTree = true;
  } else {
    nearestNeighborRecurse(currNode->right, key, bpq);
    isLeftTree = false;
  }
  
  if (bpq.size() < bpq.maxSize() ||
      gsKDTreeTraits<KeyType>::fabs(key, currPoint, currLevel%gsKDTreeTraits<KeyType>::size()) < bpq.worst()) {
    // Recursively search the other half of the tree if necessary
    if (isLeftTree) nearestNeighborRecurse(currNode->right, key, bpq);
    else nearestNeighborRecurse(currNode->left, key, bpq);
  }
}

template <class KeyType, class ValueType>
ValueType gsKDTree<KeyType, ValueType>::kNNValue(const KeyType& key, std::size_t k) const
{
  // BPQ with maximum size k
  gsBoundedPriorityQueue<ValueType> bpq(k); 
  if (empty()) throw std::out_of_range("gsKDTree is empty");

  // Recursively search the kd-tree with pruning
  nearestNeighborRecurse(root_, key, bpq);

  // Ensure finite values; non-standard 'distance' functions can be
  // used to exclude data elements that are close to the given key but
  // on the 'wrong' side of the hyperplane. This allows to exclude
  // nearest neighbours that are, e.g., smaller than the given key.
  if (!std::isfinite(bpq.best()))
      throw std::out_of_range("gsKDTree does not contain finite value");
  
  // Count occurrences of all ValueType in the kNN set
  std::unordered_map<ValueType, int> counter;
  while (!bpq.empty()) {
    ++counter[bpq.dequeueMin()];
  }

  // Return the most frequent element in the kNN set
  ValueType result;
  int cnt = -1;
  for (const auto &p : counter) {
    if (p.second > cnt) {
      result = p.first;
      cnt = p.second;
    }
  }
  return result;
}

template <class KeyType, class ValueType>
ValueType& gsKDTree<KeyType, ValueType>::kNNValue(const KeyType& key, std::size_t k)
{
  // BPQ with maximum size k
  gsBoundedPriorityQueue<ValueType*> bpq(k);
  if (empty())
    throw std::out_of_range("gsKDTree is empty");

  // Recursively search the kd-tree with pruning
  nearestNeighborRecurse(root_, key, bpq);
  
  // Ensure finite values; non-standard 'distance' functions can be
  // used to exclude data elements that are close to the given key but
  // on the 'wrong' side of the hyperplane. This allows to exclude
  // nearest neighbours that are, e.g., smaller than the given key.
  if (!std::isfinite(bpq.best()))
      throw std::out_of_range("gsKDTree does not contain finite value");
  
  // Count occurrences of all ValueType in the kNN set
  std::unordered_map<ValueType*, int> counter;
  while (!bpq.empty()) {
    ++counter[bpq.dequeueMin()];
  }

  // Return the most frequent element in the kNN set
  ValueType* result = nullptr;
  int cnt = -1;
  for (const auto &p : counter) {
    if (p.second > cnt) {
      result = p.first;
      cnt = p.second;
    }
  }
  return *result;
}

template <class KeyType, class ValueType>
void gsKDTree<KeyType, ValueType>::print(std::ostream& os) const
{
  os << "KD-tree: size= " << size() << ", dimension= " << dimension() << ".\n";
}
  
} //namespace gismo
