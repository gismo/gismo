/** @file gsKNNTree.h

    @brief Provides declaration of kd-tree interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Junjie Dong (junjiedong.umich@gmail.com),
    minor adjustments to G+Smo by M. Möller

    The original source code can be found at
    https://github.com/junjiedong/KNNTree
*/

#pragma once

#include <stdexcept>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm>

#include <gsUtils/gsKDTree.h>
#include <gsUtils/gsBoundedPriorityQueue.h>

namespace gismo
{

template <class T>
struct gsKNNTreeTraits
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
class gsKNNTree {

    struct NodeData { // Node datatype
        KeyType point;
        int level;  // level of the node in the tree, starts at 0 for the root
        ValueType value;
        NodeData(const KeyType& _key, int _level, const ValueType& _value=ValueType()):
        point(_key), level(_level), value(_value) {}
    };

    typedef gsKDTree<real_t,NodeData> Node;
    Node * m_root; /// Root node of the gsKDTree

public:

    // Constructs an empty gsKNNTree.
    gsKNNTree() : m_root(NULL) { }

    // Efficiently build a balanced kd-tree from a large set of data
    gsKNNTree(std::vector<std::pair<KeyType, ValueType> >& data)
    {
        m_root = buildTree(data.begin(), data.end(), 0);
    }

    // Frees up all the dynamically allocated resources
    ~gsKNNTree() { clear(); }

    // Frees up all the dynamically allocated resources
    void clear() { delete m_root; }

    const Node & root() const { return *m_root; }

    // Deep-copies the contents of another gsKNNTree into this one.
    gsKNNTree(const gsKNNTree& other)
    {
        m_root = new Node(*other.m_root);
    }

    // Deep-copies the contents of another gsKNNTree into this one.
    gsKNNTree& operator=(const gsKNNTree& other)
    {
        if (this != &other) // make sure we don't self-assign
            m_root = new Node(*other.m_root);
        return *this;
    }

    // Returns the dimension of the data stored in this gsKNNTree.
    std::size_t dimension() const
    { return gsKNNTreeTraits<KeyType>::size(); }

    // Returns the number of elements in the kd-tree.
    std::size_t size() const { return m_root->numNodes(); }

    // Returns true if this gsKNNTree is empty and false otherwise.
    bool empty() const { return(m_root == nullptr); }

    // Returns true if the specified key is contained in the gsKNNTree.
    bool contains(const KeyType& key) const
    {
        auto node = findNode(m_root, key);
        return node != NULL && node->point == key;
    }

    /*
     * Inserts the data with the given key into the gsKNNTree,
     * associating it with the specified value. If another data element
     * with the same key already existed in the tree, the new value will
     * overwrite the existing one.
     */
    void insert(const KeyType& key, const ValueType& value=ValueType());

    /*
     * Returns a reference to the value associated with the data stored
     * under the given key in the gsKNNTree. If the key does not exist,
     * then it is added to the gsKNNTree using the default value of
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
     * gsKNNTree nearest to the data element associated with the given
     * key and returns the most common value associated with those data
     * elements. In the event of a tie, one of the most frequent value
     * will be chosen.
     */
    ValueType kNNValue(const KeyType& key, std::size_t k) const;

    /*
     * Given a key and an integer k, finds the k data elements in the
     * gsKNNTree nearest to the data element associated with the given
     * key and returns a reference to the most common value associated
     * with those data elements. In the event of a tie, one of the most
     * frequent value will be chosen.
     */
    ValueType& kNNValue(const KeyType& key, std::size_t k);

    /// Prints the object as a string.
    void print(std::ostream &os) const;

    /// Print (as string) operator
    friend std::ostream &operator<<(std::ostream &os, const gsKNNTree &obj)
    {
        obj.print(os);
        return os;
    }

private:

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

}; // class gsKNNTree

template <class KeyType, class ValueType>
typename gsKNNTree<KeyType, ValueType>::Node*
gsKNNTree<KeyType, ValueType>::buildTree(typename std::vector<std::pair<KeyType, ValueType> >::iterator start,
                                         typename std::vector<std::pair<KeyType, ValueType>>::iterator  end,
                                         int currLevel)
{
    if (start >= end) return NULL; // empty tree

    int axis = currLevel % gsKNNTreeTraits<KeyType>::size(); // the axis to split on
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

    Node* newNode = new Node(axis, mid->first[axis],new NodeData(mid->first, currLevel, mid->second));

    newNode->left = buildTree(start, mid, currLevel + 1);
    if (nullptr == newNode->left) newNode->axis = -1;
    else newNode->left->parent = newNode;

    newNode->right = buildTree(mid + 1, end, currLevel + 1);
    if (nullptr == newNode->right) newNode->axis = -1;
    else newNode->right->parent = newNode;

    return newNode;
}

template <class KeyType, class ValueType>
typename gsKNNTree<KeyType, ValueType>::Node*
gsKNNTree<KeyType, ValueType>::findNode(typename gsKNNTree<KeyType, ValueType>::Node* currNode,
                                        const KeyType& key) const
{
    if (currNode == NULL || currNode->data().point == key) return currNode;

    const KeyType& currPoint = currNode->data().point; // point contains the position..
    // NOTE: split position is not used, it is equal to: currNode->data().point[currNode->axis]
    if (gsKNNTreeTraits<KeyType>::islhalf(key, currPoint, (currNode->axis+1)%gsKNNTreeTraits<KeyType>::size()))
    {
        // recurse to the left side
        return currNode->left == NULL ? currNode : findNode(currNode->left, key);
    } else {
        // recurse to the right side
        return currNode->right == NULL ? currNode : findNode(currNode->right, key);
    }
}

template <class KeyType, class ValueType>
void gsKNNTree<KeyType, ValueType>::insert(const KeyType& key, const ValueType& value)
{
    auto targetNode = findNode(m_root, key);
    if (targetNode == NULL) { // this means the tree is empty
        m_root = new Node();
        m_root->data = new NodeData(key, 0, value);
    } else {
        if (targetNode->data().point == key) { // key is already in the tree, simply update its value
            targetNode->data().value = value;
        } else { // construct a new node and insert it to the right place (child of targetNode)
            int currLevel = targetNode->data().level;
            Node* newNode = new Node ( (currLevel + 1)%currLevel%gsKNNTreeTraits<KeyType>::size(), 000, new NodeData(key, currLevel + 1, value) );

            if (gsKNNTreeTraits<KeyType>::islhalf(key, targetNode->data().point, currLevel%gsKNNTreeTraits<KeyType>::size())) {
                targetNode->left = newNode;
            } else {
                targetNode->right = newNode;
            }
        }
    }
}

template <class KeyType, class ValueType>
const ValueType& gsKNNTree<KeyType, ValueType>::at(const KeyType& key) const
{
    auto node = findNode(m_root, key);
    if (node == NULL || node->data().point != key) {
        throw std::out_of_range("Key not found in gsKNNTree");
    } else {
        return node->data().value;
    }
}

template <class KeyType, class ValueType>
ValueType& gsKNNTree<KeyType, ValueType>::at(const KeyType& key)
{
    const gsKNNTree<KeyType, ValueType>& constThis = *this;
    return const_cast<ValueType&>(constThis.at(key));
}

template <class KeyType, class ValueType>
ValueType& gsKNNTree<KeyType, ValueType>::operator[](const KeyType& key)
{
    auto node = findNode(m_root, key);
    if (node != NULL && node->data().point == key) { // key is already in the tree
        return node->data().value;
    } else { // insert key with default ValueType value, and return reference to the new ValueType
        insert(key);
        if (node == NULL) return m_root->data().value; // the new node is the root
        else return (node->left != NULL && node->left->data().point == key) ? node->left->data().value: node->right->data().value;
    }
}

template <class KeyType, class ValueType>
void gsKNNTree<KeyType, ValueType>::nearestNeighborRecurse(const typename gsKNNTree<KeyType, ValueType>::Node* currNode,
                                                           const KeyType& key,
                                                           gsBoundedPriorityQueue<ValueType>& bpq) const
{
    if (currNode == NULL) return;
    const KeyType& currPoint = currNode->data().point;

    // Add the current point to the BPQ if it is closer to 'key' that some point in the BPQ
    bpq.enqueue(currNode->data().value, gsKNNTreeTraits<KeyType>::distance(key, currPoint));

    // Recursively search the half of the tree that contains Point 'key'
    int currLevel = currNode->level;
    bool isLeftTree;
    if (gsKNNTreeTraits<KeyType>::islhalf(key, currPoint, currLevel%gsKNNTreeTraits<KeyType>::size())) {
        nearestNeighborRecurse(currNode->left, key, bpq);
        isLeftTree = true;
    } else {
        nearestNeighborRecurse(currNode->right, key, bpq);
        isLeftTree = false;
    }

    if (bpq.size() < bpq.maxSize() ||
        gsKNNTreeTraits<KeyType>::fabs(key, currPoint, currLevel%gsKNNTreeTraits<KeyType>::size()) < bpq.worst()) {
        // Recursively search the other half of the tree if necessary
        if (isLeftTree) nearestNeighborRecurse(currNode->right, key, bpq);
        else nearestNeighborRecurse(currNode->left, key, bpq);
    }
}

template <class KeyType, class ValueType>
void gsKNNTree<KeyType, ValueType>::nearestNeighborRecurse(const typename gsKNNTree<KeyType, ValueType>::Node* currNode,
                                                           const KeyType& key,
                                                           gsBoundedPriorityQueue<ValueType*>& bpq) const
{
    if (currNode == NULL) return;
    const KeyType& currPoint = currNode->data().point;

    // Add the current point to the BPQ if it is closer to 'key' that some point in the BPQ
    bpq.enqueue(const_cast<ValueType*>(&(currNode->data().value)), gsKNNTreeTraits<KeyType>::distance(key, currPoint));

    // Recursively search the half of the tree that contains Point 'key'
    int currLevel = currNode->level;
    bool isLeftTree;
    if (gsKNNTreeTraits<KeyType>::islhalf(key, currPoint, currLevel%gsKNNTreeTraits<KeyType>::size())) {
        nearestNeighborRecurse(currNode->left, key, bpq);
        isLeftTree = true;
    } else {
        nearestNeighborRecurse(currNode->right, key, bpq);
        isLeftTree = false;
    }

    if (bpq.size() < bpq.maxSize() ||
        gsKNNTreeTraits<KeyType>::fabs(key, currPoint, currLevel%gsKNNTreeTraits<KeyType>::size()) < bpq.worst()) {
        // Recursively search the other half of the tree if necessary
        if (isLeftTree) nearestNeighborRecurse(currNode->right, key, bpq);
        else nearestNeighborRecurse(currNode->left, key, bpq);
    }
}

template <class KeyType, class ValueType>
ValueType gsKNNTree<KeyType, ValueType>::kNNValue(const KeyType& key, std::size_t k) const
{
    // BPQ with maximum size k
    gsBoundedPriorityQueue<ValueType> bpq(k);
    if (empty()) throw std::out_of_range("gsKNNTree is empty");

    // Recursively search the kd-tree with pruning
    nearestNeighborRecurse(m_root, key, bpq);

    // Ensure finite values; non-standard 'distance' functions can be
    // used to exclude data elements that are close to the given key but
    // on the 'wrong' side of the hyperplane. This allows to exclude
    // nearest neighbours that are, e.g., smaller than the given key.
    if (!math::isfinite(bpq.best()))
        throw std::out_of_range("gsKNNTree does not contain finite value");

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
ValueType& gsKNNTree<KeyType, ValueType>::kNNValue(const KeyType& key, std::size_t k)
{
    // BPQ with maximum size k
    gsBoundedPriorityQueue<ValueType*> bpq(k);
    if (empty())
        throw std::out_of_range("gsKNNTree is empty");

    // Recursively search the kd-tree with pruning
    nearestNeighborRecurse(m_root, key, bpq);

    // Ensure finite values; non-standard 'distance' functions can be
    // used to exclude data elements that are close to the given key but
    // on the 'wrong' side of the hyperplane. This allows to exclude
    // nearest neighbours that are, e.g., smaller than the given key.
    if (!math::isfinite(bpq.best()))
        throw std::out_of_range("gsKNNTree does not contain finite value");

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
void gsKNNTree<KeyType, ValueType>::print(std::ostream& os) const
{
    os << "kNN-tree: size= " << size() << ", dimension= " << dimension() << ".\n";
}

} //namespace gismo
