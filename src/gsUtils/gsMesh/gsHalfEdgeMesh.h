/** @file gsHalfEdgeMesh.h

    @brief Provides declaration of the gsHalfEdgeMesh class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl
*/

#pragma once

#include <gsUtils/gsMesh/gsMesh.h>
#include <queue>

namespace gismo
{
enum triangleVertexIndex
{
    first, second, third, error
};

/**
 * @brief gsHalfEdgeMesh is a gsMesh implementation that handels Halfedges
 *
 * According to Floater's algorithm the N points of the triangle mesh are given by x1, ..., xN where x1, ..., xn are inner points and x(n+1), ..., xN are boundary points.
 * In MeshInfo the points are stored without special order and then two vectors containig indices are constructed for easily accessing point xi s. t. xi is inner vertex for i<=n and boundary vertex for i > n.
 * In the following the 'vertexIndex' is used for the indices of the ordered vertices 1, ..., N and by 'internVertexIndex' the indices 0,...,N-1 of the unordered m_vertices vector is meant, s. t. *m_vertices[internVertexIndex] = x(vertexIndex).
 *
 * There are getter functions for the number of vertices (N), triangles, halfedges, inner vertices (n) and boundary vertices plus one can get a vertex with particular vertex index or the index from a particular vertex.
 * The vertices of the triangle mesh get ordered, s. t. the first n vertex indices have correspondent inner vertices and the last N-n vertex indices correspond to the boundary vertices ordered counter-clockwise.
 * */
template <class T>
class gsHalfEdgeMesh : public gsMesh<>
{
public:
    /**
 * @brief Class that maintains boundary of triangle mesh.
 *
 * A Boundary class object is given by a chain of the boundary halfedges.
 * The halfedges are ordered counter clockwise.
 *
 * An object of the class can be constructed by a vector of all unordered halfedges of the mesh.
 * There are methods to get the number of halfedges, number of vertices, length, halfedge lengths, first and last halfedge, vertex indices, and distances between two vertices.
 *
 * A boundary can be printed.
 *
 * */
    class Boundary
    {
    public:
        /**
     * @brief Class that maintains chain of halfedges.
     *
     * The halfedges of the chain are stored in a list, in order to easily insert new halfedges at the beginning and at the end of the list.
     *
     * There are functions for testing whether the chain is empty or closed.
     * Using getter functions the number of halfedges or vertices, as well as the length and the lengths of all halfedges can be returned.
     * Furthermore the vertex indices ordered like they appear in the chain and the first and last halfedge can be returned.
     * The distance between two vertices can be calculated.
     * It can be questioned wheter a vertex is contained in the chain and wheter another halfedge can be appended at the beginning or end.
     *
     * A chain can be printed by printing its halfedges.
     *
     * */
        class Chain
        {
        public:
            /**
             * @brief Class that maintains directed halfedges in any dimension.
             *
             * The class Halfedge represents a halfedge in any dimension given by its origin and end point indices and its length.
             * A halfedge can be constructed with desired origin and end point and length.
             * Two halfedges can be compared by == and !=.
             * There are getter function for origin, end and length.
             * There are functions to test wheter a halfedge can be appended to another one before or afterwards and wheter a halfedge is a twin of another.
             *
             * Halfedges can be outputted like (origin--end: length).
             * */
            class Halfedge
            {

            public:
                /**
                 * @brief Default constructor
                 * The default constructor sets indices (m_origin, m_end) and length (m_length) to 0.
                 * */
                Halfedge()
                    : m_origin(0), m_end(0), m_length(0)
                {}

                /**
                 * @brief Copy constructor
                 * */
                Halfedge(const Halfedge &halfedge)
                {
                    m_origin = halfedge.m_origin;
                    m_end = halfedge.m_end;
                    m_length = halfedge.m_length;
                }

                /**
                 * @brief Assignment operator
                 * */
                Halfedge &operator=(const Halfedge &rhs)
                {
                    m_origin = rhs.m_origin;
                    m_end = rhs.m_end;
                    m_length = rhs.m_length;
                    return *this;
                }

                /**
                 * @brief Constructor
                 * This constructor sets the origin- and end-point indices as well as length to preferred values.

                 * @param[in] origin #index of origin vertex
                 * @param[in] end #index of end vertex
                 * @param[in] length #length of the halfedge
                 * */
                Halfedge(const std::size_t origin, const std::size_t end, const double length)
                {
                    if (length < -1e-8)
                        std::cerr << "Warning: [" << __PRETTY_FUNCTION__
                                  << "] Origin and end must be indices > 0 and length should be positiv or 0. One of the values is not correct:"
                                  << std::endl << "origin: " << origin << std::endl << "end: " << end << std::endl
                                  << "length: "
                                  << length;
                    m_origin = origin;
                    m_end = end;
                    m_length = length;
                }

                /**
               * @brief Comparison equal operator
               * The comarison operator equals
               *  TRUE if origin and end point indices equal and the length difference of the two halfedges is < 1e-8
               *  FALSE otherwise
               * */
                bool operator==(const Halfedge &rhs) const
                {
                    return (m_origin == rhs.m_origin && m_end == rhs.m_end && (m_length - rhs.m_length) < 1e-8);
                }

                /**
                 * @brief Comparison not equal operator
                 * */
                bool operator!=(const Halfedge &rhs) const
                {
                    return !((*this) == rhs);
                }

                /**
                 * @brief Get origin vertex index
                 * @return index of origin vertex
                 * */
                std::size_t getOrigin() const
                {
                    return m_origin;
                }

                /**
                 * @brief Get end vertex index
                 * @return index of end vertex
                 * */
                std::size_t getEnd() const
                {
                    return m_end;
                }

                /**
                 * @brief Get length of halfedge
                 * @return length of halfedge
                 * */
                double getLength() const
                {
                    return m_length;
                }

                /**
                 * @brief Tells if halfedge can be added at end.
                 * @param[in] nextHalfedge #halfedge which we want to know about if can be added at end
                 * @return TRUE if halfedge can be added at end and FALSE otherwise.
                 * */
                bool isPrev(const Halfedge &nextHalfedge) const
                {
                    return (m_end == nextHalfedge.m_origin);
                }

                /**
                 * @brief Tells if halfedge can be added at beginning.
                 * @param[in] previousHalfedge #halfedge which we want to know about if can be added at beginning
                 * @return TRUE if halfedge can be added at beginning and FALSE otherwise.
                 * */
                bool isNext(const Halfedge &previousHalfedge) const
                {
                    return (m_origin == previousHalfedge.m_end);
                }

                /**
                 * @brief Tells if halfedge is twin.
                 * @param[in] halfedge #halfedge which we want to know about if it is a twin
                 * @return TRUE if halfedge is a twin and FALSE otherwise.
                 * */
                bool isTwin(const Halfedge &halfedge) const
                {
                    return (m_origin == halfedge.m_end && m_end == halfedge.m_origin
                        && (m_length - halfedge.m_length) < 1e-8);
                }

            private:
                std::size_t m_origin; ///< index of origin vertex
                std::size_t m_end; ///< index of end vertex
                double m_length; ///< length of halfedge
            };

            /**
             * @brief Default constructor
             * */
            Chain()
            {}

            /**
             * @brief Copy constructor
             * */
            Chain(const Chain &chain)
            {
                m_chainedHalfedges = chain.m_chainedHalfedges;
            }

            /**
             * @brief Assignment operator
             * */
            Chain &operator=(const Chain &rhs)
            {
                m_chainedHalfedges = rhs.m_chainedHalfedges;
                return *this;
            }

            /**
             * @brief Tells whether chain is empty or not
             *
             * The chain is empty if there are no halfedges stored in the list m_chainedHalfedges yet. This is tested with size() operator of lists.
             * It is returned
             *  TRUE if chain is empty.
             *  FALSE otherwise.
             *
             * @return bool value
             * */
            bool isEmpty() const
            {
                return m_chainedHalfedges.empty();
            }

            /**
             * @brief Tells whether chain is closed or not
             *
             * The chain is closed if the origin point of the first halfedge and the end point of the last halfedge equal each other.
             * It is returned
             *  TRUE if the chain is closed or empty.
             *  FALSE otherwise.
             * In case the chain is empty, e. g. there are no halfedges stored yet, a warning is printed.
             *
             * @return bool value
             * */
            bool isClosed() const
            {
                if (m_chainedHalfedges.empty())
                {
                    std::cerr << "Warning: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                              << std::endl;
                    return true;
                }
                return (m_chainedHalfedges.front().getOrigin() == m_chainedHalfedges.back().getEnd());
            }

            /**
             * @brief Get number of vertices
             *
             * The number of vertices equals the number of halfedges if the chain is closed.
             * Otherwise the number of halfedges has to be increased by 1 to obtain the number of vertices.
             *
             * @return number of vertices
             * */
            std::size_t getNumberOfVertices() const
            {
                if (this->isClosed())
                    return m_chainedHalfedges.size();
                else
                    return m_chainedHalfedges.size() + 1;
            }

            /**
             * @brief Get length of the chain.
             *
             * The length of the chain is obtained by adding the lengths of all halfedges.
             *
             * @return length
             * */
            double getLength() const
            {
                double length = 0;
                for (typename std::list<Halfedge>::const_iterator it = m_chainedHalfedges.begin();
                     it != m_chainedHalfedges.end(); ++it)
                {
                    length += it->getLength();
                }
                return length;
            }

            /**
             * @brief Get vector of halfedge lengths
             *
             * The list of halfedges is traversed and all halfedge lengths are stored in a vector, s. t. the order is maintained.
             * E. g. returnVector[i] stores the length of the (i+1)-th halfedge of the chain.
             *
             * @return vector of halfedge lengths
             * */
            const std::vector<double> getHalfedgeLengths() const
            {
                if (this->isEmpty())
                {
                    std::cerr << "Warning: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                              << std::endl;
                }
                std::vector<double> lengths;
                for (typename std::list<Halfedge>::const_iterator it = m_chainedHalfedges.begin();
                     it != m_chainedHalfedges.end(); ++it)
                {
                    lengths.push_back(it->getLength());
                }
                return lengths;
            }

            /**
             * @brief Get first halfedge
             *
             * The first halfedge stored in the list m_chainedHalfedges is returned.
             * If the chain is empty an error message is printed.
             *
             * @return first halfedge
             * */
            const Halfedge &getFirstHalfedge() const
            {
                if (this->isEmpty())
                {
                    std::cerr << "Error: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                              << std::endl;
                }
                return m_chainedHalfedges.front();
            }

            /**
             * @brief Get last halfedge
             *
             * The last halfedge stored in the list m_chainedHalfedges is returned.
             * If the chain is empty an error message is printed.
             *
             * @return last halfedge
             * */
            const Halfedge &getLastHalfedge() const
            {
                if (this->isEmpty())
                {
                    std::cerr << "Error: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                              << std::endl;
                }
                return m_chainedHalfedges.back();
            }

            /**
             * @brief Get list of vertex indices
             *
             * The list m_chainedHalfedges is traversed and every vertex is stored in a list once.
             * This is done maintaining the order of the vertices.
             * Firstly every origin vertex index is stored and in case the list is not closed, the end vertex index of the last halfedge is stored, too.
             *
             * If the list m_chainedHalfedges is still empty, an empty list is returned and a warning is printed.
             *
             * @return list of vertex indices
             * */
            const std::list<std::size_t> getVertexIndices() const
            {
                if (this->isEmpty())
                {
                    std::cerr << "Warning: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                              << std::endl;
                }
                std::list<std::size_t> vertexIndices;
                for (typename std::list<Halfedge>::const_iterator it = m_chainedHalfedges.begin();
                     it != m_chainedHalfedges.end(); ++it)
                {
                    vertexIndices.push_back(it->getOrigin());
                }
                if (m_chainedHalfedges.size() == 1 || !(this->isClosed()))
                {
                    vertexIndices.push_back(m_chainedHalfedges.back().getEnd());
                }
                return vertexIndices;
            }

            /**
             * @brief Get shortest distance between vertices
             *
             * The shortest distance between two vertices i and j on a closed chain is calculated and returned.
             *
             * Just in case the chain is closed, it is checked wheter the distance between the two given vertices in a clockwise or counterclockwise direction is shorter.
             * Otherwise the distance from smaller to greater number of vertex is returned automatically.
             *
             * If the chain is empty a warning is printed and 0 is returned.
             * If either of the vertex indices is > number of vertices in the chain or < 1 an error message is printed and 0 is returned, too.
             *
             * @param[in] i #number of the first vertex in the chain
             * @param[in] j #number of the second vertex in the chain
             * @return (shortest) distance between vertices
             * */
            double getShortestDistanceBetween(std::size_t i, std::size_t j) const
            {
                if (this->isEmpty())
                {
                    std::cout << "Warning: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                              << std::endl;
                    return 0;
                }
                if (i > this->getNumberOfVertices() || j > this->getNumberOfVertices() || i < 1 || j < 1)
                {
                    std::cout << "Error: [" << __PRETTY_FUNCTION__ << "] FirstIndex: " << i << " and second index: "
                              << j
                              << " must be positiv integers smaller than the number of points of the chain, which is "
                              << this->getNumberOfVertices() << "." << std::endl;
                    return 0;
                }
                if (i > j)
                    std::swap(i, j);//myFunctions::orderIntegers(i, j);
                double distance = 0;
                std::vector<double> l = this->getHalfedgeLengths();
                for (std::size_t z = i - 1; z < j - 1; z++)
                {
                    distance += l[z];
                }
                if (this->isClosed() && (this->getLength() - distance < distance - 1e-8))
                {
                    distance = this->getLength() - distance;
                }
                return distance;
            }

            /**
             * @brief Get distance between vertices
             *
             * The distance between two vertices i and j in a particular direction is calculated.
             * If the chain is empty a warning is printed and 0 is returned.
             * If either of the vertices is > number of vertices in the chain or < 1 an error message is printed and 0 is returned, too.
             *
             * If i < j the chain simply is traversed from i-th to j-th vertex and halfedge lengths are summed up.
             * Otherwise the remaining length of the chain is returned, provided that the chain is closed.
             * If i > j and chain is not closed, the distance from i to j is returned and a warning is printed.
             *
             * @param[in] i #number of the vertex in the chain, where the length should start, e. g. i=1 for first vertex
             * @param[in] j #number of the vertex in the chain, where the length should end, e. g. j=4 for fourth vertex
             *
             * @return distance between vertices, e.g. between first and fourth chain vertex
             **/
            double getDistanceBetween(std::size_t i, std::size_t j) const
            {
                if (this->isEmpty())
                {
                    std::cout << "Warning: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                              << std::endl;
                    return 0;
                }
                if (i > this->getNumberOfVertices() || j > this->getNumberOfVertices() || i < 1 || j < 1)
                {
                    std::cout << "Error: [" << __PRETTY_FUNCTION__ << "] FirstIndex: " << i << " and second index: "
                              << j
                              << " must be positiv integers smaller than the number of points of the chain, which is "
                              << this->getNumberOfVertices() << "." << std::endl;
                    return 0;
                }
                bool ordered = (i < j);
                if (!ordered)
                    std::swap(i, j);//myFunctions::orderIntegers(i, j);
                double distance = 0;
                std::vector<double> l = this->getHalfedgeLengths();
                for (std::size_t z = i - 1; z < j - 1; z++)
                {
                    distance += l[z];
                }
                if (this->isClosed() && !ordered)
                {
                    distance = this->getLength() - distance;
                }
                else if (!ordered)
                    std::cout << "Warning: [" << __PRETTY_FUNCTION__
                              << "] The chain is supposed to be closed in case the input is not ordered." << std::endl;
                return distance;
            }

            /**
             * @brief Tells if vertex is contained in chain
             *
             * If the chain is empty, a warning is printed and FALSE is returned.
             *
             * The list m_chainedHalfedges is traversed and every origin vertex index of the halfedges is checked for equality with index.
             * If the input vertex index is found TRUE is returned, FALSE otherwise.
             * For not closed chains, the end vertex index for the last halfedge is checked, too.
             *
             * @param[in] vertexIndex vertex index of the searched point
             * @return bool value
             **/
            bool isVertexContained(const std::size_t &vertexIndex) const
            {
                if (this->isEmpty())
                {
                    std::cerr << "Warning: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                              << std::endl;
                    return false;
                }
                for (typename std::list<Halfedge>::const_iterator it = m_chainedHalfedges.begin();
                     it != m_chainedHalfedges.end(); ++it)
                {
                    if (it->getOrigin() == vertexIndex)
                        return true;
                }
                if (!(this->isClosed()) && this->getLastHalfedge().getEnd() == vertexIndex) //here ! was added
                    return true;
                return false;
            }

            /**
             * @brief Tells if halfedge is appendable at beginning
             *
             * For empty chains TRUE is returned.
             * Otherwise the first halfedge in the list m_chainedHalfedges is tested with isNext(previousHalfedge) and obtained bool value is returned.
             *
             * @param[in] previousHalfedge halfedge which is tested to append at beginning
             * @return bool value
             * */
            bool isAppendableAsPrev(const Halfedge &previousHalfedge) const
            {
                if (m_chainedHalfedges.empty())
                    return true;
                else
                    return m_chainedHalfedges.front().isNext(previousHalfedge);
            }

            /**
             * @brief Tells if halfedge is appendable at end
             *
             * For empty chains TRUE is returned.
             * Otherwise the last halfedge in the list m_chainedHalfedges is tested with isPrev(nextHalfedge) and obtained bool value is returned.
             *
             * @param[in] nextHalfedge halfedge which is tested to append at end
             * @return bool value
             * */
            bool isAppendableAsNext(const Halfedge &nextHalfedge) const
            {
                if (m_chainedHalfedges.empty())
                    return true;
                else
                    return m_chainedHalfedges.back().isPrev(nextHalfedge);
            }

            /**
             * @brief Appends halfedge at beginning of chain if possible
             *
             * The method tests whether halfedge is appendable at the beginning using isAppendableAsPrev(prevHalfedge) and appends the halfedge if possible.
             * Otherwise a warning is printed.
             *
             * @param[in] prevHalfedge halfedge that should be appended at beginning
             * @param[out] m_chainedHalfedges chain with appended prevHalfedge
             * */
            void appendPrevHalfedge(const Halfedge &prevHalfedge)
            {
                if (!this->isAppendableAsPrev(prevHalfedge))
                {
                    std::cout << "Warning: [" << __PRETTY_FUNCTION__
                              << "] This halfedge is not appendable at the beginning." << std::endl;
                    std::cout << "The first halfedge of the chain has origin " << this->getFirstHalfedge().getOrigin()
                              << " and prevHalfedge has the end " << prevHalfedge.getEnd() << "." << std::endl;
                }
                else
                {
                    m_chainedHalfedges.push_front(prevHalfedge);
                }
            }

            /**
             * @brief Appends halfedge at end of chain if possible
             *
             * The method tests whether halfedge is appendable at the end using isAppendableAsNext(nextHalfedge) and appends the halfedge if possible.
             * Otherwise a warning is printed.
             *
             * @param[in] nextHalfedge halfedge that should be appended at end
             * @param[out] m_chainedHalfedges chain with appended nextHalfedge
             * */
            void appendNextHalfedge(const Halfedge &nextHalfedge)
            {
                if (!isAppendableAsNext(nextHalfedge))
                {
                    std::cout << "Warning: [" << __PRETTY_FUNCTION__ << "] This halfedge is not appendable at the end."
                              << std::endl;
                    std::cout << "The last halfedge of the chain has end " << this->getLastHalfedge().getEnd()
                              << " and nextHalfedge has the origin " << nextHalfedge.getOrigin() << "." << std::endl;
                }
                else
                {
                    m_chainedHalfedges.push_back(nextHalfedge);
                }
            }

        private:
            std::list<Halfedge> m_chainedHalfedges; ///< list of halfedges - todo? maps<int,maps<int,double>>

        };

        /**
         * @brief Default constructor
         **/
        Boundary()
        {}

        /**
         * @brief Copy constructor
         **/
        Boundary(const Boundary &boundary)
        {
            m_boundary = boundary.m_boundary;
        }

        /**
         * @brief Assignment operator
         **/
        Boundary &operator=(const Boundary &rhs)
        {
            m_boundary = rhs.m_boundary;
            return *this;
        }

        /**
         * @brief Constructor
         *
         * Boundary is constructed from given halfedges.
         * First all halfedges that do not have a twin halfedge contained in the input vector are found.
         * Then the first halfedge is added to the chain.
         * Step by step the algorithm tries to append a halfedge, non-fitting ones are stored in a queue temporary, until a suitable halfedge is found.
         * All halfedges are put back together and this procedure repeats until the chain is closed AND all halfedges are appended.
         *
         * Otherwise a error message is printed, input is not suitable.
         *
         * @param[in] halfedges vector of halfedges of the triangle mesh
         **/
        Boundary(const std::vector<typename Chain::Halfedge> &halfedges)
        {
            std::list<typename Chain::Halfedge> unsortedNonTwinHalfedges = findNonTwinHalfedges(halfedges);
            m_boundary.appendNextHalfedge(unsortedNonTwinHalfedges.front());
            unsortedNonTwinHalfedges.pop_front();
            std::queue<typename Chain::Halfedge> nonFittingHalfedges;
            while (!unsortedNonTwinHalfedges.empty())
            {
                if (m_boundary.isAppendableAsNext(unsortedNonTwinHalfedges.front()))
                {
                    m_boundary.appendNextHalfedge(unsortedNonTwinHalfedges.front());
                    unsortedNonTwinHalfedges.pop_front();
                    while (!nonFittingHalfedges.empty())
                    {
                        unsortedNonTwinHalfedges.push_back(nonFittingHalfedges.front());
                        nonFittingHalfedges.pop();
                    }
                }
                else if (m_boundary.isAppendableAsPrev(unsortedNonTwinHalfedges.front()))
                {
                    m_boundary.appendPrevHalfedge(unsortedNonTwinHalfedges.front());
                    unsortedNonTwinHalfedges.pop_front();
                    while (!nonFittingHalfedges.empty())
                    {
                        unsortedNonTwinHalfedges.push_back(nonFittingHalfedges.front());
                        nonFittingHalfedges.pop();
                    }
                }
                else
                {
                    nonFittingHalfedges.push(unsortedNonTwinHalfedges.front());
                    unsortedNonTwinHalfedges.pop_front();
                }
            }
            if (!m_boundary.isClosed())
                std::cout << "Warning: [" << __PRETTY_FUNCTION__
                          << "] Boundary is not closed although it should be. End points are: " << std::endl
                          << m_boundary.getFirstHalfedge().getOrigin() << std::endl << " and "
                          << m_boundary.getLastHalfedge().getEnd() << std::endl;
        }

        /**
         * @brief Get number of vertices
         *
         * This getter function returns the number of vertices of the triangle mesh.
         *
         * @return number of vertices
         */
        std::size_t getNumberOfVertices() const
        {
            return m_boundary.getNumberOfVertices();
        }

        /**
         * @brief Get length
         *
         * This function returns the length of the boundary, e.g. the sum of all the halfedge lengths.
         *
         * @return length
         */
        double getLength() const
        {
            return m_boundary.getLength();
        }

        /**
         * @brief Get halfedge lengths
         *
         * This function returns a vector of the ordered halfedge lengths, e. g. vector[i] = length of (i+1)-th halfedge
         *
         * @return vector of halfedges
         */
        const std::vector<double> getHalfedgeLengths() const
        {
            return m_boundary.getHalfedgeLengths();
        }

        /**
         * @brief Get list of vertex indices in the chain.
         *
         * This method returns a list of all counter-clockwise ordered indices of the boundary.
         *
         * @return list of vertex indices
         * */
        const std::list<std::size_t> getVertexIndices() const
        {
            return m_boundary.getVertexIndices();
        }

        /**
         * @brief Get distance between vertices
         *
         * The distance between i-th and j-th vertex of the boundary chain is returned.
         * In case the chain is closed, it is checked wheter the distance in a clockwise or counterclockwise direction is shorter, and this one is returned.
         *
         * @param[in] i #number of the first vertex
         * @param[in] j #number of the second vertex
         *
         * @return (shortest) distance between vertices
         * */
        double getShortestDistanceBetween(const std::size_t &i, const std::size_t &j) const
        {
            return m_boundary.getShortestDistanceBetween(i, j);
        }

        /**
         * @brief Get distance between vertices
         *
         * The distance between i-th and j-th vertex of the boundary chain is returned.
         * In case i > j the distance between j-th and i-th vertex is returned, provided a closed chain.
         * Otherwise a warning is printed
         *
         * @param[in] i #number of the first vertex
         * @param[in] j #number of the second vertex
         *
         * @return (shortest) distance between i-th and j-th vertex
         * */
        double getDistanceBetween(const std::size_t &i, const std::size_t &j) const
        {
            return m_boundary.getDistanceBetween(i, j);
        }

        /**
         * @brief Tells if vertex is contained in boundary chain.
         *
         * @param[in] vertexIndex #index (extern) of the searched point
         *
         * @return TRUE if it is contained and FALSE otherwise
         * */
        bool isVertexContained(const std::size_t &internVertexIndex) const
        {
            return m_boundary.isVertexContained(internVertexIndex);
        }

    private:
        /**
         * @brief Finds halfedges without twin halfedge
         *
         * This private method takes a vector of unordered halfedges and finds the halfedges that do not have a twin halfedge contained in the same vector.
         *
         * @param[in] allHalfedges const std::vector<Halfedge>& - vector of all halfedges in the mesh
         *
         * @return list of non-twin halfedges (boundary halfedges)
         */
        const std::list<typename Chain::Halfedge> findNonTwinHalfedges(const std::vector<typename Chain::Halfedge> &allHalfedges)
        {
            std::queue<typename Chain::Halfedge> queue0;
            std::queue<typename Chain::Halfedge> queue1;
            bool actualQueue = 0;
            std::list<typename Chain::Halfedge> nonTwinHalfedges;
            for (std::size_t i = 0; i < allHalfedges.size(); ++i)
            {
                queue0.push(allHalfedges[i]);
            }
            while (!(queue0.empty() && queue1.empty()))
            {
                if (actualQueue == 0)
                {
                    nonTwinHalfedges.push_back(queue0.front());
                    queue0.pop();
                    while (!queue0.empty())
                    {
                        if (nonTwinHalfedges.back().isTwin(queue0.front()))
                        {
                            queue0.pop();
                            nonTwinHalfedges.pop_back();
                            while (!queue0.empty())
                            {
                                queue1.push(queue0.front());
                                queue0.pop();
                            }
                        }
                        else
                        {
                            queue1.push(queue0.front());
                            queue0.pop();
                        }
                    }
                    actualQueue = 1;
                }
                else if (actualQueue == 1)
                {
                    nonTwinHalfedges.push_back(queue1.front());
                    queue1.pop();
                    while (!queue1.empty())
                    {
                        if (nonTwinHalfedges.back().isTwin(queue1.front()))
                        {
                            queue1.pop();
                            nonTwinHalfedges.pop_back();
                            while (!queue1.empty())
                            {
                                queue0.push(queue1.front());
                                queue1.pop();
                            }
                        }
                        else
                        {
                            queue0.push(queue1.front());
                            queue1.pop();
                        }
                    }
                    actualQueue = 0;
                }
            }
            return nonTwinHalfedges;
        }

        Chain m_boundary; ///< boundary chain
    };

    /**
     * @brief Default constructor
     * The default constructor sets the number of inner vertices to 0.
     * */
    //gsHalfEdgeMesh() : m_n(0) { }

    /**
     * @brief Constructor
     * This constructor uses the STL-file named 'filename' for construction.
     *
     * @param[in] filename const std::string& filename of the STL-file containing the triangle mesh
     * */
    gsHalfEdgeMesh(const gsMesh<> &mesh);

    /**
     * @brief Get number of vertices
     * The number of vertices of the triangle mesh in the STL-file is returned.
     *
     * @return number of vertices
     */
    std::size_t getNumberOfVertices() const;

    /**
     * @brief Get number of triangles
     * The number of triangles of the triangle mesh in the STL-file is returned.
     *
     * @return number of triangles
     */
    std::size_t getNumberOfTriangles() const;

    /**
     * @brief Get number of inner vertices
     * The number of inner vertices of the triangle mesh in the STL-file is returned.
     *
     * @return number of inner vertices
     */
    std::size_t getNumberOfInnerVertices() const;

    /**
     * @brief Get number of boundary vertices
     * The number of boundary vertices of the triangle mesh in the STL-file is returned.
     *
     * @return number of boundary vertices
     */
    std::size_t getNumberOfBoundaryVertices() const;

    /**
     * @brief Get vertex
     * The vertex with index 'vertexIndex' is returned.
     *
     * @param[in] vertexIndex const int - index of the vertex that should be returned
     * @return vertex with index 'vertexIndex'
     */
    const gsMesh<>::gsVertexHandle &getVertex(const std::size_t vertexIndex) const;

    /**
     * @brief Get vertex index
     * The vertex index of a three-dimensional vertex from the triangle mesh is returned.
     *
     * @param[in] vertex const ESS_IO::IO_Vertex& - three-dimensional vertex from the triangle mesh
     */
    std::size_t getVertexIndex(const gsMesh<>::gsVertexHandle &vertex) const;

    /**
     * @brief Get vertex index for firts, second or third vertex of triangle
     * Returns the vertex index of a vertex with local vertex index 1,2 or 3 of the triangle with number 'triangleIndex'.
     *
     * @param[in] localVertexIndex const int - local vertex index (1,2,3) of the vertex
     * @param[in] triangleIndex const int - number of the triangle
     * @return vertex index
     */
    std::size_t getGlobalVertexIndex(const std::size_t localVertexIndex, const std::size_t triangleIndex) const;

    /**
     * @brief Get length of the boundary of the triangle mesh.
     * The length of the boundary of the triangle mesh is returned.
     *
     * @return length of boundary
     */
    double getBoundaryLength() const;

    /**
     * @brief Get boundary part lengths between corners
     * A vector containing the numbers of the boundary corners serves as input like [3,10,34,59] for the 3., 10., 34. and 59. boundary vertex being the corners.
     * The lengths of the parts between the corners, e. g. between 3. and 10. boundary vertex as first part length and so on, are returned in a vector of doubles.
     * In case one of the corner entries is > number of boundary vertices, an error message is printed.
     * Mostly there will be 4 corners, as the parameter domain usually is [0,1]x[0,1], but the input can have more than 4 corners too.
     *
     * @param[in] corners std::vector<int>& - vector of boundary point numbers for corners, [1,2,3,4] stands for 1., 2., 3., 4. boundary vertex serve as corners
     * @return vector with part lengths of the boundary
     */
    std::vector<double> getCornerLengths(/*const*/ std::vector<int> &corners) const;
    //std::vector<double> getBoundaryPartLengths(const std::vector<int>& corners) const;

    /**
     * @brief Get chord lengths of boundary
     * A vector storing the lengths of the halfedges of the boundary is returned.
     * The lengths are ordered in the same way they occur in the boundary, meaning returnVector[i] and returnVector[i+1] hold lengths of consecutive halfedges for all i.
     *
     * @return vector of lengths
     */
    const std::vector<double> getBoundaryChordLengths() const;

    /**
     * @brief Get distance between vertices
     * The distance between i-th and j-th boundary vertex is returned.
     * In is checked wheter the distance in a clockwise or counterclockwise direction is shorter, and this one is returned.
     *
     * @param[in] i #number of the first boundary vertex
     * @param[in] j #number of the second boundary vertex
     * @return (shortest) distance between vertices
     * */
    double getShortestBoundaryDistanceBetween(std::size_t i, std::size_t j) const;

    /**
     * @brief Get halfedge length
     * The length of the halfedge with origin and end vertex is returned.
     * In case one of the input vertex indices is greater than the number of vertices, a error message is printed. Nothing is returned then.
     *
     * @param[in] originVertexIndex const int - vertex index of origin vertex
     * @param[in] endVertexIndex const int - vertex index of origin vertex
     * @return length of halfedge
     */
    double getHalfedgeLength(const std::size_t originVertexIndex, const std::size_t endVertexIndex) const;

    /**
     * @brief Returns queue of all opposite halfedges of vertex
     * The opposite halfedge of a point in a triangle is meant to be the halfedge lying opposite of the point, e. g. the halfedge in the triangle not containing the point.
     * Therefore all halfedges of triangles containing the vertex are stored in the return queue.
     * Usually the function is used for inner points. By using the second optional input bool value to 0 it can be used for boundary points too.
     * A warning is printed if the vertex is a boundary vertex although the optional bool value was not set to 0.
     * An error message is printed if the vertex index > number of vertices.
     *
     * @param[in] vertexIndex const int - vertex index
     * @param[in] innerVertex bool - optional bool value, should be set to 0 if vertex is boundary vertex
     * @return queue of opposite halfedges
     */
    const std::queue<typename Boundary::Chain::Halfedge>
    getOppositeHalfedges(const std::size_t vertexIndex, const bool innerVertex = 1) const;

    /**
     * @brief Returns if a vertex is contained in triangle
     * The integers
     *  0 for not contained
     *  1 for being vertex1
     *  2 for being vertex2
     *  3 for being vertex3
     * in triangle are returned.
     * If vertex index > number of vertices a warning is printed and 0 is returned.
     * If number of triangle > number of all triangles a warning is printed and 0 is returned.
     *
     * @param[in] vertexIndex int - vertex index
     * @param[in] triangleIndex int - number of triangle
     * @return 0,1,2 or 3 for [not contained], [vertex1], [vertex2], [vertex3]
     */
    triangleVertexIndex isTriangleVertex(std::size_t vertexIndex, std::size_t triangleIndex) const;

private:
    /**
     * @brief Tells whether a vertex is a boundary vertex or not.
     * The boundary is constructed and it is tested whether the point with index 'internVertexIndex' stored in m_vertices is a boundary vertex.
     * In case the internVertexIndex > N-1 a warning is printed, because that is not supposed to happen.
     *
     * @param[in] internVertexIndex const int - intern vertex index
     * @return bool value
     */
    bool isBoundaryVertex(const std::size_t internVertexIndex) const;

    /**
     * @brief Get intern index of the vertex stored in disordered vertex vector
     * The index i of the intern vector of vertices, s. t. vertices[i] = vertex, is returned.
     * In case the vertex is not stored in the vector, -1 is returned.
     *
     * @param[in] vertex const ESS_IO::IO_Vertex& - vertex that one wants to know the index of
     * @return int intern index
     */
    std::size_t getInternVertexIndex(const gsMesh<>::gsVertexHandle &vertex) const;

    /**
     * @brief Returns halfedge of triangle
     * For a given triangle the first, second or third halfedge are constructed and returned.
     * The first halfedge is meant to be the halfedge from vertex2 to vertex1,
     * the second one is meant to be the one from vertex3 to vertex2 and
     * the third one is meant to be the one from vertex1 to vertex3.
     * The halfedges are constructed in the opposite direction because the boundary is supposed to have counterclockwise direction, if the triangles are given in clockwise direction.
     * In case the number of the halfedge is not 1, 2 or 3 a warning is printed and the first halfedge is returned.
     *
     * @param[in] triangle const ESS_IO::IO_Triangle& - triangle
     * @param[in] numberOfHalfedge int - number of the desired halfedge
     * @return numberofHalfedge-th halfedge of triangle
     */
    const typename Boundary::Chain::Halfedge
    getInternHalfedge(const gsMesh<>::gsFaceHandle &triangle, std::size_t numberOfHalfedge) const;

    /**
     * @brief Creates ordering for vertices
     * The vertices in m_vertices are not ordered. Therefore two index vectors m_sorting and m_inverseSorting are created.
     * These vectors are supposed to be used like
     *  xi = *(vertices[m_sorting[i-1])
     *  vertexIndex(vertex) = m_inverseSorting[internVertexIndex]
     * For construction of the index vectors all vertices stored in m_vertices are passed through.
     * The intern vertex indices of the non-boundary vertices are stored in the first n m_sorting entries. The other way round the vertex indices 1,...,n are stored in the corresponding m_inverseSorting entries.
     * At last the boundary is traversed and for every vertex the intern vertex index is stored in the last N-n m_sorting entries as well as the vertex index is stored in the corresponding m_inverseSorting entry.
     *
     * @param[out] m_sorting vector storing the intern vertex indices at place of vertex index
     * @param[out] m_inverseSorting vector storing the vertex indices at place of intern vertex indices
     */
    void sortVertices();

    std::vector<typename Boundary::Chain::Halfedge> m_halfedges; ///< vector of halfedges
    Boundary m_boundary; ///< boundary of the mesh
    std::size_t m_n; ///< number of inner vertices in the mesh
    std::vector<std::size_t>
        m_inverseSorting; ///< vector of indices s. t. m_inverseSorting[internVertexIndex] = vertexIndex
    std::vector<std::size_t>
        m_sorting; ///< vector that stores the internVertexIndices s. t. m_sorting[vertexIndex-1] = internVertexIndex

};//class gsHalfEdgeMesh

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHalfEdgeMesh.hpp)
#endif