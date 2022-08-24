/** @file gsParametrization.hpp

    @brief Provides implementation gsParametrization class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl

*/

#include <gsIO/gsOptionList.h>
#include <gsModeling/gsLineSegment.h>

namespace gismo
{

template<class T>
bool gsParametrization<T>::rangeCheck(const std::vector<index_t> &corners, const size_t minimum, const size_t maximum)
{
    for (std::vector<index_t>::const_iterator it = corners.begin(); it != corners.end(); it++)
    {
        if ((size_t)*it < minimum || (size_t)*it > maximum)
        { return false; }
    }
    return true;
}

template<class T>
gsOptionList gsParametrization<T>::defaultOptions()
{
    gsOptionList opt;
    opt.addInt("boundaryMethod", "boundary methodes: {1:chords, 2:corners, 3:smallest, 4:restrict, 5:opposite, 6:distributed}", 4);
    opt.addInt("parametrizationMethod", "parametrization methods: {1:shape, 2:uniform, 3:distance}", 1);
    std::vector<index_t> corners;
    opt.addMultiInt("corners", "vector for corners", corners);
    opt.addReal("range", "in case of restrict or opposite", 0.1);
    opt.addInt("number", "number of corners, in case of corners", 4);
    opt.addReal("precision", "precision to calculate", 1E-8);
    return opt;
}

template<class T>
gsParametrization<T>::gsParametrization(gsMesh<T> &mesh, const gsOptionList & list) : m_mesh(mesh)
{
    m_options.update(list, gsOptionList::addIfUnknown);
}

template<class T>
void gsParametrization<T>::calculate(const size_t boundaryMethod,
                                     const size_t paraMethod,
                                     const std::vector<index_t> &cornersInput,
                                     const T rangeInput,
                                     const size_t numberInput)
{
    GISMO_ASSERT(boundaryMethod >= 1 && boundaryMethod <= 6,
                 "The boundary method " << boundaryMethod << " is not valid.");
    GISMO_ASSERT(paraMethod >= 1 && paraMethod <= 3, "The parametrization method " << paraMethod << " is not valid.");
    size_t n = m_mesh.getNumberOfInnerVertices();
    size_t N = m_mesh.getNumberOfVertices();
    size_t B = m_mesh.getNumberOfBoundaryVertices();
    Neighbourhood neighbourhood(m_mesh, paraMethod);

    T w = 0;
    std::vector<T> halfedgeLengths = m_mesh.getBoundaryChordLengths();
    std::vector<index_t> corners;
    std::vector<T> lengths;

    switch (boundaryMethod)
    {
        case 1:
            m_parameterPoints.reserve(n + B);
            for (size_t i = 1; i <= n + 1; i++)
            {
                m_parameterPoints.push_back(Point2D(0, 0, i));
            }
            for (size_t i = 0; i < B - 1; i++)
            {
                w += halfedgeLengths[i] * (1. / m_mesh.getBoundaryLength()) * 4;
                m_parameterPoints.push_back(Neighbourhood::findPointOnBoundary(w, n + i + 2));
            }
            break;
        case 2:
            corners = cornersInput;
        case 3:
        case 4:
        case 5:
        case 6: // N
            if (boundaryMethod != 2)
                corners = neighbourhood.getBoundaryCorners(boundaryMethod, rangeInput, numberInput);

            m_parameterPoints.reserve(N);
            for (size_t i = 1; i <= N; i++)
            {
                m_parameterPoints.push_back(Point2D(0, 0, i));
            }

            lengths = m_mesh.getCornerLengths(corners);
            m_parameterPoints[n + corners[0] - 1] = Point2D(0, 0, n + corners[0]);

            for (size_t i = corners[0] + 1; i < corners[0] + B; i++)
            {
                w += halfedgeLengths[(i - 2) % B]
                    / findLengthOfPositionPart(i > B ? i - B : i, B, corners, lengths);
                m_parameterPoints[(n + i - 1) > N - 1 ? n + i - 1 - B : n + i - 1] =
                    Neighbourhood::findPointOnBoundary(w, n + i > N ? n + i - B : n + i);
            }
            break;
        default:
            GISMO_ERROR("boundaryMethod not valid: " << boundaryMethod);
    }

    constructAndSolveEquationSystem(neighbourhood, n, N);
}

template<class T>
void gsParametrization<T>::constructAndSolveEquationSystem(const Neighbourhood &neighbourhood,
                                                           const size_t n,
                                                           const size_t N)
{
    gsMatrix<T> A;
    A.resize(n, n);
    std::vector<T> lambdas;
    gsVector<T> b1(n), b2(n);
    b1.setZero(); b2.setZero();

    for (size_t i = 0; i < n; i++)
    {
        lambdas = neighbourhood.getLambdas(i);
        for (size_t j = 0; j < n; j++)
        {
            A(i, j) = ( i==j ? T(1) : -lambdas[j] );
        }

        for (size_t j = n; j < N; j++)
        {
            b1(i) += (lambdas[j]) * (m_parameterPoints[j][0]);
            b2(i) += (lambdas[j]) * (m_parameterPoints[j][1]);
        }
    }

    gsVector<T> u(n), v(n);
    Eigen::PartialPivLU<typename gsMatrix<T>::Base> LU = A.partialPivLu();
    u = LU.solve(b1);
    v = LU.solve(b2);

    for (size_t i = 0; i < n; i++)
        m_parameterPoints[i] << u(i), v(i);
}

template<class T>
const typename gsParametrization<T>::Point2D &gsParametrization<T>::getParameterPoint(size_t vertexIndex) const
{
    return m_parameterPoints[vertexIndex - 1];
}

template<class T>
gsMatrix<T> gsParametrization<T>::createUVmatrix()
{
    gsMatrix<T> m(2, m_mesh.getNumberOfVertices());
    for (size_t i = 1; i <= m_mesh.getNumberOfVertices(); i++)
    {
        m.col(i - 1) << this->getParameterPoint(i)[0], this->getParameterPoint(i)[1];
    }
    return m;
}

template<class T>
gsMatrix<T> gsParametrization<T>::createXYZmatrix()
{
    gsMatrix<T> m(3, m_mesh.getNumberOfVertices());
    for (size_t i = 1; i <= m_mesh.getNumberOfVertices(); i++)
    {
        m.col(i - 1) << m_mesh.getVertex(i)->x(), m_mesh.getVertex(i)->y(), m_mesh.getVertex(i)->z();
    }
    return m;
}

template<class T>
gsMesh<T> gsParametrization<T>::createFlatMesh()
{
    gsMesh<T> mesh;
    mesh.reserve(3 * m_mesh.getNumberOfTriangles(), m_mesh.getNumberOfTriangles(), 0);
    for (size_t i = 0; i < m_mesh.getNumberOfTriangles(); i++)
    {
        typename gsMesh<T>::VertexHandle v[3];
        for (size_t j = 1; j <= 3; ++j)
        {
            v[j - 1] = mesh.addVertex(getParameterPoint(m_mesh.getGlobalVertexIndex(j, i))[0],
                                      getParameterPoint(m_mesh.getGlobalVertexIndex(j, i))[1]);
        }
        mesh.addFace(v[0], v[1], v[2]);
    }
    return mesh.cleanMesh();
}

template<class T>
gsParametrization<T>& gsParametrization<T>::setOptions(const gsOptionList& list)
{
    m_options.update(list, gsOptionList::ignoreIfUnknwon);
    return *this;
}

template<class T>
gsParametrization<T>& gsParametrization<T>::compute()
{
    calculate(m_options.getInt("boundaryMethod"),
              m_options.getInt("parametrizationMethod"),
              m_options.getMultiInt("corners"),
              m_options.getReal("range"),
              m_options.getInt("number"));

    return *this;
}

template<class T>
T gsParametrization<T>::findLengthOfPositionPart(const size_t position,
                                                      const size_t numberOfPositions,
                                                      const std::vector<index_t> &bounds,
                                                      const std::vector<T> &lengths)
{
    GISMO_UNUSED(numberOfPositions);
    GISMO_ASSERT(1 <= position && position <= numberOfPositions, "The position " << position
                 << " is not a valid input. There are only " << numberOfPositions << " possible positions.");
    GISMO_ASSERT(rangeCheck(bounds, 1, numberOfPositions), "The bounds are not a valid input. They have to be out of the possible positions, which only are "
                 << numberOfPositions << ". ");
    size_t numberOfBounds = bounds.size();
    size_t s = lengths.size();
    if (position > (size_t)bounds[numberOfBounds - 1] || position <= (size_t)bounds[0])
        return lengths[s - 1];
    for (size_t i = 0; i < numberOfBounds; i++)
    {
        if (position - (size_t)bounds[0] + 1 > (size_t)bounds[i] - (size_t)bounds[0] + 1
            && position - (size_t)bounds[0] + 1 <= (size_t)bounds[(i + 1) % numberOfBounds] - (size_t)bounds[0] + 1)
            return lengths[i];
    }
    return 0;
}


//******************************************************************************************
//******************************* nested class Neighbourhood *******************************
//******************************************************************************************
template<class T>
gsParametrization<T>::Neighbourhood::Neighbourhood(const gsHalfEdgeMesh<T> & meshInfo, const size_t parametrizationMethod)  : m_basicInfos(meshInfo)
{
    m_localParametrizations.reserve(meshInfo.getNumberOfInnerVertices());
    for(size_t i=1; i <= meshInfo.getNumberOfInnerVertices(); i++)
    {
        m_localParametrizations.push_back(LocalParametrization(meshInfo, LocalNeighbourhood(meshInfo, i), parametrizationMethod));
    }

    m_localBoundaryNeighbourhoods.reserve(meshInfo.getNumberOfVertices() - meshInfo.getNumberOfInnerVertices());
    for(size_t i=meshInfo.getNumberOfInnerVertices()+1; i<= meshInfo.getNumberOfVertices(); i++)
    {
        m_localBoundaryNeighbourhoods.push_back(LocalNeighbourhood(meshInfo, i,0));
    }
}

template<class T>
const std::vector<T>& gsParametrization<T>::Neighbourhood::getLambdas(const size_t i) const
{
    return m_localParametrizations[i].getLambdas();
}

template<class T>
const std::vector<index_t> gsParametrization<T>::Neighbourhood::getBoundaryCorners(const size_t method, const T range, const size_t number) const
{
    std::vector<std::pair<T , size_t> > angles;
    std::vector<index_t> corners;
    angles.reserve(m_localBoundaryNeighbourhoods.size());
    for(typename std::vector<LocalNeighbourhood>::const_iterator it=m_localBoundaryNeighbourhoods.begin(); it!=m_localBoundaryNeighbourhoods.end(); it++)
    {
        angles.push_back(std::pair<T , size_t>(it->getInnerAngle(), it->getVertexIndex() - m_basicInfos.getNumberOfInnerVertices()));
    }
    std::sort(angles.begin(), angles.end());
    if(method == 3)
    {
        this->takeCornersWithSmallestAngles(4, angles, corners);
        std::sort(corners.begin(), corners.end());
        gsDebug << "According to the method 'smallest inner angles' the following corners were chosen:\n";
        for(std::vector<index_t>::iterator it=corners.begin(); it!=corners.end(); it++)
        {
            gsDebug << (*it) << "\n";
        }
    }
    else if(method == 5)
    {
        searchAreas(range, angles, corners);
        gsDebug << "According to the method 'nearly opposite corners' the following corners were chosen:\n";
        for(size_t i=0; i<corners.size(); i++)
        {
            gsDebug << corners[i] << "\n";
        }
    }
    else if(method == 4)
    {
        bool flag = true;
        corners.reserve(4);
        corners.push_back(angles.front().second);
        angles.erase(angles.begin());
        while(corners.size() < 4)
        {
            flag = true;
            for(std::vector<index_t>::iterator it=corners.begin(); it!=corners.end(); it++)
            {
                if(m_basicInfos.getShortestBoundaryDistanceBetween(angles.front().second, *it) < range*m_basicInfos.getBoundaryLength())
                    flag = false;
            }
            if(flag)
                corners.push_back(angles.front().second);
            angles.erase(angles.begin());
        }
        std::sort(corners.begin(), corners.end());
        for(size_t i=0; i<corners.size(); i++)
        {
            gsDebug << corners[i] << "\n";
        }
    }
    else if(method == 6)
    {
        T oldDifference = 0;
        T newDifference = 0;
        std::vector<index_t> newCorners;
        std::vector<T> lengths;
        angles.erase(angles.begin()+number, angles.end());
        gsDebug << "Angles:\n";
        for(size_t i=0; i<angles.size(); i++)
        {
            gsDebug << angles[i].first << ", " << angles[i].second << "\n";
        }
        newCorners.reserve((angles.size()*(angles.size()-1)*(angles.size()-2)*(angles.size()-3))/6);
        corners.reserve((angles.size()*(angles.size()-1)*(angles.size()-2)*(angles.size()-3))/6);
        for(size_t i=0; i<angles.size(); i++)   // n
        {
            for(size_t j=i+1; j<angles.size(); j++) // * (n-1)/2
            {
                for(size_t k=j+1; k<angles.size(); k++) // * (n-2)/3
                {
                    for(size_t l=k+1; l<angles.size(); l++) // * (n-3)/4
                    {
                        newCorners.push_back(angles[i].second);
                        newCorners.push_back(angles[j].second);
                        newCorners.push_back(angles[k].second);
                        newCorners.push_back(angles[l].second);
                        std::sort(newCorners.begin(), newCorners.end());
                        lengths = m_basicInfos.getCornerLengths(newCorners);
                        std::sort(lengths.begin(), lengths.end());
                        newDifference = math::abs(lengths[0] - lengths[3]);
                        if(oldDifference == 0 || newDifference < oldDifference)
                        {
                            corners.erase(corners.begin(), corners.end());
                            corners.push_back(angles[i].second);
                            corners.push_back(angles[j].second);
                            corners.push_back(angles[k].second);
                            corners.push_back(angles[l].second);
                            std::sort(corners.begin(), corners.end());
                        }
                        newCorners.erase(newCorners.begin(), newCorners.end());
                        oldDifference = newDifference;
                    }
                }
            }
        }
        gsDebug << "According to the method 'evenly distributed corners' the following corners were chosen:\n";
        for(size_t i=0; i<corners.size(); i++)
        {
            gsDebug << corners[i] << "\n";
        }
    }
    return corners;
}

template<class T>
const typename gsParametrization<T>::Point2D gsParametrization<T>::Neighbourhood::findPointOnBoundary(const T w, size_t vertexIndex)
{
    GISMO_ASSERT(0 <= w && w <= 4, "Wrong value for w.");
    if(0 <= w && w <=1)
        return Point2D(w,0, vertexIndex);
    else if(1<w && w<=2)
        return Point2D(1,w-1, vertexIndex);
    else if(2<w && w<=3)
        return Point2D(1-w+2,1, vertexIndex);
    else if(3<w && w<=4)
        return Point2D(0,1-w+3, vertexIndex);
    return Point2D();
}

//*****************************************************************************************************
//*****************************************************************************************************
//*******************THE******INTERN******FUNCTIONS******ARE******NOW******FOLLOWING*******************
//*****************************************************************************************************
//*****************************************************************************************************

template<class T>
void gsParametrization<T>::Neighbourhood::takeCornersWithSmallestAngles(size_t number, std::vector<std::pair<T , size_t> >& sortedAngles, std::vector<index_t>& corners) const
{
    sortedAngles.erase(sortedAngles.begin()+number, sortedAngles.end());

    corners.clear();
    corners.reserve(sortedAngles.size());
    for(typename std::vector<std::pair<T, size_t> >::iterator it=sortedAngles.begin(); it!=sortedAngles.end(); it++)
        corners.push_back(it->second);
}

template<class T>
std::vector<T> gsParametrization<T>::Neighbourhood::midpoints(const size_t numberOfCorners, const T length) const
{
    std::vector<T> midpoints;
    midpoints.reserve(numberOfCorners-1);
    T n = 1./numberOfCorners;
    for(size_t i=1; i<numberOfCorners; i++)
    {
        midpoints.push_back(i*length*n);
    }
    return midpoints;
}

template<class T>
void gsParametrization<T>::Neighbourhood::searchAreas(const T range, std::vector<std::pair<T, size_t> >& sortedAngles, std::vector<index_t>& corners) const
{
    T l = m_basicInfos.getBoundaryLength();
    std::vector<T> h = m_basicInfos.getBoundaryChordLengths();
    this->takeCornersWithSmallestAngles(1,sortedAngles, corners);
    std::vector<std::vector<std::pair<T , size_t> > > areas;
    areas.reserve(3);
    for(size_t i=0; i<3; i++)
    {
        areas.push_back(std::vector<std::pair<T , size_t> >());
    }
    std::vector<T> midpoints = this->midpoints(4, l);

    T walkAlong = 0;
    for(size_t i=0; i<h.size(); i++)
    {
        walkAlong += h[(corners[0]+i-1)%h.size()];
        for(int j = 2; j>=0; j--)
        {
            if(math::abs(walkAlong-midpoints[j]) <= l*range)
            {
                areas[j].push_back(std::pair<T , size_t>(m_localBoundaryNeighbourhoods[(corners[0]+i)%(h.size())].getInnerAngle(), (corners[0]+i)%h.size() + 1));
                break;
            }
        }
    }
    std::sort(areas[0].begin(), areas[0].end());
    std::sort(areas[1].begin(), areas[1].end());
    std::sort(areas[2].begin(), areas[2].end());
    bool smaller = false;
    //corners.reserve(3);
    for(size_t i=0; i<areas[0].size(); i++)
    {
        if(areas[0][i].second > (size_t)corners[0] || areas[0][i].second < (size_t)corners[0])
        {
            corners.push_back(areas[0][i].second);
            if(areas[0][i].second < (size_t)corners[0])
            {
                smaller = true;
            }
            break;
        }
    }
    for(size_t i=0; i<areas[1].size(); i++)
    {
        if(smaller)
        {
            if(areas[1][i].second > (size_t)corners[1] && areas[1][i].second < (size_t)corners[0])
            {
                corners.push_back(areas[1][i].second);
                break;
            }
        }
        else if(areas[1][i].second > (size_t)corners[1] || areas[1][i].second < (size_t)corners[0])
        {
            corners.push_back(areas[1][i].second);
            if(areas[1][i].second < (size_t)corners[0])
            {
                smaller = true;
            }
            break;
        }
    }
    for(size_t i=0; i<areas[2].size(); i++)
    {
        if(smaller)
        {
            if(areas[2][i].second > (size_t)corners[2] && areas[2][i].second < (size_t)corners[0])
            {
                corners.push_back(areas[2][i].second);
                break;
            }
        }
        else if(areas[2][i].second > (size_t)corners[2] || areas[2][i].second < (size_t)corners[0])
        {
            corners.push_back(areas[2][i].second);
            break;
        }
    }
}

//*******************************************************************************************
//**************************** nested class LocalParametrization ****************************
//*******************************************************************************************

template<class T>
gsParametrization<T>::LocalParametrization::LocalParametrization(const gsHalfEdgeMesh<T>& meshInfo, const LocalNeighbourhood& localNeighbourhood, const size_t parametrizationMethod)
{
    m_vertexIndex = localNeighbourhood.getVertexIndex();
    std::list<size_t> indices = localNeighbourhood.getVertexIndicesOfNeighbours();
    size_t d = localNeighbourhood.getNumberOfNeighbours();
    switch (parametrizationMethod)
    {
        case 1:
        {
            std::list<T> angles = localNeighbourhood.getAngles();
            VectorType points;
            T theta = 0;
            T nextAngle = 0;
            for(typename std::list<T>::iterator it = angles.begin(); it!=angles.end(); ++it)
            {
                theta += *it;
            }
            Point2D p(0, 0, 0);
            T length = (*meshInfo.getVertex(indices.front()) - *meshInfo.getVertex(m_vertexIndex)).norm();
            Point2D nextPoint(length, 0, indices.front());
            points.reserve(indices.size());
            points.push_back(nextPoint);
            gsVector<T> actualVector = nextPoint - p;
            indices.pop_front();
            T thetaInv = 1./theta;
            gsVector<T> nextVector;
            while(!indices.empty())
            {
                length = (*meshInfo.getVertex(indices.front()) - *meshInfo.getVertex(m_vertexIndex)).norm();
                //length =  (meshInfo.getVertex(indices.front()) - meshInfo.getVertex(m_vertexIndex) ).norm();
                nextAngle = angles.front()*thetaInv * 2 * EIGEN_PI;
                nextVector = (Eigen::Rotation2D<T>(nextAngle).operator*(actualVector).normalized()*length) + p;
                nextPoint = Point2D(nextVector[0], nextVector[1], indices.front());
                points.push_back(nextPoint);
                actualVector = nextPoint - p;
                angles.pop_front();
                indices.pop_front();
            }
            calculateLambdas(meshInfo.getNumberOfVertices(), points);
        }
            break;
        case 2:
            m_lambdas.reserve(meshInfo.getNumberOfVertices());
            for(size_t j=1; j <= meshInfo.getNumberOfVertices(); j++)
            {
                m_lambdas.push_back(0); // Lambda(m_vertexIndex, j, 0)
            }
            while(!indices.empty())
            {
                m_lambdas[indices.front()-1] += (1./d);
                indices.pop_front();
            }
            break;
        case 3:
        {
            std::list<T> neighbourDistances = localNeighbourhood.getNeighbourDistances();
            T sumOfDistances = 0;
            for(typename std::list<T>::iterator it = neighbourDistances.begin(); it != neighbourDistances.end(); it++)
            {
                sumOfDistances += *it;
            }
            T sumOfDistancesInv = 1./sumOfDistances;
            m_lambdas.reserve(meshInfo.getNumberOfVertices());
            for(size_t j=1; j <= meshInfo.getNumberOfVertices(); j++)
            {
                m_lambdas.push_back(0); //Lambda(m_vertexIndex, j, 0)
            }
            for(typename std::list<T>::iterator it = neighbourDistances.begin(); it != neighbourDistances.end(); it++)
            {
                m_lambdas[indices.front()-1] += ((*it)*sumOfDistancesInv);
                indices.pop_front();
            }
        }
            break;
        default:
            GISMO_ERROR("parametrizationMethod not valid: " << parametrizationMethod);
    }
}

template<class T>
const std::vector<T>& gsParametrization<T>::LocalParametrization::getLambdas() const
{
    return m_lambdas;
}

//*****************************************************************************************************
//*****************************************************************************************************
//*******************THE******INTERN******FUNCTIONS******ARE******NOW******FOLLOWING*******************
//*****************************************************************************************************
//*****************************************************************************************************

template<class T>
void gsParametrization<T>::LocalParametrization::calculateLambdas(const size_t N, VectorType& points)
{
    m_lambdas.reserve(N);
    for(size_t j=1; j <= N; j++)
    {
        m_lambdas.push_back(0); //Lambda(m_vertexIndex, j, 0)
    }
    Point2D p(0, 0, 0);
    size_t d = points.size();
    std::vector<T> my(d, 0);
    size_t l=1;
    size_t steps = 0;
    //size_t checkOption = 0;
    for(typename VectorType::const_iterator it=points.begin(); it != points.end(); it++)
    {
        gsLineSegment<2,T> actualLine(p, *it);
        for(size_t i=1; i < d-1; i++)
        {
            if(l+i == d)
                steps = d - 1;
            else
                steps = (l+i)%d - 1;
            //checkoption is set to another number, in case mu is negativ
            if(actualLine.intersectSegment(*(points.begin()+steps), *(points.begin()+(steps+1)%d)/*, checkOption*/))
            {
                //BarycentricCoordinates b(p, *it, *(points.begin()+steps), *(points.begin()+(steps+1)%d));
                // calculating Barycentric Coordinates
                gsMatrix<T, 3, 3> matrix;
                matrix.topRows(2).col(0) = *it;
                matrix.topRows(2).col(1) = *(points.begin()+steps);
                matrix.topRows(2).col(2) = *(points.begin()+(steps+1)%d);
                matrix.row(2).setOnes();

                gsVector3d<T> vector3d;
                vector3d << p, 1;
                gsVector3d<T> delta = matrix.partialPivLu().solve(vector3d);
                my[l-1] = delta(0);
                my[steps] = delta(1);
                my[(steps + 1)%d] = delta(2);
                break;
            }
        }
        for(size_t k = 1; k <= d; k++)
        {
            m_lambdas[points[k-1].getVertexIndex()-1] += (my[k-1]);
        }
        std::fill(my.begin(), my.end(), 0);
        l++;
    }
    for(typename std::vector<T>::iterator it=m_lambdas.begin(); it != m_lambdas.end(); it++)
    {
        *it /= d;
    }
    for(typename std::vector<T>::iterator it=m_lambdas.begin(); it != m_lambdas.end(); it++)
    {
        if(*it < 0)
            gsInfo << *it << "\n";
    }
}

//*******************************************************************************************
//***************************** nested class LocalNeighbourhood *****************************
//*******************************************************************************************

template<class T>
gsParametrization<T>::LocalNeighbourhood::LocalNeighbourhood(const gsHalfEdgeMesh<T>& meshInfo, const size_t vertexIndex, const bool innerVertex)
{
    GISMO_ASSERT(!((innerVertex && vertexIndex > meshInfo.getNumberOfInnerVertices()) || vertexIndex < 1),
                 "Vertex with index " << vertexIndex << " does either not exist (< 1) or is not an inner vertex (> "
                 << meshInfo.getNumberOfInnerVertices() << ").");

    gsVertex<T> tmp;
    m_vertexIndex = vertexIndex;
    std::queue<typename gsHalfEdgeMesh<T>::Halfedge>
        allHalfedges = meshInfo.getOppositeHalfedges(m_vertexIndex, innerVertex);
    std::queue<typename gsHalfEdgeMesh<T>::Halfedge> nonFittingHalfedges;
    m_neighbours.appendNextHalfedge(allHalfedges.front());
    tmp = *meshInfo.getVertex(allHalfedges.front().getOrigin()) - *meshInfo.getVertex(m_vertexIndex);
    m_angles.push_back(tmp.angle((*meshInfo.getVertex(allHalfedges.front().getEnd())
                               - *meshInfo.getVertex(vertexIndex))));
    m_neighbourDistances.push_back(allHalfedges.front().getLength());
    allHalfedges.pop();
    while (!allHalfedges.empty())
    {
        if (m_neighbours.isAppendableAsNext(allHalfedges.front()))
        {
            m_neighbours.appendNextHalfedge(allHalfedges.front());
            tmp = *meshInfo.getVertex(allHalfedges.front().getOrigin()) - *meshInfo.getVertex(m_vertexIndex);
            m_angles
                .push_back(tmp.angle((*meshInfo.getVertex(allHalfedges.front().getEnd())
                                   - *meshInfo.getVertex(m_vertexIndex))));
            m_neighbourDistances.push_back(allHalfedges.front().getLength());
            allHalfedges.pop();
            while (!nonFittingHalfedges.empty())
            {
                allHalfedges.push(nonFittingHalfedges.front());
                nonFittingHalfedges.pop();
            }
        }
        else if (m_neighbours.isAppendableAsPrev(allHalfedges.front()))
        {
            m_neighbours.appendPrevHalfedge(allHalfedges.front());
            tmp = *meshInfo.getVertex(allHalfedges.front().getOrigin()) - *meshInfo.getVertex(m_vertexIndex);
            m_angles
                .push_front(tmp.angle((*meshInfo.getVertex(allHalfedges.front().getEnd())
                                    - *meshInfo.getVertex(m_vertexIndex))));
            m_neighbourDistances.push_back(allHalfedges.front().getLength());
            allHalfedges.pop();
            while (!nonFittingHalfedges.empty())
            {
                allHalfedges.push(nonFittingHalfedges.front());
                nonFittingHalfedges.pop();
            }
        }
        else
        {
            nonFittingHalfedges.push(allHalfedges.front());
            allHalfedges.pop();
        }
    }
}

template<class T>
size_t gsParametrization<T>::LocalNeighbourhood::getVertexIndex() const
{
    return m_vertexIndex;
}

template<class T>
size_t gsParametrization<T>::LocalNeighbourhood::getNumberOfNeighbours() const
{
    return m_neighbours.getNumberOfVertices();
}

template<class T>
const std::list<size_t> gsParametrization<T>::LocalNeighbourhood::getVertexIndicesOfNeighbours() const
{
    return m_neighbours.getVertexIndices();
}

template<class T>
const std::list<T>& gsParametrization<T>::LocalNeighbourhood::getAngles() const
{
    return m_angles;
}

template<class T>
T gsParametrization<T>::LocalNeighbourhood::getInnerAngle() const
{
    T angle = 0;
    for(typename std::list<T>::const_iterator it=m_angles.begin(); it!=m_angles.end(); it++)
    {
        angle += (*it);
    }
    return angle;
}

template<class T>
std::list<T> gsParametrization<T>::LocalNeighbourhood::getNeighbourDistances() const
{
    return m_neighbourDistances;
}

} // namespace gismo
