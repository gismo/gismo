#pragma once

#include <gsDomain/gsTensorDomain.h>
#include <gsDomain/gsIndexSubDomain.h>
#include <numeric>

namespace gismo
{

template <class T>
gsCompositeDomainIterator<T>::gsCompositeDomainIterator(index_t _id) : Base(_id) { }

template <class T>
gsCompositeDomainIterator<T>::gsCompositeDomainIterator(const domainContainer& _dom)
    : Base(), m_domains(_dom)
{
    GISMO_ASSERT(!m_domains.empty(), "Empty..");
    m_curDomain = m_domains.begin();

    m_offset.reserve(m_domains.size()+1);
    m_offset.push_back(0);
    for( auto & sd : m_domains )
        m_offset.push_back(m_offset.back()+sd->numElements());
    m_curOffset = m_offset.begin();
    m_cur = (*m_curDomain)->beginAll();
    this->m_pside.patch = m_cur.patch();
}

template <class T>
gsCompositeDomainIterator<T>::gsCompositeDomainIterator(const gsCompositeDomainIterator & other) = default;

template <class T>
typename gsCompositeDomainIterator<T>::domainIter gsCompositeDomainIterator<T>::clone() const { return domainIter(new gsCompositeDomainIterator(*this)); }

template <class T>
gsCompositeDomainIterator<T>::~gsCompositeDomainIterator() { }

template <class T>
index_t gsCompositeDomainIterator<T>::patch() const
{
    if (m_curDomain == m_domains.end())
        return -1;
    return m_cur.patch();
}

template <class T>
size_t gsCompositeDomainIterator<T>::id() const
{
    if (m_curDomain == m_domains.end())
        return m_offset.back();
    return *m_curOffset + m_cur.id();
}

template <class T>
index_t gsCompositeDomainIterator<T>::subdomainIndex() const
{
    if (m_curDomain == m_domains.end())
        return -1;
    return std::distance(m_domains.begin(), m_curDomain);
}

template <class T>
size_t gsCompositeDomainIterator<T>::localId() const { return m_cur.localId(); }

template <class T>
void gsCompositeDomainIterator<T>::next()
{
    if (m_curDomain == m_domains.end()) return;
    
    ++m_cur; // Always increment the child iterator first
    
    // Check if the incremented iterator is still within the current child domain
    if (m_cur != (*m_curDomain)->endAll()) // Using comparison with end()
    {
        this->m_pside.patch = m_cur.patch(); // Update patch ID from the current element
    }
    else // m_cur is no longer good, meaning it reached the end of the child domain
    {
        // Advance to the next child domain
        ++m_curOffset;
        ++m_curDomain;
        if (m_curDomain != m_domains.end())
        {
            m_cur = (*m_curDomain)->beginAll(); // Reinitialize m_cur for the new child domain
            this->m_pside.patch = m_cur.patch(); // Update patch ID from the first element of the new child domain
        } else {
            this->m_pside.patch = -1; // End of iteration
        }
    }
}

template <class T>
void gsCompositeDomainIterator<T>::next(index_t increment)
{
    GISMO_ASSERT(increment >= 0, "Negative increment is not supported.");
    if (increment == 0)
        return;
        
    if (m_curDomain == m_domains.end()) return;

    const size_t target_global_pos = (*m_curOffset + m_cur.id()) + increment;

    const size_t total_elements = m_offset.back();
    if (target_global_pos >= total_elements)
    {
        // Move iterator to the end
        m_curDomain = m_domains.end();
        m_curOffset = m_offset.end() - 1;
        this->m_pside.patch = -1; // End of iteration
        return;
    }

    // Find which subdomain the target position falls into
    // std::upper_bound takes a range [first, last)
    // Here, m_offset.end() should be correct as it includes the sentinel total_elements
    auto it_offset = std::upper_bound(m_offset.begin(), m_offset.end(), target_global_pos);
    --it_offset; // Decrement to get the start offset of the correct subdomain

    // Update current domain and offset
    m_curOffset = it_offset;
    m_curDomain = m_domains.begin() + std::distance(m_offset.begin(), it_offset);
    
    // Adjust local position within the current subdomain
    index_t local_pos = target_global_pos - *m_curOffset;
    m_cur = (*m_curDomain)->beginAll(); // Reinitialize m_cur for the new child domain
    m_cur += local_pos; // Advance m_cur to the exact local position

    this->m_pside.patch = m_cur.patch(); // Update patch ID after m_cur is fully advanced
}

template <class T>
gsVector<T> gsCompositeDomainIterator<T>::lowerCorner() const
{ return m_cur.lowerCorner(); }

template <class T>
gsVector<T> gsCompositeDomainIterator<T>::upperCorner() const
{ return m_cur.upperCorner(); }


template<class T>
gsCompositeDomain<T>::gsCompositeDomain(index_t patchId)
    : Base(patchId), m_domains(), m_topology(nullptr)
{
}

template<class T>
gsCompositeDomain<T>::gsCompositeDomain(const gsMultiBasis<T> & multiBasis)
    : Base(-1), m_domains(multiBasis.nPieces()), m_topology(&multiBasis.topology())
{
    for (index_t i = 0; i != multiBasis.nPieces(); ++i)
    {
        m_domains[i] = multiBasis.basis(i).domain();
        m_domains[i]->setPatchId(i);
    }
}

template<class T>
gsCompositeDomain<T>::gsCompositeDomain(const gsMultiPatch<T> & mp)
    : Base(-1), m_domains(mp.nPieces()), m_topology(&mp.topology())
{
    for (index_t i = 0; i != mp.nPieces(); ++i)
    {
        m_domains[i] = mp.patch(i).basis().domain();
        m_domains[i]->setPatchId(i);
    }
}

template<class T>
void gsCompositeDomain<T>::addDomain(Ptr domain)
{
    m_domains.push_back(domain);
}

template<class T>
void gsCompositeDomain<T>::addDomain(Ptr domain, index_t patchID)
{
    domain->setPatchId(patchID);
    m_domains.push_back(domain);
}

template<class T>
void gsCompositeDomain<T>::keepAlive(Ptr domain)
{
    m_dependencies.push_back(domain);
}

template<class T>
typename gsCompositeDomain<T>::Ptr gsCompositeDomain<T>::subdomain(index_t k) const { return m_domains[k]; }

template<class T>
index_t gsCompositeDomain<T>::getGlobalPatchId(index_t k) const { return m_domains[k]->patchId(); }

template<class T>
gsCompositeDomain<T> gsCompositeDomain<T>::patchDomain(index_t patchId) const
{
    gsCompositeDomain<T> result;
    for (const auto& domain : m_domains)
    {
        if (domain->patchId() == patchId)
            result.addDomain(domain);
    }
    return result;
}

template<class T>
size_t gsCompositeDomain<T>::nPieces() const { return m_domains.size(); }

template<class T>
const typename gsCompositeDomain<T>::domainContainer & gsCompositeDomain<T>::subdomains() const { return m_domains;}

template<class T>
typename gsCompositeDomain<T>::iterator gsCompositeDomain<T>::beginAll() const
{
    if (m_domains.empty())
        return iterator(new gsDomainIteratorEnd<T>(0));
    return iterator(new gsCompositeDomainIterator<T>(m_domains));
}

template<class T>
size_t gsCompositeDomain<T>::numElements() const
{
    size_t sz = 0;
    for (const auto & d : m_domains)
        sz += d->numElements();
    return sz;
}

template<class T>
short_t gsCompositeDomain<T>::degree(short_t i) const
{
    GISMO_ASSERT(!m_domains.empty(), "Empty composite domain.");
    short_t result = m_domains[0]->degree(i);
    for (size_t p = 1; p < m_domains.size(); ++p)
        if (m_domains[p]->degree(i) > result )
            result = m_domains[p]->degree(i);
    return result;
}

template<class T>
short_t gsCompositeDomain<T>::dim() const { return m_domains.front()->dim(); }

template<class T>
gsMatrix<T> gsCompositeDomain<T>::boundingBox() const
{
    GISMO_NO_IMPLEMENTATION;
}

template<class T>
gsMesh<T> gsCompositeDomain<T>::mesh() const
{
    GISMO_NO_IMPLEMENTATION;
}

template<class T>
std::ostream &gsCompositeDomain<T>::print(std::ostream &os) const
{
    os << "Domain container with " << m_domains.size() << " domains.";
    return os;
}

template<class T>
const gsBoxTopology & gsCompositeDomain<T>::topology() const { return *m_topology; }

// Implementation of decompose method
template<class T>
typename gsDomain<T>::Ptr
gsCompositeDomain<T>::decompose(size_t npieces) const
{
    auto result = memory::make_shared(new gsCompositeDomain<T>());
    // Get number of patches
    const size_t npatches = m_domains.size();
    if (npatches == 0 || npieces == 0)
        return result;
    // Get the number of elements per patch
    std::vector<size_t> elements_per_patch(m_domains.size(), 0);
    size_t total_elements = 0;
    for (size_t i = 0; i < m_domains.size(); ++i)
    {
        elements_per_patch[i] = m_domains[i]->numElements();
        total_elements += elements_per_patch[i];
    }

    // if the number of pieces is smaller than the number of patches, we merge patches
    if (npatches > npieces)
    {
        gsDebug<<"Merging " << npatches << " patches into " << npieces << " pieces.\n";
        // Karmarkar Karp merging of patches until we reach npieces
        std::vector<std::vector<size_t>> patchIndices_per_piece(npieces);
        std::vector<size_t> piece_sizes(npieces, 0);
        for (size_t i = 0; i < npatches; ++i)
        {
            // Find the piece with the minimum size
            size_t minIdx = 0;
            size_t min_size = std::numeric_limits<size_t>::max();
            for (size_t j = 0; j < npieces; ++j)
            {
                if (piece_sizes[j] < min_size)
                {
                    min_size = piece_sizes[j];
                    minIdx = j;
                }
            }
            // Assign patch i to piece minIdx
            patchIndices_per_piece[minIdx].push_back(i);
            piece_sizes[minIdx] += elements_per_patch[i];
        }

        // Now create the pieces
        for (size_t i = 0; i < npieces; ++i)
        {
            GISMO_ASSERT(!patchIndices_per_piece[i].empty(), "Piece has no patches assigned.");

            typename gsCompositeDomain<T>::Ptr piece_domain = memory::make_shared(new gsCompositeDomain<T>());
            for (size_t patchIdx : patchIndices_per_piece[i])
            {
                dynamic_cast<gsCompositeDomain<T>*>(piece_domain.get())->addDomain(m_domains[patchIdx]);
                // piece_domain->addDomain(m_domains[patchIdx]);
            }
            result->addDomain(piece_domain);
        }
    }
    // If the number of pieces is greater than or equal to the number of patches, we decompose each patch individually
    else // if (npatches <= npieces)
    {
        gsDebug<<"Decomposing " << npatches << " patches into " << npieces << " pieces.\n";
        // Strategy: Decompose each patch into a number of pieces such that every piece has roughly `elementsPerPiece` elements.
        // We use the Karmarkar-Karp (KK) algorithm to distribute the pieces among the patches.
        std::vector<size_t> pieces_per_patch(npatches, 1);
        size_t current_total_pieces = npatches;
        while (current_total_pieces < npieces)
        {
            // Find the patch with the maximum number of elements per assigned piece
            size_t maxIdx = 0;
            real_t max_ratio = -1.0;
            for (size_t i = 0; i < npatches; ++i)
            {
                real_t ratio = (real_t)elements_per_patch[i] / pieces_per_patch[i];
                if (ratio > max_ratio)
                {
                    max_ratio = ratio;
                    maxIdx = i;
                }
            }
            pieces_per_patch[maxIdx]++;
            current_total_pieces++;
        }
        gsDebug<<gsAsVector<size_t>(pieces_per_patch)<<"\n";
        gsDebug<<"Total pieces after KK distribution: " << current_total_pieces << "\n";
        // Sum the pieces_per_patch to verify it matches npieces
        size_t verify_total = std::accumulate(pieces_per_patch.begin(), pieces_per_patch.end(), 0);
        gsDebugVar(verify_total);
        // Now decompose each patch accordingly and collect the results
        for (size_t i = 0; i < npatches; ++i)
        {
            if (pieces_per_patch[i] == 1)
            {
                result->addDomain(m_domains[i]);
            }
            else if (pieces_per_patch[i] > 1)
            {
                typename gsDomain<T>::Ptr decomposed_patch = m_domains[i]->decompose(pieces_per_patch[i]);
                if (dynamic_cast<gsCompositeDomain<T>*>(decomposed_patch.get()) == nullptr)
                {
                    GISMO_ERROR("Decomposition of patch did not return a composite domain.");
                }
                else
                {
                    gsDebugVar(dynamic_cast<gsCompositeDomain<T>*>(decomposed_patch.get())->nPieces());
                    gsDebugVar(pieces_per_patch[i]);

                    gsCompositeDomain<T>* decomposed_composite =
                        dynamic_cast<gsCompositeDomain<T>*>(decomposed_patch.get());
                    for (size_t j = 0; j < decomposed_composite->nPieces(); ++j)
                    {
                        auto sub = decomposed_composite->subdomain(j);
                        sub->setPatchId(m_domains[i]->patchId());
                        result->addDomain(sub);
                    }
                }
            }
            else // pieces_per_patch[i] == 0
            {
                // This should not happen due to the KK algorithm above
                GISMO_ERROR("Invalid number of pieces assigned to patch.");
            }
        }
    }
    return result;
}


} // namespace gismo
