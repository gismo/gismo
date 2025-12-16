#include <gsDomain/gsTensorDomain.h>
#include <gsDomain/gsIndexSubDomain.h>
#include <numeric>

namespace gismo
{

// Implementation of decompose method
template<class T>
typename gsDomain<T>::Ptr
gsCompositeDomain<T>::decompose(index_t npieces, decompositionStrategy strategy) const
{
    switch (strategy)
    {
        case decompositionStrategy::tensor:
        {
            const index_t npatches = m_domains.size();
            
            if (npatches == 0)
                return memory::make_shared(new gsCompositeDomain<T>());
            
            if (npatches == npieces) 
                return memory::make_shared(new gsCompositeDomain<T>(*this));
            
            if (npatches < npieces) 
            {
                auto result = memory::make_shared(new gsCompositeDomain<T>());
                
                std::vector<index_t> piecesPerPatch(npatches);
                index_t basePieces = npieces / npatches;
                index_t remainder = npieces % npatches;
                
                for (index_t i = 0; i < npatches; ++i) 
                {
                    piecesPerPatch[i] = basePieces + (i < remainder ? 1 : 0);
                }
                
                for (index_t i = 0; i < npatches; ++i) 
                {
                    auto decomposed = m_domains[i]->decompose(piecesPerPatch[i]);
                    
                    if (auto composite = std::dynamic_pointer_cast<gsCompositeDomain<T>>(decomposed))
                    {
                        for (size_t j = 0; j < composite->nPieces(); ++j)
                        {
                            auto sub = composite->subdomain(j);
                            sub->setPatchId(m_domains[i]->patchId());
                            result->addDomain(sub);
                        }
                    } else if (decomposed)
                    {
                        decomposed->setPatchId(m_domains[i]->patchId());
                        result->addDomain(decomposed);
                    }
                }
                
                return result;
            }
            
            // npatches > npieces
            auto final_result = memory::make_shared(new gsCompositeDomain<T>());
            
            if (npieces == 0) return final_result; // Handle npieces = 0 explicitly

            std::vector<typename gsDomain<T>::Ptr> temp_subdomains(npieces);
            for (index_t i = 0; i < npieces; ++i) {
                temp_subdomains[i] = memory::make_shared(new gsCompositeDomain<T>());
            }

            std::vector<size_t> current_sizes(npieces, 0); // Tracks num elements in each temp_subdomain
            
            for (index_t i = 0; i < npatches; ++i) 
            {
                index_t minIdx = 0;
                for (index_t j = 1; j < npieces; ++j) 
                {
                    if (current_sizes[j] < current_sizes[minIdx]) 
                    {
                        minIdx = j;
                    }
                }
                
                // Safely cast and add domain
                std::dynamic_pointer_cast<gsCompositeDomain<T>>(temp_subdomains[minIdx])->addDomain(m_domains[i]);
                current_sizes[minIdx] += m_domains[i]->numElements();
            }

            // Now, add the merged subdomains to the final result
            for (const auto& subdomain_ptr : temp_subdomains) {
                final_result->addDomain(subdomain_ptr);
            }
            return final_result;
        }

                case decompositionStrategy::localOptimalBalancing:
        {
            const size_t npatches = m_domains.size();
            // If requested pieces are less than available patches, use the tensor strategy to merge patches.
            // The tensor strategy handles merging into npieces composite subdomains while preserving elements.
            if (npieces < npatches) {
                return this->decompose(npieces, decompositionStrategy::tensor);
            }

            if (npatches == 0)
                return memory::make_shared(new gsCompositeDomain<T>());

            size_t total_elements = this->numElements();
            if (total_elements == 0)
                return memory::make_shared(new gsCompositeDomain<T>());
            
            auto result = memory::make_shared(new gsCompositeDomain<T>());
            
            std::vector<index_t> pieces_per_patch(npatches);
            size_t assigned_pieces = 0;
            for(size_t i = 0; i < npatches; ++i)
            {
                double ratio = (double)m_domains[i]->numElements() / total_elements;
                pieces_per_patch[i] = std::max((index_t)1, (index_t)round(ratio * npieces));
                assigned_pieces += pieces_per_patch[i];
            }
            
            index_t diff = npieces - assigned_pieces;
            while(diff != 0)
            {
                if (diff > 0)
                {
                    double max_rem = -1;
                    index_t best_patch = -1;
                    for(size_t i = 0; i < npatches; ++i)
                    {
                        double ratio = (double)m_domains[i]->numElements() / total_elements;
                        double rem = (ratio * npieces) - pieces_per_patch[i];
                        if (rem > max_rem) {
                            max_rem = rem;
                            best_patch = i;
                        }
                    }
                    pieces_per_patch[best_patch]++;
                    diff--;
                }
                else
                {
                    double min_rem = 2; // Some value larger than any possible remainder
                    index_t best_patch = -1;
                     for(size_t i = 0; i < npatches; ++i)
                    {
                        if (pieces_per_patch[i] > 1) { // Only decrement if > 1
                            double ratio = (double)m_domains[i]->numElements() / total_elements;
                            double rem = (ratio * npieces) - (pieces_per_patch[i] - 1); // Check effect of decrement
                            if (rem < min_rem) {
                                min_rem = rem;
                                best_patch = i;
                            }
                        }
                    }
                    if (best_patch != -1) {
                        pieces_per_patch[best_patch]--;
                        diff++;
                    } else {
                        break;
                    }
                }
            }


            for(size_t i = 0; i < npatches; ++i)
            {
                const auto& patch_domain = m_domains[i];
                size_t num_elems_in_patch = patch_domain->numElements();
                if (num_elems_in_patch == 0) continue;

                std::vector<index_t> all_elem_indices;
                all_elem_indices.reserve(num_elems_in_patch);
                for(auto it = patch_domain->beginAll(); it != patch_domain->endAll(); ++it)
                {
                    all_elem_indices.push_back(it.localId());
                }

                size_t num_pieces_for_patch = pieces_per_patch[i];
                size_t base_size = num_elems_in_patch / num_pieces_for_patch;
                size_t remainder = num_elems_in_patch % num_pieces_for_patch;
                
                size_t current_pos = 0;
                for(size_t j = 0; j < num_pieces_for_patch; ++j)
                {
                    size_t sub_size = base_size + (j < remainder ? 1 : 0);
                    if (sub_size > 0)
                    {
                        std::vector<index_t> sub_indices(all_elem_indices.begin() + current_pos, 
                                                         all_elem_indices.begin() + current_pos + sub_size);
                        auto sub_domain = memory::make_shared(new gsIndexSubDomain<T>(*patch_domain, std::move(sub_indices)));
                        sub_domain->setPatchId(patch_domain->patchId());
                        result->addDomain(sub_domain);
                        current_pos += sub_size;
                    }
                }
            }
            return result;
        }

        case decompositionStrategy::optimalBalancing:
        {
            const size_t npatches = m_domains.size();
            if (npatches == 0)
                return memory::make_shared(new gsCompositeDomain<T>());

            size_t total_elements = this->numElements();
            if (total_elements == 0)
                return memory::make_shared(new gsCompositeDomain<T>());

            if (npieces > total_elements) {
                npieces = total_elements;
            }

            std::vector<std::pair<index_t, index_t>> all_elements;
            all_elements.reserve(total_elements);
            for(const auto& patch_domain : m_domains)
            {
                for(auto it = patch_domain->beginAll(); it != patch_domain->endAll(); ++it)
                {
                    all_elements.emplace_back(patch_domain->patchId(), it.id());
                }
            }

            size_t base_size = total_elements / npieces;
            size_t remainder = total_elements % npieces;
            auto result = memory::make_shared(new gsCompositeDomain<T>());
            
            size_t current_pos = 0;
            for(size_t i = 0; i < npieces; ++i)
            {
                size_t chunk_size = base_size + (i < remainder ? 1 : 0);
                if (chunk_size == 0) continue;

                std::map<index_t, std::vector<index_t>> patch_to_elements;
                for(size_t j = 0; j < chunk_size; ++j)
                {
                    const auto& elem = all_elements[current_pos + j];
                    patch_to_elements[elem.first].push_back(elem.second);
                }

                auto sub_composite = memory::make_shared(new gsCompositeDomain<T>());
                for(const auto& pair : patch_to_elements)
                {
                    const index_t patch_id = pair.first;
                    const auto& indices = pair.second;

                    Ptr parent_domain = nullptr;
                    for(const auto& d : m_domains)
                        if (d->patchId() == patch_id)
                            parent_domain = d;
                    
                    if (!parent_domain) {
                        GISMO_ERROR("Original patch not found for ID!");
                        return typename gsDomain<T>::Ptr(); // Return null pointer on error
                    }
                    
                    auto index_subdomain = memory::make_shared(new gsIndexSubDomain<T>(*parent_domain, indices));
                    index_subdomain->setPatchId(patch_id);
                    sub_composite->addDomain(index_subdomain);
                }
                result->addDomain(sub_composite);

                current_pos += chunk_size;
            }
            return result;
        }
    }
    
    return typename gsDomain<T>::Ptr();
}

} // namespace gismo
