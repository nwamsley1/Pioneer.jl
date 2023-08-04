isotopes = sort(collect(zip(getIsotopes(best_psms[:,:sequence], best_psms[:,:precursor_idx], QRoots(6), 6, 2), best_psms[:,:RT])), by = x->last(x))
