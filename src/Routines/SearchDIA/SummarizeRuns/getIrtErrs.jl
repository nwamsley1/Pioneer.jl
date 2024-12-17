function get_irt_errs(
    fwhms::Dictionary{Int64, 
                        @NamedTuple{
                            median_fwhm::Float32,
                            mad_fwhm::Float32
                        }},
    prec_to_irt::Dictionary{UInt32, 
    @NamedTuple{best_prob::Float32, 
                best_ms_file_idx::UInt32, 
                best_scan_idx::UInt32, 
                best_irt::Float32, 
                mean_irt::Union{Missing, Float32}, 
                var_irt::Union{Missing, Float32}, 
                n::Union{Missing, UInt16}, 
                mz::Float32}}
    ,
    params::FirstPassSearchParameters
)
    #Get upper bound on peak fwhm. Use median + n*standard_deviation
    #estimate standard deviation by the median absolute deviation. 
    #n is a user-defined paramter. 
    fwhms = map(x->x[:median_fwhm] + params.fwhm_nstd*x[:mad_fwhm],
    fwhms)
    println("fwhms $fwhms")
    #Get variance in irt of apex accross runs. Only consider precursor identified below q-value threshold
    #in more than two runs .
    irt_std = median(
                skipmissing(map(x-> (x[:n] > 2) ? sqrt(x[:var_irt]/(x[:n] - 1)) : missing, prec_to_irt))
                )
    println("irt_std $irt_std")
    #Number of standard deviations to cover 
    irt_std *= params.irt_nstd
    #dictionary maping file name to irt tolerance. 
    return map(x->x+irt_std, fwhms)
end