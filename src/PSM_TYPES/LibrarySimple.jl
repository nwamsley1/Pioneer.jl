struct LibSimplePSM{H,L<:Real} <: PSM
    #H is "high precision"
    #L is "low precision"

    #Ion Count Statistics
    best_rank::UInt8 #Highest ranking predicted framgent that was observed
    topn::UInt8 #How many of the topN predicted fragments were observed. 

    b_count::UInt8
    y_count::UInt8
    isotope_count::UInt8

    #Basic Metrics 
    error::H

    #Spectral Simmilarity
    scribe::L
    city_block::L
    spectral_contrast::L
    matched_ratio::L
    entropy_score::L

    #Non-scores/Labels
    precursor_idx::UInt32
    scan_idx::UInt32
end