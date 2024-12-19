function getFragIndices(base_type::Char, index::UInt8, sequence_length::UInt8)
    #This ion types run from the N-terminus to the right
    if (base_type=='a')|(base_type=='b')|(base_type=='c')
        return one(UInt8), index
    elseif base_type=='p' #precursor ion covers the entire sequence 
        return one(UInt8), sequence_length
    #Other ion types run from the C-terminus to the left 
    else
        return sequence_length - index + one(UInt8), sequence_length
    end
end

"""
    countSulfurs(plain_sequence::AbstractString, 
                    mods_iterator::Base.RegexMatchIterator,
                    mods_to_sufur_diff::Dict{String, Int8})::Int8

Counts the number of Sulfurs in the sequence (assumex M/S are the only Sulfur containing characters). 
Adds or subtracts any sulfurs contributed by the modifications. 

### Input
-'plain_sequence::AbstractString' -- Plaine peptide sequence
-'mods_iterator::Base.RegexMatchIterator' -- Iteratres through sequence modifications that might contribute/subtract sulfurs
-'mods_to_sufur_diff::Dict{String, Int8}' -- maps modification name to the number of sulfurs gained or lost 

### Output
-'::UInt8' Number of Sulfurs in the precursor (if mods_to_sulfur_diff is incorrectly specified there 
might be an error because UInt8 doesn't allow for negative values)

### Examples
julia> countSulfurs(
           "LLDYRTIMHDENK",
           Base.RegexMatchIterator(r"(?<=\\().*?(?=\\))", "1(7,M,Oxidation)(13,K,AnExampleMod)", false),
           Dict("AnExampleMod" => Int8(2))
       )
3
"""
function countSulfurs!(seq_idx_to_sulfur::Vector{UInt8},
                        plain_sequence::AbstractString, 
                        mods_iterator::Base.RegexMatchIterator,
                        mods_to_sulfur_diff::Dict{String, Int8})::Int8
    sulfur_count = zero(UInt8)
    for (idx, aa) in enumerate(plain_sequence)
        has_sulfur = (aa=='C')|(aa=='M')
        seq_idx_to_sulfur[idx] = UInt8(has_sulfur)
        sulfur_count += has_sulfur
    end

    for mod in mods_iterator
        mod_string = getModName(mod.match)
        if haskey(mods_to_sulfur_diff, mod_string)
            n_sulfur = mods_to_sulfur_diff[mod_string]
            seq_idx_to_sulfur[getModIndex(mod.match)] += n_sulfur
            sulfur_count += n_sulfur
        end
    end
    return sulfur_count
end

function fillIsotopeMods!(seq_idx_to_iso_mod_mass::Vector{Float32},
                        mods_iterator::Base.RegexMatchIterator,
                        iso_mod_to_mass::Dict{String, Float32})

    for mod in mods_iterator
        mod_string = getModName(mod.match)
        if haskey(iso_mod_to_mass, mod_string)
            seq_idx_to_iso_mod_mass[getModIndex(mod.match)] += iso_mod_to_mass[mod_string]
        end
    end
    return nothing
end

function applyIsotopeMod(charge::UInt8, seq_idx_to_iso_mod::Vector{Float32}, start_idx::UInt8, stop_idx::UInt8)::Float32
    mod_mass = zero(Float32)
    for i in range(start_idx, stop_idx)
        mod_mass += seq_idx_to_iso_mod[i]
    end
    return mod_mass/charge
end


function parseKoinaFragments(
    precursor_table::Arrow.Table,
    fragment_table::Arrow.Table,
    FragAnnotationType::FragAnnotation,
    ion_annotation_set::Set{String},
    frag_name_to_idx::Dict{String, UInt16},
    precursor_batch_size::Int64,
    immonium_data_path::String,
    out_dir::String,
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32},
    model_type::KoinaModelType)
    
    #Make output paths and clear any previous entry
    fragment_table_path = joinpath(out_dir, "fragments_table.arrow")
    fragment_indices_path = joinpath(out_dir, "prec_to_frag.arrow")
    rm(fragment_table_path, force = true)
    rm(fragment_indices_path, force = true)

    #Get unique ion types and associated data
    immonium_to_sulfur_count =getImmoniumToSulfurDict(immonium_data_path)
    ion_annotation_to_features_dict = Dict{String, PioneerFragAnnotation}()
    FragAnnotationType = typeof(FragAnnotationType)
    for ion in collect(ion_annotation_set)
        frag_annotation = FragAnnotationType(ion)
        ion_annotation_to_features_dict[ion] = parseAnnotation(frag_annotation,
        immonium_to_sulfur_count=immonium_to_sulfur_count)
    end
    #initialize counters
    precursor_idx, frag_idx = one(UInt32), one(UInt64)
    total_precursors = length(precursor_table[1])
    #writ eht output in batches. 
    #`precursor_batch_size` determines the number of precursors in each batch
    pbar = ProgressBar(total=ceil(Int64, total_precursors/precursor_batch_size))
    while precursor_idx <= total_precursors
        #returns the first precursor_idx and frag_idx of the next batch 
        precursor_idx, frag_idx = appendPioneerLibBatch(
            UInt32(precursor_idx), #first precursor in the batch
            UInt64(frag_idx),#first fragment in the batch
            fragment_table_path,
            fragment_indices_path,
            precursor_batch_size,
            precursor_table[:sequence],
            precursor_table[:mods],
            precursor_table[:isotope_mods],
            fragment_table[:intensities],
            fragment_table[:mz],
            fragment_table[:annotation],
            fragment_table[:precursor_idx],
            ion_annotation_to_features_dict,
            frag_name_to_idx,
            mods_to_sulfur_diff,
            iso_mod_to_mass
        )
        update(pbar)
    end
    Arrow.append(
        fragment_indices_path,
        DataFrame((start_idx = UInt64[frag_idx - 1]))
    )
    return ion_annotation_to_features_dict
end

function appendPioneerLibBatch(
    first_prec_idx::UInt32,
    first_frag_idx::UInt64,
    fragment_table_path::String,
    fragment_indices_path::String,
    n_precursors_in_batch::Int64,
    precursor_sequence::AbstractVector{String},
    precursor_structural_mods::AbstractVector{Union{Missing, String}},
    precursor_isotope_mods::AbstractVector{Union{Missing, String}},
    frag_intensities::AbstractVector{Float32},
    frag_mzs::AbstractVector{Float32},
    frag_annotations::AbstractVector{String},
    precursor_idx::AbstractVector{UInt32},
    ion_annotation_to_data_dict::Dict{String, PioneerFragAnnotation},
    frag_name_to_idx::Dict{String, UInt16},
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32})
    #=
    Batch size is determined by a fixed number of precursors or end of the file.
    Get the number of fragments to pre-allocate for the batch of precursors 
    =#
    n = first_frag_idx
    pid = first_prec_idx
    precs_encountered = 1
    while n <= length(precursor_idx)
        if pid != precursor_idx[n]
            if precs_encountered == n_precursors_in_batch
                n = n - 1
                pid = precursor_idx[n]
                break
            end
            pid = precursor_idx[n]
            precs_encountered += 1
        end
        n += 1
    end
    n = min(n, length(precursor_idx))
    n_frags = n - first_frag_idx + 1
    n_precursors = precs_encountered
    prec_frags = Vector{PioneerFrag}(undef, n_frags)
    pid_to_frag_idxs = Vector{UInt64}(undef, n_precursors)
    prec_sulfur_count = Vector{UInt8}(undef, n_precursors)
    prec_mods_iterator = nothing
    seq_idx_to_sulfur = zeros(UInt8, 255)
    seq_idx_to_iso_mod = zeros(Float32, 255)
    #prec_frags = Vector{Vector{PioneerFrag}}(undef, n_precursors)
    last_pid = zero(UInt32)
    precursor_sulfur_count, precursor_length = zero(UInt8), zero(UInt8)
    frag_idx_start, frag_idx_stop = first_frag_idx, first_frag_idx
    batch_pid = 0
    for frag_idx in range(first_frag_idx, n)
        #Prec id for the fragment 
        frag_annotation = frag_annotations[frag_idx]
        pid = precursor_idx[frag_idx]
    
        #New precursor encountered. Update precursor counstants
        if pid != last_pid
            batch_pid += 1
            frag_idx_start, frag_idx_stop = UInt64(frag_idx), UInt64(frag_idx)
            last_pid = pid
            prec_mods_iterator = parseMods(precursor_structural_mods[pid])
            prec_sulfur_count[batch_pid] = countSulfurs!(
                seq_idx_to_sulfur,
                precursor_sequence[pid],
                prec_mods_iterator,
                mods_to_sulfur_diff
            )
            iso_mods_iterator = parseMods(precursor_isotope_mods[pid])
            #reset the isotope mods vector 
            fill!(seq_idx_to_iso_mod, zero(Float32))

            fillIsotopeMods!(
                seq_idx_to_iso_mod,
                iso_mods_iterator,
                iso_mod_to_mass)

            precursor_sulfur_count = prec_sulfur_count[batch_pid]
            precursor_length = UInt8(length(precursor_sequence[pid]))
        end
        #Update fragment index range for the current precursor
        pid_to_frag_idxs[batch_pid] = frag_idx_start#range(frag_idx_start, frag_idx_stop)
        #Get data for the fragment type 
        frag_data = nothing
        frag_name = nothing
        if startswith(frag_annotation, "Int")
            frag_name = "Int/"*parseInternalIon(frag_annotation)
        else
            frag_name = frag_annotation
        end
        frag_data = ion_annotation_to_data_dict[frag_name]
        #Frag mz and the predicted intensity 
        frag_mz, frag_intensity = frag_mzs[frag_idx], frag_intensities[frag_idx]
        #Start and stop indices of the fragment within the peptide sequence 
        start_idx, stop_idx = getFragIndices(frag_data.base_type, frag_data.frag_index, precursor_length)
        
        sulfur_count = zero(UInt8)
        if frag_data.immonium==false
            for i in range(start_idx, stop_idx)
                sulfur_count += seq_idx_to_sulfur[i]
            end
        end
        if frag_data.internal
            start_idx = parse(UInt8, first(match(r"/([0-9]+)", frag_annotation).captures))+one(UInt8)
            internal_ion_seq = (match(r"(?<=/)[A-Za-z]+(?=[^A-Za-z])", frag_annotation).match)
           
            internal_ion_length = length(internal_ion_seq)
            stop_idx = UInt8(start_idx + internal_ion_length - one(UInt8))
        end
        #Count the number of sulfurs in the fragment 
        sulfur_count += frag_data.sulfur_diff
        #Adjust the fragment m/z to account for isotope labeling mods 
        frag_mz += applyIsotopeMod(frag_data.charge, 
                        seq_idx_to_iso_mod,
                        start_idx, stop_idx)

        prec_frags[frag_idx-first_frag_idx+1] = PioneerFrag(
            frag_mz,
            Float16(frag_intensity),
            frag_name_to_idx[frag_name],
            frag_data.base_type == 'y',
            frag_data.base_type == 'b',
            frag_data.base_type == 'p',
            ((frag_data.base_type == 'a') | (frag_data.base_type == 'x')  | (frag_data.base_type == 'c') | (frag_data.base_type == 'z')),
            frag_data.neutral_diff,
            frag_data.frag_index,
            frag_data.charge,
            frag_data.isotope,
            frag_data.internal,
            frag_data.immonium,
            (start_idx, stop_idx),
            min(sulfur_count, precursor_sulfur_count) #Calculate number of sulfurs in the fragment 
        )
        frag_idx_stop += one(UInt64)
    end
    Arrow.append(
        fragment_table_path,
        DataFrame(prec_frags)
    )
    Arrow.append(
        fragment_indices_path,
        DataFrame((start_idx = pid_to_frag_idxs))
    )
    return first_prec_idx + n_precursors_in_batch, first_frag_idx + n_frags
end

function parseKoinaFragments(
    precursor_table::Arrow.Table,
    fragment_table::Arrow.Table,
    FragAnnotationType::FragAnnotation,
    ion_annotation_set::Set{String},
    frag_name_to_idx::Dict{String, UInt16},
    precursor_batch_size::Int64,
    immonium_data_path::String,
    out_dir::String,
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32},
    model_type::SplineCoefficientModel)
    
    #Make output paths and clear any previous entry
    fragment_table_path = joinpath(out_dir, "fragments_table.arrow")
    fragment_indices_path = joinpath(out_dir, "prec_to_frag.arrow")
    rm(fragment_table_path, force = true)
    rm(fragment_indices_path, force = true)

    #Get unique ion types and associated data
    immonium_to_sulfur_count =getImmoniumToSulfurDict(immonium_data_path)
    ion_annotation_to_features_dict = Dict{String, PioneerFragAnnotation}()
    FragAnnotationType = typeof(FragAnnotationType)
    for ion in collect(ion_annotation_set)
        frag_annotation = FragAnnotationType(ion)
        ion_annotation_to_features_dict[ion] = parseAnnotation(frag_annotation,
        immonium_to_sulfur_count=immonium_to_sulfur_count)
    end
    #initialize counters
    precursor_idx, frag_idx = one(UInt32), one(UInt64)
    total_precursors = length(precursor_table[1])
    #writ eht output in batches. 
    #`precursor_batch_size` determines the number of precursors in each batch
    pbar = ProgressBar(total=ceil(Int64, total_precursors/precursor_batch_size))
    while precursor_idx <= total_precursors
        #returns the first precursor_idx and frag_idx of the next batch 
        precursor_idx, frag_idx = appendPioneerLibBatch(
            UInt32(precursor_idx), #first precursor in the batch
            UInt64(frag_idx),#first fragment in the batch
            fragment_table_path,
            fragment_indices_path,
            precursor_batch_size,
            precursor_table[:sequence],
            precursor_table[:mods],
            precursor_table[:isotope_mods],
            fragment_table[:coefficients],
            fragment_table[:intensities],
            fragment_table[:mz],
            fragment_table[:annotation],
            fragment_table[:precursor_idx],
            ion_annotation_to_features_dict,
            frag_name_to_idx,
            mods_to_sulfur_diff,
            iso_mod_to_mass
        )
        update(pbar)
    end
    Arrow.append(
        fragment_indices_path,
        DataFrame((start_idx = UInt64[frag_idx - 1]))
    )
    return ion_annotation_to_features_dict
end

function appendPioneerLibBatch(
    first_prec_idx::UInt32,
    first_frag_idx::UInt64,
    fragment_table_path::String,
    fragment_indices_path::String,
    n_precursors_in_batch::Int64,
    precursor_sequence::AbstractVector{String},
    precursor_structural_mods::AbstractVector{Union{Missing, String}},
    precursor_isotope_mods::AbstractVector{Union{Missing, String}},
    frag_coefficients::AbstractVector{NTuple{N, Float32}},
    frag_intensities::AbstractVector{Float32},
    frag_mzs::AbstractVector{Float32},
    frag_annotations::AbstractVector{String},
    precursor_idx::AbstractVector{UInt32},
    ion_annotation_to_data_dict::Dict{String, PioneerFragAnnotation},
    frag_name_to_idx::Dict{String, UInt16},
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32}) where {N}
    #=
    Batch size is determined by a fixed number of precursors or end of the file.
    Get the number of fragments to pre-allocate for the batch of precursors 
    =#
    n = first_frag_idx
    pid = first_prec_idx
    precs_encountered = 1
    while n <= length(precursor_idx)
        if pid != precursor_idx[n]
            if precs_encountered == n_precursors_in_batch
                n = n - 1
                pid = precursor_idx[n]
                break
            end
            pid = precursor_idx[n]
            precs_encountered += 1
        end
        n += 1
    end
    n = min(n, length(precursor_idx))
    n_frags = n - first_frag_idx + 1
    n_precursors = precs_encountered
    prec_frags = Vector{PioneerSplineFrag{N}}(undef, n_frags)
    pid_to_frag_idxs = Vector{UInt64}(undef, n_precursors)
    prec_sulfur_count = Vector{UInt8}(undef, n_precursors)
    prec_mods_iterator = nothing
    seq_idx_to_sulfur = zeros(UInt8, 255)
    seq_idx_to_iso_mod = zeros(Float32, 255)
    #prec_frags = Vector{Vector{PioneerFrag}}(undef, n_precursors)
    last_pid = zero(UInt32)
    precursor_sulfur_count, precursor_length = zero(UInt8), zero(UInt8)
    frag_idx_start, frag_idx_stop = first_frag_idx, first_frag_idx
    batch_pid = 0
    for frag_idx in range(first_frag_idx, n)
        #Prec id for the fragment 
        frag_annotation = frag_annotations[frag_idx]
        pid = precursor_idx[frag_idx]
    
        #New precursor encountered. Update precursor counstants
        if pid != last_pid
            batch_pid += 1
            frag_idx_start, frag_idx_stop = UInt64(frag_idx), UInt64(frag_idx)
            last_pid = pid
            prec_mods_iterator = parseMods(precursor_structural_mods[pid])
            prec_sulfur_count[batch_pid] = countSulfurs!(
                seq_idx_to_sulfur,
                precursor_sequence[pid],
                prec_mods_iterator,
                mods_to_sulfur_diff
            )
            iso_mods_iterator = parseMods(precursor_isotope_mods[pid])
            #reset the isotope mods vector 
            fill!(seq_idx_to_iso_mod, zero(Float32))

            fillIsotopeMods!(
                seq_idx_to_iso_mod,
                iso_mods_iterator,
                iso_mod_to_mass)

            precursor_sulfur_count = prec_sulfur_count[batch_pid]
            precursor_length = UInt8(length(precursor_sequence[pid]))
        end
        #Update fragment index range for the current precursor
        pid_to_frag_idxs[batch_pid] = frag_idx_start#range(frag_idx_start, frag_idx_stop)
        #Get data for the fragment type 
        frag_data = nothing
        frag_name = nothing
        if startswith(frag_annotation, "Int")
            frag_name = "Int/"*parseInternalIon(frag_annotation)
        else
            frag_name = frag_annotation
        end
        frag_data = ion_annotation_to_data_dict[frag_name]
        #Frag mz and the predicted intensity 
        frag_mz, frag_intensity, frag_coef = frag_mzs[frag_idx], frag_intensities[frag_idx], frag_coefficients[frag_idx]
        #Start and stop indices of the fragment within the peptide sequence 
        start_idx, stop_idx = getFragIndices(frag_data.base_type, frag_data.frag_index, precursor_length)
        
        sulfur_count = zero(UInt8)
        if frag_data.immonium==false
            for i in range(start_idx, stop_idx)
                try
                sulfur_count += seq_idx_to_sulfur[i]
                catch e
                    println("frag_name $frag_name")
                    println("frag_data $frag_data")
                    println("precursor_length $precursor_length")
                    println("sequence ", precursor_sequence[pid])
                    throw(e)
                end
            end
        end
        if frag_data.internal
            start_idx = parse(UInt8, first(match(r"/([0-9]+)", frag_annotation).captures))+one(UInt8)
            internal_ion_seq = (match(r"(?<=/)[A-Za-z]+(?=[^A-Za-z])", frag_annotation).match)
           
            internal_ion_length = length(internal_ion_seq)
            stop_idx = UInt8(start_idx + internal_ion_length - one(UInt8))
        end
        #Count the number of sulfurs in the fragment 
        sulfur_count += frag_data.sulfur_diff
        #Adjust the fragment m/z to account for isotope labeling mods 
        frag_mz += applyIsotopeMod(frag_data.charge, 
                        seq_idx_to_iso_mod,
                        start_idx, stop_idx)

        prec_frags[frag_idx-first_frag_idx+1] = PioneerSplineFrag(
            frag_mz,
            frag_coef,
            Float16(frag_intensity),
            frag_name_to_idx[frag_name],
            frag_data.base_type == 'y',
            frag_data.base_type == 'b',
            frag_data.base_type == 'p',
            ((frag_data.base_type == 'a') | (frag_data.base_type == 'x')  | (frag_data.base_type == 'c') | (frag_data.base_type == 'z')),
            frag_data.neutral_diff,
            frag_data.frag_index,
            frag_data.charge,
            frag_data.isotope,
            frag_data.internal,
            frag_data.immonium,
            (start_idx, stop_idx),
            min(sulfur_count, precursor_sulfur_count) #Calculate number of sulfurs in the fragment 
        )
        frag_idx_stop += one(UInt64)
    end
    Arrow.append(
        fragment_table_path,
        DataFrame(prec_frags)
    )
    Arrow.append(
        fragment_indices_path,
        DataFrame((start_idx = pid_to_frag_idxs))
    )
    return first_prec_idx + n_precursors_in_batch, first_frag_idx + n_frags
end
