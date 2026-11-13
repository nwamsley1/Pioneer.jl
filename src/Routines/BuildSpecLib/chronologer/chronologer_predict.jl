# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# src/chronologer/chronologer_predict.jl

"""
    predict_retention_times(chronologer_in_path, chronologer_out_path;
                            rt_model = DEFAULT_RT_MODEL)

Predict retention times for peptides through Koina, with the model named by
`rt_model` (a key of `RT_MODEL_CONFIGS`).

Parameters:
- chronologer_in_path::String: Path to Arrow file containing peptide data.
                               Must have a 'koina_sequence' column.
- chronologer_out_path::String: Where the same table is written with an `rt`
                                column of predictions.
"""
function predict_retention_times(chronologer_in_path::String, chronologer_out_path::String;
                                 rt_model::String = DEFAULT_RT_MODEL)
    try
        chronologer_table = DataFrame(Tables.columntable(Arrow.Table(chronologer_in_path)))
        predictions = predict_rt_koina(chronologer_table; rt_model = rt_model)
        chronologer_table[!, :rt] = predictions
        Arrow.write(chronologer_out_path, chronologer_table)
        return
    catch e
        @user_warn "Retention time prediction with $(rt_model) failed through Koina." exception=e
        rethrow(e)
    end
end

"""
Helper function to predict RTs using Koina service.
"""
function predict_rt_koina(chronologer_table::DataFrame;
                          rt_model::String = DEFAULT_RT_MODEL)::Vector{Float32}
    haskey(RT_MODEL_CONFIGS, rt_model) || error(
        "Unknown rt_model '$rt_model'. Valid: $(join(sort(collect(keys(RT_MODEL_CONFIGS))), ", "))")
    model = RetentionTimeModel(rt_model)

    # Prepare batches
    batches = prepare_koina_batch(
        model,
        chronologer_table,
        batch_size=1000
    )

    # Make requests
    results = make_koina_batch_requests(
        batches,
        KOINA_URLS[rt_model]
    )
    
    # Parse results
    rt_predictions = Float32[]
    for result in results
        batch_result = parse_koina_batch(model, result)
        append!(rt_predictions, batch_result.fragments.rt)
    end
    
    return rt_predictions
end

"""
    predict_ion_mobility(in_path, out_path, im_model)

Read the precursor table at `in_path`, add `ccs` (Å², from the Koina model
`im_model`) and `inv_ion_mobility` (1/K0, Vs/cm²) columns, and write it to
`out_path`. Written to a new file rather than in place for the same reason
`predict_retention_times` is: the input Arrow file may still be mmap-locked.
"""
function predict_ion_mobility(in_path::String, out_path::String, im_model::String)
    table = DataFrame(Tables.columntable(Arrow.Table(in_path)))
    ccs = predict_ccs_koina(table, im_model)
    table[!, :ccs] = ccs
    table[!, :inv_ion_mobility] =
        ccs_to_inv_ion_mobility.(ccs, table.precursor_charge, table.mz)
    Arrow.write(out_path, table)
    return
end

"""
Helper function to predict CCS values using a Koina ion-mobility model.
"""
function predict_ccs_koina(table::DataFrame, im_model::String)::Vector{Float32}
    model = IonMobilityModel(im_model)
    if hasproperty(table, :sequence) && hasproperty(table, :mods)
        cfg = IM_MODEL_CONFIGS[im_model]
        dropped = Dict{String, Int}()
        seqs = [im_koina_sequence(s, m, cfg, dropped) for (s, m) in zip(table.sequence, table.mods)]
        isempty(dropped) || @user_warn "Ion-mobility model $im_model cannot encode " *
            join(("$k ($v precursors)" for (k, v) in sort!(collect(dropped))), ", ") *
            "; their residues are predicted as unmodified for CCS. Fragments and retention times keep them."
        table = DataFrame(koina_sequence = seqs, precursor_charge = table.precursor_charge)
    end
    batches = prepare_koina_batch(model, table, batch_size=1000)
    results = make_koina_batch_requests(batches, KOINA_URLS[im_model])
    ccs = Float32[]
    for result in results
        append!(ccs, parse_koina_batch(model, result).fragments.ccs)
    end
    length(ccs) == nrow(table) || error(
        "ion-mobility model $im_model returned $(length(ccs)) values for $(nrow(table)) precursors")
    return ccs
end

"""
    im_koina_sequence(sequence, mods, cfg, dropped) -> String

The Koina sequence for an ion-mobility model from a precursor's `sequence` and
`mods` ("(1,n,Unimod:1)(7,M,Unimod:35)"): each modification `cfg.supported_mods`
covers at its site is written as `[UNIMOD:n]`, N-terminal ones as a ProForma
prefix when `cfg.nterm_prefix`; any other is left off, so its residue is
predicted as unmodified, and counted in `dropped` ("Unimod:2062 on C" => n).
"""
function im_koina_sequence(sequence::AbstractString, mods, cfg, dropped::Dict{String, Int})
    ismissing(mods) && return String(sequence)
    prefix = ""
    on_residue = [String[] for _ in 1:length(sequence)]
    for m in parseMods(mods)
        name = getModName(m.match)
        site = only(split(m.match, ',')[2])
        id = unimod_id(name)
        if id === nothing || !occursin(site, get(cfg.supported_mods, id, ""))
            key = "$name on $(_site_label(site))"
            dropped[key] = get(dropped, key, 0) + 1
            continue
        end
        if site == 'n' && cfg.nterm_prefix
            prefix *= "[UNIMOD:$id]"
        else
            push!(on_residue[getModIndex(m.match)], "[UNIMOD:$id]")
        end
    end
    io = IOBuffer()
    isempty(prefix) || print(io, prefix, '-')
    for (i, aa) in enumerate(sequence)
        print(io, aa)
        foreach(t -> print(io, t), on_residue[i])
    end
    return String(take!(io))
end

"""
    ccs_to_inv_ion_mobility(ccs, charge, mz)

Mason–Schamp conversion as used by AlphaPeptDeep for Bruker timsTOF data
(N2 drift gas, 28 Da): `1/K0 = CCS · sqrt(μ) / (z · 1059.62245)`, where
`μ = M·28/(M+28)` is the reduced mass of the ion (`M = mz·z`) and N2.
"""
function ccs_to_inv_ion_mobility(ccs::Real, charge::Integer, mz::Real)::Float32
    M = Float64(mz) * charge
    μ = M * 28.0 / (M + 28.0)
    return Float32(Float64(ccs) * sqrt(μ) / (charge * 1059.62245))
end

