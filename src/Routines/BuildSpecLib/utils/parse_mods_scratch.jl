
sequence = "PEPMTIDME"

# Create the variable modifications vector
var_mods = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}()
push!(var_mods, (p=r"M", r="Unimod:35"))
# Get the matches
matches = matchVarMods(sequence, var_mods)
countVarModCombinations(matches, 2)
var_mods_string = Vector{String}(undef, countVarModCombinations(matches, 2))
fillVarModStrings!(var_mods_string, matches, "", 2)

# Example usage:
sequence = "ILSISADIETIGEILK"
mods = "(1,tag6)(15,tag6)"
mod_masses = zeros(Float32, length(sequence))
mod_to_mass = Dict("tag6" => 229.163f0)

# Strip mods and store masses
fillModMasses!(mod_masses, sequence, mods, mod_to_mass)

# Get mass offset for full peptide
full_mass_offset = getMassOffset(1, length(sequence), mod_masses)

# Get mass offset for a b5 ion (first 5 amino acids)
b5_mass_offset = getMassOffset(1, 5, mod_masses)

# Get mass offset for a y5 ion (last 5 amino acids)
y5_mass_offset = getMassOffset(length(sequence)-4, length(sequence), mod_masses)


# Example usage:
sequence = "ILSISADIETIGEILK"
mods = "(1,tag6)(8,Unimod:35)(15,tag6)"
mod_masses = zeros(Float32, length(sequence))
mod_to_mass = Dict("tag6" => 229.163f0, "Unimod:35" => 15.995f0)

# First fill all modification masses
fillModMasses!(mod_masses, sequence, mods, mod_to_mass)
println("After fillModMasses!: ", mod_masses)

# Then strip only tag6 modifications
stripMods!(mod_masses, sequence, mods, mod_to_mass, Set(["tag6"]))
println("After stripMods!: ", mod_masses)

# Mass offset calculations will now exclude tag6 but include Unimod:35
full_mass_offset = getMassOffset(1, length(sequence), mod_masses)
println("Full mass offset (should only include Unimod:35): ", full_mass_offset)


# Example usage:
# Example usage:
sequence = "I(tag6)LSISADI(Unimod:35)ETIGEILK(tag6)"
mods = ""
mod_masses = zeros(Float32, 255)
mod_to_mass = Dict("tag6" => 229.163f0, "Unimod:35" => 15.995f0)
parseMods!(mods, sequence)
#mods should now be "(1,I,tag6)(8,I,Unimod:35)(15,K,tag6)"
#format is position, amino acid, modification name

# Then strip only tag6 modifications
stripped_sequence, mods_string = stripMods!(mod_masses, sequence, mods, mod_to_mass, Set(["tag6"]))
println("After stripMods!: ", mod_masses[1]) # should be -229.163 (negative because the mod is being stripped)
println("After stripMods!: ", mod_masses[8]) # should be 15.995
println("After stripMods!: ", mod_masses[15]) # should be -229.163 (negative because the mod is being stripped)
println("stripped_seuqnece $stripped_sequence") #should be "ILSISADI(Unimod:35)ETIGEILK"
println("mods_string $mods_string") #should be "(8,I,Unimod:35)"


# Example usage and test
sequence = "I(tag6)LSISADI(Unimod:35)ETIGEILK(tag6)"
mods = ""
mod_masses = zeros(Float32, 255)
mod_to_mass = Dict("tag6" => 229.163f0, "Unimod:35" => 15.995f0)

# First parse the mods
mods = parseMods!(sequence)
println("Parsed mods: ", mods)

# Then strip tag6 modifications
stripped_sequence, remaining_mods = stripMods!(mod_masses, sequence, mod_to_mass, Set(["tag6"]))
println("Stripped sequence: ", stripped_sequence)
println("Remaining mods: ", remaining_mods)
println("Mod masses at positions 1, 8, 15: ", mod_masses[1], ", ", mod_masses[8], ", ", mod_masses[15])



iso_mod_masses = Dict(
    "tag6" => [("d0" = 0.0f0, "d4" = 4.0f0, "d8" = 8.0f0)]
    )

    Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}()
mod_key = "tag6"
mod_channel = @NamedTuple((:channel, :mass), ("d0", 0.0f0))
iso_mod_masses = zeros(Float32, 255)
structural_mods = "(1,I,tag6)(5,L,tag5)(6,M,Unimod:53)(16,K,tag6)"
isotopic_mods = "(2, d0)" # refers to the second mod in "structural" mods. This is the 'd0' channel of tag5
#now add d0 channel for tag6
addIsoMods!(isotopic_mods, structural_mods, mod_key, mod_channel)
println("After addIsoMods! isotopic_mods: ", isotopic_mods) #should be "(1, d0)(2, d0)(4, d0)"

# Test the function
let
    mod_key = "tag6"
    mod_channel = (channel="d0", mass=0.0f0)
    structural_mods = "(1,I,tag6)(5,L,tag5)(6,M,Unimod:53)(16,K,tag6)"
    isotopic_mods = "(2, d0)"  # existing iso mod for tag5
    
    result = addIsoMods!(isotopic_mods, structural_mods, mod_key, mod_channel)
    println("Original isotopic_mods: ", isotopic_mods)
    println("After addIsoMods!: ", result)
end

iso_mods_dict = Dict(
    "tag6" => Dict(
        "d0" => 0.0f0,
        "d4" => 4.0f0,
        "d8" => 8.0f0
    ),
    "tag5" => Dict(
        "d0" => 0.0f0,
        "d4" => 4.0f0,
        "d8" => 8.0f0
    )
)
structural_mods = "(1,I,tag6)(5,L,tag5)(6,M,Unimod:53)(16,K,tag6)"
isotopic_mods = "(1, d0)(2, d4)(4, d8)"
iso_mod_masses = zeros(Float32, 255)
getIsoModMasses!(iso_mod_masses, structural_mods, isotopic_mods, iso_mods_dict)
println("After getIsoModMasses!: ", 
iso_mod_masses[1], ", ", # d0 for tag6 should be 0.0
iso_mod_masses[5], ", ", # d4 for tag5 should be 4.0
iso_mod_masses[16], ", ", # d8 for tag6 should be 8.0
iso_mod_masses[7]) # no iso mod so should be zero )

# Test the function
let
    iso_mods_dict = Dict(
        "tag6" => Dict(
            "d0" => 0.0f0,
            "d4" => 4.0f0,
            "d8" => 8.0f0
        ),
        "tag5" => Dict(
            "d0" => 0.0f0,
            "d4" => 4.0f0,
            "d8" => 8.0f0
        )
    )
    structural_mods = "(1,I,tag6)(5,L,tag5)(6,M,Unimod:53)(16,K,tag6)"
    isotopic_mods = "(1, d0)(5, d4)(16, d8)"
    iso_mod_masses = zeros(Float32, 255)
    
    getIsoModMasses!(iso_mod_masses, structural_mods, isotopic_mods, iso_mods_dict)
    println("Masses at positions:")
    println("pos 1 (tag6, d0): ", iso_mod_masses[1])  # Should be 0.0
    println("pos 5 (tag5, d4): ", iso_mod_masses[5])  # Should be 4.0
    println("pos 16 (tag6, d8): ", iso_mod_masses[16])  # Should be 8.0
    println("pos 6 (no iso mod): ", iso_mod_masses[6])  # Should be 0.0
end

sequence = "PEPTIDE"
structural_mods = "(1,P,Phospho)(4,T,MyFavMod)(7,E,Acetyl)"

sequence, structural_mods = reverseSequence(sequence::String, structural_mods::String)
@assert sequence = "DITPEPE" #sequence is reversed except for the last AA
@assert structural_mods = "(1,E,Phospho)(3,T,MyFavMod)(7,Acetyl)" #The modification on the T gets swapped to the new position but mods on first and last position stay on the same position 

# Test the function
let
    sequence = "PEPTIDE"
    structural_mods = "(1,P,Phospho)(4,T,MyFavMod)(7,E,Acetyl)"
    
    rev_seq, rev_mods = reverseSequence(sequence, structural_mods)
    
    @assert rev_seq == "EDITPEP" "Sequence reversal failed"
    @assert rev_mods == "(1,E,Phospho)(4,T,MyFavMod)(7,P,Acetyl)" "Modification reversal failed"
    println("Reversal tests passed!")
end


sequence = "PEPTIDE"
mods = "(1,n,mymod-nterm)(4,T,MyMod)(7,E,Acetyl)"
rev_seq, rev_mods = reverseSequence(sequence, mods)
# rev_seq  = "EDITPEP"
# rev_mods = "(1,E,Phospho)(4,T,MyMod)(7,P,Acetyl)"

let
    Random.seed!(1844)
    sequence = "PEPTIDE"
    mods = "(1,n,mymod-nterm)(1,P,mymod)(4,T,Phospho)(7,E,Acetyl)(7,c,mymod-cterm)"
    shuffled_seq, shuffled_mods = shuffleSequence(sequence, mods)
    
    # Test that:
    # 1. Seed set so sequence should always be DPTPIEE
    @assert shuffled_seq == "DPTPIEE"
    @assert shuffled_mods == "(1,n,mymod-nterm)(2,P,mymod)(3,T,Phospho)(7,E,Acetyl)(7,c,mymod-cterm)"

end

srows = (test_lib[!,:sequence].=="ETEELHHDR").&(test_lib[!,:isotopic_mods].=="(1, d8)")
test_lib[srows,[:is_decoy,:precursor_idx,:sequence,:prec_charge,:structural_mods,:isotopic_mods,:prec_mz,:frag_mz,:frag_type,:frag_series_number,:library_intensity]]
kl[kl[!,:PeptideSequence].=="ETEELHHDR",[:ModifiedPeptide,:PeptideSequence,:PrecursorCharge,:PrecursorMz,:ProductMz,:FragmentType,:FragmentCharge,:FragmentSeriesNumber]]




test_df = DataFrame(Arrow.Table("/Users/nathanwamsley/temp/test.poin/precursors_table.arrow"))
test_df[test_df[!,:sequence].=="ETEELHHDR",[:precursor_idx,:mz]]
kl[kl[!,:PeptideSequence].=="ETEELHHDR",[:ModifiedPeptide,:PeptideSequence,:PrecursorCharge,:PrecursorMz,:ProductMz,:FragmentType,:FragmentCharge,:FragmentSeriesNumber]]


#=
if abspath(PROGRAM_FILE) == @__FILE__
    using DataFrames
    
    # Test setup
    struct BasicEmpiricalLibrary
        libdf::DataFrame
    end
    
    # Create test data with multiple fragments per precursor
    test_df = DataFrame(
        sequence = repeat(["PEPTIDE", "PROTEIN"], inner=3),
        structural_mods = repeat(["(2,P,Oxidation)", ""], inner=3),
        precursor_idx = repeat([1, 2], inner=3),
        fragment_idx = repeat(1:3, outer=2),
        is_decoy = falses(6)
    )
    test_lib = BasicEmpiricalLibrary(test_df)
    
    # Simple reverseSequence implementation for testing
    function reverseSequence(seq::String, mods::String)
        return reverse(seq), mods  # Simplified for testing
    end
    
    # Generate decoys
    getRevDecoys!(test_lib)
    
    # Verify results
    result_df = test_lib.libdf
    
    # Basic checks
    @assert size(result_df, 1) == 12 "Should have original + decoy rows"
    @assert all(result_df[7:end, :is_decoy]) "New rows should be marked as decoys"
    
    # Check precursor consistency
    for precursor_idx in unique(result_df.precursor_idx)
        decoy_rows = result_df[result_df.is_decoy .& (result_df.precursor_idx .== precursor_idx), :]
        if !isempty(decoy_rows)
            # All fragments of same precursor should have same sequence
            @assert length(unique(decoy_rows.sequence)) == 1 "All fragments of precursor $precursor_idx should have same sequence"
            # Sequence should be reversed
            orig_seq = first(result_df[.!result_df.is_decoy .& (result_df.precursor_idx .== precursor_idx), :sequence])
            @assert first(decoy_rows.sequence) == reverse(orig_seq) "Decoy sequence should be reversed"
        end
    end
    
    println("All tests passed!")
end
=#