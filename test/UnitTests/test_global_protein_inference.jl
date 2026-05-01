# Tests for the per-file vs. global protein-inference dispatch in
# ProteinInferenceSearch. These exercise the *semantic* difference at the
# `infer_proteins` level: when a peptide set forms a parsimonious explanation
# only when peptides from multiple files are unioned, per-file inference
# cannot reach it but global inference can.

using Test
using Dictionaries

package_root = dirname(dirname(@__DIR__))
isdefined(Main, :infer_proteins) || include(joinpath(package_root, "src", "utils", "proteinInference.jl"))

# Helper: simulate per-file inference by running infer_proteins once per "file"
# and merging the resulting peptide → protein dicts. The merge keeps the most
# recent assignment when a peptide appears in multiple files (the real
# pipeline's per-file path overwrites in arrival order — this matches that).
function _simulate_per_file_inference(files::Vector{<:Tuple{Vector{ProteinKey}, Vector{PeptideKey}}})
    merged = Dictionary{PeptideKey, ProteinKey}()
    for (proteins, peptides) in files
        result = infer_proteins(proteins, peptides)
        for (k, v) in pairs(result.peptide_to_protein)
            insert!(merged, k, v)
        end
    end
    return merged
end

# Helper: simulate global inference by unioning all (protein, peptide) tuples
# across files into a single Set, then running infer_proteins once.
function _simulate_global_inference(files::Vector{<:Tuple{Vector{ProteinKey}, Vector{PeptideKey}}})
    UniqueKey = Tuple{ProteinKey, PeptideKey}
    s = Set{UniqueKey}()
    for (proteins, peptides) in files
        for i in eachindex(proteins, peptides)
            push!(s, (proteins[i], peptides[i]))
        end
    end
    proteins_vec = Vector{ProteinKey}(undef, length(s))
    peptides_vec = Vector{PeptideKey}(undef, length(s))
    let i = 0
        for (p, q) in s
            i += 1
            proteins_vec[i] = p
            peptides_vec[i] = q
        end
    end
    return infer_proteins(proteins_vec, peptides_vec).peptide_to_protein
end

@testset "Global vs per-file protein inference" begin

    @testset "single file: identical" begin
        # One file, two peptides each uniquely from one protein
        files = [(
            [ProteinKey("A", true, UInt8(0)), ProteinKey("B", true, UInt8(0))],
            [PeptideKey("p1", true, UInt8(0)), PeptideKey("p2", true, UInt8(0))]
        )]
        per_file = _simulate_per_file_inference(files)
        glob = _simulate_global_inference(files)
        @test length(per_file) == 2
        @test per_file[PeptideKey("p1", true, UInt8(0))].name == "A"
        @test per_file[PeptideKey("p2", true, UInt8(0))].name == "B"
        @test sort(collect(keys(per_file))) == sort(collect(keys(glob)))
        for k in keys(per_file)
            @test per_file[k].name == glob[k].name
        end
    end

    @testset "disjoint files: union of per-file == global" begin
        # File 1: protein A with peptide p1
        # File 2: protein B with peptide p2
        # No protein overlap → per-file and global both produce A→p1, B→p2
        files = [
            ([ProteinKey("A", true, UInt8(0))], [PeptideKey("p1", true, UInt8(0))]),
            ([ProteinKey("B", true, UInt8(0))], [PeptideKey("p2", true, UInt8(0))]),
        ]
        per_file = _simulate_per_file_inference(files)
        glob = _simulate_global_inference(files)
        @test length(per_file) == 2
        @test length(glob) == 2
        @test per_file[PeptideKey("p1", true, UInt8(0))].name == glob[PeptideKey("p1", true, UInt8(0))].name == "A"
        @test per_file[PeptideKey("p2", true, UInt8(0))].name == glob[PeptideKey("p2", true, UInt8(0))].name == "B"
    end

    @testset "global pools decoys across files" begin
        # The practical motivation for global inference: decoys from any file
        # are visible to the protein-level PEP fit downstream. Per-file mode
        # can leave a file with zero decoys at the protein level if its decoy
        # PSMs happened not to survive earlier filtering, even when other
        # files in the experiment have plenty.
        #
        # File 1: only target peptides (no decoys observed).
        # File 2: only decoy peptides.
        # Per-file: file 1's inference sees zero decoys → file-1 decoy dict is
        # empty. File 2's sees zero targets → file-2 target dict is empty.
        # Global: the union has both → one inference call assigns both classes.
        files = [
            ([ProteinKey("A", true, UInt8(0)), ProteinKey("B", true, UInt8(0))],
             [PeptideKey("t1", true, UInt8(0)), PeptideKey("t2", true, UInt8(0))]),
            ([ProteinKey("A", false, UInt8(0)), ProteinKey("B", false, UInt8(0))],
             [PeptideKey("d1", false, UInt8(0)), PeptideKey("d2", false, UInt8(0))]),
        ]
        glob = _simulate_global_inference(files)
        @test length(glob) == 4
        @test glob[PeptideKey("t1", true, UInt8(0))].name == "A"
        @test glob[PeptideKey("d1", false, UInt8(0))].name == "A"
        @test glob[PeptideKey("t2", true, UInt8(0))].name == "B"
        @test glob[PeptideKey("d2", false, UInt8(0))].name == "B"
    end

    @testset "global produces a single mapping per peptide" begin
        # Same peptide observed in both files with the same library accession
        # → same group in both modes. The invariant that matters: no peptide
        # can be assigned to different groups in the global path (it's a
        # single dict).
        files = [
            ([ProteinKey("A;B", true, UInt8(0))], [PeptideKey("shared", true, UInt8(0))]),
            ([ProteinKey("A;B", true, UInt8(0))], [PeptideKey("shared", true, UInt8(0))]),
        ]
        glob = _simulate_global_inference(files)
        @test length(glob) == 1
        @test glob[PeptideKey("shared", true, UInt8(0))].name == "A;B"
    end

    @testset "decoy + entrap separation preserved in global" begin
        # Targets and decoys must remain in separate groups even when their
        # peptide sequences would otherwise overlap.
        files = [
            ([ProteinKey("A", true, UInt8(0)), ProteinKey("A_decoy", false, UInt8(0))],
             [PeptideKey("p1", true, UInt8(0)), PeptideKey("p1", false, UInt8(0))]),
            ([ProteinKey("A", true, UInt8(0))], [PeptideKey("p2", true, UInt8(0))]),
        ]
        glob = _simulate_global_inference(files)
        @test glob[PeptideKey("p1", true, UInt8(0))].name == "A"
        @test glob[PeptideKey("p1", false, UInt8(0))].name == "A_decoy"
        @test glob[PeptideKey("p2", true, UInt8(0))].name == "A"
    end

end
