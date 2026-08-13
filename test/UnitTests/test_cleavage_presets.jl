# Cleavage-rule regexes for the enzyme presets, checked against digests worked
# out by hand.
#
# `digest_sequence` ends a peptide at the *start* of a regex match, so the first
# element of a pattern is the P1 residue (cleaved after) and anything following
# is context the match consumes -- which is why the caller iterates with
# `overlap = true`. Enzymes that cut N-terminal to a residue are expressed with
# a lookahead instead, as `digest_fasta` already does for Asp-N.
#
# Sequences below are built so every rule has a motif to act on, and the
# expected peptides are written out rather than derived, so a pattern that
# changes meaning has to disagree with a human to pass.
#
# Run standalone: julia --project=. test/UnitTests/test_cleavage_presets.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
end

"""Every fully cleaved peptide: no length filtering, no missed cleavages."""
digest_all(sequence, pattern; missed = 0) =
    first(Pioneer.digest_sequence(sequence, Regex(pattern), 100, 1, missed))

@testset "cleavage presets" begin

# ═══════════════════════════════════════════════════════════════════════════════
# Trypsin family, Lys-N, Asp-N, Glu-C
#
#   A  K  G  R  S  T  K  D  E  K  P  A  R  Q
#   1  2  3  4  5  6  7  8  9 10 11 12 13 14
#
# K2 and K7 are followed by ordinary residues; K10 is followed by proline; R4
# and R13 are not. The final Q means no rule has to act at the C-terminus,
# which is handled separately from the regex.
# ═══════════════════════════════════════════════════════════════════════════════

SEQ = "AKGRSTKDEKPARQ"

@testset "Trypsin — after K/R, not before P" begin
    # K10 is blocked by P11, so DEKPAR stays whole.
    @test digest_all(SEQ, "[KR][^P|\$]") == ["AK", "GR", "STK", "DEKPAR", "Q"]
end

@testset "Trypsin/P — after K/R, always" begin
    # The one difference from Trypsin: K10 now cuts, splitting DEKPAR.
    @test digest_all(SEQ, "[KR][^|\$]") == ["AK", "GR", "STK", "DEK", "PAR", "Q"]
end

@testset "Lys-C — after K, not before P" begin
    # R4 and R13 are ignored; K10 is still blocked.
    @test digest_all(SEQ, "K[^P|\$]") == ["AK", "GRSTK", "DEKPARQ"]
end

@testset "Lys-C/P — after K, always" begin
    @test digest_all(SEQ, "K[^|\$]") == ["AK", "GRSTK", "DEK", "PARQ"]
end

@testset "Arg-C — after R, not before P" begin
    @test digest_all(SEQ, "R[^P|\$]") == ["AKGR", "STKDEKPAR", "Q"]
end

@testset "Lys-N — before K" begin
    # Cuts after A1, T6 and E9, the residues preceding each K, so every peptide
    # after the first begins with K.
    @test digest_all(SEQ, "[^K](?=K)") == ["A", "KGRST", "KDE", "KPARQ"]
end

@testset "Asp-N — before D" begin
    # One D, at 8: the cut lands after K7.
    @test digest_all(SEQ, "[^D](?=[D])") == ["AKGRSTK", "DEKPARQ"]
end

@testset "Glu-C, ammonium bicarbonate — after E, not before P" begin
    @test digest_all(SEQ, "E[^P|\$]") == ["AKGRSTKDE", "KPARQ"]
end

@testset "Glu-C, phosphate — after D/E, not before P" begin
    # D8 and E9 are adjacent, so E comes out as a single residue.
    @test digest_all(SEQ, "[DE][^P|\$]") == ["AKGRSTKD", "E", "KPARQ"]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Asp-N with and without glutamate, where the two differ
#
#   G  G  S  E  A  A  D  Q  Q
#   1  2  3  4  5  6  7  8  9
# ═══════════════════════════════════════════════════════════════════════════════

SEQ_DE = "GGSEAADQQ"

@testset "Asp-N ignores glutamate" begin
    @test digest_all(SEQ_DE, "[^D](?=[D])") == ["GGSEAA", "DQQ"]
end

@testset "Asp-N + Glu also cuts before E" begin
    # The extra cut is after S3, giving the E4 peptide its own start.
    @test digest_all(SEQ_DE, "[^DE](?=[DE])") == ["GGS", "EAA", "DQQ"]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Chymotrypsin and CNBr
#
#   A  F  P  G  W  S  Y  L  M  T  F  D
#   1  2  3  4  5  6  7  8  9 10 11 12
#
# F2 sits before a proline, so a rule that respects proline must leave it alone.
# ═══════════════════════════════════════════════════════════════════════════════

SEQ_AROM = "AFPGWSYLMTFD"

@testset "Chymotrypsin — after F/W/Y, not before P" begin
    # F2 blocked by P3; W5, Y7 and F11 all cut.
    @test digest_all(SEQ_AROM, "[FWY][^P|\$]") == ["AFPGW", "SY", "LMTF", "D"]
end

@testset "Chymotrypsin, broad — L and M as well" begin
    @test digest_all(SEQ_AROM, "[FWYLM][^P|\$]") == ["AFPGW", "SY", "L", "M", "TF", "D"]
end

@testset "CNBr — after M" begin
    @test digest_all(SEQ_AROM, "M[^|\$]") == ["AFPGWSYLM", "TFD"]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Properties that hold across the set
# ═══════════════════════════════════════════════════════════════════════════════

@testset "a sequence ending in the P1 residue is not double-counted" begin
    # The C-terminal peptide is emitted by an explicit call after the regex
    # loop, so a pattern able to match the final residue could plausibly add it
    # twice. It does not: the second attempt spans zero residues and is
    # rejected by the length test. Both spellings are therefore safe.
    for pattern in ("[KR][^P|\$]", "[KR](?!P)")
        peptides = digest_all("AAAKBBBK", pattern)
        @test peptides == ["AAAK", "BBBK"]
        @test length(peptides) == length(unique(peptides))
    end
end

@testset "consumed context does not hide the next site" begin
    # Adjacent cleavage sites: the trailing [^P|$] consumes R4 while matching
    # K3, and only `overlap = true` lets R4 match in its own right.
    @test digest_all("AAKRBB", "[KR][^P|\$]") == ["AAK", "R", "BB"]
end

@testset "missed cleavages extend rather than replace" begin
    # One missed cleavage keeps every fully cleaved peptide and adds each
    # adjacent pair.
    peptides = digest_all(SEQ, "[KR][^P|\$]"; missed = 1)
    for fully in ["AK", "GR", "STK", "DEKPAR", "Q"]
        @test fully in peptides
    end
    for joined in ["AKGR", "GRSTK", "STKDEKPAR", "DEKPARQ"]
        @test joined in peptides
    end
end

@testset "Trypsin and Trypsin/P differ only at proline" begin
    # No proline anywhere: the two rules must agree exactly.
    without_proline = "AAKBBRCCKDDR"
    @test digest_all(without_proline, "[KR][^P|\$]") ==
          digest_all(without_proline, "[KR][^|\$]")
end

end # top-level testset
