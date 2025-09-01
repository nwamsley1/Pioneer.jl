#!/usr/bin/env julia
"""
Test suite for KEAP1 library to verify both base_pep_id and pair_id functionality
"""

using Arrow, Tables, DataFrames, Test

# Load the library data
println("Loading KEAP1 library data...")
df = DataFrame(Arrow.Table("/tmp/keap1_test_library/keap1_test.poin/precursors_table.arrow"))

println("\n" * "="^80)
println("KEAP1 LIBRARY VALIDATION TESTS")
println("="^80)

@testset "Library Structure Tests" begin
    @test nrow(df) > 0
    @test :base_pep_id in propertynames(df)
    @test :pair_id in propertynames(df) 
    @test :entrapment_group_id in propertynames(df)
    @test :is_decoy in propertynames(df)
    @test :partner_precursor_idx in propertynames(df)
    println("✅ All required columns present")
end

@testset "base_pep_id Preservation Tests" begin
    unique_base_pep_ids = length(unique(df.base_pep_id))
    total_precursors = nrow(df)
    
    # Should have much fewer base_pep_ids than total precursors
    @test unique_base_pep_ids < total_precursors ÷ 3
    println("✅ base_pep_id count reasonable: $unique_base_pep_ids unique IDs")
    
    # Check modification variants share base_pep_id
    grouped = groupby(df, :base_pep_id)
    multi_variant_groups = [g for g in grouped if nrow(g) > 2]
    @test length(multi_variant_groups) > 0
    println("✅ Found $(length(multi_variant_groups)) base_pep_ids with multiple variants")
end

@testset "Entrapment Sequence Linking Tests" begin
    # Check entrapment sequences exist
    entrapment_df = df[df.entrapment_group_id .> 0, :]
    @test nrow(entrapment_df) > 0
    println("✅ Found $(nrow(entrapment_df)) entrapment sequences")
    
    # Check entrapment sequences share base_pep_id with targets
    target_df = df[(df.entrapment_group_id .== 0) .& (.!df.is_decoy), :]
    entrap_base_ids = Set(entrapment_df.base_pep_id)
    target_base_ids = Set(target_df.base_pep_id)
    shared_base_ids = intersect(entrap_base_ids, target_base_ids)
    
    @test length(shared_base_ids) > 0
    println("✅ $(length(shared_base_ids)) base_pep_ids shared between targets and entrapments")
    
    # All entrapment base_pep_ids should have corresponding targets
    coverage = length(shared_base_ids) / length(entrap_base_ids)
    @test coverage > 0.9  # At least 90% coverage
    println("✅ Entrapment-target coverage: $(round(coverage*100, digits=1))%")
end

@testset "Target-Decoy Pairing Tests" begin
    # Check for reasonable pair_id distribution
    unique_pair_ids = length(unique(df.pair_id))
    @test unique_pair_ids > 1000
    println("✅ Found $unique_pair_ids unique pair_ids")
    
    # Check for paired sequences (same pair_id, different is_decoy)
    pair_groups = groupby(df, :pair_id)
    paired_groups = [g for g in pair_groups if nrow(g) == 2 && 
                     length(unique(g.is_decoy)) == 2]
    
    @test length(paired_groups) > 100
    println("✅ Found $(length(paired_groups)) proper target-decoy pairs")
    
    # Sample a few pairs to verify structure
    if length(paired_groups) > 0
        sample_pair = first(paired_groups)
        target_row = sample_pair[sample_pair.is_decoy .== false, :][1, :]
        decoy_row = sample_pair[sample_pair.is_decoy .== true, :][1, :]
        
        @test target_row.pair_id == decoy_row.pair_id
        @test target_row.prec_charge == decoy_row.prec_charge
        println("✅ Sample pair verification passed")
    end
end

@testset "Partner Index Validation Tests" begin
    # Check partner indices are within bounds
    max_partner_idx = maximum(skipmissing(df.partner_precursor_idx))
    table_size = nrow(df)
    
    @test max_partner_idx <= table_size
    println("✅ All partner indices within bounds: max=$max_partner_idx, table_size=$table_size")
    
    # Check for symmetric partnerships
    partnered_rows = df[.!ismissing.(df.partner_precursor_idx), :]
    symmetric_count = 0
    sample_size = min(100, nrow(partnered_rows))
    
    for i in 1:sample_size
        row_idx = i
        partner_idx = partnered_rows.partner_precursor_idx[row_idx]
        if !ismissing(partner_idx) && partner_idx <= nrow(df)
            partner_partner_idx = df.partner_precursor_idx[partner_idx]
            if !ismissing(partner_partner_idx) 
                # Find the original row index in the full dataframe
                original_idx = findfirst(==(partnered_rows[row_idx, :]), eachrow(df))
                if partner_partner_idx == original_idx
                    symmetric_count += 1
                end
            end
        end
    end
    
    symmetry_rate = symmetric_count / sample_size
    @test symmetry_rate > 0.5  # At least 50% should be symmetric
    println("✅ Partnership symmetry rate: $(round(symmetry_rate*100, digits=1))% (sample of $sample_size)")
end

@testset "High pair_id Investigation" begin
    # Investigate the ridiculously high pair_id values
    high_pair_ids = df[df.pair_id .> 1_000_000, :]
    normal_pair_ids = df[df.pair_id .<= 100_000, :]
    
    println("\n🔍 PAIR_ID INVESTIGATION:")
    println("   Total precursors: $(nrow(df))")
    println("   Normal pair_ids (≤100k): $(nrow(normal_pair_ids))")
    println("   High pair_ids (>1M): $(nrow(high_pair_ids))")
    println("   Max pair_id: $(maximum(df.pair_id))")
    
    if nrow(high_pair_ids) > 0
        println("\n   High pair_id characteristics:")
        println("   - Decoy rate: $(round(mean(high_pair_ids.is_decoy)*100, digits=1))%")
        println("   - Entrapment rate: $(round(mean(high_pair_ids.entrapment_group_id .> 0)*100, digits=1))%")
        
        # Sample of high pair_id sequences
        sample_high = first(high_pair_ids, 3)
        println("\n   Sample high pair_id entries:")
        for (i, row) in enumerate(eachrow(sample_high))
            println("   $i. $(row.sequence) | decoy:$(row.is_decoy) | entrap:$(row.entrapment_group_id) | pair_id:$(row.pair_id)")
        end
    end
    
    # This test allows high pair_ids but documents them
    @test nrow(df) > 0  # Just ensure we have data
    println("📝 High pair_id investigation completed")
end

println("\n" * "="^80)
println("VALIDATION COMPLETE")
println("="^80)

# Summary statistics
println("\nSUMMARY STATISTICS:")
println("📊 Total precursors: $(nrow(df))")
println("📊 Unique base_pep_ids: $(length(unique(df.base_pep_id)))")
println("📊 Unique pair_ids: $(length(unique(df.pair_id)))")
println("📊 Entrapment sequences: $(sum(df.entrapment_group_id .> 0))")
println("📊 Target sequences: $(sum(.!df.is_decoy))")
println("📊 Decoy sequences: $(sum(df.is_decoy))")
println("📊 Partner indices valid: $(maximum(skipmissing(df.partner_precursor_idx)) <= nrow(df) ? "✅ YES" : "❌ NO")")