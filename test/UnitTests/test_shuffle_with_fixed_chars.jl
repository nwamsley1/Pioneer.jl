using Test
using Random

# For testing just the function without dependencies
function shuffle_fast_with_positions_and_fixed_chars!(
    s::String, 
    positions::Vector{UInt8}, 
    fixed_chars::Set{Char},
    fixed_positions::Vector{Int},
    movable_positions::Vector{Int},
    temp_positions::Vector{UInt8}
)
    ss = sizeof(s)
    l = length(s)
    
    # Count fixed and movable positions
    n_fixed = 0
    n_movable = 0
    
    # Create indices vector for byte positions
    v = Vector{Int}(undef, l)
    i = 1
    for j in 1:l
        v[j] = i
        i = nextind(s, i)
    end
    
    # Identify fixed and movable positions
    p = pointer(s)
    for j in 1:l
        c = Char(unsafe_load(p, v[j]))
        if j == l || c in fixed_chars  # Last position or fixed character
            n_fixed += 1
            fixed_positions[n_fixed] = j
        else
            n_movable += 1
            movable_positions[n_movable] = j
        end
    end
    
    # Generate permutation only for movable positions
    if n_movable > 0
        perm = randperm(n_movable)
        
        # Copy positions for movable characters to temp storage
        for idx in 1:n_movable
            temp_positions[idx] = positions[movable_positions[idx]]
        end
        
        # Apply permutation to positions
        for (new_idx, old_idx) in enumerate(perm)
            positions[movable_positions[new_idx]] = temp_positions[old_idx]
        end
        
        # Build the output string
        u = Vector{UInt8}(undef, ss)
        
        # Fill in the shuffled string
        for j in 1:l
            # Check if j is in fixed_positions (up to n_fixed)
            is_fixed = false
            for k in 1:n_fixed
                if fixed_positions[k] == j
                    is_fixed = true
                    break
                end
            end
            
            if is_fixed
                # Keep fixed characters in place
                u[v[j]] = unsafe_load(p, v[j])
            else
                # Find position in movable_positions
                idx = 0
                for k in 1:n_movable
                    if movable_positions[k] == j
                        idx = k
                        break
                    end
                end
                source_pos = movable_positions[perm[idx]]
                u[v[j]] = unsafe_load(p, v[source_pos])
            end
        end
        
        return String(u)
    else
        # All positions are fixed, return original string
        return s
    end
end

@testset "shuffle_fast_with_positions_and_fixed_chars!" begin
    
    @testset "Basic functionality - multiple runs" begin
        # Set seed for reproducibility
        Random.seed!(12345)
        
        sequence = "ABCDEFGHIJK"
        fixed_chars = Set(['C', 'F', 'I'])
        
        # Pre-allocate vectors
        fixed_positions = Vector{Int}(undef, 255)
        movable_positions = Vector{Int}(undef, 255)
        temp_positions = Vector{UInt8}(undef, 255)
        
        # Run multiple times to ensure consistency
        for run in 1:100
            positions = collect(UInt8, 1:length(sequence))
            
            shuffled = shuffle_fast_with_positions_and_fixed_chars!(
                sequence, positions, fixed_chars,
                fixed_positions, movable_positions, temp_positions
            )
            
            # Check that fixed characters are in their original positions
            @test shuffled[3] == 'C'  # Position 3
            @test shuffled[6] == 'F'  # Position 6
            @test shuffled[9] == 'I'  # Position 9
            @test shuffled[end] == 'K'  # Last position always fixed
            
            # Check that all original characters are present
            @test sort(collect(shuffled)) == sort(collect(sequence))
        end
    end
    
    @testset "Position tracking - multiple runs" begin
        Random.seed!(54321)
        
        sequence = "PEPTIDEK"
        fixed_chars = Set(['P', 'K'])
        
        fixed_positions = Vector{Int}(undef, 255)
        movable_positions = Vector{Int}(undef, 255)
        temp_positions = Vector{UInt8}(undef, 255)
        
        for run in 1:100
            positions = collect(UInt8, 1:length(sequence))
            
            shuffled = shuffle_fast_with_positions_and_fixed_chars!(
                sequence, positions, fixed_chars,
                fixed_positions, movable_positions, temp_positions
            )
            
            # Verify position tracking
            for (idx, char) in enumerate(shuffled)
                original_idx = positions[idx]
                @test sequence[original_idx] == char
            end
            
            # Fixed positions should map to themselves
            @test positions[1] == 1  # First P
            @test positions[3] == 3  # Second P  
            @test positions[8] == 8  # K (last position)
        end
    end
    
    @testset "Edge cases - multiple runs" begin
        Random.seed!(11111)
        
        fixed_positions = Vector{Int}(undef, 255)
        movable_positions = Vector{Int}(undef, 255)
        temp_positions = Vector{UInt8}(undef, 255)
        
        # Test with empty sequence
        for run in 1:10
            sequence = ""
            fixed_chars = Set(['A'])
            positions = UInt8[]
            
            shuffled = shuffle_fast_with_positions_and_fixed_chars!(
                sequence, positions, fixed_chars,
                fixed_positions, movable_positions, temp_positions
            )
            @test shuffled == ""
        end
        
        # Test with single character
        for run in 1:10
            sequence = "A"
            positions = UInt8[1]
            fixed_chars = Set(['A'])
            shuffled = shuffle_fast_with_positions_and_fixed_chars!(
                sequence, positions, fixed_chars,
                fixed_positions, movable_positions, temp_positions
            )
            @test shuffled == "A"
        end
        
        # Test with all characters fixed
        for run in 1:50
            sequence = "AAAA"
            fixed_chars = Set(['A'])
            positions = collect(UInt8, 1:4)
            shuffled = shuffle_fast_with_positions_and_fixed_chars!(
                sequence, positions, fixed_chars,
                fixed_positions, movable_positions, temp_positions
            )
            @test shuffled == "AAAA"
        end
        
        # Test with no fixed characters (except last)
        for run in 1:100
            sequence = "ABCDE"
            fixed_chars = Set{Char}()
            positions = collect(UInt8, 1:5)
            shuffled = shuffle_fast_with_positions_and_fixed_chars!(
                sequence, positions, fixed_chars,
                fixed_positions, movable_positions, temp_positions
            )
            @test shuffled[end] == 'E'  # Last always fixed
            @test sort(collect(shuffled)) == sort(collect(sequence))
        end
    end
    
    @testset "Multiple runs produce different results" begin
        Random.seed!(99999)
        
        sequence = "ABCDEFGHIJK"
        fixed_chars = Set(['C'])
        
        fixed_positions = Vector{Int}(undef, 255)
        movable_positions = Vector{Int}(undef, 255)
        temp_positions = Vector{UInt8}(undef, 255)
        
        results = Set{String}()
        for _ in 1:50  # More runs to ensure we see variation
            positions = collect(UInt8, 1:length(sequence))
            
            shuffled = shuffle_fast_with_positions_and_fixed_chars!(
                sequence, positions, fixed_chars,
                fixed_positions, movable_positions, temp_positions
            )
            push!(results, shuffled)
        end
        
        # Should get multiple different shuffles (with 9 movable positions, we expect many variations)
        @test length(results) > 10
        
        # Verify all results are valid
        for result in results
            @test result[3] == 'C'  # Fixed char
            @test result[end] == 'K'  # Last char
            @test sort(collect(result)) == sort(collect(sequence))
        end
    end
    
    @testset "Complex peptide sequences - multiple runs" begin
        Random.seed!(7777)
        
        sequence = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
        fixed_chars = Set(['R', 'K', 'H', 'P'])
        
        fixed_positions = Vector{Int}(undef, 255)
        movable_positions = Vector{Int}(undef, 255)
        temp_positions = Vector{UInt8}(undef, 255)
        
        for run in 1:100
            positions = collect(UInt8, 1:length(sequence))
            
            shuffled = shuffle_fast_with_positions_and_fixed_chars!(
                sequence, positions, fixed_chars,
                fixed_positions, movable_positions, temp_positions
            )
            
            # Check all R, K, H, P remain in place
            for (idx, char) in enumerate(sequence)
                if char in fixed_chars
                    @test shuffled[idx] == char
                end
            end
            
            # Last character always fixed
            @test shuffled[end] == sequence[end]
            
            # All characters preserved
            @test sort(collect(shuffled)) == sort(collect(sequence))
        end
    end
    
    @testset "Pre-allocated vector reuse - multiple runs" begin
        Random.seed!(33333)
        
        fixed_positions = Vector{Int}(undef, 255)
        movable_positions = Vector{Int}(undef, 255)
        temp_positions = Vector{UInt8}(undef, 255)
        fixed_chars = Set(['X'])
        
        # Run multiple times with different sequences
        sequences = ["ABCXDEF", "GHIXJKL", "MNOXPQR"]
        
        # Run each sequence multiple times
        for seq in sequences
            for run in 1:50
                positions = collect(UInt8, 1:length(seq))
                shuffled = shuffle_fast_with_positions_and_fixed_chars!(
                    seq, positions, fixed_chars,
                    fixed_positions, movable_positions, temp_positions
                )
                
                # Verify X stays in place
                for (idx, char) in enumerate(seq)
                    if char == 'X'
                        @test shuffled[idx] == 'X'
                    end
                end
                
                # Last char fixed
                @test shuffled[end] == seq[end]
                @test sort(collect(shuffled)) == sort(collect(seq))
            end
        end
    end
    
    @testset "Peptide-specific sequences - multiple runs" begin
        Random.seed!(44444)
        
        # Test with peptide sequences containing fixed amino acids
        sequence = "MKTVRQERLKSIVR"
        fixed_chars = Set(['R', 'K'])  # Common enzymatic cleavage sites
        
        fixed_positions = Vector{Int}(undef, 255)
        movable_positions = Vector{Int}(undef, 255)
        temp_positions = Vector{UInt8}(undef, 255)
        
        for run in 1:50
            positions = collect(UInt8, 1:length(sequence))
            
            shuffled = shuffle_fast_with_positions_and_fixed_chars!(
                sequence, positions, fixed_chars,
                fixed_positions, movable_positions, temp_positions
            )
            
            # Check that K and R stay in their original positions
            # Original sequence: MKTVRQERLKSIVR
            # K at position 2, R at positions 5, 8, 14
            # K at position 10
            @test shuffled[2] == 'K'   # K at position 2
            @test shuffled[5] == 'R'   # R at position 5
            @test shuffled[8] == 'R'   # R at position 8
            @test shuffled[10] == 'K'  # K at position 10
            @test shuffled[14] == 'R'  # R at position 14 (last)
            @test sort(collect(shuffled)) == sort(collect(sequence))
        end
    end
    
    @testset "Stress test - many fixed characters" begin
        Random.seed!(55555)
        
        # Test with many fixed characters
        sequence = "RKPHDERKPHDERKPHDERKPHDE"
        fixed_chars = Set(['R', 'K', 'P', 'H', 'D'])  # Most chars are fixed
        
        fixed_positions = Vector{Int}(undef, 255)
        movable_positions = Vector{Int}(undef, 255)
        temp_positions = Vector{UInt8}(undef, 255)
        
        for run in 1:100
            positions = collect(UInt8, 1:length(sequence))
            
            shuffled = shuffle_fast_with_positions_and_fixed_chars!(
                sequence, positions, fixed_chars,
                fixed_positions, movable_positions, temp_positions
            )
            
            # Verify all fixed characters remain in place
            for (idx, char) in enumerate(sequence)
                if char in fixed_chars
                    @test shuffled[idx] == char
                end
            end
            
            @test shuffled[end] == 'E'  # Last char
            @test sort(collect(shuffled)) == sort(collect(sequence))
        end
    end
end

# Run tests
@testset "Performance characteristics" begin
    Random.seed!(66666)
    
    sequence = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    fixed_chars = Set(['A', 'M', 'Z'])
    
    # Pre-allocate all vectors
    positions = collect(UInt8, 1:length(sequence))
    fixed_positions = Vector{Int}(undef, 255)
    movable_positions = Vector{Int}(undef, 255)
    temp_positions = Vector{UInt8}(undef, 255)
    
    # Warm up
    for _ in 1:10
        shuffle_fast_with_positions_and_fixed_chars!(
            sequence, positions, fixed_chars,
            fixed_positions, movable_positions, temp_positions
        )
    end
    
    # Measure allocations (should be minimal - just for output string and randperm)
    allocs = @allocated shuffle_fast_with_positions_and_fixed_chars!(
        sequence, positions, fixed_chars,
        fixed_positions, movable_positions, temp_positions
    )
    
    @test allocs < 1000  # Should be much less than this
end