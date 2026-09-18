@testset "Attach Pass-1 predictions while merging folds" begin
    function write_scored_fold_fixture(dir, name, ids; infold=:scores)
        path = joinpath(dir, name * ".arrow")
        n = length(ids)
        main = DataFrame(
            precursor_idx=UInt32.(ids),
            scan_idx=UInt32.(100 .+ ids),
            target=iseven.(ids),
            cv_fold=fill(UInt8(endswith(name, "fold1")), n),
            trace_prob=fill(-1f0, n),
            mbr_recovered=trues(n),
        )
        predictions = DataFrame(
            precursor_idx=main.precursor_idx,
            scan_idx=main.scan_idx,
            trace_prob_prepass=Float64.(ids) ./ 10,
        )
        if infold !== :absent
            predictions[!, :trace_prob_infold] = infold === :nan ?
                fill(NaN32, n) : Float64.(ids) ./ 20
        end
        auxiliary = Int16.(10 .+ ids)
        auxiliary_path = path * ".features.sidecar.arrow"
        pass1_path = path * Pioneer.PASS1_SIDECAR_SUFFIX
        Arrow.write(path, main)
        Arrow.write(auxiliary_path, (; auxiliary, trace_prob_prepass=fill(-99f0, n)))
        Arrow.write(pass1_path, predictions)

        expected = copy(main)
        expected[!, :decoy] = .!main.target
        expected[!, :trace_prob_prepass] = Float32.(predictions.trace_prob_prepass)
        if infold !== :absent
            expected[!, :trace_prob_infold] = Float32.(predictions.trace_prob_infold)
        end
        expected[!, :trace_prob] = copy(expected.trace_prob_prepass)
        expected[!, :mbr_recovered] = falses(n)
        expected[!, :auxiliary] = auxiliary
        return (; path, pass1_path, expected, inputs=[path, auxiliary_path, pass1_path])
    end

    function check_merged_table(path, expected)
        actual = DataFrame(Arrow.Table(path))
        @test names(actual) == names(expected)
        @test all(name -> isequal(actual[!, name], expected[!, name]), names(expected))
        @test all(name -> eltype(actual[!, name]) === eltype(expected[!, name]), names(expected))
        @test Pioneer.PSMFileReference(path).row_count == nrow(expected)
    end

    for infold in (:scores, :nan, :absent)
        @testset "In-fold scores: $infold" begin
            mktempdir() do dir
                fold0 = write_scored_fold_fixture(dir, "run_fold0", [4, 1]; infold)
                fold1 = write_scored_fold_fixture(dir, "run_fold1", [6, 2]; infold)
                inputs = vcat(fold0.inputs, fold1.inputs)
                original_bytes = Dict(path => read(path) for path in inputs)
                merged_path = joinpath(dir, "run.arrow")

                result = Pioneer._merge_scored_folds!([fold0.path, fold1.path], merged_path)

                @test result.rows == 4
                @test Set(result.cleanup_paths) == Set(inputs)
                @test all(path -> read(path) == original_bytes[path], inputs)
                check_merged_table(merged_path, vcat(fold0.expected, fold1.expected))
            end
        end
    end

    @testset "Missing and empty folds" begin
        mktempdir() do dir
            absent_path = joinpath(dir, "absent_fold.arrow")
            merged_path = joinpath(dir, "run.arrow")
            @test isnothing(Pioneer._merge_scored_folds!([absent_path], merged_path))
            @test !isfile(merged_path)

            fold1 = write_scored_fold_fixture(dir, "run_fold1", [3, 1])
            result = Pioneer._merge_scored_folds!([absent_path, fold1.path], merged_path)
            @test result.rows == 2
            @test Set(result.cleanup_paths) == Set(fold1.inputs)
            check_merged_table(merged_path, fold1.expected)

            empty_fold = write_scored_fold_fixture(dir, "run_fold0", Int[])
            result = Pioneer._merge_scored_folds!([empty_fold.path, fold1.path], merged_path)
            @test result.rows == 2
            @test Set(result.cleanup_paths) == Set(vcat(empty_fold.inputs, fold1.inputs))
            check_merged_table(merged_path, fold1.expected)

            result = Pioneer._merge_scored_folds!([empty_fold.path], merged_path)
            @test result.rows == 0
            @test Set(result.cleanup_paths) == Set(empty_fold.inputs)
            check_merged_table(merged_path, empty_fold.expected)
        end
    end

    for mismatch in (:row_count, :precursor_idx, :scan_idx, :missing_predictions)
        @testset "Reject $mismatch without altering inputs or output" begin
            mktempdir() do dir
                valid_fold = write_scored_fold_fixture(dir, "run_fold0", [6, 3])
                fold = write_scored_fold_fixture(dir, "run_fold1", [4, 1])
                predictions = DataFrame(Arrow.Table(fold.pass1_path); copycols=true)
                if mismatch === :missing_predictions
                    rm(fold.pass1_path)
                else
                    if mismatch === :row_count
                        predictions = predictions[1:1, :]
                    else
                        reverse!(predictions[!, mismatch])
                    end
                    Pioneer.writeArrow(fold.pass1_path, predictions)
                end
                inputs = filter(isfile, vcat(valid_fold.inputs, fold.inputs))
                original_bytes = Dict(path => read(path) for path in inputs)
                merged_path = joinpath(dir, "run.arrow")
                Arrow.write(merged_path, (; marker=[123]))
                original_output = read(merged_path)

                @test_throws ErrorException Pioneer._merge_scored_folds!(
                    [valid_fold.path, fold.path], merged_path,
                )
                @test read(merged_path) == original_output
                @test all(path -> read(path) == original_bytes[path], inputs)
            end
        end
    end

    @testset "Write failure retains source folds and predictions" begin
        mktempdir() do dir
            fold = write_scored_fold_fixture(dir, "run_fold0", [4, 1])
            original_bytes = Dict(path => read(path) for path in fold.inputs)
            merged_path = joinpath(dir, "missing_directory", "run.arrow")

            @test_throws Exception Pioneer._merge_scored_folds!([fold.path], merged_path)
            @test !isfile(merged_path)
            @test all(path -> read(path) == original_bytes[path], fold.inputs)
        end
    end
end
