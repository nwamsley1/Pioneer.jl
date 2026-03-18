using PProf, CodecZlib

function analyze_profile(filepath)
    PB = PProf.ProtoBuf

    data = transcode(GzipDecompressor, read(filepath))
    decoder = PB.ProtoDecoder(IOBuffer(data))
    prof = PB.decode(decoder, PProf.perftools.profiles.Profile)

    strings = prof.string_table
    loc_map = Dict(l.id => l for l in prof.location)
    funcs = getproperty(prof, Symbol("#function"))
    func_map = Dict(f.id => f for f in funcs)

    self_counts = Dict{String, Int64}()
    inclusive_counts = Dict{String, Int64}()
    total_samples = 0

    for sample in prof.sample
        val = sample.value[1]
        total_samples += val

        if !isempty(sample.location_id)
            loc = loc_map[sample.location_id[1]]
            if !isempty(loc.line)
                fname = strings[func_map[loc.line[1].function_id].name + 1]
                self_counts[fname] = get(self_counts, fname, 0) + val
            end
        end

        seen = Set{String}()
        for loc_id in sample.location_id
            loc = loc_map[loc_id]
            for line in loc.line
                fname = strings[func_map[line.function_id].name + 1]
                if !(fname in seen)
                    push!(seen, fname)
                    inclusive_counts[fname] = get(inclusive_counts, fname, 0) + val
                end
            end
        end
    end

    println("Total samples: $total_samples\n")

    println("TOP 50 BY SELF (exclusive) TIME:")
    println("="^100)
    sorted_self = sort(collect(self_counts), by=x->-x[2])
    for (i, (name, count)) in enumerate(sorted_self[1:min(50, length(sorted_self))])
        pct = round(100.0 * count / total_samples, digits=1)
        name_trunc = length(name) > 65 ? name[1:65] : name
        println("  $(lpad(i, 3)). $(rpad(name_trunc, 65)) $(lpad(count, 8))  $(lpad(pct, 5))%")
    end

    println("\n\nTOP 50 BY INCLUSIVE TIME:")
    println("="^100)
    sorted_inc = sort(collect(inclusive_counts), by=x->-x[2])
    for (i, (name, count)) in enumerate(sorted_inc[1:min(50, length(sorted_inc))])
        pct = round(100.0 * count / total_samples, digits=1)
        name_trunc = length(name) > 65 ? name[1:65] : name
        println("  $(lpad(i, 3)). $(rpad(name_trunc, 65)) $(lpad(count, 8))  $(lpad(pct, 5))%")
    end
end

analyze_profile(ARGS[1])
