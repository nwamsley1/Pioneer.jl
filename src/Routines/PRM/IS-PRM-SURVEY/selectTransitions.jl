"""
    selectTransitionsPRM(window_center::Float32, precursorList::Vector{Precursor}, params)

Given a SORTED vector of `Precursor` in order os ascending MZ, gets the subset of `Precursor` with MZ within the tolerance.
The tolerance is specified based on the `window_center` and `params[:lower_tol]` and `params[:upper_tol]` . Every routine implementing a method `SearchRaw` should
implement a `selectTransitions` method. 

### Input

- `window_center::Float32` -- The isolation window center for an MS2 scan. 
- `precursorList::Vector{Precursor}` -- List of possible `Precursor` 
- `MS_TABLE::Arrow.Table` -- Search parameters. A named tuple that must have fields [:upper_tol] and [:lower_tol]

### Output
- Returns Vector{Precursor} which is a subset of `precursorList` that satisfies the constraints. This output is 
also sorted by MZ just like `precursorLit`

### Notes

### Examples 

"""
function selectTransitions(window_center::T, 
                                ptable::PrecursorTable, 
                                right_precursor_tolerance::U,
                                left_precursor_tolerance::U,
                                #transition_charges::Vector{UInt8},
                                #transition_isotopes::Vector{UInt8},
                                #b_start::Int64,
                                #y_start::Int64,
                                #fragment_match_ppm::Float32
                                ) where {T,U<:Real}

    transitions = Vector{Transition}();

    for prec_id in precursorRangeQuery(ptable, window_center, left_precursor_tolerance, right_precursor_tolerance)
        if isassigned(getPrecIDToTransitions(ptable), prec_id)
            append!(transitions, getTransitions(ptable, prec_id))
        end
    end
    sort!(transitions, by=x->getMZ(x))

    transitions
end