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
function selectTransitionsPRM(window_center::Float32, 
                                precursorList::Vector{Precursor}, 
                                prec_id_to_transitions::Dictionary{UInt32, Vector{Transition}},
                                right_precursor_tolerance::Float32,
                                left_precursor_tolerance::Float32,
                                transition_charges::Vector{UInt8},
                                transition_isotopes::Vector{UInt8},
                                b_start::Int64,
                                y_start::Int64,
                                fragment_match_ppm::Float32)

    function getPrecursors(window_center::Float32, precursorList::Vector{Precursor}, left_precursor_tolerance::Float32, right_precursor_tolerance::Float32)
        l_bnd, u_bnd = window_center - left_precursor_tolerance, window_center + right_precursor_tolerance
        start, stop = searchsortedfirst(precursorList, l_bnd,lt=(t,x)->getMZ(t)<x), searchsortedlast(precursorList, u_bnd,lt=(x,t)->getMZ(t)>x)
        return @view(precursorList[start:stop])
    end

    transitions = Vector{Transition}();

    for precursor in getPrecursors(window_center, precursorList, left_precursor_tolerance, right_precursor_tolerance)
        for charge in transition_charges, isotope in transition_isotopes
            if !isassigned(prec_id_to_transitions, getPrecID(precursor))
                insert!(prec_id_to_transitions, 
                        getPrecID(precursor),
                        getTransitions(precursor, 
                                                charge = charge, 
                                                isotope = isotope, 
                                                b_start = b_start, 
                                                y_start = y_start,
                                                ppm = fragment_match_ppm)
                        ) 
            end
        end
        append!(transitions, prec_id_to_transitions[getPrecID(precursor)])
    end
    sort!(transitions, by=x->getMZ(x))

    transitions
end