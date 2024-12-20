const MASSES = Dict{String, Float64}(
    "H" => 1.007825035,
    "C" => 12.0,
    "N" => 14.003074,
    "O" => 15.994915,
    "P" => 30.973762,
    "S" => 31.972071,
    "H2O" => 18.010565,
    "NH3" => 17.026549,
    "CO" => 27.994915,
    "CO2" => 43.989830,
    "H3PO4" => 97.976896
)

"""
Parse a molecular formula string into a mass difference.
"""
function parse_modification_mass(formula::String)::Float64
    total_mass = 0.0
    current_number = ""
    current_element = ""
    sign = 1.0
    
    if startswith(formula, '-')
        sign = -1.0
        formula = formula[2:end]
    elseif startswith(formula, '+')
        formula = formula[2:end]
    end
    
    for char in formula
        if isdigit(char)
            current_number *= char
        elseif isuppercase(char)
            if !isempty(current_element)
                count = isempty(current_number) ? 1 : parse(Int, current_number)
                total_mass += MASSES[current_element] * count
                current_number = ""
            end
            current_element = string(char)
        elseif islowercase(char)
            current_element *= char
        end
    end
    
    # Handle last element
    if !isempty(current_element)
        count = isempty(current_number) ? 1 : parse(Int, current_number)
        total_mass += MASSES[current_element] * count
    end
    
    return sign * total_mass
end

"""
Count sulfur atoms in modification formula.
"""
function count_sulfur_atoms(formula::String)::Int8
    matches = eachmatch(r"S([0-9]*)", formula)
    count = 0
    for m in matches
        count += isempty(m.captures[1]) ? 1 : parse(Int8, m.captures[1])
    end
    return sign(formula) * count
end

"""
Get modification indices from sequence.
"""
function get_modification_indices(sequence::String, aa_pattern::Regex)::Vector{Int}
    return [m.offset for m in eachmatch(aa_pattern, sequence)]
end

"""
Apply modifications to sequence.
"""
function apply_modifications(sequence::String, mods::String)::Tuple{String,Float64}
    if ismissing(mods)
        return (sequence, 0.0)
    end
    
    total_mass_shift = 0.0
    modified_sequence = sequence
    
    for mod in eachmatch(r"\((\d+),(\w+),([^)]+)\)", mods)
        index = parse(Int, mod.captures[1])
        aa = mod.captures[2]
        mod_name = mod.captures[3]
        
        # Insert modification marker and calculate mass shift
        modified_sequence = modified_sequence[1:index-1] * 
                          aa * "[" * mod_name * "]" * 
                          modified_sequence[index+1:end]
                          
        total_mass_shift += parse_modification_mass(mod_name)
    end
    
    return (modified_sequence, total_mass_shift)
end
