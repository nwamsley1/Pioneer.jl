struct PeptideMod
    position::UInt8
    aa::Char
    mod_name::String
end

getPosition(mod::PeptideMod) = mod.position
getAA(mod::PeptideMod) = mod.aa
getModName(mod::PeptideMod) = mod.mod_name
getModString(mod::PeptideMod) = join(['(', string(mod.position), ',', mod.aa, ',', mod.mod_name, ')'])
getModString(mods::Vector{PeptideMod}) = join([getModString(mod) for mod in mods])
getModString(mods::Missing) = missing
"""
Defines a custom sort order for PeptideMod objects.
Sorts by position, then by mod_name, then by aa.
"""
function Base.isless(a::PeptideMod, b::PeptideMod)
    # Compare positions first
    if a.position != b.position
        return a.position < b.position
    end
    
    # If positions are equal, compare mod_name
    if a.mod_name != b.mod_name
        return a.mod_name < b.mod_name
    end
    
    # If positions and mod_name are equal, compare aa
    return a.aa < b.aa
end