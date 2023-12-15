isoxml = readxml("./data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml")
isosplines = root(isoxml)

for node in eachelement(isosplines)
    if hasnode(node)
        if node.name == "model"
        for subnode in eachelement(node)
            println(node.name)
            println(node.content)
        end
        break
        end
    end
end
collect(eachelement(isosplines))[100]

xdoc = parse_file("./data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml")

i = 0
for c in root(xdoc)["model"]
    println(attributes_dict(c))
    println(attributes_dict(c["knots"][1]))
    println(content(c["knots"][1]))

    println(attributes_dict(c["coefficients"][1]))
    println(content(c["coefficients"][1]))

end