include("src/Utils/isotopeSplines.jl")
iso_splines = parseIsoXML("./data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml");


#Use exact isotope composition from expasy.org 
#https://education.expasy.org/student_projects/isotopident/cgi-bin/iso.pl

pep_a = "PEPTIDEPEPTIDEMMM"
pep_mass_a = 1973.831f0
pep_sulfur_count = 4 #1 means zero for for 3 sulfurs must be 4
exact_iso_abundances = Float32[0.953933484,1.000000000,0.707123653,0.373977793,0.162691088,0.060424297]
spline_iso_abundances = [iso_splines.splines[pep_sulfur_count][x](pep_mass_a) for x in range(1, 6)]
spline_iso_abundances = spline_iso_abundances./maximum(spline_iso_abundances)
percent_deviance = 100.0*(exact_iso_abundances .- spline_iso_abundances)./exact_iso_abundances

pep_a = "PEPTIDEPEPTIDE"
pep_mass_a = 1580.709f0
pep_sulfur_count = 1
exact_iso_abundances = Float32[1.000000000,0.840226254,0.406760020,0.144049251,0.041139503,0.009968661]
spline_iso_abundances = [iso_splines.splines[pep_sulfur_count][x](pep_mass_a) for x in range(1, 6)]
spline_iso_abundances = spline_iso_abundances./maximum(spline_iso_abundances)
percent_deviance = 100.0*(exact_iso_abundances .- spline_iso_abundances)./exact_iso_abundances
 
pep_a = "VLSMSNLR"
pep_mass_a = 918.496f0
pep_sulfur_count = 2
exact_iso_abundances = Float32[1.000000000,0.493068635,0.187476610,0.052138016,0.011652472,0.002172227]
spline_iso_abundances = [iso_splines.splines[pep_sulfur_count][x](pep_mass_a) for x in range(1, 6)]
spline_iso_abundances = spline_iso_abundances./maximum(spline_iso_abundances)
percent_deviance = 100.0*(exact_iso_abundances .- spline_iso_abundances)./exact_iso_abundances

pep_a = "NNSPDGNTDSSALDCYNPMTNQWSPCAPMSVPRNR"
pep_mass_a = 3952.662f0
pep_sulfur_count = 5
exact_iso_abundances = Float32[0.411367694,0.852177261,1.000000000,0.851991166,0.581595404,0.335087683]
spline_iso_abundances = [iso_splines.splines[pep_sulfur_count][x](pep_mass_a) for x in range(1, 6)]
spline_iso_abundances = spline_iso_abundances./maximum(spline_iso_abundances)
percent_deviance = 100.0*(exact_iso_abundances .- spline_iso_abundances)./exact_iso_abundances