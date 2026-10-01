# Add Koina ion-mobility predictions (ccs, inv_ion_mobility) to an existing .poin library's precursors_table.arrow,
# in place, with the same functions BuildSpecLib's im_model step uses. Used to give the committed E. coli test library
# (test/integration/ecoli_lib.poin) the 1/K0 column that searching the timsTOF .d fixture needs.
# Usage: julia --project=<Pioneer> add_im_to_lib.jl <lib.poin> [alphapept_ccs|im2deep]
using Pioneer, Arrow, DataFrames, Tables
using Pioneer: predict_ccs_koina, ccs_to_inv_ion_mobility
lib = ARGS[1]
im_model = length(ARGS) >= 2 ? ARGS[2] : "alphapept_ccs"
p = joinpath(lib, "precursors_table.arrow")
t = DataFrame(Tables.columntable(Arrow.Table(p)))
println(nrow(t), " precursors; columns: ", join(names(t), ", "))
hasproperty(t, :koina_sequence) || error("precursors_table.arrow has no koina_sequence column")
if hasproperty(t, :inv_ion_mobility)
    println("library already has inv_ion_mobility; nothing to do"); exit(0)
end
t[!, :precursor_charge] = t.prec_charge          # column name prepare_koina_batch reads
t0 = time()
ccs = predict_ccs_koina(t, im_model)
println("Koina ", im_model, ": ", length(ccs), " predictions in ", round(time() - t0, digits = 1), " s")
select!(t, Not(:precursor_charge))
t[!, :ccs] = Float32.(ccs)
t[!, :inv_ion_mobility] = Float32.(ccs_to_inv_ion_mobility.(ccs, t.prec_charge, t.mz))
tmp = p * ".new"
Arrow.write(tmp, t)
mv(tmp, p; force = true)
chk = Arrow.Table(p)
println("wrote ", p, ": ", length(chk.inv_ion_mobility), " rows; 1/K0 range ",
        extrema(chk.inv_ion_mobility), "; ccs range ", extrema(chk.ccs))
