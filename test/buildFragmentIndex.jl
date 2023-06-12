file_path = "/Users/n.t.wamsley/RIS_temp/HAMAD_MAY23/mouse_SIL_List/UP000000589_10090.fasta.gz"
fixed_mods = [(p=r"C", r="C[Carb]")]
var_mods = [(p=r"(K$)", r="[Hlys]"), (p=r"(R$)", r="[Harg]"), (p=r"(M)", r="[MOx]")]
test_table = PrecursorTable()
testFASTA = [
FastaEntry("A", "", "PEPTIDEK", false),
FastaEntry("A", "", "PEPTIDEK", false),
FastaEntry("A", "", "DRAGRACER", false),
FastaEntry("B", "", "PEPTIDEK", false),
FastaEntry("B", "", "DENTIST", false),
FastaEntry("C", "", "DENTISTRY", false)]