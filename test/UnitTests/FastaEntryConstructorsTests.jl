# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

using Test

@testset "FastaEntry constructors (compat)" begin
    # Inputs as AbstractString (SubString) and Integer for entrapment id
    acc = SubString("sp|P12345|TEST", 1:12)
    desc = SubString("some description", 1:5)
    gene = SubString("GENE1", 1:5)
    prot = SubString("PROT1", 1:5)
    org = SubString("human", 1:5)
    protm = SubString("human", 1:5)
    seq = SubString("PEPTIDEK", 1:8)

    # 15-arg canonical constructor
    e1 = FastaEntry(acc, desc, gene, prot, org, protm, seq,
                    UInt32(1), missing, missing, UInt8(0),
                    UInt32(0), UInt32(0), 0, false)
    @test get_id(e1) == String(acc)
    @test get_description(e1) == String(desc)
    @test get_gene(e1) == String(gene)
    @test get_protein(e1) == String(prot)
    @test get_organism(e1) == String(org)
    @test get_proteome(e1) == String(protm)
    @test get_sequence(e1) == String(seq)
    @test get_entrapment_pair_id(e1) == 0
    @test get_base_pep_id(e1) == 0
    @test get_charge(e1) == 0

    # 16-arg compatibility constructor (base_prec_id ignored)
    e2 = FastaEntry(acc, desc, gene, prot, org, protm, seq,
                    UInt32(1), missing, missing, UInt8(0),
                    UInt32(1), UInt32(2), UInt32(999), 1, false)
    @test get_base_pep_id(e2) == 2
    @test get_entrapment_pair_id(e2) == 1
end

