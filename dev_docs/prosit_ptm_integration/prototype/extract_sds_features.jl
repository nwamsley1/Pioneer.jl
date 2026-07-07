# Site-determining-ion feature extractor on the COMBINED search results (prototype).
# S_vs_targets per passing C>1 instance: target vs sibling targets; decoy vs ALL targets incl parent.
# Distinguishing = fragments whose cumulative phospho count differs; support matched to observed peaks.
# Usage: julia --project=Pioneer extract_sds_features.jl <results_dir> <lib.poin> [n_sample] [out.arrow]
using Arrow, DataFrames, Statistics, Combinatorics, Random
using DataFrames: combine, groupby

const PROTON=1.0072764668; const H2O=18.0105646863; const PHOS=79.96633; const CAM=57.02146
const AA = Dict{Char,Float64}('A'=>71.03711,'R'=>156.10111,'N'=>114.04293,'D'=>115.02694,'C'=>103.00919,
 'E'=>129.04259,'Q'=>128.05858,'G'=>57.02146,'H'=>137.05891,'I'=>113.08406,'L'=>113.08406,'K'=>128.09496,
 'M'=>131.04049,'F'=>147.06841,'P'=>97.05276,'S'=>87.03203,'T'=>101.04768,'W'=>186.07931,'Y'=>163.06333,'V'=>99.06841)

phpos(m)  = Int[parse(Int,x.captures[1]) for x in eachmatch(r"\((\d+),[A-Z],Unimod:21\)",String(m))]
campos(m) = Int[parse(Int,x.captures[1]) for x in eachmatch(r"\((\d+),[A-Z],Unimod:4\)", String(m))]

@inline function frag_mz(res::Vector{Float64}, ph::Vector{Int}, cm::Vector{Int}, T::Symbol, i::Int, z::Int, L::Int)
    if T===:b
        m=0.0; @inbounds for k in 1:i; m+=res[k]; end
        m += PHOS*count(<=(i),ph) + CAM*count(<=(i),cm)
    else
        m=H2O; @inbounds for k in (L-i+1):L; m+=res[k]; end
        m += PHOS*count(>=(L-i+1),ph) + CAM*count(>=(L-i+1),cm)
    end
    (m + z*PROTON)/z
end

@inline function matched_int(mz, omz, oint, tol)
    d = mz*tol*1e-6; lo=searchsortedfirst(omz,mz-d); hi=searchsortedlast(omz,mz+d)
    b=0.0; @inbounds for k in lo:hi; b=max(b,Float64(oint[k])); end; b
end

# s(inst, comp): distinguishing fragments (cumulative phospho count differs), support at each m/z.
function s_pair_rd(res, phI, phC, cm, L, omz, oint, tol; charges=(1,2))
    supI=0.0; supC=0.0; nd=0
    for T in (:b,:y), i in 1:(L-1)
        cI = T===:b ? count(<=(i),phI) : count(>=(L-i+1),phI)
        cC = T===:b ? count(<=(i),phC) : count(>=(L-i+1),phC)
        cI==cC && continue
        for z in charges
            supI += matched_int(frag_mz(res,phI,cm,T,i,z,L),omz,oint,tol)
            supC += matched_int(frag_mz(res,phC,cm,T,i,z,L),omz,oint,tol); nd+=1
        end
    end
    cap=supI+supC; (s = cap==0 ? NaN : supI/cap, nd=nd, cap=cap)
end
# fraction of the instance's own b/y (z=1,2) fragments matched -- sanity that it's really in the scan
function self_frac(res, ph, cm, L, omz, oint, tol; charges=(1,2))
    m=0; t=0
    for T in (:b,:y), i in 1:(L-1), z in charges
        t+=1; matched_int(frag_mz(res,ph,cm,T,i,z,L),omz,oint,tol)>0 && (m+=1)
    end
    m/max(t,1)
end

function main()
    resdir=ARGS[1]; lib=ARGS[2]
    nsample = length(ARGS)>=3 ? parse(Int,ARGS[3]) : 4000
    outp = length(ARGS)>=4 ? ARGS[4] : joinpath(resdir,"sds_features.arrow")
    tol=20.0; arrow_dir="/Users/n.t.wamsley/prosit_phospho_test/coon_dilution/arrow_out"

    lp=DataFrame(Arrow.Table(joinpath(lib,"precursors_table.arrow"))); n=nrow(lp)
    seqs=String.(lp.sequence); mods=String.(coalesce.(lp.structural_mods,"")); chg=collect(lp.prec_charge)
    ild=collect(lp.is_loc_decoy); isd=collect(lp.is_decoy)
    ph=[phpos(m) for m in mods]; cm=[campos(m) for m in mods]
    nSTY=[count(c->c=='S'||c=='T'||c=='Y',s) for s in seqs]
    Cn=[(0<length(ph[i])<=nSTY[i]) ? binomial(nSTY[i],length(ph[i])) : 1 for i in 1:n]
    grpT=Dict{Tuple{String,UInt8},Vector{Int}}()
    for i in 1:n; (isd[i]||ild[i]) && continue; push!(get!(grpT,(seqs[i],chg[i]),Int[]),i); end

    d=DataFrame(Arrow.Table(joinpath(resdir,"precursors_long.arrow")))
    d.ild=[ild[Int(i)] for i in d.precursor_idx]; d.C=[Cn[Int(i)] for i in d.precursor_idx]
    d=d[(d.target.==true).&(d.qval.<=0.01f0).&(coalesce.(d.global_qval,1f0).<=0.01f0).&(d.C.>1),:]
    d.frac=Float64.(coalesce.(d.iso_weight_fraction_at_scan,0.0))
    d.gsz = Int.(coalesce.(d.iso_group_size_at_scan,1))
    pf=combine(groupby(d,:precursor_idx)) do g
        j=argmax(g.frac); (scan=UInt32(g.scan_idx[j]),file=String(g.file_name[j]),ild=g.ild[j],C=g.C[j],
                           frac=Float64(g.frac[j]), gsz=g.gsz[j])
    end
    pf=pf[shuffle(MersenneTwister(1), 1:nrow(pf)), :]        # unbiased sample
    pf=pf[1:min(nsample,nrow(pf)),:]
    println("scoring ",nrow(pf)," instances (",count(pf.ild)," decoys)")

    out=DataFrame(precursor_idx=UInt32[],is_loc_decoy=Bool[],C=Int[],S_vs_targets=Float64[],
                  n_comp=Int[],self_frac=Float64[],frac=Float64[],gsz=Int[])
    for gi in groupby(pf,:file)
        af=joinpath(arrow_dir,"20240306_NML_15min_2Th_250ngYeast_"*gi.file[1]*".arrow")
        isfile(af) || (println("  missing ",af); continue)
        t=Arrow.Table(af); MZ=t.mz_array; IN=t.intensity_array
        for row in eachrow(gi)
            pid=Int(row.precursor_idx); sc=Int(row.scan); s=seqs[pid]; L=length(s); res=[AA[c] for c in s]
            omz=Float64.(MZ[sc]); oint=Float64.(IN[sc])
            p=sortperm(omz); omz=omz[p]; oint=oint[p]
            comps = row.ild ? grpT[(s,chg[pid])] : filter(!=(pid), grpT[(s,chg[pid])])
            smin=Inf; anyd=false
            for c in comps
                r=s_pair_rd(res, ph[pid], ph[c], cm[pid], L, omz, oint, tol)
                isnan(r.s) && continue
                anyd=true; smin=min(smin,r.s)
            end
            sf=self_frac(res, ph[pid], cm[pid], L, omz, oint, tol)
            push!(out,(UInt32(pid),row.ild,row.C, anyd ? smin : NaN, length(comps), sf, row.frac, row.gsz))
        end
    end
    Arrow.write(outp, out)
    println("wrote ",nrow(out)," rows -> ",outp)
    # ---- sanity + separation ----
    o=out[.!isnan.(out.S_vs_targets),:]
    for lbl in (false,true)
        s=o[o.is_loc_decoy.==lbl,:]
        println(lbl ? "DECOY " : "TARGET", " n=",nrow(s)," S median=",round(median(s.S_vs_targets),digits=3),
                " mean=",round(mean(s.S_vs_targets),digits=3)," frac(S>=0.9)=",round(count(>=(0.9),s.S_vs_targets)/max(nrow(s),1),digits=3))
    end
    cb(C)= C==2 ? "2" : C==3 ? "3" : C<=6 ? "4-6" : "7+"
    o.cb=cb.(o.C)
    println("\nS_vs_targets by C stratum (target median | decoy median):")
    for b in ["2","3","4-6","7+"]
        tt=o[(o.cb.==b).&(.!o.is_loc_decoy),:]; dd=o[(o.cb.==b).&(o.is_loc_decoy),:]
        println("  C",b,": target=",isempty(tt) ? "-" : round(median(tt.S_vs_targets),digits=2),
                "  decoy=",isempty(dd) ? "-" : round(median(dd.S_vs_targets),digits=2),
                "  (nt=",nrow(tt),",nd=",nrow(dd),")")
    end
    println("\nFLR = decoy/target at S_vs_targets cutoff (this SAMPLE; not FDR-scaled):")
    R=count(!,o.is_loc_decoy); D=count(o.is_loc_decoy)
    for tc in (0.0,0.5,0.75,0.9,0.95)
        s=o[o.S_vs_targets.>=tc,:]; r=count(!,s.is_loc_decoy); d2=count(s.is_loc_decoy)
        println("  S>=",tc,": target IDs=",r,"  decoy=",d2,"  D/R=",round(100*d2/max(r,1),digits=1),"%")
    end
end
main()
