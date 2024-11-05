function getCVFolds(
                prec_ids::AbstractVector{UInt32},
                protein_groups::AbstractVector{String}
                )
    uniq_pgs = unique(protein_groups) #Unique protein groups
    pg_to_cv_fold = Dictionary{String, UInt8}() #Protein group to cross-validation fold
    #Map protein groups to cv folds 
    for pg in uniq_pgs
        insert!(pg_to_cv_fold,
        pg,
        rand(UInt8[0, 1])
        )
    end
    #Now map pid's to cv folds 
    pid_to_cv_fold = Dictionary{UInt32, UInt8}()
    for pid in prec_ids
        insert!(
        pid_to_cv_fold,
        pid,
        pg_to_cv_fold[protein_groups[pid]]
        )
    end
    return pid_to_cv_fold
end
