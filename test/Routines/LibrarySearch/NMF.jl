@testset "NMF.jl" begin

    #Show that our sparse NMF gives the same answers using coordinate Descent

    X = Float32[2 3 3 1 4 4 3]
    W = Float32[1 1 1]
    H = Matrix(Float32[1 0 0; 1 1 0; 1 1 0; 0 1 0; 0 1 1; 0 1 1; 0 0 1]')

    NMF.solve!(NMF.CoordinateDescent{Float32}(
                                            maxiter = 100,
                                            update_H = false,
                                            tol = 1.0e-8,
                                            α = 0.0,
                                            regularization = :transformation,
                                            l₁ratio=1.0

                ), X, W, H)

    @test sum(W.-[2, 1, 3]).<1e-7



    X = Float32[2 3 3 1 4 4 3][:]
    H = sparse(Float32[1 0 0; 1 1 0; 1 1 0; 0 1 0; 0 1 1; 0 1 1; 0 0 1]')

    #Prove that our solver gets the same answer. 
    @test sum(sparseNMF(H, sparse(H'), X, tol = Float32(1e-8)).-W).<1e-7

    ##########
    #Try again but set the lasso coefficients
    ##########
    X = Float32[2 3 3 1 4 4 3]
    W = Float32[1 1 1]
    H = Matrix(Float32[1 0 0; 1 1 0; 1 1 0; 0 1 0; 0 1 1; 0 1 1; 0 0 1]')

    NMF.solve!(NMF.CoordinateDescent{Float32}(
                                            maxiter = 100,
                                            update_H = false,
                                            tol = 1.0e-8,
                                            α = 10.0,
                                            regularization = :transformation,
                                            l₁ratio=1.0

                ), X, W, H)

    #@test sum(W.-[2, 1, 3]).<1e-7



    X = Float32[2 3 3 1 4 4 3][:]
    H = sparse(Float32[1 0 0; 1 1 0; 1 1 0; 0 1 0; 0 1 1; 0 1 1; 0 0 1]')
    sparseNMF(H, sparse(H'), X, λ=Float32(10.0), tol = Float32(1e-8))
    #Prove that our solver gets the same answer. 
    @test abs(sum(sparseNMF(H, sparse(H'), X, λ=Float32(10.0/sqrt(size(H)[1])), tol = Float32(1e-8)).-W)).<1e-7

end