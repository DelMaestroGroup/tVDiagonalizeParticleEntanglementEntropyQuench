""" Definitions of the unit tests run test by running "runtests.jl" """
# provide additional print output during tests
const VERBOSE = false

# import modules used in tests  
push!(LOAD_PATH, joinpath(dirname(@__FILE__), "..", "src"))
using Test

using tVVpDiagonalize
using ArgParse
using IntFermionicbasis
using Printf
using LinearAlgebra
using OutputFileHandler
using Formatting
using QuenchCalculation

test_save_location = "./test/"

# test helper functions
function test_approx(got, truth, atol, verbose=true)
    passed = isapprox(got, truth, atol=atol)
    if verbose
        if passed
            print("\nPASS: ")
        else
            print("\nFAIL: ")
        end
        print("got: ", got, ", truth: ", truth, '\n')
    end
    @test passed
end

# define test cases
@time @testset "Test Entropy and Pair Correlations" begin

    """
    Runs test cases for spatial entanglement entropy and particle entanglement entropy
    comparing with S1_spatial_truth, S2_spatial_truth, S1_particle_truth, S2_particle_truth
    allowing for a tolerance tol.
    """
    function tests_entropies(t::Float64, V::Float64, Vp::Float64, M::Int64, N::Int64, n::Int64, l::Int64, time::Float64, S1_spatial_truth::Float64, S2_spatial_truth::Float64, S1_particle_truth::Float64, S2_particle_truth::Float64, boundary::BdryCond=PBC, tol::Float64=0.0, load_offdiag::Bool=false, save_offdiag::Bool=false; verbose=VERBOSE)
        V0 = 0.0
        Vp0 = 0.0
        Δt = 0.1
        time_num = 2

        # compute symmetry cycles 
        Cycle_leaders, Cycle_sizes, NumOfCycles = symmetry_cycles_q0R1PH1(M, N)
        # Compute structure matrix 
        AmatrixStructure = PE_StructureMatrix(M, N, Cycle_leaders, n)
        # get equilibrium ground state only if not loading states 
        if save_offdiag
            save_offdiag_sparse_Block_Diagonal_Hamiltonian_q0R1PH1(M, N, Cycle_leaders, Cycle_sizes, NumOfCycles, t, "./test")
        end
        # compute ground state Ψ and devide by cycle size below to get Ψ_coeff (only use single variable Ψ_coeff here)
        Ψ0, HRank = ground_state(M, N, Cycle_leaders, Cycle_sizes, NumOfCycles, t, V0, Vp0, boundary, load_offdiag, "./test")


        # state saving caluclation
        file_handler = FileOutputHandler(true)
        handler_name = "states"
        Ψt_output = @sprintf "psioft_%02d_%02d_%+5.3f_%+5.3f_%+5.3f_%+5.3f_%6.4f_%06.3f_%06.3f.dat" M N V0 V Vp0 Vp Δt time - 0.1 time
        # we handle writing Ψ with a different function and only grab the file using the handler
        data_to_str = x -> nothing
        # add to handler, open later
        path = joinpath("./test/", format(Ψt_output, V))
        add!(file_handler, path, data_to_str, handler_name)

        time_range = [time - 0.1, time]
        quench_state_prep(
            file_handler, M, N, V0, Vp0, V, Vp, t, time_num,
            Cycle_leaders, Cycle_sizes, NumOfCycles,
            Ψ0, time_range;
            load_offdiagonal=load_offdiag,
            ham_offdiag_storage_path="./test")
        close(file_handler.files_dict["states"])
        # load the states from disk
        times, Ψts, HRank = load_states(open(file_handler.paths_dict["states"], "r"),
            M, N, Δt, V0, V, Vp0, Vp, time_range, NumOfCycles, time_num)

        Ψ_coeff = zeros(ComplexF64, HRank)
        it = 2
        Ψ_coeff .= Ψts[:, it] ./ sqrt(dot(Ψts[:, it], Ψts[:, it]))
        for j = 1:HRank
            Ψ_coeff[j] = Ψ_coeff[j] / sqrt(Cycle_sizes[j])
        end
        # 3.1 particle entanglement  
        s_particle = particle_entropy_Ts(M, N, n, Ψ_coeff, false, AmatrixStructure)
        s_spatial = spatial_entropy(M, N, l, Ψ_coeff, Cycle_leaders)


        # test against production code result 
        if verbose
            print('\n')
        end
        test_approx(s_spatial[1], S1_spatial_truth, tol, verbose)
        test_approx(s_spatial[2], S2_spatial_truth, tol, verbose)
        test_approx(s_particle[1], S1_particle_truth, tol, verbose)
        test_approx(s_particle[2], S2_particle_truth, tol, verbose)
        if verbose
            print('\n')
        end
    end

    testcases = Vector()

    # test case 1:
    params = (
        t=1.0,
        V=1.3,
        Vp=0.4,
        M=12,
        N=6,
        n=1,
        l=6,
        time=3.4,
        S1_spatial_truth=1.2150198462585413,
        S2_spatial_truth=0.8684233157776676,
        S1_particle_truth=0.13803439950396768,
        S2_particle_truth=0.0694216265492471,
        boundary=PBC,
        tol=0.0)
    push!(testcases, params)
    # test case 2 
    params = (
        t=1.0,
        V=0.3,
        Vp=1.3,
        M=14,
        N=7,
        n=3,
        l=7,
        time=1.4,
        S1_spatial_truth=1.4908521461808377,
        S2_spatial_truth=1.1094295769295803,
        S1_particle_truth=0.36124135391973633,
        S2_particle_truth=0.15310766224311756,
        boundary=PBC,
        tol=0.0)
    push!(testcases, params)
    # test case 3
    params = (
        t=1.0,
        V=0.003,
        Vp=0.0,
        M=14,
        N=7,
        n=2,
        l=7,
        time=0.1,
        S1_spatial_truth=1.2261438388174952,
        S2_spatial_truth=0.9783657471747936,
        S1_particle_truth=5.165909851001516e-7,
        S2_particle_truth=5.4591921649915776e-8,
        boundary=PBC,
        tol=0.0)
    push!(testcases, params)
    # test case 4
    params = (
        t=1.0,
        V=100.0,
        Vp=0.0,
        M=14,
        N=7,
        n=5,
        l=7,
        time=0.3,
        S1_spatial_truth=1.8593409386588757,
        S2_spatial_truth=1.4653906531261338,
        S1_particle_truth=1.0698388451335137,
        S2_particle_truth=0.8798080908374009,
        boundary=PBC,
        tol=0.0)
    push!(testcases, params)
    # test case 5
    params = (
        t=1.0,
        V=0.003,
        Vp=0.0,
        M=14,
        N=7,
        n=1,
        l=7,
        time=10.0,
        S1_spatial_truth=1.2261456519273655,
        S2_spatial_truth=0.9783245045831165,
        S1_particle_truth=4.715608179717279e-6,
        S2_particle_truth=5.996319938361694e-7,
        boundary=PBC,
        tol=0.0,
        load_offdiag=false,
        save_offdiag=true)
    push!(testcases, params)
    # test case 6
    params = (
        t=0.1,
        V=0.3,
        Vp=0.0,
        M=14,
        N=7,
        n=1,
        l=7,
        time=2.1,
        S1_spatial_truth=1.3203934228333842,
        S2_spatial_truth=1.0336533957099077,
        S1_particle_truth=0.23391771915478965,
        S2_particle_truth=0.12771956247939875,
        boundary=PBC,
        tol=0.0,
        load_offdiag=true,
        save_offdiag=false)
    push!(testcases, params)
    # test case 7
    params = (
        t=1.0,
        V=0.003,
        Vp=0.0,
        M=14,
        N=7,
        n=1,
        l=7,
        time=10.0,
        S1_spatial_truth=1.2261456519273655,
        S2_spatial_truth=0.9783245045831165,
        S1_particle_truth=4.715608179717279e-6,
        S2_particle_truth=5.996319938361694e-7,
        boundary=PBC,
        tol=0.0,
        load_offdiag=true,
        save_offdiag=false)
    push!(testcases, params)
    # test case 8
    params = (
        t=2.3,
        V=0.3,
        Vp=1.3,
        M=14,
        N=7,
        n=3,
        l=7,
        time=2.0,
        S1_spatial_truth=1.301591509038949,
        S2_spatial_truth=1.0163356357115163,
        S1_particle_truth=0.09492374989390351,
        S2_particle_truth=0.030402812888033814,
        boundary=PBC,
        tol=0.0)
    push!(testcases, params)

    # run tests 
    @testset "Test case $i" for i in 1:length(testcases)
        tests_entropies(testcases[i]...)
    end
end

# define test cases
@time @testset "Test Serial Number" begin

    """
    Runs test cases for the serial_num function.
    """
    basis = Fermionsbasis(12, 6)
    @test serial_num(12, 6, 63) == 1
    @test serial_num(basis, 63) == 1
    @test serial_num(12, 6, 95) == 2
    @test serial_num(basis, 95) == 2
    @test serial_num(12, 6, 119) == 4
    @test serial_num(basis, 119) == 4
    basis = Fermionsbasis(30, 15)
    @test serial_num(30, 15, 1073680384) == 155117517
    @test serial_num(basis, 1073680384) == 155117517
    @test serial_num(30, 15, 65471) == 10
    @test serial_num(basis, 65471) == 10
end