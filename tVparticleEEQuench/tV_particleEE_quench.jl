# Entanglement entropy of the t-V-V` Model at half-filling - dependence on interaction strength V

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "src"))

using Pkg

Pkg.add("Formatting")
Pkg.add("Arpack")
Pkg.add("ArgParse")
Pkg.add("ProgressBars")
Pkg.add("BenchmarkTools")
Pkg.add("MKL")


using MKL
using tVVpDiagonalize
using ArgParse
using IntFermionicbasis
using Printf
using ProgressBars
using Dates
using Formatting
using QuenchCalculation
using OutputFileHandler
using Utils

# This file contains specific file output operations which we moved
# out of the main script to improve readability
include("./src/file_operations.jl")
# This file contains specific parameter checks which we moved
# out of the main script to improve readability
# In this file you find the list of accepted parameters or by
# $julia typing tV_particleEE_quench.jl --help
# into a terminal.
include("./src/parameter_operations.jl")



"""Runs entanglement calculation with quench for a range of interaction strenths V based on the input parameters.

"""
function main()
    # _____________1_Parameter_Setup________________
    c = parse_commandline()
    M, N, t, V0, V_array, Vp0, Vp, boundary, Asize, ℓsize, time_range, time_num, Δt = unpack_parameters(c)

    # _____________2_Output_Setup___________________
    file_handler = create_file_output_handler_and_write_headers(c, M, N, t, V0, V_array, Vp0, Vp, Asize, ℓsize, time_range, Δt)

    # _____________3_Calculation______________________ 
    # compute symmetry cycles 
    Cycle_leaders, Cycle_sizes, NumOfCycles = symmetry_cycles_q0R1PH1(M, N)
    # need structure matrix only when computing entanglement not when saving states only
    if ~c[:save_states]
        # Compute structure matrix 
        AmatrixStructure = PE_StructureMatrix(M, N, Cycle_leaders, Asize)
    end
    # get equilibrium ground state only if not loading states
    if ~c[:load_states]
        # save off-diagonal terms once for V0=0=V0'
        if ~c[:skip_hoffdiag_saving]
            save_offdiag_sparse_Block_Diagonal_Hamiltonian_q0R1PH1(M, N, Cycle_leaders, Cycle_sizes, NumOfCycles, t, c[:tmp])
        end
        # compute ground state Ψ and devide by cycle size below to get Ψ_coeff (only use single variable Ψ_coeff here)
        Ψ0, HRank = ground_state(M, N, Cycle_leaders, Cycle_sizes, NumOfCycles, t, V0, Vp0, boundary, ~c[:skip_hoffdiag_loading], c[:tmp])
    end

    for V in ProgressBar(V_array)
        if c[:save_states]
            # state saving caluclation
            quench_state_prep(
                file_handler, M, N, V0, Vp0, V, Vp, t, time_num,
                Cycle_leaders, Cycle_sizes, NumOfCycles,
                Ψ0, time_range;
                load_offdiagonal=~c[:skip_hoffdiag_loading],
                ham_offdiag_storage_path=c[:tmp])
        else
            # file headers 
            write_headers_to_files(file_handler, V, c, M, N, t, V0, Vp0, Vp, boundary, Asize, ℓsize, time_range, Δt)
            # load the states from disk
            times, Ψts, HRank = load_states(
                get_file(file_handler, "states"; open_as="r", replace_list=[V]),
                M, N, Δt, V0, V, Vp0, Vp, time_range, NumOfCycles, time_num)
            # perform actual calculation
            quench_entanglement_calc(
                file_handler, M, N, V, Vp, Asize, ℓsize, AmatrixStructure,
                Cycle_leaders, Cycle_sizes, HRank, Ψts, times;
                calc_obdm=c[:obdm],
                calc_spatial=c[:spatial],
                calc_g2=c[:g2])
            # close files
            close_file_handler(file_handler, c, V)
        end

    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

