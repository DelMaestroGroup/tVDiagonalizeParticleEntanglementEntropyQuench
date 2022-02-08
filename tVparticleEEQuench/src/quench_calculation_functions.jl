struct FileHeader
    M::Int64
    N::Int64
    time_num:: Int64
    basis_num:: Int64
    V0::Float64
    V::Float64
    Vp0::Float64
    Vp::Float64
    time_min::Float64
    time_max::Float64
end

"""
Compute the time evolved states and save them to file.
"""
function quench_state_prep(
        file_handler::FileOutputHandler, 
        M::Int64,
        N::Int64,
        V0::Float64,
        Vp0::Float64,
        V::Float64,
        Vp::Float64, 
        t::Float64,
        time_num::Int64,
        Cycle_leaders::Vector{Int64},
        Cycle_sizes:: Vector{Int64},
        NumOfCycles::Int64,
        Ψ0::Vector{ComplexF64},
        times::Vector{Float64};
        load_offdiagonal::Bool=false,
        ham_offdiag_storage_path::String="./") 

    EigenEnergie=ComplexF64
    #  wave function in terms of the Symmetry basis at time t
    Ψt=zeros(ComplexF64, NumOfCycles, time_num)
    TimeEvolutionFactor=zeros(ComplexF64, time_num)
     
    H,HRank = sparse_Block_Diagonal_Hamiltonian_q0R1PH1(M, N ,Cycle_leaders, Cycle_sizes, NumOfCycles, t, V, Vp, load_offdiagonal, ham_offdiag_storage_path) 
    EigenEnergies,VV = eigen(Symmetric(Matrix(H))) #full(H)
    H= Nothing
    for i_HRank = 1:HRank
        EigenEnergie= EigenEnergies[i_HRank]
        Ψn = dot(VV[:, i_HRank],Ψ0)
        TimeEvolutionFactor .= exp.(-(0.0+1.0im) * times * EigenEnergie)*Ψn
        Ψt += kron(VV[:, i_HRank],transpose(TimeEvolutionFactor))
    end
    VV= Nothing
    EigenEnergies= Nothing 
     
    # save states
    Ψtf = get_file(file_handler,"states";open_as="w",replace_list=[V])
    file_header= FileHeader(M, N, time_num, HRank, V0, V, Vp0, Vp, times[1], times[end])
    write(Ψtf, file_header.M, file_header.N, file_header.time_num, file_header.basis_num, file_header.V0, file_header.V, file_header.Vp0, file_header.Vp ,file_header.time_min, file_header.time_max,Ψt)
    flush(Ψtf)
    close(Ψtf)

    return Ψt
end

"""
Load time evolved states from file and compute entanglement.
"""
function quench_entanglement_calc(
        file_handler::FileOutputHandler,
        M::Int64,
        N::Int64,
        V::Float64,
        Vp::Float64,
        Asize::Int64,
        ℓsize::Int64,
        AmatrixStructure::AbstractArray,
        Cycle_leaders::Vector{Int64},
        CycleSize::Vector{Int64},
        HRank::Int64,
        Ψts::Matrix{ComplexF64},
        times::Vector{Float64};
        calc_obdm::Bool=false,
        calc_spatial::Bool=false,
        calc_g2::Bool=false,
        out_folder_obdm::String="./")

    Ψ_coeff = zeros(ComplexF64, HRank)
    if calc_obdm && Asize==1
        obdm = zeros(Float64,M,length(times))
    end

    # time loop
    for (it,time) in enumerate(times)
        Ψ_coeff .= Ψts[:,it]./sqrt(dot(Ψts[:,it],Ψts[:,it])) 
        for j=1:HRank
            Ψ_coeff[j]=Ψ_coeff[j]/sqrt(CycleSize[j])
        end
        # 3.1 particle entanglement 
        if calc_obdm && Asize==1  
            s_particle, obdm[:,it] = particle_entropy_Ts(M, N, Asize, Ψ_coeff, true, AmatrixStructure) 
        else
            s_particle = particle_entropy_Ts(M, N, Asize, Ψ_coeff, false, AmatrixStructure) 
        end 
        # save to file
        write(file_handler,"particleEE",time,s_particle; replace_list=[V]) 

        # 3.2 calculate spatial entanglement 
        if calc_spatial
            s_spatial = spatial_entropy(M, N, ℓsize, Ψ_coeff, Cycle_leaders) 
            # save to file
            write(file_handler,"spatialEE",time,s_spatial; replace_list=[V])  
        end

        # 3.3 pair correlation function
        if calc_g2
            g2 = pair_correlation(M, N, Ψ_coeff, Cycle_leaders)  
            # save to file 
            write(file_handler,"g2",time,g2; replace_list=[V])   
        end 

    end # end of time loop 
    if calc_obdm && Asize==1 
        # save obdm to file (one file for each V)
        save_obdm(file_handler,obdm,M,V,times) 
    end
 
    return nothing
end


function load_states(Ψtf::IOStream,M::Int64,N::Int64,Δt::Float64,V0::Float64,V::Float64,Vp0::Float64,Vp::Float64,time_range::Vector{Float64},NumOfCycles::Int64,time_num::Int64)
    file_header1 =zeros(Int64,4)
    file_header2 =zeros(Float64,6) 
    read!(Ψtf,file_header1)
    read!(Ψtf,file_header2)
    M_f=file_header1[1]
    N_f=file_header1[2]
    time_num_f=file_header1[3]
    basis_num_f=file_header1[4]
    V0_f=file_header2[1]
    V_f=file_header2[2]
    Vp0_f=file_header2[3]
    Vp_f=file_header2[4]
    time_min_f=file_header2[5]
    time_max_f=file_header2[6]
    time_range_f = lin_range(time_min_f, time_max_f, time_num_f)
    #println(length(time_range_f))
    if length(time_range_f) > 1
        Δt_f = time_range_f[2]-time_range_f[1]
    else
        Δt_f = time_range_f[1]
    end
    #println(Δt_f)
    if (M_f!=M) || (N_f!=N) || (abs(Δt_f-Δt)> 1.0E-12)||(abs(V0_f- V0)> 1.0E-12) ||(abs(V_f- V)> 1.0E-12)||(abs(Vp0_f- Vp0)> 1.0E-12)||(abs(Vp_f- Vp)> 1.0E-12) ||((time_range_f[1]- time_range[1])> 1.0E-12) ||(time_range_f[end]- time_range[end]< -1.0E-12)
        println("the file of states is not compatible with the input parameters" )
        println("M=",M," N=",N," V0=",V0," V=",V," Vp0=",Vp0," Vp=",Vp," Δt=",Δt," time_max=", time_range[end]," time_min=",time_min)
        println("M_f=",M_f," N_f=",N_f," V0_f=",V0_f," V_f=",V_f," Vp0_f=",Vp0_f," Vp_f=",Vp_f," Δt_f=",Δt_f," time_max_f=", time_max_f," time_min_f=",time_min_f)
        exit(1)
    end
    Ψt = zeros(ComplexF64, NumOfCycles, time_num)
    Ψ=zeros(ComplexF64, NumOfCycles)
    it=0
    for (it_f, time_f) in enumerate(time_range_f)
        read!(Ψtf, Ψ)
        if (abs(time_f- time_range[it+1])< 1.0E-12)
            it+=1
            Ψt[:,it]=Ψ[:]
            if it== time_num
                break
            end
        end
    end
   close(Ψtf)
   HRank = basis_num_f

   return time_range_f, Ψt, HRank
end