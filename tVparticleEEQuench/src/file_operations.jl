"""For readability of the main script, the file setup
steps are moved to this .jl file that is included
into the main script"""

function create_file_output_handler_and_write_headers(
        c::Dict{Symbol,Any},
        M::Int64,
        N::Int64,
        t::Float64,
        V0::Float64,
        V_array::Vector{Float64},
        Vp0::Float64,
        Vp::Float64, 
        Asize::Int64,
        ℓsize::Int64,
        time_range::Vector{Float64},
        Δt::Float64)

    if c[:out] === nothing
        out_folder = "./"
    else
        out_folder = c[:out]
    end 
    if c[:tmp] === nothing
        tmp_folder = "./"
    else
        tmp_folder = c[:tmp]
    end 
    if c[:out_obdm] === nothing
        out_folder_obdm = out_folder
    else
        out_folder_obdm = c[:out_obdm]
    end
    calculation_label = @sprintf "M%02d_N%02d_t%+5.3f_Vp%+5.3f_Vp0%+5.3f_V{:+5.3f}_V0%+5.3f_dt%6.4f_tstart%06.3f_tendf_%06.3f" M N t Vp Vp0 V0 Δt time_range[1] time_range[end]
    # create file output handler
    file_handler = FileOutputHandler(~c[:no_flush])
    # 2.1. output of particle entanglement (pe_01)
    if ~c[:save_states]
        handler_name = "particleEE"
        hname_format = "_V{:5.3f}"
        # function to convert data to string
        out_str_pe_01 = (t,data)->@sprintf "%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E\n" t data...
        # add file series to handler (files will be opend when first written to)
        for V in V_array
            path = joinpath(out_folder,@sprintf "particle_entanglement_n%02d_%s.dat" Asize format(calculation_label,V))  
            add!(file_handler,path,out_str_pe_01,handler_name;hname_format=hname_format,replace_list=[V]) 
        end
    end
    # 2.2. output of spatial entanglement (se_02)
    if c[:spatial] && ~c[:save_states]
        handler_name = "spatialEE"
        hname_format = "_V{:5.3f}"  
        # function to convert data to string
        out_str_se_02 = (t,data)->@sprintf "%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E\n" t data...
        # add to handler, open later
        for V in V_array
            path = joinpath(out_folder,@sprintf "spatial_entanglement_l%02d_%s.dat" ℓsize format(calculation_label,V))
            add!(file_handler,path,out_str_se_02,handler_name;hname_format=hname_format,replace_list=[V]) 
        end
    end

    # 2.3. output of pair correlations (pcf_03)
    if c[:g2] && ~c[:save_states]
        handler_name = "g2"
        hname_format = "_V{:5.3f}"  
        # function to convert data to string
        function out_str_pcf_03(t::Float64,data::Vector{Float64})::String
            str = @sprintf "%24.12E" t
            for x=1:M
                str = @sprintf "%s%24.12E" str data[x] 
            end
            str = @sprintf "%s\n" str
            return str
        end
        # add to handler, open later
        for V in V_array
            path = joinpath(out_folder,@sprintf "pair_correlations_g2_n%02d_%s.dat" Asize format(calculation_label,V))
            add!(file_handler,path,out_str_pcf_03,handler_name;hname_format=hname_format,replace_list=[V]) 
        end
    end
    # 2.4 state output file
    if c[:save_states] || c[:load_states]
        handler_name = "states"
        hname_format = "_V{:5.3f}"  
        if c[:states_file] === nothing
            if (c[:load_states])
                if c[:ftime_max] === nothing
                     ftime_max= time_range[end]
                else
                     ftime_max= c[:ftime_max]
                end 
                if c[:ftime_min] === nothing
                     ftime_min= time_range[1]
                else
                     ftime_min= c[:ftime_min]
                end
            elseif (c[:save_states])
                ftime_max= time_range[end]
                ftime_min= time_range[1]
            end
            Ψt_output = @sprintf "psioft_%02d_%02d_%+5.3f_%+5.3f_{:+5.3f}_%+5.3f_%6.4f_%06.3f_%06.3f.dat" M N V0 Vp0 Vp Δt ftime_min ftime_max
        else
            Ψt_output=c[:states_file]
        end
        # we handle writing Ψ with a different function and only grab the file using the handler
        data_to_str = x->nothing
        # add to handler, open later
        for V in V_array 
            path = joinpath(out_folder,format(Ψt_output,V))
            add!(file_handler,path,data_to_str,handler_name;hname_format=hname_format,replace_list=[V]) 
        end
    end
    # 2.5 obdm
    if c[:obdm]
        handler_name = "obdm"
        hname_format = "_V{:5.3f}"  
        obdm_name = @sprintf "obdm_%02d_%02d_{:+5.3f}_%+5.3f_tstart%+5.3f_tend%+5.3f_tstep%+5.3f.dat" M N Vp time_range[1] time_range[end] (time_range[2]-time_range[1])
        # we handle writing obdm with a different function and only grab the file using the handler
        data_to_str = x->nothing
        # add to handler, open later
        for V in V_array  
            path = joinpath(out_folder_obdm,format(obdm_name,V))
            add!(file_handler,path,data_to_str,handler_name;hname_format=hname_format,replace_list=[V])  
        end
    end
    return file_handler
end

function write_headers_to_files(
        file_handler::FileOutputHandler,
        V::Float64,
        c::Dict{Symbol,Any},
        M::Int64,
        N::Int64,
        t::Float64,
        V0::Float64, 
        Vp0::Float64,
        Vp::Float64,
        boundary::BdryCond,
        Asize::Int64,
        ℓsize::Int64,
        time_range::Vector{Float64},
        Δt::Float64)

    # 2.1 particleEE
    write_str(file_handler,"particleEE", "# M=$(M), N=$(N), V0=$(V0), Vp0=$(Vp0), Vp=$(Vp), V=$(V), t=$(t), n=$(Asize), tstart=$(c[:time_min]), tstop=$(c[:time_max]), tstep=$(c[:time_step]), $(boundary)\n";replace_list=[V])
    write_str(file_handler,"particleEE", "# start time $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n";replace_list=[V])
    write_str(file_handler,"particleEE",(@sprintf "#%24s#%24s#%24s%24s#%24s#%24s#%24s#%24s#%24s#%24s#%24s#%24s\n" "t" "S₁(n=$(Asize))" "S₂(n=$(Asize))" "S₃(n=$(Asize))" "S₄(n=$(Asize))" "S₅(n=$(Asize))" "S₆(n=$(Asize))" "S₇(n=$(Asize))" "S₈(n=$(Asize))" "S₉(n=$(Asize))" "S₁₀(n=$(Asize))" "S₀₋₅(n=$(Asize))");replace_list=[V])
    # 2.2 spatialEE
    if c[:spatial] 
    write_str(file_handler,"spatialEE", "# M=$(M), N=$(N), V0=$(V0), Vp0=$(Vp0), Vp=$(Vp), V=$(V), t=$(t), l=$(ℓsize), tstart=$(c[:time_min]), tstop=$(c[:time_max]), tstep=$(c[:time_step]), $(boundary)\n";replace_list=[V])
    write_str(file_handler,"spatialEE", "# start time $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n";replace_list=[V])
    write_str(file_handler,"spatialEE",(@sprintf "#%24s#%24s#%24s%24s#%24s#%24s#%24s#%24s#%24s#%24s#%24s#%24s\n" "V" "S₁(n=$(ℓsize))" "S₂(n=$(ℓsize))" "S₃(n=$(ℓsize))" "S₄(n=$(ℓsize))" "S₅(n=$(ℓsize))" "S₆(n=$(ℓsize))" "S₇(n=$(ℓsize))" "S₈(n=$(ℓsize))" "S₉(n=$(ℓsize))" "S₁₀(n=$(ℓsize))" "S₀₋₅(n=$(ℓsize))");replace_list=[V])
    end
    # 2.3 g2
    if c[:g2] 
    write_str(file_handler,"g2", "# M=$(M), N=$(N), V0=$(V0), Vp0=$(Vp0), Vp=$(Vp), V=$(V), t=$(t), n=$(Asize), tstart=$(c[:time_min]), tstop=$(c[:time_max]), tstep=$(c[:time_step]), $(boundary)\n";replace_list=[V])
    write_str(file_handler,"g2", "# start time $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n" ;replace_list=[V])
    write_str(file_handler,"g2",(@sprintf "%24s" "V") ;replace_list=[V] )
    for x=0:M-1
        write_str(file_handler,"g2",(@sprintf "%24d" x) ;replace_list=[V] )
    end 
    write_str(file_handler,"g2","\n";replace_list=[V] ) 
    end

    return nothing
end

function close_file_handler(file_handler::FileOutputHandler,V::Float64)
    for handler_name in ["particleEE","spatialEE","g2"]
        write_str(file_handler,handler_name, "\n\n Calculation finished at  $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))"; replace_list=[V])  
        close(file_handler,handler_name; replace_list=[V]) 
    end
end