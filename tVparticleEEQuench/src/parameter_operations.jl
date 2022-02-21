"""For readability of the main script, the parameter checking
steps are moved to this .jl file that is included
into the main script"""

# ------------------------------------------------------------------------------
function parse_commandline()
    s = ArgParseSettings()
    s.autofix_names = true
    @add_arg_table s begin
        "M"
            help = "number of sites"
            arg_type = Int
            required = true
        "N"
            help = "number of particles"
            arg_type = Int
            required = true
        "--out"
            metavar = "FOLDER"
            help = "path to output folder, default ./"
        "--tmp"
            metavar = "FOLDER"
            help = "folder for hamiltonian storage location"
        "--out-states"
            metavar = "FOLDER"
            help = "folder for states storage location, default is --out"
        "--out-obdm" 
            metavar = "FOLDER"
            help = "folder for obdm storage location (if --obdm provided)"
        "--g2"
            help = "output the pair correlation function ⟨ρ_iρ_0⟩"
            action = :store_true
        "--save-states"
            help = "save the time-evolved state to disk"
            action = :store_true
        "--load-states"
            help = "load the time-evolved state from disk"
            action = :store_true
        "--states-file"
            metavar = "FILE"
            help = "path to I/O file for states"
        "--obdm"
            help = "output the spatial dependence of the OBDM"
            action = :store_true
        "--spatial"
            help = "output the spatial entanglement entropy for ℓ = M/2"
            action = :store_true
        "--skip_hoffdiag_saving"
            help = "do not save offdiagonal terms of Hamiltonian for V=0" 
            action = :store_true
        "--skip_hoffdiag_loading"
            help = "do not load offdiagonal terms of Hamiltonian from V=0 case" 
            action = :store_true
        "--no-flush"
            help = "do not flush write buffer to output files in after computation for each V" 
            action = :store_true 

    end
    add_arg_group(s, "boundary conditions")
    @add_arg_table s begin
        "--pbc"
            help = "periodic boundary conditions (default)"
            arg_type = BdryCond
            action = :store_const
            dest_name = "boundary"
            constant = PBC
            default = PBC
        "--obc"
            help = "open boundary conditions"
            arg_type = BdryCond
            action = :store_const
            dest_name = "boundary"
            constant = OBC
            default = PBC
    end
    add_arg_group(s, "tV parameters")
    @add_arg_table s begin 
        "--V0"
            metavar = "V0"
            help = "initial V"
            arg_type = Float64
            default = 0.0
        "--Vp0"
            metavar = "Vp0"
            help = "initial Vp"
            arg_type = Float64
            default = 0.0
        "--V-start"
            metavar = "V_start"
            help = "start V"
            arg_type = Float64
            default = -2.0
        "--V-end"
            metavar = "V_end"
            help = "end V"
            arg_type = Float64
            default = 2.0
        "--V-num"
            metavar = "V_num"
            help = "num of V values"
            arg_type = Int64
            default = 100 
        "--V-list"
            metavar = "V_num"
            nargs = '*' 
            help = "multiple V values used as V_array and ignore V-start, V-end, V-step, and --logspace"
            arg_type = Float64 
            required = false
        "--Vp"
            metavar = "Vp"
            help = "final Vp"
            arg_type = Float64
            default = 0.0
        "--logspace" 
            help = "logarithmic spacing of V values around 0"
            action = :store_true 
        "--t"
            metavar = "t"
            help = "t value"
            arg_type = Float64
            default = 1.0
    end
    add_arg_group(s, "time parameters")
    @add_arg_table s begin
    "--time-min"
        metavar = "time"
        help = "minimum time"
        arg_type = Float64
        default = 0.0
    "--time-max"
        metavar = "time"
        help = "maximum time"
        arg_type = Float64
        default = 5.0
    "--time-step"
        metavar = "time"
        help = "time step"
        arg_type = Float64
        default = 0.1
    "--time-num"
        metavar = "N"
        help = "number of time"
        arg_type = Int
    "--time-log"
        help = "use logarithmic scale for time"
        action = :store_true
    "--ftime-min"
        metavar = "time"
        help = "file minimum time"
        arg_type = Float64
    "--ftime-max"
        metavar = "time"
        help = "file maximum time"
        arg_type = Float64

end
    add_arg_group(s, "entanglement entropy")
    @add_arg_table s begin
        "--ee"
            metavar = "ℓ"
            help = "compute all EEs with partition size ℓ"
            arg_type = Int
            required = true
    end

        return parse_args(s, as_symbols=true)
end


function unpack_parameters(c::Dict{Symbol,Any})
    # Number of sites
    M = c[:M]
    # Number of particles
    N = c[:N]
    # Hopping
    t = c[:t]
    # Boundary conditions
    boundary = c[:boundary]
    # Size of region A
    Asize = c[:ee]
    ℓsize = div(M, 2)
    # Initial V
    V0 = c[:V0]
    Vp0 = c[:Vp0]
    # Initial time
    time_min = c[:time_min]
    # Interaction paramers V, V'  
    if   ~(c[:V_list]===nothing) && length(c[:V_list])> 0
        V_array = c[:V_list] 
    elseif c[:logspace]
        V_array = log_range(c[:V_start],c[:V_end],c[:V_num]) 
    else
        V_array = lin_range(c[:V_start],c[:V_end],c[:V_num])
    end
    Vp = c[:Vp]
    # check parameters
    if M!=2N
        println("Not at half-filling: the number of sites =", M," and the number of particles =",N )
        exit(1)
    end
    if (c[:save_states] && c[:load_states]) || (~c[:save_states] && ~c[:load_states])
        println("exactly one option should be true: --save-states or --load-states. exit()")
        exit(1)
    end 
    if c[:time_step] === nothing
        if c[:time_num] === nothing
            time_range = c[:time_min]:0.5:c[:time_max]
            time_num=length(time_range) 
            Δt = 0.5
        else
            if c[:time_log]
                time_range = log_range(c[:time_min], c[:time_max], c[:time_num])
            else
                time_range = lin_range(c[:time_min], c[:time_max], c[:time_num])
                if length(time_range) > 1
                    Δt = time_range[2]-time_range[1]
                else
                    Δt = time_range[1]
                end
            end
            time_num=c[:time_num]
        end
    else
        if c[:time_num] === nothing
            time_range = c[:time_min]:c[:time_step]:c[:time_max]
            time_num=length(time_range) 
            Δt = c[:time_step]
        else
            println("--time-step and --time-num may not both be supplied")
            exit(1)
        end
    end

    return M,N,t,V0,V_array,Vp0,Vp,boundary,Asize,ℓsize,collect(time_range),time_num,Δt 
end
