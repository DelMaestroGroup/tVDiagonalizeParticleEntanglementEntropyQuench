# tVDiagonalizeParticleEntanglementEntropyQuench
Fork of the [code](https://github.com/DelMaestroGroup/tVDiagonalizeTimeEvaluationQuench) for ED time evolution in the t-V model after a quantum quench with optimizations 

### Example
- setup output folders
    ```
    mkdir -p ./out
    mkdir -p ./tmp
    ```
- save states for all times
    ```
    julia ./tV_particleEE_quench.jl --out "./out/" --tmp "./tmp/" --spatial --V-start -3.0 --V-end 3.0 --V-num 5 --time-min 0.0 --time-max 40 --time-step 0.1 --save-states --ee 1 10 5
    ```
- load states and compute entropy
    ```
    julia ./tV_particleEE_quench.jl --out "./out/" --tmp "./tmp/" --spatial --V-start -3.0 --V-end 3.0 --V-num 5 --time-min 0.0 --time-max 40 --time-step 0.1 --load-states --ee 1 10 5
    ```

### Required julia packages:
- ArgParse
-  Arpack  
-  MKL
-  ProgressBars 

### Unit tests
```
    julia ./test/runtests.jl
```

### Usage
```
  tV_particleEE_quench.jl [--out FOLDER] [--tmp FOLDER]
                        [--out-states FOLDER] [--out-obdm FOLDER]
                        [--g2] [--save-states] [--load-states]
                        [--states-file FILE] [--obdm] [--spatial]
                        [--skip-hoffdiag-saving]
                        [--skip-hoffdiag-loading] [--no-flush] [--pbc]
                        [--obc] [--V0 V0] [--Vp0 Vp0]
                        [--V-start V_start] [--V-end V_end]
                        [--V-num V_num] [--V-list [V_num...]]
                        [--Vp Vp] [--logspace] [--t t]
                        [--time-min time] [--time-max time]
                        [--time-step time] [--time-num N] [--time-log]
                        [--ftime-min time] [--ftime-max time] --ee ℓ
                        [-h] M N

positional arguments:
        M                     number of sites (type: Int64)
        N                     number of particles (type: Int64)

optional arguments:
        --out FOLDER          path to output folder, default ./
        --tmp FOLDER          folder for hamiltonian storage location
        --out-states FOLDER   folder for states storage location, default is
                                --out
        --out-obdm FOLDER     folder for obdm storage location (if --obdm
                                provided)
        --g2                  output the pair correlation function ⟨ρ_iρ_0⟩
        --save-states         save the time-evolved state to disk
        --load-states         load the time-evolved state from disk
        --states-file FILE    path to I/O file for states
        --obdm                output the spatial dependence of the OBDM
        --spatial             output the spatial entanglement entropy for ℓ
                                = M/2
        --skip-hoffdiag-saving
                                do not save offdiagonal terms of Hamiltonian
                                for V=0
        --skip-hoffdiag-loading
                                do not load offdiagonal terms of Hamiltonian
                                from V=0 case
        --no-flush            do not flush write buffer to output files in
                                after computation for each V
        -h, --help            show this help message and exit

boundary conditions:
        --pbc                 periodic boundary conditions (default)
        --obc                 open boundary conditions

tV parameters:
        --V0 V0               initial V (type: Float64, default: 0.0)
        --Vp0 Vp0             initial Vp (type: Float64, default: 0.0)
        --V-start V_start     start V (type: Float64, default: -2.0)
        --V-end V_end         end V (type: Float64, default: 2.0)
        --V-num V_num         num of V values (type: Int64, default: 100)
        --V-list [V_num...]   multiple V values used as V_array and ignore
                                V-start, V-end, V-step, and --logspace (type:
                                Float64)
        --Vp Vp               final Vp (type: Float64, default: 0.0)
        --logspace            logarithmic spacing of V values around 0
        --t t                 t value (type: Float64, default: 1.0)

time parameters:
        --time-min time       minimum time (type: Float64, default: 0.0)
        --time-max time       maximum time (type: Float64, default: 5.0)
        --time-step time      time step (type: Float64, default: 0.1)
        --time-num N          number of time (type: Int64)
        --time-log            use logarithmic scale for time
        --ftime-min time      file minimum time (type: Float64)
        --ftime-max time      file maximum time (type: Float64)

entanglement entropy:
        --ee ℓ                compute all EEs with partition size ℓ (type:
                                Int64)
```