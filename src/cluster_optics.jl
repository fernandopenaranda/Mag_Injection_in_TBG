job_id = parse(Int, ARGS[1])
folder = string(ARGS[2])
PID = string(ARGS[3])
rs = string(ARGS[4]) # vp or qah

job_id += 1
include("HF_Optics.jl")
include("degen_decision.jl")

evals = 20000
omegas = collect(0:0.15:40) 


p = read_struct_data(folder)
s = provide_folder_get_observables_atlas_reshuffled(folder, rs)
s = discontinuities_fix(s, rs)

self_consistent_shift_map_correlated(folder,job_id-1 ,:REAL, :x, :x, :x, p,  [s.mus[job_id]], [diagm.(s.ofmats)[job_id]], [s.ns[job_id]], omegas,
    rs = "m1", gs = rs, η = 0.5, evals = evals)

self_consistent_shift_map_correlated(folder, job_id-1,  :REAL, :y, :y, :y, p, [s.mus[job_id]], [diagm.(s.ofmats)[job_id]], [s.ns[job_id]], omegas,
    rs = "m1", gs = rs,  η = 0.5, evals = evals)

self_consistent_injection_map_correlated(folder, job_id-1, :REAL, :x, :x, :x, p,  [s.mus[job_id]], [diagm.(s.ofmats)[job_id]], [s.ns[job_id]], omegas,
    rs = "m1",gs = rs,  η = 0.5, evals = evals)

self_consistent_injection_map_correlated(folder, job_id-1, :REAL, :y, :y, :y, p,  [s.mus[job_id]], [diagm.(s.ofmats)[job_id]], [s.ns[job_id]], omegas,
    rs = "m1",gs = rs,  η = 0.5, evals = evals)

job_id -= 1
str = pwd() * "/slurm-" * string(PID) * "." * string(job_id)
str_dest = pwd() * "/Data" * "/" * string(folder) * "/" * string(job_id) * "/"
mv(str * ".out", str_dest * "output_opt.out",force=true)
mv(str * ".err", str_dest * "error_opt.err", force=true)

