job_id = parse(Int, ARGS[1])
jobs_num = parse(Int, ARGS[2])
PID = ARGS[3]
U = parse(Int, ARGS[4])
J = parse(Int, ARGS[5])

n = range(0, 4, step = 4/jobs_num)[job_id+1]
include("HF_Optics.jl")

p = ParamsHF(paramsHF(1.05, 1,  1), # custom presets
    μ = 0,
    nf = 1,
    λ = 100,
    M = 5,
    v = -7e3,
    vp = 2e3,
    γ =-30,
    sigmaz = 0,
    sigmazlayerz = 2.5,
    U1 = U,
    U2 = 0,
    J = J,
    VP =true,
    twovalleystwospins= true)

print("filling: ", n)

# tol, conv_iterations, rand_calls, int_evals, points 100 100 1000 100
config = SelfConsistency_config_random(1e-6, 15, 20, 8000, 200)
b = guided_sweep(:β, p, config, [n])
save_sbp(PID, job_id, p, b[1], b[2], b[3], [n])

str = pwd() * "/slurm-" * string(PID) * "." * string(job_id)
str_dest = pwd() * "/Data" * "/" * string(PID) * "/" * string(job_id) * "/"
mv(str * ".out", str_dest * "output.out")
mv(str * ".err", str_dest * "error.err")