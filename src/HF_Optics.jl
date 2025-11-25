using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using Unitful
using Roots
using PhysicalConstants
using PhysicalConstants.CODATA2018
using CairoMakie
using Dates
using GridInterpolations
using OptimizationOptimJL
using ForwardDiff
using IntervalRootFinding
using Interpolations
using JET
using Profile
using Arpack
using Cubature, HCubature 
using Base.Threads
using Distributed
using Dierckx
using DataFrames, CSV, Random
using PolygonOps, Statistics

const s0 = SA[1 0; 0 1]
const σ0 = SA[1 0; 0 1]
const σx = SA[0 1; 1 0]
const σy = SA[0 -1im; 1im 0]
const σz = SA[1 0; 0 -1]

const k_B = (PhysicalConstants.CODATA2018.k_B |> u"eV/mK").val
const ħ = PhysicalConstants.CODATA2018.ħ
const e = PhysicalConstants.CODATA2018.e
const C = ((e^3 / ħ^2) |> u"μA/V^2/s").val
const C_cd = ((e^2/ħ) |> u"μA/V").val
const ħ_ev_s = (ħ |> u"eV*s").val

struct My_optim{T}
    of::Vector{T}  
    e::T
    average_sampling_diff::T
    status::Symbol         
end

struct SelfConsistency_config{T}
    global_int_evals::T
    global_points::T
    global_max_evals::T
    global_penalty::T
    local_int_evals::T
    local_points::T
    local_max_evals::T
    local_penalty::T
    max_count::T
end

struct SelfConsistency_config_random{T}
    tol::Float64
    conv_iterations::T
    rand_calls::T
    int_evals::T
    points::T
end

include("hf_presets.jl")
include("optics_module.jl")
include("hf_model.jl")
include("hf_twovalleys_model.jl")
include("hf_velocities.jl")
include("hf_observables.jl")
include("hf_plots.jl")
include("optics_maps.jl")
include("hf_optics_wrapper.jl")
include("optic_maps_correlated.jl")
include("hf_selfconsistency.jl")
include("sc_penalties.jl")
include("optim_algorithms.jl")
include("model_projectors.jl")
include("rz_operator_ansatz.jl")
include("save_data.jl")
include("alphabeta_competition.jl")
include("sublattice_masses.jl")
include("figures_functions.jl")