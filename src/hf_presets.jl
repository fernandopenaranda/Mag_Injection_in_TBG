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
using Roots
using ProgressMeter
using Interpolations
using JET
using Profile

const a0 = 2.46 # Å

@with_kw struct ParamsHF # [in eV]
    μ::Float64
    θ::Float64
    ν::Int
    nmax::Int
    pointsk::Int
    v::Float64
    vp::Float64
    vpp::Float64
    M::Float64
    γ::Float64
    λ::Float64
    U1::Float64
    U2::Float64
    n::Float64 # occupation of the flat bands. Deactivated by ρG
    nf::Float64 # occupation of the f electrons for a given n
    sigmaz::Float64
    sigmazlayerz::Float64
    layerz::Float64
    vafek::Float64
    c2ybreakingmass::Float64
    J::Float64
    VP::Bool
    twovalleys::Bool
    twovalleystwospins::Bool
    KIVC::Bool
end


# ks, Gs, Ms are adimensional when required will be divided by const a0 in Å 
# Units: meV, Å. You have to multiply the velocities by a0
paramsHF(θ::Float64, nmax, ν = 1; kpoints = 30) = 
    # ParamsHF(0.0, θ, ν, nmax, kpoints, -4.303*1e3, 1.623*1e3, 3.697, -24.75, 
    #     110#=lambda(θ)=#, 100, 0.0, 0.0, 0.0)
     ParamsHF(0.0, θ, ν, nmax, kpoints, -5*1e3, 1.623*1e3, -0.0332, 3.697, -24.75, 
        110#=lambda(θ)=#, 100, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, false, false, false, false) 


lambda(θ) = a0 * 0.3375/(2pi*θ/(2*360)) #in Å

function paramsHF!(p::ParamsHF, field_name::Symbol, new_value)
    field_names = fieldnames(typeof(p))
    if field_name in field_names
        p.field_name = new_value
    else
        error("Field '$field_name' not found in struct.")
    end
    return p
end

""" returns the indices n1 and n2 that define rectangular mesh of reciprocal vectors M1 and M2"""
function gs_indices(p::ParamsHF; ϵ = 1e-10) 
# gs_indices(nmax::Int64; ϵ = 1e-10)
    nmax = p.nmax
    num_indices = (2 * nmax + 1) * (2 * nmax + 1)
    g1, g2 = bravais_vectors(p)
    indices = Vector{Tuple{Int64, Int64}}(undef, num_indices)
    index = 1
    rad = nmax * norm(g1) + ϵ
    aux = 0.0
    for n1 in -nmax:nmax
        for n2 in -nmax:nmax
            Gpoint = norm(n1*g1 + n2*g2)
            if Gpoint ≤ rad
                indices[index] = (n1, n2)
                index += 1
            else nothing end
        end
    end
    return indices[1:index-1]
end


""" returns the G in Gmesh with indices (i, j) """
function gs_vector(g1, g2, index_pair::Tuple) 
    n1g1, n2g2 = index_pair .* [g1, g2]
    return [n1g1[1]+n2g2[1], n1g1[2]+n2g2[2]]
end

function gs_vector!(g::Vector, g1, g2, index_pair)
    n1g1, n2g2 = index_pair .* [g1, g2]
    g .= [n1g1[1]+n2g2[1], n1g1[2]+n2g2[2]]
end

# Geometry
bravais_vectors(p::ParamsHF) = bravais_vectors(p.θ, 2π*SA[1.0, 1/√3]/a0, 2π*SA[-1.0, 1/√3]/a0)
bravais_vectors(θ::Float64, b1::SVector{2, Float64}, b2::SVector{2, Float64}) = g1(θ, b1), g2(θ, b2)
g1(θ::Float64, b1::SVector{2, Float64}) = rot_bot(θ) * b1 .- rot_top(θ) * b1
g2(θ::Float64, b2::SVector{2, Float64}) = rot_bot(θ) * b2 .- rot_top(θ) * b2
κ(G1, G2) = 1/3 * (2*G1+G2)
m(G1, G2) = 1/2 * (G1+G2)
rot_top(θ) = SA[cos(θ*π/360) -sin(θ*π/360); sin(θ*π/360) cos(θ*π/360)]
rot_bot(θ) = rot_top(-θ)

#### AUXILIARY
ham_matrix_size(p::ParamsHF) = (p.twovalleys == true ? 2 : 1) *(p.twovalleystwospins == true ? 4 : 1) * ham_matrix_size(gs_indices(p))
ham_matrix_size(indices::Vector{Tuple{Int64, Int64}}) = 2 + 4length(indices)    
# 2 are the f orbitals and 4length. Sublattice degof freedom