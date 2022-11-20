using BSON: @save, @load
using Statistics
using DelimitedFiles
using Printf
using MDToolbox
using StatsBase 
using LinearAlgebra
using Random

include("../src/afm.jl")

nq = 576
qs = readdlm("../data/quaternion/QUATERNION_LIST_$(nq)_Orient")
models = load_model("../data/t1r/cluster.pdb");
nmodel = size(models, 1)

seed = MersenneTwister(123); ############################################
sigma_noise = 3  #################################
nframe = 100
probe_radius = 25 #########################
sharpness = 10
param = AfmizeConfig(sharpness * (pi / 180.0),
    probe_radius, 
    MDToolbox.Point2D(-250, -200), 
    MDToolbox.Point2D(250, 200), 
    MDToolbox.Point2D(6.25, 6.25), 
    MDToolbox.defaultParameters())
T_rot = Matrix{Float64}(I, nq, nq)
@load "../data/t1r/t1r.bson" T pi_i

run(`mkdir -p data/test_case/radius_$(probe_radius)_realization02/sharpness_$(sharpness)/q_$(nq)`)

Threads.@threads for init_q in 1:nq
    imodel_array = rand_generate(seed, nframe, T, pi_i)

    p0_q = zeros(Float64, nq)
    p0_q[init_q] = 1.0
    iq_array = rand_generate(seed, nframe, T_rot, p0_q);
    dxdy_array = randn(seed, Float64, nframe, 2);

    afms = []
    for iframe in 1:nframe
        imodel = imodel_array[iframe]
        model = models[imodel, :]
        iq = iq_array[iframe]
        q = qs[iq, :]
        model = MDToolbox.rotate(model, q)
        model.xyz[:, 1:3:end] .+=  dxdy_array[iframe, 1]
        model.xyz[:, 2:3:end] .+=  dxdy_array[iframe, 2]
        afm = MDToolbox.afmize(model, param)
        h, w = size(afm)
        afm .+= randn(seed, h, w) * sigma_noise
        push!(afms, afm)
    end

    @save "data/test_case/radius_$(probe_radius)_realization02/sharpness_$(sharpness)/q_$(nq)/iq_$(init_q)_noise_$(sigma_noise)_nframe_$(nframe).bson" afms param imodel_array dxdy_array
end