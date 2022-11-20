using BSON: @save, @load
using Statistics
using DelimitedFiles
using Printf
using MDToolbox
using StatsBase 
using LinearAlgebra
using Random

pred_radius = 20

arg_id = parse(Int, ARGS[1])
#arg_id = 1

include("../../src/afm.jl")

nq = 576
qs = readdlm("../../data/quaternion/QUATERNION_LIST_$(nq)_Orient")
models = load_model("../../data/t1r/cluster.pdb");
nmodel = size(models, 1)

sigma_noise = 3
nframe = 100
test_radius = 25
sharpness = 10
param = AfmizeConfig(sharpness * (pi / 180.0),
    test_radius, 
    MDToolbox.Point2D(-250, -200), 
    MDToolbox.Point2D(250, 200), 
    MDToolbox.Point2D(6.25, 6.25), 
    MDToolbox.defaultParameters())
T_rot = Matrix{Float64}(I, nq, nq)

@load "../../data/t1r/t1r.bson" T pi_i

true_T = deepcopy(T)
true_pi_i = deepcopy(pi_i);

@load "../data/viterbi/test_radius_$(test_radius)/pred_radius_$(pred_radius)/sharpness_$(sharpness)/q_$(nq)/T_baumwelch.bson" T learned_T learned_pi_i obs_list log_emission ids nid

seed = MersenneTwister(arg_id);
extracted_qs = sample(seed, 1:nq, nmodel, replace=false)

function make_merged_T(ids, T, p, nmodel)
    ret_T = zeros(nmodel, nmodel)
    ret_p = zeros(nmodel)
    for i in 1:length(ids)
        i_state = ids[i][1]
        for j in 1:length(ids)
            j_state = ids[j][1]
            ret_T[i_state, j_state] += T[i, j]
        end

        ret_p[i_state] += p[i]
    end

    for i in 1:nmodel
        if sum(ret_T[i, :]) <= 1e-100
            continue 
        end
        ret_T[i, :] ./= sum(ret_T[i, :])
    end

    return ret_T, ret_p
end

learned_T2, learned_pi2 = make_merged_T(ids, learned_T, learned_pi_i, nmodel)

for (itest, init_q) in enumerate(extracted_qs)
    imodel_array = rand_generate(seed, nframe, learned_T2, learned_pi2)

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

    #@save "data/test_case/radius_$(probe_radius)/sharpness_$(sharpness)/q_$(nq)/iq_$(init_q)_noise_$(sigma_noise)_nframe_$(nframe).bson" afms param imodel_array dxdy_array
    @save "arg_id_$(arg_id)_iq_$(init_q).bson" afms param imodel_array dxdy_array
end

using Distributed
addprocs(8)

@everywhere using Printf
@everywhere using DelimitedFiles
@everywhere using Random
@everywhere using Statistics
@everywhere using ProgressMeter
@everywhere using Base.Threads
@everywhere using MDToolbox
@everywhere using BSON: @save, @load
@everywhere include("../../src/afm.jl")

param = [AfmizeConfig(10.0 * (pi / 180),
    r, 
    MDToolbox.Point2D(-250, -200), 
    MDToolbox.Point2D(250, 200), 
    MDToolbox.Point2D(6.25, 6.25), ################################
    MDToolbox.defaultParameters())
    for r in [pred_radius]]

for (itest, init_q) in enumerate(extracted_qs)
    @load "arg_id_$(arg_id)_iq_$(init_q).bson" afms
    r = getposterior_parallel_use_local_afm(models, afms, qs, sharpness, param)
    @save "arg_id_$(arg_id)_iq_$(init_q)_likelihood.bson" param r
end

rmprocs(2:9)

rs = []
for init_q in extracted_qs
    @load "arg_id_$(arg_id)_iq_$(init_q)_likelihood.bson" param r
    push!(rs, r)
end

ids = get_valid_ids(rs, nmodel, nq)
nid = size(ids, 1)
T = rand(nid, nid)
pi_i = ones(nid) ./ nid
expanded_init_T = make_expanded_transition(T, ids, qs)
expanded_init_pi_i = make_expanded_pi_i(pi_i, ids)
obs_list, log_emission = make_expanded_emission(ids, rs)

learned_T, learned_pi_i, _ = baumwelch(obs_list, deepcopy(expanded_init_T), expanded_init_pi_i, log_emission)
@show learned_pi_i

@save "arg_id_$(arg_id)_baumwelch.bson" true_T true_pi_i T learned_T learned_pi_i obs_list log_emission ids nid

