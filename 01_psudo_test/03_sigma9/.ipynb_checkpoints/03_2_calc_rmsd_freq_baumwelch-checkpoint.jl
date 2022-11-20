using BSON: @save, @load
using Statistics
using DelimitedFiles
using Printf
using MDToolbox
using StatsBase 
using LinearAlgebra
using Random
using Plots

include("../../src/afm.jl")

function calc_rmsd_freq_baumwelch()
    nq = 576
    qs = readdlm("../../data/quaternion/QUATERNION_LIST_$(nq)_Orient")
    models = load_model("../../data/t1r/cluster.pdb");
    nmodel = size(models, 1)
    seed = MersenneTwister(334)

    test_radius = 25
    pred_radii = [15, 18, 20, 25, 30, 32, 35]
    #pred_radii =          [20, 25, 30, 32, 35]
    sigma_noise = 9
    sharpness = 10
    nframe = 100
    extracted_qs = sample(seed, 1:576, 50, replace=false)

    for pred_radius in pred_radii
        run(`mkdir -p ../data/viterbi/test_radius_$(test_radius)_sigma_$(sigma_noise)/pred_radius_$(pred_radius)/sharpness_$(sharpness)/q_$(nq)`)
        rs = []
        for iq in extracted_qs
            @load "../data/result/test_radius_$(test_radius)_sigma_$(sigma_noise)/pred_radius_$(pred_radius)/sharpness_$(sharpness)/q_$(nq)/iq_$(iq)_noise_$(sigma_noise)_nframe_$(nframe).bson" params r
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
        
        rmsds_arr = zeros(nq, nframe)

        for (itest, iq) in enumerate(extracted_qs)
            @load "../data/test_case/radius_$(test_radius)_sigma_$(sigma_noise)/sharpness_$(sharpness)/q_$(nq)/iq_$(iq)_noise_$(sigma_noise)_nframe_$(nframe).bson" afms param imodel_array dxdy_array
            @load "../data/result/test_radius_$(test_radius)_sigma_$(sigma_noise)/pred_radius_$(pred_radius)/sharpness_$(sharpness)/q_$(nq)/iq_$(iq)_noise_$(sigma_noise)_nframe_$(nframe).bson" params r

            estimated_state = viterbi(obs_list[itest], learned_T, learned_pi_i, log_emission)
            est_imodel_arr, est_iq_arr = devide_ids(estimated_state, ids)
            rmsds = calc_rmsds_of_movie(models, qs, imodel_array, [iq for i in 1:nframe], est_imodel_arr, est_iq_arr)

            rmsds_arr[iq, :] = rmsds
        end

        @save "../data/viterbi/test_radius_$(test_radius)_sigma_$(sigma_noise)/pred_radius_$(pred_radius)/sharpness_$(sharpness)/q_$(nq)/rmsds_baumwelch.bson" rmsds_arr
        @save "../data/viterbi/test_radius_$(test_radius)_sigma_$(sigma_noise)/pred_radius_$(pred_radius)/sharpness_$(sharpness)/q_$(nq)/T_baumwelch.bson" T learned_T learned_pi_i obs_list log_emission ids nid
    end 
end

calc_rmsd_freq_baumwelch()
