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

function calc_rmsd_freq_rand()
    nq = 576
    qs = readdlm("../../data/quaternion/QUATERNION_LIST_$(nq)_Orient")
    models = load_model("../../data/t1r/cluster.pdb");
    nmodel = size(models, 1)

    test_radius = 25
    pred_radii = [15, 18, 20, 25, 30, 32, 35]
    sigma_noise = 6
    sharpness = 10
    nframe = 100

    for pred_radius in pred_radii
        run(`mkdir -p ../data/viterbi/test_radius_$(test_radius)_sigma_$(sigma_noise)/pred_radius_$(pred_radius)/sharpness_$(sharpness)/q_$(nq)`)
        rmsds_arr = zeros(nq, nframe)
        for iq in 1:nq
            @load "../data/test_case/radius_$(test_radius)_sigma_$(sigma_noise)/sharpness_$(sharpness)/q_$(nq)/iq_$(iq)_noise_$(sigma_noise)_nframe_$(nframe).bson" afms param imodel_array dxdy_array
            @load "../data/result/test_radius_$(test_radius)_sigma_$(sigma_noise)/pred_radius_$(pred_radius)/sharpness_$(sharpness)/q_$(nq)/iq_$(iq)_noise_$(sigma_noise)_nframe_$(nframe).bson" params r

            ids = get_valid_ids([r], nmodel, nq)
            nid = size(ids, 1)

            expanded_T = rand(nid, nid)
            for i in 1:nid
                expanded_T[i, :] ./= sum(expanded_T[i, :])
            end
            expanded_pi_i = rand(nid)
            expanded_pi_i ./= sum(expanded_pi_i)
            
            obs_list, log_emission = make_expanded_emission(ids, [r])
            estimated_state = viterbi(obs_list[1], expanded_T, expanded_pi_i, log_emission)
            est_imodel_arr, est_iq_arr = devide_ids(estimated_state, ids)
            rmsds = calc_rmsds_of_movie(models, qs, imodel_array, [iq for i in 1:nframe], est_imodel_arr, est_iq_arr)

            rmsds_arr[iq, :] = rmsds
        end

        @save "../data/viterbi/test_radius_$(test_radius)_sigma_$(sigma_noise)/pred_radius_$(pred_radius)/sharpness_$(sharpness)/q_$(nq)/rmsds_rand.bson" rmsds_arr
    end 
end

calc_rmsd_freq_rand()