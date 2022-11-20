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

function calc_rmsd_freq_origin()
    nq = 576
    qs = readdlm("../../data/quaternion/QUATERNION_LIST_$(nq)_Orient")
    models = load_model("../../data/t1r/cluster.pdb");
    nmodel = size(models, 1)

    test_radius = 25
    pred_radii = [15, 18, 20, 25, 30, 32, 35]
    sigma_noise = 3
    sharpness = 10
    nframe = 100

    for pred_radius in pred_radii
        rmsds_arr = zeros(nq, nframe)
        for iq in 1:nq
            @load "../data/test_case/radius_$(test_radius)_realization03/sharpness_$(sharpness)/q_$(nq)/iq_$(iq)_noise_$(sigma_noise)_nframe_$(nframe).bson" afms param imodel_array dxdy_array
            @load "../data/result/test_radius_$(test_radius)_realization03/pred_radius_$(pred_radius)/sharpness_$(sharpness)/q_$(nq)/iq_$(iq)_noise_$(sigma_noise)_nframe_$(nframe).bson" params r

            est_imodel_arr, est_iq_arr, est_dxdy_arr = get_max_structues(r)
            rmsds = calc_rmsds_of_movie(models, qs, imodel_array, [iq for i in 1:nframe], est_imodel_arr, est_iq_arr)
            rmsds_arr[iq, :] = rmsds
        end

        @save "../data/viterbi/test_radius_$(test_radius)_realization03/pred_radius_$(pred_radius)/sharpness_$(sharpness)/q_$(nq)/rmsds_origin.bson" rmsds_arr
    end
end

calc_rmsd_freq_origin()