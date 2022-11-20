using Distributed
@everywhere using Printf
@everywhere using DelimitedFiles
@everywhere using Random
@everywhere using Statistics
@everywhere using ProgressMeter
@everywhere using Base.Threads
@everywhere using MDToolbox
@everywhere using BSON: @save, @load

@everywhere include("../src/afm.jl")

function analyze()
    nq = 576
    qs = readdlm("../data/quaternion/QUATERNION_LIST_$(nq)_Orient")
    models = load_model("../data/t1r/cluster.pdb");

    test_radius = 25
    pred_radii = [15, 18, 20, 25, 30, 32, 35]
    #pred_radii = [20]
    sigma_noise = 3
    nframe = 100
    sharpness = 10
    for pred_radius in pred_radii
        run(`mkdir -p data/result/test_radius_$(test_radius)/pred_radius_$(pred_radius)/sharpness_$(sharpness)/q_$(nq)`)

        params = [AfmizeConfig(10.0 * (pi / 180),
            r, 
            MDToolbox.Point2D(-250, -200), 
            MDToolbox.Point2D(250, 200), 
            MDToolbox.Point2D(6.25, 6.25), ################################
            MDToolbox.defaultParameters())
            for r in [pred_radius]]

        for iq in 1:nq
            @load "data/test_case/radius_$(test_radius)/sharpness_$(sharpness)/q_$(nq)/iq_$(iq)_noise_$(sigma_noise)_nframe_$(nframe).bson" afms

            # all models with parallel map
            r = getposterior_parallel_use_local_afm(models, afms, qs, sharpness, params)

            @save "data/result/test_radius_$(test_radius)/pred_radius_$(pred_radius)/sharpness_$(sharpness)/q_$(nq)/iq_$(iq)_noise_$(sigma_noise)_nframe_$(nframe).bson" params r
        end
    end
end

@time analyze()