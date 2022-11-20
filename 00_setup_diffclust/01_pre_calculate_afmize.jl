using BSON: @save, @load
using Statistics
using DelimitedFiles
using Printf
using MDToolbox
using StatsBase 

include("../src/afm.jl")

function pre_calculate_afmize()
    nq = 576
    # nq = 4608
    qs = readdlm("../data/quaternion/QUATERNION_LIST_$(nq)_Orient")
    models = load_model("../data/t1r/cluster_2.pdb")

    nmodel = size(models, 1)
    radii = [15, 18, 20, 25, 30, 32, 35]
    # radii = [31, 33, 35, 37, 39, 41, 43, 45]
    # radii = [28, 42]
    # radii = [25]
    # radii = [40, 42, 45]
    sharpness = 10

    for radius in radii
        for imodel in 1:1:nmodel
            command = `mkdir -p ../data/calc_afms/radius_$(radius)/sharpness_$(sharpness)/q_$(nq)/model_diffclust_$(imodel)`
            run(command)
        end
    end

    params = [AfmizeConfig(sharpness * (pi / 180.0),
        r, 
        MDToolbox.Point2D(-250, -200), 
        MDToolbox.Point2D(250, 200), 
        MDToolbox.Point2D(6.25, 6.25), 
        MDToolbox.defaultParameters())
        for r in radii]
    nparam = length(params)

    Threads.@threads for imodel in 1:nmodel
        model = models[imodel, :]
        for iq in 1:nq
            q = qs[iq, :]
            model_rotated = MDToolbox.rotate(model, q)
            ### loop over afmize parameters
            for iparam in 1:nparam
                if isfile("../data/calc_afms/radius_$(radii[iparam])/sharpness_$(sharpness)/q_$(nq)/model_diffclust_$(imodel)/$(iq).bson")
                    continue
                end
                param = params[iparam]
                calculated = MDToolbox.afmize(model_rotated, param)
                @save "../data/calc_afms/radius_$(radii[iparam])/sharpness_$(sharpness)/q_$(nq)/model_diffclust_$(imodel)/$(iq).bson" calculated
            end
        end
    end
end

pre_calculate_afmize()
