using BSON: @save, @load
using Statistics
using DelimitedFiles
using Printf
using MDToolbox
using StatsBase

function quate_dist(q1, q2)
    # return min(sum(q1 .+ q2) .^ 2, sum(q1 .- q2) .^ 2)
    # return acos(max(-1, min(1, 2.0 * sum(q1.*q2).^2 - 1.0)))
    return 1.0 - sum(q1.*q2)^2
end

function load_model(pdb_file)
    models = readpdb(pdb_file);
    for iatom = 1:models.natom
        models.atomname[iatom] = models.resname[iatom]
    end
    MDToolbox.decenter!(models)    
    
    return models
end

function rand_sample(seed, p)
  p_cum = cumsum(p) ./ sum(p)
  r = rand(seed)
  for i = 1:length(p_cum)
    if r <= p_cum[i]
      return i
    end
  end
end

function rand_generate(seed, nframe, T, pi_i)
  states = zeros(typeof(nframe), nframe)

  states[1] = rand_sample(seed, pi_i)
  for iframe = 2:nframe
    states[iframe] = rand_sample(seed, T[states[iframe-1], :])
  end
  return states
end

function my_logprob_eachmodel(afm_array, imodel::Int64, q_array::Matrix{Float64}, sharpness, param_array)
    nafm = length(afm_array)
    nq = size(q_array, 1)
    nparam = length(param_array)
    logprob_array = zeros(eltype(afm_array[1]), nafm, nq, nparam)
    dxdy_array = Array{Tuple{Int,Int}}(undef, nafm, nq, nparam)
    dxdy = Array{Tuple{Int,Int}}(undef, nafm)
    x_center = ceil(Int64, (size(afm_array[1],2)/2)+1.0)
    y_center = ceil(Int64, (size(afm_array[1],1)/2)+1.0)
  
    ### loop over rotations
    for iq in 1:nq
        ### loop over afmize parameters
        for iparam in 1:nparam
            radius = Int(param_array[iparam].probeRadius)
            @load "../data/calc_afms/radius_$(radius)/sharpness_$(sharpness)/q_$(nq)/model_diffclust_$(imodel)/$(iq).bson"  calculated
            for iafm in 1:nafm
                observed = afm_array[iafm]
                logprob = MDToolbox.calcLogProb(observed, calculated)
                logprob_array[iafm, iq, iparam] = logsumexp(logprob)
                dx = argmax(logprob)[2] - x_center
                dy = argmax(logprob)[1] - y_center
                dxdy_array[iafm, iq, iparam] = (dx, dy)
            end
        end
    end

    for iafm in 1:nafm
        ind = argmax(logprob_array[iafm, :, :])
        dxdy[iafm] = dxdy_array[iafm, :, :][ind]
    end

    return (logprob_array=logprob_array, dxdy=dxdy)
end


function getposterior_parallel_use_local_afm(models::TrjArray, afm_array, q_array::Matrix{Float64}, sharpness, param_array)
    nmodel = models.nframe
    nafm = length(afm_array)
    nq = size(q_array, 1)
    nparam = length(param_array)

    p = pmap(x -> my_logprob_eachmodel(afm_array, x, q_array, sharpness, param_array), 1:nmodel)

    logprob_all = []
    logprob_model = []
    logprob_q = []
    logprob_param = []
    dxdy_best = []
    for iafm = 1:nafm
        pp = zeros(Float64, nmodel, nq, nparam)
        for imodel = 1:nmodel
            for iq = 1:nq
                for iparam = 1:nparam
                    pp[imodel, iq, iparam] = p[imodel].logprob_array[iafm, iq, iparam]
                end
            end
        end

        push!(logprob_all, pp)

        t = zeros(Float64, nmodel)
        for imodel = 1:nmodel
            t[imodel] = logsumexp(pp[imodel, :, :][:])
        end
        push!(logprob_model, t)

        imax = argmax(t)
        push!(dxdy_best, p[imax].dxdy[iafm])

        t = zeros(Float64, nq)
        for iq = 1:nq
            t[iq] = logsumexp(pp[:, iq, :][:])
        end
        push!(logprob_q, t)

        t = zeros(Float64, nparam)
        for iparam = 1:nparam
            t[iparam] = logsumexp(pp[:, :, iparam][:])
        end
        push!(logprob_param, t)
    end

    return (all=logprob_all, 
            model=logprob_model, 
            q=logprob_q, 
            param=logprob_param,
            dxdy=dxdy_best)
end

function make_emission_max(results, nmodel, nq)
    emission_max = zeros(nmodel, nq)
    rs = deepcopy(results)
    test_case_num = size(rs, 1)
    
    for icase in 1:test_case_num
        nframe = size(rs[icase].all, 1)
        for iframe = 1:nframe
            rs[icase].all[iframe] .-= maximum(rs[icase].all[iframe])
            rs[icase].all[iframe] = exp.(rs[icase].all[iframe])
            emission = reshape(sum(rs[icase].all[iframe], dims = 3), nmodel, nq)
            emission_max .= max.(emission_max, emission)
        end
    end
    
    return emission_max
end

function get_valid_ids(results, nmodel, nq)
    emission_max = make_emission_max(results, nmodel, nq)
    
    arr = reshape(deepcopy(emission_max), nmodel * nq)
    sort!(arr, rev=true)
    cou = 0
    border = 1
    for ele in arr
        if ele == 1
            continue
        end
        cou += 1
        if cou >= 30 || ele < 1e-300
            break
        end
        border = ele
    end
    
    ids = findall(x -> x >= border, emission_max)
    return ids
end

function make_expanded_transition(T, ids, qs)
    nid = size(ids, 1)
    expanded_T = zeros(nid, nid)
    for i in 1:nid
        for j in 1:nid
            if quate_dist(qs[ids[i][2], :], qs[ids[j][2], :]) > 0.1
                continue
            end
            
            expanded_T[i, j] = T[ids[i][1], ids[j][1]]
        end
        expanded_T[i, :] ./= sum(expanded_T[i, :])
    end
    
    return expanded_T
end

function make_expanded_pi_i(pi_i, ids)
    nid = size(ids, 1)
    return [pi_i[ids[i][1]] for i in 1:nid]
end

function logsumexp_all(x)
    max_x = maximum(x)
    exp_x = exp.(x .- max_x)
    return log(sum(exp_x)) + max_x
end

function logsumexp_rows(x)
    max_x = maximum(x, dims=1)
    exp_x = exp.(x .- max_x)
    return log.(sum(exp_x, dims=1)) .+ max_x
end

function make_expanded_emission(ids, results)
    rs = deepcopy(results)
    nid = size(ids, 1)
    test_case_num = size(rs, 1)
    frame_sum = 0
    for icase in 1:test_case_num
        frame_sum += size(rs[icase].all, 1)
    end
    
    count = 1
    obs_list = []
    #log_emission = ones(nid, frame_sum) .* -10
    log_emission = ones(nid, frame_sum) .* -10000
    
    for icase in 1:test_case_num
        obs = Array{Int64, 1}()
        nframe = size(rs[icase].all, 1)
        for iframe in 1:nframe
            push!(obs, count)
            #rs[icase].all[iframe] .-= maximum(rs[icase].all[iframe])
            #rs[icase].all[iframe] = exp.(rs[icase].all[iframe])
            rs[icase].all[iframe] .= rs[icase].all[iframe] .- logsumexp_all(rs[icase].all[iframe])
            #rs[icase].all[iframe] .= exp.(rs[icase].all[iframe])
            for (i, id) in enumerate(ids)
                tmp = 0.0
                for iparam in 1:size(rs[icase].all[iframe], 3)
                    if typeof(id) == CartesianIndex{2}
                        tmp += exp(rs[icase].all[iframe][id, iparam])
                    else
                        tmp += exp(rs[icase].all[iframe][id..., iparam])
                    end
                end
                if tmp > 1e-100
                    log_emission[i, count] = log(tmp)
                end
            end

            #emission[:, count] ./= sum(emission[:, count])
            count += 1
        end
        
        push!(obs_list, obs)
    end

    #@show minimum(log_emission, dims=1)
    return obs_list, log_emission
end

function devide_ids(state_arr, ids)
    imodel_arr = []
    iq_arr = []
    for state in state_arr
        push!(imodel_arr, ids[state][1])
        push!(iq_arr, ids[state][2])
    end
    
    return imodel_arr, iq_arr
end

function calc_rmsds_of_movie(models, qs, imodel_arr_a, iq_arr_a, imodel_arr_b, iq_arr_b)
    rmsds = []
    nframe = size(imodel_arr_a, 1)
    for iframe in 1:nframe
        a_model = deepcopy(models[imodel_arr_a[iframe], :])
        MDToolbox.rotate!(a_model, qs[iq_arr_a[iframe], :])
        b_model = deepcopy(models[imodel_arr_b[iframe], :])
        MDToolbox.rotate!(b_model, qs[iq_arr_b[iframe], :])
        push!(rmsds, compute_rmsd(a_model, b_model)[1])
    end
    
    return rmsds
end

function get_max_structues(result)
    nframe = size(result.all, 1)
    imodel_arr = []
    iq_arr = []
    dxdy_arr = []
    
    for iframe in 1:nframe
        id = argmax(result.all[iframe])
        push!(imodel_arr, id[1])
        push!(iq_arr, id[2])    
        push!(dxdy_arr, result.dxdy[iframe])
    end
    
    return imodel_arr, iq_arr, dxdy_arr
end

function viterbi(observation, T, pi_i, log_emission)
    nframe = size(observation, 1)
    nstate = size(T, 1)
    P = zeros(eltype(T), nstate, nframe)
    I = zeros(eltype(T), nstate, nframe)
    state_estimated = zeros(eltype(observation), nframe)

    # initialization
    P[:, 1] .= log.(pi_i) .+ log_emission[:, observation[1]]
    I[:, 1] .= zeros(eltype(T), nstate)

    # argmax forward
    Z = zeros(eltype(T), nstate, nstate)
    for t = 2:nframe
        Z .= P[:, t-1] .+ log.(T)
        I[:, t] .= getindex.(argmax(Z, dims=1), 1)[:]
        P[:, t] .= maximum(Z, dims=1)[:] .+ log_emission[:, observation[t]]
    end

    # termination
    P_star = maximum(P[:, nframe])
    state_estimated[nframe] = argmax(P[:, nframe])
    #@show P

    # decoding
    for t = (nframe-1):-1:1
        state_estimated[t] = I[state_estimated[t+1], t+1]
    end

    return state_estimated
end

function load_afm_from_asd(asd_path)
    asd = readasd(asd_path)
    afms = []
    for i = 1:size(asd.frames, 1)
        push!(afms, asd.frames[i].data)
    end
    
    return afms
end

function transitionmatrix(C; TOLERANCE=10^(-8), verbose=false)
  c = Matrix{Float64}(C)
  nstate = size(c, 1)

  c_sym = c + c'
  x = c_sym

  c_rs = sum(c, dims=2)
  x_rs = sum(x, dims=2)

  delta = 10 * TOLERANCE
  logL_old = rand(nstate, nstate)
  logL = rand(nstate, nstate)
  count_iteration = 0

  x = zeros(Float64, nstate, nstate)
  x_new = zeros(Float64, nstate, nstate)

  while delta > TOLERANCE
    count_iteration = count_iteration + 1
    logL_old = deepcopy(logL)
  
    # fixed-point method
    for i = 1:nstate
      for j = 1:nstate
        denom = 0
        if !iszero(x_rs[i]) denom += (c_rs[i] / x_rs[i]) end
        if !iszero(x_rs[j]) denom += (c_rs[j] / x_rs[j]) end
        if !iszero(denom) x_new[i, j] = (c[i, j] + c[j, i]) / denom end
      end
    end
    
    # update
    x_rs .= sum(x_new, dims=2)
    x .= x_new
    for i = 1:nstate
      for j = 1:nstate
        if !iszero(x[i, j]) & !iszero(x_rs[i])
          logL[i, j] = c[i, j] * log(x[i, j] / x_rs[i]);
        end
      end
    end
        
    delta = sum(abs.(logL_old .- logL)) / (nstate ^ 2)
        
    if verbose & (mod(count_iteration, 10) == 0)
      Printf.@printf("%d iteration  LogLikelihood = %8.5e  delta = %8.5e  tolerance = %8.5e\n", count_iteration, sum(logL), delta, TOLERANCE)
    end
  end
  
  pi_i = x_rs ./ sum(x_rs)
  pi_i = pi_i[:]
  t = x ./ x_rs

  return t, pi_i
end

function baumwelch_forward(data_list, T, pi_i, emission)
    ndata = length(data_list)
    nstate = length(T[1, :])
    logL = zeros(Float64, ndata)
    alpha_list = []
    factor_list = []
    for idata = 1:ndata
        data = data_list[idata]
        nframe = length(data)
        alpha  = zeros(Float64, (nframe, nstate))
        factor = zeros(Float64, nframe)
        alpha[1, :] = pi_i.*emission[:, data[1]]
        factor[1] = sum(alpha[1, :])

        alpha[1, :] = alpha[1, :]./factor[1]
        for iframe = 2:nframe
            alpha[iframe, :] = sum(alpha[iframe-1, :] .* T, dims=1)' .* emission[:, data[iframe]]
            factor[iframe] = sum(alpha[iframe, :])
            alpha[iframe, :] = alpha[iframe, :]./factor[iframe]
        end
        logL[idata] = sum(log.(factor))
        push!(alpha_list, alpha)
        push!(factor_list, factor)
    end
    logL, alpha_list, factor_list
end

function baumwelch_backward(data_list, factor_list, T, pi_i, emission)
    ndata = length(data_list)
    nstate = length(T[1, :])
    logL = zeros(Float64, ndata)
    beta_list = []
    for idata = 1:ndata
        data   = data_list[idata]
        factor = factor_list[idata]
        nframe = length(data)
        beta   = zeros(Float64, (nframe, nstate))
        beta[nframe, :] .= 1.0
        for iframe = (nframe-1):-1:1
            beta[iframe, :] = sum((T .* (emission[:, data[iframe+1]] .* beta[iframe+1, :])'), dims=2) ./ factor[iframe+1]
        end
        logL[idata] = sum(log.(factor))
        push!(beta_list, beta)
    end
    logL, beta_list
end

function baumwelch(data_list, T0, pi_i0, log_emission0; TOLERANCE = 10.0^(-8), MAXITERATION=200)
    ## setup
    check_convergence = Inf64
    count_iteration = 1
    #if not isinstance(data_list, list):
    #    data_list = [data_list]
    ndata = length(data_list)
    logL_old = ones(Float64, ndata)
    nobs = length(log_emission0[1, :])
    
    max_emission = maximum(log_emission0, dims=1)
    nstate = length(T0[1, :])
    T = similar(T0)
    log_emission0_scaled = log_emission0 .- maximum(log_emission0, dims=1)
    log_emission0_scaled .= max.(log_emission0_scaled, -10.0)
    log_emission0_scaled .= log_emission0_scaled .- logsumexp_rows(log_emission0_scaled)
    #@show minimum(log_emission0_scaled, dims=1)
    emission0_scaled = exp.(log_emission0_scaled)
    pi_i = similar(pi_i0)
    while (check_convergence > TOLERANCE) && (count_iteration <= MAXITERATION)
        @show count_iteration
        ## E-step
        logL, alpha_list, factor_list = baumwelch_forward(data_list, T0, pi_i0, emission0_scaled)
        #if count_iteration == 1
        #    @show alpha_list
        #end
        #print("1"); println(logL)
        logL2, beta_list = baumwelch_backward(data_list, factor_list, T0, pi_i0, emission0_scaled)
        #print("2"); println(logL2)
        log_alpha_list = []
        for a in alpha_list
            push!(log_alpha_list, log.(a))
        end
        log_beta_list = []
        for b in beta_list
            push!(log_beta_list, log.(b))
        end
        log_T0 = log.(T0)
        #log_emission0 = log.(emission0)

        ## M-step
        # pi
        # pi = np.zeros(nstate, dtype=np.float64)
        # log_gamma_list = []
        # for idata in range(ndata):
        #     log_gamma_list.append(log_alpha_list[idata] + log_beta_list[idata])
        #     pi = pi + np.exp(log_gamma_list[idata][0, :])
        # pi = pi/np.sum(pi)
        # pi_i = pi_i0
        # emission
        # emission = np.zeros((nstate, nobs), dtype=np.float64)
        # for idata in range(ndata):
        #     data = data_list[idata]
        #     for istate in range(nstate):
        #         for iobs in range(nobs):
        #             id = (data == iobs)
        #             if np.any(id):
        #                 emission[istate, iobs] = emission[istate, iobs] + np.sum(np.exp(log_gamma_list[idata][id, istate]))
        # emission[np.isnan(emission)] = 0.0
        # emission = emission / np.sum(emission, axis=1)[:, None]
        #emission = emission0
        # T
        T = zeros(Float64, (nstate, nstate))
        for idata = 1:ndata
          data = data_list[idata]
          nframe = length(data)
          for iframe = 2:nframe
            #log_xi = bsxfun(@plus, log_alpha{idata}(iframe-1, :)', log_beta{idata}(iframe, :));
            log_xi = log_alpha_list[idata][iframe-1, :] .+ log_beta_list[idata][iframe, :]'
            #T = T .+ exp(bsxfun(@plus, log_xi, log_emission0(:, data(iframe))') + log_T0)./factor{idata}(iframe);
            T = T .+ exp.((log_xi .+ log_emission0_scaled[:, data[iframe]]') .+ log_T0) ./ factor_list[idata][iframe]
          end
        end
        #T[np.isnan(T)] = 0.0
        for i in 1:size(T, 1)
            if iszero(sum(T[i, :]))
                T[i, i] = 1
                continue
            end
            T[i, :] = T[i, :] ./ sum(T[i, :])
        end

        ## reversible T
        T, pi_i = transitionmatrix(T)

        ## Check convergence
        count_iteration += 1
        check_convergence = sum(abs.(logL_old .- logL)) / (nstate ^ 2)
        if mod(count_iteration, 100) == 0
            Printf.@printf("%d iteration: delta = %e  tolerance = %e\n" , count_iteration, check_convergence, TOLERANCE)
        end
        logL_old = logL
        pi_i0 = pi_i
        #emission0 = emission
        T0 = T
    end
    T, pi_i, log_emission0
end
