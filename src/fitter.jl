function g_quad(n::Vector{T} where T<:Real)
	nparams = length(n)
	quadvec = ones(Int(quad_n_params(nparams)))
	ind = 1
	for i in 1:nparams
		quadvec[i+1] = n[i]
		for j in i+1:nparams
			quadvec[nparams+1+ind] = n[i]*n[j]
			ind = ind+1
		end
		quadvec[end-nparams+i] = n[i]*n[i]
	end
	return quadvec
end

function g_cubic(n::Vector{T} where T<:Real)
	nparams = length(n)
	#cubicvec = ones(Int(cubic_n_params(nparams)))
    T = eltype(n)
    #cubicvec = StaticArrays.SizedArray{Tuple{Int(cubic_n_params(nparams)),}}( ones(T,Int(cubic_n_params(nparams))) )
    cubicvec = ones(T,Int(cubic_n_params(nparams)))
    #cubicvec[1:quad_n_params(length(n))] = g_quad(n)
	indquad = quad_n_params(length(n))
    indij=1
    indi2j=1
    indijk=1
    indstartijk=indquad+nparams*(nparams-1)
	for i in 1:nparams
        cubicvec[i+1] = n[i]
		for j in 1:nparams
            if j > i 
                cubicvec[nparams+1+indij] = n[i]*n[j]
                indij += 1
                for k in j+1:nparams
                    cubicvec[indstartijk+indijk] = n[i]*n[j]*n[k]
                    indijk += 1
                end
            end
            if j != i
                cubicvec[indquad+indi2j] = n[i]^2*n[j]
                indi2j += 1
            end
		end
        cubicvec[end-nparams+i] = n[i]*n[i]*n[i]
        cubicvec[indquad-nparams+i] = n[i]*n[i]
	end
	return cubicvec
end

#= function g_cubic_test(n::Vector{<:String})
	nparams = length(n)
	cubicvec = Vector{String}(undef,Int(cubic_n_params(nparams)))
    #cubicvec[1:quad_n_params(length(n))] = g_quad(n)
	indquad = quad_n_params(length(n))
    indij=1
    indi2j=1
    indijk=1
    indstartijk=indquad+nparams*(nparams-1)
    cubicvec[1] = "1"
	for i in 1:nparams
        cubicvec[i+1] = n[i]
		for j in 1:nparams
            if j > i 
                cubicvec[nparams+1+indij] = string(n[i]," * ",n[j])
                indij += 1
                for k in j+1:nparams
                    cubicvec[indstartijk+indijk] = string(n[i]," * ",n[j]," * ",n[k])
                    indijk += 1
                end
            end
            if j != i
                cubicvec[indquad+indi2j] = string(n[i],"^2"," * ",n[j])
                indi2j += 1
            end
		end
        cubicvec[end-nparams+i] = string(n[i],"^3")
        cubicvec[indquad-nparams+i] = string(n[i],"^2")
	end
	return cubicvec
end =#

function g_quad(n::AbstractArray{<:Real,2})
	return mapslices(g_quad,n,dims=2)
end

function g_quad(n::Vector{Vector{T}} where T<:Real)
	return broadcast(g_quad,n)
end

function g_cubic(n::AbstractArray{<:Real,2})
	return mapslices(g_cubic,n,dims=2)
end

function g_cubic(n::Vector{Vector{T}} where T<:Real)
	return broadcast(g_cubic,n)
end

function g(n,p)
    m = g_quad(n)
    m*p
end

function g_err(n,covar)
    m = g_quad(n)
    m'*covar*m
end

function g_cubic(n,p)
    m = g_cubic(n)
    m*p
end

function g_cubic(n::Vector{T} where T<:Real,p::Vector{T} where T<:Real)
    m = g_cubic(n)
    m'*p
end

function g_err_cubic(n,covar)
    m = g_cubic(n)
    m'*covar*m
end

function poly_n_params(n,P)
    N = 1
    for i in 1:n
        fact = 1
        for j in 1:i
            fact = fact*j
        end
        S=1
        for k in 0:i-1
            S = S*(P+k)
        end
        N += S/fact
    end
    return Int(N)
end

function quad_n_params(n)
    N = Int(1+3/2*n+n^2/2)
    return N
end

function cubic_n_params(n)
    N = Int(1+11/6*n+n^2+1/6*n^3)
    return N
end

function clean_binvals_paramaters(values,parameters,ibin)
    binval = values[:,2+ibin]
    if sum(broadcast(isnan,binval))+sum(binval.==Inf) > 0
        print(string("Attention ",string(sum(broadcast(isnan,binval))+sum(binval.==Inf))," NaN/Inf Values found and deleted."))
        nanloc = findall(i->isnan(i),binval)
        infloc = findall(i->abs(i)==Inf,binval)
        okloc = filter(i->i ∉ vcat(nanloc,infloc),collect(1:length(binval)))

        binval = binval[okloc]
        parameters = parameters[okloc,:]
    end
    (binval,parameters)
end

function fit_a_bin(g,parameters,values,ibin,p0=1;cubic=false)
    if p0 == 1
        p0 = [1.0 for i in 1:quad_n_params(length(parameters[1,:]))]
        if cubic
            p0 = [1.0 for i in 1:cubic_n_params(length(parameters[1,:]))]
        end
    end
    binval,parameters = clean_binvals_paramaters(values,parameters,ibin)
    if cubic
        fit_res = curve_fit(g_cubic,parameters,binval,p0,autodiff=:forwarddiff)
    else
        fit_res = curve_fit(g,parameters,binval,p0,autodiff=:forwarddiff)
    end
    return fit_res
end

function chisq(x,y)
    chi_square = 0
    for i in 1:length(x)
        chi_square += ((x[i]-y[i])^2/y[i])
    end
    return chi_square
end

function chisqerr(x,y,sy)
    chi_square = 0
    for i in 1:length(x)
        chi_square += ((x[i]-y[i])^2/sy[i])
    end
    return chi_square
end

# function do_chi_square(container::Vector{Observablecontainer},obsname::String,params_clean,i_bin,fit_params,err,α=0.05)
#     fit_values = g(params_clean,fit_params)
#     mc_data = get_bin_value_for_observable(container,obsname)[:,2+i_bin]
#     chi_square = chisq(fit_values,mc_data)
#     chi_square_err = chisqerr(fit_values,mc_data,err)
#     dof = length(mc_data)-length(fit_params)
#     chi_hypo=0
#     if dof > 0
#         chi_hypo = cquantile(Chisq(dof),1-α)
#     end
#     return [chi_square_err/dof,chi_square_err,dof,chi_hypo,chi_square,(length(mc_data)-1)]
# end

function do_chi_square(container::Observablecontainer,params_clean,i_bin,fit_params,err,α=0.05;cubic=false)
    if cubic
        fit_values = g_cubic(params_clean,fit_params)
    else
        fit_values = g(params_clean,fit_params)
    end
    mc_data = get_bin_value_for_observable(container)[:,2+i_bin]
    chi_square = chisq(fit_values,mc_data)
    chi_square_err = chisqerr(fit_values,mc_data,err)
    dof = length(mc_data)-length(fit_params)
    chi_hypo=0
    if dof > 0
        chi_hypo = cquantile(Chisq(dof),1-α)
    end
    return [chi_square_err/dof,chi_square_err,dof,chi_hypo,chi_square,(length(mc_data)-1)]
end

# function calculate_fit_for_observables(container::Vector{Observablecontainer},obsname::String; par_mask=[true for i in 1:length(container[1].hc[1].parname)])
#     values = get_bin_value_for_observable(container,obsname)
#     sumw2 = get_sumw2_for_obervable(container,obsname)
#     binedges = get_bin_edges_for_observables(container,obsname)
#     binlen = binedges[2:end] .- binedges[1:end-1]
#     allparams = values[:,2]
#     params = [allparams[i][par_mask,:] for i in 1:length(allparams)]
#     n_params = Int(length(params[1])/2)
#     n_values = length(params)
#     n_bins = length(values[1,:])-2
#     params_clean = Array{Float64,1}[]
#     for i in params
#         subparams = []
#         for j in 1:n_params
#             push!(subparams,i[n_params+j])
#         end
#         push!(params_clean,subparams)
#     end
#     params_clean = transpose(reshape(vcat(params_clean...),n_params,n_values))
#     fit_coefficiants = []
#     chi_test_res = []
#     for i_bin in 1:n_bins
#         coeffs = fit_a_bin(g,params_clean,values,i_bin)
#         push!(fit_coefficiants, coeffs)
#         push!(chi_test_res,do_chi_square(container,obsname,params_clean,i_bin,coeffs,sumw2[:,2+i_bin]/binlen[i_bin]))
#     end
#     fit_coefficiants = permutedims(reshape(hcat(fit_coefficiants...), (length(fit_coefficiants[1]), length(fit_coefficiants))))
#     chi_test_res = permutedims(reshape(hcat(chi_test_res...), (length(chi_test_res[1]), length(chi_test_res))))
#     return [fit_coefficiants,chi_test_res]
# end

function calculate_pull(g,params,coeffs,covar,yval,yerr)
    fitval = g(params,coeffs)
    fiterr = g_err(param,covar)
    (fitval-yval)/sqrt(yerr^2-fiterr^2)
end

function calculate_fit_for_observables(container::Observablecontainer; par_mask=[true for i in 1:length(container.hc[1].parname)],cubic=false)
    values = get_bin_value_for_observable(container)
    sumw2 = get_sumw2_for_obervable(container)
    binedges = get_bin_edges_for_observables(container)
    binlen = binedges[2:end] .- binedges[1:end-1]
    allparams = values[:,2]
    params = [allparams[i][par_mask,:] for i in 1:length(allparams)]
    n_params = Int(length(params[1])/2)
    n_values = length(params)
    n_bins = length(values[1,:])-2
    params_clean = Array{Float64,1}[]
    for i in params
        subparams = []
        for j in 1:n_params
            push!(subparams,i[n_params+j])
        end
        push!(params_clean,subparams)
    end
    params_clean = transpose(reshape(vcat(params_clean...),n_params,n_values))
    fit_coefficiants = []
    chi_test_res = []
    covvars = []
    for i_bin in 1:n_bins
        fit_res = fit_a_bin(g,params_clean,values,i_bin,cubic=cubic)
        covar = estimate_covar(fit_res)
        push!(fit_coefficiants, fit_res.param)
        push!(chi_test_res,do_chi_square(container,params_clean,i_bin,fit_res.param,sumw2[:,2+i_bin]/binlen[i_bin],cubic=cubic))
        push!(covvars, estimate_covar(fit_res))
    end
    fit_coefficiants = permutedims(reshape(hcat(fit_coefficiants...), (length(fit_coefficiants[1]), length(fit_coefficiants))))
    chi_test_res = permutedims(reshape(hcat(chi_test_res...), (length(chi_test_res[1]), length(chi_test_res))))
    return [fit_coefficiants,chi_test_res,covvars]
end

function do_ipol_prof(container::Vector{Observablecontainer};do_chi_sq=true,LIST="/ceph/groups/e4/users/slacagnina/overH70222/longlist_raw.txt")
    observables = get_observables_from_container(container)
    binedges = []
    coefficiants = []

    params = get_clean_params_from_mc(container)
    parsets = g_cubic(params)
    pinv_pars = pinv(reshape(hcat(parsets...),length(parsets[1,:]),length(parsets[:,1])))
    pinv_pars = pinv(parsets)

    goodobs = []
    for i in 1:length(observables)
        if LIST != ""
            if getweight(LIST,string(observables[i]),weight=0) == 0
                continue
            end
        end
        push!(goodobs,observables[i])
        values = get_bin_value_for_observable(container[i])
        sumw2 = get_sumw2_for_obervable(container[i])
        n_bins = length(values[1,:])-2

        coeff_binres = []
        chisq_binres = []
        for i_bin in 1:n_bins
            binvals, = clean_binvals_paramaters(values,params,i_bin)
            coeffs = pinv_pars*binvals
            chi_sq = 1
            if do_chi_sq
                chi_sq = do_chi_square(container[i],params,i_bin,coeffs,sumw2[:,2+i_bin],cubic=true)
            end 
            push!(coeff_binres,coeffs)
            push!(chisq_binres,chi_sq)
        end
        coeff_binres=permutedims(reshape(hcat(coeff_binres...),length(coeff_binres[1]),length(coeff_binres)))
        chisq_binres=permutedims(reshape(hcat(chisq_binres...),length(chisq_binres[1]),length(chisq_binres)))

        push!(binedges,get_bin_edges_for_observables(container[i]))
        push!(coefficiants,[coeff_binres,chisq_binres])

        print(string("Calculating Observable Number ",i," out of ",length(observables)," "))
        print(string("Done with ",100*i/length(observables), "% \n"))
    end
    hcat(goodobs,binedges,coefficiants)
end

function do_ipol(container::Vector{Observablecontainer};par_mask=[true for i in 1:length(container[1].hc[1].parname)],parallel=true,prof=false,cubic=false)
    if parallel
        do_ipol_par(container::Vector{Observablecontainer},par_mask=par_mask,cubic=cubic)
    elseif prof
        do_ipol_prof(container::Vector{Observablecontainer})
    else
        do_ipol_nonpar(container::Vector{Observablecontainer},par_mask=par_mask,cubic=cubic)
    end
end

function do_ipol_nonpar(container::Vector{Observablecontainer};par_mask=[true for i in 1:length(container[1].hc[1].parname)],cubic=false)
    observables = get_observables_from_container(container)
    binedges = Vector{Vector{Float64}}(undef,length(observables))
    coefficiants = Vector{Vector{Matrix{Float64}}}(undef,length(observables))
    coefficiants = Vector{Any}(undef,length(observables))
    for i in 1:length(observables)
        binedges[i] = get_bin_edges_for_observables(container[i])
        coefficiants[i] = calculate_fit_for_observables(container[i],par_mask=par_mask,cubic=cubic)
        print(string("Calculating Observable Number ",i," out of ",length(observables)," "))
        print(string("Done with ",100*i/length(observables), "% \n"))
    end
    hcat(observables,binedges,coefficiants)
end

function do_ipol_par(container::Vector{Observablecontainer};par_mask=[true for i in 1:length(container[1].hc[1].parname)],cubic=false)
    observables = get_observables_from_container(container)
    binedges = Vector{Vector{Float64}}(undef,length(observables))
    coefficiants = Vector{Vector{Matrix{Float64}}}(undef,length(observables))
    result = Folds.collect( do_ipol_step(i,observables,container,par_mask=par_mask,cubic=cubic) for i in 1:length(observables) )
    resultform = Matrix{Any}(undef,length(result),3)
    for i in 1:length(result)
    	resultform[i,1] = result[i][1][i]
    	resultform[i,2] = result[i][2]
    	resultform[i,3] = result[i][3]
    end
    return resultform
end

function do_ipol_step(i::Int,observables::Vector{String},container::Vector{Observablecontainer};par_mask=[true for i in 1:length(container[1].hc[1].parname)],cubic=false)
    println(string("Calculating Observable Number ",i," out of ",length(observables)," "))
    binedges = get_bin_edges_for_observables(container[i])
    coefficiants = calculate_fit_for_observables(container[i],par_mask=par_mask,cubic=cubic)
    println(string("Done with ",100*i/length(observables), "%"))
    return [observables,binedges,coefficiants]
end

function format_ipol_par(ipol)
    ipolform = Matrix{Any}(undef,length(ipol),3)
    for i in 1:length(ipol)
    	ipolform[i,1] = ipol[i][1][i]
    	ipolform[i,2] = ipol[i][2]
    	ipolform[i,3] = ipol[i][3]
    end
    return ipolform
end

function save_my_ipol(ipol,output="./myipol")
    #save(string(output,".jld"),"ipol",ipol) old JLD1
    JLD2.@save string(output,".jld2") ipol
end

function read_my_ipol(input="./myipol")
    #load(string(input,".jld"))["ipol"]
    JLD2.@load string(input,".jld2") ipol
end
