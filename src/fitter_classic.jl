
function calculate_fit_for_observables_classic(container::Vector{Observablecontainer},obsname::String)
    values = get_bin_value_for_observable(container,obsname)
    sumw2 = get_sumw2_for_obervable(container,obsname)
    names = values[:,1]
    params = values[:,2]
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
    for i_bin in 1:n_bins
        coeffs = fit_a_bin(g,params_clean,values,i_bin)
        push!(fit_coefficiants, coeffs)
        push!(chi_test_res,do_chi_square(container,obsname,params_clean,i_bin,coeffs,sumw2[:,2+i_bin]))
    end
    fit_coefficiants = permutedims(reshape(hcat(fit_coefficiants...), (length(fit_coefficiants[1]), length(fit_coefficiants))))
    chi_test_res = permutedims(reshape(hcat(chi_test_res...), (length(chi_test_res[1]), length(chi_test_res))))
    return [fit_coefficiants,chi_test_res]
end

function do_ipol_classic(container::Vector{Observablecontainer};parallel=true,prof=false)
    if parallel
        do_ipol_par_classic(container::Vector{Observablecontainer})
    elseif prof
        do_ipol_prof(container::Vector{Observablecontainer})
    else
        do_ipol_nonpar_classic(container::Vector{Observablecontainer})
    end
end

function do_ipol_nonpar_classic(container::Vector{Observablecontainer})
    observables = get_observables_from_container(container)
    binedges = Vector{Vector{Float64}}(undef,length(observables))
    coefficiants = Vector{Vector{Matrix{Float64}}}(undef,length(observables))
    for i in 1:length(observables)
        binedges[i] = get_bin_edges_for_observables(container,observables[i])
        coefficiants[i] = calculate_fit_for_observables_classic(container,observables[i])
        print(string("Calculating Observable Number ",i," out of ",length(observables)," "))
        print(string("Done with ",100*i/length(observables), "% \n"))
    end
    hcat(observables,binedges,coefficiants)
end

function do_ipol_par_classic(container::Vector{Observablecontainer})
    observables = get_observables_from_container(container)
    binedges = Vector{Vector{Float64}}(undef,length(observables))
    coefficiants = Vector{Vector{Matrix{Float64}}}(undef,length(observables))
    result = Folds.collect( do_ipol_step_classic(i,observables,container) for i in 1:length(observables) )
    resultform = Matrix{Any}(undef,length(result),3)
    for i in 1:length(result)
    	resultform[i,1] = result[i][1][i]
    	resultform[i,2] = result[i][2]
    	resultform[i,3] = result[i][3]
    end
    return resultform
end

function do_ipol_step_classic(i::Int,observables::Vector{String},container::Vector{Observablecontainer})
    binedges = get_bin_edges_for_observables(container,observables[i])
    coefficiants = calculate_fit_for_observables_classic(container,observables[i])
    print(string("Calculating Observable Number ",i," out of ",length(observables)," "))
    print(string("Done with ",100*i/length(observables), "% \n"))
    return [observables,binedges,coefficiants]
end
