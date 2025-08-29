
@with_kw struct Histcontainer
    h::Histogram
    sumw2::Vector{Float64}
    name::String
    parname::Vector{String}
    parval::Vector{Float64}
    parval_ID::Int64=NaN
end

@with_kw struct Observablecontainer
    hc::Vector{Histcontainer}
    name::String
    bins::Vector{Float64}
    nbins = length(bins)-1
end

# function Histcontainer(h::Histogram,sumw2::Vector{Float64},name::String,parname::Vector{String},parval::Vector{Float64},parval_ID::Int64)
#     Histcontainer(h,sumw2,name,parname,parval,parval_ID)
# end

function Observablecontainer(hc::Vector{Histcontainer},name::String,bins::Vector{Float64})
    Observablecontainer(hc,name,bins,length(bins)-1)
end


function load_Histcontainer_from_file(yoda_file::String,params_file::String,file_ID::Int=NaN)
    file_par = get_paramas_from_file(params_file)
    par_names = file_par[:,1]
    par_vals = file_par[:,2]
    #println(file_ID)
    #println("yoda_file: ", yoda_file)
    file_hists = get_all_histograms(yoda_file)
    file_hists_sumw2 = get_all_sumw2(yoda_file)
    n_hists  = length(file_hists[:,1])
    container = Vector{Histcontainer}(undef,n_hists)
    for i_hist in 1:n_hists
        @assert file_hists_sumw2[i_hist,1] == file_hists[i_hist,1]
        container[i_hist] = Histcontainer(file_hists[i_hist,2],file_hists_sumw2[i_hist,2],file_hists[i_hist,1],par_names,par_vals,file_ID)
    end
    container
end

function sort_Histcontainer_by_observable(container::Vector{Vector{Histcontainer}})
    number_of_observables =  length(container[1])
    number_of_samples = length(container)
    obs_container = Vector{Observablecontainer}(undef,number_of_observables)
    for i_observable in 1:number_of_observables
        observable_hiscontainer = Vector{Histcontainer}(undef,number_of_samples)
        observable_name = container[1][i_observable].name
        bins = container[1][i_observable].h.edges[1]
        for i_sample in 1:number_of_samples
            obs_index = locate_observable(container[i_sample],observable_name)
            observable_hiscontainer[i_sample] = container[i_sample][obs_index]
        end
        obs_container[i_observable] = Observablecontainer(observable_hiscontainer,observable_name,bins)
    end
    obs_container
end


function load_Histcontainer_from_folder(mc_folder::String; mc_numbers::Vector{String}=Vector{String}([]), mode_csv::Bool=false, parallel::Bool=true)
    if mc_numbers == []
        mc_numbers = cd(readdir,mc_folder)
    end
    container = Histcontainer[]
    container = Vector{Histcontainer}(undef,length(mc_numbers))

    yoda_files = Vector{String}(undef,length(mc_numbers))
    param_files = Vector{String}(undef,length(mc_numbers))
    mc_IDs = Vector{Int64}(undef,length(mc_numbers))
    i_nums = Vector{String}(undef,length(mc_numbers))

    for i in 1:length(mc_numbers)
        i_num = mc_numbers[i]
        i_nums[i] = i_num
        mc_IDs[i] = parse(Int64,i_num)
        files = cd(readdir,string(mc_folder,i_num))
        println(i_num)
        if mode_csv == true
            yoda_files[i] = files[occursin.(Ref(".csv"),files)][1]
        else
            yoda_files[i] = files[occursin.(Ref(".yoda"),files)][1]
        end
        param_files[i] = files[occursin.(Ref("params.dat"),files)][1]
        if parallel == true
            continue
        else
            tempcont = load_Histcontainer_from_file(string(mc_folder,i_num,"/",yoda_file),string(mc_folder,i_num,"/",param_file),mc_ID)
            push!(container,tempcont)
        end
    end
    if parallel == true
        container = Folds.collect( load_Histcontainer_from_file(string(mc_folder,i_nums[i],"/",yoda_files[i]),string(mc_folder,i_nums[i],"/",param_files[i]),mc_IDs[i]) for i in 1:length(mc_numbers) )
    end
    container
end

function get_observables_from_container_full(container::Vector{Observablecontainer})
    [obsc.name for obsc in container]
end

function get_observables_from_container(container::Vector{Observablecontainer})
    [obsc.name for obsc in container if !occursin("RAW",obsc.name)]
end

function give_bin_i(hc::Histcontainer,i::Int64)
    if i > length(hc.h.weights)
        return 999999
    else
        hc.h.weights[i]
    end
end

function get_paramas_from_mc(container::Vector{Observablecontainer})
    nsam = length(container[1].hc)
    npar = length(container[1].hc[1].parname)

    m_all = Matrix{Any}(undef,(nsam,2))
    for i in 1:nsam
        m = Matrix{Any}(undef,(npar,2))
        m[:,1] = container[1].hc[i].parname
        m[:,2] = container[1].hc[i].parval
        m_all[i,1] = container[1].hc[i].parval_ID
        m_all[i,2] = m
    end
    m_all
end

function get_paramas_from_mc(container::Observablecontainer)
    nsam = length(container.hc)
    npar = length(container.hc[1].parname)

    m_all = Matrix{Any}(undef,(nsam,2))
    for i in 1:nsam
        m = Matrix{Any}(undef,(npar,2))
        m[:,1] = container.hc[i].parname
        m[:,2] = container.hc[i].parval
        m_all[i,1] = container.hc[i].parval_ID
        m_all[i,2] = m
    end
    m_all
end

function convert_get_paramas_from_mc(container::Vector{Observablecontainer},digits::Int64=3)
    pars = get_paramas_from_mc(container)
    nrow = length(pars[:,1])
    for i in 1:nrow
        num = pars[i,1]
        pars[i,1] = string("0000",num)[end-digits+1:end]
    end
    pars
end

function convert_get_paramas_from_mc(container::Observablecontainer,digits::Int64=3)
    pars = get_paramas_from_mc(container)
    nrow = length(pars[:,1])
    for i in 1:nrow
        num = pars[i,1]
        pars[i,1] = string("0000",num)[end-digits+1:end]
    end
    pars
end

function make_clean_params(unclean_params)
    n_params = Int(length(unclean_params[1])/2)
    n_values = length(unclean_params)
    params_clean = Array{Float64,1}[]
    for i in unclean_params
        subparams = []
        for j in 1:n_params
            push!(subparams,i[n_params+j])
        end
        push!(params_clean,subparams)
    end
    params_clean = transpose(reshape(vcat(params_clean...),n_params,n_values))
end

function get_clean_params_from_mc(container::Vector{Observablecontainer})
    mc_params = get_paramas_from_mc(container::Vector{Observablecontainer})
    make_clean_params(mc_params[:,2])
end

function locate_observable(container::Vector{Histcontainer},name::String)
    n = length(container)
    index = Vector{Int}(undef,n)
    names = Vector{String}(undef,n)
    for i in 1:length(container)
        if occursin(name,container[i].name)
            return i
        end
    end
end

function locate_observable(container::Vector{Observablecontainer},name::String)
    n = length(container)
    for i in 1:length(container)
        if occursin(name,container[i].name)
            return i
        end
    end
end

function locate_mcID(container::Vector{Histcontainer},name::String)
    n = length(container)
    for i in 1:length(container)
        ID = parse(Int64,name)
        digits = 0
        if ID == 1
                digits = 1
            elseif ID == 0
                digits = 1
            else
                digits = Int(floor(log10(ID)+1))
        end

        if isequal(name[end-digits+1:end],string(container[i].parval_ID))
            return i
        end
    end
end

function locate_mcID(container::Vector{Histcontainer},ID::Int64)
    n = length(container)
    for i in 1:length(container)
        if ID == container[i].parval_ID
            return i
        end
    end
end

function get_bin_edges_for_observables(container::Vector{Observablecontainer},obsname::String)
    container[locate_observable(container,obsname)].bins
end

function get_bin_edges_for_observables(container::Observablecontainer)
    container.bins
end

function get_bin_value_for_observable(container::Vector{Observablecontainer},obsname::String)
    mc_params = convert_get_paramas_from_mc(container)
    obscontainer = container[locate_observable(container,obsname)]
    n_bins = 0
    hist_vs_params = []
    for i in mc_params[:,1]
        hists = obscontainer.hc
        hist = hists[locate_mcID(hists,i)].h
        n_bins = length(hist.weights)
        hist_vs_params = vcat(hist_vs_params,hist.weights)
    end
    hist_vs_params = reshape(hist_vs_params,(n_bins,length(mc_params[:,1])))
    hcat(mc_params,transpose(hist_vs_params))
end

function get_bin_value_for_observable(container::Observablecontainer)
    mc_params = convert_get_paramas_from_mc(container)
    n_bins = length(container.bins)-1
    ids = [container.hc[i].parval_ID for i in 1:length(container.hc)]
    #sp = sortperm(ids)
    #hist_vs_params = reshape(vcat([container.hc[sp][i].h.weights for i in 1:length(ids)]...),(n_bins,length(mc_params[:,1])))
    hist_vs_params = reshape(vcat([container.hc[i].h.weights for i in 1:length(ids)]...),(n_bins,length(mc_params[:,1])))
    hcat(mc_params,transpose(hist_vs_params))
end

function get_sumw2_for_obervable(container::Vector{Observablecontainer},obsname::String)
    mc_params = convert_get_paramas_from_mc(container)
    obscontainer = container[locate_observable(container,obsname)]
    n_params = length(mc_params[:,1])
    n_bins = 0
    sumw2_vs_params = []
    for i in mc_params[:,1]
        hists = obscontainer.hc
        hist = hists[locate_mcID(hists,i)].sumw2
        n_bins = length(hist)
        sumw2_vs_params = vcat(sumw2_vs_params,hist)
    end
    sumw2_vs_params = reshape(sumw2_vs_params,(n_bins,length(mc_params[:,1])))
    hcat(mc_params,transpose(sumw2_vs_params))
end

function get_sumw2_for_obervable(container::Observablecontainer)
    mc_params = convert_get_paramas_from_mc(container)
    n_bins = length(container.bins)-1
    ids = [container.hc[i].parval_ID for i in 1:length(container.hc)]
    sp = sortperm(ids)
    hist_vs_params = reshape(vcat([container.hc[sp][i].sumw2 for i in 1:length(ids)]...),(n_bins,length(mc_params[:,1])))
    hcat(mc_params,transpose(hist_vs_params))
end

function filter_container_by_ranges(container::Vector{Vector{Histcontainer}}, parameter_ranges::Vector{Vector{Float64}})
    parameter_ranges = hcat(parameter_ranges...)
    for j in 1:length(parameter_ranges[:,1])
        print(parameter_ranges[j,:])
        container = [container[i] for i in 1:length(container) if container[i][1].parval[j] > parameter_ranges[j,1]]
        container = [container[i] for i in 1:length(container) if container[i][1].parval[j] < parameter_ranges[j,2]]
        print(length(container))
    end
    return container
end

function filter_container_by_ranges(container::Vector{Vector{Histcontainer}},parameter_ranges::Vector{Uniform{Float64}})
    filter_container_by_ranges(container,[[i.a for i in parameter_ranges],[i.b for i in parameter_ranges]])
end