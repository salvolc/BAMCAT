#################################
### Plotting the tuned samples###
#################################
function plot_overview_1D_2D(
    samples::DensitySampleVector,
    name::String,
    parameter_names::Vector{String},
    parameter_symbols::Vector{Expr}
    )

    thispath = pwd()
    if(!(name in readdir()))
    	mkpath(name)
    end
    cd(name)
    if(!("2D" in readdir()))
    	mkpath("2D")
    end
    cd(thispath)

    Plots.plot(samples)

    Plots.savefig(string("./",name,"/results.pdf"))
    Plots.savefig(string("./",name,"/results.png"))

    for i in 1:length(parameter_names)
    	Plots.plot(samples,parameter_symbols[i],xlabel=parameter_names[i],ylabel=string("P(",parameter_names[i],")"))
    	Plots.savefig(string("./",name,"/",i,"_",parameter_names[i],"results.pdf"))
    	Plots.savefig(string("./",name,"/",i,"_",parameter_names[i],"results.png"))

    	for j in i:length(parameter_names)
    		Plots.plot(samples,(parameter_symbols[i],parameter_symbols[j]),xlabel=parameter_names[i],ylabel=parameter_names[j])
    		Plots.savefig(string("./",name,"/2D/",i,"_",j,"_",parameter_names[i],"_",parameter_names[j],"results.pdf"))
    		Plots.savefig(string("./",name,"/2D/",i,"_",j,"_",parameter_names[i],"_",parameter_names[j],"results.png"))
    	end
    end
end

function smival_for_samples(samples,parameter_names,p,bins,atol)
    syms =  keys(samples[1].v)
	smival = Vector{Any}(undef,length(parameter_names))
	for i in 1:length(parameter_names)
		try get_smallest_interval_edges(samples,syms[i],p,bins=bins,atol=atol)
			smival[i] = get_smallest_interval_edges(samples,syms[i],p,bins=bins,atol=atol)
		catch
			smival[i] = (lower = 0.0, upper = 0.0)
		end
	end
    return smival
end

function print_results(
    samples::DensitySampleVector,
    name::String,
    parameter_names::Vector{String};
    digits=3,
    p=0.68,
    bins=250,
    atol = 0.01,
    )

    if(!(name in readdir()))
    	mkpath(name)
    end

    f = open(string(name,"/results.txt"),"w")
    write(f,string("Mode: ",mode(samples),"\n"))
    write(f,string("Mean: ",mean(samples),"\n"))
    write(f,string("Var: ",var(samples),"\n"))
    write(f,"\n\n")

    header = [" ",parameter_names...]
    results_table=[["Mode",round.(collect(mode(samples)[1]),digits=digits)...],["Mean",round.(collect(mean(samples)[1]),digits=digits)...],["Var",round.(collect(var(samples)[1]),digits=digits)...],["Std",round.(sqrt.(abs.(collect(var(samples)[1]))),digits=digits)...]]
    results_table = permutedims(hcat(results_table...))
    pretty_table(f,results_table,header)
    write(f,"\n\n")
    header = [" ","Mode","Mean","Var","Std"]
    results_table=[parameter_names,[round.(collect(mode(samples)[1]),digits=digits)...],[round.(collect(mean(samples)[1]),digits=digits)...],[round.(collect(var(samples)[1]),digits=digits)...],[round.(sqrt.(abs.(collect(var(samples)[1]))),digits=digits)...]]
    results_table = (hcat(results_table...))
    pretty_table(f,results_table,header)

    #smival = [get_smallest_interval_edges(samples,i,p,bins=bins,atol=atol) for i in 1:length(parameter_names)]
    smival = smival_for_samples(samples,parameter_names,p,bins,atol)

    smival_for = [[[round(i.lower[j],digits=digits),round(i.upper[j],digits=digits)] for j in 1:length(i.lower)] for i in smival]

    header = [" ","Mode","Mean","Var","Std","smallest $(Int(p*100))% interval"]
    results_table=[parameter_names,[round.(collect(mode(samples)[1]),digits=digits)...],[round.(collect(mean(samples)[1]),digits=digits)...],[round.(collect(var(samples)[1]),digits=digits)...],[round.(sqrt.(abs.(collect(var(samples)[1]))),digits=digits)...],smival_for]
    results_table = (hcat(results_table...))
    pretty_table(f,results_table,header)

    header = [" ","Marginal Mode","Mean","smallest $(Int(p*100))% interval"]
    results_table=[parameter_names,[round.(collect(bat_marginalmode(samples).result[1]),digits=digits)...],[round.(collect(mean(samples)[1]),digits=digits)...],smival_for]
    results_table = (hcat(results_table...))
    pretty_table(f,results_table,header)
    
    close(f)
end

##########Nice Plots##########
function plot_overview_1D_2D_pyplot(
    samples::DensitySampleVector,
    name::String,
    parameter_names::Vector{String},
    parameter_names_y::Vector{String},
    parameter_symbols
    )

    thispath = pwd()
    if(!(name in readdir()))
    	mkpath(name)
    end
    cd(name)
    if(!("2D" in readdir()))
    	mkpath("2D")
    end
    cd(thispath)

    pyplot()
    Plots.plot(samples,vsel=parameter_symbols,vsel_label=parameter_names,vsel_label_y=parameter_names_y,globalmode=true)
    #Plots.plot(samples,vsel=collect(1:length(parameter_names)),vsel_label=parameter_names,vsel_label_y=parameter_names_y,globalmode=true)
    #Plots.plot!(size=(2500,1500),xtickfontsize=13,ytickfontsize=13,xguidefontsize=13,yguidefontsize=13)
    Plots.plot!(size=(3000,2000),xtickfontsize=14,ytickfontsize=14,xguidefontsize=14,yguidefontsize=14)

    Plots.savefig(string(name,"/results.pdf"))
    Plots.savefig(string(name,"/results.png"))

    for i in 1:length(parameter_names)
    	Plots.plot(samples,parameter_symbols[i],xlabel=parameter_names[i],ylabel=string("P(",parameter_names[i],")"))
    	Plots.savefig(string(name,"/",i,"_",parameter_names[i],"results.pdf"))
    	Plots.savefig(string(name,"/",i,"_",parameter_names[i],"results.png"))

    	for j in i:length(parameter_names)
    		Plots.plot(samples,(parameter_symbols[i],parameter_symbols[j]),xlabel=parameter_names[i],ylabel=parameter_names[j])
    		Plots.savefig(string(name,"/2D/",i,"_",j,"_",parameter_names[i],"_",parameter_names[j],"results.pdf"))
    		Plots.savefig(string(name,"/2D/",i,"_",j,"_",parameter_names[i],"_",parameter_names[j],"results.png"))
    	end
    end
end

##########Compare Plots##########
function plot_comparison_pyplot(
    samples_A::DensitySampleVector,
    samples_B::DensitySampleVector,
    name::String,
    parameter_names::Vector{String},
    parameter_names_y::Vector{String},
    parameter_symbols::Vector{Expr},
    label_A="H7",
    label_B="TheP8I"
    )

    if(!(name in readdir()))
    	mkpath(name)
    end
    cd(name)
    if(!("2D" in readdir()))
    	mkpath("2D")
    end
    cd("../")

    for i in 1:length(parameter_names)
    	plot(samples_A,i,seriestype=:stephist,linecolor=:blue,marginalmode=false,globalmode=true,label=label_A)
    	plot!(samples_B,i,seriestype=:stephist,linecolor=:red,marginalmode=false,globalmode=true,label=label_B)
    	plot!(xlabel=parameter_names[i],ylabel=parameter_names_y[i])
    	Plots.savefig(string("./",names,"/",i,"_",parameter_names[i],"results_comp.pdf"))
    	Plots.savefig(string("./",names,"/",i,"_",parameter_names[i],"results_comp.png"))
    end
end

function plot_contour_2D(
        ref_sample,
        all_samples,
        tpl::Tuple{Symbol,Symbol},
        labels,#::Vector{String},
        all_colors,
        RESPath::String,
        parameter_names::Vector{String},
        parameter_symbols::Vector{Symbol};
        interval=[0.9],
        alpha=0.1,
        smoothing=0.5,
    )

    @assert(length(all_colors[1]) == length(interval))
    @assert(length(parameter_names) == length(parameter_symbols))

    if(!(RESPath in readdir()))
    	mkpath(RESPath)
    end

    modeval = mode(ref_sample)
    indx = [findall(x->x==tpl[1],fieldnames(modeval[])), findall(x->x==tpl[2],fieldnames(modeval[]))]
    parameter_names_tpl = [parameter_names[findall(x->x==tpl[1],parameter_symbols)[1]],parameter_names[findall(x->x==tpl[2],parameter_symbols)[1]]]

    modeval = [modeval[1][indx[1][1]],modeval[1][indx[2][1]]]
    
    plt = plot(ref_sample, tpl, xlabel=parameter_names_tpl[1], ylabel=parameter_names_tpl[2],st=:smallest_intervals_contourf,smoothing=smoothing,color=[:red1,:yellow1,:green1],globalmode=true,marginalmode=false)
    #plt = plot()
    if smoothing == 0
        plt = plot(ref_sample, tpl, xlabel=parameter_names_tpl[1], ylabel=parameter_names_tpl[2])
    end
    for i in 1:length(all_samples)
        samples = all_samples[i]
        plot!([modeval[1]], [modeval[2]], st=:shape, lw=0, color=all_colors[i], label=labels[i])
        plot!(samples, tpl, intervals=interval, colors=all_colors[i],   lw=1.2, st=:smallest_intervals_contourf, alpha=alpha, bins=100, smoothing=1, marginalmode=false)#, xlims=(-1,1), ylims=(-1,1))
    end
    #plot!()
    plot!(size=(600,400),xtickfontsize=13,ytickfontsize=13,xguidefontsize=13,yguidefontsize=13,xlabel=parameter_names_tpl[1], ylabel=parameter_names_tpl[2],legend=:topleft,framestyle=:box)
    savefig(string(RESPath,"contour_$(parameter_names_tpl[1])_$(parameter_names_tpl[2]).png"))
    savefig(string(RESPath,"contour_$(parameter_names_tpl[1])_$(parameter_names_tpl[2]).pdf"))
end

####################################################
####This Plots all fits obs wise min max ranges ####
####################################################
function plot_fit_ranges(ipol, sorted; folder::String="fit_per_obs_filled",cubic=false)
    names = ipol[:,1]
    params = get_clean_params_from_mc(sorted)
    for iobs in 1:length(names)
        obs = names[iobs]
        if occursin("RAW",obs)
            continue
        end
        binvalues = get_bin_value_for_observable(sorted,obs)
        binedges = get_bin_edges_for_observables(sorted,obs)
        binerrs = get_sumw2_for_obervable(sorted,obs)
        nparams = length(params[1,:])
        nbins = length(binedges)-1
        nparcombos = length(params[:,1])

        parval = params[:,:]
        coeffs = ipol[:,3][iobs][1]
        fitvals = [g(parval,coe)[1] for coe in eachrow(coeffs)]
        fit_minvals = fitvals
        fit_maxvals = fitvals
        orgval = Vector{Float64}(binvalues[1,3:end])
        mc_minvals = orgval
        mc_maxvals = orgval


        for iparcombo in 1:nparcombos
            orgval = Vector{Float64}(binvalues[iparcombo,3:end])
            if sum(Vector{Float64}(binerrs[iparcombo,3:end]) .< 0) > 0
                continue
            end
            orgvalerr = sqrt.(Vector{Float64}(binerrs[iparcombo,3:end]))
            parval = params[:,:]
            coeffs = ipol[:,3][iobs][1]
            if cubic
                fitvals = [g_cubic(parvals,coe)[iparcombo] for coe in eachrow(coeffs)]
            else
                fitvals = [g(parval,coe)[iparcombo] for coe in eachrow(coeffs)]
            end
            fit_minvals = [fitvals[i] < fit_minvals[i] ? fitvals[i] : fit_minvals[i] for i in 1:nbins]
            fit_maxvals = [fitvals[i] > fit_maxvals[i] ? fitvals[i] : fit_maxvals[i] for i in 1:nbins]

            mc_minvals = [fitvals[i] < mc_minvals[i] ? orgval[i] : mc_minvals[i] for i in 1:nbins]
            mc_maxvals = [fitvals[i] > mc_maxvals[i] ? orgval[i] : mc_maxvals[i] for i in 1:nbins]

        end

        for iparcombo in 1:nparcombos
            orgval = Vector{Float64}(binvalues[iparcombo,3:end])
            if sum(Vector{Float64}(binerrs[iparcombo,3:end]) .< 0) > 0
                continue
            end
            orgvalerr = sqrt.(Vector{Float64}(binerrs[iparcombo,3:end]))
            parval = params[:,:]
            coeffs = ipol[:,3][iobs][1]
            if cubic
                fitvals = [g_cubic(parvals,coe)[iparcombo] for coe in eachrow(coeffs)]
            else
                fitvals = [g(parval,coe)[iparcombo] for coe in eachrow(coeffs)]
            end            
            fit_midval = fit_minvals.+((fit_maxvals .- fit_minvals)/2)
            fithist = fit(Histogram,binedges[1:end-1],FrequencyWeights(fit_midval),binedges)
            fithistmin = fit(Histogram,binedges[1:end-1],FrequencyWeights(fit_minvals),binedges)
            orghist = fit(Histogram,binedges[1:end-1],Weights(orgval),binedges)
            binmid = binedges[1:end-1]+((binedges[2:end]-binedges[1:end-1])/2)
            plot(orghist,seriestype=:step,label="MC",color="red")
            plot!(fithistmin,fill_between=mc_maxvals,linecolor=:blue,color=:blue,alpha=0.5,label="Fit")
            plot!(fithist,seriestype=:step,color=:grey,label="")
            xlabel!("$(obs)")
            ylabel!("Entries (normalized)")
            title!("$(obs), Par: $(iparcombo)/$(nparcombos)")
            if(!(folder in readdir()))
                mkpath(folder)
            end
            cd(folder)
            if(!(obs in readdir()))
                mkpath(obs[2:end])
            end
            cd("../")
            savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$(iparcombo)-$(nparcombos).pdf"))
            savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$(iparcombo)-$(nparcombos).png"))
        end
    end
end

function plot_fit_ranges_data(ipol,sorted,rivet_data_path; folder::String="fit_per_obs_filled_data",par_mask=[true for i in 1:length(sorted[1].hc[1].parname)],LIST="/ceph/groups/e4/users/slacagnina/overH70222/longlist_raw.txt",cubic=false)
    names = ipol[:,1]
    params = get_clean_params_from_mc(sorted)
    params = params[:,par_mask]
    prepath = pwd()
    mkpath(folder)
    for iobs in 1:length(names)
        obs = names[iobs]
        if occursin("RAW",obs)
            continue
        end
        if LIST != ""
            if getweight(LIST,string(obs),weight=0) == 0
                continue
            end
        end
        binvalues = get_bin_value_for_observable(sorted,obs)
        binedges = get_bin_edges_for_observables(sorted,obs)
        binerrs = get_sumw2_for_obervable(sorted,obs)
        nbins = length(binedges)-1
        nparcombos = length(params[:,1])
        
        if !(string(split(ipol[:,1][iobs],"/")[2], ".yoda") in readdir(rivet_data_path))
            print("Skipped $obs, it is not in the Rivet Data Path")
            continue
        end

        data_histos = get_all_ref_histograms_to_names(
        string(rivet_data_path, split(ipol[:,1][iobs],"/")[2], ".yoda"),
        obs
        )
        data_sumw2 = get_all_ref_sumw2_to_names(
        string(rivet_data_path, split(ipol[:,1][iobs],"/")[2], ".yoda"),
        obs
        )

        parval = params[:,:]
        coeffs = ipol[:,3][iobs][1]
        fitvals_per_comb=[]
        if cubic
            fitvals_per_comb = reshape(vcat([g_cubic(parval,coe) for coe in eachrow(coeffs)]...),(nparcombos,nbins))
        else
            fitvals_per_comb = reshape(vcat([g(parval,coe) for coe in eachrow(coeffs)]...),(nparcombos,nbins))
        end
        mc_minvals = [minimum(i) for i in eachcol(binvalues[:,3:end])]
        mc_maxvals = [maximum(i) for i in eachcol(binvalues[:,3:end])]


        fit_minvals=[minimum(i) for i in eachcol(fitvals_per_comb)]
        fit_maxvals=[maximum(i) for i in eachcol(fitvals_per_comb)]

        # fitvals = [g(parval,coe)[1] for coe in eachrow(coeffs)]
        # fit_minvals = fitvals
        # fit_maxvals = fitvals
        # orgval = Vector{Float64}(binvalues[1,3:end])
        # mc_minvals = orgval
        # mc_maxvals = orgval


        # for iparcombo in 1:nparcombos
        #     orgval = Vector{Float64}(binvalues[iparcombo,3:end])
        #     if sum(Vector{Float64}(binerrs[iparcombo,3:end]) .< 0) > 0
        #         continue
        #     end
        #     orgvalerr = sqrt.(Vector{Float64}(binerrs[iparcombo,3:end]))
        #     parval = params[:,:]
        #     coeffs = ipol[:,3][iobs][1]
            # if cubic
            #     fitvals = [g_cubic(parvals,coe)[iparcombo] for coe in eachrow(coeffs)]
            # else
            #     fitvals = [g(parval,coe)[iparcombo] for coe in eachrow(coeffs)]
            # end
        #     fit_minvals = [fitvals[i] < fit_minvals[i] ? fitvals[i] : fit_minvals[i] for i in 1:nbins]
        #     fit_maxvals = [fitvals[i] > fit_maxvals[i] ? fitvals[i] : fit_maxvals[i] for i in 1:nbins]

        #     mc_minvals = [fitvals[i] < mc_minvals[i] ? orgval[i] : mc_minvals[i] for i in 1:nbins]
        #     mc_maxvals = [fitvals[i] > mc_maxvals[i] ? orgval[i] : mc_maxvals[i] for i in 1:nbins]

        # end
        
        datavals = data_histos[2].weights
        dataerr = data_sumw2[2]
        datahist = Histogram(binedges,datavals)

        binmid = binedges[1:end-1]+((binedges[2:end]-binedges[1:end-1])/2)

        fithistmin = Histogram(binedges,fit_minvals)
        mchistmin = Histogram(binedges,mc_minvals)

        
        cd(folder)
        if (!(obs in readdir()))
            mkpath(split(obs,"/")[2])
            mkpath(string(split(obs,"/")[2],"/CSV/"))
        end
        cd(prepath)

        df = DataFrame(xmin=collect(datahist.edges[1])[1:end-1],xmax=collect(datahist.edges[1])[2:end],
        yval_data = datahist.weights,
        yerr_data = dataerr,
        yval_mcmax = mc_minvals,
        yval_mcmin = mc_maxvals,
        yval_fitmax = fit_maxvals,
        yval_fitmin = fit_minvals
        )
        CSV.write(string(folder,"/",split(obs,"/")[2],"/CSV/",split(obs,"/")[3],".csv"),df)

        plot(datahist,seriestype=:step,label="Data",color="black")
        plot!(binmid,datavals,yerr=dataerr,linewidth=0)
        plot!(fithistmin,fill_between=fit_maxvals,linecolor=:blue,color=:blue,alpha=0.3,label="Fit")
        plot!(mchistmin,fill_between=mc_maxvals,linecolor=:red,color=:red,alpha=0.3,label="MC")
        xlabel!("$(obs)")
        ylabel!("Entries")
        savefig(string(folder,"/",split(obs,"/")[2],"/",split(obs,"/")[3],".pdf"))
        savefig(string(folder,"/",split(obs,"/")[2],"/",split(obs,"/")[3],".png"))

        #yerr = sqrt.(dataerr)
        #for i in 1:length(yerr)
        #    while (datavals[i] - yerr[i]) < 0
        #        yerr[i] = yerr[i]/10
        #    end
        #end

        #plot(datahist,seriestype=:step,label="Data",color="black",yaxis=:log)
        #plot!(binmid,datavals,yerr=sqrt.(dataerr)./10,linewidth=0)
        #plot!(binmid,datavals,yerr=yerr,linewidth=0)
        #plot!(fithistmin,fill_between=fit_maxvals,linecolor=:blue,color=:blue,alpha=0.3,label="Fit")
        #plot!(mchistmin,fill_between=mc_maxvals,linecolor=:red,color=:red,alpha=0.3,label="MC")
        #xlabel!("$(obs)")
        #ylabel!("Entries")
        #savefig(string(folder,"/",split(obs,"/")[2],"/","ylog_",split(obs,"/")[3],".pdf"))
        #savefig(string(folder,"/",split(obs,"/")[2],"/","ylog_",split(obs,"/")[3],".png"))
    end
end

function plot_fit_ranges_data_ylog(
    ipol,
    sorted,
    rivet_data_path; 
    folder::String="fit_per_obs_filled_data", 
    obsnum::Int64=-1,
    par_mask=[true for i in 1:length(sorted[1].hc[1].parname)],
    LIST="/ceph/groups/e4/users/slacagnina/overH70222/longlist_raw.txt",
    cubic=false,
    PDF=false,
    )
    names = ipol[:,1]
    params = get_clean_params_from_mc(sorted)
    params = params[:,par_mask]
    prepath = pwd()
    mkpath(folder)
    for iobs in 1:length(names)
        if obsnum != -1
            if iobs != obsnum
                continue
            end
        end
        obs = names[iobs]
        if occursin("RAW",obs)
            continue
        end
        if LIST != ""
            if getweight(LIST,string(obs),weight=0) == 0
                continue
            end
        end
        binvalues = get_bin_value_for_observable(sorted,obs)
        binedges = get_bin_edges_for_observables(sorted,obs)
        binerrs = get_sumw2_for_obervable(sorted,obs)
        nbins = length(binedges)-1
        nparcombos = length(params[:,1])
        
        if !(string(split(ipol[:,1][iobs],"/")[2], ".yoda") in readdir(rivet_data_path))
            print("Skipped $obs, it is not in the Rivet Data Path")
            continue
        end

        data_histos = get_all_ref_histograms_to_names(
        string(rivet_data_path, split(ipol[:,1][iobs],"/")[2], ".yoda"),
        obs
        )
        data_sumw2 = get_all_ref_sumw2_to_names(
        string(rivet_data_path, split(ipol[:,1][iobs],"/")[2], ".yoda"),
        obs
        )

        parval = params[:,:]
        coeffs = ipol[:,3][iobs][1]
        fitvals_per_comb=[]
        if cubic
            fitvals_per_comb = reshape(vcat([g_cubic(parval,coe) for coe in eachrow(coeffs)]...),(nparcombos,nbins))
        else
            fitvals_per_comb = reshape(vcat([g(parval,coe) for coe in eachrow(coeffs)]...),(nparcombos,nbins))
        end

        mc_minvals = [minimum(i) for i in eachcol(binvalues[:,3:end])]
        mc_maxvals = [maximum(i) for i in eachcol(binvalues[:,3:end])]

        fit_minvals=[minimum(i) for i in eachcol(fitvals_per_comb)]
        fit_maxvals=[maximum(i) for i in eachcol(fitvals_per_comb)]

        
        datavals = data_histos[2].weights
        binl = binedges[2:end] .- binedges[1:end-1]
        dataerr = data_sumw2[2]
        yerr = dataerr
        for i in 1:length(yerr)
            while (datavals[i] - yerr[i]) < 0
                println("Errors had to be reduced foe $obs in bin number $i")
                yerr[i] = yerr[i]/10
                if datavals[i] < 0 
                    datavals[i] = 0
                    yerr[i] = 0
                    println("Errors had to be set to 0 for $obs in bin number $i")
                    break 
                end
            end
        end
        
        binmid = binedges[1:end-1]+((binedges[2:end]-binedges[1:end-1])/2)
        mask = ones(length(datavals))
        for i in 1:length(datavals)
            if fit_minvals[i] <= 0
                if length(fit_minvals[fit_minvals .> 0]) > 0
                    fit_minvals[i] = minimum(fit_minvals[fit_minvals .> 0])*0.9
                else
                    fit_minvals[i] = 0.01 .* fit_maxvals[i]
                end
            end
            if datavals[i] <= 0
                datavals[i] = minimum(datavals[datavals .> 0])*0.9
            end
            if mc_minvals[i] <= 0
                if length(mc_minvals[mc_minvals .> 0]) > 0
                    mc_minvals[i] = minimum(mc_minvals[mc_minvals .> 0])*0.9
                else
                    mc_minvals[i] = 0.01 .* fit_maxvals[i]
                end
            end
        end
        
        fithistmin = Histogram(binedges,fit_minvals)
        mchistmin = Histogram(binedges,mc_minvals)
        datahist = Histogram(binedges,datavals)
        
        cd(folder)
        if (!(obs in readdir()))
            mkpath(split(obs,"/")[2])
            mkpath(string(split(obs,"/")[2],"/CSV/"))
        end
        cd(prepath)
        
        df = DataFrame(xmin=collect(datahist.edges[1])[1:end-1],xmax=collect(datahist.edges[1])[2:end],
        yval_data = datahist.weights,
        yerr_data = dataerr,
        yval_mcmax = mc_minvals,
        yval_mcmin = mc_maxvals,
        yval_fitmax = fit_maxvals,
        yval_fitmin = fit_minvals
        )
        CSV.write(string(folder,"/",split(obs,"/")[2],"/CSV/ylog_",split(obs,"/")[3],".csv"),df)

        obstitle = obs
        if LIST != ""
            obstitle = getobsname(LIST,obs)
        end

        plt = scatter(binmid,datavals,yerr=yerr,markershape=:circle,ms=2.5,lw=1,label="",color="black",yaxis=:log,yminorticks=4,xminorticks=2,minorgrid=true,borderstyle=:box,marker_legend_size=1)#linecolor=invisible())
        scatter!(binmid,datavals,label="Data",ms=1,msw=0,legendmarkersize=2,color=:black)#,thickness_scaling=1.4)
        plot!(fithistmin,fill_between=fit_maxvals,linecolor=:blue,color=:blue,alpha=0.3,label="Fit")
        plot!(mchistmin,fill_between=mc_maxvals,linecolor=:red,color=:red,alpha=0.3,label="MC")
        ylabel!(L"\frac{1}{\sigma}\frac{d\sigma}{dS}")
        #ylabel!(L"1/\sigma \cdot d \sigma/dS")
        xlabel!("$(obstitle)")
        copy_ticks_log(plt,plt[1],xminorticks=2,yminorticks=4)
        plot!(size=(450,300),top_margin=2Plots.mm)
        if PDF
            savefig(string(folder,"/",split(obs,"/")[2],"/","ylog_",split(obs,"/")[3],".pdf"))
        end
        savefig(string(folder,"/",split(obs,"/")[2],"/","ylog_",split(obs,"/")[3],".png"))

        
        rfitmax = fit_maxvals ./ datavals
        rmcmax = mc_maxvals ./ datavals
        rfitmin = fit_minvals ./ datavals
        rmcmin = mc_minvals ./ datavals
        rfithistmin = Histogram(binedges,rfitmin)
        rmchistmin = Histogram(binedges,rmcmin)

        plt = scatter(binmid,datavals,yerr=yerr,markershape=:circle,ms=3,lw=1,label="Data",color="black",yaxis=:log,yminorticks=4,xminorticks=2,minorgrid=true,xformatter=_->"",borderstyle=:box)
        plot!(fithistmin,fill_between=fit_maxvals,linecolor=:blue,color=:blue,alpha=0.3,label="Fit")
        plot!(mchistmin,fill_between=mc_maxvals,linecolor=:red,color=:red,alpha=0.3,label="MC")
        ylabel!("Entries")
        plot!(size=(510,320),top_margin=2Plots.mm,bottom_margin=-10Plots.mm)
        
        pltr = plot(binmid,ones(length(binmid)),link=:all,lw=0,yminorticks=4,xminorticks=2,minorgrid=true,borderstyle=:box,ylim=(minimum(vcat(rfitmin...,rmcmin...)),maximum(vcat(rfitmax...,rmcmax...))),color=:black,label="Data")
        hline!([1],color=:black,label="")
        plot!(rfithistmin,fill_between=fit_maxvals./ datavals,linecolor=:blue,color=:blue,alpha=0.3,label="")#,label="Fit")
        plot!(rmchistmin,fill_between=mc_maxvals./ datavals,linecolor=:red,color=:red,alpha=0.3,label="")#,label="MC")
        xlabel!("$(obstitle)")
        ylabel!("Ratio")
        plot!(top_margin=-1.9Plots.mm)
        
        l = @layout [a{0.7h};b{0.3h}]
        plot(plt,pltr,layout=l,size=(510,600),link=:x)

        if PDF
            savefig(string(folder,"/",split(obs,"/")[2],"/","ratio_ylog_",split(obs,"/")[3],".pdf"))
        end
        savefig(string(folder,"/",split(obs,"/")[2],"/","ratio_ylog_",split(obs,"/")[3],".png"))
    end
end

function plot_fit_ranges_compact(ipol, sorted; folder::String="fit_per_obs_compact",par_mask=[true for i in 1:length(sorted[1].hc[1].parname)],cubic=false,LIST="",PDF=false)
    names = ipol[:,1]
    params = get_clean_params_from_mc(sorted)
    params = params[:,par_mask]
    prepath = pwd()
    mkpath(folder)
    for iobs in 1:length(names)
        obs = names[iobs]
        if occursin("RAW",obs)
            continue
        end
        if LIST != ""
            if getweight(LIST,string(obs),weight=0) == 0
                continue
            end
        end
        binvalues = get_bin_value_for_observable(sorted,obs)
        binedges = get_bin_edges_for_observables(sorted,obs)
        binerrs = get_sumw2_for_obervable(sorted,obs)
        nparams = length(params[1,:])
        nbins = length(binedges)-1
        nparcombos = length(params[:,1])

        parval = params[:,:]
        coeffs = ipol[:,3][iobs][1]
        fitvals = [g(parval,coe)[1] for coe in eachrow(coeffs)]
        fit_minvals = fitvals
        fit_maxvals = fitvals
        orgval = Vector{Float64}(binvalues[1,3:end])
        mc_minvals = orgval
        mc_maxvals = orgval


        for iparcombo in 1:nparcombos
            orgval = Vector{Float64}(binvalues[iparcombo,3:end])
            if sum(Vector{Float64}(binerrs[iparcombo,3:end]) .< 0) > 0
                continue
            end
            orgvalerr = sqrt.(Vector{Float64}(binerrs[iparcombo,3:end]))
            parval = params[:,:]
            coeffs = ipol[:,3][iobs][1]
            if cubic
                fitvals = [g_cubic(parvals,coe)[iparcombo] for coe in eachrow(coeffs)]
            else
                fitvals = [g(parval,coe)[iparcombo] for coe in eachrow(coeffs)]
            end
            fit_minvals = [fitvals[i] < fit_minvals[i] ? fitvals[i] : fit_minvals[i] for i in 1:nbins]
            fit_maxvals = [fitvals[i] > fit_maxvals[i] ? fitvals[i] : fit_maxvals[i] for i in 1:nbins]

            mc_minvals = [fitvals[i] < mc_minvals[i] ? orgval[i] : mc_minvals[i] for i in 1:nbins]
            mc_maxvals = [fitvals[i] > mc_maxvals[i] ? orgval[i] : mc_maxvals[i] for i in 1:nbins]

        end

        binmid = binedges[1:end-1]+((binedges[2:end]-binedges[1:end-1])/2)
        mc_mid = mc_minvals .+ (0.5 .* (mc_maxvals .- mc_minvals))

        fithistmin = Histogram(binedges,fit_minvals)
        mchistmin = Histogram(binedges,mc_minvals)
        mchistmid = Histogram(binedges,mc_mid)


        cd(folder)
        if (!(obs in readdir()))
            mkpath(split(obs,"/")[2])
            mkpath(string(split(obs,"/")[2],"/CSV/"))
        end
        cd(prepath)

        df = DataFrame(xmin=collect(fithistmin.edges[1])[1:end-1],xmax=collect(fithistmin.edges[1])[2:end],
        yval_mcmid = mchistmin.weights,
        yval_mcmax = mc_minvals,
        yval_mcmin = mc_maxvals,
        yval_fitmax = fit_maxvals,
        yval_fitmin = fit_minvals
        )
        CSV.write(string(folder,"/",split(obs,"/")[2],"/CSV/",split(obs,"/")[3],".csv"),df)

        plot(mchistmid,seriestype=:step,label="Data",color="black")
        plot!(fithistmin,fill_between=fit_maxvals,linecolor=:blue,color=:blue,alpha=0.3,label="Fit")
        plot!(mchistmin,fill_between=mc_maxvals,linecolor=:red,color=:red,alpha=0.3,label="MC")

        xlabel!("$(obs)")
        ylabel!("Entries")
        if PDF
            savefig(string(folder,"/",split(obs,"/")[2],"/",split(obs,"/")[3],".pdf"))
        end
        savefig(string(folder,"/",split(obs,"/")[2],"/",split(obs,"/")[3],".png"))

        #plot(mchistmid,seriestype=:step,label="Data",color="black",yaxis=:log)
        #plot!(fithistmin,fill_between=fit_maxvals,linecolor=:blue,color=:blue,alpha=0.3,label="Fit")
        #plot!(mchistmin,fill_between=mc_maxvals,linecolor=:red,color=:red,alpha=0.3,label="MC")
        #xlabel!("$(obs)")
        #ylabel!("Entries")
        #savefig(string(folder,"/",split(obs,"/")[2],"/","ylog_",split(obs,"/")[3],".pdf"))
        #savefig(string(folder,"/",split(obs,"/")[2],"/","ylog_",split(obs,"/")[3],".png"))
        
    end
end

####################################################
####Plot chi2 to data for 2 mc samples plots #######
####################################################

function plot_tuned_chisq(ipol;
    NOMINAL_PATH::String="",
    TUNED_PATH::String="",
    RIVET_DATA_PATH="/ceph/groups/e4/users/slacagnina/Programs/new/Herwig7/share/Rivet/",
    folder="./chisq_tune/tune/",
    MCPATH = "/ceph/groups/e4/users/slacagnina/overH70222/grid_inputs_H7_v42/",
    CONSTR = "/ceph/groups/e4/users/slacagnina/overH70222/bad_data.txt",
    LIST = "/ceph/groups/e4/users/slacagnina/overH70222/longlist_raw.txt",
    NAME_NOMINAL="Nominal",
    NAME_TUNED="Tuned",
    prop_err=[],
    combine_err=true,
    plot=true,
    )

    nom_hists = get_all_histograms(NOMINAL_PATH)
    nom_sumw2 = get_all_sumw2(NOMINAL_PATH)
    nom_histos_names = nom_hists[:,1]

    tune_hists = get_all_histograms(TUNED_PATH)
    tune_sumw2 = get_all_sumw2(TUNED_PATH)
    tune_histos_names = tune_hists[:,1]

    ipol_histos_names = ipol[:,1]

    analyses_input_file = string(MCPATH, lpad(1, 5, "0"), "/params.in")
    io = open(analyses_input_file, "r")
    file_string = read(io, String)
    file_lbl = split(file_string, "\n")
    analyses_input = [
        split(file_lbl[i])[4] for
        i in 1:length(file_lbl) if occursin("Rivet:Analyses", file_lbl[i])
    ]

    data_histos = Vector{Matrix{Any}}(undef, length(analyses_input))
    data_sumw2 = Vector{Matrix{Any}}(undef, length(analyses_input))
    for i in 1:length(analyses_input)
        data_histos[i] = get_all_ref_histograms_to_names(
            string(RIVET_DATA_PATH, analyses_input[i], ".yoda"),
            ipol_histos_names
        )
        data_sumw2[i] = get_all_ref_sumw2_to_names(
            string(RIVET_DATA_PATH, analyses_input[i], ".yoda"),
            ipol_histos_names
        )
    end
    data_histos = vcat(data_histos...)
    data_sumw2 = vcat(data_sumw2...)
    data_histos_names = data_histos[:, 1]
    n_hist = length(data_histos_names)
    n_bins = [length(data_histos[i_hist, 2].weights) for i_hist = 1:n_hist]
    bad_data=split(read(open(CONSTR,"r"),String),"\n")
    chisq_nom=[]
    chisq_tune=[]
    chisq_prob_nom=[]
    chisq_prob_tune=[]
    all_names = []

    for i_hist in 1:n_hist
        n_bin = n_bins[i_hist]
        nom_index = findall(i -> i == data_histos_names[i_hist][5:end], nom_histos_names)
        tune_index = findall(i -> i == data_histos_names[i_hist][5:end], tune_histos_names)
        ipol_index = findall(i -> i == data_histos_names[i_hist][5:end], ipol[:,1])
        skip=false
        if nom_index == []
            println("Skip $(tune_histos_names[tune_index[1],1]) (i=$i_hist), not present in Nom Values")
            println("Skip $(ipol[ipol_index,1]) (i=$i_hist), not present in Nom Values")
            skip=true
        end
        if tune_index == [] 
            println("Skip $(nom_histos_names[nom_index[1],1]) (i=$i_hist), not present in Tune Values")
            println("Skip $(ipol[ipol_index,1]) (i=$i_hist), not present in Tune Values")
            skip=true
        end
        if ipol_index == []
            println("Skip $(nom_histos_names[nom_index[1],1]) (i=$i_hist), not present in ipol Values")
            println("Skip $(tune_histos_names[tune_index[1],1]) (i=$i_hist), not present in ipol Values")
            skip=true
        end
        if skip
            continue
        end
        nom_index = nom_index[1]
        tune_index = tune_index[1]
        ipol_index = ipol_index[1]
        if length(findall(x->occursin(x,nom_histos_names[nom_index]),bad_data))>0
            continue
        end
        weight = getweight(LIST,string(nom_histos_names[nom_index]),weight=0)
        if weight == 0
            continue
        end
        data_bins = data_histos[i_hist,2].weights
        nom_bins = nom_hists[nom_index,2].weights
        nom_err = nom_sumw2[nom_index,2]
        tune_bins = tune_hists[tune_index,2].weights
        tune_err = tune_sumw2[tune_index,2]
        
        binning = data_histos[i_hist,2].edges[1]
        binls = (binning[2:end] .- binning[1:end-1])
        
        nom_err = broadcast(sqrt,nom_err ./ binls)
        tune_err = broadcast(sqrt,tune_err ./ binls)
        
        if prop_err != [] && combine_err
            prop_err_index = findall(i -> i == data_histos_names[i_hist][5:end], prop_err.names)
            tune_err = sqrt.(prop_err.prop_err[prop_err_index][1] .^2 + tune_err .^2)
        elseif  prop_err != [] && !combine_err
            prop_err_index = findall(i -> i == data_histos_names[i_hist][5:end], prop_err.names)
            tune_err = prop_err.prop_err[prop_err_index][1]
        end
        
        chisq_nom_i = 0
        chisq_tune_i = 0
        chisq_nom_err_i = 0
        chisq_tune_err_i = 0
        tru_bins = 0
        
        for i_bin in 1:n_bin
            if  (data_histos[i_hist, 2].weights[i_bin] != 0 )
                #&& ccdf(Chisq(ipol[ipol_index, 3][2][i_bin,6]-1),ipol[ipol_index, 3][2][i_bin,5]) > 0.05)
                #&& ipol[ipol_index, 3][2][i_bin, 1] < 25)
                #&& !isapprox(0,ipol[ipol_index, 3][2][i_bin, 1],atol=0.00001)
                tru_bins += 1
                chisq_nom_i += (nom_bins[i_bin] - data_bins[i_bin])^2 / data_bins[i_bin]
                chisq_tune_i += (tune_bins[i_bin] - data_bins[i_bin])^2 / data_bins[i_bin]
                err_sum_dat_nom  = sqrt(nom_err[i_bin]^2 + data_sumw2[i_hist,2][i_bin]^2)
                err_sum_dat_tune = sqrt(tune_err[i_bin]^2 + data_sumw2[i_hist,2][i_bin]^2)
                #chisq_nom_err_i += (nom_bins[i_bin] - data_bins[i_bin])^2 / (err_sum_dat_nom)^2
                #chisq_tune_err_i += (tune_bins[i_bin] - data_bins[i_bin])^2 / (err_sum_dat_tune)^2
                
                if prop_err != []
                    chisq_nom_err_i += ( nom_bins[i_bin] - data_bins[i_bin] )^2 / (err_sum_dat_nom)^2
                    chisq_tune_err_i += (tune_bins[i_bin] - data_bins[i_bin])^2 / (err_sum_dat_tune)^2
                else
                    chisq_nom_err_i += (nom_bins[i_bin] - data_bins[i_bin])^2 / (data_sumw2[i_hist,2][i_bin])^2
                    chisq_tune_err_i += (tune_bins[i_bin] - data_bins[i_bin])^2 / (data_sumw2[i_hist,2][i_bin])^2
                end
            end
        end
        if tru_bins != 0
            push!(chisq_prob_nom,ccdf(Chisq(tru_bins),chisq_nom_err_i))
            push!(chisq_prob_tune,ccdf(Chisq(tru_bins),chisq_tune_err_i))  
        end
        push!(chisq_nom,chisq_nom_i)
        push!(chisq_tune,chisq_tune_i)
        push!(all_names,data_histos_names[i_hist])
    end
    
    mkpath(folder)
    
    if plot 
        maxval = maximum(vcat(chisq_nom...,chisq_tune...))
        stepsize = maxval/20
        binning = collect(0:stepsize:maxval+stepsize)
        histogram(chisq_nom,bins=binning,alpha=0.6,label="$NAME_NOMINAL \n n entries = $(length(chisq_nom)) \n mean = $(round(mean(chisq_nom),digits=3))")
        histogram!(chisq_tune,bins=binning,alpha=0.6,label="$NAME_TUNED \n n entries = $(length(chisq_tune))\n mean = $(round(mean(chisq_tune),digits=3)) ",xlabel="\$\\chi^{2}\$",ylabel="Entries")
        savefig(string(folder,"/chisq_all.pdf"))
        savefig(string(folder,"/chisq_all.png"))
        
        binning = collect(0:0.03:1)
        histogram(chisq_nom[chisq_nom .< 1],bins=binning,alpha=0.6,label="$NAME_NOMINAL \n n entries = $(length(chisq_nom[chisq_nom .< 1])) \n mean = $(round(mean(chisq_nom[chisq_nom .< 1]),digits=3))")
        histogram!(chisq_tune[chisq_tune .< 1],bins=binning,alpha=0.6,label="$NAME_TUNED \n n entries = $(length(chisq_tune[chisq_tune .< 1])) \n mean = $(round(mean(chisq_tune[chisq_tune .< 1]),digits=3)) ",xlabel="\$\\chi^{2}\$",ylabel="Entries")
        savefig(string(folder,"/chisq_all_1.pdf"))
        savefig(string(folder,"/chisq_all_1.png"))
        
        binning = collect(0:0.03:1)
        diffs = sum(chisq_nom .< 1) - sum(chisq_tune .< 1)
        sml_tune = diffs >= 0 ? chisq_tune[chisq_tune .< 1] : sort(chisq_tune)[1:sum(chisq_tune .< 1)+diffs]
        sml_nom = diffs <= 0 ? chisq_nom[chisq_nom .< 1] : sort(chisq_nom)[1:sum(chisq_nom .< 1)-diffs]
        histogram(sml_nom,bins=binning,alpha=0.6,label="$NAME_NOMINAL \n n entries = $(length(sml_nom)) \n mean = $(round(mean(sml_nom),digits=3))",color=:red)
        histogram!(sml_tune,bins=binning,alpha=0.6,label="$NAME_TUNED \n n entries = $(length(sml_tune)) \n mean = $(round(mean(sml_tune),digits=3)) ",xlabel="\$\\chi^{2}\$",ylabel="Entries",color=:blue)
        savefig(string(folder,"/chisq_all_1_adj.pdf"))
        savefig(string(folder,"/chisq_all_1_adj.png"))
        
        stepsize = maximum(chisq_nom)/20
        binning = collect(0:stepsize:maximum(chisq_nom)+stepsize)
        histogram(chisq_nom,bins=binning,label="$NAME_NOMINAL \n n entries = $(length(chisq_nom)) \n mean = $(round(mean(chisq_nom),digits=3))",xlabel="\$\\chi^{2}\$",ylabel="Entries")
        savefig(string(folder,"/chisq_nom.pdf"))
        savefig(string(folder,"/chisq_nom.png"))
        histogram(chisq_tune,bins=binning,label="$NAME_TUNED \n n entries = $(length(chisq_tune))\n mean = $(round(mean(chisq_tune),digits=3)) ",xlabel="\$\\chi^{2}\$",ylabel="Entries")
        savefig(string(folder,"/chisq_tune.pdf"))
        savefig(string(folder,"/chisq_tune.png"))
        
        binning = collect(0:0.03:1)
        histogram(chisq_nom[chisq_nom .< 1],bins=binning,label="$NAME_NOMINAL \n n entries = $(length(chisq_nom[chisq_nom .< 1])) \n mean = $(round(mean(chisq_nom[chisq_nom .< 1]),digits=3))",xlabel="\$\\chi^{2}\$",ylabel="Entries")
        savefig(string(folder,"/chisq_nom_1.pdf"))
        savefig(string(folder,"/chisq_nom_1.png"))
        histogram(chisq_tune[chisq_tune .< 1],bins=binning,label="$NAME_TUNED \n n entries = $(length(chisq_tune[chisq_tune .< 1])) \n mean = $(round(mean(chisq_tune[chisq_tune .< 1]),digits=3)) ",xlabel="\$\\chi^{2}\$",ylabel="Entries")
        savefig(string(folder,"/chisq_tune_1.pdf"))
        savefig(string(folder,"/chisq_tune_1.png"))


        binning = collect(0:1/20:1)
        histogram(chisq_prob_tune,bins=binning,alpha=0.6,label="$NAME_TUNED \n entries = $(length(chisq_prob_tune))\n mean = $(round(mean(chisq_prob_tune),digits=3)) ",xlabel="\$p\$-value",ylabel="Entries",color=:blue)
        histogram!(chisq_prob_nom,bins=binning,alpha=0.6,label="$NAME_NOMINAL \n entries = $(length(chisq_prob_nom)) \n mean = $(round(mean(chisq_prob_nom),digits=3))",color=:yellow)
        #legend("$NAME_TUNED \n entries = $(length(chisq_prob_tune))\n mean = $(round(mean(chisq_prob_tune),digits=3)) ","$NAME_NOMINAL \n entries = $(length(chisq_prob_nom)) \n mean = $(round(mean(chisq_prob_nom),digits=3))",maxrows=2)
        #l = Plots.layout(title="Legend Title", nrow=1)
        #plot!(legend=l)
        #plot!(l)
        Plots.plot!(size=(450,230),grid=true,framestyle=:box,xlim=(0.0,1.0))
        #copy_ticks(plt) xlabel="\$\\chi^{2}\$ prob"
        savefig(string(folder,"/chisq_prob_all.pdf"))
        savefig(string(folder,"/chisq_prob_all.png"))
        

        lowcutexp = -4
        binning = 10 .^ collect((lowcutexp:0.2:0))

        histogram(chisq_prob_tune[chisq_prob_tune .> 10.0^lowcutexp],bins=binning,alpha=0.6,label="$NAME_TUNED \n entries = $(length(chisq_prob_tune[chisq_prob_tune .> 10.0^lowcutexp]))\n mean = $(round(mean(chisq_prob_tune[chisq_prob_tune .> 10.0^lowcutexp]),digits=3)) ",xlabel="\$p\$-value",ylabel="Entries",xaxis=(:log,[10.0^lowcutexp,1]),color=:blue)
        histogram!(chisq_prob_nom[chisq_prob_nom .> 10.0^lowcutexp],bins=binning,alpha=0.6,label="$NAME_NOMINAL \n entries = $(length(chisq_prob_nom[chisq_prob_nom .> 10.0^lowcutexp])) \n mean = $(round(mean(chisq_prob_nom[chisq_prob_nom .> 10.0^lowcutexp]),digits=3))",color=:yellow)
        plot!(size=(450,230),grid=true,framestyle=:box,legend=:topleft)
        savefig(string(folder,"/chisq_prob_all_log.pdf"))
        savefig(string(folder,"/chisq_prob_all_log.png"))
    end
    (names=all_names,pval_nom = chisq_prob_nom, pval_tune = chisq_prob_tune)
end


function plot_corr_matrix(
    samples,
    parameter_names,
    path::String;
    name=""
    )
    thispath = pwd()
    if(!(path in readdir()))
        mkpath(path)
    end

    par_cor = cor(unshaped.(samples))
    npars = length(par_cor[1,:])
    corr = zeros(npars,npars)
    if npars == 8
        for i in 1:8
            for j in 1:8
                i_t = i
                j_t = j
                if i == 3
                    i = 8
                elseif i > 3
                    i = i-1
                end
                if j == 3 
                    j = 8
                elseif j > 3
                    j = j-1
                end
                corr[i_t,j_t] = par_cor[i,j]
            end
        end
    else
        corr = par_cor
    end

    plot()
    heatmap!(corr,
        yflip=true,
        color=:RdBu,
        xticks = (1:npars, parameter_names), 
        yticks = (1:npars, parameter_names), 
        clim = (-1,1),
        #yplotscale=:ln
    )
    for i in 1:npars
        for j in 1:npars
            annotate!(
                tuple(i, j, 
                text(string.(round.(corr',digits=2))[i,j], 10, ifelse(abs(corr[i,j]) < 0.33, :black, :white), rotation = 0, halign=:center, valign=:center)
                )
                )
        end
    end
    Plots.savefig(string(path,"/cor_table",name,".pdf"))
    Plots.savefig(string(path,"/cor_table",name,".png"))
    #plot!()
end

function plot_tuned_obs(ipol;
    NOMINAL_PATH::String="",
    TUNED_PATH::String="",
    folder="./tuned_obs/",
    MCPATH = "/ceph/groups/e4/users/slacagnina/overH70222/grid_inputs_H7_v42/",
    RIVET_DATA_PATH="/ceph/groups/e4/users/slacagnina/Programs/new/Herwig7/share/Rivet/",
    CONSTR = "/ceph/groups/e4/users/slacagnina/overH70222/bad_data.txt",
    LIST = "/ceph/groups/e4/users/slacagnina/overH70222/longlist_raw.txt",
    NAME_NOMINAL="Nominal",
    NAME_TUNED="Tuned",
    islog=false,
    prop_err_tune=[],
    prop_err_nom=[],
    combine_err=true,
    obsnum=-1,
    PDF=false,
    )
    prepath = pwd()
    mkpath(folder)

    nom_hists = get_all_histograms(NOMINAL_PATH)
    nom_sumw2 = get_all_sumw2(NOMINAL_PATH)
    nom_histos_names = nom_hists[:,1]

    tune_hists = get_all_histograms(TUNED_PATH)
    tune_sumw2 = get_all_sumw2(TUNED_PATH)
    tune_histos_names = tune_hists[:,1]

    ipol_histos_names = ipol[:,1]

    analyses_input_file = string(MCPATH, lpad(1, 5, "0"), "/params.in")
    io = open(analyses_input_file, "r")
    file_string = read(io, String)
    file_lbl = split(file_string, "\n")
    analyses_input = [
        split(file_lbl[i])[4] for
        i in 1:length(file_lbl) if occursin("Rivet:Analyses", file_lbl[i])
    ]

    data_histos = Vector{Matrix{Any}}(undef, length(analyses_input))
    data_sumw2 = Vector{Matrix{Any}}(undef, length(analyses_input))
    for i in 1:length(analyses_input)
        data_histos[i] = get_all_ref_histograms_to_names(
            string(RIVET_DATA_PATH, analyses_input[i], ".yoda"),
            ipol_histos_names
        )
        data_sumw2[i] = get_all_ref_sumw2_to_names(
            string(RIVET_DATA_PATH, analyses_input[i], ".yoda"),
            ipol_histos_names
        )
    end

    data_histos = vcat(data_histos...)
    data_sumw2 = vcat(data_sumw2...)
    data_histos_names = data_histos[:, 1]
    n_hist = length(data_histos_names)
    n_bins = [length(data_histos[i_hist, 2].weights) for i_hist = 1:n_hist]
    bad_data=split(read(open(CONSTR,"r"),String),"\n")
    chisq_nom=[]
    chisq_tune=[]
    chisq_red_nom=[]
    chisq_red_tune=[]
    chisq_red_nom_tot=[]
    chisq_red_tune_tot=[]
    chisq_red_tune_stat=[]
    chisq_prob_nom=[]
    chisq_prob_tune=[]
    bad_ps = []
    err_stats=[]

    for i_hist in 1:n_hist
        if obsnum != -1
            if obsnum != i_hist
                continue
            end
        end
        n_bin = n_bins[i_hist]
        nom_index = findall(i -> i == data_histos_names[i_hist][5:end], nom_histos_names)
        tune_index = findall(i -> i == data_histos_names[i_hist][5:end], tune_histos_names)
        ipol_index = findall(i -> i == data_histos_names[i_hist][5:end], ipol[:,1])
        obstitle = data_histos_names[i_hist][5:end]

        skip=false
        if nom_index == []
            println("Skip $(tune_histos_names[tune_index[1],1]) (i=$i_hist), not present in Nom Values")
            println("Skip $(ipol[ipol_index,1]) (i=$i_hist), not present in Nom Values")
            skip=true
        end
        if tune_index == [] 
            println("Skip $(nom_histos_names[nom_index[1],1]) (i=$i_hist), not present in Tune Values")
            println("Skip $(ipol[ipol_index,1]) (i=$i_hist), not present in Tune Values")
            skip=true
        end
        if ipol_index == []
            println("Skip $(nom_histos_names[nom_index[1],1]) (i=$i_hist), not present in ipol Values")
            println("Skip $(tune_histos_names[tune_index[1],1]) (i=$i_hist), not present in ipol Values")
            skip=true
        end
        if skip
            continue
        end
        nom_index = nom_index[1]
        tune_index = tune_index[1]
        ipol_index = ipol_index[1]
        if length(findall(x->occursin(x,nom_histos_names[nom_index]),bad_data))>0
            continue
        end
        weight = getweight(LIST,string(nom_histos_names[nom_index]),weight=0)
        if weight == 0
            continue
        end
        data_bins = data_histos[i_hist,2].weights
        data_err = data_sumw2[i_hist,2]
        nom_bins = nom_hists[nom_index,2].weights
        nom_err = nom_sumw2[nom_index,2]
        tune_bins = tune_hists[tune_index,2].weights
        tune_err = tune_sumw2[tune_index,2]

        data_histo = data_histos[i_hist,2]
        nom_histo = nom_hists[nom_index,2]
        tune_histo = tune_hists[tune_index,2]
        
        binning = data_histo.edges[1]
        binls = (binning[2:end] .- binning[1:end-1])
        binmids = binning[1:end-1] .+ (binls./2)
        
        nom_err = broadcast(sqrt,nom_err ./ binls)
        tune_err = broadcast(sqrt,tune_err ./ binls)
        tune_err_stat = tune_err

        if prop_err_tune != [] && combine_err
            prop_err_tune_index = findall(i -> i == data_histos_names[i_hist][5:end], prop_err_tune.names)
            tune_err = sqrt.(prop_err_tune.prop_err[prop_err_tune_index][1] .^2 + tune_err .^2)
        elseif  prop_err_tune != [] && !combine_err
            prop_err_tune_index = findall(i -> i == data_histos_names[i_hist][5:end], prop_err_tune.names)
            tune_err = prop_err_tune.prop_err[prop_err_tune_index][1]
        end
        if prop_err_nom != [] && combine_err
            prop_err_nom_index = findall(i -> i == data_histos_names[i_hist][5:end], prop_err_nom.names)
            nom_err = sqrt.(prop_err_nom.prop_err[prop_err_nom_index][1] .^2 + nom_err .^2)
        elseif  prop_err_nom != [] && !combine_err
            prop_err_nom_index = findall(i -> i == data_histos_names[i_hist][5:end], prop_err_nom.names)
            nom_err = prop_err_nom.prop_err[prop_err_nom_index][1]
        end
    
        chisq_nom_err_i = 0
        chisq_tune_err_i = 0
        chisq_nom_err_tot_i = 0
        chisq_tune_err_tot_i = 0
        chisq_tune_err_stat_i = 0
        tru_bins = 0

        for i_bin in 1:n_bin
            if  (data_histos[i_hist, 2].weights[i_bin] != 0 )
                tru_bins += 1
                chisq_nom_err_i += (nom_bins[i_bin] - data_bins[i_bin])^2 / (data_sumw2[i_hist,2][i_bin])^2
                chisq_tune_err_i += (tune_bins[i_bin] - data_bins[i_bin])^2 / (data_sumw2[i_hist,2][i_bin])^2
                full_nom_err = sqrt(nom_err[i_bin]^2 + data_sumw2[i_hist,2][i_bin]^2)
                chisq_nom_err_tot_i += (nom_bins[i_bin] - data_bins[i_bin])^2 / (full_nom_err)^2
                full_tune_err = sqrt(tune_err[i_bin]^2 + data_sumw2[i_hist,2][i_bin]^2)
                chisq_tune_err_tot_i += (tune_bins[i_bin] - data_bins[i_bin])^2 / (full_tune_err)^2
                full_err_stat = sqrt(tune_err_stat[i_bin]^2 + data_sumw2[i_hist,2][i_bin]^2)
                chisq_tune_err_stat_i += (tune_bins[i_bin] - data_bins[i_bin])^2 / (full_err_stat)^2
            end
        end
        if tru_bins == 0
            continue
        end

        cd(folder)
        if(!(obstitle in readdir()))
            mkpath(obstitle[2:end])
            #mkpath(string(obstitle[2:end],"/CSV/"))
        end
        cd(prepath)

        push!(chisq_prob_nom,ccdf(Chisq(tru_bins),chisq_nom_err_i))
        push!(chisq_prob_tune,ccdf(Chisq(tru_bins),chisq_tune_err_i))  
        push!(chisq_red_nom,chisq_nom_err_i/(tru_bins))
        push!(chisq_red_tune,chisq_tune_err_i/(tru_bins))
        push!(chisq_red_nom_tot,chisq_nom_err_tot_i/(tru_bins))
        push!(chisq_red_tune_tot,chisq_tune_err_tot_i/(tru_bins))
        push!(chisq_red_tune_stat,chisq_tune_err_stat_i/(tru_bins))

        push!(err_stats,[(broadcast(sqrt,nom_sumw2[tune_index,2] ./ binls)),(broadcast(sqrt,tune_sumw2[tune_index,2] ./ binls)),prop_err_tune.prop_err[prop_err_tune_index][1]])

        xobslabel = data_histos_names[i_hist][5:end]
        if LIST != ""
            xobslabel = getobsname(LIST,string(xobslabel))
        end
        println(xobslabel)
        if ccdf(Chisq(tru_bins),chisq_tune_err_i) < 10^-4
            push!(bad_ps,(obsname=obstitle,pval=ccdf(Chisq(tru_bins),chisq_tune_err_i)))
        end
        
        pltlog,pltrlog = [plot_mc_v_data(binmids=binmids,data_bins=data_bins,data_histo=data_histo,data_err=data_err,nom_bins=nom_bins,nom_histo=nom_histo,nom_err=nom_err,tune_bins=tune_bins,tune_histo=tune_histo,tune_err=tune_err,chisq_prob_nom=chisq_prob_nom,chisq_prob_tune=chisq_prob_tune,xobslabel=xobslabel,pvals=false,ylog=true,ratio=i,lbln=NAME_NOMINAL,lblt=NAME_TUNED,) for i in [false,true]]
        plt,pltr = [plot_mc_v_data(binmids=binmids,data_bins=data_bins,data_histo=data_histo,data_err=data_err,nom_bins=nom_bins,nom_histo=nom_histo,nom_err=nom_err,tune_bins=tune_bins,tune_histo=tune_histo,tune_err=tune_err,chisq_prob_nom=chisq_prob_nom,chisq_prob_tune=chisq_prob_tune,xobslabel=xobslabel,pvals=false,ylog=false,ratio=i,lbln=NAME_NOMINAL,lblt=NAME_TUNED,) for i in [false,true]]

        savefig(pltlog[1],string(folder,obstitle[2:end],"/log_MCvData.png"))
        savefig(pltlog[2],string(folder,obstitle[2:end],"/log_MCvData_2.png"))
        savefig(pltlog[3],string(folder,obstitle[2:end],"/log_MCvData_3.png"))
        savefig(pltlog[4],string(folder,obstitle[2:end],"/log_MCvData_4.png"))
        savefig(pltrlog[1],string(folder,obstitle[2:end],"/log_MCvData_R.png"))
        savefig(pltrlog[2],string(folder,obstitle[2:end],"/log_MCvData_R_2.png"))
        savefig(pltrlog[3],string(folder,obstitle[2:end],"/log_MCvData_RD.png"))
        savefig(pltrlog[4],string(folder,obstitle[2:end],"/log_MCvData_RD_2.png"))
        savefig(pltrlog[5],string(folder,obstitle[2:end],"/log_MCvData_Pull.png"))
        savefig(plt[1],string(folder,obstitle[2:end],"/MCvData.png"))
        savefig(plt[2],string(folder,obstitle[2:end],"/MCvData_2.png"))
        savefig(plt[3],string(folder,obstitle[2:end],"/MCvData_3.png"))
        savefig(plt[4],string(folder,obstitle[2:end],"/MCvData_4.png"))
        savefig(pltr[1],string(folder,obstitle[2:end],"/MCvData_R.png"))
        savefig(pltr[2],string(folder,obstitle[2:end],"/MCvData_R_2.png"))
        savefig(pltr[3],string(folder,obstitle[2:end],"/MCvData_RD.png"))
        savefig(pltr[4],string(folder,obstitle[2:end],"/MCvData_RD_2.png"))
        savefig(pltr[5],string(folder,obstitle[2:end],"/MCvData_Pull.png"))
        
        if PDF
            if occursin("MULTIPLICITIES",obstitle[2:end])
                pltr3 = plot(pltr[3],xticks=false,xlims=(Plots.xlims(pltr[3])[1]*1.001,Plots.xlims(pltr[3])[2]*0.999),legend=:topright,xlabel="   ",size=(374,340))
                pltrlog3 = plot(pltrlog[3],xticks=false,xlims=(12,28),xlabel="   ",size=(374,340))
                #xlims=(Plots.xlims(pltr[3])[1]*1.001,Plots.xlims(pltr[3])[2]*0.999)
                savefig(pltr3,string(folder,obstitle[2:end],"/MCvData_RD.pdf"))
                savefig(pltrlog3,string(folder,obstitle[2:end],"/log_MCvData_RD.pdf"))
            else
                pltr3 = plot(pltr[3],legend=:right,xlabel="   ")
                pltrlog3 = Plots.plot(pltrlog[3],xlim=(0.004,0.85),xscale=:log,size=(374,340))
                savefig(pltr3,string(folder,obstitle[2:end],"/MCvData_RD.pdf"))
                savefig(pltrlog3,string(folder,obstitle[2:end],"/log_MCvData_RD.pdf"))            
            end
        end

    end
    for i in bad_ps
        mkpath(string(folder,"bad_bins_1/"))
        cp(string(folder,i.obsname,"/MCvData.png"),string(folder,"bad_bins_1/",replace(i.obsname,"/"=>"_"),".png"),force=true)
        mkpath(string(folder,"bad_bins_2/"))
        cp(string(folder,i.obsname,"/MCvData_2.png"),string(folder,"bad_bins_2/",replace(i.obsname,"/"=>"_"),".png"),force=true)
        mkpath(string(folder,"bad_bins_1_log/"))
        if exists_str("MCvData_log.png",readdir(string(folder,i.obsname)))
            cp(string(folder,i.obsname,"/MCvData_log.png"),string(folder,"bad_bins_1_log/",replace(i.obsname,"/"=>"_"),".png"),force=true)
        end
        mkpath(string(folder,"bad_bins_2_log/"))
        if exists_str("MCvData_2_log.png",readdir(string(folder,i.obsname)))
            cp(string(folder,i.obsname,"/MCvData_2_log.png"),string(folder,"bad_bins_2_log/",replace(i.obsname,"/"=>"_"),".png"),force=true)
        end
    end
    #return err_stats
    return [chisq_red_nom,chisq_red_tune,chisq_red_nom_tot,chisq_red_tune_stat,chisq_red_tune_tot]
    #return bad_ps
end


function plot_mc_v_data(;
    binmids=[],
    data_bins=[],
    data_histo=[],
    data_err=[],
    nom_bins=[],
    nom_histo=[],
    nom_err=[],
    tune_bins=[],
    tune_histo=[],
    tune_err=[],
    chisq_prob_nom=[],
    chisq_prob_tune=[],
    xobslabel="",
    pvals=false,
    ylog = false,
    ratio = false,
    lbln="",
    lblt="",
    )

    yplotscale = :identity
    if ylog
        if  (length(findall(x->x<=1e-10,(data_bins.-data_err))) == 0 &&
            length(findall(x->x<=1e-10,(nom_bins .- nom_err))) == 0 &&
            length(findall(x->x<=1e-10,(tune_bins.-tune_err))) == 0 &&
            length(findall(x->x<=1e-10,(nom_bins .- sqrt.(nom_err.^2 .+ data_err.^2)))) == 0 &&
            length(findall(x->x<=1e-10,(tune_bins.- sqrt.(tune_err.^2 .+ data_err.^2)))) == 0)
            yplotscale = :log10
        else
            println("WARNING: Could not create log plot, values smaller 0, scale set to lin")
            yplotscale = :identity
        end
    end

    if (lbln=="") lbln="Nominal" end
    if (lbln=="") lblt = "Tuned" end
    if pvals && chisq_prob_nom != [] && chisq_prob_tune != []
        lbln = string(lbln," \n p=$(round(chisq_prob_nom[end],digits=3))")
        lblt = string(lblt," \n p=$(round(chisq_prob_tune[end],digits=3))")
    end

    plt1 = plot(data_histo,seriestype=:step,label="Data",color=:black,legend=:best,xlabel=latexstring(xobslabel),ylabel="Entries",yscale=yplotscale,lw=0)
    #plot!(binmids,data_bins,yerr=data_err,linewidth=0,color=:black,label="",yscale=yplotscale)
    scatter!(binmids,data_bins,yerr=data_err,linewidth=1,ms=3,msw=1,color=:black,label="")
    plot!(nom_histo,seriestype=:step,label=lbln,color=:red,yscale=yplotscale)
    plot!(binmids,nom_bins,linewidth=0,color=:red,label="",yscale=yplotscale,yerr=0)#nom_err)
    scatter!(binmids,nom_bins,linewidth=1,ms=0,color=:red,label="",yerr=0)#nom_err)
    plot!(tune_histo,seriestype=:step,label=lblt,color=:blue,yscale=yplotscale)
    plot!(binmids,tune_bins,linewidth=0,color=:blue,label="",yscale=yplotscale,yerr=0)#tune_err)
    scatter!(binmids,tune_bins,linewidth=1,ms=0,color=:blue,label="",yerr=0)#tune_err)

    plt2 = scatter(binmids,data_bins,lw=0,markersize=3,label="",color=:black,legend=:best,xlabel=latexstring(xobslabel),ylabel="Entries",yscale=yplotscale)
    plot!(binmids,data_bins,lw=0,yscale=yplotscale,label="Data",color=:black)
    plot!(nom_histo,seriestype=:step,label=lbln,color=:red,yscale=yplotscale)
    plot!(binmids,nom_bins,yerr=sqrt.(nom_err.^2 .+ data_err.^2),linewidth=0,color=:red,label="",yscale=yplotscale)
    scatter!(binmids,nom_bins,yerr=sqrt.(nom_err.^2 .+ data_err.^2),linewidth=1,ms=0,color=:red,label="")
    plot!(tune_histo,seriestype=:step,label=lblt,color=:blue,yscale=yplotscale)
    plot!(binmids,tune_bins,yerr=sqrt.(tune_err.^2 .+ data_err.^2),linewidth=0,color=:blue,label="",yscale=yplotscale)
    scatter!(binmids,tune_bins,yerr=sqrt.(tune_err.^2 .+ data_err.^2),linewidth=1,ms=0,color=:blue,label="")

    plt3 = scatter(binmids,data_bins,lw=0,markersize=3,label="",color=:black,legend=:best,xlabel=latexstring(xobslabel),ylabel="Entries",yscale=yplotscale)
    plot!(binmids,data_bins,lw=0,yscale=yplotscale,label="Data",color=:black)
    plot!(nom_histo,seriestype=:step,label=lbln,color=:red,yscale=yplotscale)
    plot!(Histogram(data_histo.edges[1],nom_bins .- sqrt.(nom_err.^2 .+ data_err.^2)),fill_between=nom_bins .+ sqrt.(nom_err.^2 .+ data_err.^2),lw=0,color=:red,alpha=0.3,label="")
    plot!(tune_histo,seriestype=:step,label=lblt,color=:blue,yscale=yplotscale)
    plot!(Histogram(data_histo.edges[1],tune_bins .- sqrt.(tune_err.^2 .+ data_err.^2)),fill_between=tune_bins .+ sqrt.(tune_err.^2 .+ data_err.^2),lw=0,color=:blue,alpha=0.3,label="")
    
    if !occursin("multiplicity",xobslabel) nom_histo.edges[1][1]=0.003;tune_histo.edges[1][1]=0.003;data_histo.edges[1][1]=0.003 end
    plt4 = Plots.scatter(binmids,data_bins,lw=0,markersize=2,label="",color=:black,legend=:best,xlabel=latexstring(xobslabel),ylabel="Entries",yscale=yplotscale)
    Plots.plot!([NaN], [NaN], color=:black, linewidth=1, label="Data")
    Plots.plot!(binmids,data_bins,linecolor=invisible(),yscale=yplotscale,label="",color=:black)
    Plots.plot!(nom_histo,seriestype=:step,label=lbln,color=:red,yscale=yplotscale)
    Plots.plot!(Histogram(data_histo.edges[1],nom_bins .- sqrt.(nom_err.^2)),fill_between=nom_bins .+ sqrt.(nom_err.^2),lw=0,color=:red,alpha=0.3,label="")
    Plots.plot!(tune_histo,seriestype=:step,label=lblt,color=:blue,yscale=yplotscale)
    Plots.plot!(Histogram(data_histo.edges[1],tune_bins .- sqrt.(tune_err.^2)),fill_between=tune_bins .+ sqrt.(tune_err.^2),lw=0,color=:blue,alpha=0.3,label="")
    Plots.scatter!(binmids,data_bins,yerr=data_err,linewidth=1,ms=2,msw=1,color=:black,label="")
    if occursin("multiplicity",xobslabel) Plots.plot!(xlabel="",ylabel=latexstring(xobslabel),ylims=(0.15,0.25),legend=:topright) end
    if !occursin("multiplicity",xobslabel) Plots.plot!(ylabel=L"1/\sigma dS/d\sigma",ylims=(0.003,45)) end
    #Plots.ylims()[2]*0+22)  \frac{1}{\sigma}\frac{d\sigma}{dS} ylims=(0.08,Plots.ylims()[2]*0+0.25)
    
    if ratio == false
        return (plt1 = plt1,plt2 = plt2,plt3=plt3,plt4=plt4)
    end

    plot!(plt1,xformatter=_->"",borderstyle=:box,minorgrid=true,size=(510,320),top_margin=2Plots.mm,bottom_margin=-10Plots.mm,legend=:best)
    plot!(plt2,xformatter=_->"",borderstyle=:box,minorgrid=true,size=(510,320),top_margin=2Plots.mm,bottom_margin=-10Plots.mm,legend=:best)
    plot!(plt3,xformatter=_->"",borderstyle=:box,minorgrid=true,size=(510,320),top_margin=2Plots.mm,bottom_margin=-10Plots.mm,legend=:best)
    Plots.plot!(plt4,xformatter=_->"",borderstyle=:box,minorgrid=true,size=(510,320),top_margin=2Plots.mm,bottom_margin=-10Plots.mm,legend=:best)

    dtratio = ratio_histogram(data_histo,tune_histo)
    dtratiovals = dtratio.weights
    dtratioerrs = [i.err for i in measurement.(data_bins,data_err) ./ tune_bins]
    fitratioerrs = [i.err for i in measurement.(tune_bins,tune_err) ./ tune_bins]
    fitratioerrshist = Histogram(data_histo.edges[1],1 .- fitratioerrs)
    binl = (data_histo.edges[1][2:end] .- data_histo.edges[1][1:end-1])./2
    ntratio = ratio_histogram(nom_histo,tune_histo)
    ntratiovals = ntratio.weights
    ntratioerrs = [i.err for i in measurement.(nom_bins,nom_err) ./ tune_bins]

    limval = maximum(vcat(
        abs.(1 .- (dtratiovals .+ dtratioerrs))...,
        abs.(1 .- (dtratiovals .- dtratioerrs))...,        
        abs.(1 .- ntratio.weights)...
    )) * 1.10

    pltr1 = plot(binmids,ones(length(binmids)),link=:all,lw=0,yminorticks=4,xminorticks=2,minorgrid=true,borderstyle=:box,ylim=(1-limval,1+limval),color=:blue,label="")
    hline!([1],color=:blue,label="")
    plot!(fitratioerrshist,fill_between=1 .+ fitratioerrs,linecolor=:blue,color=:blue,alpha=0.3,label="",lw=0)
    scatter!(binmids,dtratiovals,color=:black,label="",xerr=0,yerr=dtratioerrs,lw=0.5,ms=3)
    scatter!(binmids,ntratio.weights,color=:red,label="",xerr=0,yerr=0,lw=0,ms=3,msw=0)
    xlabel!(latexstring(xobslabel))
    ylabel!("Ratio")
    plot!(top_margin=-1.9Plots.mm)
    l = @layout [a{0.7h};b{0.3h}]
    pltr1 = plot(plt1,pltr1,layout=l,size=(510,600),link=:x)
    
    limval = maximum(vcat(
        abs.(1 .- (dtratiovals .+ dtratioerrs))...,
        abs.(1 .- (dtratiovals .- dtratioerrs))...,        
        abs.(1 .- (ntratiovals .+ ntratioerrs))...,
        abs.(1 .- (ntratiovals .- ntratioerrs))...,
        abs.(1 .- ntratio.weights)...
    )) * 1.10
    
    pltr2 = plot(binmids,ones(length(binmids)),link=:all,lw=0,yminorticks=4,xminorticks=2,minorgrid=true,borderstyle=:box,ylim=(1-limval,1+limval),color=:blue,label="")
    hline!([1],color=:blue,label="")
    plot!(fitratioerrshist,fill_between=1 .+ fitratioerrs,linecolor=:blue,color=:blue,alpha=0.3,label="",lw=0)
    scatter!(binmids,dtratiovals,color=:black,label="",xerr=0,yerr=dtratioerrs,lw=0.5,ms=3)
    scatter!(binmids,ntratio.weights,color=:black,label="",xerr=0,yerr=ntratioerrs,lw=1,ms=2.5,msw=1)
    scatter!(binmids,ntratio.weights,color=:red,label="",xerr=0,yerr=ntratioerrs,lw=0,ms=3,msw=0)
    xlabel!(latexstring(xobslabel))
    ylabel!("Ratio")
    plot!(top_margin=-1.9Plots.mm)
    l = @layout [a{0.7h};b{0.3h}]
    pltr2 = plot(plt3,pltr2,layout=l,size=(510,600),link=:x)


    tdratio = ratio_histogram(tune_histo,data_histo)
    ndratio = ratio_histogram(nom_histo,data_histo)
    tdratiovals = tdratio.weights
    ndratiovals = ndratio.weights
    tdratioerrs = [i.err for i in measurement.(tune_bins,tune_err) ./ data_bins]
    ndratioerrs = [i.err for i in measurement.(nom_bins,nom_err) ./ data_bins]
    
    fitratioerrs = [i.err for i in measurement.(data_bins,data_err) ./ data_bins]
    fitratioerrshist = Histogram(data_histo.edges[1],1 .- fitratioerrs)
    binl = (data_histo.edges[1][2:end] .- data_histo.edges[1][1:end-1])./2
    limval = maximum(vcat(
        abs.(1 .- (tdratiovals .+ tdratioerrs))...,
        abs.(1 .- (tdratiovals .- tdratioerrs))...,
        abs.(1 .- (ndratiovals .+ ndratioerrs))...,
        abs.(1 .- (ndratiovals .- ndratioerrs))...,
        abs.(1 .- tdratio.weights)...
    ))*1.20

    pltr3 = Plots.plot(binmids,ones(length(binmids)),link=:all,lw=0,yminorticks=4,xminorticks=2,minorgrid=true,borderstyle=:box,ylim=(1-limval,1+limval),color=:black,label="")
    Plots.hline!([1],color=:black,label="")
    Plots.plot!(fitratioerrshist,fill_between=1 .+ fitratioerrs,linecolor=:black,color=:black,alpha=0.3,label="",lw=0)
    Plots.scatter!(binmids,tdratiovals,color=:black,label="",xerr=0,yerr=tdratioerrs,lw=1,ms=1.5,msw=1)
    Plots.scatter!(binmids,tdratiovals,color=:blue,label="",xerr=0,yerr=tdratioerrs,lw=0,ms=2,msw=0)
    Plots.scatter!(binmids,ndratiovals,color=:black,label="",xerr=0,yerr=ndratioerrs,lw=1,ms=1.5,msw=1)
    Plots.scatter!(binmids,ndratiovals,color=:red,label="",xerr=0,yerr=ndratioerrs,lw=0,ms=2,msw=0)
    Plots.xlabel!(latexstring(xobslabel))
    Plots.ylabel!("Ratio")
    Plots.plot!(top_margin=-1.9Plots.mm,bottom_margin=-2Plots.mm)
    l = @layout [a{0.7h};b{0.3h}]
    pltr3 = Plots.plot(plt4,pltr3,layout=l,size=(374,400),link=:x,left_margin=-0Plots.mm) #340xnorm

    pltr4 = plot(binmids,ones(length(binmids)),link=:all,lw=0,yminorticks=4,xminorticks=2,minorgrid=true,borderstyle=:box,ylim=(1-limval,1+limval),color=:black,label="")
    hline!([1],color=:black,label="")
    plot!(fitratioerrshist,fill_between=1 .+ fitratioerrs,linecolor=:black,color=:black,alpha=0.3,label="",lw=0)
    scatter!(binmids,tdratiovals,color=:black,label="",xerr=0,yerr=tdratioerrs,lw=1,ms=2.5,msw=1)
    scatter!(binmids,tdratiovals,color=:blue,label="",xerr=0,yerr=tdratioerrs,lw=0,ms=3,msw=0)
    scatter!(binmids,ndratiovals,color=:red,label="",xerr=0,yerr=0,lw=0.5,ms=3,msw=0)
    xlabel!(latexstring(xobslabel))
    ylabel!("Ratio")
    plot!(top_margin=-1.9Plots.mm)
    l = @layout [a{0.7h};b{0.3h}]
    pltr4 = plot(plt1,pltr4,layout=l,size=(510,600),link=:x)

    pulltune = abs.(data_bins .- tune_bins) ./ sqrt.(data_err.^2 .+ tune_err.^2)
    pullnom = abs.(data_bins .- nom_bins) ./ sqrt.(data_err.^2 .+ nom_err.^2)
    limval = maximum(vcat(pulltune...,pullnom...))*1.05

    pltr5 = plot(binmids,ones(length(binmids)),link=:all,lw=0,yminorticks=4,xminorticks=2,minorgrid=true,borderstyle=:box,ylim=(1-limval,1+limval),color=:black,label="")
    hline!([1],color=:black,label="")
    #plot!(fitratioerrshist,fill_between=1 .+ fitratioerrs,linecolor=:black,color=:black,alpha=0.3,label="",lw=0)
    scatter!(binmids,pulltune,color=:blue,label="",xerr=0,yerr=0,lw=0.5,ms=3,msw=0)
    scatter!(binmids,pullnom,color=:red,label="",xerr=0,yerr=0,lw=0.5,ms=3,msw=0)
    xlabel!(latexstring(xobslabel))
    ylabel!("Pull")
    plot!(top_margin=-1.9Plots.mm)
    l = @layout [a{0.7h};b{0.3h}]
    pltr5 = plot(plt1,pltr5,layout=l,size=(510,600),link=:x)

    return (plt1=pltr1,plt2=pltr2,plt3=pltr3,plt4=pltr4,plt5=pltr5)
end


function plot_two_hist_with_ratio(h1,h2;
    h1err=[],h2err=[],ylog=true,h1label="",h2label="",rylim=[],xlabel="",xlim=[],xplotscale=:identity,ylim=[],yscale=:identity
    )
    h1vals = h1.weights
    h2vals = h2.weights
    binmids = h1.edges[1][1:end-1] .+ (h1.edges[1][2:end] .- h1.edges[1][1:end-1]) ./ 2
    if h1err == []
        h1err = zeros(length(h1vals))
    end
    if h2err == []
        h2err = zeros(length(h2vals))
    end
    if ylog
        if  (length(findall(x->x<=1e-10,(h1vals.-h1err))) == 0 &&
            length(findall(x->x<=1e-10,(h2vals.-h2err))) == 0)
            yplotscale = :log10
        else
            println("WARNING: Could not create log plot, values smaller 0, scale set to lin")
            yplotscale = :identity
        end
    end

    plt = plot(h1,seriestype=:step,color=:black,yscale=yplotscale)
    plot!(Histogram(h1.edges[1],h1vals .- (h1err)),fill_between=h1vals .+ (h1err),lw=0,color=:black,alpha=0.3,label=h1label)
    scatter!(binmids,h2.weights,lw=0,markersize=2,color=:blue,yscale=yplotscale,yerr=h2err,label=h2label)
    
    plot!(plt,xformatter=_->"",borderstyle=:box,minorgrid=true,size=(510,320),top_margin=2Plots.mm,bottom_margin=-7Plots.mm,legend=:best,xscale=xplotscale)
    if xlim != []
        plot!(xlim=xlim)
    end
    if ylim != []
        plot!(ylim=ylim)
    end
    if yscale != :identity
        plot!(yscale=yscale)
    end

    ratio = ratio_histogram(h2,h1)
    ratiovals = ratio.weights
    h2ratioerrs = [i.err for i in measurement.(h2vals,h2err) ./ h2vals]
    h1ratioerrs = [i.err for i in measurement.(h1vals,h1err) ./ h1vals]
    limval = maximum(vcat(
        abs.(1 .- (ratiovals .+ h2ratioerrs))...,
        abs.(1 .- (ratiovals .- h2ratioerrs))...,
        abs.(1 .- (ratiovals .+ h1ratioerrs))...,
        abs.(1 .- (ratiovals .- h1ratioerrs))...,
        abs.(1 .- ratio.weights)...
    ))*1.20

    h1errhist = Histogram(ratio.edges[1],1 .- h1ratioerrs)

    pltr = plot(binmids,ones(length(binmids)),link=:all,lw=0,yminorticks=4,xminorticks=2,minorgrid=true,borderstyle=:box,ylim=(1-limval,1+limval),color=:black,label="")
    hline!([1],color=:black,label="")
    scatter!(binmids,ratiovals,color=:black,label="",xerr=0,yerr=h2ratioerrs,lw=1,ms=1.5,msw=1)
    scatter!(binmids,ratiovals,color=:blue,label="",xerr=0,yerr=h2ratioerrs,lw=0,ms=2,msw=0)
    plot!(h1errhist,fill_between=1 .+ h1ratioerrs,linecolor=:black,color=:black,alpha=0.3,label="",lw=0,xscale=xplotscale,xlabel=xlabel)
    if rylim !=[]
        plot!(ylim=rylim)
    end
    if xlim != []
        plot!(xlim=xlim)
    end    

    l = @layout [a{0.7h};b{0.3h}]
    pltf = plot(plt,pltr,layout=l,size=(374,400),link=:x,left_margin=-0Plots.mm) #340xnorm
    return pltf
end

