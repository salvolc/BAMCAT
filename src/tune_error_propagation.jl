function propagate_error(samples,
    ipol;
    g_func=g_cubic,
    n=1e6,
    keys_all=[:p1,:p2,:p3,:p4,:p5,:p6,:p7,:p8],
    p=0.68,
    atol=0.05,
    nbins=200,
    parallel=true,
    )
    obs_names = ipol[:,1]
    prop_error = Vector{Vector{Float64}}(undef,length(obs_names))
    prop_error_hists = Vector{Vector{Histogram}}(undef,length(obs_names))
    prop_error_asy = Vector{Vector{Vector{Float64}}}(undef,length(obs_names))
    resamples_full = bat_sample(samples,OrderedResampling(nsamples=n)).result
    resamples = collect(resamples_full.v)
    keys=keys_all[1:length(resamples[1])]
    result = []
    if parallel
        result = Folds.collect(propagate_error_for_iobs(resamples,ipol,keys,g_func,iobs,p=p,atol=atol,nbins=nbins) for iobs in 1:length(obs_names))
    else
        result = [propagate_error_for_iobs(resamples,ipol,keys,g_func,iobs,p=p,atol=atol,nbins=nbins) for iobs in 1:length(obs_names)]
    end
    prop_error = [i[1] for i in result]
    prop_error_asy = [i[2] for i in result]
    prop_error_hists = [i[3] for i in result]
    #hcat(obs_names[1:4],prop_error,prop_error_asy,prop_error_hists)
    (names = obs_names,prop_err=prop_error,prop_err_int=prop_error_asy,prop_err_hists = prop_error_hists)
end

function propagate_error_for_iobs(resamples,ipol,keys,g_func,iobs;p=0.68,atol=0.05,nbins=200)
    println(ipol[:,1][iobs])
    binedges = ipol[iobs,2]
    nbins_obs = length(binedges)-1
    prop_error_iobs = Vector{Float64}(undef,nbins_obs)
    prop_error_hists = Vector{Histogram}(undef,nbins_obs)
    prop_error_asy_iobs = Vector{Vector{Float64}}(undef,nbins_obs)
    for ibin in 1:nbins_obs
        coeffs = ipol[iobs,3][1][ibin,:]
        pars = [[i[keyi] for keyi in keys] for i in resamples]
        vals = [g_func(ipars,coeffs) for ipars in pars]
        h1 = fit(Histogram,vals,nbins=nbins)
        hist_p = BAT.get_smallest_intervals(h1, [p])[1][1]
        lower, upper = BAT.get_interval_edges(hist_p, atol=atol)
        #@assert(length(lower) > 1,"ERROR: HOW DID THIS HAPPEN; ERRORS ON PROPAGATION SHOULD BE GAUSSIAN")
        if length(lower) > 1
            println("ERROR: WAS IST HIER LOS!!!")
        end
        prop_error_iobs[ibin] =  std(vals)#(upper[1]-lower[1])/2 using 
        prop_error_hists[ibin] = h1
        allul = [[lower[i],upper[i]] for i in 1:length(lower)]
        prop_error_asy_iobs[ibin] = vcat(allul...)
    end
    (prop_error_iobs,prop_error_asy_iobs,prop_error_hists)
end

function filter_err_prop(prop_err;LIST="/ceph/groups/e4/users/slacagnina/overH70222/longlist_raw.txt")
    all = hcat([prop_err[i] for i in 1:length(prop_err)]...)
    indx = findall(x->getweight(LIST,x,weight=0)!=0,prop_err.names)
    all_filt = all[indx,:]
    prop_err_filt = (names=all_filt[:,1],prop_err=all_filt[:,2],prop_err_int=all_filt[:,3],prop_err_hists=all_filt[:,4])
end

function give_plots_err_prop(prop_err_filt;LIST="/ceph/groups/e4/users/slacagnina/overH70222/longlist_raw.txt")
    plts = []
    plts_obs = []
    for iobs in 1:length(prop_err_filt.names)
        xlab = getobsname(LIST,prop_err_filt.names[iobs])
        pobs=[]
        for ibin in 1:length(prop_err_filt.prop_err[iobs])
            mid = prop_err_filt.prop_err_hists[iobs][ibin].edges[1][findmax(prop_err_filt.prop_err_hists[iobs][ibin].weights)[2]]
            plt = plot(prop_err_filt.prop_err_hists[iobs][ibin],xlim=(mid-5*prop_err_filt.prop_err[iobs][ibin],mid+5*prop_err_filt.prop_err[iobs][ibin]),label="Parametrization",title="Observable $(xlab), Bin: $(ibin)",st=:step,yformatter=:plain)
            yt = (Plots.yticks(plt)[1][1],string.(parse.(Float64,Plots.yticks(plt)[1][2])./10000))
            plot!(yticks=yt)
            push!(plts,plt)
            push!(pobs,plt)
        end
        push!(plts_obs,(names=prop_err_filt.names[iobs],plts=pobs))
    end
    (plts=plts,plts_obs=plts_obs)
end

function give_plots_err_prop_for_obs(prop_err_filt,obsname)
    plts,plts_obs = give_plots_err_prop(prop_err_filt;LIST="/ceph/groups/e4/users/slacagnina/overH70222/longlist_raw.txt")
    vcat([i.plts for i in plts_obs[give_all_str(obsname,prop_err_filt.names)]]...)
end

function plot_err_prop_for_obs(prop_err_filt,obsname;folder="./err_prop_bins/",PDF=false)
    indx = give_all_str(obsname,prop_err_filt.names)
    plts, plts_obs = give_plots_err_prop(prop_err_filt)
    for i in indx
        fname = getobsname("/ceph/groups/e4/users/slacagnina/overH70222/longlist_raw.txt",prop_err_filt.names[i])
        fname = replace(replace(replace(fname,"\\"=>""),"\$"=>""),"/"=>"")
        mkpath(string(folder,obsname))
        pltbins = plts_obs[i].plts
        nbins = length(pltbins)
        for j in 1:nbins
            plt = pltbins[j]
            plot!(pltbins[j],size=(450,250),titlefont=("Helevetica",10),framestyle=:box,ylabel=string("Entries / ",L"10^4"),xlabel="Bin value",st=:stephist,formatter=:plain)
            xl = Plots.xlims(pltbins[j])
            yl = Plots.ylims(pltbins[j])
            #annotate!([(xl[1]*1.02, yl[2] * 1.051, Plots.text(L"\times10^{4}", 11, :black, :center))])
            
            Plots.savefig(pltbins[j],string(folder,obsname,"/",fname,"_Bin_$(j)",".png"))
            if PDF
                Plots.savefig(pltbins[j],string(folder,obsname,"/",fname,"_Bin_$(j)",".pdf"))
            end
        end
    end
end