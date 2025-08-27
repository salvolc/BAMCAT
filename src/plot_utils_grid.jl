####################################################
####This Plots all fits obs wise min max ranges ####
####################################################,
function plot_fit_grid(ipol,grid::String; folder::String="fit_grid")
    gridhc = load_Histcontainer_from_folder(grid)
    gridsorted = sort_Histcontainer_by_observable(gridhc)
    plot_fit_grid(ipol,gridsorted,folder)
end

function plot_fit_grid(ipol,grid::Vector{Observablecontainer}; folder::String="fit_grid", par_mask=[true for i in 1:length(grid[1].hc[1].parname)],LIST="",cubic=false)
    prepath = pwd()
    mkpath(folder)
    hist_names = ipol[:,1]
    pars = get_clean_params_from_mc(grid)
    nperpar = Int(length(pars[:,1])/length(pars[1,:]))

    delindst = collect(1:nperpar:length(pars[:,1]))[findall(x->x==0,par_mask)]
    delinded = collect(1:nperpar:length(pars[:,1]))[findall(x->x==0,par_mask)] .+ (nperpar-1)

    delindst_sort = sort(delindst,rev=true)
    delinded_sort = sort(delinded,rev=true)

    for i in zip(delindst_sort,delinded_sort)
        pars = pars[1:end .∉ [collect(i[1]:i[2])] , 1:end]
    end
    pars = pars[:,par_mask]
    parnames = grid[1].hc[1].parname[par_mask]

    for iobs in 1:length(hist_names)
        obs = hist_names[iobs]
        if occursin("RAW",obs)
            continue
        end
        if LIST != ""
            if getweight(LIST,string(obs),weight=0) == 0
                continue
            end
        end
        binvalues = get_bin_value_for_observable(grid,obs)
        binedges = get_bin_edges_for_observables(grid,obs)
        binerrs = get_sumw2_for_obervable(grid,obs)
        for i in zip(delindst_sort,delinded_sort)
            binvalues = binvalues[1:end .∉ [collect(i[1]:i[2])] , 1:end]
            binerrs = binerrs[1:end .∉ [collect(i[1]:i[2])] , 1:end]
        end
        nparams = length(pars[1,:])
        nbins = length(binedges)-1
        nparcombos = length(pars[:,1])
        nvalsperpar = Int(length(pars[:,1])/nparams)
        parvals = pars[:,:]
        coeffs = ipol[:,3][iobs][1]
        if cubic
            fitvals = [g_cubic(parvals,coe) for coe in eachrow(coeffs)]
        else
            fitvals = [g(parvals,coe) for coe in eachrow(coeffs)]
        end

        for ibin in 1:nbins
            if sum(binerrs[:,ibin+2] .< 0) > 0
                print("Skipped $obs in bin $ibin for negative sumw2.")
                continue
            end
            binl = binedges[ibin+1]-binedges[ibin]
            for ipar in 1:nparams

                parbeg = 1+nvalsperpar*(ipar-1)
                parend = nvalsperpar*ipar
                parval = parvals[parbeg:parend,:][:,ipar]
                fitval = fitvals[ibin][parbeg:parend]
                orgval = binvalues[parbeg:parend,2+ibin]
                orgvalerr = broadcast(sqrt,binerrs[parbeg:parend,2+ibin] ./ binl) 

                c2 = chisqerr(fitval,orgval,orgvalerr)
                plot(parval,orgval,yerr=orgvalerr,color="grey",label="MC")
                plot!(parval,fitval,label="Fit \n chi2 = $c2")
                par_name = parnames[ipar]
                if occursin("/",par_name)
                    par_name = split(par_name,"/")[end]
                end
                xlabel!("$par_name")
                ylabel!("Entries")
                title!("$(obs), Bin: $(ibin)")

                cd(folder)
                if(!(obs in readdir()))
                    mkpath(obs[2:end])
                end
                cd(prepath)
                #savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$par_name-$(ibin)-2.pdf"))
                savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$par_name-$(ibin)-2.png"))
            end
        end
    end
end

function plot_fit_grid_pulls_2(ipol,grid::Vector{Observablecontainer}; folder::String="fit_grid", par_mask=[true for i in 1:length(grid[1].hc[1].parname)],LIST="",cubic=false,pdf=false,inbound=false,pro=0.10)
    prepath = pwd()
    mkpath(folder)
    hist_names = ipol[:,1]
    pars = get_clean_params_from_mc(grid)
    nperpar = Int(length(pars[:,1])/length(pars[1,:]))

    delindst = collect(1:nperpar:length(pars[:,1]))[findall(x->x==0,par_mask)]
    delinded = collect(1:nperpar:length(pars[:,1]))[findall(x->x==0,par_mask)] .+ (nperpar-1)

    delindst_sort = sort(delindst,rev=true)
    delinded_sort = sort(delinded,rev=true)

    for i in zip(delindst_sort,delinded_sort)
        pars = pars[1:end .∉ [collect(i[1]:i[2])] , 1:end]
    end
    pars = pars[:,par_mask]
    parnames = grid[1].hc[1].parname[par_mask]

    def_pars = [pars[end,1],pars[1,2:end]...]
    def_pars_vars = Array{Float64,2}(undef,length(def_pars)*2+2,length(def_pars))

    for i in 1:length(def_pars_vars[:,1])
        par_id = i%length(def_pars)+1
        sign = i > length(def_pars) ? +1 : -1
        def_pars_vars[i,:] = def_pars
        def_pars_vars[i,par_id] = (1+sign*pro)*def_pars[par_id]
    end
    def_pars_vars[end-1,:] = (1+pro) .* def_pars
    def_pars_vars[end,:] = (1-pro) .* def_pars

    if length(ipol[1,3]) < 3
        println("Ipol object does not contain covar matrix.")
    end

    for iobs in 1:length(hist_names)
        obs = hist_names[iobs]
        if occursin("RAW",obs)
            continue
        end
        if LIST != ""
            if getweight(LIST,string(obs),weight=0) == 0
                continue
            end
        end
        binvalues = get_bin_value_for_observable(grid,obs)
        binedges = get_bin_edges_for_observables(grid,obs)
        binerrs = get_sumw2_for_obervable(grid,obs)
        for i in zip(delindst_sort,delinded_sort)
            binvalues = binvalues[1:end .∉ [collect(i[1]:i[2])] , 1:end]
            binerrs = binerrs[1:end .∉ [collect(i[1]:i[2])] , 1:end]
        end
        nparams = length(pars[1,:])
        nbins = length(binedges)-1
        nparcombos = length(pars[:,1])
        nvalsperpar = Int(length(pars[:,1])/nparams)
        parvals = pars[:,:]
        coeffs = ipol[:,3][iobs][1]
        if cubic
            fitvals = [g_cubic(parvals,coe) for coe in eachrow(coeffs)]
        else
            fitvals = [g(parvals,coe) for coe in eachrow(coeffs)]
        end

        if cubic
            def_fitvals_vars = [g_cubic(def_pars_vars,coe) for coe in eachrow(coeffs)]
        else
            def_fitvals_vars = [g(def_pars_vars,coe) for coe in eachrow(coeffs)]
        end

        
        for ibin in 1:nbins
            if sum(binerrs[:,ibin+2] .< 0) > 0
                print("Skipped $obs in bin $ibin for negative sumw2.")
                continue
            end
            binl = binedges[ibin+1]-binedges[ibin]
            for ipar in 1:nparams

                def_fitval_var = def_fitvals_vars[ibin]
                def_fitval_var_ind = [argmin(def_fitval_var),argmax(def_fitval_var)]
                def_parval_low = def_pars_vars[argmin(def_fitval_var),:]
                def_parval_high = def_pars_vars[argmax(def_fitval_var),:]
                
                parbeg = 1+nvalsperpar*(ipar-1)
                parend = nvalsperpar*ipar
                parval = parvals[parbeg:parend,:][:,ipar]
                fitval = fitvals[ibin][parbeg:parend]
                orgval = binvalues[parbeg:parend,2+ibin]
                orgvalerr = broadcast(sqrt,binerrs[parbeg:parend,2+ibin] ./ binl) 
                covvar = ipol[iobs,3][3][ibin]

                def_parvals_low = Matrix{Float64}(undef,length(parval),length(def_parval_low))
                def_parvals_high = Matrix{Float64}(undef,length(parval),length(def_parval_low))

                for i in 1:length(parval)
                    def_parvals_low[i,:] = def_parval_low
                    def_parvals_high[i,:] = def_parval_high
                    def_parvals_low[i,ipar] = parval[i]
                    def_parvals_high[i,ipar] = parval[i]
                end
                
                
                parvecs = parvals[parbeg:parend,:]
                
                if cubic
                    fiterrs = [sqrt(g_err_cubic(parvecs[i,:],covvar)) for i in 1:length(parvecs[:,1])]
                    default_vars_low = g_cubic(def_parvals_low,coeffs[ibin,:])
                    default_vars_high = g_cubic(def_parvals_high,coeffs[ibin,:])
                else
                    fiterrs = [sqrt(g_err(parvecs[i,:],covvar)) for i in 1:length(parvecs[:,1])]
                    default_vars_low = g(def_parvals_low,coeffs[ibin,:])
                    default_vars_high = g(def_parvals_high,coeffs[ibin,:])
                end

                dfl = [default_vars_low[i] < default_vars_high[i] ? default_vars_low[i] : default_vars_high[i] for i in 1:length(default_vars_low)]
                dfh = [default_vars_low[i] > default_vars_high[i] ? default_vars_low[i] : default_vars_high[i] for i in 1:length(default_vars_low)]

                if inbound
                    if ipar == 1 && nparams == 6
                        parval=parval[1:end-6]
                        fitval=fitval[1:end-6]
                        fiterrs=fiterrs[1:end-6]
                        orgval=orgval[1:end-6]
                        orgvalerr=orgvalerr[1:end-6]
                        dfl=dfl[1:end-6]
                        dfh=dfh[1:end-6]
                    end
                    if ipar == 1 && nparams == 8
                        parval=parval[1:end-3]
                        fitval=fitval[1:end-3]
                        fiterrs=fiterrs[1:end-3]
                        orgval=orgval[1:end-3]
                        orgvalerr=orgvalerr[1:end-3]
                        dfl=dfl[1:end-3]
                        dfh=dfh[1:end-3]
                    end
                    if ipar == 2 && nparams == 8
                        parval=parval[1:end-5]
                        fitval=fitval[1:end-5]
                        fiterrs=fiterrs[1:end-5]
                        orgval=orgval[1:end-5]
                        orgvalerr=orgvalerr[1:end-5]
                        dfl=dfl[1:end-5]
                        dfh=dfh[1:end-5]
                    end
                end
                                
                respull = (fitval .- orgval) ./ sqrt.(abs.(orgvalerr .^ 2 .- fiterrs .^ 2))
                chisqval = (fitval .- orgval).^2 ./ ((orgvalerr .^ 2 .+ fiterrs .^ 2))
        

                c2 = chisqerr(fitval,orgval,orgvalerr)
                plot(parval,orgval,yerr=orgvalerr,color="grey",label="MC")
                plot!(parval,fitval,label="Fit")#label="Fit \n chi2 = $c2")
                plot!(parval,fitval .- fiterrs,fill_between=fitval .+ fiterrs,linecolor=:red,color=:red,alpha=0.3,label="Fit Err")
                plot!(parval,dfl,fill_between=dfh,linecolor=:blue,color=:blue,alpha=0.3,label="Default Var")
                vline!([def_pars[ipar]],label="Default",color="grey")

                par_name = parnames[ipar]
                if occursin("/",par_name)
                    par_name = split(par_name,"/")[end]
                end
                if LIST != ""
                    obstitle = getobsname(LIST,obs)
                    title!("$(obstitle), Bin: $(ibin)")
                else
                    title!("$(obs), Bin: $(ibin)")
                end
                ylabel!("Entries")
                xlabel!("$par_name")

                cd(folder)
                if(!(obs in readdir()))
                    mkpath(obs[2:end])
                    mkpath(string(obs[2:end],"/CSV/"))
                end
                cd(prepath)

                df = DataFrame(xval=parval, yval_orgval=orgval, yval_orgvalerr=orgvalerr, yval_fitval=fitval, yval_fitvalerr=fiterrs
                )
                CSV.write(string(folder,obs,"/CSV/",split(obs,"/")[end-1],"-$par_name-$(ibin)-2.csv"),df)
                

                if pdf 
                    savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$par_name-$(ibin)-2.pdf"))
                end
                savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$par_name-$(ibin)-2.png"))

                scatter(parval,respull,yerr=1,color="black",label="Pull",ms=3,msc="black")
                xlabel!("$par_name")
                ylabel!("Res Pull")
                title!("$(obs), Bin: $(ibin)")
                if pdf
                    savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$par_name-$(ibin)-Pull.pdf"))
                end
                savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$par_name-$(ibin)-Pull.png"))
                
                hchisqpull = fit(Histogram,chisqval)
                if cubic
                    funcchi = Chisq(length(chisqval)-4)
                else
                    funcchi = Chisq(length(chisqval)-3)
                end
                chiprob = ccdf(funcchi,sum(chisqval))

                plot(hchisqpull,label="chisq ,n=$(length(chisqval)) \n chiprob = $(chiprob)",seriestype=:step)
                xlabel!("$par_name")
                ylabel!("Chisq Pull")
                title!("$(obs), Bin: $(ibin)")
                if pdf
                    savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$par_name-$(ibin)-chisq.pdf"))
                end
                savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$par_name-$(ibin)-chisq.png"))
            end
        end
    end
end

function calc_fit_grid_pulls(ipol,grid::Vector{Observablecontainer}; skip_neg=false, folder::String="fit_grid", par_mask=[true for i in 1:length(grid[1].hc[1].parname)],LIST="",cubic=false,inbound=false)
    prepath = pwd()
    mkpath(folder)
    hist_names = ipol[:,1]
    pars = get_clean_params_from_mc(grid)
    nperpar = Int(length(pars[:,1])/length(pars[1,:]))

    delindst = collect(1:nperpar:length(pars[:,1]))[findall(x->x==0,par_mask)]
    delinded = collect(1:nperpar:length(pars[:,1]))[findall(x->x==0,par_mask)] .+ (nperpar-1)

    delindst_sort = sort(delindst,rev=true)
    delinded_sort = sort(delinded,rev=true)

    for i in zip(delindst_sort,delinded_sort)
        pars = pars[1:end .∉ [collect(i[1]:i[2])] , 1:end]
    end
    pars = pars[:,par_mask]
    parnames = grid[1].hc[1].parname[par_mask]

    if length(ipol[1,3]) < 3
        println("Ipol object does not contain covar matrix.")
    end

    pulls = []
    pulls_per_par = [[] for i in 1:length(pars[1,:])]
    chiprob = []
    chiprob_per_par = [[] for i in 1:length(pars[1,:])]
    chisqvals = []
    chisqvals_per_par = [[] for i in 1:length(pars[1,:])]
    iden = []

    for iobs in 1:length(hist_names)
        obs = hist_names[iobs]
        if occursin("RAW",obs)
            continue
        end
        if LIST != ""
            if getweight(LIST,string(obs),weight=0) == 0
                continue
            end
        end
        binvalues = get_bin_value_for_observable(grid,obs)
        binedges = get_bin_edges_for_observables(grid,obs)
        binerrs = get_sumw2_for_obervable(grid,obs)
        for i in zip(delindst_sort,delinded_sort)
            binvalues = binvalues[1:end .∉ [collect(i[1]:i[2])] , 1:end]
            binerrs = binerrs[1:end .∉ [collect(i[1]:i[2])] , 1:end]
        end
        nparams = length(pars[1,:])
        nbins = length(binedges)-1
        nparcombos = length(pars[:,1])
        nvalsperpar = Int(length(pars[:,1])/nparams)
        parvals = pars[:,:]
        coeffs = ipol[:,3][iobs][1]
        if cubic
            fitvals = [g_cubic(parvals,coe) for coe in eachrow(coeffs)]
        else
            fitvals = [g(parvals,coe) for coe in eachrow(coeffs)]
        end

        for ibin in 1:nbins
            if sum(binerrs[:,ibin+2] .< 0) > 0
                print("Skipped $obs in bin $ibin for negative sumw2.")
                continue
            end
            binl = binedges[ibin+1]-binedges[ibin]
            for ipar in 1:nparams

                parbeg = 1+nvalsperpar*(ipar-1)
                parend = nvalsperpar*ipar
                parval = parvals[parbeg:parend,:][:,ipar]
                fitval = fitvals[ibin][parbeg:parend]
                orgval = binvalues[parbeg:parend,2+ibin]
                orgvalerr = broadcast(sqrt,binerrs[parbeg:parend,2+ibin] ./ binl) 
                covvar = ipol[iobs,3][3][ibin]

                parvecs = parvals[parbeg:parend,:]
                if cubic
                    fiterrs = [sqrt(g_err_cubic(parvecs[i,:],covvar)) for i in 1:length(parvecs[:,1])]
                else
                    fiterrs = [sqrt(g_err(parvecs[i,:],covvar)) for i in 1:length(parvecs[:,1])]
                end

                if inbound
                    if ipar == 1 && nparams == 6
                        parval=parval[1:end-5]
                        fitval=fitval[1:end-5]
                        fiterrs=fiterrs[1:end-5]
                        orgval=orgval[1:end-5]
                        orgvalerr=orgvalerr[1:end-5]
                    end
                    if ipar == 1 && nparams == 8
                        parval=parval[1:end-3]
                        fitval=fitval[1:end-3]
                        fiterrs=fiterrs[1:end-3]
                        orgval=orgval[1:end-3]
                        orgvalerr=orgvalerr[1:end-3]
                    end
                    if ipar == 2 && nparams == 8
                        parval=parval[1:end-5]
                        fitval=fitval[1:end-5]
                        fiterrs=fiterrs[1:end-5]
                        orgval=orgval[1:end-5]
                        orgvalerr=orgvalerr[1:end-5]
                    end
                end
                
                chisqval = (fitval .- orgval).^2 ./ ((orgvalerr .^ 2 .+ fiterrs .^ 2))
                hchisqpull = fit(Histogram,chisqval)
                if cubic
                    funcchi = Chisq(length(chisqval)-4)
                else
                    funcchi = Chisq(length(chisqval)-3)
                end
                chiprobval = ccdf(funcchi,sum(chisqval))
                push!(chiprob_per_par[ipar],chiprobval)
                push!(chisqvals_per_par[ipar],chisqval)
                
                if skip_neg
                    pullsi = []
                    chisqvali = 0
                    chiprobi = 0
                    chisqvalif = []
                    for i in 1:length(orgvalerr)
                       if orgvalerr[i]^ 2 - fiterrs[i]^ 2 < 0
                            continue
                        end 
                        push!(pullsi,(fitval[i] - orgval[i]) / sqrt(abs(orgvalerr[i]^2 + fiterrs[i]^2)))
                        chisqvali += chisqval[i]
                        push!(chisqvalif,chisqval[i])
                    end
                    if length(chisqvalif)-4 > 0
                        if cubic
                            chiprobi = ccdf(Chisq(length(chisqvalif)-4),sum(chisqvalif))
                        else
                            chiprobi = ccdf(Chisq(length(chisqvalif)-3),sum(chisqvalif))
                        end
                    else
                        chiprobi = 10
                    end
                    push!(pulls,Vector{Float64}(pullsi))
                    push!(chiprob,chiprobi)
                    push!(chisqvals,Vector{Float64}(chisqvalif))
                    push!(iden,string(obs," Bin ",ibin," Par ",parnames[ipar]))

                else
                    push!(chiprob,chiprobval)
                    push!(chisqvals,chisqval)
                    push!(iden,string(obs," Bin ",ibin," Par ",parnames[ipar]))
                    respull = (fitval .- orgval) ./ sqrt.(abs.(orgvalerr .^ 2 .+ fiterrs .^ 2))
                    push!(pulls,respull)
                    push!(pulls_per_par[ipar],respull)
                end
            end
        end
    end
    return pulls,pulls_per_par,chiprob,chiprob_per_par,chisqvals,chisqvals_per_par,iden
end

function plot_all_pulls(ipol,grid::Vector{Observablecontainer};xmax=10,folder::String="fit/fit_grid_pulls",normalize=false,skip_neg=false,cubic=false,inbound=false)
    mkpath(folder)
    res = calc_fit_grid_pulls(ipol,gridsorted,skip_neg=skip_neg,cubic=cubic,inbound=inbound)
    #pulls = res[1]
    respulls = vcat(res[1]...)

    xmax=10
    hrp = fit(Histogram,respulls[abs.(respulls) .< xmax],nbins = 64)
    if normalize
        hrp = normalize(hrp,mode=:pdf)
    end
    fvals = fit(Normal,respulls[abs.(respulls) .< 5])

    plot(hrp,label="Residual Pulls",seriestype=:step)
    x = -xmax:0.1:xmax
    #plot!(x,pdf(Normal(),x)*maximum(hrp.weights)*sqrt(2*pi),label="Normal Distribution")
    plot!(x,pdf(Normal(fvals.μ,fvals.σ),x)*maximum(hrp.weights)*sqrt(2*pi*fvals.σ^2),label=string("Fitted normal: \n",L"$\mu$ = ","$(round(fvals.μ,digits=3))","\n",L"$\sigma$ = ","$(round(fvals.σ,digits=3))"))
    plot!(size=(420,2/3*400),xlabel=L"$\mathrm{p}$",ylabel="Entries") #legend=:topleft

    fname = "grid_pulls"
    if skip_neg
        fname = "grid_pulls_skip_neg"
    end

    savefig(string(folder,"/",fname,".png"))
    savefig(string(folder,"/",fname,".pdf"))

    plot(hrp,label="Residual Pulls",seriestype=:step)
    
    savefig(string(folder,"/",fname,"_pulls.png"))
    savefig(string(folder,"/",fname,"_pulls.pdf"))

    hst = fit(Histogram,Vector{Float64}(res[3]),collect(0:0.05:1.1))
    plt = plot(hst,label="Test sample",seriestype=:step,ylabel="Entries",xlabel="\$p\$-values")#,framestyle=:box)
    #\$p(\\chi^{2}_{\\mathrm{grid}})\$
    copy_ticks(plt,plt[1])
    plot!(size=(450,300))
    savefig(string(folder,"/",fname,"_respvals.png"))
    savefig(string(folder,"/",fname,"_respvals.pdf"))

    binning = collect(0:0.1:10)
    chisq_vals = vcat(res[5]...)
    hst = fit(Histogram,chisq_vals,binning)
    plot(hst,label="Residual Chisq",seriestype=:step)
    savefig(string(folder,"/",fname,"_reschisq.png"))
    savefig(string(folder,"/",fname,"_reschisq.pdf"))

    df = DataFrame(pulls=res[1],chiprob=res[3],chisqvals=res[5],iden=res[7])
    CSV.write(string(folder,"/",fname,"_pulls.csv"),df)
    df = DataFrame(pulls_per_par=res[2],chiprob_per_par=res[4],chisqvals_per_par=res[6])
    CSV.write(string(folder,"/",fname,"_pulls_per_par.csv"),df)
    df = DataFrame(bad_bins=res[7][res[3] .< 0.05])
    CSV.write(string(folder,"/",fname,"_bad_bins.csv"),df)
end

function plot_fit_grid_pulls(ipol,grid::Vector{Observablecontainer}; 
    folder::String="fit_grid", 
    par_mask=[true for i in 1:length(grid[1].hc[1].parname)],
    LIST="",
    cubic=false,
    pdf=false,
    inbound=false,
    pro=0.10,
    parnames=[],
    obsnum=[],
)
    prepath = pwd()
    mkpath(folder)
    hist_names = ipol[:,1]
    pars = get_clean_params_from_mc(grid)
    nperpar = Int(length(pars[:,1])/length(pars[1,:]))

    delindst = collect(1:nperpar:length(pars[:,1]))[findall(x->x==0,par_mask)]
    delinded = collect(1:nperpar:length(pars[:,1]))[findall(x->x==0,par_mask)] .+ (nperpar-1)

    delindst_sort = sort(delindst,rev=true)
    delinded_sort = sort(delinded,rev=true)

    for i in zip(delindst_sort,delinded_sort)
        pars = pars[1:end .∉ [collect(i[1]:i[2])] , 1:end]
    end
    pars = pars[:,par_mask]

    if parnames == []
        parnames = grid[1].hc[1].parname[par_mask]
    end

    def_pars = [pars[end,1],pars[1,2:end]...]
    def_pars_vars = Array{Float64,2}(undef,length(def_pars)*2+3,length(def_pars))

    for i in 1:length(def_pars_vars[:,1])
        par_id = i%length(def_pars)+1
        sign = i > length(def_pars) ? +1 : -1
        def_pars_vars[i,:] = def_pars
        def_pars_vars[i,par_id] = (1+sign*pro)*def_pars[par_id]
    end
    def_pars_vars[end-2,:] = (1+pro) .* def_pars
    def_pars_vars[end-1,:] = (1-pro) .* def_pars
    def_pars_vars[end,:] = (1) .* def_pars

    if length(ipol[1,3]) < 3
        println("Ipol object does not contain covar matrix.")
    end

    for iobs in 1:length(hist_names)
        obs = hist_names[iobs]
        if occursin("RAW",obs)
            continue
        end
        if LIST != ""
            if getweight(LIST,string(obs),weight=0) == 0
                continue
            end
        end
        if obsnum != []
            if obsnum != iobs 
                continue
            end
        end
        binvalues = get_bin_value_for_observable(grid,obs)
        binedges = get_bin_edges_for_observables(grid,obs)
        binerrs = get_sumw2_for_obervable(grid,obs)
        for i in zip(delindst_sort,delinded_sort)
            binvalues = binvalues[1:end .∉ [collect(i[1]:i[2])] , 1:end]
            binerrs = binerrs[1:end .∉ [collect(i[1]:i[2])] , 1:end]
        end
        nparams = length(pars[1,:])
        nbins = length(binedges)-1
        nparcombos = length(pars[:,1])
        nvalsperpar = Int(length(pars[:,1])/nparams)
        parvals = pars[:,:]
        coeffs = ipol[:,3][iobs][1]
        if cubic
            fitvals = [g_cubic(parvals,coe) for coe in eachrow(coeffs)]
        else
            fitvals = [g(parvals,coe) for coe in eachrow(coeffs)]
        end

        if cubic
            def_fitvals_vars = [g_cubic(def_pars_vars,coe) for coe in eachrow(coeffs)]
        else
            def_fitvals_vars = [g(def_pars_vars,coe) for coe in eachrow(coeffs)]
        end

        
        for ibin in 1:nbins
            if sum(binerrs[:,ibin+2] .< 0) > 0
                print("Skipped $obs in bin $ibin for negative sumw2.")
                continue
            end
            binl = binedges[ibin+1]-binedges[ibin]
            for ipar in 1:nparams

                def_fitval_var = def_fitvals_vars[ibin]
                def_fitval_var_ind = [argmin(def_fitval_var),argmax(def_fitval_var)]
                def_parval_low = def_pars_vars[argmin(def_fitval_var),:]
                def_parval_high = def_pars_vars[argmax(def_fitval_var),:]
                
                parbeg = 1+nvalsperpar*(ipar-1)
                parend = nvalsperpar*ipar
                parval = parvals[parbeg:parend,:][:,ipar]
                fitval = fitvals[ibin][parbeg:parend]
                orgval = binvalues[parbeg:parend,2+ibin]
                orgvalerr = broadcast(sqrt,binerrs[parbeg:parend,2+ibin] ./ binl) 
                covvar = ipol[iobs,3][3][ibin]

                def_parvals_low = Matrix{Float64}(undef,length(parval),length(def_parval_low))
                def_parvals_high = Matrix{Float64}(undef,length(parval),length(def_parval_low))

                #def_pars_vars
                def_pars_vars_full = Array{Float64,3}(undef,length(def_pars_vars[:,1]),length(parval),nparams)

                for i in 1:length(def_pars_vars[:,1])
                    for j in 1:length(parval)
                        def_pars_vars_full[i,j,:] = def_pars_vars[i,:]
                        def_pars_vars_full[i,j,ipar] = parval[j]
                    end
                end           
                
                parvecs = parvals[parbeg:parend,:]
                
                if cubic
                    fiterrs = [sqrt(g_err_cubic(parvecs[i,:],covvar)) for i in 1:length(parvecs[:,1])]
                    default_vars_full = [g_cubic(def_pars_vars_full[:,i,:],coeffs[ibin,:]) for i in 1:length(parval)]
                else
                    fiterrs = [sqrt(g_err(parvecs[i,:],covvar)) for i in 1:length(parvecs[:,1])]
                    default_vars_full = [g(def_pars_vars_full[:,i,:],coeffs[ibin,:]) for i in 1:length(parval)]
                end

                dfl = [minimum(default_vars_full[i]) < fitval[i] ? minimum(default_vars_full[i]) : fitval[i] for i in 1:length(parval)]
                dfh = [maximum(default_vars_full[i]) > fitval[i] ? maximum(default_vars_full[i]) : fitval[i] for i in 1:length(parval)]
                #dfh = [default_vars_low[i] > default_vars_high[i] ? default_vars_low[i] : default_vars_high[i] for i in 1:length(default_vars_low)]

                if inbound
                    if ipar == 1 && nparams == 6
                        parval=parval[1:end-5]
                        fitval=fitval[1:end-5]
                        fiterrs=fiterrs[1:end-5]
                        orgval=orgval[1:end-5]
                        orgvalerr=orgvalerr[1:end-5]
                        dfl=dfl[1:end-5]
                        dfh=dfh[1:end-5]
                    end
                    if ipar == 1 && nparams == 8
                        parval=parval[1:end-3]
                        fitval=fitval[1:end-3]
                        fiterrs=fiterrs[1:end-3]
                        orgval=orgval[1:end-3]
                        orgvalerr=orgvalerr[1:end-3]
                        dfl=dfl[1:end-3]
                        dfh=dfh[1:end-3]
                    end
                    if ipar == 2 && nparams == 8
                        parval=parval[1:end-5]
                        fitval=fitval[1:end-5]
                        fiterrs=fiterrs[1:end-5]
                        orgval=orgval[1:end-5]
                        orgvalerr=orgvalerr[1:end-5]
                        dfl=dfl[1:end-5]
                        dfh=dfh[1:end-5]
                    end
                end
                                
                respull = (fitval .- orgval) ./ sqrt.(abs.(orgvalerr .^ 2 .+ fiterrs .^ 2))
                chisqval = (fitval .- orgval).^2 ./ ((orgvalerr .^ 2 .+ fiterrs .^ 2))

                c2 = chisqerr(fitval,orgval,orgvalerr)
                plt = plot(parval,orgval,yerr=orgvalerr,color="grey",label="MC")
                plot!(parval,fitval,label="Fit")#label="Fit \n chi2 = $c2")
                plot!(parval,fitval .- fiterrs,fill_between=fitval .+ fiterrs,linecolor=:red,color=:red,alpha=0.3,label="Fit Err")
                plot!(parval,dfl,fill_between=dfh,linecolor=:blue,color=:blue,alpha=0.3,label="Default Var")
                vline!([def_pars[ipar]],label="Default",color="grey")

                par_name = parnames[ipar]
                if occursin("/",par_name)
                    par_name = split(par_name,"/")[end]
                end
                if LIST != "" && pdf == false
                    obstitle = getobsname(LIST,obs)
                    title!("$(obstitle), Bin: $(ibin)")
                else
                    #title!("$(obs), Bin: $(ibin)")
                end
                #ylabel!("Entries")
                ylabel!(L"\frac{1}{\sigma}\frac{d\sigma}{dS}")
                #ylabel!(L"1/\sigma \cdot d \sigma/dS")
                xlabel!("$par_name")

                cd(folder)
                if(!(obs in readdir()))
                    mkpath(obs[2:end])
                    mkpath(string(obs[2:end],"/CSV/"))
                end
                cd(prepath)

                df = DataFrame(xval=parval, yval_orgval=orgval, yval_orgvalerr=orgvalerr, yval_fitval=fitval, yval_fitvalerr=fiterrs
                )
                CSV.write(string(folder,obs,"/CSV/",split(obs,"/")[end-1],"-$par_name-$(ibin)-2.csv"),df)

                copy_ticks(plt,plt[1])
                plot!(size=(510,405),top_margin=6Plots.mm)
                if pdf 
                    plot!(size=(500,260),top_margin=0Plots.mm)
                    savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$par_name-$(ibin)-2.pdf"))
                end
                savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$par_name-$(ibin)-2.png"))

                scatter(parval,respull,yerr=1,color="black",label="Pull",ms=3,msc="black")
                xlabel!("$par_name")
                ylabel!("Res Pull")
                title!("$(obs), Bin: $(ibin)")
                if pdf
                    savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$par_name-$(ibin)-Pull.pdf"))
                end
                savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$par_name-$(ibin)-Pull.png"))
                
                hchisqpull = fit(Histogram,chisqval)
                if cubic
                    funcchi = Chisq(length(chisqval)-4)
                else
                    funcchi = Chisq(length(chisqval)-3)
                end
                chiprob = ccdf(funcchi,sum(chisqval))

                plot(hchisqpull,label="chisq ,n=$(length(chisqval)) \n chiprob = $(chiprob)",seriestype=:step)
                xlabel!("$par_name")
                ylabel!("Chisq Pull")
                title!("$(obs), Bin: $(ibin)")
                if pdf
                    savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$par_name-$(ibin)-chisq.pdf"))
                end
                savefig(string(folder,obs,"/",split(obs,"/")[end-1],"-$par_name-$(ibin)-chisq.png"))
            end
        end
    end
end

function plot_sensitivities(ipol,gridsorted::Vector{Observablecontainer}; folder::String="sens", par_mask=[true for i in 1:length(gridsorted[1].hc[1].parname)],LIST="",cubic=false,pdf=false,inbound=false,pro=0.10,parnames=[])

    prepath = pwd()
    mkpath(folder)
    hist_names = ipol[:,1]
    pars = get_clean_params_from_mc(gridsorted)
    nperpar = Int(length(pars[:,1])/length(pars[1,:]))

    delindst = collect(1:nperpar:length(pars[:,1]))[findall(x->x==0,par_mask)]
    delinded = collect(1:nperpar:length(pars[:,1]))[findall(x->x==0,par_mask)] .+ (nperpar-1)

    delindst_sort = sort(delindst,rev=true)
    delinded_sort = sort(delinded,rev=true)

    for i in zip(delindst_sort,delinded_sort)
        pars = pars[1:end .∉ [collect(i[1]:i[2])] , 1:end]
    end
    pars = pars[:,par_mask]
    if parnames == []
        parnames = gridsorted[1].hc[1].parname[par_mask]
    end

    def_pars = [pars[end,1],pars[1,2:end]...]
    def_pars_vars = Array{Float64,2}(undef,length(def_pars)*2+3,length(def_pars))

    for i in 1:length(def_pars_vars[:,1])
        par_id = i%length(def_pars)+1
        sign = i > length(def_pars) ? +1 : -1
        def_pars_vars[i,:] = def_pars
        def_pars_vars[i,par_id] = (1+sign*pro)*def_pars[par_id]
    end
    def_pars_vars[end-2,:] = (1+pro) .* def_pars
    def_pars_vars[end-1,:] = (1-pro) .* def_pars
    def_pars_vars[end,:] = (1) .* def_pars

    full_sens = []
    full_sens_ipol = []

    for iobs in 1:length(hist_names)

        obs = hist_names[iobs]
        if occursin("RAW",obs)
            continue
        end
        if LIST != ""
            if getweight(LIST,string(obs),weight=0) == 0
                continue
            end
        end

        binvalues = get_bin_value_for_observable(gridsorted,obs)
        binedges = get_bin_edges_for_observables(gridsorted,obs)
        binerrs = get_sumw2_for_obervable(gridsorted,obs)

        for i in zip(delindst_sort,delinded_sort)
            binvalues = binvalues[1:end .∉ [collect(i[1]:i[2])] , 1:end]
            binerrs = binerrs[1:end .∉ [collect(i[1]:i[2])] , 1:end]
        end

        nparams = length(pars[1,:])
        nbins = length(binedges)-1
        nparcombos = length(pars[:,1])
        nvalsperpar = Int(length(pars[:,1])/nparams)
        parvals = pars[:,:]
        coeffs = ipol[:,3][iobs][1]

        if cubic
            fitvals = [g_cubic(parvals,coe) for coe in eachrow(coeffs)]
        else
            fitvals = [g(parvals,coe) for coe in eachrow(coeffs)]
        end

        sens = Matrix{Float64}(undef,(nparams,nbins))
        sensipol = Matrix{Float64}(undef,(nparams,nbins))

        for ibin in 1:nbins
            if sum(binerrs[:,ibin+2] .< 0) > 0
                print("Skipped $obs in bin $ibin for negative sumw2.")
                continue
            end
            # if  binvalues 
                
            # end
            binl = binedges[ibin+1]-binedges[ibin]
            f = x->g_cubic(x,coeffs[ibin,:])
            grads = ForwardDiff.gradient(f,def_pars)
            
            for ipar in 1:nparams
                
                parbeg = 1+nvalsperpar*(ipar-1)
                parend = nvalsperpar*ipar
                parval = parvals[parbeg:parend,:][:,ipar]
                orgval = binvalues[parbeg:parend,2+ibin]

                parvalmid = findmin(abs.(parval.-def_pars[ipar]))[2]
                while parval[parvalmid] == 0.0
                    parvalmid += 1
                end
                #fitval = fitvals[ibin][parbeg:parend]

                delpar = parval[parvalmid+1]-parval[parvalmid]
                delmc = orgval[parvalmid+1]-orgval[parvalmid]
                #s1 = (delmc/orgval[parvalmid])/(delpar/parval[parvalmid]) 
                #s2 = (orgval[parvalmid+2]-orgval[parvalmid]/orgval[parvalmid])/(parval[parvalmid+2]-parval[parvalmid]/parval[parvalmid]) 
                s1 = (log(abs(orgval[parvalmid]))-log(abs(orgval[parvalmid+1]))) / (log(abs(parval[parvalmid]))-log(abs(parval[parvalmid+1])))
                s2 = (log(abs(orgval[parvalmid]))-log(abs(orgval[parvalmid+2]))) / (log(abs(parval[parvalmid]))-log(abs(parval[parvalmid+2])))

                #if occursin("Constituent",parnames[ipar])
                #    sens[ipar,ibin] = s1
                #else
                    sens[ipar,ibin] = (s1+s2)/2#(delmc/orgval[parvalmid])/(delpar/parval[parvalmid]) 
                    sensipol[ipar,ibin] = grads[ipar] * (parval[parvalmid]/orgval[parvalmid])
                #end

            end
        end
        binmids = binedges[1:end-1] .+ (binedges[2:end]-binedges[1:end-1]) ./ 2
        
        plot()
        if nbins > 1
            [plot!(binmids,sens[i,:],label=split(parnames[i],"/")[end],markershape=:auto,markersize=3,lw=1,legend = :outertopleft,size=(1000,500),margin=4Plots.mm) for i in 1:nparams]
        else
            [plot!([(binedges[1]+binedges[2])/2],sens[i,:],label=split(parnames[i],"/")[end],markershape=:auto,lw=0,legend = :outertopleft,margin=4Plots.mm) for i in 1:nparams]
        end
        if LIST != ""
            obstitle = getobsname(LIST,obs)
            #title!("$(obstitle), Bin: $(ibin)")
            xlabel!("$obstitle")
            else
            xlabel!("$obstitle")
            #title!("$(obs), Bin: $(ibin)")
        end
        plot!(size=(600,300),margin=1Plots.mm,leftmargin=-14Plots.mm,bottommargin=2Plots.mm)
        ylabel!("\n Sensitivity")

        cd(folder)
        if(!(obs in readdir()))
            mkpath(obs[2:end])
            #mkpath(string(obs[2:end],"/CSV/"))
        end
        cd(prepath)

        savefig(string(folder,obs,"/",split(obs,"/")[end-1],"_sensitivities.png"))
        savefig(string(folder,obs,"/",split(obs,"/")[end-1],"_sensitivities.pdf"))

        plot()
        if nbins > 1
            [plot!(binmids,sensipol[i,:],label=split(parnames[i],"/")[end],markershape=:auto,markersize=3,lw=1,legend = :outertopleft,size=(1000,500),margin=4Plots.mm) for i in 1:nparams]
        else
            [plot!([(binedges[1]+binedges[2])/2],sensipol[i,:],label=split(parnames[i],"/")[end],markershape=:auto,lw=0,legend = :outertopleft,margin=4Plots.mm) for i in 1:nparams]
        end
        plot!(size=(600,300),margin=1Plots.mm,leftmargin=-14Plots.mm,bottommargin=2Plots.mm)
        ylabel!("\n Sensitivity")
        if LIST != ""
            obstitle = getobsname(LIST,obs)
            #title!("$(obstitle), Bin: $(ibin)")
            xlabel!("$obstitle")
            else
            xlabel!("$obstitle")
            #title!("$(obs), Bin: $(ibin)")
        end

        savefig(string(folder,obs,"/ipol",split(obs,"/")[end-1],"_sensitivities.png"))
        savefig(string(folder,obs,"/ipol",split(obs,"/")[end-1],"_sensitivities.pdf"))

        push!(full_sens,(name=obstitle,sens=sens,code=obs))
        push!(full_sens_ipol,(name=obstitle,sens=sensipol,code=obs))
    end
    [full_sens,full_sens_ipol]
end

function print_sens_table(fullsens;fname="",header=[],baddims=[],iround=true,imean=true)
    meansens=[]
    if imean
        meansens = permutedims(reshape(vcat(vcat([[mean.([a for a in eachrow(fullsens[iobs].sens)]).*100] for iobs in 1:length(fullsens)]...)...),(length(fullsens[1].sens[:,1]),length(fullsens))))
    else
        meansens = permutedims(reshape(vcat(vcat([[var.([a for a in eachrow(fullsens[iobs].sens)]).*100] for iobs in 1:length(fullsens)]...)...),(length(fullsens[1].sens[:,1]),length(fullsens))))
    end
    if iround
        meansens = abs.(round.(Int,meansens[1:end .∉ [baddims],1:end]))
    else
        meansens = abs.(round.(meansens[1:end .∉ [baddims],1:end]))
    end
    namesens = [fullsens[iobs].name for iobs in 1:length(fullsens)][1:end .∉ [baddims],1:end]
    tablesens = hcat(namesens,meansens)
    pretty_table(tablesens,header=header)

    f=open(string(fname),"w+")
    pretty_table(f,tablesens,header=header)
    close(f)

    namesens = [fullsens[iobs].code for iobs in 1:length(fullsens)][1:end .∉ [baddims],1:end]
    tablesens = hcat(namesens,meansens)
    pretty_table(tablesens,header=header)
    fname = string(fname,".tex")
    f=open(string(fname),"w+")
    pretty_table(f,tablesens,header=header,backend=Val(:latex))
    close(f)
end