######################################
####Plotting for the interpolation####
######################################
function plot_chi2_int_per_obs(ipol;folder::String="fit_chisq")
    if(!(folder in readdir()))
        mkpath(folder)
    end
    chi_red = []
    for i in 1:length(ipol[:,3])
        chi_red = ipol[:,3][i][2][:,1]
        name = ipol[:,1][i]
        pathname_obs=string(folder,"/")

        for i in split(name,"/")[2:end-1]
            pathname_obs = string(pathname_obs,string(i),"/")
        end
        if(!(pathname_obs in readdir()))
            mkpath(pathname_obs)
        end

        name=split(name,"/")[end]
        if length(chi_red[chi_red .!= Inf]) != 0
            histogram(chi_red[chi_red .!= Inf],bins=16)
            savefig(string(pathname_obs,name,"_chisq.pdf"))
            savefig(string(pathname_obs,name,"_chisq.png"))
        end
    end
end

function plot_chi2_int(ipol;folder::String="fit_chisq")
    if(!(folder in readdir()))
        mkpath(folder)
    end
    chi_red_all=[]
    chi_red = []
    for i in 1:Int(length(ipol[:,3]))
        if occursin("RAW",ipol[i,1])
            continue
        end
        chi_red = ipol[:,3][i][2][:,1]
        name = ipol[:,1][i]
        if length(chi_red[(chi_red .!= Inf) .& (chi_red .> 0)]) != 0
            push!(chi_red_all,chi_red[(chi_red .!= Inf) .& (chi_red .> 0)])
        end
    end
    chi_red_all = vcat(chi_red_all...)
    hst = fit(Histogram,chi_red_all,nbins=64)
    plot(hst,label="n entries = $(length(chi_red_all))",xlabel="\$\\chi_{\\mathrm{red}}\$",ylabel="Entries")
    df = DataFrame(xmin=collect(hst.edges[1])[1:end-1],xmax=collect(hst.edges[1])[2:end],yval=hst.weights)
    mkpath(string("./",folder,"/CSV/"))
    CSV.write(string("./",folder,"/CSV/all.csv"),df)
    savefig(string("./",folder,"/all.pdf"))
    savefig(string("./",folder,"/all.png"))

    cut_vals=[100,70,50,10,5,1,0.1,0.01]
    for cut_vals_i in cut_vals
        numname = cut_vals_i%1 == 0 ? Int(cut_vals_i) : replace(string(cut_vals_i),"."=>"")
        hst = fit(Histogram,chi_red_all[chi_red_all .< cut_vals_i],nbins=128)
        df = DataFrame(xmin=collect(hst.edges[1])[1:end-1],mxax=collect(hst.edges[1])[2:end],yval=hst.weights)
        CSV.write(string("./",folder,"/CSV/all_sm$(numname).csv"),df)
        plot(hst,label="n entries = $(length(chi_red_all[chi_red_all .< cut_vals_i]))",xlabel=xlabel="\$\\chi_{\\mathrm{red}}\$",ylabel="Entries")
        savefig(string("./",folder,"/all_sm$(numname).pdf"))
        savefig(string("./",folder,"/all_sm$(numname).png"))
    end
end

function plot_chi2_prob_int(ipol;folder::String="fit_chisq")
    if(!(folder in readdir()))
        mkpath(folder)
    end
    chi_prob_all=[]
    chi_prob = []
    for i in 1:Int(length(ipol[:,3]))
        if occursin("RAW",ipol[i,1])
            continue
        end
        chi_prob = [ccdf(Chisq(ipol[i,3][2][j,6]-1),ipol[i,3][2][j,5]) for j in 1:length(ipol[i,3][2][:,5])]
        if length(chi_prob[(chi_prob .!= Inf) ]) != 0
            push!(chi_prob_all,chi_prob[(chi_prob .!= Inf) ])
        end
    end
    chi_prob_all = vcat(chi_prob_all...)

    hst = fit(Histogram,chi_prob_all,nbins=64)

    plot(hst,label="n entries = $(length(chi_prob_all))",xlabel="p",ylabel="Entries")
    df = DataFrame(xmin=collect(hst.edges[1])[1:end-1],xmax=collect(hst.edges[1])[2:end],yval=hst.weights)
    mkpath(string("./",folder,"/CSV/"))
    CSV.write(string("./",folder,"/CSV/all_prob.csv"),df)
    savefig(string("./",folder,"/all_prob.pdf"))
    savefig(string("./",folder,"/all_prob.png"))
    
    cut_vals=[10,1,0.1,0.01]
    for cut_vals_i in cut_vals
        numname = cut_vals_i%1 == 0 ? Int(cut_vals_i) : replace(string(cut_vals_i),"."=>"")
        hst = fit(Histogram,chi_prob_all[chi_prob_all .< cut_vals_i],nbins=128)
        df = DataFrame(xmin=collect(hst.edges[1])[1:end-1],mxax=collect(hst.edges[1])[2:end],yval=hst.weights)
        CSV.write(string("./",folder,"/CSV/all_prob_sm$(numname).csv"),df)
        plot(hst,label="n entries = $(length(chi_prob_all[chi_prob_all .< cut_vals_i]))",xlabel="p",ylabel="Entries")
        savefig(string("./",folder,"/all_prob_sm$(numname).pdf"))
        savefig(string("./",folder,"/all_prob_sm$(numname).png"))
    end
end

function plot_chi2_int_all(ipol;folder::String="fit_chisq_all")
    if(!(folder in readdir()))
        mkpath(folder)
    end
    chi_red_all=[]
    chi_red = []
    for i in 1:Int(length(ipol[:,3]))
        chi_red = ipol[:,3][i][2][:,1]
        name = ipol[:,1][i]
        if length(chi_red[chi_red .!= Inf]) != 0
            push!(chi_red_all,chi_red[chi_red .!= Inf])
        end
    end
    chi_red_all = vcat(chi_red_all...)
    histogram(chi_red_all,bins=64,label="n entries = $(length(chi_red_all))")
    savefig(string("./",folder,"/all.pdf"))
    savefig(string("./",folder,"/all.png"))
    histogram(chi_red_all[chi_red_all .< 50],bins=128,label="n entries = $(length(chi_red_all[chi_red_all .< 50]))")
    savefig(string("./",folder,"/all_sm50.pdf"))
    savefig(string("./",folder,"/all_sm50.png"))
    histogram(chi_red_all[chi_red_all .< 10],bins=128,label="n entries = $(length(chi_red_all[chi_red_all .< 10]))")
    savefig(string("./",folder,"/all_sm10.pdf"))
    savefig(string("./",folder,"/all_sm10.png"))
    histogram(chi_red_all[chi_red_all .< 0.1],bins=128,label="n entries = $(length(chi_red_all[chi_red_all .< 0.1]))")
    savefig(string("./",folder,"/all_sm01.pdf"))
    savefig(string("./",folder,"/all_sm01.png"))
end

function plot_chisqprob(ipol;name="",folder="fit",dof=0)
    PREFIX = folder
    mkpath(folder)
    nameobs = ipol[:,1]
    binedges = ipol[:,2]

    big_chisq_list = []
    [push!(big_chisq_list,ipol[ihist,3][2]) for ihist in 1:length(nameobs) if !occursin("RAW",nameobs[ihist])]

    fulllist = vcat(big_chisq_list...)

    funcchi = Chisq(Int(floor(fulllist[1,3])))

    if dof != 0
        funcchi = Chisq(Int(dof))
        name=string(name,"_dof$(dof)")
    end

    chiprob = ccdf.(funcchi,fulllist[:,2])
    
    onered = []
    zerored = []
    onered_inverse = []
    onered_inverse_nozero = []

    for i in 1:length(fulllist[:,1])
        if isapprox(fulllist[i,1],1,atol=0.2)
            push!(onered,chiprob[i])
        else
            push!(onered_inverse,chiprob[i])
            if isapprox(fulllist[i,1],0,atol=0.2)
                push!(zerored,chiprob[i])
            else
                push!(onered_inverse_nozero,chiprob[i])
            end
        end
    end

    onered = Vector{Float64}(onered)
    zerored = Vector{Float64}(zerored)
    onered_inverse = Vector{Float64}(onered_inverse)
    onered_inverse_nozero = Vector{Float64}(onered_inverse_nozero)

    bins = collect(0:0.03:1.02)

    plot(fit(Histogram,chiprob,bins),label="chiprob,n=$(length(chiprob))",seriestype=:step)
    plot!(fit(Histogram,onered,bins),label="onered,n=$(length(onered))",seriestype=:step)
    plot!(fit(Histogram,onered[onered .> 0.01],bins),label="onered > 0.01,n=$(length(onered[onered .> 0.01]))",seriestype=:step)
    plot!(fit(Histogram,onered[onered .> 0.05],bins),label="onered > 0.05,n=$(length(onered[onered .> 0.05]))",seriestype=:step)
    plot!(fit(Histogram,zerored,bins),label="zerored,n=$(length(zerored))",seriestype=:step)
    plot!(fit(Histogram,onered_inverse,bins),label="onered_inverse,n=$(length(onered_inverse))",seriestype=:step)
    plot!(fit(Histogram,onered_inverse_nozero,bins),label="onered_inverse_nozero,n=$(length(onered_inverse_nozero))",seriestype=:step)

    savefig(string(PREFIX,"/",name,"_full.png"))
    savefig(string(PREFIX,"/",name,"_full.pdf"))

    plot(fit(Histogram,onered,bins),label="onered,n=$(length(onered))",seriestype=:step)
    plot!(fit(Histogram,onered[onered .> 0.01],bins),label="onered > 0.01,n=$(length(onered[onered .> 0.01]))",seriestype=:step)
    plot!(fit(Histogram,onered[onered .> 0.05],bins),label="onered > 0.05,n=$(length(onered[onered .> 0.05]))",seriestype=:step)
    plot!(fit(Histogram,zerored,bins),label="zerored,n=$(length(zerored))",seriestype=:step)

    savefig(string(PREFIX,"/",name,"_full_better.png"))
    savefig(string(PREFIX,"/",name,"_full_better.pdf"))

    plot(fit(Histogram,chiprob,bins),label="chiprob,n=$(length(chiprob))",seriestype=:step)
    savefig(string(PREFIX,"/",name,"_all.png"))
    savefig(string(PREFIX,"/",name,"_all.pdf"))

    plot(fit(Histogram,chiprob[chiprob .> 0.01],bins),label="chiprob,n=$(length(chiprob[chiprob .> 0.01]))",seriestype=:step)
    savefig(string(PREFIX,"/",name,"_all_cut.png"))
    savefig(string(PREFIX,"/",name,"_all_cut.pdf"))
    
    plot(fit(Histogram,onered,bins),label="onered,n=$(length(onered))",seriestype=:step)
    plot!(fit(Histogram,onered[onered .> 0.01],bins),label="onered > 0.01,n=$(length(onered[onered .> 0.01]))",seriestype=:step)
    plot!(fit(Histogram,onered[onered .> 0.05],bins),label="onered > 0.05,n=$(length(onered[onered .> 0.05]))",seriestype=:step)
    savefig(string(PREFIX,"/",name,"_onered.png"))
    savefig(string(PREFIX,"/",name,"_onered.pdf"))

    plot(fit(Histogram,zerored,bins),label="zerored,n=$(length(zerored))",seriestype=:step)
    savefig(string(PREFIX,"/",name,"_zerored.png"))
    savefig(string(PREFIX,"/",name,"_zerored.pdf"))

    plot(fit(Histogram,onered_inverse,bins),label="onered_inverse,n=$(length(onered_inverse))",seriestype=:step)
    plot!(fit(Histogram,onered_inverse_nozero,bins),label="onered_inverse_nozero,n=$(length(onered_inverse_nozero))",seriestype=:step)
    savefig(string(PREFIX,"/",name,"_onred_inv.png"))
    savefig(string(PREFIX,"/",name,"_onred_inv.pdf"))

end

function plot_chisqprob_constr(ipol;name="",folder="fit",LIST="",dof=0)
    PREFIX = folder
    name = string(name,"_constr")
    mkpath(folder)
    mkpath(string(folder,"/CSV/"))
    nameobs = ipol[:,1]
    binedges = ipol[:,2]
    big_chisq_list = []
    [push!(big_chisq_list,ipol[ihist,3][2]) for ihist in 1:length(nameobs) if (!occursin("RAW",nameobs[ihist]) && getweight(LIST,nameobs[ihist],weight=0) != 0)]

    fulllist = vcat(big_chisq_list...)
    #funcchi = Chisq(Int(floor(fulllist[1,3])))
    funcchi = Chisq(Int(floor(ipol[1,3][2][1,6])))
    if dof != 0
        funcchi = Chisq(Int(dof))
        name=string(name,"_dof$(dof)")
    end
    chiprob = ccdf.(funcchi,fulllist[:,2])
    
    onered = []
    zerored = []
    onered_inverse = []
    onered_inverse_nozero = []

    for i in 1:length(fulllist[:,1])
        if isapprox(fulllist[i,1],1,atol=0.2)
            push!(onered,chiprob[i])
        else
            push!(onered_inverse,chiprob[i])
            if isapprox(fulllist[i,1],0,atol=0.2)
                push!(zerored,chiprob[i])
            else
                push!(onered_inverse_nozero,chiprob[i])
            end
        end
    end

    onered = Vector{Float64}(onered)
    zerored = Vector{Float64}(zerored)
    onered_inverse = Vector{Float64}(onered_inverse)
    onered_inverse_nozero = Vector{Float64}(onered_inverse_nozero)

    bins = collect(0:0.08:1.02)

    hst_chiprob = fit(Histogram,chiprob,bins)
    hst_onered = fit(Histogram,onered,bins)
    hst_onered001 = fit(Histogram,onered[onered .> 0.01],bins)
    hst_onered005 = fit(Histogram,onered[onered .> 0.05],bins)
    hst_zerored = fit(Histogram,zerored,bins)
    hst_onered_inverse = fit(Histogram,onered_inverse,bins)
    hst_onered_inverse_nozero = fit(Histogram,onered_inverse_nozero,bins)

    df = DataFrame(xmin=collect(hst_chiprob.edges[1])[1:end-1],xmax=collect(hst_chiprob.edges[1])[2:end], yval_chiprob = hst_chiprob.weights, yval_onered = hst_onered.weights, yval_onered001 = hst_onered001.weights, yval_onered005 = hst_onered005.weights, yval_zerored = hst_zerored.weights, yval_onered_inverse = hst_onered_inverse.weights, yval_onered_inverse_nozero = hst_onered_inverse_nozero.weights
    )
    CSV.write(string(PREFIX,"/CSV/",name,"_full.csv"),df)

    plot(hst_chiprob,label="chiprob,n=$(length(chiprob))",seriestype=:step)
    plot!(hst_onered,label="onered,n=$(length(onered))",seriestype=:step)
    plot!(hst_onered001,label="onered > 0.01,n=$(length(onered[onered .> 0.01]))",seriestype=:step)
    plot!(hst_onered005,label="onered > 0.05,n=$(length(onered[onered .> 0.05]))",seriestype=:step)
    plot!(hst_zerored,label="zerored,n=$(length(zerored))",seriestype=:step)
    plot!(hst_onered_inverse,label="onered_inverse,n=$(length(onered_inverse))",seriestype=:step)
    plot!(hst_onered_inverse_nozero,label="onered_inverse_nozero,n=$(length(onered_inverse_nozero))",seriestype=:step)
    
    savefig(string(PREFIX,"/",name,"_full.png"))
    savefig(string(PREFIX,"/",name,"_full.pdf"))

    
    df = DataFrame(xmin=collect(hst_chiprob.edges[1])[1:end-1],xmax=collect(hst_chiprob.edges[1])[2:end], yval_onered = hst_onered.weights, yval_onered001 = hst_onered001.weights, yval_onered005 = hst_onered005.weights, yval_zerored = hst_zerored.weights)
    CSV.write(string(PREFIX,"/CSV/",name,"_full_better.csv"),df)
    
    plot(hst_onered,label="onered,n=$(length(onered))",seriestype=:step)
    plot!(hst_onered001,label="onered > 0.01,n=$(length(onered[onered .> 0.01]))",seriestype=:step)
    plot!(hst_onered005,label="onered > 0.05,n=$(length(onered[onered .> 0.05]))",seriestype=:step)
    plot!(hst_zerored,label="zerored,n=$(length(zerored))",seriestype=:step)

    savefig(string(PREFIX,"/",name,"_full_better.png"))
    savefig(string(PREFIX,"/",name,"_full_better.pdf"))

    bins = collect(0:0.03:1.02)

    hst_chiprob = fit(Histogram,chiprob,bins)

    df = DataFrame(xmin=collect(hst_chiprob.edges[1])[1:end-1],xmax=collect(hst_chiprob.edges[1])[2:end], yval = hst_chiprob.weights)
    CSV.write(string(PREFIX,"/CSV/",name,"_all.csv"),df)

    plot(fit(Histogram,chiprob,bins),label="chiprob,n=$(length(chiprob))",seriestype=:step)
    savefig(string(PREFIX,"/",name,"_all.png"))
    savefig(string(PREFIX,"/",name,"_all.pdf"))

    bins = collect(0:0.08:1.02)

    hst_chiprob = fit(Histogram,chiprob[chiprob .> 0.01],bins)

    df = DataFrame(xmin=collect(hst_chiprob.edges[1])[1:end-1],xmax=collect(hst_chiprob.edges[1])[2:end], yval = hst_chiprob.weights)
    CSV.write(string(PREFIX,"/CSV/",name,"_all_cut.csv"),df)

    plot(hst_chiprob,label="chiprob,n=$(length(chiprob[chiprob .> 0.01]))",seriestype=:step)
    savefig(string(PREFIX,"/",name,"_all_cut.png"))
    savefig(string(PREFIX,"/",name,"_all_cut.pdf"))

    hst_onered = fit(Histogram,onered,bins)
    hst_onered001 = fit(Histogram,onered[onered .> 0.01],bins)
    hst_onered005 = fit(Histogram,onered[onered .> 0.05],bins)
    hst_zerored = fit(Histogram,zerored,bins)
    hst_onered_inverse = fit(Histogram,onered_inverse,bins)
    hst_onered_inverse_nozero = fit(Histogram,onered_inverse_nozero,bins)

    df = DataFrame(xmin=collect(hst_chiprob.edges[1])[1:end-1],xmax=collect(hst_chiprob.edges[1])[2:end], yval_onered = hst_onered.weights, yval_onered001 = hst_onered001.weights, yval_onered005 = hst_onered005.weights)
    CSV.write(string(PREFIX,"/CSV/",name,"_onered.csv"),df)
    
    plot(hst_onered,label="onered,n=$(length(onered))",seriestype=:step)
    plot!(hst_onered001,label="onered > 0.01,n=$(length(onered[onered .> 0.01]))",seriestype=:step)
    plot!(hst_onered005,label="onered > 0.05,n=$(length(onered[onered .> 0.05]))",seriestype=:step)
    savefig(string(PREFIX,"/",name,"_onered.png"))
    savefig(string(PREFIX,"/",name,"_onered.pdf"))

    df = DataFrame(xmin=collect(hst_chiprob.edges[1])[1:end-1],xmax=collect(hst_chiprob.edges[1])[2:end], yval = hst_zerored.weights)
    CSV.write(string(PREFIX,"/CSV/",name,"_zerored.csv"),df)

    plot(hst_zerored,label="zerored,n=$(length(zerored))",seriestype=:step)
    savefig(string(PREFIX,"/",name,"_zerored.png"))
    savefig(string(PREFIX,"/",name,"_zerored.pdf"))

    df = DataFrame(xmin=collect(hst_chiprob.edges[1])[1:end-1],xmax=collect(hst_chiprob.edges[1])[2:end], yval_onered_inverse = hst_onered_inverse.weights, yval_onered_inverse_nozero = hst_onered_inverse_nozero.weights)
    CSV.write(string(PREFIX,"/CSV/",name,"_onred_inv.csv"),df)

    plot(hst_onered_inverse,label="onered_inverse,n=$(length(onered_inverse))",seriestype=:step)
    plot!(hst_onered_inverse_nozero,label="onered_inverse_nozero,n=$(length(onered_inverse_nozero))",seriestype=:step)
    savefig(string(PREFIX,"/",name,"_onred_inv.png"))
    savefig(string(PREFIX,"/",name,"_onred_inv.pdf"))
end

function plot_chired_for_chisqprob(ipol;name="",folder="fit",xmax=15,dof=0)
    PREFIX = folder
    mkpath(folder)
    nameobs = ipol[:,1]
    binedges = ipol[:,2]

    big_chisq_list = []
    [push!(big_chisq_list,ipol[ihist,3][2]) for ihist in 1:length(nameobs) if !occursin("RAW",nameobs[ihist])]

    fulllist = vcat(big_chisq_list...)

    funcchi = Chisq(Int(floor(fulllist[1,3])))
    if dof != 0
        funcchi = Chisq(Int(dof))
        name=string(name,"_dof$(dof)")
    end
    chiprob = ccdf.(funcchi,fulllist[:,2])
    
    oneprob = []
    zeroprob01 = []
    zeroprob05 = []
    goodrest = []

    for i in 1:length(fulllist[:,1])
        if (isnan(fulllist[i,1]) || fulllist[i,1]^2 == Inf) 
            continue
        end
        if chiprob[i] <= 0.01
            push!(zeroprob01,fulllist[i,1])
        elseif chiprob[i] <= 0.05
            push!(zeroprob05,fulllist[i,1])
        elseif chiprob[i] < 0.9
            push!(goodrest,fulllist[i,1])
        else
            push!(oneprob,fulllist[i,1])
        end
    end
    oneprob = Vector{Float64}(oneprob)
    zeroprob01 = Vector{Float64}(zeroprob01)
    zeroprob05 = Vector{Float64}(zeroprob05)
    goodrest = Vector{Float64}(goodrest)

    binedg = collect(0:xmax/32:xmax)
    honeprob = fit(Histogram,vcat(oneprob,goodrest,zeroprob05,zeroprob01),binedg)
    binedg = collect(honeprob.edges[1])

    hgoodrest = fit(Histogram,vcat(goodrest,zeroprob05,zeroprob01),binedg)
    hzeroprob05 = fit(Histogram,vcat(zeroprob05,zeroprob01),binedg)
    hzeroprob01 = fit(Histogram,vcat(zeroprob01),binedg)

    plot(honeprob,label="p > 0.9,n=$(length(oneprob))",xlabel="\$\\chi^{2}_{\\text{red}}\$")
    plot!(hgoodrest,label="0.9 > p > 0.05,n=$(length(goodrest))")
    plot!(hzeroprob05,label="0.05 > p > 0.01,n=$(length(zeroprob05))")
    plot!(hzeroprob01,label="0.01 > p,n=$(length(zeroprob01))")

    savefig(string(PREFIX,"/",name,"_redchi_full.png"))
    savefig(string(PREFIX,"/",name,"_redchi_full.pdf"))

    honeprob = fit(Histogram,vcat(oneprob,goodrest,zeroprob05),binedg)
    hgoodrest = fit(Histogram,vcat(goodrest,zeroprob05),binedg)
    hzeroprob05 = fit(Histogram,vcat(zeroprob05),binedg)
    plot(honeprob,label="p > 0.9,n=$(length(oneprob))",xlabel="\$\\chi^{2}_{\\text{red}}\$")
    plot!(hgoodrest,label="0.9 > p > 0.05,n=$(length(goodrest))")
    plot!(hzeroprob05,label="0.05 > p > 0.01,n=$(length(zeroprob05))")
    
    savefig(string(PREFIX,"/",name,"_redchi_nozero.png"))
    savefig(string(PREFIX,"/",name,"_redchi_nozero.pdf"))

    honeprob = fit(Histogram,vcat(oneprob,goodrest,zeroprob05),binedg)
    hgoodrest = fit(Histogram,vcat(goodrest,zeroprob05),binedg)
    plot(honeprob,label="p > 0.9,n=$(length(oneprob))",xlabel="\$\\chi^{2}_{\\text{red}}\$")
    plot!(hgoodrest,label="0.9 > p > 0.05,n=$(length(goodrest))")
    
    savefig(string(PREFIX,"/",name,"_redchi_nozero_05.png"))
    savefig(string(PREFIX,"/",name,"_redchi_nozero_05.pdf"))
end

function plot_chired_for_chisqprob_constr(ipol;
    name="",
    folder="fit",
    LIST="",
    xmax=15,
    dof=0,
    recalc=false,
    sorted = "",
)
    PREFIX = folder
    name = string(name,"_constr")
    mkpath(folder)
    mkpath(string(folder,"/CSV/"))
    nameobs = ipol[:,1]
    binedges = ipol[:,2]
    big_chisq_list = []
    [push!(big_chisq_list,ipol[ihist,3][2]) for ihist in 1:length(nameobs) if (!occursin("RAW",nameobs[ihist]) && getweight(LIST,nameobs[ihist],weight=0) != 0)]

    fulllist = vcat(big_chisq_list...)
    if recalc
        @assert sorted != "" 
        fulllist[:,1] = recalc_chired_ipolerr(ipol,sorted,full=true,combine_err=false)
    end

    funcchi = Chisq(Int(floor(fulllist[1,3])))
    if dof != 0
        funcchi = Chisq(Int(dof))
        name=string(name,"_dof$(dof)")
    end
    chiprob = ccdf.(funcchi,fulllist[:,2])
    
    oneprob = []
    zeroprob01 = []
    zeroprob05 = []
    goodrest = []

    for i in 1:length(fulllist[:,1])
        if (isnan(fulllist[i,1]) || fulllist[i,1]^2 == Inf) 
            continue
        end
        if chiprob[i] <= 0.01
            push!(zeroprob01,fulllist[i,1])
        elseif chiprob[i] <= 0.05
            push!(zeroprob05,fulllist[i,1])
        elseif chiprob[i] < 0.9
            push!(goodrest,fulllist[i,1])
        else
            push!(oneprob,fulllist[i,1])
        end
    end
    oneprob = Vector{Float64}(oneprob)
    zeroprob01 = Vector{Float64}(zeroprob01)
    zeroprob05 = Vector{Float64}(zeroprob05)
    goodrest = Vector{Float64}(goodrest)

    binedg = collect(0:xmax/32:xmax)
    honeprob = fit(Histogram,vcat(oneprob,goodrest,zeroprob05,zeroprob01),binedg)
    binedg = collect(honeprob.edges[1])

    hgoodrest = fit(Histogram,vcat(goodrest,zeroprob05,zeroprob01),binedg)
    hzeroprob05 = fit(Histogram,vcat(zeroprob05,zeroprob01),binedg)
    hzeroprob01 = fit(Histogram,vcat(zeroprob01),binedg)

    df = DataFrame(xmin=collect(hgoodrest.edges[1])[1:end-1],xmax=collect(hgoodrest.edges[1])[2:end],yval_oneprob = honeprob.weights)
    CSV.write(string(PREFIX,"/CSV/",name,"_redchi_full_basic.csv"),df)
    plt = plot(honeprob,label="Approximation model \n n=$(sum(honeprob.weights))",xlabel="\$\\chi^{2}_{\\mathrm{red}}\$",ylabel="Entries")
    copy_ticks_lin(plt,plt[1])
    plot!(size=(450,300))
    savefig(string(PREFIX,"/",name,"_redchi_full_basic.png"))
    savefig(string(PREFIX,"/",name,"_redchi_full_basic.pdf"))
    

    df = DataFrame(xmin=collect(hgoodrest.edges[1])[1:end-1],xmax=collect(hgoodrest.edges[1])[2:end],yval_oneprob = honeprob.weights,yval_goodrest = hgoodrest.weights,yval_zeroprob05 = hzeroprob05.weights,yval_zeroprob01 = hzeroprob01.weights
    )
    CSV.write(string(PREFIX,"/CSV/",name,"_redchi_full.csv"),df)

    plot(honeprob,label="p > 0.9,n=$(length(oneprob))",xlabel="\$\\chi^{2}_{\\mathrm{red}}\$",ylabel="Entries")
    plot!(hgoodrest,label="0.9 > p > 0.05,n=$(length(goodrest))")
    plot!(hzeroprob05,label="0.05 > p > 0.01,n=$(length(zeroprob05))")
    plot!(hzeroprob01,label="0.01 > p,n=$(length(zeroprob01))")

    savefig(string(PREFIX,"/",name,"_redchi_full.png"))
    savefig(string(PREFIX,"/",name,"_redchi_full.pdf"))

    honeprob = fit(Histogram,vcat(oneprob,goodrest,zeroprob05),binedg)
    hgoodrest = fit(Histogram,vcat(goodrest,zeroprob05),binedg)
    hzeroprob05 = fit(Histogram,vcat(zeroprob05),binedg)

    df = DataFrame(xmin=collect(honeprob.edges[1])[1:end-1],xmax=collect(honeprob.edges[1])[2:end], yval_oneprob = honeprob.weights, yval_goodrest = hgoodrest.weights, yval_zeroprob05 = hzeroprob05.weights)
    CSV.write(string(PREFIX,"/CSV/",name,"_redchi_nozero.csv"),df)

    plot(honeprob,label="p > 0.9,n=$(length(oneprob))",xlabel="\$\\chi^{2}_{\\mathrm{red}}\$")
    plot!(hgoodrest,label="0.9 > p > 0.05,n=$(length(goodrest))")
    plot!(hzeroprob05,label="0.05 > p > 0.01,n=$(length(zeroprob05))")
    
    savefig(string(PREFIX,"/",name,"_redchi_nozero.png"))
    savefig(string(PREFIX,"/",name,"_redchi_nozero.pdf"))

    honeprob = fit(Histogram,vcat(oneprob,goodrest,zeroprob05),binedg)
    hgoodrest = fit(Histogram,vcat(goodrest,zeroprob05),binedg)
    df = DataFrame(xmin=collect(honeprob.edges[1])[1:end-1],xmax=collect(honeprob.edges[1])[2:end], yval_oneprob = honeprob.weights, yval_goodrest = hgoodrest.weights)
    CSV.write(string(PREFIX,"/CSV/",name,"_redchi_nozero_05.csv"),df)

    plot(honeprob,label="p > 0.9,n=$(length(oneprob))",xlabel="\$\\chi^{2}_{\\mathrm{red}}\$")
    plot!(hgoodrest,label="0.9 > p > 0.05,n=$(length(goodrest))")
    
    savefig(string(PREFIX,"/",name,"_redchi_nozero_05.png"))
    savefig(string(PREFIX,"/",name,"_redchi_nozero_05.pdf"))
end

function recalc_pulls(ipol,sorted;LIST="",skip_neg=true,old=false)
    nameobs = ipol[:,1]
    binedges = ipol[:,2]
    params = get_clean_params_from_mc(sorted)
    pulls = []
    #@assert length(nameobs) == length(sorted)
    for iobs in 1:length(nameobs)
        if occursin("RAW",nameobs[iobs])
            continue
        end
        if LIST != ""
            if old 
                if getweight(LIST,string(iobs),weight=0,numbersafeguard=false) == 0
                    continue
                end
            else
                if getweight(LIST,nameobs[iobs],weight=0,numbersafeguard=false) == 0
                    continue
                end
            end
        end
        binvalues = get_bin_value_for_observable(sorted[iobs])
        binedges = get_bin_edges_for_observables(sorted[iobs])
        binls = binedges[2:end].-binedges[1:end-1]
        binerrs = get_sumw2_for_obervable(sorted[iobs]) #this is sumw2/binl need /binl and sqrt
        binerrs = binerrs[:,3:end]
        nbins = length(binls)

        pulls_per_obs = []

        for ibin in 1:nbins
            covvar = ipol[iobs,3][3][ibin]
            coeffs = ipol[iobs,3][1][ibin,:]
            
            fitvals = g_cubic(params,coeffs)
            mcvals = binvalues[:,2+ibin]

            mcerrs = broadcast(sqrt, binerrs[:,ibin] ./ binls[ibin])
            fiterrs = [sqrt(g_err_cubic(params[i,:],covvar)) for i in 1:length(params[:,1])]


            if skip_neg
                pullsi = []
                for i in 1:length(mcerrs)
                    if mcerrs[i]^ 2 - fiterrs[i]^ 2 < 0
                        continue
                    end 
                    #push!(pullsi,(mcvals[i] - fitvals[i]) / sqrt((mcerrs[i]^2 - fiterrs[i]^2)))
                    push!(pullsi,(mcvals[i] - fitvals[i]) / sqrt((mcerrs[i]^2 - fiterrs[i]^2)))
                end
                push!(pulls_per_obs,Vector{Float64}(pullsi))
                #push!(iden,string(obs," Bin ",ibin," Par ",parnames[ipar]))
            else
                pull = (mcvals .- fitvals) ./ sqrt.(abs.(mcerrs.^2 .- fiterrs.^2))
                push!(pulls_per_obs,pull)
            end            
        end
        push!(pulls,pulls_per_obs)        
    end
    pulls
end

function plot_recalc_pulls(ipol,sorted;folder="recalc_pulls",skip_neg=false,old=false)
    mkpath(folder)
    pulls = recalc_pulls(ipol,sorted,skip_neg=skip_neg,old=old)
    pulls_flat = vcat(vcat(pulls...)...)

    xmax=10
    hrp = fit(Histogram,pulls_flat[abs.(pulls_flat) .< xmax],nbins = 64)
    # if normalize
    #     hrp = normalize(hrp,mode=:pdf)
    # end
    fvals = fit(Normal,pulls_flat[abs.(pulls_flat) .< 5])

    plot(hrp,label="Residual Pulls",seriestype=:step)
    x = -xmax:0.1:xmax
    #plot!(x,pdf(Normal(),x)*maximum(hrp.weights)*sqrt(2*pi),label="Normal Distribution")
    plot!(x,pdf(Normal(fvals.μ,fvals.σ),x)*maximum(hrp.weights)*sqrt(2*pi*fvals.σ^2),label=string("Fitted normal: \n",L"$\mu$ = ","$(round(fvals.μ,digits=3))","\n",L"$\sigma$ = ","$(round(fvals.σ,digits=3))"))
    plot!(size=(420,2/3*400),xlabel=L"$\mathrm{p}$",ylabel="Entries") #legend=:topleft

    fname = "ipol_pulls"
    if skip_neg
        fname = "ipol_pulls_skip_neg"
    end

    savefig(string(folder,"/",fname,".png"))
    savefig(string(folder,"/",fname,".pdf"))

    plot(hrp,label="Residual Pulls",seriestype=:step)
    
    savefig(string(folder,"/",fname,"_pulls.png"))
    savefig(string(folder,"/",fname,"_pulls.pdf"))
end

function recalc_chired_ipolerr(ipol,sorted;
    #folder="chired_ipolerr",
    LIST="",
    combine_err = false,
    full = false,
    filter=true,
    add_rel_err=0.0,
)
    nameobs = ipol[:,1]
    binedges = ipol[:,2]
    @assert length(ipol[1,3]) > 2 #error matrix on ipol[iobs,3][3][ibin]
    params = get_clean_params_from_mc(sorted)
    dof = length(params[:,1])-length(ipol[1,3][1][1,:]) 
    redchi = []
    pvals = []

    for iobs in 1:length(nameobs)
        if occursin("RAW",nameobs[iobs])
            continue
        end
        if LIST != ""
            #if getweight(LIST,string(iobs),weight=0) == 0
            if getweight(LIST,nameobs[iobs],weight=0) == 0
                println(string(nameobs[iobs],"\t",string(iobs)))
                continue
            end
        end
        binvalues = get_bin_value_for_observable(sorted[iobs])
        binedges = get_bin_edges_for_observables(sorted[iobs])
        binls = binedges[2:end].-binedges[1:end-1]
        binerrs = get_sumw2_for_obervable(sorted[iobs]) #this is sumw2/binl need /binl and sqrt
        binerrs = binerrs[:,3:end]
        nbins = length(binls)

        for ibin in 1:nbins
            covvar = ipol[iobs,3][3][ibin]
            coeffs = ipol[iobs,3][1][ibin,:]
            
            fitvals = g_cubic(params,coeffs)
            mcvals = binvalues[:,2+ibin]

            mcerrs = broadcast(sqrt, binerrs[:,ibin] ./ binls[ibin])
            fiterrs = [sqrt(g_err_cubic(params[i,:],covvar)) for i in 1:length(params[:,1])]
            
            adderr = 0
            if add_rel_err > 0.0001
                adderr = mcvals*add_rel_err
            end
            sumerr2 = mcerrs .^ 2 #.+ fiterrs .^2
            if combine_err
                sumerr2 = mcerrs .^ 2 .+ fiterrs .^2
            end
            if adderr != 0
                sumerr2 = mcerrs .^ 2 .+ adderr .^2
            end
            
            if filter
                relerr_entr = (mcerrs./mcvals).^(-2)
                msk = (relerr_entr .> 10)# .& (relerr_entr .< 1000) 
                fitvals=fitvals[msk]
                mcvals=mcvals[msk]
                sumerr2=sumerr2[msk]
                if length(mcvals) == 0
                    continue
                end
            end

            chi2 = sum((fitvals - mcvals).^2 ./ sumerr2)

            pval = ccdf(Chisq(dof),chi2)
            #if typeof(chi2) != typeof(1.0) print(typeof(chi2)) end
            if ((!isnan(chi2/dof) && !isinf(chi2/dof)) || full) push!(redchi,chi2/dof) end
            if ((!isnan(chi2/dof) && !isinf(chi2/dof)) || full) push!(pvals,pval) end
        end
    end
    redchi = Vector{Float64}(redchi)
    return redchi,pvals
end
