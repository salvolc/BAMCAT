function get_paramas_from_file(filename)
    io = open(filename, "r")
    file_string = read(io,String)
    file_lbl = split(file_string,"\n")
    name=[]
    value=[]
    for j in 1:length(file_lbl)
        line = file_lbl[j]
        if line != ""
            push!(name,split(line," ")[1])
            push!(value,parse(Float64,split(line," ")[2]))
        end
    end
    hcat(name,value)
end

function execute(cmd::Cmd)
  out = Pipe()
  err = Pipe()

  process = run(pipeline(ignorestatus(cmd), stdout=out, stderr=err))
  close(out.in)
  close(err.in)

  (
    stdout = String(read(out)),
    stderr = String(read(err)),
    code = process.exitcode
  )
end

function generate_inputs_for_grid(parameter_ranges;
    analysis_config_path::String="./",
    input_grid_path::String="input_grid",
    n_per_par::Int=5
    )

    npars = length(parameter_ranges)
    thispath = pwd()
    execute(`mkdir -p $input_grid_path`)
    full_grid_path = string(thispath,"/",input_grid_path)
    cd(analysis_config_path)
    full_anaysis_path = pwd()

    NAME          = execute(`python -c "from analysis_config import *; print(NAME)"`).stdout[1:end-1]
    CONTAINER     = execute(`python -c "from analysis_config import *; print(CONTAINER)"`).stdout[1:end-1]
    CACHELOCATION = execute(`python -c "from analysis_config import *; print(CACHELOCATION)"`).stdout[1:end-1]
    ENERGIES      = execute(`python -c "from analysis_config import *;  [print(en)  for en in energies]"`).stdout[1:end-1]
    PREFIXES      = split(replace(replace(execute(`python -c "from analysis_config import *; print(prefixes)"`).stdout[2:end-2],"'"=>""),","=>""))
    #VARIABLES     = split(replace(replace(execute(`python -c "from analysis_config import *; print(variables)"`).stdout[2:end-2],"'"=>""),","=>""))
    VARIABLES_RAW = split(replace(execute(`python -c "from analysis_config import *; print(variables)"`).stdout[2:end-2],"'"=>""),"],")

    VARIABLES=[]
    for i in 1:length(VARIABLES_RAW)
        push!(VARIABLES,Vector{String}(split(replace(replace(replace(VARIABLES_RAW[i],"["=>""),"]"=>"")," "=>""),",")))
    end

    low = [parameter_ranges[i].a for i in 1:npars]
    high = [parameter_ranges[i].b for i in 1:npars]
    ranges = Matrix{Float64}(undef,npars,n_per_par)

    for i in 1:npars
        ranges[i,:] = collect(range(low[i],stop=high[i],length=n_per_par))
    end
    values = Vector{Float64}(undef,npars)

    cd(full_grid_path)

    for i in 1:npars
        for j in 1:n_per_par
            n_mc = (i-1)*n_per_par+(j-1)
            mc_name = string(lpad(n_mc, 5, "0"))
            execute(`mkdir -p $mc_name`)
            execute(`cp $full_anaysis_path/bin/do.sh $mc_name`)
            execute(`chmod +x $mc_name/do.sh`)

            open(string(mc_name,"/cache.txt"), "w") do io
                write(io, string(CACHELOCATION,"Herwig-cache_",NAME,".tar.gz"))
            end
            open(string(mc_name,"/container.txt"), "w") do io
                write(io, string(CONTAINER))
            end
            open(string(mc_name,"/energies.txt"), "w") do io
                write(io, string(ENERGIES))
            end
            open(string(mc_name,"/seed.txt"), "w") do io
                write(io, string(floor(Int,rand()*1e8)))
            end

            values = ranges[:,trunc(Int,(n_per_par/2)+1)]
            values[i] = ranges[i,j]

            io = open(string(mc_name,"/params.dat"),"w")
            for i_write in 1:npars
                write(io,string(VARIABLES[i_write][1]," ",values[i_write],"\n"))
            end
            close(io)

            io = open(string(mc_name,"/params.in"),"w")
            for i_write in 1:npars
                [write(io,string("set ",PREFIXES[i_write],VARS," ",values[i_write],"\n")) for VARS in VARIABLES[i_write]]
            end
            close(io)
            
            execute(`cat $full_anaysis_path/cards/Analysis_$NAME".in" ">>" $mc_name/params.in`)
            analyses = execute(`cat $full_anaysis_path/cards/Analysis_$NAME".in"`).stdout
            open(string(mc_name,"/params.in"), "a") do io
                write(io, string(analyses))
            end
        end
    end
    cd(thispath)
end

function run_one_mc(path::String)
    cd(path)
    run(`./do.sh`)
    cd("../")
end

function run_grid_inputs(inputpath::String="./input_grid/")
    currentpath = pwd()
    mc_nums = cd(readdir,inputpath)
    cd(inputpath)
    Folds.collect(run_one_mc(string(i)) for i in mc_nums)
    cd(currentpath)
end

function run_grid_inputs_st(inputpath::String="./input_grid/")
    currentpath = pwd()
    mc_nums = cd(readdir,inputpath)
    cd(inputpath)
    for i in mc_nums
        cd(string(i))
        run(`./do.sh`)
        cd("..")
    end
    cd(currentpath)
end

function getweight(file::String, obsname::String;weight=1,n_weigths=2,numbersafeguard=true)
    if numbersafeguard
        obs = 1
        try
            obs = parse(Float64,obsname)
        catch
        end
        @assert typeof(obs) != Float64 "Your Input was a Number is this is intended please user numbersafeguard = false as a flag for getweight()"
    end
    if weight == -1
        return 1
    end
    io = open(file,"r")
    file_lbl = readlines(io)
    file_format = [split(replace(line," "=>""),"&") for line in file_lbl if length(split(line,"&")) > 2][2:end]
    ana_format = [[i[1],i[2]] for i in file_format]
    index = findall(i -> occursin(obsname,string("/",i[1],"/",i[2],"/")), ana_format)
    if index == []
        return 0
    end
    if weight == 0
        return 1
    end
    parse(Int64,file_format[index][1][end-n_weigths+weight])
end

function getobsname(file::String, obsname::String)
    io = open(file,"r")
    file_lbl = readlines(io)
    file_format = [split(line,"&") for line in file_lbl if length(split(line,"&")) > 2]
    index = findall(i -> occursin(obsname,string("/",replace(i[1]," "=>""),"/",replace(i[2]," "=>""),"/")), file_format)
    if index == []
        return 0
    end
    file_format[index][1][3]
end

function give_category_of_obs(;obsname="",LIST="/ceph/groups/e4/users/slacagnina/overH70222/category_list.txt")
    lbl = readlines(LIST)
    for i in 1:length(lbl)
        ls = [replace(j," " => "") for j in split(lbl[i],"&")]
        if length(ls) > 2 && occursin("$(ls[1])/$(ls[2])",obsname) 
            return ls[end] 
        end 
    end
end

function chi2_for_grid(ipol,grid::Vector{Observablecontainer})
    hist_names = ipol[:,1]
    params = get_clean_params_from_mc(grid)
    pars = get_clean_params_from_mc(grid)
    chi2 = []
    chi2_marg = []
    for iobs in 1:length(hist_names)
        obs = hist_names[iobs]
        binvalues = get_bin_value_for_observable(grid,obs)
        binedges = get_bin_edges_for_observables(grid,obs)
        binerrs = get_sumw2_for_obervable(grid,obs)
        nparams = length(pars[1,:])
        nbins = length(binedges)-1
        nparcombos = length(pars[:,1])
        nvalsperpar = Int(length(pars[:,1])/nparams)
        parvals = pars[:,:]
        coeffs = ipol[:,3][iobs][1]
        fitvals = [g(parvals,coe) for coe in eachrow(coeffs)]

        for ibin in 1:nbins
            if sum(binerrs[:,ibin+2] .< 0) > 0
                print("Skipped $obs in bin $ibin for negative sumw2.")
                continue
            end
            sumchi2=0
            sumchi2err=0
            for ipar in 1:nparams
                parbeg = 1+nvalsperpar*(ipar-1)
                parend = nvalsperpar*ipar
                parval = parvals[parbeg:parend,:][:,ipar]
                fitval = fitvals[ibin][parbeg:parend]
                orgval = binvalues[parbeg:parend,2+ibin]
                orgvalerr = broadcast(sqrt,binerrs[parbeg:parend,2+ibin])

                chisq_val = chisq(orgval,fitval)
                chisqerr_val = chisqerr(orgval,fitval,orgvalerr)
                
                chisqprob = ccdf(Chisq(nvalsperpar-1),chisq_val)
                chired = chisq_val/ipol[iobs,3][2][ibin,3]

                push!(chi2,(name=obs,bin=ibin,par=ipar,chisq=chisq_val,chisqerr=chisqerr_val,chisqprob=chisqprob,chired=chired))

                sumchi2 += chisq_val
                sumchi2err += chisqerr_val
            end
            sumchi2prob = ccdf(Chisq(nvalsperpar*nparams),sumchi2)
            sumchi2red = sumchi2err/ipol[iobs,3][2][ibin,3]
            push!(chi2_marg,(name=obs,bin=ibin,chisq=sumchi2,chisqerr=sumchi2err,chisqprob=sumchi2prob,chired=sumchi2red))
        end
    end
    return chi2,chi2_marg
 end

 function collapse_mcpars(MCPATH::String)
    mc_folders = readdir(MCPATH)
    for i in mc_folders
        params_file=string(MCPATH,i,"/params.dat")
        file_par = get_paramas_from_file(params_file)
        v = file_par[:,2]
        indices = findfirst.(isequal.(unique(v)), [v])
        data = file_par[indices,:]
        io = open(string(MCPATH,i,"/params.dat"),"w")
        [write(io,string(data[i,1]," ",data[i,2],"\n")) for i in 1:length(data[:,1])]
        close(io) 
    end
end

function chisq_data_mc(MC_PATH::String,RIVET_DATA_PATH::String)
    nom_hists = get_all_histograms(string(MC_PATH,"/output.yoda"))
    nom_sumw2 = get_all_sumw2(string(MC_PATH,"/output.yoda"))
    ipol_histos_names = nom_hists[:,1]

    analyses_input_file = string(MC_PATH,"/params.in")
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

    #bad_data=split(read(open("bad_data.txt","r"),String),"\n")
    chisq_nom=[]

    for i_hist = 1:n_hist
        n_bin = n_bins[i_hist]
        ipol_index =
            findall(i -> i == data_histos_names[i_hist][5:end], ipol_histos_names)[1]
        weight = getweight(LIST,string(ipol_histos_names[ipol_index]),weight=1)
        if weight == 0
            continue
        end
        #if length(findall(x->occursin(x,ipol_histos_names[ipol_index]),bad_data))>0
        #    continue
        #end
        data_bins = data_histos[i_hist,2].weights
        nom_bins = nom_hists[ipol_index,2].weights
        nom_err = nom_sumw2[ipol_index,2]

        chisq_nom_i = 0

        for i_bin = 1:n_bin
            if  (data_histos[i_hist, 2].weights[i_bin] != 0 
                && ccdf(Chisq(ipol[ipol_index, 3][2][i_bin,6]-1),ipol[ipol_index, 3][2][i_bin,5]) > 0.05)
                #&& ipol[ipol_index, 3][2][i_bin, 1] < 25)
                #&& !isapprox(0,ipol[ipol_index, 3][2][i_bin, 1],atol=0.00001)
                chisq_nom_i += (nom_bins[i_bin] - data_bins[i_bin])^2 / data_bins[i_bin]
            end
        end

        push!(chisq_nom,chisq_nom_i)
    end
    chisq_nom
end

function chisq_data_mc_folder(MC_FOLDER::String,RIVET_DATA_PATH::String)
    folders = cd(readdir,MC_FOLDER)
    chisq_vals = Vector{Vector{Float64}}(undef,length(folders))
    for i in 1:length(folders)
        MC_PATH = string(MC_FOLDER,"/",folders[i],"/")
        chisq_vals[i] = chisq_data_mc(MC_PATH,RIVET_DATA_PATH)
    end
    chisq_vals
end

function chisq_data_mc_folder_parallel(MC_FOLDER::String,RIVET_DATA_PATH::String)
    folders = cd(readdir,MC_FOLDER)
    Folds.collect( chisq_data_mc(string(MC_FOLDER,"/",folders[i],"/"),RIVET_DATA_PATH) for i in 1:length(folders))
end

function grid_chi2(ipol,grid::Vector{Observablecontainer};par_mask=[true for i in 1:length(grid[1].hc[1].parname)])
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

    c2all = []
    for iobs in 1:length(hist_names)
        obs = hist_names[iobs]
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
        fitvals = [g(parvals,coe) for coe in eachrow(coeffs)]
        c2obs = Vector{Vector{Float64}}(undef,nbins)
        binl = binedges[2:end]-binedges[1:end-1]

        for ibin in 1:nbins
            if sum(binerrs[:,ibin+2] .< 0) > 0
                print("Skipped $obs in bin $ibin for negative sumw2.")
                continue
            end
            c2bin = Vector{Float64}(undef,nparams)
            for ipar in 1:nparams
                parbeg = 1+nvalsperpar*(ipar-1)
                parend = nvalsperpar*ipar
                parval = parvals[parbeg:parend,:][:,ipar]
                fitval = fitvals[ibin][parbeg:parend]
                orgval = binvalues[parbeg:parend,2+ibin]
                orgvalerr = broadcast(sqrt,binerrs[parbeg:parend,2+ibin]) ./ binl[ibin]
                c2 = chisqerr(fitval,orgval,orgvalerr)
                c2bin[ipar] = c2
            end
            c2obs[ibin] = c2bin
        end
        push!(c2all,c2obs)
    end
    return c2all
end

function ipol_chisqprob_lists(ipol;name="",folder="fit",dof=0)
    PREFIX = folder
    mkpath(folder)
    nameobs = ipol[:,1]
    binedges = ipol[:,2]

    funcchi = Chisq(Int(floor(ipol[1,3][2][1,3])))
    if dof != 0
        funcchi = Chisq(Int(dof))
        name=string(name,"_dof$(dof)")
    end

    f = open(string(PREFIX,"/",name,"_onered_list.txt"),"w")
    for i in 1:length(nameobs)
        if occursin("RAW",nameobs[i])
            continue
        end
        bedg = binedges[i]
        nn = nameobs[i]
        write(f,string("$nn [ "))
        for j in 1:length(bedg)-1
            if isapprox(ipol[i,3][2][j,1],1,atol=0.2)
                chiprob = ccdf(funcchi,ipol[i,3][2][j,2])
                write(f,string(round(chiprob,digits=3),"[$j], "))
            end
        end
    write(f,string(" ] \n"))
    end
    close(f)


    f = open(string(PREFIX,"/",name,"_onered_list_inverse.txt"),"w")
    for i in 1:length(nameobs)
        if occursin("RAW",nameobs[i])
            continue
        end
        bedg = binedges[i]
        nn = nameobs[i]
        write(f,string("$nn [ "))
        for j in 1:length(bedg)-1
            if !isapprox(ipol[i,3][2][j,1],1,atol=0.2)
                chiprob = ccdf(funcchi,ipol[i,3][2][j,2])
                write(f,string(round(chiprob,digits=3),"[$j], "))
            end
        end
    write(f,string(" ] \n"))
    end
    close(f)

    f = open(string(PREFIX,"/",name,"_onered_list_inverse_nozero.txt"),"w")
    for i in 1:length(nameobs)
        if occursin("RAW",nameobs[i])
            continue
        end
        bedg = binedges[i]
        nn = nameobs[i]
        write(f,string("$nn [ "))
        for j in 1:length(bedg)-1
            if !isapprox(ipol[i,3][2][j,1],1,atol=0.2) && !isapprox(ipol[i,3][2][j,1],0,atol=0.2)
                chiprob = ccdf(funcchi,ipol[i,3][2][j,2])
                write(f,string(round(chiprob,digits=3),"[$j], "))
            end
        end
    write(f,string(" ] \n"))
    end
    close(f)


    f = open(string(PREFIX,"/",name,"_zerored_list.txt"),"w")
    for i in 1:length(nameobs)
        if occursin("RAW",nameobs[i])
            continue
        end
        bedg = binedges[i]
        nn = nameobs[i]
        write(f,string("$nn [ "))
        for j in 1:length(bedg)-1
            if ipol[i,3][2][j,1] < 0.8 #isapprox(ipol[i,3][2][j,1],0,atol=0.2)
                chiprob = ccdf(funcchi,ipol[i,3][2][j,2])
                write(f,string(round(chiprob,digits=3),"[$j], "))
            end
        end
    write(f,string(" ] \n"))
    end
    close(f)
end

function ipol_chisqprob_lists_constr(ipol;name="",folder="fit",LIST="/ceph/groups/e4/users/slacagnina/overH70222/longlist_raw.txt",dof=0)
    PREFIX = folder
    name=string(name,"_constr")
    mkpath(folder)
    nameobs = ipol[:,1]
    binedges = ipol[:,2]

    funcchi = Chisq(Int(floor(ipol[1,3][2][1,3])))
    #funcchi = Chisq(Int(floor(ipol[1,3][2][1,6])))
    if dof != 0
        funcchi = Chisq(Int(dof))
        name=string(name,"_dof$(dof)")
    end

    f = open(string(PREFIX,"/",name,"_onered_list.txt"),"w")
    for i in 1:length(nameobs)
        if (occursin("RAW",nameobs[i]) || getweight(LIST,nameobs[i],weight=0) == 0)
            continue
        end
        bedg = binedges[i]
        nn = nameobs[i]
        write(f,string("$nn [ "))
        for j in 1:length(bedg)-1
            if isapprox(ipol[i,3][2][j,1],1,atol=0.2)
                chiprob = ccdf(funcchi,ipol[i,3][2][j,2])
                write(f,string(round(chiprob,digits=3),"[$j], "))
            end
        end
    write(f,string(" ] \n"))
    end
    close(f)


    f = open(string(PREFIX,"/",name,"_onered_list_inverse.txt"),"w")
    for i in 1:length(nameobs)
        if (occursin("RAW",nameobs[i]) || getweight(LIST,nameobs[i],weight=0) == 0)
            continue
        end
        bedg = binedges[i]
        nn = nameobs[i]
        write(f,string("$nn [ "))
        for j in 1:length(bedg)-1
            if !isapprox(ipol[i,3][2][j,1],1,atol=0.2)
                chiprob = ccdf(funcchi,ipol[i,3][2][j,2])
                write(f,string(round(chiprob,digits=3),"[$j], "))
            end
        end
    write(f,string(" ] \n"))
    end
    close(f)

    f = open(string(PREFIX,"/",name,"_onered_list_inverse_nozero.txt"),"w")
    for i in 1:length(nameobs)
        if (occursin("RAW",nameobs[i]) || getweight(LIST,nameobs[i],weight=0) == 0)
            continue
        end
        bedg = binedges[i]
        nn = nameobs[i]
        write(f,string("$nn [ "))
        for j in 1:length(bedg)-1
            if !isapprox(ipol[i,3][2][j,1],1,atol=0.2) && !isapprox(ipol[i,3][2][j,1],0,atol=0.2)
                chiprob = ccdf(funcchi,ipol[i,3][2][j,2])
                write(f,string(round(chiprob,digits=3),"[$j], "))
            end
        end
    write(f,string(" ] \n"))
    end
    close(f)


    f = open(string(PREFIX,"/",name,"_zerored_list.txt"),"w")
    for i in 1:length(nameobs)
        if (occursin("RAW",nameobs[i]) || getweight(LIST,nameobs[i],weight=0) == 0)
            continue
        end
        bedg = binedges[i]
        nn = nameobs[i]
        write(f,string("$nn [ "))
        for j in 1:length(bedg)-1
            if ipol[i,3][2][j,1] < 0.8 #isapprox(ipol[i,3][2][j,1],0,atol=0.2)
                chiprob = ccdf(funcchi,ipol[i,3][2][j,2])
                write(f,string(round(chiprob,digits=3),"[$j], "))
            end
        end
    write(f,string(" ] \n"))
    end
    close(f)
end

select_chainid(c...) = x -> x.info.chainid in reduce(vcat, c)

@userplot LogPosteriorPlot

@recipe function f(
    lpp::LogPosteriorPlot; 
    chainids=unique(lpp.args[1].info.chainid),
)
    samples = lpp.args[1]
        
    for cid in chainids
        suc = filter(select_chainid(cid), samples)
    
        @series begin
            seriestype --> :line
            label --> "Chain ID: $cid"
            xguide --> "Steps"
            yguide --> "log(posterior)"
            suc.logd
        end
    end
end

function read_hepdata_csv(csv_file_path)
    alllines = readlines(csv_file_path)
    headerindx = findall(x->occursin("sys",x) || occursin("error",x),alllines)
    alldf = Vector{DataFrame}(undef,length(headerindx))
    for i in 1:length(alldf)
        stop = headerindx[i]
        while alllines[stop+1] != ""
            stop += 1
        end
        csf = CSV.read(IOBuffer(join(alllines[headerindx[i]:stop+1],"\n")),DataFrame)
        alldf[i] = DataFrame(csf)
    end
    deleteat!(alldf,findall(x->nrow(x)==0,alldf))
end

function read_hepdata_csv(csv_folder_path,csv_file_name)
    read_hepdata_csv(string(csv_folder_path,csv_file_name))
end

function read_hepdata_csv_exp(csv_folder_path)
    table_list = readdir(csv_folder_path)
    all_tables = Vector{Vector{DataFrame}}(undef,length(table_list))
    for i in 1:length(table_list)
        all_tables[i] = read_hepdata_csv(string(csv_folder_path,"/Table$(i).csv"))
    end
    all_tables
end

function ratio_histogram(h1,h2)
    @assert isapprox(h1.edges[1],h2.edges[1])
    binning = h1.edges[1]
    h1w = h1.weights 
    h2w = h2.weights
    rtw = h1w ./ h2w
    #fit(Histogram,rtw,binning)
    Histogram(binning,rtw)
end

function twiny(sp::Plots.Subplot)
    sp[:top_margin] = max(sp[:top_margin], 30Plots.px)
    plot!(sp.plt, inset = (sp[:subplot_index], bbox(0,0,1,1)))
    twinsp = sp.plt.subplots[end]
    twinsp[:xaxis][:mirror] = true
    twinsp[:background_color_inside] = RGBA{Float64}(0,0,0,0)
    Plots.link_axes!(sp[:yaxis], twinsp[:yaxis])
    twinsp
end
twiny(plt::Plots.Plot = current()) = twiny(plt[1])


function copy_ticks(sp::Plots.Subplot)
    ptx = twinx(sp)
    plot!(ptx,xlims=xlims(plt),ylims=ylims(plt),xformatter=_->"",yformatter=_->"")
    pty = twiny(sp)
    plot!(pty,xlims=xlims(plt),ylims=ylims(plt),xformatter=_->"",yformatter=_->"")
end

function copy_ticks(plt::Plots.Plot,sp::Plots.Subplot)
    ptx = twinx(sp)
    plot!(ptx,xlims=xlims(plt),ylims=ylims(plt),xformatter=_->"",yformatter=_->"")
    pty = twiny(sp)
    plot!(pty,xlims=xlims(plt),ylims=ylims(plt),xformatter=_->"",yformatter=_->"")
end
copy_ticks(plt::Plots.Plot = current()) = copy_ticks(plt[1])

function exists(a,b)
    length(findall(x->x==a,b)) > 0
end

function exists_str(a,b)
    length(findall(x->occursin(a,x),b)) > 0
end

function give_all_str(a,b)
    findall(x->occursin(a,x),b)
end


function copy_ticks_log(plt::Plots.Plot,sp::Plots.Subplot;yminorticks=0,xminorticks=0)
    ptx = twinx(sp)
    plot!(ptx,xlims=xlims(plt),ylims=ylims(plt),xformatter=_->"",yformatter=_->"",yscale=:log10,yminorticks=yminorticks,xminorticks=xminorticks)
    pty = twiny(sp)
    plot!(pty,xlims=xlims(plt),ylims=ylims(plt),xformatter=_->"",yformatter=_->"",yscale=:log10,yminorticks=yminorticks,xminorticks=xminorticks)
end

function copy_ticks_lin(plt::Plots.Plot,sp::Plots.Subplot;yminorticks=0,xminorticks=0)
    ptx = twinx(sp)
    plot!(ptx,xlims=xlims(plt),ylims=ylims(plt),xformatter=_->"",yformatter=_->"",yscale=:lin,yminorticks=yminorticks,xminorticks=xminorticks)
    pty = twiny(sp)
    plot!(pty,xlims=xlims(plt),ylims=ylims(plt),xformatter=_->"",yformatter=_->"",yscale=:lin,yminorticks=yminorticks,xminorticks=xminorticks)
end