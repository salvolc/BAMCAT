
function give_histogram_names(file_lbl)
    histogram_names = []
    for i in 1:length(file_lbl)
        line = file_lbl[i]
        if occursin("HISTO1D",line) && occursin("BEGIN",line)
            push!(histogram_names,split(line)[3])
        end
    end
    return histogram_names
end

function give_histogram(i,file_lbl,histogram_names=give_histogram_names(file_lbl))
    histblock=false
    datablock=false
    data = []

    for j in 1:length(file_lbl)
        line = file_lbl[j]

        if occursin("BEGIN",line)
            if histogram_names[i] == split(line)[3]
                histblock = true
            end
        end

        if histblock && occursin("# xlow",line)
            datablock=true
            continue
        end

        if occursin("END",line)
            histblock=false
            datablock=false
        end

        if datablock
            push!(data,parse.(Float64,split(line,"\t")))
        end
    end
    data = hcat(data...)
    binl = data[2,:] .- data[1,:]
    Histogram(vcat(data[1,:],data[2,length(data[1,:])]),data[3,:]./binl)
end

function give_histogram_csv(i,file_lbl,histogram_names=give_histogram_names(file_lbl))
    histblock=false
    datablock=false
    data = []

    for j in 1:length(file_lbl)
        line = file_lbl[j]

        if occursin("BEGIN",line)
            if histogram_names[i] == split(line)[3]
                histblock = true
                datablock = true
                continue
            end
        end

        if occursin("END",line)
            histblock=false
            datablock=false
        end

        if datablock
            push!(data,parse.(Float64,split(line," ")))
        end
    end
    data = hcat(data...)
    edge = data[1,1:end]
    values = data[2,1:end-1]
    Histogram(edge,values)
end


function get_all_histograms(filename)
    io = open(filename, "r")
    file_string = read(io,String)
    file_lbl = split(file_string,"\n")
    names = give_histogram_names(file_lbl)
    if endswith(filename,".csv")
        histograms = [give_histogram_csv(i,file_lbl,names) for i in 1:length(names)]
    else
        histograms = [give_histogram(i,file_lbl,names) for i in 1:length(names)]
    end
    return hcat(names,histograms)
end


function give_sumw2(i,file_lbl,histogram_names=give_histogram_names(file_lbl),unweighted=false)
    histblock=false
    datablock=false
    data = []
    for j in 1:length(file_lbl)
        line = file_lbl[j]
        if occursin("BEGIN",line)
            if histogram_names[i] == split(line)[3]
                histblock = true
            end
        end
        if histblock && occursin("# xlow",line)
            datablock=true
            continue
        end
        if occursin("END",line)
            histblock=false
            datablock=false
        end
        if datablock
            push!(data,parse.(Float64,split(line,"\t")))
        end
    end
    data = hcat(data...)
    binl = data[2,:] .- data[1,:]
    data[4,:]./binl
end

function give_sumw2_csv(i,file_lbl,histogram_names=give_histogram_names(file_lbl),unweighted=false)
    # temp solution because no sumw2 in pythia output
    histblock=false
    datablock=false
    data = []
    for j in 1:length(file_lbl)
        line = file_lbl[j]
        if occursin("BEGIN",line)
            if histogram_names[i] == split(line)[3]
                histblock = true
                datablock = true
                continue
            end
        end
        if occursin("END",line)
            histblock=false
            datablock=false
        end
        if datablock
            push!(data,parse.(Float64,split(line," ")))
        end
    end
    data = hcat(data...)
    data[3,1:end-1]
end

function get_all_sumw2(filename)
    io = open(filename, "r")
    file_string = read(io,String)
    file_lbl = split(file_string,"\n")
    names = give_histogram_names(file_lbl)
    if endswith(filename,".csv")
        histograms = [give_sumw2_csv(i,file_lbl,names) for i in 1:length(names)]
    else
        histograms = [give_sumw2(i,file_lbl,names) for i in 1:length(names)]
    end
    return hcat(names,histograms)
end


function give_ref_histogram_names(file_lbl)
    histogram_names = []
    for i in 1:length(file_lbl)
        line = file_lbl[i]
        if occursin("YODA_SCATTER2D_V2",line) && occursin("BEGIN",line)
            push!(histogram_names,split(line)[3])
        end
    end
    return histogram_names
end

function give_ref_histogram(i,file_lbl,histogram_names=give_ref_histogram_names(file_lbl))
    histblock=false
    datablock=false
    data = []

    for j in 1:length(file_lbl)
        line = file_lbl[j]

        if occursin("BEGIN",line)
            if histogram_names[i] == split(line)[3]
                histblock = true
            end
        end

        if histblock && occursin("# xval",line)
            datablock=true
            continue
        end

        if occursin("END",line)
            histblock=false
            datablock=false
        end

        if datablock
            push!(data,parse.(Float64,split(line,"\t")))
        end
    end
    data = hcat(data...)
    Histogram(vcat(data[1,:].-data[2,:],data[1,end]+data[3,end]),data[4,:])
end


function get_all_ref_histograms_to_mc(filename,filename_mc)
    io_mc = open(filename_mc, "r")
    file_string_mc = read(io_mc,String)
    file_lbl_mc = split(file_string_mc,"\n")
    names_mc = give_histogram_names(file_lbl_mc)
    get_all_ref_histograms_to_names(filename,names_mc)
end

function get_all_ref_histograms_to_names(filename,names_mc)
    io = open(filename, "r")
    file_string = read(io,String)
    file_lbl = split(file_string,"\n")
    names = give_ref_histogram_names(file_lbl)

    ind = []
    for n in 1:length(names)
        if(sum(occursin.(Ref(string(split(names[n],"/")[end-1],"/",split(names[n],"/")[end])),names_mc)) > 0)
            push!(ind,n)
        end
    end

    histograms = [give_ref_histogram(i,file_lbl,names) for i in ind]
    names = [names[i] for i in ind]
    return hcat(names,histograms)
end



function give_ref_sumw2(i,file_lbl,histogram_names=give_ref_histogram_names(file_lbl))
    histblock=false
    datablock=false
    data = []

    for j in 1:length(file_lbl)
        line = file_lbl[j]

        if occursin("BEGIN",line)
            if histogram_names[i] == split(line)[3]
                histblock = true
            end
        end

        if histblock && occursin("# xval",line)
            datablock=true
            continue
        end

        if occursin("END",line)
            histblock=false
            datablock=false
        end

        if datablock
            push!(data,parse.(Float64,split(line,"\t")))
        end
    end
    data = hcat(data...)
    return data[5,:]
end

function get_all_ref_sumw2_to_mc(filename,filename_mc)
    io_mc = open(filename_mc, "r")
    file_string_mc = read(io_mc,String)
    file_lbl_mc = split(file_string_mc,"\n")
    names_mc = give_histogram_names(file_lbl_mc)
    get_all_ref_sumw2_to_names(filename,names_mc)
end

function get_all_ref_sumw2_to_names(filename,names_mc)
    io = open(filename, "r")
    file_string = read(io,String)
    file_lbl = split(file_string,"\n")
    names = give_ref_histogram_names(file_lbl)

    ind = []
    for n in 1:length(names)
        if(sum(occursin.(Ref(string(split(names[n],"/")[end-1],"/",split(names[n],"/")[end])),names_mc)) > 0)
            push!(ind,n)
        end
    end

    histograms = [give_ref_sumw2(i,file_lbl,names) for i in ind]
    names = [names[i] for i in ind]
    return hcat(names,histograms)
end
