# function g_pro_bin(a,n)
#     p = n
#     n = a.p
#     nparams = 0
#     param = Vector{Float64}(undef,length(a.p))
#     if length(n[1,:]) == 1
#         nparams = length(n)
#         param = n
#     else
#         nparams = length(n[1,:])
#         for i in 1:nparams
#             push!(param,n[:,i])
#         end
#     end

#     m = ones(length(param[1]))
#     if length(param[1]) == 1
#         m = [1.]
#     end

#     for i in 1:nparams
#         push!(m,param[i])
#     end
#     for i in 1:nparams
#         for j in i+1:nparams
#                 push!(m,param[i].*param[j])
#         end
#     end
#     for i in 1:nparams
#         push!(m,param[i].*param[i])
#     end

#     hcat(m...)*p
# end

# function g_pro_bin(a,n)
#     a = collect(a)
#     if length(a) == 8
#         if a[3] > a[2]/2
#             return 0
#         end
#     end
#     m = g_cubic(a)
#     m'*n
# end

function g_pro_bin(a,n)
    a = collect(a)
    m = g_quad(a)
    m'*n
end

function g_pro_bin(a,n)
    #pars = StaticArrays.SVector{length(collect(a))}(collect(a))
    pars = collect(a)
    if length(a) == 8
        #pars = StaticArrays.SVector{length(collect(a))}(collect(a)[[1,2,8,3,4,5,6,7]])
        pars = collect(a)[[1,2,8,3,4,5,6,7]]
    end
    m = g_cubic(pars)
    m'*n
end

function make_per_bin_function(coeffs)
    return x->g_pro_bin(x,coeffs)[1]
end

function test_l()
    for i in 1:500
        g_pro_bin((p=[1,2,3,4,5,6,7],),ones(36))
    end
end