

function _AudzeEgliasObjective!(dim::Continous,dist,LHC)
    n,d = size(LHC) 
    dist .= 0.0

    # Squared l-2 norm of distances between all (unique) points
    l = 0
    for i = 2:n
        for j = 1:i-1
            l += 1
            for k = 1:d
                dist[l] += (LHC[i,k]-LHC[j,k])^2
            end
            dist[l] = 1/dist[l]
        end
    end
    output = 1/(sum(dist))
end

function _AudzeEgliasObjective!(dim::Categorical,dist,LHC)
    n,d = size(LHC) 
    dist .= 0.0

    # Squared l-2 norm of distances between all points of the same category
    l = 0

    for i = 2:size(LHC,1)
        for j = 1:i-1
            l += 1
            for k = 1:2
                dist[l] += (LHC[i,k]-LHC[j,k])^2
            end
            dist[l] = 1/dist[l]
        end
    end
    output2 = 1/(sum(dist))
end



"""
    function AudzeEgliasObjective!(dist,LHC::T) where T <: AbstractArray
Return the scalar which should be maximized when using the Audze-Eglias
distance as the objective function. Note this is the inverse of the typical
Audze-Eglias distance which normally is minimized.
"""
function AudzeEgliasObjective!(dist,LHC::T; dims::Array{V,1} =[Continous() for i in 1:size(LHC,2)],
                                            weights::Array{Float64,1}=ones(Float64,length(dims)).*1/length(dims),
                                            ) where T <: AbstractArray where V <: LHCDimension
    
    out = 0.0

    #Compute the objetive function among all points
    continousDimInds = findall(x->x==Continous(),dims)
    if !isempty(continousDimInds)
        out += _AudzeEgliasObjective!(Continous(),dist,LHC)*
                                        sum(weights[continousDimInds])
    end
    

    categoricalDimInds = findall(x->typeof(x)==Categorical,dims)
    for i in categoricalDimInds
        for j = 1:dims[i].levels
            subLHC = @view LHC[LHC[:,i] .== j,:] 
            out += _AudzeEgliasObjective!(dims[i],dist,subLHC)*weights[i]
        end
    end

    return out
end


"""
    function AudzeEgliasObjective(LHC::T) where T <: AbstractArray
Same as AudzeEgliasObjective!(dist,LHC::Array) but creating a new distance array.
"""
function AudzeEgliasObjective(LHC::T;   dims::Array{V,1}=[Continous() for i in 1:size(LHC,2)],
                                        weights::Array{Float64,1}=ones(Float64,length(dims)).*1/length(dims),
                                        ) where T <: AbstractArray where V <: LHCDimension

    n = size(LHC,1) 
    dist = zeros(Float64,Int(n*(n-1)*0.5))

    return AudzeEgliasObjective!(dist,LHC; dims=dims, weights=weights)

end
