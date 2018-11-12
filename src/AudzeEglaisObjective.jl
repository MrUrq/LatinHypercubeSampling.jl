

function _AudzeEglaisObjective!(dim::Continuous,dist,LHC)
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

function _AudzeEglaisObjective!(dim::Categorical,dist,LHC)
    n,d = size(LHC) 
    dist .= 0.0

    # Squared l-2 norm of distances between all points of the same category
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
    if output == Inf
        output = 0
    end
    output    
end



"""
    function AudzeEglaisObjective!(dist,LHC::T) where T <: AbstractArray
Return the scalar which should be maximized when using the Audze-Eglais
distance as the objective function. Note this is the inverse of the typical
Audze-Eglais distance which normally is minimized.
"""
function AudzeEglaisObjective!(dist,LHC::T; dims::Array{V,1} =[Continuous() for i in 1:size(LHC,2)],
                                            interSampleWeight::Float64=1.0,
                                            ) where T <: AbstractArray where V <: LHCDimension
    
    out = 0.0

    #Compute the objective function among all points
    out += _AudzeEglaisObjective!(Continuous(),dist,LHC)*interSampleWeight
    

    categoricalDimInds = findall(x->typeof(x)==Categorical,dims)
    for i in categoricalDimInds
        for j = 1:dims[i].levels
            subLHC = @view LHC[LHC[:,i] .== j,:] 
            out += _AudzeEglaisObjective!(dims[i],dist,subLHC)*dims[i].weight
        end
    end

    return out
end


"""
    function AudzeEglaisObjective(LHC::T) where T <: AbstractArray
Same as AudzeEglaisObjective!(dist,LHC::Array) but creating a new distance array.
"""
function AudzeEglaisObjective(LHC::T;   dims::Array{V,1}=[Continuous() for i in 1:size(LHC,2)],
                                        interSampleWeight::Float64=1.0,
                                        ) where T <: AbstractArray where V <: LHCDimension

    n = size(LHC,1) 
    dist = zeros(Float64,Int(n*(n-1)*0.5))

    return AudzeEglaisObjective!(dist,LHC; dims=dims, interSampleWeight=interSampleWeight)

end
