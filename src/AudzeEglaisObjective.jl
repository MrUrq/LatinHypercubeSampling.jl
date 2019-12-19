@inline function _AudzeEglaisDist(LHC)
    n,d = size(LHC) 
    dist = 0.0
    dist_tmp = 0.0

    # Squared l-2 norm of distances between all (unique) points
    for i = 2:n
        for j = 1:i-1
            for k = 1:d
                @inbounds dist_tmp += (LHC[i,k]-LHC[j,k])^2
            end
            @inbounds dist += 1/dist_tmp
            dist_tmp = 0.0
        end
    end
    output = 1/dist
    return output
end

function _AudzeEglaisObjective(dim::Continuous,LHC)
    output = _AudzeEglaisDist(LHC)    
    return output
end

function _AudzeEglaisObjective(dim::Categorical,LHC)
    output = _AudzeEglaisDist(LHC)     
    output == Inf ? 0 : output    
end

"""
    function AudzeEglaisObjective!(LHC::T) where T <: AbstractArray
Return the scalar which should be maximized when using the Audze-Eglais
distance as the objective function. Note this is the inverse of the typical
Audze-Eglais distance which normally is minimized.
"""
function AudzeEglaisObjective(LHC::T; dims::Array{V,1} =[Continuous() for i in 1:size(LHC,2)],
                                            interSampleWeight::Float64=1.0,
                                            ) where T <: AbstractArray where V <: LHCDimension
    
    out = 0.0

    #Compute the objective function among all points
    out += _AudzeEglaisObjective(Continuous(),LHC)*interSampleWeight

    #Compute the objective function within each categorical dimension
    categoricalDimInds = findall(x->typeof(x)==Categorical,dims)
    for i in categoricalDimInds
        for j = 1:dims[i].levels
            subLHC = @view LHC[LHC[:,i] .== j,:] 
            out += _AudzeEglaisObjective(dims[i],subLHC)*dims[i].weight
        end
    end

    return out
end