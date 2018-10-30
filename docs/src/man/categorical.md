# Categorical Latin Hypercube Sampling Plan

Categorical Latin Hypercube plans allows one to mix discrete and continuous data 
in the same plan. 
 
## Example
Say we have two continuous dimensions as well as one on/off, discrete, dimension. 
These can be included in the same sampling plan with

```julia-repl
julia> numPoints = 100
julia> dims = [Continuous(),Continuous(),Categorical(2)]
julia> initialSample = randomLHC(numPoints,dims)
julia> X = LHCoptim!(initialSample,gens;dims=dims)[1]
```


!!! note

    This is no longer strictly a Latin Hypercube because of the introduction of 
    categorical values.

The objective function is altered to include the separation within each plane of the 
categorical values in addition to the separation between all points. The weights
 for each plane can be supplied by the user. 

Large emphasis can be put on keeping the separation within each plane by
 increasing its weight. This is similar to doing a separate LHC for each categorical
 dimension. 
```julia-repl
julia> weights = [1,1,1000]
julia> julia> X = LHCoptim!(initialSample,gens;dims=dims,weights=weights)[1]
```

```@setup x
using PlotlyJS, LatinHypercubeSampling # hide
numPoints = 100 # hide
weights = [1.0,1.0,1000.0] # hide
dims = [Continuous(),Continuous(),Categorical(2)]  # hide
initialSample = randomLHC(numPoints,dims) # hide 
X = LHCoptim!(initialSample,50;dims=dims,weights=weights)[1] # hide 

function plotlhc(X,titletext) # hide
    x1 = X[X[:,3].==1,:] # hide
    x2 = X[X[:,3].==2,:] # hide
    trace1 = scatter(;x=x1[:,1], y=x1[:,2], # hide
                        mode="markers", name="Category 1", # hide
                        marker_size=12) # hide

    trace2 = scatter(;x=x2[:,1], y=x2[:,2], # hide
                        mode="markers", name="Category 2", # hide
                        marker_size=12) # hide

    
    data = [trace1, trace2] # hide
    layout = Layout(height=650, # hide
                    width=740, # hide
                    title=titletext, # hide
                    xaxis=attr(title="Continuous dim. 1"), # hide
                    yaxis=attr(title="Continuous dim. 2"), # hide
                    margin=attr(l=100, r=30, b=50, t=90), # hide
                                ) # hide
    plot(data, layout) # hide
end # hide

p = plotlhc(X,"Promote in-plane separation") # hide
pkgpath = abspath(joinpath(dirname(Base.find_package("LatinHypercubeSampling")), "..")) # hide
savedir = joinpath(pkgpath,"docs","src","assets","example1.html") # hide
PlotlyJS.savehtml(p,savedir,:embed) # hide
```
```@raw html
    <iframe src="../assets/example1.html" height="765" width="765" frameborder="0" seamless="seamless" scrolling="no"></iframe>
```



Similarly we can turn of the separation within planes entirely with 
```julia-repl
julia> weights = [1.0,1.0,0.0]
julia> julia> X = LHCoptim!(initialSample,gens;dims=dims,weights=weights)[1]
```
```@setup x
weights = [1.0,1.0,0.0] # hide
X = LHCoptim!(initialSample,50;dims=dims,weights=weights)[1] # hide 


p = plotlhc(X,"Promote inter sample separation") # hide
pkgpath = abspath(joinpath(dirname(Base.find_package("LatinHypercubeSampling")), "..")) # hide
savedir = joinpath(pkgpath,"docs","src","assets","example2.html") # hide
PlotlyJS.savehtml(p,savedir,:embed) # hide
```
```@raw html
    <iframe src="../assets/example2.html" height="765" width="765" frameborder="0" seamless="seamless" scrolling="no"></iframe>
```