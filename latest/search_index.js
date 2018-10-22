var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#LatinHypercubeSampling-1",
    "page": "Home",
    "title": "LatinHypercubeSampling",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "LatinHypercubeSampling is a Julia package for the creation of optimised Latin Hypercube (LHC) Sampling Plans. The genetic optimisation algorithm is largely based on the work by Bates et al. [1]. The package includes additional functionality for the creation of an optimised subset of an existing plan, as well  as the inclusion of discrete parameters."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The package is registered and can be installed with Pkg.add.julia> Pkg.add(\"LatinHypercubeSampling\")"
},

{
    "location": "index.html#Reference-1",
    "page": "Home",
    "title": "Reference",
    "category": "section",
    "text": "[1]: Stuart Bates, Johann Sienz, and Vassili Toropov. \"Formulation of the Optimal Latin Hypercube Design of Experiments Using a Permutation Genetic Algorithm\", 45th AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics & Materials Conference, Structures, Structural Dynamics, and Materials and Co-located Conferences, () https://doi.org/10.2514/6.2004-2011"
},

{
    "location": "man/lhcoptim.html#",
    "page": "Optimised Latin Hypercube Sampling Plan",
    "title": "Optimised Latin Hypercube Sampling Plan",
    "category": "page",
    "text": ""
},

{
    "location": "man/lhcoptim.html#Optimised-Latin-Hypercube-Sampling-Plan-1",
    "page": "Optimised Latin Hypercube Sampling Plan",
    "title": "Optimised Latin Hypercube Sampling Plan",
    "category": "section",
    "text": "Create an optimised Latin Hypercube Sampling Plan using a genetic based optimisation algorithm. The objective function is the inverse of the Audze-Eglais function defined astextmax  U = textmax sum_p=1^P sum_q=p+1^P L^2_pqwhere L^2_pq is the square of the Euclidean norm.note: Note\nThis package maximises the inverse of the Audze-Eglais objective function."
},

{
    "location": "man/lhcoptim.html#LatinHypercubeSampling.LHCoptim-Tuple{Int64,Int64,Any}",
    "page": "Optimised Latin Hypercube Sampling Plan",
    "title": "LatinHypercubeSampling.LHCoptim",
    "category": "method",
    "text": "function LHCoptim(n::Int,d::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8)\n\nProduce an optimized Latin Hyper Cube with d dimensions and n sample points. Optimization is run for gens generations.\n\n\n\n\n\n"
},

{
    "location": "man/lhcoptim.html#Function-1",
    "page": "Optimised Latin Hypercube Sampling Plan",
    "title": "Function",
    "category": "section",
    "text": "LHCoptim(n::Int,d::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8)"
},

{
    "location": "man/lhcoptim.html#Example-1",
    "page": "Optimised Latin Hypercube Sampling Plan",
    "title": "Example",
    "category": "section",
    "text": "The LHCoptim function run for many generations to create an optimised 120 point plan in 2 dimensions.julia> LHCoptim(120,2,gens)(Image: )"
},

{
    "location": "man/sublhcoptim.html#",
    "page": "Optimised subset of LHC Sampling Plan",
    "title": "Optimised subset of LHC Sampling Plan",
    "category": "page",
    "text": ""
},

{
    "location": "man/sublhcoptim.html#Optimised-subset-of-LHC-Sampling-Plan-1",
    "page": "Optimised subset of LHC Sampling Plan",
    "title": "Optimised subset of LHC Sampling Plan",
    "category": "section",
    "text": "Generate an optimised subset of an existing plan. The optimisation of the subset is based on a genetic algorithm."
},

{
    "location": "man/sublhcoptim.html#LatinHypercubeSampling.subLHCoptim-Tuple{Int64,Int64,Any}",
    "page": "Optimised subset of LHC Sampling Plan",
    "title": "LatinHypercubeSampling.subLHCoptim",
    "category": "method",
    "text": "function subLHCoptim(X,n::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8)\n\nProduce an optimized Latin Hyper Cube with n sample points from a subset of points in X. Optimization is run for gens generations.\n\n\n\n\n\n"
},

{
    "location": "man/sublhcoptim.html#LatinHypercubeSampling.subLHCindex-Tuple{Any,Any}",
    "page": "Optimised subset of LHC Sampling Plan",
    "title": "LatinHypercubeSampling.subLHCindex",
    "category": "method",
    "text": "function subLHCindex(X,Xsub)\n\nIndex in the large LHC to get the subLHC.\n\n\n\n\n\n"
},

{
    "location": "man/sublhcoptim.html#Functions-1",
    "page": "Optimised subset of LHC Sampling Plan",
    "title": "Functions",
    "category": "section",
    "text": "subLHCoptim(n::Int,d::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8)subLHCindex(X,Xsub)"
},

{
    "location": "man/sublhcoptim.html#Example-1",
    "page": "Optimised subset of LHC Sampling Plan",
    "title": "Example",
    "category": "section",
    "text": "Create an optimised subset LHC plan from an existing plan with 120 points in 2 dimensions.(Image: )julia> subLHCoptim(X,Xsub)(Image: )The indices of the subset in the larger plan can be extracted withjulia> subLHCindex(X,Xsub)"
},

{
    "location": "man/categorical.html#",
    "page": "Categorical Latin Hypercube Sampling Plan",
    "title": "Categorical Latin Hypercube Sampling Plan",
    "category": "page",
    "text": ""
},

{
    "location": "man/categorical.html#Categorical-Latin-Hypercube-Sampling-Plan-1",
    "page": "Categorical Latin Hypercube Sampling Plan",
    "title": "Categorical Latin Hypercube Sampling Plan",
    "category": "section",
    "text": "Categorical Latin Hypercube plans allows one to mix discrete and continous data  in the same plan. "
},

{
    "location": "man/categorical.html#Example-1",
    "page": "Categorical Latin Hypercube Sampling Plan",
    "title": "Example",
    "category": "section",
    "text": "Say we have two continous dimensions as well as one on/off, discrete, dimension.  These can be included in the same sampling plan withjulia> numPoints = 100\njulia> dims = [Continous(),Continous(),Categorical(2)]\njulia> initialSample = randomLHC(numPoints,dims)\njulia> X = LHCoptim!(initialSample,gens;dims=dims)[1]note: Note\nThis is no longer strictly a Latin Hypercube because of the introduction of  categorical values.The objective function is altered to include the separation within each plane of the  categorical values in addition to the separation between all points. The weights  for each plane can be supplied by the user. Large emphasis can be put on keeping the separation within each plane by  increasing its weight. This is similar to doing a separate LHC for each categorical  dimension. julia> weights = [1,1,1000]\njulia> julia> X = LHCoptim!(initialSample,gens;dims=dims,weights=weights)[1]using PlotlyJS, LatinHypercubeSampling # hide\nnumPoints = 100 # hide\nweights = [1,1,1000] # hide\ndims = [Continous(),Continous(),Categorical(2)]  # hide\ninitialSample = randomLHC(numPoints,dims) # hide \nX = LHCoptim!(initialSample,50;dims=dims,weights=weights)[1] # hide \n\nfunction plotlhc(X,titletext) # hide\n    x1 = X[X[:,3].==1,:] # hide\n    x2 = X[X[:,3].==2,:] # hide\n    trace1 = scatter(;x=x1[:,1], y=x1[:,2], # hide\n                        mode=\"markers\", name=\"Category 1\", # hide\n                        marker_size=12) # hide\n\n    trace2 = scatter(;x=x2[:,1], y=x2[:,2], # hide\n                        mode=\"markers\", name=\"Category 2\", # hide\n                        marker_size=12) # hide\n\n    \n    data = [trace1, trace2] # hide\n    layout = Layout(height=650, # hide\n                    width=740, # hide\n                    title=titletext, # hide\n                    xaxis=attr(title=\"Continous dim. 1\"), # hide\n                    yaxis=attr(title=\"Continous dim. 2\"), # hide\n                    margin=attr(l=100, r=30, b=50, t=90), # hide\n                                ) # hide\n    plot(data, layout) # hide\nend # hide\n\np = plotlhc(X,\"Promote in-plane separation\") # hide\npkgpath = abspath(joinpath(dirname(Base.find_package(\"LatinHypercubeSampling\")), \"..\")) # hide\nsavedir = joinpath(pkgpath,\"docs\",\"src\",\"assets\",\"example1.html\") # hide\nPlotlyJS.savehtml(p,savedir,:embed) # hide    <iframe src=\"../assets/example1.html\" height=\"765\" width=\"765\" frameborder=\"0\" seamless=\"seamless\" scrolling=\"no\"></iframe>Similarly we can turn of the separation within planes entirely with julia> weights = [1,1,0]\njulia> julia> X = LHCoptim!(initialSample,gens;dims=dims,weights=weights)[1]weights = [1,1,0] # hide\nX = LHCoptim!(initialSample,50;dims=dims,weights=weights)[1] # hide \n\n\np = plotlhc(X,\"Promote inter sample separation\") # hide\npkgpath = abspath(joinpath(dirname(Base.find_package(\"LatinHypercubeSampling\")), \"..\")) # hide\nsavedir = joinpath(pkgpath,\"docs\",\"src\",\"assets\",\"example2.html\") # hide\nPlotlyJS.savehtml(p,savedir,:embed) # hide    <iframe src=\"../assets/example2.html\" height=\"765\" width=\"765\" frameborder=\"0\" seamless=\"seamless\" scrolling=\"no\"></iframe>"
},

]}
