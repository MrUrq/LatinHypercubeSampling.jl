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
    "text": "LatinHypercubeSampling is a Julia package for the creation of optimised Latin Hypercube (LHC) Sampling Plans. The genetic optimisation algorithm is largely based on the work by Bates et al. [1]. The package includes additional functionality for the creation of an optimised subset of an existing plan."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The package is registered in METADATA.jl and can be installed with Pkg.add.julia> Pkg.add(\"LatinHypercubeSampling\")"
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
    "text": "function LHCoptim(n::Int,d::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8)\n\nProduce an optimized Latin Hyper Cube with d dimensions and n sample points. Optimization is run for gens generations.\n\n\n\n"
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
    "text": "The LHCoptim function run for many generations to create an optimised 120 point plan in 2 dimensions.julia> LHCoptim(120,2,gens)(Image: alt text)"
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
    "text": "function subLHCoptim(X,n::Int,gens;popsize::Int=100,ntour::Int=2,ptour=0.8)\n\nProduce an optimized Latin Hyper Cube with n sample points from a subset of points in X. Optimization is run for gens generations.\n\n\n\n"
},

{
    "location": "man/sublhcoptim.html#LatinHypercubeSampling.subLHCindex-Tuple{Any,Any}",
    "page": "Optimised subset of LHC Sampling Plan",
    "title": "LatinHypercubeSampling.subLHCindex",
    "category": "method",
    "text": "function subLHCindex(X,Xsub)\n\nIndex in the large LHC to get the subLHC.\n\n\n\n"
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
    "text": "Create an optimised subset LHC plan from an existing plan with 120 points in 2 dimensions.(Image: alt text)julia> subLHCoptim(X,Xsub)(Image: alt text)The indices of the subset in the larger plan can be extracted withjulia> subLHCindex(X,Xsub)"
},

{
    "location": "man/refinelhc.html#",
    "page": "Refine LHC Sampling Plan",
    "title": "Refine LHC Sampling Plan",
    "category": "page",
    "text": ""
},

{
    "location": "man/refinelhc.html#Refine-LHC-Sampling-Plan-1",
    "page": "Refine LHC Sampling Plan",
    "title": "Refine LHC Sampling Plan",
    "category": "section",
    "text": "Refine existing plan through mutation only. Can be useful for fast improvement to existing plan.note: Note\nLikely to converge to local optima."
},

{
    "location": "man/refinelhc.html#LatinHypercubeSampling.refineLHCoptim-Tuple{Any,Any}",
    "page": "Refine LHC Sampling Plan",
    "title": "LatinHypercubeSampling.refineLHCoptim",
    "category": "method",
    "text": "function refineLHCoptim(X,gens;popsize::Int=100)\n\nRefine an existing plan by mutation only.\n\n\n\n"
},

{
    "location": "man/refinelhc.html#Function-1",
    "page": "Refine LHC Sampling Plan",
    "title": "Function",
    "category": "section",
    "text": "refineLHCoptim(X,gens;popsize::Int=100)"
},

]}
