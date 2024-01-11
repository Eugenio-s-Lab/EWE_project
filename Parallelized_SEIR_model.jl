using Distributions
using Random
using LinearAlgebra
using StatsBase
using Base.Threads
using DelimitedFiles
using PyCall

function p_from_rate(x::Array,delta_t)
    # Compute the probability of transition from the rate in the case of the multi-patch model
    
    return -exp.(-x .* delta_t) .+ 1
end

function p_from_rate_single_patch(x::Float64,delta_t)
    # Compute the probability of transition from the rate in the case of the single-patch model
    
    return -exp(-x * delta_t) + 1
end

function force_infection(map::Array, I::Array, beta::Float64)
    # Compute the force of infection for every tile of the multi-patch model
    
    return map * I * beta
end

 function SEIRS_multi_patch!(E::Array,I::Array,R::Array,T,delta_t,N::Array,parameters, map::Array)
#= T is the total observation time (days)
   delta_t is the time-step (days)
   N is the array of total population in study in each tile
   parameters is a dictionary of values for the parameters defining the mathematical process 
   map is the colocation map used =#

    p_mu=p_from_rate_single_patch(parameters["mu"],delta_t)
    p_alpha=p_from_rate_single_patch(parameters["alpha"],delta_t)
    p_omega=p_from_rate_single_patch(parameters["omega"],delta_t)
    
    for t = range(1,convert(Int32,round(T/delta_t))-1)
        p_beta=p_from_rate(force_infection(map[1+convert(Int32,floor((t-1)/convert(Int32,round(7/dt)))),:,:].*transpose(N),I[t,:]./N,parameters["beta"]),delta_t)
	    for i = range(1,size(E,2))
            E_new=rand(Binomial(N[i] - E[t,i] - I[t,i] - R[t,i],p_beta[i]))
            I_new=rand(Binomial(E[t,i],p_alpha))
            R_new=rand(Binomial(I[t,i],p_mu))
            S_new=rand(Binomial(R[t,i],p_omega))
            E[t+1,i]=E[t,i]+E_new-I_new
            I[t+1,i]=I[t,i]+I_new-R_new
            R[t+1,i]=R[t,i]+R_new-S_new
        end
        #if ((sum(I[t+1,:])+sum(E[t+1,:]))==0)
         #   break
        #end
    end
end

#function SEIRS_ODE_multi!(du,u,params,t)
    #= ODE version of the multi-patch stochastic model:
       p[1] is the beta rate (S-->E),
       p[2] is the alpha rate (E-->I),
       p[3] is the mu rate (I-->R),
       p[4] is the omega rate (R-->S),
       N is the total population,
       map is the colocation map between tiles =#
    
    #N,map,p=params
    #force_of_infection=force_infection(map.*N,u[2,:]./N,p[1])    
    #Tiles=size(map,1)
    #for i = range(1,Tiles) 
    #    S=N[i]-u[1,i]-u[2,i]-u[3,i]
    #    du[1,i] = force_of_infection[i]*S - p[2]*u[1,i]
    #    du[2,i] = p[2]*u[1,i] - p[3]*u[2,i]
    #    du[3,i] = p[3]*u[2,i] - p[4]*u[3,i]
    #end
#end 

starting=time()

Iterations_multi=parse(Int32,ARGS[1])

include(ARGS[2])
address=address_colocation
dir=readdir(address)
coloc_3d=zeros(Float64,size(dir,1)-1,110,110)
coloc_baseline=zeros(Float64,110,110)
for i in 1:size(dir,1)
    if (endswith(dir[i],"_baseline.txt"))
        coloc=readdlm("$address/$(dir[i])", ',', '\n')
        coloc=replace(coloc, NaN => 0)
        coloc_baseline.=coloc
    elseif (endswith(dir[i],".txt"))
        coloc = readdlm("$address/$(dir[i])", ',', '\n')
        coloc=replace(coloc, NaN => 0)
        coloc_3d[i,:,:].=coloc
    end
end
N=Int32.(readdlm(address_population, ',', '\n'))
 
open_file = open(dictionary_address, "r")
dict_str = read(open_file, String)
close(open_file)

# Parse the dictionary string using Python's ast module
python = pyimport("ast")
dict_data = python.literal_eval(dict_str)


println("Reading part successful")

params = include("params.jl")
index_seeding=dict_data[params["Seeding_department"]]
time_seeding=1+(Date(params["time_seeding"])-Date("2023-03-27")).value
params["R_local"]=parse(Float64,ARGS[3])
params["beta"]=(params["R_local"]*params["mu"])/(coloc_baseline[index_seeding,index_seeding]*N[index_seeding])
const T=size(coloc_3d,1)*7
const dt=1
const Tiles = size(coloc_3d,2)
const T_size = convert(Int32, round(T/dt))

@threads for i = 1:Iterations_multi
    R = zeros(Int32, T_size, Tiles)
    E = zeros(Int32, T_size, Tiles)
    I = zeros(Int32, T_size, Tiles)
    I[time_seeding,index_seeding]=1
    SEIRS_multi_patch!(E, I, R, T, dt, N, params, coloc_3d)
    integer_part=convert(Int32,floor(params["R_local"]))
    decimal_part=10*round(params["R_local"]-floor(params["R_local"]),digits=2)
    decimal_part_10=convert(Int32,floor(decimal_part))
    decimal_part_100=convert(Int32,10*round((decimal_part-floor(decimal_part)),digits=1))
    file_path = "$output_address/I+E_output_$(i)_R_local_$(integer_part)_$(decimal_part_10)"*"$(decimal_part_100).txt"
    file_path_1 = "$output_address/R_output_$(i)_R_local_$(integer_part)_$(decimal_part_10)"*"$(decimal_part_100).txt"
    writedlm(file_path, I.+E , ',')
    writedlm(file_path_1, R , ',')
    println("Writing part $i successful with R local $(params["R_local"])")
end

ending=time()
println(ending-starting)