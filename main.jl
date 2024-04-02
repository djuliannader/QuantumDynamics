
push!(LOAD_PATH, pwd())
import body
import norm
import wigner
import tunneling
using Plots


println("\r Quantum-Dynamics 1D ")
println("\r Finite Difference Method")


# Reading data from input file
#------------------------------------
open("input.dat") do f
 K1=readline(f)
 K2=readline(f)
 L = parse(Float64, K2)
 K3=readline(f)
 K4=readline(f)
 N = parse(Int64, K4)
 K5=readline(f)
 K6=readline(f)
 k = parse(Int64, K6)
 K7=readline(f)
 K8=readline(f)
 K9=readline(f)
 K10=readline(f)


 tt=r"([0-9])" 
 tpl = [parse(Int64,t.match) for t in eachmatch(tt, K10)]
 
#------------------------------------

# printing information 
println("------------------------------------------------")
println("Space from -L to L with L=",L)
println("Partitioning the space in N=",N, " subintervals")

@time begin

#------------------------------------------------------------------------------------
#------ Calculations of the stationary state starts ---------------------------------
#------------------------------------------------------------------------------------


# calling function which performs the diagonalizacion
r=body.diagonalization(L,N,k)


# printing results
# calling function which normalize the wave function and priting results
 println("-------------    results of the stationary (k=",k,") state     ---------------")
 wf=norm.normalizing(r[2],2L/N)
 x=[-L+(2L/N)*i for  i in 1:(N-1)]
 open("wavefunction_stationary.dat","w") do io
     for i in 1:length(x)
       println(io,x[i]," ",round(real(wf[i]),digits=16))
     end
 end
 println("See file wavefunction_stationary.dat for the wave function ")
 # calling routine to calculate wigner function
  wig=wigner.wignerf(wf,L,N,"wigner_stationary.dat")
 println("Ground state energy = ",r[1][1])
 println("Energy of the ",k," state = ",r[1][k])
 # calling function which calculate the derivative of the wave function at the center of coordinates
 der=norm.derivative2(wf,2L/N)
 println("Second derivative at the center of coordinates: ",der)


# calling function which calculate the classical turning points
 tpoints=tunneling.turnningpoints(wf,L,N,r[1][k])
 println("Classical turning points: ",tpoints)

# calling function which transmision coefficient
 if K10=="True"
   tcoef=tunneling.wkbt(wf,tpoints[tpl[1]],tpoints[tpl[2]],L,N,r[1][k])
   println("WKB transmission coeficient: ",tcoef)
 end

println("----------------------------------------------------------------")


#------------------------------------------------------------------------------------
#------ Calculations of the stationary state ends   ---------------------------------
#------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------
#------ Calculations of the dynamics of coherent state starts -----------------------
#------------------------------------------------------------------------------------


end


end





