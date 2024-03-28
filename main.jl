
push!(LOAD_PATH, pwd())
import body
import norm
import wigner
import tunneling
using Plots


println("\r Gross-Pitaeskii 1D ")
println("\r Selfconsistent method")


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
 K11=readline(f)
 K12=readline(f)


 tt=r"([0-9])" 
 tpl = [parse(Int64,t.match) for t in eachmatch(tt, K12)]
 
#------------------------------------

# printing information 
println("------------------------------------------------")
println("Space from -L to L with L=",L)
println("Partitioning the space in N=",N, " subintervals")
println("------------------------------------------------")

# calling function which performs selfconsistent method
@time begin
r=body.diagonalization(L,N,k)
end

# printing results
# calling function which normalize the wave function
 wf=norm.normalizing(r[2],2L/N)
 println("-----------------------------------------------")
 println("Ground state energy = ",r[1][1])
 # calling function which calculate the derivative of the wave function at the center of coordinates
 der=norm.derivative2(wf,2L/N)
 println("Second derivative at the center of coordinates: ",der)


# calling function which calculate the classical turning points
 tpoints=tunneling.turnningpoints(wf,L,N,r[1][1])
 println("Classical turning points: ",tpoints)

# calling function which transmision coefficient
 if K10=="True"
   tcoef=tunneling.wkbt(wf,tpoints[tpl[1]],tpoints[tpl[2]],L,N,r[1][1])
   println("WKB transmission coeficient: ",tcoef)
 end

# calling routine to calculate wigner function
  wig=wigner.wignerf(wf,L,N)


# printing the wave function
if K12=="True"
   println("-------------")
   println("u:",wf)
   println("-------------")
 end



x=[-L+(2L/N)*i for  i in 1:(N-1)]
plot(x,wf,title="wave function")
 xlabel!("x")
 ylabel!("y(x)")
end



