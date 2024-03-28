push!(LOAD_PATH, pwd())
module tunneling
export turnningpoints
export wkbt
import potential



function turnningpoints(u::Vector{Float64},LL,NN,Ener)
         b=0
         d=2*LL/NN
         x=[-LL+i*d for i in 1:(NN-1)]
         difaa=abs(potential.V(x[1])+b*u[1]*u[1]-Ener)
	 difa=abs(potential.V(x[2])+b*u[2]*u[2]-Ener)
	 tplist=Vector{Float64}()
	 for i in 3:(NN-1)
	     dif=abs(potential.V(x[i])+b*u[i]*u[i]-Ener)
	     if (difa<dif) && (difa<difaa)
	      	tp=x[i-1]
		append!(tplist,tp)
	     end
	    difaa=difa
	    difa=dif
	 end
	 return tplist
	 end



function wkbt(u::Vector{Float64},tp1,tp2,LL,NN,Ener)
	d=2*LL/NN
	b=0
	ni=trunc(Int,((tp1-(-LL))/d))
	nf=trunc(Int,((tp2-(-LL))/d))
	x=[-LL+i*d for i in 1:(NN-1)]
	t1=[((2*abs(Ener-(potential.V(x[i])+b*u[i]*u[i])))^(1/2))*d for i in ni:nf]
        suma=sum(t1)
	T=exp(-2*suma)
        return T
	end


end

