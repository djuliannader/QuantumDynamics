module body
push!(LOAD_PATH, pwd())
using LinearAlgebra
export diagonalization
import potential
import norm

function diagonalization(l,n,k)
     h=2*l/n
     hbar=10/10
     d=[-2/(h*h) for i in 1:(n-1)]
     du=[1/(h*h) for i in 1:(n-2)]
     A1=Array(Tridiagonal(du, d, du))
     da2=[potential.V(-l+i*h) for i in 1:(n-1)]
     A2=Array(Diagonal(da2))
     H=(-hbar^2/2)*A1+A2
     HV=eigvecs(H)
     wf=[HV[i,k] for i in 1:(n-1)]
     ev=eigvals(H)
     return [ev,wf]
end

end