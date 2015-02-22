# Programming Assignment 01
# James Folberth
# Spring 2015

# the sparse module is included in Base.  It uses CSC

# You can run PA1Driver by doing
# julia> include("prog01.jl")
# julia> PA1Driver(Lx, Ly, sigma)
#
# For example, you could do
# julia> PA1Driver(2,3,1.)
# 
# Be sure that the third argument to PA1Driver is a float by appending
# a period, or doing 1e0

# Note that I'm running Julia 0.4
# Here are the specifics.  I don't know if v0.3 is vastly different from v0.4
# julia> versioninfo()
# Julia Version 0.4.0-dev+2677
# Commit 42e45a9 (2015-01-13 17:14 UTC)
# Platform Info:
#   System: Linux (x86_64-unknown-linux-gnu)
#   CPU: Intel(R) Core(TM) i5-2467M CPU @ 1.60GHz
#   WORD_SIZE: 64
#   BLAS: libmkl_rt
#   LAPACK: libmkl_rt
#   LIBM: libimf
#   LLVM: libLLVM-3.5.0

function ModelProblem1D(L::Int, sigma::Float64 = 1e0) 
   m = 2^L-1; n = m;
   h = 1e0/float64(m+1)

   # Build matrix in COO format and then pass to Julia's sparse to form CSC

   row = Array(Int, 3*n-2)
   col = Array(Int, 3*n-2)
   val = Array(Float64, 3*n-2)

   row[1:2] = [1,1]
   col[1:2] = [1,2]
   val[1:2] = [2e0/h^2+sigma, -1e0/h^2]

   for i in 2:n-1
      row[3*i-3:3*i-1] = [i,i,i]
      col[3*i-3:3*i-1] = [i-1,i,i+1]
      val[3*i-3:3*i-1] = [-1e0/h^2, 2e0/h^2+sigma, -1e0/h^2]
   end

   row[3*n-3:3*n-2] = [n,n]
   col[3*n-3:3*n-2] = [n-1,n]
   val[3*n-3:3*n-2] = [-1e0/h^2, 2e0/h^2+sigma]

   # convert to built-in CSC format
   return sparse(row, col, val, m, n)
end

function ModelProblem2D(Lx::Int, Ly::Int, sigma::Float64 = 1e0)
   m = 2^Lx-1; n = 2^Ly-1;
   Ax = ModelProblem1D(Lx,0e0)
   Ix = speye(Float64, m)
   Ay = ModelProblem1D(Ly,0e0)
   Iy = speye(Float64, n)

   # kron respects sparse structure
   #return kron(Ax,Iy)+kron(Ix,Ay)+sigma*speye(Float64,m*n) # constant x
   return kron(Iy,Ax)+kron(Ay,Ix)+sigma*speye(Float64,m*n) # constant y
end

function SparseMatVec(A::SparseMatrixCSC,x::Array{Float64,1})
   # CSC mat-vec
   A.n == length(x) || error("dimension mismatch")
   b = zeros(Float64,A.m) # need to zero prior to loop, unlike CSR
   for c in 1:A.n # length(A.colptr) == size(A,2)
      for ri in A.colptr[c]:(A.colptr[c+1]-1)
         @inbounds b[A.rowval[ri]] += A.nzval[ri]*x[c] # stop bounds checking
      end
   end
   return b
end

function PA1Driver(Lx::Int, Ly::Int, sigma::Float64)
# You can run PA1Driver by doing
# julia> include("prog01.jl")
# julia> PA1Driver(Lx, Ly, sigma)
#
# For example, you could do
# julia> PA1Driver(2,3,1.)

   A1 = ModelProblem1D(Lx,sigma)
   A2 = ModelProblem2D(Lx,Ly,sigma)

   x = ones(Float64,size(A2,2))
   b2 = SparseMatVec(A2,x)

   println("Rows of matrix files are column pointer, row index, and value")

   # write to files
   fA1 = open("A1.txt","w+")
   if isopen(fA1)
      print(fA1,A1.colptr,"\n\n")
      print(fA1,A1.rowval,"\n\n")
      print(fA1,A1.nzval)
      close(fA1)
      println("1D model problem A matrix written to 'A1.txt'")
   end

   fA2 = open("A2.txt","w+")
   if isopen(fA2)
      print(fA2,A2.colptr,"\n\n")
      print(fA2,A2.rowval,"\n\n")
      print(fA2,A2.nzval,"\n\n")
      close(fA2)
      println("2D model problem A matrix written to 'A2.txt'")
   end

   fb2 = open("b2.txt","w+")
   if isopen(fb2)
      print(fb2,b2)
      close(fb2)
      println("SparseMatVec test output written to 'b2.txt'")
   end

end

# These are just some routines to test wheter or not I goofed something up
function test()
   @time A = ModelProblem2D(10,12,1e0)
   #x = ones(Float64,size(A,2))

   #@time b = SparseMatVec(A,x)
   #@time b_control = A*x # use Julia's matvec
   #println("||b-A*x||_Inf = ",norm(b-b_control,Inf))

   Lx = 8
   A = ModelProblem1D(Lx,1e0)
   x = linspace(0,1,2^Lx+1)
   f = exp(-x).*(pi^2*sin(pi*x)+2.*pi*cos(pi*x))
   f[1]=0.; f[end]=0.
   v = A\f[2:end-1]
   v = [0, v, 0]

   plot(x,sin(pi*x).*exp(-x)-v)

   #Lx = 6
   #Ly = 8
   #A = ModelProblem2D(Lx,Ly,1e0)
   #x = linspace(0,1,2^Lx+1)
   #y = linspace(0,1,2^Ly+1)

   #u = Array(Float64, (length(x)-2, length(y)-2))
   #f = Array(Float64, (length(x)-2, length(y)-2))
   ##@time begin
   #for j = 2:length(y)-1
   #   for i = 2:length(x)-1
   #      # XXX @simd can't inline sin/cos/exp, so won't get SIMD to work
   #      @inbounds u[i-1,j-1] = sin(pi*x[i])*sin(pi*y[j])*exp(-x[i]^2+y[j]^2)
   #      @inbounds f[i-1,j-1] = exp(-x[i]^2+y[j]^2)*(
   #         -4.*pi*y[j]*cos(pi*y[j])*sin(pi*x[i]) + (
   #         4*pi*x[i]*cos(pi*x[i]) + (1+2*pi^2-4*x[i]^2-4*y[j]^2)*sin(pi*x[i]))
   #         *sin(pi*y[j]))
   #   end
   #end
   ##end
   #  
   #v = A\reshape(f,(length(f),1)) # constant y
   #v = reshape(v, (length(x)-2, length(y)-2))

   #error = u-v
   #display(maxabs(error[:]))

   return
end

test()
#PA1Driver(2,2,1e0)
