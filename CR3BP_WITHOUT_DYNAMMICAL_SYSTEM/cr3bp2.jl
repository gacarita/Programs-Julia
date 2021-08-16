using DifferentialEquations, Plots, Symbolics,LinearAlgebra


@inline @inbounds function vy0(C,p,v)
   """
   Calculate Vy0 from jacobian value
   input : (C,p,v) ...
            C = Jacobi cte
            v = x,y,vx
            p = mass ratio parameters mu
   """
       #cartesian coordinates
       x, y, x̂= v

       #distance vectors
       r1 = ((x+μ)^2 + y^2)^(1/2)
       r2 = ((1-x-μ)^2 + y^2)^(1/2)
       #jacobi integral - constant C.
       ŷ = (x^2 + y^2 - x̂^2 + 2*((1-μ)/r1 + μ/r2) - C)^(1/2)
       return ŷ
end

# Planar Circular Restricted 3body problem
# Equations of motion in co-rotating frame:

global μ = 0.1

# Initial Conditions

x0 = 1.983
y0 = 0
vx0 = 0
vyu = -vy0(0,μ,[x0,y0,vx0]) #sign negative for retrograde orbits
u0 = [x0, 0, 0,vyu] # vector with initial conditions

# Time Span
t_end = 2*pi*5e+2
tspan=(0.0,t_end)

global eigenvaluee = 0

@inline @inbounds function cr3bp_rhs(du,u,μ,t)
   #cartesian coordinates
   x, y, x̂, ŷ = u

   #distance vectors
   r1 = ((x + μ)^2 + y^2)^(1/2)
   r2 = ((-1 + x + μ)^2 + y^2)^(1/2)

   #equations of motion
   du[1] = x̂ # F1
   du[2] = ŷ # F2
   du[3] =  2ŷ  + x - (1 - μ)*(μ + x)/(r1^3) - μ*(-1 + μ + x)/(r2^3) # F3
   du[4] = -2x̂  + y - (1 - μ)*y/(r1^3)       - μ*y/(r2^3) # F4
   eigenvalues = f_jac(u,μ,t)
   global eigenvaluee += log((eigenvalues[1]))
end


global J=zeros((4,4))


@inline @inbounds function f_jac(u,μ,t)
  x, y, x̂, ŷ = u

  #J=zeros(Float64, (4,4))

  J[1,1] = 0
  J[1,2] = 0
  J[1,3] = 1
  J[1,4] = 0

  J[2,1] = 0
  J[2,2] = 0
  J[2,3] = 0
  J[2,4] = 1

  J[3,1] = 1 + (μ - 1)*(((y^2 + (x + μ)^2)^1.5)^-1) + 1.5μ*(x + μ - 1)*(2x + 2μ - 2)*((y^2 + (x + μ - 1)^2)^0.5)*(((y^2 + (x + μ - 1)^2)^1.5)^-2) - (μ*(((y^2 + (x + μ - 1)^2)^1.5)^-1)) - (1.5(-x - μ)*(2x + 2μ)*(1 - μ)*((y^2 + (x + μ)^2)^0.5)*(((y^2 + (x + μ)^2)^1.5)^-2))
  J[3,2] = 3.0y*μ*(x + μ - 1)*((y^2 + (x + μ - 1)^2)^0.5)*(((y^2 + (x + μ - 1)^2)^1.5)^-2) - (3.0y*(-x - μ)*(1 - μ)*((y^2 + (x + μ)^2)^0.5)*(((y^2 + (x + μ)^2)^1.5)^-2))
  J[3,3] = 0
  J[3,4] = 2

  J[4,1] = 1.5y*μ*(2x + 2μ - 2)*((y^2 + (x + μ - 1)^2)^0.5)*(((y^2 + (x + μ - 1)^2)^1.5)^-2) + 1.5y*(2x + 2μ)*(1 - μ)*((y^2 + (x + μ)^2)^0.5)*(((y^2 + (x + μ)^2)^1.5)^-2)
  J[4,2] = 1 + (μ - 1)*(((y^2 + (x + μ)^2)^1.5)^-1) + 3.0μ*(y^2)*((y^2 + (x + μ - 1)^2)^0.5)*(((y^2 + (x + μ - 1)^2)^1.5)^-2) + 3.0(y^2)*(1 - μ)*((y^2 + (x + μ)^2)^0.5)*(((y^2 + (x + μ)^2)^1.5)^-2) - (μ*(((y^2 + (x + μ - 1)^2)^1.5)^-1))
  J[4,3] = -2
  J[4,4] = 0
  JJ = eigvals(J)
  return JJ
  #nothing
end


# solver

prob = ODEProblem(cr3bp_rhs,u0,tspan,μ)
solution2 = solve(prob,Feagin10(), reltol=1e-12, abstol=1e-12,saveat=2*pi)

p2 = scatter(solution2,vars=(1,2),title="integrator Feagin10()",
     xaxis="x (AU)",yaxis="y (AU)")
p1 = scatter(solution2,vars=(3,4),title="integrator Feagin10()",
          xaxis="x (AU)",yaxis="y (AU)")

plot(p2,p1)
#println(eigenvalue/t_end)
println(eigenvalue/length(solution2[1]))
