### A Pluto.jl notebook ###
# v0.14.5

using LinearAlgebra, DifferentialEquations, Plots, DataFrames, CSV

# Rossler System

global lyapunov_list_x = []
global lyapunov_list_y = []
global lyapunov_list_z = []

global lya_sum = 0

@inline @inbounds function loop(du,u,p,t)
	# Defferential Equations

	x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz = u
	r = sqrt(xx^2 + yy^2 +zz^2)
	append!(lyapunov_list_x,log((abs(xx+yy+zz)))/t)
	append!(lyapunov_list_y,log(abs((yx+yy+yz))/t))
	append!(lyapunov_list_z,log(abs(((zx+zy+zz)/1))/(t)))
	xr = sqrt(xx^2 + yx^2 + zx^2)
	yr = sqrt(yx^2 + yy^2 + yz^2)
	zr = sqrt(zx^2 + zy^2 + zz^2)

	xx = xx/xr
	xy = xy/xr
	xz = xz/xr


	yx = xx/yr
	yy = xy/yr
	yz = xz/yr


	zx = zx/xr
	zy = zy/xr
	zz = zz/xr

	a,b,c = p

	du[1] = -y -z
	du[2] = x +a*y
	du[3] = b + z*(x-c)

	# Variational Equations
	du[4] = -xy-xz
	du[5] = xx + a*xy
	du[6] = z*xx + (x-c)*xy
	du[7] = -yy -yz
	du[8] = yx + a*yy
	du[9] = z*yx + (x-c)*yz
	du[10] = -zy -zz
	du[11] = zx + a*zy
	du[12] = z*zx + (x-c)*zz
end

a= 0.1
b= 0.1
c= 18.0
p = [a,b,c]
u0 = [0,1,2,1,0,0,0,1,0,0,0,1]
tspan =(0.0,850)

prob = ODEProblem(loop,u0,tspan,p)

solution = solve(prob,Feagin10(), reltol=1e-10, abstol=1e-10)

df = DataFrame(solution)

p0 = scatter(solution,vars=(1,2),
     xaxis="x",yaxis="y",zaxis="z")

p1 = scatter(solution,vars=(1,2,3),
     xaxis="x",yaxis="y",zaxis="z")
#p2 = plot(lyapunov_list_x,xlim=(10000,20000),ylim=(-2,2))
p2 = plot(lyapunov_list_x,ylim=(-0.025,0.025))
p3 = plot(lyapunov_list_y)
p4 = plot(lyapunov_list_z)



CSV.write("df.csv",df)
#print(lyapunov_list)
print(length(lyapunov_list_x))
plot(p0,p1,p2,p3,p4)
