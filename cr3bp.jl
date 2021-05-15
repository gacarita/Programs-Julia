### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 1b497e02-b4df-11eb-2dcd-4f4ac0a22336
using DynamicalSystems, Plots # also exports relevant StaticArrays names

# ╔═╡ 4cf2fb21-30c3-4a15-8705-4217047f3c7a
begin
	# Planar Circular Restricted 3body problem
	# Equations of motion:
	
	@inline @inbounds function loop(u,p,t)
		#cartesian coordinates
		x, y, x̂, ŷ = u
		μ = p
		#distance vectors	
		r1 = ((x + μ)^2 + y^2)^(1/2)
		r2 = ((-1 + x + μ)^2 + y^2)^(1/2)
		
		#equations of motion
		du1 = x̂ # F1
		du2 = ŷ # F2
		du3 =  2ŷ  + x - (1 - μ)*(μ + x)/(r1^3) - μ*(-1 + μ + x)/(r2^3) # F3
		du4 = -2x̂  + y - (1 - μ)*y/(r1^3)       - μ*y/(r2^3) # F4
		return SVector{4}(du1,du2,du3,du4)
	end
	
	@inline @inbounds function loop_jac(u,p,t)
		#cartesian coordinates
		x, y, x̂, ŷ = u
		μ = p
		#distance vectors	
		r1 = ((x+μ)^2 + y^2)^(1/2)
		r2 = ((1-x-μ)^2 + y^2)^(1/2)
		dr1x = (x+μ)/r1
		dr2x = (1-x-μ)/r2
		dr1y = y/r1
		dr2y = y/r2
		
		JF3x = 1 -x*dr1x - r1 - x*dr2x - r2 
		JF4y = 1- y*dr1y - r1 - y*dr2y - r2
		
		#Jacobian Matrix
		
		J = @SMatrix [0  0  1  0;
		0  0  0  1;
		JF3x  0  0  2;
		0  JF4y -2  0]
		
		return J
	end
end

# ╔═╡ 1a300c74-bdc0-4953-85f6-e43055028ffc
@inline @inbounds function vy0(C,p,v)
		#cartesian coordinates
		x, y, x̂= v
		μ = p
		#distance vectors	
		r1 = ((x+μ)^2 + y^2)^(1/2)
		r2 = ((1-x-μ)^2 + y^2)^(1/2)
		
		ŷ = (x^2 + y^2 - x̂^2 + 2*((1-μ)/r1 + μ/r2) - C)^(1/2)
		
		return ŷ
	end

# ╔═╡ 3702611c-91e7-4277-8a3a-f7eb0862892c
begin
	x0 = 1.983
	y0 = 0
	vx0 = 0
	vyu = -vy0(0,0.1,[x0,y0,vx0])#sign negative for retrograde orbits
	u0 = [x0, 0, 0,vyu]
	
end

# ╔═╡ 45989493-fcf5-41d8-9427-7ae9246c3852
ds = ContinuousDynamicalSystem(loop,u0,0.1, loop_jac)

# ╔═╡ de85306f-043f-47aa-813e-36aba1695c37
tr = trajectory(ds, 1000.0; dt = 0.01, abstol = 1e-14, reltol = 1e-14)

# ╔═╡ fc8e582b-7281-4c6a-a6e7-325cf02dca00
scatter(tr[:,1],tr[:,2])

# ╔═╡ 94752d46-3c0b-4330-8e52-ebae31e4647b
begin
	plane = (3, 0)
	psos = poincaresos(ds, plane, 20000; u0 = u0,direction = 1,rootkw = (xrtol = 1e-9, atol = 1e-9))
	scatter(psos[:,1], psos[:,2], ms=0.5,title = "Poincaré Plot!")
end

# ╔═╡ 4053df00-aba0-4a32-8493-51b8f1bc7756


# ╔═╡ Cell order:
# ╠═1b497e02-b4df-11eb-2dcd-4f4ac0a22336
# ╠═4cf2fb21-30c3-4a15-8705-4217047f3c7a
# ╠═1a300c74-bdc0-4953-85f6-e43055028ffc
# ╠═3702611c-91e7-4277-8a3a-f7eb0862892c
# ╠═45989493-fcf5-41d8-9427-7ae9246c3852
# ╠═de85306f-043f-47aa-813e-36aba1695c37
# ╠═fc8e582b-7281-4c6a-a6e7-325cf02dca00
# ╠═94752d46-3c0b-4330-8e52-ebae31e4647b
# ╠═4053df00-aba0-4a32-8493-51b8f1bc7756
