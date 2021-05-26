### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 45df12de-1560-4e86-a8ad-246a6fc29158
using DynamicalSystems, Plots # also exports relevant StaticArrays names

# ╔═╡ 646d9821-3730-401c-b95d-1aa90716c744
# Lorenz system
# Equations of motion:
@inline @inbounds function loop(u, p, t)
	σ = p[1]; ρ = p[2]; β = p[3]
	du1 = σ*(u[2]-u[1])
	du2 = u[1]*(ρ-u[3]) - u[2]
	du3 = u[1]*u[2] - β*u[3]
	return SVector{3}(du1, du2, du3)
end

# ╔═╡ 4901065e-5ddf-4ada-bba6-66c7068e804c
@inline @inbounds function loop_jac(u, p, t)
	σ, ρ, β = p
	J = @SMatrix [-σ  σ  0;
	ρ - u[3]  (-1)  (-u[1]);
	u[2]   u[1]  -β]
	return J
end	


# ╔═╡ 1ff35071-5db8-43b9-9fc7-0a88ca1c8f8f
ds = ContinuousDynamicalSystem(loop, rand(3), [10.0, 28.0, 8/3], loop_jac)

# ╔═╡ eae61860-d594-41f9-836a-e485f8272d98
tr = trajectory(ds, 100.0; dt = 0.001, abstol = 1e-14, reltol = 1e-14)

# ╔═╡ bef4a800-afcd-4b41-ab4e-c530dc64751a
plot(tr[:,1],tr[:,2])

# ╔═╡ Cell order:
# ╠═45df12de-1560-4e86-a8ad-246a6fc29158
# ╠═646d9821-3730-401c-b95d-1aa90716c744
# ╠═4901065e-5ddf-4ada-bba6-66c7068e804c
# ╠═1ff35071-5db8-43b9-9fc7-0a88ca1c8f8f
# ╠═eae61860-d594-41f9-836a-e485f8272d98
# ╠═bef4a800-afcd-4b41-ab4e-c530dc64751a
