using Symbolics,SymbolicUtils, Latexify

@variables x̂ ŷ x y μ
a = x̂
b = ŷ
c = 2ŷ  + x - (1 - μ)*(μ + x)/(((x + μ)^2 + y^2)^(3/2)) - μ*(-1 + μ + x)/(((-1 + x + μ)^2 + y^2)^(3/2))
d = -2x̂  + y - (1 - μ)*y/(((x + μ)^2 + y^2)^(3/2)) - μ*y/(((-1 + x + μ)^2 + y^2)^(3/2))
#r1 = ((x + μ)^2 + y^2)^(1/2)
#r2 = ((-1 + x + μ)^2 + y^2)^(1/2)

#Symbolics.jacobian([x̂, ŷ, 2ŷ  + x - (1 - μ)*(μ + x)/(((x + μ)^2 + y^2)^(3/2)) - μ*(-1 + μ + x)/(((-1 + x + μ)^2 + y^2)^(3/2)),
#                        -2x̂  + y - (1 - μ)*y/(((x + μ)^2 + y^2)^(3/2))
#                        - μ*y/(((-1 + x + μ)^2 + y^2)^(3/2))], [x, y, x̂, ŷ])

#A = Symbolics.jacobian([x̂,ŷ,a,b],[x,y,x̂,ŷ])

A = Symbolics.jacobian([d],[y])
C = simplify.(A)
#latexify(A)
#println(C)
