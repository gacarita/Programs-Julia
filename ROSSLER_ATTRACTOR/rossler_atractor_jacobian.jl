### A Pluto.jl notebook ###
# v0.14.5

using Symbolics,SymbolicUtils, Latexify

@variables x y z a b c

# Rossler System

#a=0.2
#b=0.2
#c=5.7

dxdt = -y-z
dydt = x+a*y
dzdt = b+z*(x-c)



Symbolics.jacobian([dxdt,dydt,dzdt],[x,y,z])
