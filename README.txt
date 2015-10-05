The continuous wavelet implementation after:

 Continuous Wavelet Transform Library
 Copyright (C) 2005 Stepan V.Karpenko <carp@mail.ru>

was re implemented to learn the programming in Julia: "http:julialang.org".
Julia is a high-level, high-performance dynamic programming language with
a easy to understand syntax. Julia includes a large mathematical library and
is complement by a large number of packages. The best of Julia is that it can
be used in an IJulia like Python in IPython. The just in time compiling makes 
it easy to get high calculation speed. 

The continuous wavelet transformation is sometimes useful to analyze time 
depended data like the Fourier transformation. But continuous wavelets
can be used even when the frequency of the signal is not constant. I have 
implemented the real part of the wavelet functions in Julia. The main functions
are the:

basic implementation 
cwt(wname::ASCIIString, x::Array{Float64,1}, hight::Int64,width::Int64,
    da::Float64=0) ===> Array{Float64,2}, Array{Float64,1} 

and the fast implementation 
cwtf(wname::ASCIIString, x::Array{Float64,1}, hight::Int64,width::Int64,
da::Float64=0) ===> Array{Float64,2}, Array{Float64,1} 

The wavelet names (wname) can be find in the global dictionary D. The
x is the input data and height, width are the describe the size of the
output array.  The last parameter da is used to control the range of
the shape parameter a. 

Some examples can be found in the testing part: 

test(wname)
testr(wname)
testg10(wname)
testsweep(wname)

Here are som examples how to use the wavelets.jl:

using PyPlot
include("wavelets.jl")

@time z,x=testsweep("Morlet")
t=1:1000
plot(t,x)                         # shows the sweep signal
imshow(z)                         # shows the cwt of the signal
plot(vec(z[30,:]))                # show parts of the cwt

Please compare the calculation time between cwt and cwtf. 

I hope this piece of code is useful for some,

Ralf Wieland. 
