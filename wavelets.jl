# a collection of wavelets according to 
#  Continuous Wavelet Transform Library
#  Copyright (C) 2005 Stepan V.Karpenko <carp@mail.ru>

#  Haar. Real part (Time domain)


const TINY = nextfloat(0.0)
function HAARreal(x::Float64, a::Float64, b::Float64)
    a= a==0.0 ? TINY : a
    x = (x-b)/a
    if x >= -0.5 && x < 0
	return 1.0
    end
    if x >= 0.0 && x < 0.5
        return -1.0
    end
    return 0.0
end
	 
# French Hat. Real part (Time domain) 
function FRHATreal(x::Float64, a::Float64, b::Float64)
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    if abs(x) <= 1.0/3.0                    
        return  1.0
    end
    if abs(x) >  1.0/3.0 && abs(x) <= 1.0  
        return -1.0/2.0
    end
    return 0.0;
end

# Splash. Real part (Time domain) 
function SPLASHreal(x::Float64, a::Float64, b::Float64)
    Fb = 0.2 # bandwidth parameter 
    c = sqrt( 2.0/(Fb*Fb*Fb) )
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    return c * x * exp( -abs(x)/Fb )
end

# Mexican Hat. Real part (Time domain) 
function MEXHATreal(x::Float64, a::Float64, b::Float64)
    # c = 2 / ( sqrt(3) * pi^(1/4) ) 
    c = 0.8673250705840776
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    return c * (1.0 - x2) * exp(-x2/2)
    end

# Poisson. Real part (Time domain)  
function POISSONreal(x::Float64, a::Float64, b::Float64)
    # c = 2.0 / sqrt(pi) */
    c = 1.128379167095513
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    dnm = 1.0 + x2
    dnm = dnm==0 ? TINY : dnm
    return c * (1.0 - x2) / ( dnm*dnm )
    end

# Morlet. Real part (Time domain) 
function MORLETreal(x::Float64, a::Float64, b::Float64)
    Fb = 2.0 # bandwidth parameter 
    Fc = 0.8 # wavelet center frequency 
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    return exp(-x2/Fb) * cos(2*pi*Fc*x)
    end

# Complex Morlet. Real part (Time domain) *
function CMORLETreal(x::Float64, a::Float64, b::Float64)
    Fb = 2.0 # bandwidth parameter 
    Fc = 0.8 # wavelet center frequency 
    c = 1.0 / sqrt(pi*Fb)
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    return c * exp(-x2/Fb) * cos(2*pi*Fc*x)
    end

# Complex Shannon wavelet. Real part (Time domain) 
function SHANNONreal(x::Float64, a::Float64, b::Float64)
    # c = 1 / Integral( f(x)^2, x=-Inf..+Inf )^(1/2) 
    c = 2.2360679774997897  # L2 norm 
    Fb = 2.0                # bandwidth parameter 
    Fc = 0.8                # wavelet center frequency 
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x= x == 0.0 ? TINY : x
    return c * sin(pi*Fb*x) * cos(2*pi*Fc*x) / ( sqrt(Fb) * pi * x )
    end

# Complex Shannon wavelet with exponential decay. Real part (Time domain) 
function ESHANNONreal(x::Float64, a::Float64, b::Float64)
    # c = 1 / Integral( f(x)^2, x=-Inf..+Inf )^(1/2) 
    c =  5.9605667650728377
    Fb = 1.0 # bandwidth parameter 
    Fc = 1.0 # wavelet center frequency 
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x= x == 0.0 ? TINY : x
    return c*sin(pi*Fb*x)*cos(2*pi*Fc*x)/(sqrt(Fb)*pi*x*exp(Fb*abs(x)))
    end

# Gaussian 1. Real part (Time domain) 
function GAUSS1real(x::Float64, a::Float64, b::Float64)
    # c = sqrt(2) / pi^(1/4) 
    c = 1.062251932027197
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    return -c * x * exp(-x2/2)
    end

function GAUSS2real(x::Float64, a::Float64, b::Float64)
    # c = 2 / ( sqrt(3) * pi^(1/4) ) *
    c = 0.8673250705840776
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    return c * exp(-x2/2) * (x2 - 1.0)
    end

function GAUSS3real(x::Float64, a::Float64, b::Float64)
    # c = 2 * sqrt(2) / ( sqrt(15) * pi^(1/4) ) 
    c = 0.5485445389623982
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    return -c * x * exp(-x2/2) * (x2 - 3.0)
    end

function GAUSS4real(x::Float64, a::Float64, b::Float64)
    # c = 4 / ( sqrt(105) * pi^(1/4) ) 
    c = 0.2932093894547376
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    x4 = x2 * x2
    return c * exp(-x2/2) * (x4 - 6.0*x2 + 3.0)
    end

function GAUSS5real(x::Float64, a::Float64, b::Float64)
    # c = 4 * sqrt(2) / ( sqrt(945) * pi^(1/4) ) 
    c = 0.1382202317273416
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    x4 = x2 * x2
    return -c * x * exp(-x2/2) * (x4 - 10.0*x2 + 15.0)
    end

function GAUSS6real(x::Float64, a::Float64, b::Float64)
    # c = 8 / ( sqrt(10395) * pi^(1/4) ) *
    c = 0.05893730483821539
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    x4 = x2 * x2
    x6 = x4 * x2
    return c * exp(-x2/2) * (x6 - 15.0*x4 + 45.0*x2 - 15.0)
    end

function GAUSS7real(x::Float64, a::Float64, b::Float64)
    # c = 8 * sqrt(2) / ( sqrt(135135) * pi^(1/4) ) 
    c = 0.02311711288066360;
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    x4 = x2 * x2
    x6 = x4 * x2
    return -c * x * exp(-x2/2) * (x6 - 21.0*x4 + 105.0*x2 - 105.0)
    end

function GAUSS8real(x::Float64, a::Float64, b::Float64)
    # c = 16 / ( sqrt(2027025) * pi^(1/4) ) 
    c = 0.008441176126088454;
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    x4 = x2 * x2
    x6 = x4 * x2
    x8 = x4 * x4
    return c * exp(-x2/2) * (x8 - 28.0*x6 + 210.0*x4 - 420.0*x2 + 105.0)
    end

function GAUSS9real(x::Float64, a::Float64, b::Float64)
    # c = 16 * sqrt(2) / ( sqrt(34459425) * pi^(1/4) ) 
    c = 0.002895299525125788;
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    x4 = x2 * x2
    x6 = x4 * x2
    x8 = x4 * x4
    return -c * x * exp(-x2/2) * (x8 - 36.0*x6 + 378.0*x4 - 1260.0*x2 + 945.0)
    end

function GAUSS10real(x::Float64, a::Float64, b::Float64)
    # c = 32 / ( sqrt(654729075) * pi^(1/4) )
    c = 0.0009393592071302543;
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    x4 = x2 * x2
    x6 = x4 * x2
    x8 = x4 * x4
    x10 = x8 * x2
    return c*exp(-x2/2)*(x10-45.0*x8 +630.0*x6-3150.0*x4+4725.0*x2-945.0)
    end

# Paul 1. Real part (Time domain) 
function PAUL1real(x::Float64, a::Float64, b::Float64)
    # c = 2 * sqrt(2) / sqrt(pi) 
    c = 1.595769121605731
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    x4 = x2 * x2
    return -c * x / (x4 + 2.0*x2 + 1.0)
    end

function PAUL2real(x::Float64, a::Float64, b::Float64)
    # c = 2 * sqrt(2) / sqrt(3*pi) 
    c = 0.9213177319235614;
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    x4 = x2 * x2
    x6 = x4 * x2
    return c * (3.0*x2 - 1.0) / (x6 + 3.0*x4 + 3.0*x2 + 1.0)
    end

function PAUL3real(x::Float64, a::Float64, b::Float64)
    # c = 16 / sqrt(5*pi) 
    c = 4.037012035232256
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    x4 = x2 * x2
    x6 = x4 * x2
    x8 = x4 * x4
    return -c * x * (x2 - 1.0) / (x8 + 4.0*x6 + 6.0*x4 + 4.0*x2 + 1.0)
    end

function PAUL4real(x::Float64, a::Float64, b::Float64)
    # c = 8 * sqrt(2) / sqrt(35*pi) 
    c = 1.078936850151577;
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    x4 = x2 * x2
    x6 = x4 * x2
    x8 = x4 * x4
    x10 = x8 * x2
    return c * (5.0*x4 - 10.0*x2 + 1.0) / 
           (x10 + 5.0*x8 + 10.0*x6 + 10.0*x4 + 5.0*x2 + 1.0)
    end

function PAUL5real(x::Float64, a::Float64, b::Float64)
    # c = 32 / ( 3 * sqrt(7*pi) ) 
    c = 2.274598598644513;
    a= a==0.0 ? TINY : a
    x = (x - b) / a
    x2 = x * x
    x4 = x2 * x2
    x6 = x4 * x2
    x8 = x4 * x4
    x10 = x8 * x2
    x12 = x6 * x6
    return -c * x * (3.0*x4 - 10.0*x2 + 3.0) /
             (x12 + 6.0*x10 + 15.0*x8 + 20.0*x6 + 15.0*x4 + 6.0*x2 + 1.0)
    end

# generic function using function pointer: y=W(PAUL4real,x,1.0,0.0)
function W(f, x, a::Float64, b::Float64)
    if typeof(x) == Float64
        return f(x,a,b)
    end
    if typeof(x) == Array{Float64,1}
        A=[f(x[i],a,b) for i in 1:length(x)]
        return A
    end
    return NaN
 end

# dict for all defined wavelets

type Wavelet
    esl
    esh
    f 
end

global D=Dict()
D["Haar"]=Wavelet(-1.0,1.0,HAARreal)
D["FrHat"]=Wavelet(-1.5,1.5,FRHATreal)
D["Splash"]=Wavelet(-2.0,2.0,SPLASHreal)
D["MexHat"]=Wavelet(-5.0,5.0,MEXHATreal)
D["Poisson"]=Wavelet(-15.0,15.0,POISSONreal)
D["Morlet"]=Wavelet(-4.0,4.0,MORLETreal)
D["CMorlet"]=Wavelet(-4.0,4.0,CMORLETreal)
D["Shannon"]=Wavelet(-10.0,10.0,SHANNONreal)
D["EShannon"]=Wavelet(-4.0,4.0,ESHANNONreal)
D["Wave"]=Wavelet(-5.0,5.0,GAUSS1real)
D["Gauss1"]=Wavelet(-5.0,5.0,GAUSS1real)
D["Gauss2"]=Wavelet(-5.0,5.0,GAUSS2real)
D["Gauss3"]=Wavelet(-5.0,5.0,GAUSS3real)
D["Gauss4"]=Wavelet(-5.0,5.0,GAUSS4real)
D["Gauss5"]=Wavelet(-5.0,5.0,GAUSS5real)
D["Gauss6"]=Wavelet(-5.0,5.0,GAUSS6real)
D["Gauss7"]=Wavelet(-5.0,5.0,GAUSS7real)
D["Gauss8"]=Wavelet(-5.0,5.0,GAUSS8real)
D["Gauss9"]=Wavelet(-5.0,5.0,GAUSS9real)
D["Gauss10"]=Wavelet(-5.0,5.0,GAUSS10real)
D["Paul1"]=Wavelet(-12.0,12.0,PAUL1real)
D["Paul2"]=Wavelet(-6.0,6.0,PAUL2real)
D["Paul3"]=Wavelet(-4.0,4.0,PAUL3real)
D["Paul4"]=Wavelet(-3.0,3.0,PAUL4real)
D["Paul5"]=Wavelet(-3.0,3.0,PAUL5real)

# high level funtions for wavelet transformation
# conv_line convolves a data array x with a W of the same length

function conv_line(wname, x, a)
    # a is the size of window in pixels of the length(x)
    if typeof(x) != Array{Float64,1}
        return NaN
    end
    if haskey(D,wname) == false
        return NaN
    end
    f=D[wname].f
    b=0.0 # a/2.0
    iter=length(x)  
    A=zeros(iter)
    for i = 1:iter          # iterate over all elements of A
        for j = 1:iter      # convolution step 
            A[i]+=f(convert(Float64,j),a,b)*x[j]
        end
        A[i] /= sqrt(a)
        b+=1
    end
    return A
end

# main function to produce an image of the wavlet transformation
# inputs: wname, x
#
# parameters:
#   width of the resulting image in pixels
#   height of the resulting image in pixels
#   check:
#       width<=length(x) && hight<=length(x)
#
# wavelet parameters in pixels:
#   a = 1   start with 1 as default
#  da = floor(length(x)/hight) a= a==0 ? 1 : a
#  da = parameter
#   b = 0   start with no shift
#  db = floor(length(x)/width) b= b==0 ? 1 : b

function cwt(wname::ASCIIString,
             x::Array{Float64,1},
             height::Int64,width::Int64,
             da::Float64=0)
    # perform all checks
    if haskey(D,wname) == false
        return NaN
    end
    l::Int64=length(x)
    if width>l || height>l
        return NaN
    end
    # calculate the paramters for iteration
    if da==0.0
        da= floor(l/hight)
        da = da==0.0 ? 1.0 : da
    end
    db = floor(l/width)
    db = db==0.0 ? 1.0 : db
    a0=1.0
    b0=0.0
    println("da=",da," db=",db)
    f=D[wname].f
    # define the output array and fill it with zero
    A = zeros(height,width)
    # main loop
    a=a0
    for i::Int32=1:height          # loop for a
        b=b0
        for j::Int32=1:width      # loop for b
            for k::Float64=1.0:l  # Perform convolution 
                A[i,j] += f(k,a,b)*x[k]
            end
            A[i,j] /= sqrt(a)
            b += db
         end
        a += da
    end
    return A
end

function cwtf(wname::ASCIIString,
             x::Array{Float64,1},
             height::Int64,width::Int64,
             da::Float64=0)
    # perform all checks
    if haskey(D,wname) == false
        return NaN
    end
    l::Int64=length(x)
    if width>l || height>l
        return NaN
    end
    # calculate the paramters for iteration
    if da==0.0
        da= floor(l/height)
        da = da==0.0 ? 1.0 : da
    end
    db = floor(l/width)
    db = db==0.0 ? 1.0 : db
    a0=1.0
    b0=0.0
    println("da=",da," db=",db)
    f=D[wname].f
    # define the output array and fill it with zero
    A = zeros(height,width)
    # main loop
    a=a0
    l3=3*l
    l2=int(l3/2)
    FX=zeros(l3)                   # define an array for convolution
    for i::Int32=1:height          # loop for a
        b=b0
        for k1::Float64=1.0:l3 # fill F
                FX[k1]=f(k1,a,l3/2.0)
        end
        for j::Int32=1:width      # loop for b
            for k::Int32=1:l  # Perform convolution 
                A[i,j] += FX[(l2-b+k)]*x[k]
            end
            A[i,j] /= sqrt(a)
            b += db
         end
        a += da
    end
    return A
end
    
# simple test function returns an 250x250 Array{Float64,2}
# can be visualized usind PyPlot
# imshow(z)

function test(wname)
    t=linspace(0,10*pi,1000)
    y=sin(t)
    y+=sin(3*t)
    y+=sin(5*t)
    y+=sin(7*t)
    y+=sin(9*t)
    z=cwt(wname,y,250,250,2.0)
    return z,y
end

function testr(wname)
    y=zeros(1000)
    y[100:200]=1
    y[300:400]=1
    y[500:600]=1
    y[700:800]=1
    y[900:1000]=1
    z=cwt(wname,y,250,250,2.0)
    return z,y
end

function testg10(wname)
    y = [GAUSS2real(t,50.0,500.0) for t=1.0:1000]
    z=cwt(wname,y,250,250,2.0)
    return z,y
end

function testsweep(wname)
    sweep=[sin((1.-x/1001.0)^2*pi*20) for x=1:1000]
    z=cwtf(wname,sweep,500,500,2.0)
    return z,sweep
end
