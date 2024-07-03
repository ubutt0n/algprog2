struct Polynomial{T}
    coeff::Vector{T}
    function Polynomial{T}(coeff)where T
        k = findfirst(x -> !iszero(x), coeff)
        return new(coeff[k:end])
    end
end
ord(p::Polynomial) = length(p.coeff) - 1
Base. +(p::Polynomial{T}, q::Polynomial{T}) where T = begin
    n, m = length(p.coeff), length(q.coeff)
    if n > m
        r = copy(p.coeff)
        r[n-m+1:end] .+= q.coeff
    else
        r = copy(q.coeff)
        r[m-n+1:end] .+= p.coeff
    end
    return Polynomial{T}(r)
end

Base. +(p::Polynomial{T}, q::T) where T = p + Polynomial{T}([q])
Base. +(p::T, q::Polynomial{T}) where T = q + Polynomial{T}([p])

Base. -(p::Polynomial{T}) where T = Polynomial{T}(-p.coeff)
Base. -(p::Polynomial{T}, q::Polynomial{T}) where T = p + (-q)
Base. -(p::Polynomial{T}, q::T) where T = p + Polynomial{T}([-q])
Base. -(p::T, q::Polynomial{T}) where T = Polynomial{T}([p]) - q

Base. *(p::Polynomial{T}, q::Polynomial{T}) where T = begin
    res = zeros(T, ord(p) + ord(q) + 1)
    p_new = []
    append!(p_new, p.coeff)
    append!(p_new, zeros(Int, ord(q)))

    q_new = []
    append!(q_new, q.coeff)
    append!(q_new, zeros(Int, ord(p)))

    for i in range(0, (ord(p) + ord(q)))
        sum = 0
        for s in range(0, i)
            sum = sum + (p_new[s+1]*q_new[i-s+1])
        end
        res[i+1] = sum
    end
    return Polynomial{T}(res)
end

Base. *(p::Polynomial{T}, q::T) where T = begin
    res = []
    for i in p.coeff
        append!(res, i*q)
    end
    return Polynomial{T}(res)
end
Base. *(p::T, q::Polynomial{T}) where T = begin
    res = []
    for i in q.coeff
        append!(res, i*p)
    end
    return Polynomial{T}(res)
end

#=Base. divrem(p::Polynomial{T}, q::Polynomial{T}) where T = begin
    r = Polynomial{AbstractFloat}(p.coeff)
    q_2 = Polynomial{AbstractFloat}(q.coeff)
    d_coeff = zeros(AbstractFloat, ord(p)-ord(q)+1)
    d_index = firstindex(d_coeff)
    ord_r = ord(r)
    while ord(r) >= ord(q_2)
        d_coeff[d_index] = r.coeff[begin]/q_2.coeff[begin]
        c = r.coeff[begin]/q_2.coeff[begin]
        if r.coeff == (c*q_2).coeff
            return [Polynomial{T}(d_coeff), 0]
        end
        r = r - c*q_2
        @assert ord(r) < ord_r
        d_index += ord_r - ord(r)
        ord_r = ord(r)
    end

    return [Polynomial{T}(d_coeff), r]
end=#

Base. copy(p::Polynomial{T}) where T = Polynomial{T}(p.coeff)

Base. divrem(p::Polynomial{T}, q::Polynomial{T}) where T = begin
    r = copy(p)
    d_coeff = zeros(T, ord(p)-ord(q)+1)
    d_index = firstindex(d_coeff)
    ord_r = ord(r)
    while ord(r) >= ord(q)
        d_coeff[d_index] = r.coeff[begin]/q.coeff[begin]
        c = T(r.coeff[begin]/q.coeff[begin])
        r = r - c*q
        @assert ord(r) < ord_r
        d_index += ord_r - ord(r)
        ord_r = ord(r)
    end
    d = Polynomial{T}(d_coeff)
    return d, r
end

Base. ÷(p::Polynomial{T}, q::Polynomial{T}) where T = divrem(p::Polynomial{T}, q::Polynomial{T})[1]
Base. %(p::Polynomial{T}, q::Polynomial{T}) where T = divrem(p::Polynomial{T}, q::Polynomial{T})[2]
(p::Polynomial)(x) = begin 
    res = 0
    for i in range(0, ord(p))
        res = res*x + p.coeff[i+1]
    end
    return res
end

valdiff(p::Polynomial, x) = begin
    res = 0
    res_diff = 0
    for i in range(0, ord(p))
        res_diff = res_diff*x + res
        res = res*x + p.coeff[i+1]
    end
    return res, res_diff
end
Base. display(p::Polynomial) = begin
    res = ""
    for i in range(0, ord(p))
        if i != ord(p)
            res = res*"$(p.coeff[i+1])x^$(ord(p)-i) + "
        else
            res = res*"$(p.coeff[i+1])"
        end
    end
    print(res)
end
struct Dual{T} <: Number
    a::T
    b::T
end
real(x::Dual) = x.a
imag(x::Dual) = x.b
conj(x::Dual{T}) where T = Dual{T}(x.a, -x.b)

zero(::Type{Dual{T}}) where T = Dual{T}(Base.zero(T), Base.zero(T))
one(::Type{Dual{T}}) where T = Dual{T}(Base.one(T), Base.zero(T))
eps(::Type{Dual{T}}) where T = Dual{T}(Base.zero(T), Base.one(T))

Base. +(x::Dual{T}, y::Dual{T}) where T = Dual{T}(x.a+y.a, x.b+y.b)
Base. +(x::Dual{T}, y::T) where T = Dual{T}(x.a+y, x.b)
Base. +(x::T, y::Dual{T}) where T = y + x
Base. +(x::Dual{T}, y::Number) where T = Dual{T}(x.a+y, x.b)
Base. +(x::Number, y::Dual{T}) where T = y + x

Base. -(x::Dual{T}) where T = Dual{T}(-x.a, -x.b)
Base. -(x::Dual{T}, y::Dual{T}) where T = x + (-y)
Base. -(x::Dual{T}, y::T) where T = Dual{T}(x.a - y, x.b)
Base. -(x::T, y::Dual{T}) where T = -y + x
Base. -(x::Dual{T}, y::Number) where T = Dual{T}(x.a - y, x.b)
Base. -(x::Number, y::Dual{T}) where T = -y + x

Base. *(x::Dual{T}, y::Dual{T}) where T = Dual{T}(x.a*y.a, x.a*y.b + x.b*y.a)
Base. *(x::Dual{T}, y::T) where T = Dual{T}(x.a*y, x.b*y)
Base. *(x::T, y::Dual{T}) where T = Dual{T}(y.a*x, y.b*x)
Base. *(x::Dual{T}, y::Number) where T = Dual{T}(x.a*y, x.b*y)
Base. *(x::Number, y::Dual{T}) where T = Dual{T}(y.a*x, y.b*x)

Base. /(x::Dual{T}, y::T) where T = Dual{AbstractFloat}(x.a/y, x.b/y)
Base. /(x::Dual{T}, y::Number) where T = Dual{AbstractFloat}(x.a/y, x.b/y)
Base. /(x::Dual{T}, y::Dual{T}) where T = (x * conj(y))/(y.a*y.a)
Base. /(x::T, y::Dual{T}) where T = Dual{T}(x, 0) / y
Base. /(x::Number, y::Dual{T}) where T = Dual{T}(x, 0) / y

Base. ^(x::Dual{T}, y::Integer) where T = begin
    res = Dual{T}(x.a, x.b)
    if y == 0
        return Dual{T}(1, 0)
    end
    for _ in range(1, y-1)
        res = res * x
    end
    return res
end 

#(a + bϵ)^(c + dϵ) = aᶜ + ϵ(b(ca^(c-1)) + d(a^c * ln(a))
Base. ^(x::Dual{T}, y::Dual{T}) where T = Dual{AbstractFloat}(x.a^y.a, (x.b*(y.a*(x.a^(y.a-1))) + y.b*((x.a^y.a) * log(x.a))))
Base. ^(x::Dual{T}, y::T) where T = x^Dual{T}(y, 0)
Base. ^(x::Dual{T}, y::Number) where T = x^Dual{T}(y, 0)
Base. ^(x::T, y::Dual{T}) where T = Dual{T}(x, 0)^y
Base. ^(x::Number, y::Dual{T}) where T = Dual{T}(x, 0)^y

Base. sqrt(x::Dual{T}) where T = Dual{AbstractFloat}(sqrt(x.a), (x.b/(2*sqrt(x.a))))

Base. sin(x::Dual{T}) where T = Dual{T}(sin(x.a), x.b*cos(x.a))
Base. cos(x::Dual{T}) where T = Dual{T}(cos(x.a), x.b*(-sin(x.a)))
Base. tan(x::Dual{T}) where T = Dual{T}(tan(x.a), x.b*(1/(cos(x.a)^2)))
Base. cot(x::Dual{T}) where T = Dual{T}(cot(x.a), x.b*(-1/(sin(x.a)^2)))
Base. asin(x::Dual{T}) where T = Dual{T}(asin(x.a), x.b*(1/sqrt(1 - x.a^2)))
Base. acos(x::Dual{T}) where T = Dual{T}(acos(x.a), x.b*(-1/sqrt(1 - x.a^2)))
Base. atan(x::Dual{T}) where T = Dual{T}(atan(x.a), x.b*(1/(1 + x.a^2)))
Base. acot(x::Dual{T}) where T = Dual{T}(acot(x.a), x.b*(-1/(1 + x.a^2)))
Base. exp(x::Dual{T}) where T = Dual{T}(exp(x.a), x.b*exp(x.a))
Base. log(x::Dual{T}) where T = Dual{T}(log(x.a), x.b*(1/x.a))
Base. log2(x::Dual{T}) where T = Dual{T}(log2(x.a), x.b*(1/(x.a*log(2))))
Base. log10(x::Dual{T}) where T = Dual{T}(log10(x.a), x.b*(1/(x.a*log(10))))
Base. log(a::AbstractFloat, x::Dual{T}) where T = Dual{T}(log(a, x.a), x.b*(1/(x.a*log(a))))
Base. sqrt(x::Dual{T}) where T = Dual{T}(sqrt(x.a), x.b*(1/(2*sqrt(x.a))))
valdiff(f::Function, x) = begin
    a = f(Dual{AbstractFloat}(x, 1))
    return [real(a), imag(a)]
end
valdiff(p::Polynomial, x, ::Type{Dual}) = begin
    return valdiff(x -> p(x), x)
end