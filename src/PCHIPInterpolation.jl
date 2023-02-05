module PCHIPInterpolation
#
# Simple PCHIP implementation following Fritsch and Carlson, SIAM J. NUMER. ANAL. 17 (1980) 238-246.
# Derivative calculated in a fashion similar to SciPy's PchipInterpolate
#

export Interpolator, integrate

using ArgCheck: @argcheck
using RecipesBase: @recipe, @series, RecipesBase


struct Interpolator{Xs,Ys,Ds}
    xs::Xs
    ys::Ys
    ds::Ds

    function Interpolator(xs::AbstractVector, ys::AbstractVector)
        @argcheck length(xs) ≥ 2
        @argcheck length(xs) == length(ys) DimensionMismatch
        foldl(xs) do a,b
            @argcheck a < b "xs must be strictly increasing"
            return b
        end

        ds = _initial_ds_scipy(xs, ys)
        new{typeof(xs),typeof(ys),typeof(ds)}(deepcopy(xs), deepcopy(ys), ds)
    end

    function Interpolator(xs::AbstractVector, ys::AbstractVector, _ds::AbstractVector)
        @argcheck length(xs) ≥ 2
        @argcheck length(xs) == length(ys) == length(_ds) DimensionMismatch
        foldl(xs) do a,b
            @argcheck a < b "xs must be strictly increasing"
            return b
        end

        new{typeof(xs),typeof(ys),typeof(_ds)}(deepcopy(xs), deepcopy(ys), deepcopy(_ds))
    end
end

_ϕ(t) = 3t^2 - 2t^3
_ψ(t) = t^3 - t^2

@inline function _value_with_index(pchip::Interpolator, x::Number, i)
    x1 = @inbounds pchip.xs[i]
    x2 = @inbounds pchip.xs[i+1]
    @assert x1 ≤ x ≤ x2
    h = x2 - x1

    y1 = pchip.ys[i]
    y2 = pchip.ys[i+1]

    d1 = pchip.ds[i]
    d2 = pchip.ds[i+1]

    return (y1 * _ϕ((x2-x)/h)
           + y2 * _ϕ((x-x1)/h)
           - d1*h * _ψ((x2-x)/h)
           + d2*h * _ψ((x-x1)/h))
end

@inline (pchip::Interpolator)(x::Number) = _value_with_index(pchip, x, _findindex(pchip, x))

function _integrate_segment(pchip::Interpolator, i, x1=nothing, x2=nothing)
    if isnothing(x1)
        x1 = @inbounds pchip.xs[i]
        y1 = pchip.ys[i]
    else
        y1 = _value_with_index(pchip, x1, i)
    end

    if isnothing(x2)
        x2 = @inbounds pchip.xs[i+1]
        y2 = pchip.ys[i+1]
    else
        y2 = _value_with_index(pchip, x2, i)
    end

    return (x2 - x1)/6 * (y1 + 4*_value_with_index(pchip, (x1 + x2)/2, i) + y2)  # Simpson's rule
end

function integrate(pchip::Interpolator, a::Number, b::Number)
    if b < a
        return -integrate(pchip, b, a)
    end

    i = _findindex(pchip, a)
    j = _findindex(pchip, b)

    if i == j
        return _integrate_segment(pchip, i, a, b)
    end

    integral = _integrate_segment(pchip, i, a, nothing) + _integrate_segment(pchip, j, nothing, b)

    for k in range(i+1, stop=j-1)
        integral += _integrate_segment(pchip, k)
    end

    return integral
end


@inline _findindex(pchip::Interpolator{<:AbstractRange}, x) = _findindex_base(pchip.xs, x) # Base binary search has an overload for ranges
@inline _findindex(pchip::Interpolator, x) = _findindex_custom(pchip.xs, x) # Otherwise, the custom search is preferred


@inline function _findindex_base(xs, x) # Generic binary search from Julia Base
    i = searchsortedlast(xs, x)

    if i < firstindex(xs)
        throw(DomainError(x, "Below interpolation range"))
    end

    if i == lastindex(xs)
        if x != @inbounds xs[i]
            throw(DomainError(x, "Above interpolation range"))
        end
        i -= 1 # Treat right endpoint as part of rightmost interval
    end

    return i
end

@inline function _findindex_custom(xs, x) # SciPy-like binary search from SimplePCHIP
    imin = firstindex(xs)

    if x < @inbounds xs[imin]
        throw(DomainError(x, "Below interpolation range"))
    end

    imax = lastindex(xs)
    xmax = @inbounds xs[imax]

    if x > xmax
        throw(DomainError(x, "Above interpolation range"))
    elseif x == xmax
        return imax - 1 # Treat right endpoint as part of rightmost interval
    end

    i = imin + (imax - imin + 1)÷2
    while imin < imax
        if x < @inbounds xs[i]
            imax = i - 1
        elseif x >= @inbounds xs[i+1]
            imin = i + 1
        else
            break
        end
        i = imin + (imax - imin + 1)÷2
    end

    return i
end


"Similar to how SciPy's PCHIP does it"
function _initial_ds_scipy(xs::AbstractVector, ys::AbstractVector)
    h(i) = xs[i+1]-xs[i]
    Δ(i) = (ys[i+1]-ys[i]) / h(i)

    N = length(xs)
    ds = similar(ys./xs)
    if N == 2
        ds[:] .= Δ(1)
    else
        Δl = Δ(1)
        hl = h(1)
        for i ∈ 2:N-1
            Δr = Δ(i)
            hr = h(i)
            if sign(Δl) != sign(Δr) || Δl ≈ 0.0 || Δr ≈ 0.0
                ds[i] = 0.0
            else
                wl = 2hl + hr
                wr = hl + 2hr
                ds[i] = (wl + wr) / (wl/Δl + wr/Δr)
            end
            Δl = Δr
            hl = hr
        end
        ds[1] = _edge_derivative(h(1), h(2), Δ(1), Δ(2))
        ds[N] = _edge_derivative(h(N-1), h(N-2), Δ(N-1), Δ(N-2))
    end
    ds
end

function _edge_derivative(h1, h2, Δ1, Δ2)
    d = ((2h1 + h2)*Δ1 - h2*Δ2) / (h1 + h2)
    if sign(d) != sign(Δ1)
        d = 0.0
    elseif sign(Δ1) != sign(Δ2)  &&  abs(d) > abs(3Δ1)
        d = 3Δ1
    end
    d
end

@recipe function f(pchip::Interpolator; markershape=:none) \
    @series begin
        markershape := :none
        plotdensity = clamp(10length(pchip.xs), 1000, 100000)
        x = range(pchip.xs[1], stop=pchip.xs[end], length=plotdensity)
        return x, pchip.(x)
    end
    if markershape !== :none
        @series begin
            seriestype := :scatter
            primary := false
            return pchip.xs, pchip.ys
        end
    end
    return nothing
end

end  # module
