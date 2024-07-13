module PCHIPInterpolation

export Interpolator, integrate


_is_strictly_increasing(xs::AbstractVector) =
    all(a < b for (a, b) in zip(xs[begin:end-1], xs[begin+1:end]))
_is_strictly_increasing(xs::AbstractRange) = step(xs) > 0

function _pchip_edge_derivative(h1, h2, Δ1, Δ2)
    d = ((2h1 + h2) * Δ1 - h2 * Δ2) / (h1 + h2)
    if sign(d) != sign(Δ1)
        d = zero(d)
    elseif sign(Δ1) != sign(Δ2) && abs(d) > abs(3Δ1)
        d = 3Δ1
    end
    return d
end

function _pchip_ds_scipy(xs::AbstractVector, ys::AbstractVector)
    h(i) = xs[i+1] - xs[i]
    Δ(i) = (ys[i+1] - ys[i]) / h(i)

    length(ys) != length(xs) && throw(DimensionMismatch)
    ds = similar(ys ./ xs)

    is = eachindex(xs, ys, ds)

    if length(ds) == 2
        ds[:] .= Δ(first(is))
    else
        Δl = Δ(first(is))
        hl = h(first(is))
        for i in is[begin+1:end-1]
            Δr = Δ(i)
            hr = h(i)
            if sign(Δl) != sign(Δr) || Δl ≈ zero(Δl) || Δr ≈ zero(Δr)
                ds[i] = zero(eltype(ds))
            else
                wl = 2hr + hl
                wr = hr + 2hl
                ds[i] = (wl + wr) / (wl / Δl + wr / Δr)
            end
            Δl = Δr
            hl = hr
        end
        ds[begin] = _pchip_edge_derivative(
            h(first(is)),
            h(first(is) + 1),
            Δ(first(is)),
            Δ(first(is) + 1),
        )
        ds[end] = _pchip_edge_derivative(
            h(last(is) - 1),
            h(last(is) - 2),
            Δ(last(is) - 1),
            Δ(last(is) - 2),
        )
    end
    return ds
end

struct Interpolator{Xs,Ys,Ds}
    xs::Xs
    ys::Ys
    ds::Ds
    extrapolate::Bool

    function Interpolator(xs::AbstractVector, ys::AbstractVector; extrapolate::Bool = false)
        length(eachindex(xs, ys)) ≥ 2 ||
            throw(ArgumentError("inputs must have at least 2 elements"))
        _is_strictly_increasing(xs) ||
            throw(ArgumentError("xs must be strictly increasing"))

        ds = _pchip_ds_scipy(xs, ys)

        new{typeof(xs),typeof(ys),typeof(ds)}(xs, ys, ds, extrapolate)
    end

    function Interpolator(
        xs::AbstractVector,
        ys::AbstractVector,
        ds::AbstractVector;
        extrapolate::Bool = false,
    )
        length(eachindex(xs, ys, ds)) ≥ 2 ||
            throw(ArgumentError("inputs must have at least 2 elements"))
        _is_strictly_increasing(xs) ||
            throw(ArgumentError("xs must be strictly increasing"))

        new{typeof(xs),typeof(ys),typeof(ds)}(xs, ys, ds, extrapolate)
    end
end


@inline function _findinterval_base(xs::AbstractVector, x) # Generic binary search from Julia Base
    @boundscheck length(xs) ≥ 2 || throw(BoundsError(xs))

    VERSION < v"1.6.2" && isnan(x) && return lastindex(xs) - 1

    i = searchsortedlast(xs, x)

    if i == lastindex(xs) && x == @inbounds xs[i]
        i -= 1 # Treat right endpoint as part of rightmost interval
    end

    return i
end

@inline function _findinterval_custom(xs::AbstractVector, x) # SciPy-like binary search from SimplePCHIP
    @boundscheck length(xs) ≥ 2 || throw(BoundsError(xs))

    imin = firstindex(xs)

    if x < @inbounds xs[imin]
        return imin - 1
    end

    imax = lastindex(xs)
    xmax = @inbounds xs[imax]

    if x > xmax
        return imax
    elseif x == xmax
        return imax - 1 # Treat right endpoint as part of rightmost interval
    end

    i = imin + (imax - imin + 1) ÷ 2
    while imin < imax
        if x < @inbounds xs[i]
            imax = i - 1
        elseif x >= @inbounds xs[i+1]
            imin = i + 1
        else
            break
        end
        i = imin + (imax - imin + 1) ÷ 2
    end

    return i
end

@inline _findinterval(itp::Interpolator{<:AbstractRange}, x) = _findinterval_base(itp.xs, x) # Base binary search has an overload for ranges
@inline _findinterval(itp::Interpolator, x) = _findinterval_custom(itp.xs, x) # Otherwise, the custom search is preferred


@inline _x(::Interpolator, x) = x
@inline _x(itp::Interpolator, x, _) = _x(itp, x)
@inline _x(itp::Interpolator, ::Val{:begin}, i) = @inbounds itp.xs[i]
@inline _x(itp::Interpolator, ::Val{:end}, i) = @inbounds itp.xs[i+1]

@inline _evaluate(itp::Interpolator, ::Val{:begin}, i) = @inbounds itp.ys[i]
@inline _evaluate(itp::Interpolator, ::Val{:end}, i) = @inbounds itp.ys[i+1]
@inline _derivative(itp::Interpolator, ::Val{:begin}, i) = @inbounds itp.ds[i]
@inline _derivative(itp::Interpolator, ::Val{:end}, i) = @inbounds itp.ds[i+1]


@inline _ϕ(t) = 3t^2 - 2t^3
@inline _ψ(t) = t^3 - t^2

function _evaluate(itp::Interpolator, x, i)
    if itp.extrapolate
        if i < firstindex(itp.xs)
            i = firstindex(itp.xs)
        elseif i >= lastindex(itp.xs)
            i = lastindex(itp.xs) - 1
        end
    elseif i < firstindex(itp.xs) || i >= lastindex(itp.xs)
        return float(eltype(itp.xs))(NaN)
    end

    x1 = _x(itp, Val(:begin), i)
    x2 = _x(itp, Val(:end), i)
    h = x2 - x1

    y1 = _evaluate(itp, Val(:begin), i)
    y2 = _evaluate(itp, Val(:end), i)

    d1 = _derivative(itp, Val(:begin), i)
    d2 = _derivative(itp, Val(:end), i)

    return (
        y1 * _ϕ((x2 - x) / h) + y2 * _ϕ((x - x1) / h) - d1 * h * _ψ((x2 - x) / h) +
        d2 * h * _ψ((x - x1) / h)
    )
end

@inline _evaluate(itp::Interpolator, x) = _evaluate(itp, x, _findinterval(itp, x))

@inline (itp::Interpolator)(x::Number) = _evaluate(itp, x)


@inline function _integrate(itp::Interpolator, a, b, i)
    if !itp.extrapolate && (i < firstindex(itp.xs) || i >= lastindex(itp.xs))
        return float(eltype(itp.xs))(NaN)
    end

    a_ = _x(itp, a, i)
    b_ = _x(itp, b, i)
    return (b_ - a_) / 6 * (
        _evaluate(itp, a, i) + 4 * _evaluate(itp, (a_ + b_) / 2, i) + _evaluate(itp, b, i)
    ) # Simpson's rule
end

@inline function _integrate(itp::Interpolator, a, b, i, j)
    if i == j
        return _integrate(itp, a, b, i)
    end

    integral = _integrate(itp, a, Val(:end), i)
    for k = i+1:j-1
        integral += _integrate(itp, Val(:begin), Val(:end), k)
    end
    integral += _integrate(itp, Val(:begin), b, j)

    return integral
end

@inline function _integrate(itp::Interpolator, a, b)
    if b < a
        return -_integrate(itp, b, a)
    end

    return _integrate(itp, a, b, _findinterval(itp, a), _findinterval(itp, b))
end

@inline integrate(itp::Interpolator, a::Number, b::Number) = _integrate(itp, a, b)


if !isdefined(Base, :get_extension) # For compatibility with Julia < 1.9
    include("../ext/PCHIPInterpolationRecipesBaseExt.jl")
end

end
