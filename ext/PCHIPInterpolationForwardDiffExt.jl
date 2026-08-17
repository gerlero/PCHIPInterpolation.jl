module PCHIPInterpolationForwardDiffExt

using PCHIPInterpolation: Interpolator
import PCHIPInterpolation: _findinterval
using ForwardDiff: Dual

@inline _findinterval(itp::Interpolator{<:AbstractRange}, x::Dual) =
    _findinterval(itp, x.value)
@inline _findinterval(itp::Interpolator, x::Dual) = _findinterval(itp, x.value)

end
