module PCHIPInterpolationRecipesBaseExt

using PCHIPInterpolation
using RecipesBase

@recipe function _(itp::Interpolator; markershape=:none)
    @series begin
        markershape := :none
        plotdensity = clamp(10 * length(itp.xs), 1000, 100000)
        if itp.extrapolate
            Δxs = maximum(diff(itp.xs)) / 2
            x = range(first(itp.xs) - Δxs, last(itp.xs) + Δxs, length=plotdensity)
        else
            x = range(first(itp.xs), last(itp.xs), length=plotdensity)
        end
        return x, itp.(x)
    end
    if markershape !== :none
        @series begin
            seriestype := :scatter
            primary := false
            return itp.xs, itp.ys
        end
    end
    return nothing
end

end
