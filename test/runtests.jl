using PCHIPInterpolation
using Test

using ForwardDiff: derivative
using Plots: plot
using QuadGK: quadgk
using OffsetArrays: OffsetVector


function is_monotone(ys)
    direction = first(ys) ≈ last(ys) ? 0 : sign(last(ys) - first(ys))
    return all(
        ya ≈ yb || sign(yb - ya) == direction for
        (ya, yb) in zip(ys[begin:end-1], ys[begin+1:end])
    )
end

function is_piecewise_monotone(itp::Interpolator, N = 10000)
    for i in eachindex(itp.xs)[begin:end-1]
        xsi = range(itp.xs[i], stop = itp.xs[i+1], length = N)
        ysi = itp.(xsi)

        @assert first(ysi) == itp.ys[i]
        @assert last(ysi) == itp.ys[i+1]

        !is_monotone(ysi) && return false
    end
    return true
end

@testset "PCHIPInterpolation" begin

    @testset "helper functions" begin
        @test is_monotone([0.0, 0.0, 0.0])
        @test is_monotone([0.0, 1.0, 10.0])
        @test is_monotone([10.0, 1.0, 0.0])
        @test !is_monotone([0.0, 1.0, 0.0])
        @test !is_monotone([1.0, 0.0, 1.0])
    end

    @testset "interval search" begin
        @testset "valid xs" begin
            xs = 1.0:3.0
            for xs in (xs, collect(xs))
                @assert length(xs) == 3
                @assert PCHIPInterpolation._is_strictly_increasing(xs)
                for search in (
                    x -> (@inferred PCHIPInterpolation._findinterval_base(xs, x)),
                    x -> (@inferred PCHIPInterpolation._findinterval_custom(xs, x)),
                )
                    @test search(0.0) == 0
                    @test search(1.0 - eps(1.0)) == 0
                    @test search(1.0) == 1
                    @test search(1.0 + eps(1.0)) == 1
                    @test search(2.0 - eps(2.0)) == 1
                    @test search(2.0) == 2
                    @test search(2.0 + eps(2.0)) == 2
                    @test search(3.0 - eps(3.0)) == 2
                    @test search(3.0) == 2
                    @test search(3.0 + eps(3.0)) == 3
                    @test search(4.0) == 3
                end
            end
        end

        @testset "xs too short" begin
            for xs in (1.0:1.0, collect(1.0:1.0), 1.0:0.0, collect(1.0:0.0))
                @assert length(xs) < 2
                for search in (
                    x -> (@inferred PCHIPInterpolation._findinterval_base(xs, x)),
                    x -> (@inferred PCHIPInterpolation._findinterval_custom(xs, x)),
                )
                    for x in (-1.0, 0.0, 1.0, 2.0, NaN)
                        @test_throws BoundsError search(x) # Must throw (no valid intervals exist)
                    end
                end
            end
        end
    end

    @testset "Interpolator" begin
        @testset "PCHIP" begin
            xs = [0.0, 1.0]
            ys = [1.0, 1.0]
            @test is_piecewise_monotone(@inferred Interpolator(xs, ys))

            xs = [0.0, 1.0]
            ys = [1.0, 0.0]
            @test is_piecewise_monotone(@inferred Interpolator(xs, ys))

            xs = [0.0, 1.0]
            ys = [0.0, 1.0]
            @test is_piecewise_monotone(@inferred Interpolator(xs, ys))

            xs = [0.0, 1.0, 2.0, 100.0]
            ys = [0.0, 0.0, 1.0, 1.0]
            @test is_piecewise_monotone(@inferred Interpolator(xs, ys))

            xs = [1.0, 2.0, 5.0, 6.0, 10.0, 22.0]
            ys = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
            @test is_piecewise_monotone(@inferred Interpolator(xs, ys))

            xs = [0.0, 10.0, 10000.0]
            ys = [-100.0, 1.2, 1.1]
            @test is_piecewise_monotone(@inferred Interpolator(xs, ys))

            xs = range(0, stop = 2π, length = 5000)
            ys = sin.(xs)
            @test is_piecewise_monotone(@inferred Interpolator(xs, ys))
            @test is_piecewise_monotone(@inferred Interpolator(collect(xs), ys))

            # Example data from Fritsch and Carlson, SIAM J. NUMER. ANAL. 17 (1980) 238-246.
            xs = [0.0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15]
            ys = [10.0, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85]
            @test is_piecewise_monotone(@inferred Interpolator(xs, ys))

            xs = [0.0, 1.2, 2.0, 5.0, 10.0, 11.0]
            ys = [2.0, 2.1, 1.0, 0.0, 0.0, 3.0]
            itp = @inferred Interpolator(xs, ys)
            @test itp.(xs) == ys
            @inferred itp(xs[begin])

            xs = 1:10
            ys = 2:11
            itp = @inferred Interpolator(xs, ys)
            @test itp.(xs) == ys
            @inferred itp(xs[begin])

            xs = 2:101
            ys = collect(3:102)
            itp = @inferred Interpolator(xs, ys)
            @test itp.(xs) == ys
            @inferred itp(xs[begin])
            itp2 = @inferred Interpolator(collect(xs), ys)
            xs2 = [2.73, 3.14, 10, 15, 50, 60, 85]
            @test itp2.(xs2) == itp.(xs2)
            @inferred itp(xs2[begin])
        end

        @testset "custom derivatives" begin
            xs = 1:4
            ys = xs
            ds = -ones(length(ys))
            itp = @inferred Interpolator(xs, ys, ds)
            @test !is_piecewise_monotone(itp)
        end

        @testset "continuity" begin
            xs = [0.0, 1.2, 2.0, 5.0, 10.0, 11.0]
            ys = [2.0, 2.1, 1.0, 0.0, 0.0, 3.0]
            itp = @inferred Interpolator(xs, ys)

            for x in xs[begin+1:end-1]
                y = @inferred itp(x)
                @test (@inferred itp(x - eps(x))) ≈ y atol = 1e-12
                @test (@inferred itp(x + eps(x))) ≈ y atol = 1e-12
            end
        end

        @testset "bad arguments" begin
            @testset "xs too short" begin
                @test_throws ArgumentError Interpolator([1.0], [4.0])
                @test_throws ArgumentError Interpolator([], [])
            end

            @testset "mismatched lengths" begin
                @test_throws DimensionMismatch Interpolator(
                    [1.0, 2.0, 3.0, 4.0],
                    [4.0, 3.0, 2.0],
                )
                @test_throws DimensionMismatch Interpolator(
                    [1.0, 2.0, 3.0],
                    [4.0, 3.0, 2.0, 1.0],
                )
                @test_throws DimensionMismatch Interpolator(
                    [1.0, 2.0, 3.0],
                    [1.0, 2.0, 3.0],
                    [1, 1, 1, 1],
                )
            end

            @testset "non-increasing xs" begin
                @test_throws ArgumentError Interpolator(
                    [1.0, 2.0, 3.0, 3.0],
                    [4.0, 3.0, 2.0, 1.0],
                )
                @test_throws ArgumentError Interpolator(
                    [1.0, 2.0, 3.0, 2.0],
                    [4.0, 3.0, 2.0, 1.0],
                    [1, 1, 1, 1],
                )
            end
        end
    end

    @testset "integrate" begin
        @testset "zero" begin
            itp = @inferred Interpolator([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
            @test (@inferred integrate(itp, 1, 3)) == 0
            @test (@inferred integrate(itp, 3, 1)) == 0
            @test (@inferred integrate(itp, 1, 2)) == 0
            @test (@inferred integrate(itp, 2, 3)) == 0
            @test (@inferred integrate(itp, 1.25, 2)) == 0
            @test (@inferred integrate(itp, 1, 2.25)) == 0
            @test (@inferred integrate(itp, 1.25, 2.25)) == 0
        end

        @testset "linear" begin
            itp = @inferred Interpolator([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0])
            @test integrate(itp, 1, 4) == 3 * 3 / 2 + 3
            @test integrate(itp, 1, 2) + integrate(itp, 2, 4) == integrate(itp, 1, 4)
            @test integrate(itp, 1, 2.75) + integrate(itp, 2.75, 4) == integrate(itp, 1, 4)
            @test integrate(itp, 3, 1) == -integrate(itp, 1, 3)
        end

        @testset "QuadGK" begin
            # Example data from Fritsch and Carlson, SIAM J. NUMER. ANAL. 17 (1980) 238-246.
            xs = [0.0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15]
            ys = [10.0, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85]
            itp = @inferred Interpolator(xs, ys)
            for a = 0:1.5:15, b = 0:1.5:15
                integral, err = quadgk(itp, a, b)
                @test (@inferred integrate(itp, a, b)) ≈ integral atol = err
            end
        end

        @testset "ForwardDiff" begin
            xs = [0.0, 1.2, 2.0, 5.0, 10.0, 11.0]
            ys = [2.0, 2.1, 1.0, 0.0, 0.0, 3.0]
            itp = @inferred Interpolator(xs, ys)
            for x in [0, 1, 3.14, 5, 9, 11]
                @test derivative(b -> integrate(itp, 0, b), x) ≈ itp(x)
                @test derivative(b -> integrate(itp, 1.1, b), x) ≈ itp(x)
                @test derivative(b -> integrate(itp, 11, b), x) ≈ itp(x)
            end
        end
    end

    @testset "ForwardDiff" begin
        @testset "zero" begin
            itp = @inferred Interpolator([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
            @test (@inferred derivative(itp, 1)) == 0
            @test (@inferred derivative(itp, 2)) == 0
            @test (@inferred derivative(itp, 3)) == 0
            @test (@inferred derivative(itp, 1.25)) == 0
            @test (@inferred derivative(itp, 2.25)) == 0
        end

        @testset "linear" begin
            itp = @inferred Interpolator([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0])
            @test (@inferred derivative(itp, 1)) == -1
            @test (@inferred derivative(itp, 2)) == -1
            @test (@inferred derivative(itp, 3)) == -1
            @test (@inferred derivative(itp, 4)) == -1
            @test (@inferred derivative(itp, 2.75)) == -1
        end

        @testset "derivative continuity" begin
            xs = [0.0, 1.2, 2.0, 5.0, 10.0, 11.0]
            ys = [2.0, 2.1, 1.0, 0.0, 0.0, 3.0]
            itp = Interpolator(xs, ys)

            for x in xs[begin+1:end-1]
                d = @inferred derivative(itp, x)
                @test (@inferred derivative(itp, x - eps(x))) ≈ d atol = 1e-12
                @test (@inferred derivative(itp, x + eps(x))) ≈ d atol = 1e-12
            end
        end

        @testset "custom derivatives" begin
            xs = 1:4
            ys = xs
            ds = -ones(length(ys))
            itp = @inferred Interpolator(xs, ys, ds)
            @test all(derivative(itp, x) == d for (x, d) in zip(xs, ds))
        end
    end

    @testset "out of domain" begin
        itp = @inferred Interpolator([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0])
        @test itp(1) == 4
        @test itp(4) == 1
        @test isnan(itp(1 - 1e-6))
        @test isnan(itp(4 + 1e-6))
        @test isnan(integrate(itp, 1 - 1e-6, 4 + 1e-6))
        @test isnan(integrate(itp, 1 - 1e-6, 3))
        @test isnan(integrate(itp, 2, 4 + 1e-6))
    end

    @testset "NaN propagation" begin
        itp = Interpolator(1.0:4.0, [4.0, 3.0, 2.0, 1.0])
        @test isnan(@inferred itp(NaN))
        @test isnan(@inferred integrate(itp, 1, NaN))
        @test isnan(@inferred integrate(itp, NaN, 4))
        @test isnan(@inferred integrate(itp, NaN, NaN))
        @test isnan(@inferred derivative(itp, NaN))

        itp = Interpolator(collect(1.0:4.0), [4.0, 3.0, 2.0, 1.0])
        @test isnan(@inferred itp(NaN))
        @test isnan(@inferred integrate(itp, 1, NaN))
        @test isnan(@inferred integrate(itp, NaN, 4))
        @test isnan(@inferred integrate(itp, NaN, NaN))
        @test isnan(@inferred derivative(itp, NaN))
    end

    @testset "plot recipe" begin
        xs = [0.0, 1.2, 2.0, 5.0, 10.0, 11.0]
        ys = [2.0, 2.1, 1.0, 0.0, 0.0, 3.0]
        itp = @inferred Interpolator(xs, ys)
        plot(itp)
        plot(itp, markershape = :auto)
    end

    @testset "OffsetArrays" begin
        xs = [0.0, 1.2, 2.0, 5.0, 10.0, 11.0]
        ys = [2.0, 2.1, 1.0, 0.0, 0.0, 3.0]

        oxs = OffsetVector(xs, -1)
        oys = OffsetVector(ys, -1)
        @assert firstindex(oxs) == firstindex(oys) == 0

        @testset "valid" begin
            oitp = Interpolator(oxs, oys)
            itp = Interpolator(xs, ys)

            @test oitp(0) == itp(0)
            @test oitp(1) == itp(1)
            @test oitp(3.0) == itp(3.0)
            @test oitp(5.0) == itp(5.0)
            @test oitp(11.0) == itp(11.0)
            @test integrate(oitp, 1.5, 6.5) == integrate(itp, 1.5, 6.5)
            @test integrate(oitp, 5, 11) == integrate(itp, 5, 11)
            @test integrate(oitp, 0, 11) == integrate(itp, 0, 11)
        end

        @testset "out of domain" begin
            oitp = Interpolator(oxs, oys)

            @test isnan(oitp(0.0 - eps(0.0)))
            @test isnan(oitp(11.0 + eps(11.0)))
            @test isnan(integrate(oitp, 0.0 - eps(0.0), 1.0))
            @test isnan(integrate(oitp, 0.0, 11.0 + eps(11.0)))
        end

        @testset "incompatible arguments" begin
            @test_throws DimensionMismatch Interpolator(oxs, ys)
            @test_throws DimensionMismatch Interpolator(xs, oys)
            @test_throws DimensionMismatch Interpolator(0:lastindex(oxs), oys)
            @test_throws DimensionMismatch Interpolator(oys, zeros(lastindex(oys)))
            @test_throws DimensionMismatch Interpolator(oxs, oys, zeros(length(oxs)))
            @test_throws DimensionMismatch Interpolator(oxs, oys, zeros(lastindex(oxs)))
        end
    end

    @testset "regression tests" begin
        @testset "gh-8" begin
            # Regression test for https://github.com/gerlero/PCHIPInterpolation.jl/issues/8
            x = collect(range(0, stop = 1, length = 24))
            itp = Interpolator(x, x)
            @test itp(1) == 1
            x = collect(range(0, stop = 1, length = 25))
            itp = Interpolator(x, x)
            @test itp(1) == 1
            x = collect(range(0, stop = 1, length = 100))
            itp = Interpolator(x, x)
            @test itp(1) == 1
        end

        @testset "gh-31" begin
            # Test correctness of _pchip_ds_scipy implementation following https://github.com/gerlero/PCHIPInterpolation.jl/issues/31
            # Tests that the interpolated values are the same as SciPy
            xs = [0.0, 1.2, 2.0, 5.0, 10.0, 11.0]
            ys = [2.0, 2.1, 1.0, 0.0, 0.0, 3.0]

            xps = [0.6, 1.6, 3.5, 7.5, 10.5]
            yps = [2.08750000, 1.61081474, 0.27194471, 0.00000000, 1.06250000]

            itp = @inferred Interpolator(xs, ys)
            @test itp.(xps) ≈ yps atol = 1e-8
        end
    end
end
