using CSV
using DataFrames
using DifferentialMobilityAnalyzers
using Underscores
using Distributions
using Memoize
using Dates
using Random
using LinearAlgebra
using Lazy
using RegularizationTools
using Memoize

import Base.|>
|>(args...) = args[end](args[1:end-1]...)

@memoize function initializeDMAs(Dd, q)
    t, p = 295.15, 1e5
    r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
    Λ₁ = DMAconfig(t, p, q.sa1, q.sh1, r₁, r₂, l, 0.0, :-, 10, :cylindrical)
    Λ₂ = DMAconfig(t, p, q.sa2, q.sh2, r₁, r₂, l, 0.0, :-, 10, :cylindrical)
    δ₁ = setupDMA(Λ₁, dtoz(Λ₁, 500e-9), dtoz(Λ₁, 10e-9), 120)
    δ₂ = setupDMA(Λ₂, dtoz(Λ₂, 5.0 * Dd), dtoz(Λ₂, 0.8 * Dd), 60)
    Λ₁, Λ₂, δ₁, δ₂
end

function poisson_noise(Qcpc, N; t = 2.0)
    Q = Qcpc * 1e6
    c = N * Q * t

    map(c) do i
        f = rand(Poisson(i))
        f[1] / (Q * t)
    end
end

@memoize function TDMAFunctions(Dd, q)
    Λ₁, Λ₂, δ₁, δ₂ = initializeDMAs(Dd, q)

    O(k) = mapfoldl(zs -> (δ₂.Ω(Λ₂, δ₁.Z, zs / k, k) .* δ₂.Tl(Λ₂, δ₁.Z, k))', vcat, δ₂.Z)
    T₁(zˢ, k) = δ₁.Ω(Λ₁, δ₁.Z, zˢ / k, k) .* δ₁.Tc(k, δ₁.Dp) .* δ₁.Tl(Λ₁, δ₁.Dp)
    DMA₁(𝕟, zˢ, gf) = @_ map((gf ⋅ (T₁(zˢ, _) * 𝕟)), 1:3)
    DMA₂(𝕟, k) = O(k) * 𝕟

    function TDMA(𝕟, zˢ, gf)
        ℕ = @> DMA₁(𝕟, zˢ, gf)
        r = mapreduce(k -> DMA₂(ℕ[k].N, k), +, 1:3)
        return SizeDistribution([[]], δ₂.De, δ₂.Dp, δ₂.ΔlnD, r ./ δ₂.ΔlnD, r, :TDMA)
    end

    model(𝕟, P, Dd, gf) = sum(@_ map(P[_] * TDMA(𝕟, dtoz(Λ₁, Dd), gf[_]), 1:length(P)))

    return Λ₁, Λ₂, δ₁, δ₂, T₁, DMA₁, DMA₂, model
end

@memoize function getMatrix(𝕟, δ₂, Dd, model)
    mgf = δ₂.Dp ./ (Dd * 1e9)
    λ = @_ map(model(𝕟, 1.0, Dd, _), mgf)
    𝐁 = hcat((@_ map(_.N, λ))...)

    return 𝐁
end


function psd()
    lpm = 1.66666e-5
    Dd = 100e-9
    Ax = [[3500.0, 70.0, 1.6], [900.0, 130.0, 1.4]]
    q = (sh1 = 5.0 * lpm, sa1 = 1.0 * lpm, sh2 = 5.0 * lpm, sa2 = 1.0 * lpm)
    Λ₁, Λ₂, δ₁, δ₂, T₁, DMA₁, DMA₂, model = TDMAFunctions(Dd, q)
    𝕟 = DMALognormalDistribution(Ax, δ₁)
    N = 𝕟.N
    Random.seed!(2000)
    N = poisson_noise(q.sa1, 𝕟.N; t = 2.0)
    S = N ./ 𝕟.ΔlnD
    return SizeDistribution(
        Ax,
        reverse(𝕟.De),
        reverse(𝕟.Dp),
        reverse(𝕟.ΔlnD),
        reverse(S),
        reverse(N),
        :noisy,
    )
end

function list_to_df(file, ts, ℕ; bins = false)
    tint = Dates.value.(ts)
    𝕟 = ℕ[1]
    Dmin = minimum(𝕟.Dp)
    Dmax = maximum(𝕟.Dp)

    theDp = 𝕟.Dp
    upDp = 𝕟.De[2:end]
    lowDp = 𝕟.De[1:end-1]

    function single_line(j)
        df = DataFrame(File = file, timeISO8601 = ts[j], timeInt64 = tint[j])
        theN = @_ map(_ < 0 ? 0 : _, ℕ[j].N)
        map(i -> df[!, Symbol("DMA$i")] = [theN[i]], 1:length(𝕟.Dp))
        return df
    end
    df = map(single_line, 1:length(ts))... |> vcat

    dummy1 = DataFrame(
        File = "Lower bin bounds",
        timeISO8601 = DateTime(2020, 1, 1, 0, 0, 0),
        timeInt64 = Dates.value(DateTime(2020, 1, 1, 0, 0, 0)),
    )
    map(i -> dummy1[!, Symbol("DMA$i")] = [lowDp[i]], 1:length(𝕟.Dp))

    dummy2 = DataFrame(
        File = "Midpoints",
        timeISO8601 = DateTime(2020, 1, 1, 0, 0, 0),
        timeInt64 = Dates.value(DateTime(2020, 1, 1, 0, 0, 0)),
    )
    map(i -> dummy2[!, Symbol("DMA$i")] = [theDp[i]], 1:length(𝕟.Dp))

    dummy3 = DataFrame(
        File = "Upper bin bounds",
        timeISO8601 = DateTime(2020, 1, 1, 0, 0, 0),
        timeInt64 = Dates.value(DateTime(2020, 1, 1, 0, 0, 0)),
    )
    map(i -> dummy3[!, Symbol("DMA$i")] = [upDp[i]], 1:length(𝕟.Dp))

    if bins == true
        return [dummy1; dummy2; dummy3; df]
    else
        return df
    end
end

function list_to_htdma_df(file, ts, Dd, ℕ; bins = false)
    tint = Dates.value.(ts)

    function single_line(j)
        𝕟 = ℕ[j]
        theDp = 𝕟.Dp
        upDp = 𝕟.De[2:end]
        lowDp = 𝕟.De[1:end-1]
        df = DataFrame(
            File = file,
            timeISO8601 = ts[j],
            timeInt64 = tint[j],
            Dd = Dd[j] * 1e9,
        )
        theN = @_ map(_ < 0 ? 0 : _, ℕ[j].N)
        map(i -> df[!, Symbol("DMA$i")] = [theN[i]], 1:length(𝕟.Dp))

        dummy1 = DataFrame(
            File = "Lower bin bounds",
            timeISO8601 = DateTime(2020, 1, 1, 0, 0, 0),
            timeInt64 = Dates.value(DateTime(2020, 1, 1, 0, 0, 0)),
            Dd = 0.0,
        )
        map(i -> dummy1[!, Symbol("DMA$i")] = [lowDp[i]], 1:length(𝕟.Dp))

        dummy2 = DataFrame(
            File = "Midpoints",
            timeISO8601 = DateTime(2020, 1, 1, 0, 0, 0),
            timeInt64 = Dates.value(DateTime(2020, 1, 1, 0, 0, 0)),
            Dd = 0.0,
        )
        map(i -> dummy2[!, Symbol("DMA$i")] = [theDp[i]], 1:length(𝕟.Dp))

        dummy3 = DataFrame(
            File = "Upper bin bounds",
            timeISO8601 = DateTime(2020, 1, 1, 0, 0, 0),
            timeInt64 = Dates.value(DateTime(2020, 1, 1, 0, 0, 0)),
            Dd = 0.0,
        )
        map(i -> dummy3[!, Symbol("DMA$i")] = [upDp[i]], 1:length(𝕟.Dp))

        return [dummy1; dummy2; dummy3; df]
    end
    df = map(single_line, 1:length(ts))... |> vcat

end

function htdma(Dd, i)
    lpm = 1.66666e-5
    Ax = [[3500.0, 70.0, 1.6], [900.0, 130.0, 1.4]]
    q = (sh1 = 5.0 * lpm, sa1 = 0.6 * lpm, sh2 = 5.0 * lpm, sa2 = 1.0 * lpm)
    Λ₁, Λ₂, δ₁, δ₂, T₁, DMA₁, DMA₂, model = TDMAFunctions(Dd, q)
    𝕟 = DMALognormalDistribution(Ax, δ₁)

    𝕣 = model(𝕟, 1.0, Dd, 1.3)
    gf = 𝕣.Dp ./ (Dd * 1e9)

    Normalize(x) = x ./ sum(x)
    Random.seed!(2000)
    if i == 1
        f = @> (0.7 * pdf.(Normal(1.3, 0.1), gf) + pdf.(Normal(1.7, 0.25), gf)) Normalize
    elseif i == 2
        kk = argmin(abs.(1.3 .- gf))
        f = @> zeros(length(gf)) setindex!(1.0, kk)
    else
        kk = argmin(abs.(1.2 .- gf))
        ll = argmin(abs.(1.6 .- gf))
        f = @> zeros(length(gf)) setindex!(0.6, kk) setindex!(0.4, ll)
    end
    𝐁 = getMatrix(𝕟, δ₂, Dd, model)
    N = 𝐁 * f
    N = poisson_noise(q.sa2, N; t = 2.0)

    SizeDistribution(
        Ax,
        reverse(δ₂.De),
        reverse(δ₂.Dp),
        reverse(δ₂.ΔlnD),
        reverse(N ./ δ₂.ΔlnD),
        reverse(N),
        :noisy,
    )
end

function main()
    ts = DateTime(2022, 2, 1, 0, 0, 0):Minute(15):DateTime(2022, 2, 1, 2, 15, 0) |> collect
    ℕ = mapfoldl(_ -> psd(), vcat, 1:10)
    list_to_df("none", ts, ℕ; bins = true) |> CSV.write("smpsgenerated.csv")

    Dd = [50, 100, 150, 50, 100, 150, 50, 100, 150, 50] * 1e-9
    i = [1, 2, 3, 1, 2, 3, 1, 2, 3, 1]
    𝕄 = @>> map(htdma, Dd, i) foldl(vcat)
    list_to_htdma_df("none", ts, Dd, 𝕄; bins = true) |> CSV.write("htdmagenerated.csv")
end

main()
