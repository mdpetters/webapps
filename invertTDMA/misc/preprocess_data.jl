using NetCDF
using DifferentialMobilityAnalyzers
using DataFrames
using CSV
using Underscores
using Dates

import Base.|>
|>(args...) = args[end](args[1:end-1]...)

function loadHTDMA(file)
    basetime = ncread(file, "base_time")
    bt = DateTime(1970, 1, 1, 0, 0, 0) + Dates.Second(convert(Int, basetime[1]))
    offset = ncread(file, "time_offset")
    timestamp = bt .+ Dates.Second.(offset)
    dNdlog10D = ncread(file, "aerosol_concentration")
    Dp = ncread(file, "bin_center")
    ΔD = ncread(file, "bin_width")
    Dd = ncread(file, "dry_diameter_setting")
    RH = ncread(file, "humid_rh")
    Nt = ncread(file, "total_concentration")
    Q = ncread(file, "sample_flow")
    gf = hcat((@_ map(Dp[:, _] ./ Dd[_], 1:length(Dd)))...)

    De1 = Dp .- ΔD / 2
    De2 = Dp .+ ΔD / 2
    Δlog10D = log10.(De2 ./ De1)
    ΔlnD = log.(De2 ./ De1)
    N = dNdlog10D .* Δlog10D

    ℕ = map(1:length(timestamp)) do i
        De = [De1[:, i][1:end]; De2[end, i]]
        SizeDistribution(
            [[]],
            (De),
            (Dp[:, i]),
            (ΔlnD[:, i]),
            (N[:, i] ./ ΔlnD[:, i]),
            (N[:, i]),
            :sgpaoshtdama,
        )
    end

    return timestamp, ℕ, Dd
end

function list_to_htdma_df(file, ts, Dd, ℕ; bins = false)
    tint = Dates.value.(ts)

    function single_line(j)
        𝕟 = ℕ[j]

        theDp = 𝕟.Dp
        upDp = 𝕟.De[2:end]
        lowDp = 𝕟.De[1:end-1]

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

        df = DataFrame(File = file, timeISO8601 = ts[j], timeInt64 = tint[j], Dd = Dd[j])
        theN = @_ map(_ < 0 ? 0 : _, ℕ[j].N)
        map(i -> df[!, Symbol("DMA$i")] = [theN[i]], 1:length(𝕟.Dp))

        return [dummy1; dummy2; dummy3; df]
    end
    df = map(single_line, 1:length(ts))... |> vcat
end

function loadSMPS(file)
    basetime = ncread(file, "base_time")
    bt = DateTime(1970, 1, 1, 0, 0, 0) + Dates.Second(convert(Int, basetime[1]))
    offset = ncread(file, "time_offset")
    timestamp = bt .+ Dates.Second.(offset)
    status_flag = ncread(file, "status_flag")
    mpsd = ncread(file, "number_size_distribution")
    Dp = ncread(file, "diameter_midpoint")
    Nt = ncread(file, "total_concentration")
    Qsh = ncread(file, "sheath_flow")

    if Qsh[1] ≠ 5.0
        throw("Error: wrong sheath flow")
    end

    De = sqrt.(Dp[1:end-1] .* Dp[2:end])
    Dp = Dp[2:end-1]
    Δlog10D = log10.(De[2:end] ./ De[1:end-1])
    ΔlnD = log.(De[2:end] ./ De[1:end-1])

    tpsd = mpsd[2:end-1, :]
    ℕ = map(1:length(timestamp)) do i
        SizeDistribution(
            [[]],
            (De[65:end]),
            (Dp[65:end]),
            (ΔlnD[65:end]),
            (tpsd[65:end, i] .* Δlog10D[65:end] ./ ΔlnD[65:end]),
            (tpsd[65:end, i] .* Δlog10D[65:end]),
            :sgpaossmps,
        )
    end
    return timestamp, ℕ
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

file = "sgpaoshtdma/sgpaoshtdmaE13.a1.20200214.000231.custom.cdf"
ts, ℕ, Dd = loadHTDMA(file)
df = list_to_htdma_df(file, ts, Dd, ℕ)
df |> CSV.write("htdmasgp.csv")

file1 = "sgpaossmps/sgpaossmpsE13.a1.20200214.000000.nc"
ts1, ℕ1 = loadSMPS(file1)

res = map(1:length(ts)) do j
    tts = ts[j]
    ii = argmin(abs.(ts1 .- tts))
    ts1[ii], ℕ1[ii]
end

ts1 = map(x -> x[1], res)
ℕ1 = map(x -> x[2], res)
df = list_to_df(file1, ts1, ℕ1; bins = true)
df |> CSV.write("smpssgp.csv")
df[4, 5:end] |> Vector |> sum
