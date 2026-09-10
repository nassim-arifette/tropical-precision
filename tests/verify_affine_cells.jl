# Separate rational verification of cell witnesses and exact pointwise supports.
# This does not establish the universal cell bound or the asymptotic theorem.
if !isdefined(@__MODULE__, :check_combination)
    include("verify_witnesses.jl")
end

function cells_support(x,w,k,sign)
    radius = sign == 1 ? 1-x[k] : x[k]
    bounds = [w[j] >= 0 ? 1-x[j] : x[j] for j in eachindex(x)]
    nodes = Set([0//BigInt(1),radius])
    for j in eachindex(x)
        j != k && push!(nodes,min(radius,bounds[j]))
    end
    maximum(sign*w[k]*s + sum(abs(w[j])*min(s,bounds[j])
                             for j in eachindex(x) if j != k) for s in nodes)
end

function cells_deficit(x,w)
    values = [min(0,sum(w))-dotq(w,x)]
    for j in eachindex(x)
        x[j] > minimum(x) && push!(values,cells_support(x,w,j,1))
        x[j] < maximum(x) && push!(values,cells_support(x,w,j,-1))
    end
    minimum(values)
end

function cell_verify(case)
    n = case["dimension"]
    n >= 3 || error("Cell dimension")
    interval = vector(case["interval"])
    length(interval) == 2 || error("Cell endpoints")
    a,b = interval
    h = b-a
    0 <= a < b <= 1 && h <= min(a,1-b) || error("Interior-cell hypothesis")
    G,T = vector.(case["generators"]),vector.(case["context"])
    x,w = vector(case["target"]),vector(case["affine_coefficients"])
    length(G) == case["generator_count"] || error("Generator count")
    all(g -> length(g) == 2n,G) || error("Generator dimension")
    length(x) == length(w) == n || error("Point/direction dimension")
    all(v -> a <= v <= b,x) || error("Target outside cell")
    all(t -> length(t) == n && all(v -> 0 <= v <= 1,t),T) || error("Context outside invariant")
    sum(abs.(w)) == 1 || error("Affine normalization")
    check_slice(G,n)
    check_combination(vcat(G,lift.(T)),lift(x),vector(case["coefficients"]))
    concrete = minimum(vcat([min(0,sum(w))],[dotq(w,t) for t in T]))
    value = dotq(w,x)
    gap = concrete-value
    concrete == rational(case["minimum_concrete_union"]) || error("Concrete minimum")
    value == rational(case["target_value"]) || error("Target value")
    gap == rational(case["certified_gap"]) || error("Gap")
    formula = n == 3 ? h*max(a,1-b)/(2-h) : h/2
    gap == formula || error("Not an attaining cell witness")
    cells_deficit(x,w) == gap || error("Pointwise sector deficit")
    gap
end

function cell_negative_checks(original)
    mutations = [
        c -> (c["coefficients"][1] = "1/2"),
        c -> (c["context"][1][1] = "2"),
        c -> (c["affine_coefficients"][1] = "2"),
        c -> (c["certified_gap"] = "99"),
        c -> (c["target"][1] = "7/13"),
        c -> (c["interval"][1] = "0"),
        c -> popfirst!(c["generators"]),
        c -> (c["generator_count"] += 1),
    ]
    for mutation in mutations
        bad = deepcopy(original)
        mutation(bad)
        rejected = false
        try
            cell_verify(bad)
        catch
            rejected = true
        end
        rejected || error("Accepted a malformed cell witness")
    end
    length(mutations)
end

function run_affine_cells(input,output)
    abspath(input) != abspath(output) || error("Receipt would overwrite input")
    receipt = Dict{String,Any}(
        "complete"=>false, "time_utc"=>string(now(UTC)), "julia"=>string(VERSION),
        "arithmetic"=>"Rational{BigInt}",
        "source_sha256"=>Dict("tests/"*name => bytes2hex(sha256(read(joinpath(@__DIR__,name))))
                             for name in ["verify_affine_cells.jl","verify_witnesses.jl"]),
        "scope"=>"Cell-attaining witnesses and exact pointwise supports; not a universal upper-bound or asymptotic proof")
    mkpath(dirname(output))
    open(io -> TOML.print(io,receipt),output,"w")
    try
        receipt["input_sha256"] = bytes2hex(sha256(read(input)))
        data = TOML.parsefile(input)
        get(data,"complete",false) || error("Incomplete input")
        data["schema"] == "affine-cell-controls-v1" || error("Wrong input schema")
        cases = data["witnesses"]
        isempty(cases) && error("No witnesses")
        gaps = cell_verify.(cases)
        receipt["witnesses_verified"] = length(gaps)
        receipt["gaps"] = string.(gaps)
        receipt["malformed_inputs_rejected"] = cell_negative_checks(first(cases))
        receipt["complete"] = true
    catch exc
        receipt["complete"] = false
        receipt["error"] = sprint(showerror,exc)
        rethrow()
    finally
        open(io -> TOML.print(io,receipt),output,"w")
    end
    println("VERIFIED ",receipt["witnesses_verified"]," cell witnesses and eight malformed controls")
    receipt
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    directory = normpath(joinpath(@__DIR__,"..","results","theory"))
    run_affine_cells(joinpath(directory,"affine_cells.toml"),joinpath(directory,"affine_cells_julia.toml"))
end
