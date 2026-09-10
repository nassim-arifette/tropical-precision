# Independently check original/reduced deficits and explicit contextual witnesses.
# Does not prove the universal reduction or perform scalar global optimization.
if !isdefined(@__MODULE__, :cells_deficit)
    include("verify_affine_cells.jl")
end

function expected_reduction_branch(x,w)
    A = sum(v for v in w if v > 0)
    C = -sum(v for v in w if v < 0)
    A > 0 && C > 0 || error("Both signs required")
    u = sum(v*t for (t,v) in zip(x,w) if v > 0)/A
    v = -sum(v*t for (t,v) in zip(x,w) if v < 0)/C
    u < v || error("Sign averages not ordered")
    p,k = count(>(0),w),count(<(0),w)
    if v <= 1//2
        side = "left_"
    elseif u >= 1//2
        side = "right_"
        p,k = k,p
    else
        return "central"
    end
    side * (p == 1 ? "positive_singleton" : k == 1 ? "negative_singleton" : "both")
end

function two_level_verify(case)
    n = case["dimension"]
    n >= 3 || error("Dimension")
    interval = vector(case["interval"])
    length(interval) == 2 || error("Interval length")
    a,b = interval
    0 <= a < b <= 1 || error("Cell interval")
    original = vector(case["original_target"])
    weights = vector(case["original_weights"])
    x,w = vector(case["target"]),vector(case["affine_coefficients"])
    length(original) == length(weights) == length(x) == length(w) == n || error("Dimensions")
    all(t -> a <= t <= b,original) && all(t -> a <= t <= b,x) || error("Point outside cell")
    0 < sum(abs.(weights)) <= 1 || error("Original normalization")
    sum(abs.(w)) == 1 && all(!iszero,w) || error("Reduced normalization")
    length(unique(x)) == 2 || error("Not two levels")
    for sign in (-1,1)
        group = [w[j] for j in eachindex(w) if sign*w[j] > 0]
        !isempty(group) && length(unique(group)) == 1 || error("Sign weights not equal")
    end
    all(x[j] == (w[j] > 0 ? minimum(x) : maximum(x)) for j in eachindex(x)) ||
        error("Sign/level alignment")
    source = cells_deficit(original,weights)
    source > 0 && source == rational(case["original_deficit"]) || error("Original deficit")
    case["branch"] == expected_reduction_branch(original,weights) || error("Reduction branch")
    delta = cells_deficit(x,w)
    delta == rational(case["certified_gap"]) && delta >= source || error("Domination")

    G,T = vector.(case["generators"]),vector.(case["context"])
    length(G) == case["generator_count"] || error("Generator count")
    all(g -> length(g) == 2n,G) || error("Generator dimension")
    all(t -> length(t) == n && all(v -> 0 <= v <= 1,t),T) || error("Context invariant")
    check_slice(G,n)
    partition = sort(unique([0//BigInt(1),a,b,1//BigInt(1)]))
    expected = [lift(fill(0//BigInt(1),n)),lift(fill(1//BigInt(1),n))]
    append!(expected,[vcat(fill(partition[i],n),fill(-partition[i+1],n)) for i in 1:length(partition)-1])
    G == expected || error("Canonical cell representative")
    check_combination(vcat(G,lift.(T)),lift(x),vector(case["coefficients"]))
    concrete = minimum(vcat([min(0,sum(w))],[dotq(w,t) for t in T]))
    value = dotq(w,x)
    concrete == rational(case["minimum_concrete_union"]) || error("Concrete minimum")
    value == rational(case["target_value"]) || error("Target value")
    concrete-value == delta || error("Contextual gap")
    delta
end

function two_level_negative_checks(original)
    mutations = [
        c -> (c["original_deficit"] = "99"),
        c -> (c["original_target"][1] = "-1"),
        c -> (c["original_weights"][1] = "2"),
        c -> (c["branch"] = "unknown"),
        c -> (c["coefficients"][1] = "1/2"),
        c -> (c["context"][1][1] = "2"),
        c -> (c["affine_coefficients"][1] = "2"),
        c -> (c["certified_gap"] = "99"),
        c -> (c["target"][1] = "7/13"),
        c -> (c["interval"][2] = "0"),
        c -> popfirst!(c["generators"]),
        c -> (c["generator_count"] += 1),
    ]
    for mutation in mutations
        bad = deepcopy(original)
        mutation(bad)
        rejected = false
        try
            two_level_verify(bad)
        catch
            rejected = true
        end
        rejected || error("Accepted a malformed reduction witness")
    end
    length(mutations)
end

function run_affine_two_level(input,output)
    abspath(input) != abspath(output) || error("Receipt would overwrite input")
    receipt = Dict{String,Any}(
        "complete"=>false, "time_utc"=>string(now(UTC)), "julia"=>string(VERSION),
        "arithmetic"=>"Rational{BigInt}",
        "source_sha256"=>Dict("tests/"*name => bytes2hex(sha256(read(joinpath(@__DIR__,name))))
                             for name in ["verify_affine_two_level.jl","verify_affine_cells.jl","verify_witnesses.jl"]),
        "scope"=>"Exact source/reduced deficits, domination and contextual witnesses; not a universal proof or scalar global optimization")
    mkpath(dirname(output))
    open(io -> TOML.print(io,receipt),output,"w")
    try
        receipt["input_sha256"] = bytes2hex(sha256(read(input)))
        data = TOML.parsefile(input)
        get(data,"complete",false) || error("Incomplete input")
        data["schema"] == "affine-two-level-controls-v1" || error("Input schema")
        cases = data["witnesses"]
        isempty(cases) && error("No witnesses")
        gaps = two_level_verify.(cases)
        receipt["witnesses_verified"] = length(gaps)
        receipt["gaps"] = string.(gaps)
        receipt["malformed_inputs_rejected"] = two_level_negative_checks(first(cases))
        receipt["complete"] = true
    catch exc
        receipt["complete"] = false
        receipt["error"] = sprint(showerror,exc)
        rethrow()
    finally
        open(io -> TOML.print(io,receipt),output,"w")
    end
    println("VERIFIED ",receipt["witnesses_verified"]," reduction witnesses and twelve malformed controls")
    receipt
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    directory = normpath(joinpath(@__DIR__,"..","results","theory"))
    run_affine_two_level(joinpath(directory,"affine_two_level.toml"),joinpath(directory,"affine_two_level_julia.toml"))
end
