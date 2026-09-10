# Independent certificate verifier: standard-library TOML/SHA only.
# It checks rational equations and an interval-cover criterion for the complete
# diagonal slice; it does not recompute Python's sector/context optimization.
using TOML, SHA, Dates

function rational(s)
    parts = split(string(s), '/')
    length(parts) == 1 && return parse(BigInt, parts[1]) // BigInt(1)
    length(parts) == 2 || error("Malformed rational")
    parse(BigInt, parts[1]) // parse(BigInt, parts[2])
end
vector(xs) = rational.(xs)
dotq(a, b) = sum(a .* b)
lift(x) = vcat(x, -x)

function check_slice(G, n)
    # For this certificate class, equal coordinate groups force every consistent
    # point onto the diagonal. Bounds and endpoints restrict it to [0,1].
    all(g -> all(==(g[1]), g[1:n]) && all(==(g[n+1]), g[n+1:2n]), G) ||
        error("Unsupported nonuniform generator: no slice certificate supplied")
    all(g -> 0 <= g[1] <= 1 && -1 <= g[n+1] <= 0, G) || error("Invalid bounds")
    lift(zeros(Rational{BigInt}, n)) in G || error("Missing zero endpoint")
    lift(ones(Rational{BigInt}, n)) in G || error("Missing one endpoint")
    intervals = sort([(g[1], -g[n+1]) for g in G if g[1] <= -g[n+1]])
    reached = BigInt(0)//BigInt(1)
    for (a,b) in intervals
        a > reached && error("Gap in diagonal anchor coverage")
        reached = max(reached, b)
    end
    reached == 1 || error("Incomplete coverage")
    # Any covered t is reconstructed by its anchor, E(0)-t, E(1)-(1-t).
end

function verify(case)
    n = case["dimension"]
    G = vector.(case["generators"])
    T = vector.(case["context"])
    z = vector(case["target"])
    w = vector(case["affine_coefficients"])
    coeff = vector(case["coefficients"])
    length(G) == case["generator_count"] || error("Generator count")
    n >= 3 || error("Wrong dimensional claim")
    all(g -> length(g) == 2n, G) || error("Generator dimension")
    all(t -> length(t) == n && all(v -> 0 <= v <= 1, t), T) || error("Context box")
    length(w) == length(z) == n || error("Vector dimension")
    sum(abs.(w)) == 1 || error("Wrong Lipschitz normalization")
    check_slice(G,n)
    operands = vcat(G, lift.(T))
    length(coeff) == length(operands) || error("Coefficient count")
    maximum(coeff) == 0 || error("Normalization")
    reconstructed = [maximum(g[j]+c for (g,c) in zip(operands,coeff)) for j in 1:2n]
    reconstructed == lift(z) || error("Invalid reconstruction")
    true_minimum = minimum(vcat([min(0,sum(w))], [dotq(w,t) for t in T]))
    target_value = dotq(w,z)
    gap = true_minimum-target_value
    true_minimum == rational(case["minimum_concrete_union"]) || error("Concrete minimum")
    target_value == rational(case["target_value"]) || error("Target value")
    gap == rational(case["certified_gap"]) || error("Gap")
    gap >= 1//(32*(length(G)-2)) || error("Claimed rate lower bound")
    return gap
end

function rejection_checks(original)
    mutations = [
        c -> (c["coefficients"][1] = "1/2"),
        c -> (c["context"][1][1] = "2"),
        c -> (c["affine_coefficients"][1] = "2"),
        c -> (c["target"][1] = "7/13"),
        c -> (c["certified_gap"] = "99"),
        c -> (c["generator_count"] += 1),
        c -> (c["generators"] = c["generators"][1:2]),
    ]
    for change in mutations
        bad = deepcopy(original)
        change(bad)
        rejected = false
        try
            verify(bad)
        catch
            rejected = true
        end
        rejected || error("Accepted a deliberately malformed witness")
    end
    return length(mutations)
end

function check_combination(G, target, coefficients)
    length(coefficients) == length(G) || error("Coefficient count")
    maximum(coefficients) == 0 || error("Combination normalization")
    all(g -> length(g) == length(target), G) || error("Operand dimension")
    [maximum(g[j]+c for (g,c) in zip(G,coefficients))
     for j in eachindex(target)] == target || error("Containment certificate")
end

function verify_partition(case)
    n = case["dimension"]
    n >= 1 || error("Invalid dimension")
    G = vector.(case["generators"])
    partition = vector(case["partition"])
    output = vector.(case["output_generators"])
    coefficients = vector.(case["coefficients"])
    all(g -> length(g) == 2n, G) || error("Generator dimension")
    length(G) >= 3 || error("Endpoint generator lower bound")
    first(partition) == 0 && last(partition) == 1 || error("Partition endpoints")
    all(diff(partition) .> 0) || error("Partition must be strict")
    q = length(partition)-1
    q >= 1 || error("Empty partition")
    length(output) == q+2 <= length(G) || error("Generator budget")
    length(coefficients) == length(output) || error("Containment count")
    expected = [lift(fill(0//BigInt(1),n)), lift(fill(1//BigInt(1),n))]
    append!(expected, [vcat(fill(partition[i],n), fill(-partition[i+1],n)) for i in 1:q])
    output == expected || error("Output is not canonical")
    # These exact equalities certify Q subset P for every ambient point.
    for (g, coeff) in zip(output, coefficients)
        check_combination(G,g,coeff)
    end
    # Canonical bounds and coordinate groups, together with strict full coverage,
    # certify the entire output slice, not just membership of sampled points.
    indices = case["source_indices"]
    length(indices) == q || error("Source interval count")
    for i in 1:q
        k = indices[i]+1  # Python export uses zero-based source indices.
        1 <= k <= length(G) || error("Source index")
        all(G[k] .<= output[i+2]) || error("Cell is outside its anchor interval")
    end
    anchors = [(maximum(g[1:n]), minimum(-g[n+1:2n])) for g in G]
    count(ab -> max(0,ab[1]) < min(1,ab[2]), anchors) <= length(G)-2 ||
        error("Endpoint interval count")
    scope = case["input_scope"]
    if scope == "exact_one_uniform_block"
        (all(g -> all(==(g[1]),g[1:n]),G) ||
         all(g -> all(==(g[n+1]),g[n+1:2n]),G)) ||
            error("Input exact-slice class is unsupported")
        lo = -minimum(maximum(g[n+j] for g in G) for j in 1:n)
        hi = minimum(maximum(g[j] for g in G) for j in 1:n)
        for (a,b) in anchors
            l,u = max(a,lo),min(b,hi)
            (l > u || 0 <= l <= u <= 1) || error("Input slice exceeds diagonal segment")
        end
    elseif scope == "inexact_enclosure"
        point = vector(case["inexact_point"])
        length(point) == n || error("Inexact point dimension")
        !(all(==(point[1]),point) && 0 <= point[1] <= 1) ||
            error("Purported inexact point belongs to diagonal")
        check_combination(G,lift(point),vector(case["inexact_coefficients"]))
    else
        error("Unknown input scope")
    end
    return length(output)
end

function partition_rejection_checks(original)
    mutations = [
        c -> (c["coefficients"][1][1] = "1/2"),
        c -> (c["output_generators"][3][1] = "1/7"),
        c -> (c["partition"][end] = "3/4"),
        c -> (c["source_indices"][1] = length(c["generators"])),
        c -> pop!(c["output_generators"]),
        c -> (c["input_scope"] = "unrecognized"),
        c -> (c["generators"][end][end] = "2"),
    ]
    for change in mutations
        bad = deepcopy(original)
        change(bad)
        rejected = false
        try
            verify_partition(bad)
        catch
            rejected = true
        end
        rejected || error("Accepted a malformed partition certificate")
    end
    length(mutations)
end

function run_verification(input, output; negative_check=nothing)
    abspath(input) != abspath(output) || error("Receipt must not overwrite its input")
    receipt = Dict{String,Any}(
        "complete"=>false, "time_utc"=>string(now(UTC)), "julia"=>string(VERSION),
        "arithmetic"=>"Rational{BigInt}",
        "verifier_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "scope"=>"Finite rational containment/slice certificates, not universal theorem verification")
    mkpath(dirname(output))
    # Clear a former successful receipt before input I/O or mathematical work.
    open(io -> TOML.print(io,receipt), output, "w")
    try
        receipt["input_sha256"] = bytes2hex(sha256(read(input)))
        data = TOML.parsefile(input)
        get(data,"complete",false) || error("Source export is incomplete")
        if data["schema"] == "bounded-context-witnesses-v1"
            cases = data["witnesses"]
            gaps = verify.(cases)
            receipt["witnesses_verified"] = length(gaps)
            receipt["gaps"] = string.(gaps)
            check = isnothing(negative_check) ? rejection_checks : negative_check
        elseif data["schema"] == "diagonal-partition-certificates-v1"
            cases = data["cases"]
            counts = verify_partition.(cases)
            receipt["partitions_verified"] = length(counts)
            receipt["generator_certificates_verified"] = sum(counts)
            check = isnothing(negative_check) ? partition_rejection_checks : negative_check
        else
            error("Unknown certificate schema")
        end
        receipt["malformed_witnesses_rejected"] = check(first(cases))
        receipt["complete"] = true
        println("VERIFIED ",length(cases)," exact certificate cases")
    catch exc
        receipt["complete"] = false
        receipt["error"] = sprint(showerror,exc)
        rethrow()
    finally
        open(io -> TOML.print(io,receipt), output, "w")
    end
    return receipt
end

function main()
    root = normpath(joinpath(@__DIR__, ".."))
    input = length(ARGS) >= 1 ? abspath(ARGS[1]) : joinpath(root,"results","theory","exact_checks.toml")
    output = length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(root,"results","theory","julia_verification.toml")
    run_verification(input,output)
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
