# Independent rational polynomial checks for the sharp three-generator theorem.
# A polynomial with coordinate degrees d_i is identically zero if it vanishes
# on a tensor product of d_i+1 distinct nodes. No floating arithmetic is used.
using TOML, SHA, Dates

function upper_algebra(data)
    get(data,"complete",false) || error("Incomplete source")
    data["schema"] == "affine-upper-algebra-v1" || error("Wrong schema")
    names = ["weighted_denominators", "mixed_cross_multiplication",
             "mixed_square_decomposition", "positive_remainder",
             "derivative_numerator", "golden_value"]
    degrees = [[3,3,3], [2,3], [2,4], [2], [3], [3]]
    data["identity_names"] == names || error("Identity inventory")
    data["degree_bounds"] == degrees || error("Polynomial degree bounds")
    data["mixed_remainder_lower_bound"] == "7/128" || error("Positive remainder")
    data["optimum_coefficients"] == ["-3","5"] || error("Algebraic optimum")
    residuals = [
        (a,b,c) -> c*((a+b)*c+(a-b)^2) -
                   (a+b)*(c^2-(a-b)^2) - (a+b+c)*(a-b)^2,
        (a,k) -> (1-k)^2*a*(1-2a+k) - (1+k)^2*(a-k)*(1-2a) -
                   k*(8a^2-(k^2+6k+5)*a+(1+k)^2),
        (a,k) -> (8a^2-(k^2+6k+5)*a+(1+k)^2) -
                   8*(a-(k^2+6k+5)/16)^2 -
                   (1+k)^2*(1-(5+k)^2/32),
        k -> 1-(5+k)^2/32-7//BigInt(128) -
                   (1//BigInt(2)-k)*(21//BigInt(2)+k)/32,
        r -> (2r-3r^2)*(1+r)-(r^2-r^3)-2r*(1-r-r^2),
        r -> r^2-r^3-(5r-3)*(1+r)+(r+3)*(r^2+r-1),
    ]
    counts = Int[]
    for (degree,residual) in zip(degrees,residuals)
        # Negative nodes differ from the Python implementation's positive grid.
        axes = [collect(-BigInt(d):BigInt(0)) .// BigInt(1) for d in degree]
        count = 0
        for node in Iterators.product(axes...)
            residual(node...) == 0 || error("Nonzero polynomial residual")
            count += 1
        end
        push!(counts,count)
    end
    counts == data["evaluations"] || error("Unisolvent evaluation count")
    return sum(counts)
end

function upper_negative_checks(original)
    mutations = [
        d -> (d["complete"] = false),
        d -> pop!(d["identity_names"]),
        d -> (d["degree_bounds"][1][1] = 2),
        d -> (d["mixed_remainder_lower_bound"] = "8/128"),
        d -> (d["optimum_coefficients"][1] = "-2"),
    ]
    for change in mutations
        bad = deepcopy(original)
        change(bad)
        failed = false
        try
            upper_algebra(bad)
        catch
            failed = true
        end
        failed || error("Accepted malformed algebra receipt")
    end
    length(mutations)
end

function verify_affine_upper(input, output)
    abspath(input) != abspath(output) || error("Receipt would overwrite input")
    receipt = Dict{String,Any}(
        "complete"=>false, "julia"=>string(VERSION), "time_utc"=>string(now(UTC)),
        "arithmetic"=>"Rational{BigInt}",
        "verifier_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "scope"=>"Six polynomial identities by exact unisolvence; not the whole geometric proof")
    mkpath(dirname(output))
    open(io -> TOML.print(io,receipt),output,"w")
    try
        receipt["input_sha256"] = bytes2hex(sha256(read(input)))
        data = TOML.parsefile(input)
        receipt["polynomial_evaluations"] = upper_algebra(data)
        receipt["identities_verified"] = 6
        receipt["malformed_inputs_rejected"] = upper_negative_checks(data)
        receipt["complete"] = true
    catch exc
        receipt["complete"] = false
        receipt["error"] = sprint(showerror,exc)
        rethrow()
    finally
        open(io -> TOML.print(io,receipt),output,"w")
    end
    println("VERIFIED six affine upper-bound identities and five malformed controls")
    receipt
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    directory = normpath(joinpath(@__DIR__,"..","results","theory"))
    verify_affine_upper(joinpath(directory,"affine_upper.toml"),
                        joinpath(directory,"affine_upper_julia.toml"))
end
