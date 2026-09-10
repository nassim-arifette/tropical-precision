# Exercise failure paths through the actual verifier, including an error after
# all positive checks. No production input or receipt is modified.
include("verify_witnesses.jl")

root = normpath(joinpath(@__DIR__,".."))
destination = joinpath(root,"results","theory","verifier_failures.toml")
write(destination,"complete = false\n")
scratch = joinpath(root,".build","verifier-failure-checks")
mkpath(scratch)
valid = joinpath(root,"results","theory","exact_checks.toml")
bad_toml = joinpath(scratch,"malformed.toml")
write(bad_toml,"[invalid\n")
bad_witness = joinpath(scratch,"bad-witness.toml")
data = TOML.parsefile(valid)
data["witnesses"][1]["coefficients"][1] = "1/2"
open(io -> TOML.print(io,data),bad_witness,"w")
conditions = [
    ("missing_input", joinpath(scratch,"absent-input.toml"), nothing),
    ("malformed_toml", bad_toml, nothing),
    ("invalid_certificate", bad_witness, nothing),
    ("failed_negative_control", valid, c -> error("Injected negative-control failure")),
]
passed = String[]
for (name,input,negative) in conditions
    output = joinpath(scratch,name*".toml")
    write(output,"complete = true\n")
    failed = false
    try
        run_verification(input,output; negative_check=negative)
    catch
        failed = true
    end
    failed || error("Failure injection did not fail: "*name)
    result = TOML.parsefile(output)
    result["complete"] == false || error("Successful stale receipt: "*name)
    haskey(result,"error") || error("Missing diagnostic: "*name)
    push!(passed,name)
end
receipt = Dict(
    "complete"=>true, "julia"=>string(VERSION), "checks_passed"=>passed,
    "verifier_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"verify_witnesses.jl")))),
    "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
    "input_sha256"=>bytes2hex(sha256(read(valid))),
)
open(io -> TOML.print(io,receipt),destination,"w")
println("Verified ",length(passed)," failure paths leave incomplete receipts")
