# Exercise the actual D9 wrapper, including failure after mathematical checks.
if !isdefined(@__MODULE__, :verify_affine_certified)
    include("verify_affine_certified.jl")
end

function check_affine_certified_failures()
    root = normpath(joinpath(@__DIR__,".."))
    destination = joinpath(root,"results","theory","affine_certified_failures.toml")
    write(destination,"complete = false\n")
    valid = joinpath(root,"results","theory","affine_certified.toml")
    passed = String[]
    mktempdir() do scratch
        malformed = joinpath(scratch,"malformed.toml")
        write(malformed,"[invalid\n")
        bad_certificate = joinpath(scratch,"bad-certificate.toml")
        data = TOML.parsefile(valid)
        data["controls"][1]["upper"] = "0"
        open(io->TOML.print(io,data),bad_certificate,"w")
        closed_output = IOBuffer()
        close(closed_output)
        for (name,input,status_io) in [
            ("missing_input",joinpath(scratch,"absent.toml"),devnull),
            ("malformed_toml",malformed,devnull),
            ("invalid_certificate",bad_certificate,devnull),
            ("reporting_after_all_checks",valid,closed_output),
        ]
            output = joinpath(scratch,name*".toml")
            write(output,"complete = true\n")
            failed = false
            try
                verify_affine_certified(input,output; status_io=status_io)
            catch
                failed = true
            end
            failed || error("Failure injection did not fail: "*name)
            receipt = TOML.parsefile(output)
            receipt["complete"] == false || error("Stale successful receipt: "*name)
            haskey(receipt,"error") || error("Missing failure diagnostic: "*name)
            if name == "reporting_after_all_checks"
                receipt["scalar_intervals_verified"] == 21052 || error("Failure was too early")
                receipt["malformed_inputs_rejected"] == 19 || error("Negative controls did not finish")
            end
            push!(passed,name)
        end
    end
    receipt = Dict("complete"=>true,"julia"=>string(VERSION),"checks_passed"=>passed,
        "input_sha256"=>bytes2hex(sha256(read(valid))),
        "source_sha256"=>Dict("tests/"*name=>bytes2hex(sha256(read(joinpath(@__DIR__,name))))
            for name in ["check_affine_certified_failures.jl","verify_affine_certified.jl",
                         "verify_affine_cells.jl","verify_witnesses.jl"]))
    open(io->TOML.print(io,receipt),destination,"w")
    println("Verified ",length(passed)," D9 failure paths leave incomplete receipts")
    receipt
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    check_affine_certified_failures()
end
