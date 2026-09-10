# One Julia process for all positive, negative and failure-path checks.
include("check_verifier_failures.jl")
for (input, output) in [
    ("exact_checks", "julia_verification"),
    ("affine_profiles", "affine_profiles_julia"),
    ("partition_reduction", "partition_reduction_julia"),
    ("affine_dimensions", "affine_dimensions_julia"),
]
    run_verification(joinpath(root,"results","theory",input*".toml"),
                     joinpath(root,"results","theory",output*".toml"))
end
include("verify_affine_upper.jl")
verify_affine_upper(joinpath(root,"results","theory","affine_upper.toml"),
                    joinpath(root,"results","theory","affine_upper_julia.toml"))
include("verify_affine_cells.jl")
run_affine_cells(joinpath(root,"results","theory","affine_cells.toml"),
                    joinpath(root,"results","theory","affine_cells_julia.toml"))
include("verify_affine_two_level.jl")
run_affine_two_level(joinpath(root,"results","theory","affine_two_level.toml"),
                        joinpath(root,"results","theory","affine_two_level_julia.toml"))
include("verify_affine_certified.jl")
verify_affine_certified(joinpath(root,"results","theory","affine_certified.toml"),
                       joinpath(root,"results","theory","affine_certified_julia.toml"))
include("check_affine_certified_failures.jl")
check_affine_certified_failures()
include("verify_affine_boundary.jl")
verify_affine_boundary(joinpath(root,"results","theory","affine_boundary.toml"),
                      joinpath(root,"results","theory","affine_boundary_julia.toml"))
