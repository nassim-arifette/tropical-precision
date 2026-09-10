# Exact independent checks of D10 algebra and root intervals, with D8/D9
# geometry/scalar certificates as controls. The universal theorem is written.
if !isdefined(@__MODULE__, :certified_partition)
    include("verify_affine_certified.jl")
end

boundary_B(d,c)=(2(d+1)-(5d+4)*c+2d*c^2)/(2(d+1)*(1-c)^2)
function boundary_F(d,c)
    c==0 || c==1//BigInt(2) ? 0//BigInt(1) :
        c/(1+(1-c)/((d+2c)/(d+2))+c/(1-2c))
end
boundary_P(d,q,c)=(1-2c)*(2(d+1)-d*c+(q-2)*c*(d+2c))-(d+1)*(1-c)^2

function boundary_bracket(record,minimum)
    r=vector(record["root_interval"])
    length(r)==2 || error("Root bracket dimension")
    lo,hi=r
    minimum<=lo<=hi<=1//BigInt(2) || error("Root bracket domain")
    hi-lo<=1//(BigInt(2)^80) || error("Root bracket tolerance")
    lo,hi
end

function check_boundary_root(record)
    n=record["dimension"]
    n isa Integer && n>=3 || error("Boundary dimension")
    d=BigInt(n)-3
    b=rational(record["endpoint"])
    0<=b<=1//BigInt(2) || error("Boundary endpoint")
    lo,hi=boundary_bracket(record,d==0 ? 2//BigInt(5) : 3//BigInt(8))
    if b==0
        lo==hi==1//BigInt(2) || error("Zero-cell root")
    else
        boundary_B(d,lo)>=b>=boundary_B(d,hi) || error("Boundary root signs")
    end
    balanced=d//BigInt(3d+1)*b
    expected=[max(balanced,boundary_F(d,hi)),max(balanced,boundary_F(d,lo))]
    vector(record["loss_interval"])==expected || error("Boundary loss formula")
    length(record["certificates"])==1 || error("Missing boundary control")
    cert=only(record["certificates"])
    cert["dimension"]==n && vector(cert["interval"])==[0,b] || error("Wrong boundary control")
    leaves=certified_cell(cert)
    rational(cert["lower"])<=last(expected) && first(expected)<=rational(cert["upper"]) || error("Boundary formula refuted by D9")
    leaves
end

function check_budget_root(record)
    n,M=record["dimension"],record["generator_budget"]
    n isa Integer && n>=4 && M isa Integer && M>=4 || error("Budget/dimension domain")
    d,q=BigInt(n)-3,BigInt(M)-2
    lo,hi=boundary_bracket(record,3//BigInt(8))
    boundary_P(d,q,lo)>=0>=boundary_P(d,q,hi) || error("Cubic root signs")
    balanced=1//BigInt(2)/(M-1+1//BigInt(d))
    rational(record["balanced_value"])==balanced || error("Balanced law")
    lower,upper=max(balanced,boundary_F(d,hi)),max(balanced,boundary_F(d,lo))
    vector(record["loss_interval"])==[lower,upper] || error("Budget loss interval")
    vector(record["boundary_interval"])==[1//BigInt(2)-(q-2)*upper,1//BigInt(2)-(q-2)*lower] || error("Boundary widths")
    vector(record["second_order_interval"])==[BigInt(M)^2*(lower-1//(2BigInt(M))),BigInt(M)^2*(upper-1//(2BigInt(M)))] || error("Second-order values")
    leaves=0
    length(record["certificates"])<=1 || error("Ambiguous partition controls")
    for cert in record["certificates"]
        cert["dimension"]==n && cert["generator_budget"]==M || error("Wrong partition control")
        leaves+=certified_partition(cert)
        rational(cert["lower"])<=upper && lower<=rational(cert["upper"]) || error("Budget formula refuted by D9")
    end
    leaves
end

function check_transition_root(record)
    n=record["dimension"]
    n isa Integer && n>=4 || error("Transition dimension")
    d=BigInt(n)-3
    lo,hi=boundary_bracket(record,0//BigInt(1))
    switch(c)=(3d+1)*c*(d+2c)-d*(2(d+1)-d*c)
    switch(lo)<=0<=switch(hi) || error("Transition root signs")
    blo,bhi=boundary_B(d,hi),boundary_B(d,lo)
    vector(record["boundary_interval"])==[blo,bhi] || error("Transition boundary")
    kap=d//BigInt(3d+1)
    mlo,mhi=4+(1/(2bhi)-1)/kap,4+(1/(2blo)-1)/kap
    vector(record["budget_threshold_interval"])==[mlo,mhi] || error("Transition budget")
    floor(BigInt,mlo)==floor(BigInt,mhi)==record["last_balanced_budget"] || error("Integer transition")
end

function boundary_algebra(data)
    names=["k_derivative","B_derivative","F_derivative","cubic_expansion",
           "boundary_quadratic","switch_polynomial","endpoint_gap",
           "dimension_four_threshold","dimension_four_M9_root","dimension_four_M9_value"]
    degrees=[[2,2],[1,2],[1,3],[1,1,3],[1,1,2],[2,2],[2],[2],[3],[3]]
    residuals=[
        (d,c)->(d+4c)*(2d+2-d*c)+d*(d*c+2c^2)-(2d*(d+1)+8(d+1)*c-2d*c^2),
        (d,c)->(-5d-4+4d*c)*(1-c)+2*(2d+2-(5d+4)*c+2d*c^2)+d+(d+4)*c,
        (d,c)->(d+4c-4d*c-12c^2)*(1-c)+2*(d*c+2c^2)*(1-2c)-d*(1-3c)-4c*(1-3c+c^2),
        (d,q,c)->boundary_P(d,q,c)-(d+1+((q-5)*d-2)*c+(2q-5)*(1-d)*c^2-4*(q-2)*c^3),
        (d,b,c)->2*(d+1)*b*(1-c)^2-(1-2c)*(2*(d+1)-d*c)-(2*(b*(d+1)-d)*c^2+(5d+4-4b*(d+1))*c+2*(d+1)*(b-1)),
        (d,c)->(3d+1)*c*(d+2c)-d*(2*(d+1)-d*c)-(2*(3d+1)*c^2+d*(4d+1)*c-2d*(d+1)),
        d->(d+1)*(3d+1)-d*(3d+4)-1,
        c->(13//BigInt(2)+4c)*(20-41c)-(48-84c)+(41//BigInt(2))*(8c^2+5c-4),
        c->boundary_P(1,7,c)-(2-20c^3),
        c->c*(1+2c)*(1-2c)-(c-2//BigInt(5))+4*(c^3-1//BigInt(10)),
    ]
    length(data)==length(names) || error("Algebra inventory")
    total=0
    for (record,name,degree,residual) in zip(data,names,degrees,residuals)
        record["name"]==name && record["degree_bounds"]==degree || error("Identity/degree mismatch")
        axes=[collect(-BigInt(d):BigInt(0)).//BigInt(7) for d in degree]
        count=0
        for node in Iterators.product(axes...)
            residual(node...)==0 || error("Nonzero identity: "*name)
            count+=1
        end
        record["evaluations"]==count || error("Algebra evaluation count")
        total+=count
    end
    total
end

function boundary_negative_checks(data)
    count=0
    function rejects(verify,original,change)
        bad=deepcopy(original)
        change(bad)
        try
            verify(bad)
        catch
            return
        end
        error("Accepted malformed D10 input")
    end
    for change in [c->(c["dimension"]=2),c->(c["endpoint"]="1"),
                   c->(c["root_interval"][1]="0"),c->(c["loss_interval"][2]="0"),
                   c->(c["certificates"][1]["dimension"]=4),c->empty!(c["certificates"])]
        rejects(check_boundary_root,first(data["boundaries"]),change);count+=1
    end
    for change in [c->(c["dimension"]=3),c->(c["generator_budget"]=3),
                   c->(c["root_interval"][1]="0"),c->(c["loss_interval"][1]="0"),
                   c->(c["balanced_value"]="0"),c->(c["boundary_interval"][1]="0"),
                   c->(c["second_order_interval"][1]="99"),
                   c->(c["certificates"][1]["generator_budget"]+=1)]
        rejects(check_budget_root,first(data["budgets"]),change);count+=1
    end
    for change in [c->(c["dimension"]=3),c->(c["root_interval"][1]="0"),
                   c->(c["budget_threshold_interval"][1]="0"),c->(c["last_balanced_budget"]+=1)]
        rejects(check_transition_root,first(data["transitions"]),change);count+=1
    end
    for change in [c->(c[1]["degree_bounds"][1]=1),c->(c[1]["name"]="unknown")]
        rejects(boundary_algebra,data["algebra"],change);count+=1
    end
    count
end

function verify_affine_boundary(input,output; status_io=stdout)
    abspath(input)!=abspath(output) || error("Receipt would overwrite input")
    receipt=Dict{String,Any}("complete"=>false,"julia"=>string(VERSION),"time_utc"=>string(now(UTC)),
        "arithmetic"=>"Rational{BigInt}",
        "source_sha256"=>Dict("tests/"*name=>bytes2hex(sha256(read(joinpath(@__DIR__,name))))
            for name in ["verify_affine_boundary.jl","verify_affine_certified.jl","verify_affine_cells.jl","verify_witnesses.jl"]),
        "scope"=>"D10 algebra/root checks and independent D8/D9 control certificates; universal theorem is not mechanized")
    mkpath(dirname(output))
    open(io->TOML.print(io,receipt),output,"w")
    try
        receipt["input_sha256"]=bytes2hex(sha256(read(input)))
        data=TOML.parsefile(input)
        get(data,"complete",false) && data["schema"]=="affine-boundary-budget-v1" || error("Incomplete/wrong source")
        all(family->!isempty(data[family]),["boundaries","budgets","transitions"]) || error("Missing root family")
        receipt["algebra_evaluations"]=boundary_algebra(data["algebra"])
        leaves=sum(check_boundary_root.(data["boundaries"]))+sum(check_budget_root.(data["budgets"]))
        check_transition_root.(data["transitions"])
        receipt["boundary_cells_verified"]=length(data["boundaries"])
        receipt["budget_roots_verified"]=length(data["budgets"])
        receipt["transition_roots_verified"]=length(data["transitions"])
        receipt["global_partition_controls_verified"]=sum(length(r["certificates"]) for r in data["budgets"])
        receipt["scalar_intervals_verified"]=leaves
        receipt["malformed_inputs_rejected"]=boundary_negative_checks(data)
        println(status_io,"VERIFIED D10: ",receipt["boundary_cells_verified"]," boundary cells / ",
                receipt["budget_roots_verified"]," budgets / ",leaves," D9 scalar intervals")
        receipt["complete"]=true
    catch exc
        receipt["complete"]=false
        receipt["error"]=sprint(showerror,exc)
        rethrow()
    finally
        open(io->TOML.print(io,receipt),output,"w")
    end
    receipt
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    directory=normpath(joinpath(@__DIR__,"..","results","theory"))
    verify_affine_boundary(joinpath(directory,"affine_boundary.toml"),joinpath(directory,"affine_boundary_julia.toml"))
end
