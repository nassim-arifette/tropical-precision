# Exact global scalar covers and primal/dual partition certificates.
# Mathematical guarantees are conditional on the written D8/D9 theorems.
if !isdefined(@__MODULE__, :cells_deficit)
    include("verify_affine_cells.jl")
end

function certified_specs(n,a,b)
    specs=NamedTuple[]
    constants=Rational{BigInt}[]
    half=1//BigInt(2)
    for (aa,bb) in [(a,min(b,half)),(1-b,min(1-a,half))]
        aa>=bb && continue
        for kind in [:positive,:negative]
            push!(specs,(kind=kind,a=aa,b=bb,k=0,lo=0//BigInt(1),hi=half))
        end
        n>=4 && push!(constants,min((bb-aa)/2,(n-3)*bb/(3n-8)))
    end
    if a<half<b
        for k in 1:n-1
            p=n-k
            push!(specs,(kind=:central,a=a,b=b,k=k,
                         lo=max(0//BigInt(1),1-k//BigInt(2)),
                         hi=min(1//BigInt(1),p//BigInt(2))))
        end
    end
    specs,constants
end

function certified_value(n,s,t)
    s.lo<=t<=s.hi || error("Parameter outside branch")
    (t==s.lo || t==s.hi) && return 0//BigInt(1)
    if s.kind==:positive
        A,C=t,1-t
        alpha=1-2A
        base=A-C*(1-s.b)
        return max(0,min(base-A*s.a,alpha*base/(A+alpha)))
    elseif s.kind==:negative
        c,A=t,1-t
        alpha,gamma=1-2A/(n-1),1-2c
        return max(0,min(c*s.b-A*s.a,alpha*c*s.b/(alpha+A),
                         gamma*(c/A-s.a),c/(1+A/alpha+c/gamma)))
    else
        A,C=t,1-t
        p,k=n-s.k,s.k
        alpha,beta,m=1-2A/p,1-2C/k,min(A,C)
        return max(0,min(alpha/2,beta/2,m-A*s.a-C*(1-s.b),
                         (m-C*(1-s.b))/(1+A/alpha),
                         (m-A*s.a)/(1+C/beta),m/(1+A/alpha+C/beta)))
    end
end

function certified_witness(record,n,a,b,lower)
    record["dimension"]==n && vector(record["interval"])==[a,b] || error("Witness cell")
    G,T=vector.(record["generators"]),vector.(record["context"])
    x,w=vector(record["target"]),vector(record["affine_coefficients"])
    length(x)==length(w)==n || error("Witness dimension")
    length(G)==record["generator_count"] || error("Witness generator count")
    all(g->length(g)==2n,G) || error("Generator dimension")
    all(t->length(t)==n && all(v->0<=v<=1,t),T) || error("Context invariant")
    all(v->a<=v<=b,x) || error("Target outside cell")
    sum(abs.(w))==1 || error("Observation normalization")
    knots=sort(unique([0//BigInt(1),a,b,1//BigInt(1)]))
    expected=[lift(fill(0//BigInt(1),n)),lift(fill(1//BigInt(1),n))]
    append!(expected,[vcat(fill(knots[i],n),fill(-knots[i+1],n)) for i in 1:length(knots)-1])
    G==expected || error("Cell representative")
    check_slice(G,n)
    check_combination(vcat(G,lift.(T)),lift(x),vector(record["coefficients"]))
    concrete=minimum(vcat([min(0,sum(w))],[dotq(w,t) for t in T]))
    value=dotq(w,x)
    concrete==rational(record["minimum_concrete_union"]) || error("Concrete minimum")
    value==rational(record["target_value"]) || error("Target value")
    concrete-value==rational(record["certified_gap"])==lower || error("Lower witness gap")
    cells_deficit(x,w)==lower || error("Pointwise support disagreement")
end

function certified_cell(record)
    n=record["dimension"]
    n isa Integer && n>=1 || error("Cell dimension")
    endpoints=vector(record["interval"])
    length(endpoints)==2 || error("Cell interval")
    a,b=endpoints
    0<=a<=b<=1 || error("Cell endpoints")
    lower,upper=rational(record["lower"]),rational(record["upper"])
    0<=lower<=upper || error("Inconsistent cell interval")
    method=record["method"]
    if method=="zero"
        (n<=2 || a==b) && lower==upper==0 || error("False zero cell")
        return 0
    end
    n>=3 && a<b || error("Positive-width certificate class")
    if lower>0
        haskey(record,"witness") || error("Missing lower witness")
        certified_witness(record["witness"],n,a,b,lower)
    end
    if method=="interior"
        b-a<=min(a,1-b) || error("Interior hypothesis missing")
        expected=n==3 ? (b-a)*max(a,1-b)/(2-b+a) : (b-a)/2
        lower==upper==expected || error("Interior formula")
        return 0
    end
    method=="scalar" || error("Unrecognized certificate method")
    tolerance=rational(record["tolerance"])
    0<tolerance && upper-lower<=tolerance || error("Cell tolerance")
    specs,constants=certified_specs(n,a,b)
    all(c->c<=upper,constants) || error("Omitted balanced branch")
    leaves=record["leaves"]
    isempty(leaves) && error("No branch cover")
    all(leaf->leaf["branch"] isa Integer && 0<=leaf["branch"]<length(specs),leaves) || error("Unknown branch")
    for (index,spec) in enumerate(specs)
        group=sort([leaf for leaf in leaves if leaf["branch"]==index-1],by=leaf->rational(leaf["lo"]))
        isempty(group) && error("Branch omitted")
        reached=spec.lo
        for leaf in group
            lo,hi=rational(leaf["lo"]),rational(leaf["hi"])
            lo==reached && lo<hi<=spec.hi || error("Gap, overlap or reversed cover")
            value=certified_value(n,spec,(lo+hi)/2)
            value==rational(leaf["midpoint_value"]) || error("Midpoint value")
            value+hi-lo<=upper || error("Uncertified scalar upper bound")
            reached=hi
        end
        reached==spec.hi || error("Incomplete branch domain")
    end
    record["evaluations"]==2length(leaves)-length(specs) || error("Bisection count")
    length(leaves)
end

function certified_partition(record)
    n,M=record["dimension"],record["generator_budget"]
    M isa Integer && M>=3 || error("Generator budget")
    knots=vector(record["partition"])
    length(knots)==M-1 || error("Lower certificate needs exactly M-2 cells")
    first(knots)==0 && last(knots)==1 && all(diff(knots).>0) || error("Partition")
    cells=record["cells"]
    length(cells)==length(knots)-1 || error("Cell count")
    leaves=0
    for (i,cell) in enumerate(cells)
        cell["dimension"]==n && vector(cell["interval"])==knots[i:i+1] || error("Partition cell mismatch")
        leaves+=certified_cell(cell)
    end
    rational(record["lower"])==minimum(rational(c["lower"]) for c in cells) || error("Dual lower bound")
    rational(record["upper"])==maximum(rational(c["upper"]) for c in cells) || error("Primal upper bound")
    leaves
end

function certified_algebra(data)
    data["identities"]==3 && data["evaluations"]==69 || error("M4 algebra inventory")
    # Distinct rational nodes from the Python study; residual degrees <=3.
    for c in (-BigInt(11):BigInt(11)).//BigInt(13)
        (1-2c)-(1-c)^2/2==-(c^2+2c-1)/2 || error("Root crossing identity")
        (2c-6c^2)*(1-c)+2*(c^2-2c^3)==2c*(1-3c+c^2) || error("Derivative identity")
        (1-3c+c^2)-(2-5c)==c^2+2c-1 || error("Derivative sign identity")
    end
    lo,hi=vector(data["root_interval"])
    2//BigInt(5)<lo<hi<1//BigInt(2) || error("Root isolation domain")
    lo^2+2lo-1<0<hi^2+2hi-1 || error("Root isolation signs")
    vector(data["loss_interval"])==[lo^2/2,hi^2/2] || error("Irrational loss enclosure")
    1//BigInt(24)<lo^2/2 || error("Positive-singleton comparison")
end

function certified_negative_checks(cell,partition,algebra)
    mutations=[
        c->(c["upper"]="0"),
        c->(c["lower"]="1"),
        c->(c["method"]="zero"),
        c->(c["method"]="interior"),
        c->popfirst!(c["leaves"]),
        c->(c["leaves"][1]["branch"]=999),
        c->(c["leaves"][1]["midpoint_value"]="99"),
        c->(c["leaves"][1]["hi"]=c["leaves"][1]["lo"]),
        c->delete!(c,"witness"),
        c->(c["witness"]["context"][1][1]="2"),
        c->(c["witness"]["coefficients"][1]="1/2"),
        c->(c["witness"]["affine_coefficients"][1]="2"),
    ]
    function rejects(verify,original,change)
        bad=deepcopy(original)
        change(bad)
        try
            verify(bad)
        catch
            return true
        end
        error("Accepted malformed global certificate")
    end
    for change in mutations
        rejects(certified_cell,cell,change)
    end
    partition_mutations=[c->(c["generator_budget"]+=1),c->(c["lower"]="1"),
                         c->(c["upper"]="0"),c->pop!(c["cells"]),
                         c->(c["partition"][end]="3/4")]
    for change in partition_mutations
        rejects(certified_partition,partition,change)
    end
    rejects(certified_algebra,algebra,c->(c["root_interval"][1]="0"))
    rejects(certified_algebra,algebra,c->(c["loss_interval"][1]="0"))
    length(mutations)+length(partition_mutations)+2
end

function verify_affine_certified(input,output; status_io=stdout)
    abspath(input)!=abspath(output) || error("Receipt would overwrite input")
    receipt=Dict{String,Any}("complete"=>false,"time_utc"=>string(now(UTC)),
        "julia"=>string(VERSION),"arithmetic"=>"Rational{BigInt}",
        "source_sha256"=>Dict("tests/"*name=>bytes2hex(sha256(read(joinpath(@__DIR__,name))))
                             for name in ["verify_affine_certified.jl","verify_affine_cells.jl","verify_witnesses.jl"]),
        "scope"=>"Global scalar covers and dual partition intervals conditional on D8/D9; not mechanized theorem verification")
    mkpath(dirname(output))
    open(io->TOML.print(io,receipt),output,"w")
    try
        receipt["input_sha256"]=bytes2hex(sha256(read(input)))
        data=TOML.parsefile(input)
        get(data,"complete",false) && data["schema"]=="affine-certified-minimax-v1" || error("Incomplete/wrong input")
        isempty(data["controls"]) && error("No cell controls")
        isempty(data["partitions"]) && error("No finite-budget certificates")
        leaves=sum(certified_cell.(data["controls"]))
        leaves+=sum(certified_partition.(data["partitions"]))
        certified_algebra(data["algebra"])
        receipt["cell_controls_verified"]=length(data["controls"])
        receipt["partitions_verified"]=length(data["partitions"])
        receipt["scalar_intervals_verified"]=leaves
        receipt["algebra_evaluations"]=69
        receipt["malformed_inputs_rejected"]=certified_negative_checks(first(data["controls"]),first(data["partitions"]),data["algebra"])
        println(status_io,"VERIFIED global certificates: ",receipt["partitions_verified"],
                " partitions / ",receipt["scalar_intervals_verified"]," scalar intervals")
        receipt["complete"]=true
    catch exc
        receipt["error"]=sprint(showerror,exc)
        rethrow()
    finally
        open(io->TOML.print(io,receipt),output,"w")
    end
    receipt
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    directory=normpath(joinpath(@__DIR__,"..","results","theory"))
    verify_affine_certified(joinpath(directory,"affine_certified.toml"),joinpath(directory,"affine_certified_julia.toml"))
end
