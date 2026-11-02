# JumpProcesses.jl (via Catalyst.jl) on the frozen SSA workloads W1-W4.
#
# Usage: julia --project=validation/julia validation/julia/bench.jl MODELS.jl OUT.json
# MODELS.jl is written by validation/compare_julia.py from the CRN-craft
# models; it defines `WORKLOADS`, a vector of (label, reaction network,
# initial counts, species names in CRN-craft order) tuples. Each workload runs 64 trajectories over t in
# [0, 100], saving at 0 and 100, serially, for several aggregators. The first
# solve per aggregator is reported apart (it includes compilation); then
# REPEATS solves are timed. 512 trajectories with Direct are kept for the
# equivalence checks.

using Catalyst, JumpProcesses, JSON3, Random

const HORIZON = 100.0
const PATHS = 64
const CHECK_PATHS = 512
const REPEATS = parse(Int, get(ENV, "REPEATS", "5"))
const AGGREGATORS = [
    ("Direct", Direct()),
    ("SortingDirect", SortingDirect()),
    ("RSSA", RSSA()),
    ("RSSACR", RSSACR()),
    ("DirectCR", DirectCR()),
]

include(ARGS[1])

function jump_problem(rn, u0, aggregator, seed)
    tspan = (0.0, HORIZON)
    rng = Xoshiro(seed)
    if isdefined(Catalyst, :JumpInputs)  # Catalyst 14.2 and later
        inputs = JumpInputs(rn, u0, tspan)
        return JumpProblem(inputs, aggregator; save_positions = (false, false), rng = rng)
    else
        dprob = DiscreteProblem(rn, u0, tspan)
        return JumpProblem(rn, dprob, aggregator; save_positions = (false, false), rng = rng)
    end
end

function ensemble(jprob, trajectories)
    eprob = EnsembleProblem(jprob)
    return solve(eprob, SSAStepper(), EnsembleSerial();
        trajectories = trajectories, saveat = [0.0, HORIZON])
end

report = Dict{String, Any}(
    "julia" => string(VERSION),
    "packages" => Dict(
        "Catalyst" => string(pkgversion(Catalyst)),
        "JumpProcesses" => string(pkgversion(JumpProcesses)),
    ),
    "threads" => Threads.nthreads(),
    "workloads" => Dict{String, Any}(),
)

for (label, rn, u0, names) in WORKLOADS
    entry = Dict{String, Any}()
    # Catalyst orders species by appearance; map them to CRN-craft's order.
    order = [replace(string(s), "(t)" => "") for s in species(rn)]
    index = Dict(n => findfirst(==(n), order) for n in names)
    for (name, aggregator) in AGGREGATORS
        build_first = @elapsed (jprob = jump_problem(rn, u0, aggregator, 4242))
        solve_first = @elapsed ensemble(jprob, PATHS)
        times = Float64[]
        for repeat in 1:REPEATS
            jprob = jump_problem(rn, u0, aggregator, 4242 + repeat)
            push!(times, @elapsed ensemble(jprob, PATHS))
        end
        entry[name] = Dict(
            "first_build_seconds" => build_first,
            "first_solve_seconds" => solve_first,
            "runs" => times,
        )
        println(label, " ", name, " first ", round(build_first + solve_first; digits = 2),
            " s, runs ", round.(times; digits = 4))
        flush(stdout)
    end
    # Equivalence check: final counts of 512 Direct trajectories, in the
    # species order of the CRN-craft model.
    sol = ensemble(jump_problem(rn, u0, Direct(), 991), CHECK_PATHS)
    entry["final_counts"] = Dict(
        n => [Int(s.u[end][index[n]]) for s in sol.u] for n in names)
    report["workloads"][label] = entry
end

open(ARGS[2], "w") do io
    JSON3.pretty(io, report)
end
