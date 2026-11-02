# Isolated reference benchmark; this is not part of the CRN-craft engine.
const LOAD_START = time_ns()
using Catalyst, JumpProcesses, JSON3, Random, Statistics, LinearAlgebra
using Catalyst: SciMLBase
BLAS.set_num_threads(1)
const LOAD_SECONDS = (time_ns() - LOAD_START) / 1e9
const MODELS = Dict{String,Any}()
const MODEL_INPUTS = Dict{String,Any}()
const PROBLEMS = Dict{Tuple{String,String},Any}()
const ORDER = Dict{String,Vector{Int}}()
const BUILD_SECONDS = Dict{String,Float64}()
const ALGORITHMS = Dict(
    "Direct" => Direct(), "SortingDirect" => SortingDirect(),
    "DirectCR" => DirectCR(), "NRM" => NRM(), "RSSACR" => RSSACR(),
)

function load_models(path)
    include(path) # trusted generated DSL, emitted from the four frozen fixtures
    for (name, spec) in MODEL_INPUTS
        rn = MODELS[name]
        actual = [first(split(string(s), "(")) for s in species(rn)]
        inds = [findfirst(==(s), actual) for s in spec.names]
        any(isnothing, inds) && error("species mapping failed: " * name)
        ORDER[name] = Int[i for i in inds]
        for (label, algorithm) in ALGORITHMS
            t = time_ns()
            u0 = [Symbol(s) => Int(x) for (s, x) in zip(spec.names, spec.initial)]
            PROBLEMS[(name, label)] = JumpProblem(rn, u0, (0.0, 100.0);
                aggregator=algorithm, save_positions=(false, false))
            BUILD_SECONDS[name * "/" * label] = (time_ns() - t) / 1e9
        end
    end
end

function batch(prob, indices, n::Int, seed::Int, keep::Bool)
    rng = Xoshiro(seed)
    ns = length(indices)
    means = zeros(ns)
    m2 = zeros(ns)
    final = keep ? zeros(Int, n, ns) : zeros(Int, 0, 0)
    for i in 1:n
        sol = solve(prob, SSAStepper(); rng=rng, saveat=[0.0, 100.0])
        SciMLBase.successful_retcode(sol) || error("incomplete Julia trajectory")
        sol.t == [0.0, 100.0] || error("unexpected saved times")
        u = sol.u[end]
        all(>=(0), u) || error("negative population")
        for j in 1:ns
            x = u[indices[j]]
            x isa Integer || error("noninteger population")
            delta = Float64(x) - means[j]
            means[j] += delta / i
            m2[j] += delta * (Float64(x) - means[j])
            keep && (final[i, j] = x)
        end
    end
    return Dict{String,Any}("mean" => means, "variance" => m2 / (n - 1),
        "final_counts" => keep ? [collect(final[i, :]) for i in 1:n] : nothing)
end

function probe(name, states)
    rn = MODELS[name]
    sp = species(rn)
    inds = ORDER[name]
    laws = [jumpratelaw(r) for r in reactions(rn)]
    rates = []
    for values in states
        replacements = Dict(sp[inds[j]] => Int(x) for (j, x) in enumerate(values))
        push!(rates, [Float64(Catalyst.Symbolics.value(
            Catalyst.Symbolics.substitute(law, replacements))) for law in laws])
    end
    mat = Matrix(netstoichmat(rn))[inds, :]
    changes = [collect(mat[:, j]) for j in axes(mat, 2)]
    return Dict("rates" => rates, "changes" => changes,
        "reactions" => length(reactions(rn)), "species" => length(sp))
end

function handle(req)
    cmd = String(req.command)
    if cmd == "load"
        t = time_ns()
        load_models(String(req.path))
        return Dict("load_seconds" => LOAD_SECONDS,
            "model_setup_seconds" => (time_ns() - t) / 1e9,
            "problem_setup_seconds" => BUILD_SECONDS,
            "versions" => Dict("Julia" => string(VERSION),
                "Catalyst" => string(pkgversion(Catalyst)),
                "JumpProcesses" => string(pkgversion(JumpProcesses)),
                "JSON3" => string(pkgversion(JSON3))),
            "threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads())
    elseif cmd == "probe"
        return probe(String(req.name), req.states)
    elseif cmd == "run"
        key = (String(req.name), String(req.algorithm))
        # SortingDirect retains its learned order. Reset the pristine batch state
        # outside the timer; repeated batches then use identical stochastic work.
        prob = deepcopy(PROBLEMS[key])
        n, seed, keep = Int(req.paths), Int(req.seed), Bool(req.keep)
        t = time_ns()
        result = batch(prob, ORDER[key[1]], n, seed, keep)
        result["seconds"] = (time_ns() - t) / 1e9
        return result
    end
    error("unknown command")
end

for line in eachline(stdin)
    try
        result = handle(JSON3.read(line))
        println(JSON3.write(Dict("ok" => true, "result" => result)))
    catch e
        showerror(stderr, e, catch_backtrace())
        println(stderr)
        println(JSON3.write(Dict("ok" => false, "error" => sprint(showerror, e))))
    end
    flush(stdout)
    flush(stderr)
end
