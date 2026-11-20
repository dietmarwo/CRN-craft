// Study-specific reference samplers. No CRN-craft engine or ABI changes.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <array>
#include <vector>

namespace {
constexpr double log_two = 0.69314718055994530942;
constexpr double atp_decay = 0.2;
constexpr double atp_supply = 0.1;
constexpr double atp_maintenance = 0.01;
constexpr double mutant_synthesis = 1.5;
constexpr int count_limit = 1000000;

enum Mechanism { drift, size_advantage, transcription };
enum Status { horizon, collapsed, empty, count_overflow };

struct Config {
    int mechanism, population, initial_mutants, hours;
    int wild_delay, mutant_delay, wild_only;
    double half_life, mutation, initiation, step;
};

struct Endpoint {
    int wild, mutant, clones, largest, status;
    double time;
    uint64_t steps;
};

uint64_t mix(uint64_t value) {
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

class Random {
public:
    using result_type = uint64_t;
    explicit Random(uint64_t seed) : state(mix(seed)) {}
    static constexpr uint64_t min() { return 0; }
    static constexpr uint64_t max() { return UINT64_MAX; }
    uint64_t operator()() { return mix(state += UINT64_C(0x9e3779b97f4a7c15)); }
    double unit() { return ((*this)() >> 11) * 0x1.0p-53; }
    int rank(int n) {
        return std::uniform_int_distribution<int>(0, n - 1)(*this);
    }

private:
    uint64_t state;
};

struct Chance {
    double p, log_survival;

    static Chance bernoulli(double p) { return {p, std::log1p(-p)}; }
    static Chance hazard(double rate) { return {-std::expm1(-rate), -rate}; }

    int draw(Random& rng, int n, const double* zero = nullptr) const {
        if (n == 0 || p <= 0) return 0;
        if (p >= 1) return n;
        if (n == 1) return rng.unit() < p;
        if (p * n > 24 || p > 0.08) {
            return std::binomial_distribution<int>(n, p)(rng);
        }
        // Inversion of the binomial PMF. The large-mean case above avoids
        // underflow and a linear walk through a large number of successes.
        double probability = zero ? zero[n] : std::exp(n * log_survival);
        const double odds = p / (1 - p);
        double u = rng.unit();
        int k = 0;
        while (k < n && u >= probability) {
            u -= probability;
            ++k;
            probability *= odds * (n - k + 1) / k;
        }
        return k;
    }
};

// Active counts are compact. A temporary departure to a busy cohort can
// remove a label from the free pool, and completion can add that label back.
class Labels {
public:
    int new_clone(int count = 1) {
        const int id = static_cast<int>(entries.size());
        entries.push_back({count, static_cast<int>(active.size())});
        if (count) active.push_back(id);
        total += count;
        squares += static_cast<int64_t>(count) * count;
        return id;
    }

    void add(int id, int increment) {
        auto& e = entries[id];
        const int count = e.count;
        total += increment;
        e.count += increment;
        squares += static_cast<int64_t>(increment) * (count + e.count);
        if (count == 0 && e.count > 0) {
            e.index = static_cast<int>(active.size());
            active.push_back(id);
        } else if (count > 0 && e.count == 0) {
            entries[active.back()].index = e.index;
            active[e.index] = active.back();
            active.pop_back();
        }
    }

    int pick(int rank, int exclude = -1) const {
        for (int id : active) {
            if (id != exclude) {
                rank -= entries[id].count;
                if (rank < 0) return id;
            }
        }
        throw std::logic_error("label rank exceeds held copies");
    }

    int take(Random& rng) {
        const int id = pick(rng.rank(total));
        add(id, -1);
        return id;
    }

    int contested_victim(int64_t rank, int population) const {
        for (int id : active) {
            const int64_t n = entries[id].count;
            rank -= n * (population - n);
            if (rank < 0) return id;
        }
        throw std::logic_error("invalid contested replacement rank");
    }

    int count(int id) const { return entries[id].count; }

    Endpoint endpoint(int wild, int mutant, Status status, double time,
                      uint64_t steps, int first_mutant_id = 0) const {
        int kinds = 0;
        int largest = 0;
        int counted_mutants = 0;
        for (int id : active) {
            if (id >= first_mutant_id) {
                ++kinds;
                largest = std::max(largest, entries[id].count);
                counted_mutants += entries[id].count;
            }
        }
        if (counted_mutants != mutant) throw std::logic_error("lineage sum differs");
        return {wild, mutant, kinds, largest, status, time, steps};
    }

    int total = 0;
    int64_t squares = 0;

private:
    struct Entry { int count = 0; int index = 0; };
    std::vector<Entry> entries;
    std::vector<int> active;
};

Endpoint drift_aggregate(const Config& cfg, Random& rng) {
    const int n = cfg.population;
    const double delta = log_two / cfg.half_life;
    const double copies = delta * (1 - cfg.mutation) / (n - 1);
    const double mutation = delta * cfg.mutation;
    int m = cfg.initial_mutants;
    double time = 0;
    uint64_t steps = 0;
    while (m < n) {
        const int w = n - m;
        const double down = copies * m * w;
        const double sum = 2 * down + mutation * w;
        if (sum == 0) break;
        time -= std::log1p(-rng.unit()) / sum;
        if (time > cfg.hours) break;
        m += rng.unit() * sum < down ? -1 : 1;
        ++steps;
    }
    return {n - m, m, -1, -1, horizon, double(cfg.hours), steps};
}

Endpoint drift_labelled(const Config& cfg, Random& rng, Random& labels_rng) {
    const int n = cfg.population;
    const double delta = log_two / cfg.half_life;
    const double copies = delta * (1 - cfg.mutation) / (n - 1);
    const double mutation = delta * n * cfg.mutation;
    Labels labels;
    const int wild = labels.new_clone(n - cfg.initial_mutants);
    if (cfg.initial_mutants) labels.new_clone(cfg.initial_mutants);
    const int64_t pairs = int64_t(n) * n;
    double time = 0;
    uint64_t steps = 0;
    while (true) {
        // Self-replacements leave no trace. Each ordered cross-label pair
        // of genomes has the same hazard because the victim cannot parent.
        const auto cross_pairs = pairs - labels.squares;
        const double sum = mutation + copies * cross_pairs;
        if (sum == 0) break;
        time -= std::log1p(-rng.unit()) / sum;
        if (time > cfg.hours) break;
        if (rng.unit() * sum < mutation) {
            labels.take(labels_rng);
            labels.new_clone();
        } else {
            const auto rank = std::uniform_int_distribution<int64_t>(
                0, cross_pairs - 1)(labels_rng);
            const int victim = labels.contested_victim(rank, n);
            const int parent = labels.pick(labels_rng.rank(n - labels.count(victim)), victim);
            labels.add(victim, -1);
            labels.add(parent, 1);
        }
        ++steps;
    }
    const int w = labels.count(wild);
    return labels.endpoint(w, n - w, horizon, cfg.hours, steps, 1);
}

struct Context {
    explicit Context(const Config& cfg)
        : death(cfg.mechanism == size_advantage
                    ? Chance::bernoulli(log_two / cfg.half_life)
                    : Chance::hazard(log_two * cfg.step / cfg.half_life)),
          mutation(Chance::bernoulli(cfg.mutation)),
          mutant_mutation(Chance::bernoulli(cfg.wild_only ? 0 : cfg.mutation)) {
        survival[0] = 1;
        const double survival_one = 1 - death.p;
        for (int n = 1; n <= 4096; ++n) survival[n] = survival[n - 1] * survival_one;
    }

    int deaths(Random& rng, int n) const {
        return death.draw(rng, n, n <= 4096 ? survival : nullptr);
    }

    Chance death, mutation, mutant_mutation;
    double survival[4097];
};

double size_capacity(int n, double delta, double wild_delay, double initiation) {
    return n / (initiation * (2 * std::exp(-delta * wild_delay) - 1) / delta - 1);
}

struct Busy {
    int wild = 0;
    int mutant = 0;
    std::vector<int> labels;
};

template <bool labelled>
Endpoint size_cell(const Config& cfg, const Context& ctx, Random& rng, Random& label_rng) {
    const double capacity = size_capacity(cfg.population, log_two / cfg.half_life,
                                          cfg.wild_delay, cfg.initiation);
    const int buckets = std::max(cfg.wild_delay, cfg.mutant_delay) + 1;
    std::vector<Busy> queue(buckets);
    Labels free;
    int w = cfg.population - cfg.initial_mutants;
    int m = cfg.initial_mutants;
    if (labelled && m) free.new_clone(m);
    int index = 0;
    int hour = 0;
    int busy_w = 0, busy_m = 0;
    Status status = horizon;
    for (; hour < cfg.hours; ++hour) {
        auto& due = queue[index];
        const int mutated_w = ctx.mutation.draw(rng, due.wild);
        const int mutated_m = ctx.mutant_mutation.draw(rng, due.mutant);
        w += 2 * due.wild - mutated_w;
        m += 2 * due.mutant + mutated_w;
        if constexpr (labelled) {
            // Conditional on their number, mutant-template mutations have
            // a uniform subset of due parents. Return all original templates.
            for (int j = 0; j < mutated_m; ++j) {
                std::swap(due.labels[j], due.labels[j + label_rng.rank(due.mutant - j)]);
            }
            for (int j = 0; j < due.mutant; ++j) {
                free.add(due.labels[j], j < mutated_m ? 1 : 2);
            }
            due.labels.clear();
            for (int j = 0; j < mutated_w + mutated_m; ++j) free.new_clone();
        }
        due.wild = due.mutant = 0;
        w -= ctx.deaths(rng, w);
        const int free_deaths = ctx.deaths(rng, m);
        m -= free_deaths;
        if constexpr (labelled) {
            for (int j = 0; j < free_deaths; ++j) free.take(label_rng);
        }
        busy_w = busy_m = 0;
        for (auto& future : queue) {
            future.wild -= ctx.deaths(rng, future.wild);
            const int lost = ctx.deaths(rng, future.mutant);
            future.mutant -= lost;
            if constexpr (labelled) {
                for (int j = 0; j < lost; ++j) {
                    const int dead = label_rng.rank(static_cast<int>(future.labels.size()));
                    future.labels[dead] = future.labels.back();
                    future.labels.pop_back();
                }
            }
            busy_w += future.wild;
            busy_m += future.mutant;
        }
        const int total = w + m + busy_w + busy_m;
        if (!total || total > count_limit) {
            status = total == 0 ? empty : count_overflow;
            ++hour;
            break;
        }
        const auto start = Chance::bernoulli(cfg.initiation * capacity / (capacity + total));
        const int new_w = start.draw(rng, w);
        const int new_m = start.draw(rng, m);
        auto& dw = queue[(index + cfg.wild_delay) % buckets];
        auto& dm = queue[(index + cfg.mutant_delay) % buckets];
        dw.wild += new_w;
        dm.mutant += new_m;
        if constexpr (labelled) {
            for (int j = 0; j < new_m; ++j) dm.labels.push_back(free.take(label_rng));
        }
        w -= new_w;
        m -= new_m;
        busy_w += new_w;
        busy_m += new_m;
        if (++index == buckets) index = 0;
    }
    if constexpr (labelled) {
        for (const auto& future : queue) {
            for (int label : future.labels) free.add(label, 1);
        }
        return free.endpoint(w + busy_w, m + busy_m, status, hour, hour);
    }
    return {w + busy_w, m + busy_m, -1, -1, status, double(hour), uint64_t(hour)};
}

template <bool labelled>
Endpoint atp_cell(const Config& cfg, const Context& ctx, Random& rng, Random& label_rng) {
    const double delta = log_two / cfg.half_life;
    const double equilibrium = cfg.population * (atp_supply - atp_maintenance) / atp_decay;
    const double capacity = equilibrium * delta / (1 - delta);
    const double adjustment = -std::expm1(-atp_decay * cfg.step);
    const auto visits = static_cast<uint64_t>(cfg.hours / cfg.step);
    Labels labels;
    std::vector<int> faithful_parents;
    double atp = equilibrium;
    int w = cfg.population - cfg.initial_mutants;
    int m = cfg.initial_mutants;
    if (labelled && m) labels.new_clone(m);
    Status status = horizon;
    double time = cfg.hours;
    uint64_t idx = 0;
    for (; idx < visits; ++idx) {
        w -= ctx.deaths(rng, w);
        const int loss = ctx.deaths(rng, m);
        m -= loss;
        const double rate = capacity * cfg.step / (capacity + atp);
        const int births_w = Chance::hazard(rate).draw(rng, w);
        const int births_m = m ? Chance::hazard(mutant_synthesis * rate).draw(rng, m) : 0;
        const int mutated_w = ctx.mutation.draw(rng, births_w);
        const int mutated_m = ctx.mutant_mutation.draw(rng, births_m);
        w += births_w - mutated_w;
        m += births_m + mutated_w;
        if constexpr (labelled) {
            for (int j = 0; j < loss; ++j) labels.take(label_rng);
            // Only survivors present at the start of the birth draw may
            // parent a daughter, and each supplies at most one.
            faithful_parents.clear();
            for (int j = mutated_m; j < births_m; ++j) {
                faithful_parents.push_back(labels.take(label_rng));
            }
            for (int label : faithful_parents) labels.add(label, 2);
            for (int j = 0; j < mutated_w + mutated_m; ++j) labels.new_clone();
        }
        const double target = (atp_supply * w - atp_maintenance * (w + m)) / atp_decay;
        const double next = atp + (target - atp) * adjustment;
        if (w + m == 0 || w + m > count_limit || next <= 0) {
            if (!w && !m) {
                status = empty;
                time = (idx + 1) * cfg.step;
            } else if (w + m > count_limit) {
                status = count_overflow;
                time = (idx + 1) * cfg.step;
            } else {
                status = collapsed;
                time = idx * cfg.step + std::log1p(atp / -target) / atp_decay;
            }
            ++idx;
            break;
        }
        atp = next;
    }
    if constexpr (labelled) return labels.endpoint(w, m, status, time, idx);
    return {w, m, -1, -1, status, time, idx};
}

} // namespace

extern "C" int kk14_batch(const Config* cfg, int labelled, int count,
                           const uint64_t* seeds, Endpoint* out) {
    try {
        const Context ctx(*cfg);
        for (int i = 0; i < count; ++i) {
            Random demographic(seeds[i]);
            Random ancestry(seeds[i] ^ UINT64_C(0xd1b54a32d192ed03));
            if (cfg->mechanism == drift) {
                out[i] = labelled ? drift_labelled(*cfg, demographic, ancestry)
                                  : drift_aggregate(*cfg, demographic);
            } else if (cfg->mechanism == size_advantage) {
                out[i] = labelled ? size_cell<true>(*cfg, ctx, demographic, ancestry)
                                  : size_cell<false>(*cfg, ctx, demographic, ancestry);
            } else {
                out[i] = labelled ? atp_cell<true>(*cfg, ctx, demographic, ancestry)
                                  : atp_cell<false>(*cfg, ctx, demographic, ancestry);
            }
        }
    } catch (const std::exception&) {
        return 1;
    }
    return 0;
}

namespace {
using Pools = std::array<double, 4>;
using Flux = std::array<double, 2>;

struct DelayPoint {
    double hour, wild_free, wild_busy, mutant_free, mutant_busy;
};

class Delay {
public:
    Delay(double delay_w, double delay_m, double half_life, double repl,
          double step)
        : delta(log_two / half_life), initiation(repl), h(step),
          capacity(size_capacity(1000, delta, delay_w, repl)),
          lag_w(std::lround(delay_w / h)), lag_m(std::lround(delay_m / h)),
          survives_w(std::exp(-delta * delay_w)), survives_m(std::exp(-delta * delay_m)),
          history(std::max(lag_w, lag_m) + 1) {}

    Flux starts(const Pools& y) const {
        const double total = y[0] + y[1] + y[2] + y[3];
        const double rate = initiation * capacity / (capacity + total);
        return {rate * y[0], rate * y[2]};
    }

    Pools derivative(const Pools& y, Flux start, int index, bool left) const {
        const double w = survives_w * old(index - lag_w, 0, left);
        const double m = survives_m * old(index - lag_m, 1, left);
        return {2*w - start[0] - delta*y[0], start[0] - w - delta*y[1],
                2*m - start[1] - delta*y[2], start[1] - m - delta*y[3]};
    }

    Flux& at(int i) { return history[i % history.size()]; }

private:
    double old(int i, int g, bool left) const {
        // There were no initiations before t=0. At the first possible
        // completion the flux jumps; the Heun predictor uses its left limit.
        if (i < 0 || (left && i == 0)) return 0;
        return history[i % history.size()][g];
    }

    double delta, initiation, h, capacity;
    int lag_w, lag_m;
    double survives_w, survives_m;
    std::vector<Flux> history;
};
} // namespace

extern "C" int kk14_delay(double wild_delay, double mutant_delay, double half_life,
                           double initiation, double step, int n_steps,
                           int stride, DelayPoint* result) {
    try {
        Delay dde(wild_delay, mutant_delay, half_life, initiation, step);
        Pools y{999, 0, 1, 0};
        dde.at(0) = dde.starts(y);
        int count = 0;
        const auto write = [&](double hour) {
            result[count++] = {hour, y[0], y[1], y[2], y[3]};
        };
        write(0);
        for (int i = 0; i < n_steps; ++i) {
            const Pools prev = y;
            const auto d0 = dde.derivative(y, dde.at(i), i, false);
            for (int j = 0; j < 4; ++j) y[j] += step * d0[j];
            const auto d1 = dde.derivative(y, dde.starts(y), i + 1, true);
            for (int j = 0; j < 4; ++j) y[j] = prev[j] + 0.5 * step * (d0[j] + d1[j]);
            dde.at(i + 1) = dde.starts(y);
            const double wild = y[0] + y[1];
            if (wild < 1) {
                const double before = prev[0] + prev[1];
                const double portion = (before - 1) / (before - wild);
                for (int j = 0; j < 4; ++j) y[j] = prev[j] + portion * (y[j] - prev[j]);
                write(step * (i + portion));
                return count;
            }
            if ((i + 1) % stride == 0 || i + 1 == n_steps) write(step * (i + 1));
        }
        return count;
    } catch (const std::exception&) {
        return -1;
    }
}
