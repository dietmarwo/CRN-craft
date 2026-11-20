// Equation 7 on a fixed population. All kernels are local to this example.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

namespace {
constexpr double ln2 = 0.69314718055994530942;
constexpr int n_bins = 5;
constexpr int n_summary = 7;

struct Parameters {
    int copies, cap;
    double mutation, advantage, advantageous, half_days;
};

double clock(const Parameters& p) { return ln2 / p.half_days * p.copies; }

std::array<double, 3> daughters(const Parameters& p, int wild, int neutral, int advantaged) {
    const double biased = (1.0 + p.advantage) * advantaged;
    const double total = wild + neutral + biased;
    return {wild * (1 - p.mutation) / total,
            (neutral + wild * p.mutation * (1 - p.advantageous)) / total,
            (biased + wild * p.mutation * p.advantageous) / total};
}

int bin(int mutants, int total) {
    return (mutants > 0) + (20 * mutants > total) + (5 * mutants > total)
         + (2 * mutants > total);
}

// Three positive directions are (neutral++, adv unchanged), (neutral
// unchanged, adv++), and (neutral++, adv--).
constexpr std::array<int, 6> du{1,-1,0,0,1,-1};
constexpr std::array<int, 6> da{0,0,1,-1,-1,1};

std::array<double, 6> outgoing(const Parameters& p, int u, int a) {
    const int w = p.copies - u - a;
    const auto b = daughters(p, w, u, a);
    const double iw = static_cast<double>(w) / p.copies;
    const double iu = static_cast<double>(u) / p.copies;
    const double ia = static_cast<double>(a) / p.copies;
    return {b[1]*iw, b[0]*iu, b[2]*iw, b[0]*ia, b[1]*ia, b[2]*iu};
}

class Lattice {
public:
    explicit Lattice(const Parameters& p)
        : copies(p.copies), cap(p.cap), stride(copies + 3),
          space((cap + 3) * stride),
          shifts{stride,-stride,1,-1,stride-1,1-stride},
          stay(space, 0) {
        for (auto& arr : rates) arr.resize(space);
        for (int u = 0; u <= cap; ++u) {
            for (int a = 0; a <= copies - u; ++a) {
                const int here = index(u,a);
                const auto q = outgoing(p, u, a);
                stay[here] = 1;
                for (int k = 0; k < 6; ++k) {
                    rates[k][here] = q[k];
                    stay[here] -= q[k];
                }
            }
        }
    }

    int index(int u, int a) const { return (u + 1) * stride + a + 1; }

    std::vector<double> initial() const {
        std::vector<double> state(space);
        state[index(0,0)] = 1;
        return state;
    }

    void summarize(const std::vector<double>& state, double* output) const {
        std::fill(output, output+n_summary, 0);
        for (int u = 0; u <= cap; ++u) {
            for (int a = 0; a <= copies - u; ++a) {
                const double mass = state[index(u,a)];
                output[bin(a+u,copies)] += mass;
                output[6] -= std::min(0.0, mass);
            }
        }
        double total = 0;
        for (int k = 0; k < n_bins; ++k) total += output[k];
        output[5] = 1.0 - total;
    }

    void discrete(const uint64_t* times, int reports, double* out, double* final) const {
        // Transpose the stencil once so arrival coefficients are contiguous.
        std::array<std::vector<double>, 6> incoming;
        for (int k = 0; k < 6; ++k) {
            incoming[k].resize(space);
            for (int u = 0; u <= cap; ++u) {
                for (int a = 0; a <= copies-u; ++a) {
                    const int here = index(u,a);
                    incoming[k][here] = rates[k][here-shifts[k]];
                }
            }
        }
        auto p = initial();
        std::vector<double> next(space);
        uint64_t step = 0;
        for (int r = 0; r < reports; ++r) {
            for (; step < times[r]; ++step) {
                const double* from = p.data();
                double* into = next.data();
                for (int u = 0; u <= cap; ++u) {
                    const int start = index(u,0);
                    const int end = start + copies - u + 1;
                    for (int i = start; i < end; ++i) {
                        const double value =
                            stay[i]*from[i]
                            + incoming[0][i]*from[i-stride] + incoming[1][i]*from[i+stride]
                            + incoming[2][i]*from[i-1] + incoming[3][i]*from[i+1]
                            + incoming[4][i]*from[i-stride+1] + incoming[5][i]*from[i+stride-1];
                        // At most one discarded mass below 1e-280 per state and step.
                        into[i] = value >= 1e-280 ? value : 0;
                    }
                }
                p.swap(next);
            }
            summarize(p, out + n_summary*r);
        }
        if (final) {
            int k = 0;
            for (int u = 0; u <= cap; ++u) {
                for (int a = 0; a <= copies-u; ++a) final[k++] = p[index(u,a)];
            }
        }
    }

    int copies, cap, stride, space;
    std::array<int,6> shifts;
    std::array<std::vector<double>,6> rates;
    std::vector<double> stay;
};

struct Factor {
    int id;
    double lower, rhs_middle, rhs_upper, pivot_inverse, back;
};

class Direction {
public:
    Direction(const Lattice& l, int axis, double cw) {
        if (axis == 0) {
            for (int a = 0; a <= l.copies; ++a) {
                line(l, axis, cw, l.index(0,a), std::min(l.cap,l.copies-a)+1);
            }
        } else if (axis == 1) {
            for (int u = 0; u <= l.cap; ++u) {
                line(l, axis, cw, l.index(u,0), l.copies-u+1);
            }
        } else {
            for (int total = 0; total <= l.copies; ++total) {
                line(l, axis, cw, l.index(0,total), std::min(l.cap,total)+1);
            }
        }
        old.resize(factors.size()+2);
        rhs.resize(factors.size()+2);
    }

    void apply(std::vector<double>& state) {
        const int n = static_cast<int>(factors.size());
        // Padding and zero boundary coefficients separate the line solves.
        for (int i = 0; i < n; ++i) old[i+1] = state[factors[i].id];
        for (int i = 0; i < n; ++i) {
            const auto& f = factors[i];
            rhs[i+1] = (f.rhs_middle*old[i+1] + f.rhs_upper*old[i+2]
                        + f.lower*(old[i]+rhs[i])) * f.pivot_inverse;
        }
        double result = 0;
        for (int i = n-1; i >= 0; --i) {
            result = rhs[i+1] + factors[i].back*result;
            // Avoid large costs from subnormal tails without rescaling mass.
            if (std::abs(result) < 1e-280) result = 0;
            state[factors[i].id] = result;
        }
    }

private:
    void line(const Lattice& l, int axis, double cw, int start, int size) {
        const int shift = l.shifts[2*axis];
        const auto& plus = l.rates[2*axis];
        const auto& minus = l.rates[2*axis+1];
        double previous = 0;
        for (int i = 0; i < size; ++i) {
            const int id = start + i*shift;
            const double lower = i > 0 ? cw*plus[id-shift] : 0;
            const double upper = i+1 < size ? cw*minus[id+shift] : 0;
            const double exit = cw*(plus[id]+minus[id]);
            const double inv = 1 / (1 + exit - lower*previous);
            const double back = upper*inv;
            factors.push_back({id,lower,1-exit,upper,inv,back});
            previous = back;
        }
    }
    std::vector<Factor> factors;
    std::vector<double> old, rhs;
};

void split(const Parameters& p, const Lattice& l, const double* days,
           int reports, double max_substep, double* out) {
    auto state = l.initial();
    double time = 0, last_dt = 0;
    std::vector<Direction> direction;
    for (int i = 0; i < reports; ++i) {
        const double interval = days[i]-time;
        if (interval > 0) {
            const int steps = static_cast<int>(std::ceil(interval/max_substep));
            const double dt = interval/steps;
            if (dt != last_dt) {
                direction.clear();
                direction.emplace_back(l, 0, clock(p)*dt/4);
                direction.emplace_back(l, 1, clock(p)*dt/4);
                direction.emplace_back(l, 2, clock(p)*dt/2);
                last_dt = dt;
            }
            for (int k = 0; k < steps; ++k) {
                direction[0].apply(state);
                direction[1].apply(state);
                direction[2].apply(state);
                direction[1].apply(state);
                direction[0].apply(state);
            }
        }
        l.summarize(state, out + n_summary*i);
        time = days[i];
    }
}

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

    int binomial(int n, double p) {
        if (!n || p <= 0) return 0;
        if (p >= 1) return n;
        if (n == 1) return unit() < p;
        if (p*n > 24 || p > 0.08) return std::binomial_distribution<int>(n,p)(*this);
        const double odds = p/(1-p);
        double probability = std::exp(n*std::log1p(-p));
        double u = unit();
        int k = 0;
        while (k < n && u >= probability) {
            u -= probability;
            ++k;
            probability *= odds*(n-k+1)/k;
        }
        return k;
    }

private:
    uint64_t state;
};

void sample(const Parameters& p, const uint64_t* times, int reports,
            Random& rng, int* counts) {
    uint64_t step = 0;
    int i = 0, u = 0, a = 0;
    while (i < reports) {
        const auto rates = outgoing(p,u,a);
        double change = 0;
        for (double rate : rates) change += rate;
        if (change <= 0) {
            for (; i < reports; ++i) ++counts[n_bins*i+bin(u+a,p.copies)];
            return;
        }
        // The first changed step has a geometric waiting time. Its endpoint
        // belongs to that step, so a coincident observation reads the new state.
        const double wait = 1 + std::floor(std::log1p(-rng.unit())/std::log1p(-change));
        while (i < reports && (wait > times[i]-step)) {
            ++counts[n_bins*i+bin(u+a,p.copies)];
            ++i;
        }
        if (i == reports) return;
        step += static_cast<uint64_t>(wait);
        double pick = rng.unit()*change;
        for (int k=0; k<6; ++k) {
            if (pick < rates[k]) {
                u += du[k];
                a += da[k];
                break;
            }
            pick -= rates[k];
        }
    }
}

void sample_atp(const Parameters& p, const double* days, int reports, double h,
                Random& rng, int* counts, int* stops) {
    // Hourly reconstruction of Model 1. Counts at collapse are absorbed.
    constexpr double supply=0.1, loss=0.2, maintenance=0.01;
    const double delta = ln2/(24*p.half_days);
    const double eq = p.copies*(supply-maintenance)/loss;
    const double c = eq*delta/(1-delta);
    const double die = -std::expm1(-delta*h);
    const double relax = -std::expm1(-loss*h);
    double energy = eq;
    int w = p.copies, u = 0, a = 0, elapsed = 0;
    bool collapse = false;
    for (int i = 0; i < reports; ++i) {
        const int tick = static_cast<int>(std::llround(days[i]*24/h));
        for (; elapsed < tick && !collapse; ++elapsed) {
            w -= rng.binomial(w,die);
            u -= rng.binomial(u,die);
            a -= rng.binomial(a,die);
            const double rate = c/(c+energy)*h;
            const double birth = -std::expm1(-rate);
            const int ww = rng.binomial(w,birth);
            const int uu = rng.binomial(u,birth);
            const int aa = rng.binomial(a,-std::expm1(-(1+p.advantage)*rate));
            const int mutations = rng.binomial(ww,p.mutation);
            const int new_a = rng.binomial(mutations,p.advantageous);
            w += ww-mutations;
            u += uu+mutations-new_a;
            a += aa+new_a;
            const double target = (supply*w-maintenance*(w+u+a))/loss;
            energy += (target-energy)*relax;
            collapse = energy <= 0;
        }
        if (collapse) ++stops[i];
        if (w+u+a == 0) throw std::runtime_error("empty ATP cell");
        ++counts[n_bins*i+bin(u+a,w+u+a)];
    }
}
} // namespace

extern "C" int kk26_forward(const Parameters* p, const double* days, const uint64_t* ticks,
                           int reports, double max_substep, double* output, double* final) {
    try {
        const Lattice lattice(*p);
        if (max_substep > 0) split(*p, lattice, days, reports, max_substep, output);
        else lattice.discrete(ticks, reports, output, final);
    } catch (const std::exception&) {
        return 1;
    }
    return 0;
}

extern "C" int kk26_sample(const Parameters* p, const double* days, const uint64_t* ticks,
                          int reports, double atp_step, int cells, uint64_t seed,
                          int* counts, int* stops) {
    try {
        for (int i = 0; i < cells; ++i) {
            Random rng(seed+i);
            if (atp_step > 0) sample_atp(*p,days,reports,atp_step,rng,counts,stops);
            else sample(*p,ticks,reports,rng,counts);
        }
    } catch (const std::exception&) {
        return 1;
    }
    return 0;
}

extern "C" void kk26_passage(const Parameters* p, int target, const uint64_t* times,
                            int reports, double* output) {
    // Loss and takeover are absorbing. Keep their mass as part of the chain.
    std::vector<double> up(target+1), down(target+1), stay(target+1,1);
    for (int a = 1; a < target; ++a) {
        const auto q = outgoing(*p,0,a);
        up[a] = q[2];
        down[a] = q[3];
        stay[a] = 1-up[a]-down[a];
    }
    std::vector<double> state(target+1), next(target+1);
    state[1]=1;
    uint64_t step=0;
    for (int r = 0; r < reports; ++r) {
        for (; step < times[r]; ++step) {
            next[0] = state[0] + down[1]*state[1];
            for (int i = 1; i < target; ++i) {
                next[i] = up[i-1]*state[i-1]+stay[i]*state[i]+down[i+1]*state[i+1];
            }
            next[target] = state[target]+up[target-1]*state[target-1];
            state.swap(next);
        }
        double live=0;
        for (int i=1; i<target; ++i) live+=state[i];
        output[3*r]=state[0];
        output[3*r+1]=state[target];
        output[3*r+2]=live;
    }
}
