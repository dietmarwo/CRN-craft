// Example-local integrator for the four-state promoter and total mRNA count.
#include <algorithm>
#include <cmath>
#include <vector>

namespace {
constexpr int band = 4;
constexpr int stride = 2 * band + 1;

// The implicit CME matrix has a positive diagonal and is strictly
// column diagonally dominant. Its banded LU needs no row interchanges.
void solve(std::vector<double>& a, std::vector<double>& x)
{
    const int n = static_cast<int>(x.size());
    for (int col = 0; col < n; ++col) {
        double* source = &a[col * stride];
        const int stop = std::min(col + band + 1, n);
        for (int row = col + 1; row < stop; ++row) {
            double* dest = &a[row * stride + col - row + band];
            const double scale = dest[0] / source[band];
            if (scale == 0.) continue;
            x[row] -= scale * x[col];
            for (int k = 1; k < stop - col; ++k)
                dest[k] -= scale * source[k + band];
        }
    }
    for (int row = n - 1; row >= 0; --row) {
        const double* ar = &a[row * stride + band];
        for (int k = 1; k <= band && row + k < n; ++k)
            x[row] -= ar[k] * x[row + k];
        x[row] /= ar[0];
    }
}
} // namespace

extern "C" void advance(int cap, int steps, int save_every, double step, const double* par,
                        const double* forcing, double* out)
{
    const int n = 4 * (cap + 1);
    std::vector<double> previous(n, 0.), current(n, 0.), next(n), base(n * stride, 0.), matrix;
    current.front() = 1.;
    std::copy(current.begin(), current.end(), out);
    const double scale = 2. * step / 3.;
    auto event = [&](int origin, int dest, double rate) {
        base[origin * stride + band] += scale * rate;
        if (dest < n)
            base[dest * stride + origin - dest + band] -= scale * rate;
    };
    for (int m = 0; m <= cap; ++m) {
        const int j = 4 * m;
        for (int state = 0; state < 4; ++state) {
            const int i = j + state;
            base[i * stride + band] = 1.;
            event(i, i + 4, par[7 + state]);
            if (m) event(i, i - 4, m * par[11]);
        }
        event(j, j + 1, par[0]);
        event(j + 1, j + 2, par[3]);
        event(j + 2, j + 1, par[4]);
        event(j + 2, j + 3, par[5]);
        event(j + 3, j + 2, par[6]);
    }
    for (int s = 1; s <= steps; ++s) {
        matrix = base;
        const double backwards = std::max(0., par[1] - par[2] * forcing[s]);
        const double scaled = scale * backwards;
        for (int i = 1; i < n; i += 4) {
            matrix[i * stride + band] += scaled;
            matrix[(i - 1) * stride + band + 1] -= scaled;
        }
        for (int i = 0; i < n; ++i) {
            if (s == 1) {
                // Switch from the BDF2 scale to one backward-Euler startup step.
                for (int k = 0; k < stride; ++k) matrix[i * stride + k] *= 1.5;
                matrix[i * stride + band] -= .5;
                next[i] = current[i];
            } else {
                next[i] = (4. * current[i] - previous[i]) / 3.;
            }
        }
        solve(matrix, next);
        previous.swap(current);
        current.swap(next);
        if (s % save_every == 0)
            std::copy(current.begin(), current.end(), out + (s / save_every) * n);
    }
}

extern "C" void advance_mean(int steps, int save_every, double step, int spatial, const double* par,
                             const double* forcing, double* out)
{
    // The first moments close because counts only appear linearly in the rates.
    constexpr int n = 6;
    std::vector<double> previous(n, 0.), current(n, 0.), next(n), base(n * stride, 0.), matrix;
    const double scale = 2. * step / 3.;
    current[0] = 1.;
    out[0] = out[1] = 0.;
    auto event = [&](int from, int to, double rate) {
        base[from * stride + band] += scale * rate;
        base[to * stride + from - to + band] -= scale * rate;
    };
    for (int k = 0; k < n; ++k) base[k * stride + band] = 1.;
    event(0, 1, par[0]);
    event(1, 2, par[3]);
    event(2, 1, par[4]);
    event(2, 3, par[5]);
    event(3, 2, par[6]);
    for (int k = 0; k < 4; ++k)
        base[4 * stride + k] -= scale * par[7 + k];
    base[4 * stride + band] += scale * (spatial ? par[14] : par[11]);
    if (spatial) {
        event(4, 5, par[13]);
        base[5 * stride + band] += scale * par[11];
    }
    for (int s = 1; s <= steps; ++s) {
        matrix = base;
        const double reverse = scale * std::max(0., par[1] - par[2] * forcing[s]);
        matrix[stride + band] += reverse;
        matrix[band + 1] -= reverse;
        for (int i = 0; i < n; ++i) {
            if (s == 1) {
                for (int k = 0; k < stride; ++k) matrix[i * stride + k] *= 1.5;
                matrix[i * stride + band] -= .5;
                next[i] = current[i];
            } else {
                next[i] = (4. * current[i] - previous[i]) / 3.;
            }
        }
        solve(matrix, next);
        previous.swap(current);
        current.swap(next);
        if (s % save_every == 0) {
            const int at = 2 * (s / save_every);
            out[at] = current[4];
            out[at + 1] = current[5];
        }
    }
}
