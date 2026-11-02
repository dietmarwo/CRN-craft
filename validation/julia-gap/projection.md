| Workload | Fastest projected CRN | Measured speedup | Projected seconds | Fastest Julia | Julia seconds | Projected CRN / Julia |
|---|---|---|---|---|---|---|
| W1 | legacy_cached | 1.05× | 0.0749 | Julia/Direct | 0.0511 | 1.47× |
| W2 | direct_cr | 1.44× | 0.1272 | Julia/DirectCR | 0.0968 | 1.31× |
| W3 | sorting_direct | 2.33× | 0.3278 | Julia/SortingDirect | 0.2944 | 1.11× |
| W4 | rssa_cr | 1.22× | 0.1946 | Julia/RSSACR | 0.1061 | 1.83× |
| W5 | direct_cr | 1.51× | 0.0768 | Julia/DirectCR | 0.0477 | 1.61× |
| W6 | rssa_cr | 1.42× | 0.0614 | Julia/RSSACR | 0.0396 | 1.55× |
| W7 | rssa_cr | 1.27× | 0.0187 | Julia/RSSACR | 0.0133 | 1.41× |

Geometric mean projected gap: 1.46×
