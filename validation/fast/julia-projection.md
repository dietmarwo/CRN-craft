| Workload | Fastest projected CRN | Measured speedup | Projected seconds | Fastest Julia | Julia seconds | Projected CRN / Julia |
|---|---|---|---|---|---|---|
| W1 | cached_direct | 1.44× | 0.0493 | Julia/Direct | 0.0506 | 0.97× |
| W2 | cached_direct | 2.41× | 0.0849 | Julia/RSSACR | 0.0980 | 0.87× |
| W3 | cached_direct | 1.63× | 0.2322 | Julia/SortingDirect | 0.2920 | 0.80× |
| W4 | rssa_cr | 1.26× | 0.1137 | Julia/RSSACR | 0.1072 | 1.06× |
| W5 | direct_cr | 1.43× | 0.0580 | Julia/DirectCR | 0.0484 | 1.20× |
| W6 | rssa_cr | 1.27× | 0.0544 | Julia/RSSACR | 0.0388 | 1.40× |
| W7 | cached_direct | 3.19× | 0.0133 | Julia/RSSACR | 0.0138 | 0.96× |

Geometric mean projected gap: 1.02×
