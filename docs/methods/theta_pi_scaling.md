# θπ scaling — per-site vs raw window sum

## Why this matters

ANGSD's `pestPG` `tP` column is the **sum of per-site θπ across sites with
data in the window**, not a per-site density. To compare windows or
samples you must divide by `nSites` (the per-window callable-site count,
last column of pestPG).

For the 226-sample 9× cohort, `nSites` varies sample-to-sample due to
coverage dropouts. Without the division, samples with more callable sites
in a window falsely appear more diverse — a coverage artefact that
masquerades as biology.

## References

- ANGSD GitHub issue #329
- Korunes & Samuk 2021 (pixy) — the missing-data-aware denominator
- Direct trace of `thetaStat do_stat` source — `slice()` is a plain sum

## Convention in this repo

`Module 03_theta_pi/` always emits the per-site value as `theta_pi` and
preserves the raw window sum as `tP_sum` in a separate column. Downstream
analyses use `theta_pi`. The raw sum is kept for diagnostic display
(low-mappability regions, repeat-masked stretches, aggressive filter
regimes show pronounced dips in `tP_sum` that the per-site track flattens
by construction).

```
sample  chrom  window_idx  start_bp  end_bp  theta_pi  tP_sum  n_sites
                                              ↑
                                       per-site (tP / nSites)
```
