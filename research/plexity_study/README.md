# Plexity study (research, not maintained)

In-silico study of primer-dimer potential, capacity, off-target priming and (modelled) uniformity
for a tiled design, 10 to 935 targets. **This is a research artefact, not part of the tool**: it is
not linted, not in CI, and depends on internals of `design.py` (`priming_footprints`,
`cross_product_size`, `process_single_target`, `load_references`) and `scripts/design_viewer.py`,
so it may drift. The results are summarised, with all caveats, in the main README (section 4).

* `plexity_study.py`: stages `candidates`, `select` (`--repair`), `report`, `offtarget`
* `plexity_plots.py`: draws the figures; `figures/` holds the five PNGs and `plexity_study_table.csv`
* Last run against commit `7697559` (E. coli K-12, 935 designable tiles, 4 cores, ~40 min)

```bash
python research/plexity_study/plexity_study.py candidates --workdir study --n-targets 1000
python research/plexity_study/plexity_study.py select --workdir study
python research/plexity_study/plexity_study.py select --workdir study --repair
python research/plexity_study/plexity_study.py report --workdir study
python research/plexity_study/plexity_study.py offtarget --workdir study
python research/plexity_study/plexity_plots.py --workdir study --outdir research/plexity_study/figures
```

## Why keep it: predicting primer-dimer potential
The pair-occupancy model here (nearest-neighbour 3'-end dG -> `c/(c+Kd)` at a stated temperature and
primer concentration, summed into per-primer sequestration) is the most useful piece. It is a
candidate for promotion into `design.py` as a per-pool dimer report. Before that it needs: a
sensitivity analysis on the occupancy threshold (0.42% is calibrated to the tool's own -6 kcal/mol
heuristic, not to data), and validation against wet-lab dimer measurements.
