**🇺🇸 English** · [🇧🇷 Português](README.pt-BR.md)

# gmx_mmpbsa_plot

MM/PBSA plots for **gmx_MMPBSA** users (Amber force fields; an alternative to `gmx_MMPBSA_ana`).
Developed at the **LMM (Laboratório de Modelagem Molecular), UEFS**.

The script scans the subfolders, detects the files of each system and opens a menu to
choose what to plot. At the end it writes a `.txt` report with your choices and the
numerical results.

## Requirements

Python 3.8+ with `numpy`, `pandas` and `matplotlib`.

```bash
pip install numpy pandas matplotlib
```

## Expected layout

One subfolder per system (replicate, ligand, etc.), containing the `gmx_MMPBSA` output files:

```
my_project/
├── system_A/
│   ├── FINAL_RESULTS_MMPBSA.csv   # per-frame energies (gmx_MMPBSA ... -eo FINAL_RESULTS_MMPBSA.csv)
│   └── FINAL_DECOMP_MMPBSA.csv    # per-residue decomposition (-deo, with &decomp in the input)
│                                  # (or FINAL_DECOMP_MMPBSA.dat, already averaged)
├── system_B/
│   └── ...
└── gmx_mmpbsa_plot.py
```

The averages ± SD are computed from the per-frame energies (SD between frames). Not every folder
needs all files: the menu only offers what exists. Energies are already in kcal/mol (gmx_MMPBSA's
unit). The `Frame #` column is the trajectory frame, so the menu asks **how many ns one frame
represents** (`--ns-per-frame`, e.g. `0.01` for 10 ps) to build the time axis. If the CSV has both
GB and PB, choose with `--model GB|PB`. The ligand (`L:` in gmx_MMPBSA) is detected automatically.

## Usage

```bash
python gmx_mmpbsa_plot.py                 # interactive menu, current folder
python gmx_mmpbsa_plot.py /path/to/data   # another root folder
```

The first menu asks for the **language** (Português do Brasil or English US). The choice
applies to everything: menus and messages, plot texts, table headers and the `.txt`
report. In Portuguese, axis numbers use a decimal comma. The last language used is saved
in `plot_style.ini`. To skip the question: `--lang pt` or `--lang en`.
(The `--help` text is always in English.)

Then the menu asks, in this order:

1. **Systems** to include (`1,3`, `1-3` or `all`).
2. **Label** of each system. The default is the folder name; Enter keeps it.
3. **Model**: MM/GBSA (Generalized Born) or MM/PBSA (Poisson-Boltzmann). The script reads
   which one is in the CSV and suggests it; choosing one that is not in the file gives an
   error. With GBSA, the plots also include **ΔG gas (GGAS)** and **ΔG solv (GSOLV)**.
4. **Plots**: energy vs. time (total, vdW, electrostatic, polar, nonpolar; with GBSA also
   gas and solvation), averages ± SD, per-residue decomposition.
5. **Time per frame** in ns and **time range** (see "About the time axis" below).
6. **Residues**: Top 5, Top 10 or All. *All is not recommended*: with hundreds of
   residues the bar chart is unreadable, so only the sequence-position profile and the
   complete CSV are generated.
7. **Ligand**: the ligand is detected from the `L:` prefix of gmx_MMPBSA; the menu asks
   whether to keep it in the ranking. If kept, it usually dominates the top.
8. **Colors**: comma-separated hex (`#0072B2,#D55E00`) or Enter for the default palette.
9. **Organization** of the output folders (see below).
10. **Formats** (png, pdf, svg), CSV export of the data, and the interactive window.

### About the time axis

The gmx_MMPBSA CSV has only the **trajectory frame number** (`Frame #`), with no time,
so the script asks how many ns one trajectory frame represents (`--ns-per-frame`).
Example: `dt = 0.002` ps and `nstxout-compressed = 5000` give 10 ps per frame, i.e.
`0.01` ns; frame 50000 is then 500 ns. This is the time of one frame of the *original*
trajectory, not the gap between two rows of the CSV (the `interval` is already in the
frame numbers). The question only appears when an energy-vs-time plot is chosen.

### Output

By default everything goes to a new folder `MM_PBSA_plots_YYYYMMDD_HHMMSS/` created on
every run (so nothing is ever overwritten), organized by plot type:

```
MM_PBSA_plots_YYYYMMDD_HHMMSS/
├── energy_vs_time/     # energy over time (one plot per component)
├── summary/            # component averages ± SD
├── per_residue/        # per-residue decomposition (one per system)
├── tables/             # CSVs: rankings, "mean ± SD" tables, plot data
└── GB_mmpbsa_report_YYYYMMDD_HHMMSS.txt
```

The **file names inside** get the method prefix: `GB_` or `PB_`
(e.g. `GB_energy_total_vs_time.png`, `GB_mmpbsa_report_....txt`). The folder name is unchanged.

Use `-o folder_name` to choose the output folder yourself (an existing folder is reused,
and files with the same name are overwritten).

In Portuguese the folders are named `energia_vs_tempo/`, `medias/`, `por_residuo/` and
`tabelas/`. Other layouts are available: **by file format** (`png/`, `pdf/`, `svg/`,
handy when submitting) or **everything in one folder**.

### What the averages plot shows (`energy_components_mean`)

Bars with the mean of each binding-energy term, per system, read from
`FINAL_RESULTS_MMPBSA.csv`: vdW, electrostatic, polar solvation, nonpolar solvation and total
(total = sum of the four). With **MM/GBSA** the bars also include **ΔG gas** (GGAS = vdW +
electrostatic) and **ΔG solv** (GSOLV = polar + nonpolar), as usually reported in GBSA papers. The average is over the **whole trajectory** and the error bar
is the **standard deviation across frames** (the energy fluctuation), not the uncertainty
of the mean. State this in the figure caption. The entropic term (−TΔS) is not included.

## Editing labels, text and font

On the first run a `plot_style.ini` is created in the root folder. Edit it and run again:

```ini
[names]
# accepts matplotlib mathtext (comments only on their own line)
system_A = $\it{Mtb}$PS - KS1

[text_en]
energy_ylabel = ENERGY (kcal/mol)
residue_title = Per-residue energy decomposition
comp_vdw = vDW

[font]
family = Arial
size = 12
```

The texts live in two blocks, `[text_en]` and `[text_pt]`, one per language.
A `plot_style.ini` from earlier versions is upgraded automatically on the first run
(the labels you defined are kept).
The default font is **Arial**. If it is not installed (common on Linux), the script uses
Liberation Sans, which has the same metrics, and records this in the report.

For fine-tuning outside the script:

- **SVG and PDF** keep the text editable. Open them in Inkscape or Illustrator
  (for SVG the font must be installed on the machine that opens it).
- **Interactive window** (menu option): the *Figure options* button in the matplotlib
  toolbar changes title, labels, limits and colors before saving.
- **CSV** of the data: to redo the plot in Origin, Excel or Grace.

## Running without the menu

Any option skips the corresponding question. With `-y` nothing is asked: it uses the
defaults (all systems, all formats, ligand dropped), so `--plots` is required, and so are
`--time`, `--ns-per-frame` and `--top` when those plots are selected. `--model GB|PB` skips the
model question.

```bash
python gmx_mmpbsa_plot.py . --systems 1,3 --plots total,vdw,residues \
    --time 50 500 --ns-per-frame 0.01 --model GB --top 10 --ligand drop --formats png,svg --csv \
    --lang en --layout type -y
```

`python gmx_mmpbsa_plot.py --help` lists all options.

## Notes

- **Variables by model.** Polar = `EGB` (GB) or `EPB` (PB); nonpolar = `ESURF` (GB) or
  `ENPOLAR` + `EDISPER` / `ECAVITY` + `EDISPER` (PB). The script warns if the total differs
  from the sum of the four components (non-zero bonded/1-4 terms in the delta).
- **Decomposition**: besides `FINAL_DECOMP_MMPBSA.csv` (per frame), the script accepts
  `FINAL_DECOMP_MMPBSA.dat` (mean and SD already aggregated by gmx_MMPBSA). If both exist,
  the CSV is used. The `.dat` header tells whether `idecomp = 2` (1-4 terms added to vdW and
  electrostatic); mention it in the paper. Only the residues in `print_res` appear.

- The averages plot uses `FINAL_RESULTS_MMPBSA.csv`, i.e. the **whole trajectory**. The averages
  in the report for the chosen time range are computed from the per-frame CSV and may differ.
- **Per-residue decomposition and the time range.** The chosen time range also applies to
  the per-residue decomposition: with `FINAL_DECOMP_MMPBSA.csv` (per frame), mean and SD of
  each residue use only the frames inside the range (time = `Frame #` × time per frame, as in
  the energy-vs-time plots). `FINAL_DECOMP_MMPBSA.dat` is already aggregated and cannot be cut
  by time: in that case the whole trajectory is used and a warning is shown. The report states
  the range used for the residues.
- The `*_top_N_table.csv` file has the values in `mean ± SD` format.

## Authors

- Késsia Souza Santos: LMM (Laboratório de Modelagem Molecular), UEFS
