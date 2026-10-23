#!/usr/bin/env python3
"""Plot 19 signal/background stacked, N-1-style ROOT distributions.

Usage:
    python plot_19_stacks.py /path/to/root_files --outdir selection_plots
    python plot_19_stacks.py /path/to/root_files --recursive

Each ROOT file must contain TTree objects named pizero, gamma, KSzero,
pion, and kaon, each with an isSignal branch. A single ROOT file is also
accepted for backwards compatibility.
"""

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


@dataclass(frozen=True)
class PlotSpec:
    number: int
    tree: str
    variable: str
    selection: str
    nbins: int
    xmin: float
    xmax: float
    xlabel: str


def all_of(*expressions):
    """Combine ROOT selection expressions without Python-style chained cuts."""
    return " && ".join(f"({expr})" for expr in expressions)


def build_plots():
    phi = "daughterDiffOf__bo0__cm1__cmphi__bc"
    angle = "daughterAngle__bo0__cm1__bc"

    mass_window = all_of("InvM > 0.120", "InvM < 0.145")
    phi_window = all_of(f"{phi} > -1.5", f"{phi} < 1.5")
    angle_cut = f"{angle} < 1.4"
    momentum_cut = "p > 0.4"

    plots = [
        PlotSpec(1, "pizero", "InvM", all_of(phi_window, angle_cut, momentum_cut),
                 100, 0.08, 0.20, "M_{#gamma#gamma} [GeV/c^{2}]"),
        PlotSpec(2, "pizero", phi, all_of(mass_window, angle_cut, momentum_cut),
                 100, -math.pi, math.pi, "#Delta#phi [rad]"),
        PlotSpec(3, "pizero", angle, all_of(mass_window, phi_window, momentum_cut),
                 100, 0.0, math.pi, "Daughter opening angle [rad]"),
        PlotSpec(4, "pizero", "p", all_of(mass_window, phi_window, angle_cut),
                 100, 0.0, 5.0, "#pi^{0} momentum [GeV/c]"),
    ]

    gamma_energy = (
        "((clusterReg == 1 && E > 0.080) || "
        "(clusterReg == 2 && E > 0.030) || "
        "(clusterReg == 3 && E > 0.060))"
    )
    gamma_theta = "thetaInCDCAcceptance > 0.5"
    gamma_hits = "clusterNHits > 1.5"

    plots += [
        PlotSpec(5, "gamma", "clusterNHits", all_of(gamma_theta, gamma_energy),
                 100, -0.5, 99.5, "Cluster N Hits"),
        PlotSpec(6, "gamma", "theta", all_of(gamma_hits, gamma_energy),
                 100, 0.0, math.pi, "Cluster #theta [rad]"),
        PlotSpec(7, "gamma", "E", all_of(gamma_hits, gamma_theta),
                 100, 0.0, 3.0, "Cluster energy [GeV]"),
        PlotSpec(8, "KSzero", "M", "significanceOfDistance > 50",
                 100, 0.45, 0.55, "M_{K_{S}^{0}} [GeV/c^{2}]"),
        PlotSpec(9, "KSzero", "significanceOfDistance",
                 all_of("M > 0.4876", "M < 0.5076"),
                 100, 0.0, 500.0, "Significance of distance"),
    ]

    tracking_binning = {
        "dr": (100, -1.0, 1.0, "dr [cm]"),
        "dz": (100, -5.0, 5.0, "dz [cm]"),
        "nCDCHits": (80, -0.5, 79.5, "N CDC hits"),
        "nPXDHits": (6, -0.5, 5.5, "N PXD hits"),
    }

    for tree_name, pid, first_number in (
        ("pion", "pionID", 10),
        ("kaon", "kaonID", 15),
    ):
        # For each plot remove ONLY the cut on its displayed variable.
        basic_cuts = {
            "dr": "dr < 0.5",
            "dz": all_of("dz > -2", "dz < 2"),
            "nCDCHits": "nCDCHits > 20",
            "nPXDHits": "nPXDHits > 0",
            pid: f"{pid} > 0.6",
        }
        binning = {
            **tracking_binning,
            pid: (50, 0.0, 1.0, pid),
        }

        for offset, variable in enumerate(basic_cuts):
            selection = all_of(*(
                expression for name, expression in basic_cuts.items()
                if name != variable
            ))
            nbins, xmin, xmax, xlabel = binning[variable]
            plots.append(PlotSpec(first_number + offset, tree_name, variable,
                                  selection, nbins, xmin, xmax, xlabel))

    assert [p.number for p in plots] == list(range(1, 20))
    return plots


def make_plot(tree, spec, outdir, pdf_filename, is_first, is_last):
    background_name = f"hist_background_{spec.number}"
    signal_name = f"hist_signal_{spec.number}"
    args = (spec.nbins, spec.xmin, spec.xmax)
    background = ROOT.TH1D(background_name, "", *args)
    signal = ROOT.TH1D(signal_name, "", *args)
    background.Sumw2()
    signal.Sumw2()

    cut_background = all_of(spec.selection, "isSignal < 0.5")
    cut_signal = all_of(spec.selection, "isSignal >= 0.5")

    # '+': keep the pre-created binning and fill the existing histogram.
    n_background = tree.Draw(
        f"{spec.variable}>>+{background_name}", cut_background, "goff"
    )
    n_signal = tree.Draw(
        f"{spec.variable}>>+{signal_name}", cut_signal, "goff"
    )
    if n_background < 0 or n_signal < 0:
        raise RuntimeError(
            f"TTree.Draw failed: {spec.tree}.{spec.variable}. "
            "Check that all required branches exist."
        )

    background.SetFillColor(ROOT.kBlue - 9)
    background.SetLineColor(ROOT.kBlue + 1)
    signal.SetFillColor(ROOT.kRed - 9)
    signal.SetLineColor(ROOT.kRed + 1)

    stack = ROOT.THStack(
        f"stack_{spec.number}", f";{spec.xlabel};Entries / bin"
    )
    stack.Add(background)
    stack.Add(signal)

    canvas = ROOT.TCanvas(f"canvas_{spec.number}", "", 800, 600)
    canvas.SetLeftMargin(0.13)
    canvas.SetBottomMargin(0.12)
    stack.Draw("HIST")
    stack.SetMinimum(0.0)

    legend = ROOT.TLegend(0.62, 0.75, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(background, f"Background ({n_background})", "f")
    legend.AddEntry(signal, f"Signal ({n_signal})", "f")
    legend.Draw()

    output = outdir / f"{spec.number:02d}_{spec.tree}_{spec.variable}_stack.png"
    canvas.SaveAs(str(output))

    # ROOT uses '(' and ')' to open/close a multipage PDF.
    pdf_suffix = "(" if is_first else ")" if is_last else ""
    canvas.Print(str(pdf_filename) + pdf_suffix)

    print(
        f"[{spec.number:02d}/19] {spec.tree}.{spec.variable}: "
        f"background={n_background}, signal={n_signal}, "
        f"total={n_background + n_signal}"
    )


def find_input_files(input_path, recursive=False):
    """Find ROOT files in the requested directory, or accept a single file."""
    if input_path.is_file():
        if input_path.suffix.lower() != ".root":
            raise ValueError(f"Input file is not a .root file: {input_path}")
        return [input_path.resolve()]
    if not input_path.is_dir():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")

    pattern = "**/*.root" if recursive else "*.root"
    files = sorted(p.resolve() for p in input_path.glob(pattern) if p.is_file())
    if not files:
        raise FileNotFoundError(f"No .root files found in: {input_path}")
    return files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_path", type=Path,
                        help="Directory containing ROOT files, or a single ROOT file")
    parser.add_argument("--recursive", action="store_true",
                        help="Also read .root files in subdirectories")
    parser.add_argument("--outdir", type=Path, default=Path("selection_plots"),
                        help="Directory for 19 PNG files and all_stacks.pdf")
    args = parser.parse_args()

    plots = build_plots()
    input_files = find_input_files(args.input_path, args.recursive)
    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading {len(input_files)} ROOT file(s) from {args.input_path}")
    trees = {}
    for tree_name in ("pizero", "gamma", "KSzero", "pion", "kaon"):
        chain = ROOT.TChain(tree_name)
        for filepath in input_files:
            # nentries=0 forces ROOT to check the file and the named tree.
            if chain.Add(str(filepath), 0) == 0:
                raise RuntimeError(f"Missing or unreadable '{tree_name}' in {filepath}")
        trees[tree_name] = chain
        print(f"  {tree_name}: {chain.GetEntries()} entries")

    pdf_filename = args.outdir / "all_stacks.pdf"
    for i, spec in enumerate(plots):
        make_plot(trees[spec.tree], spec, args.outdir, pdf_filename,
                  is_first=(i == 0), is_last=(i == len(plots) - 1))

    print(f"Saved {len(plots)} PNG files and {pdf_filename}")


if __name__ == "__main__":
    main()
