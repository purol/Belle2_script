#!/usr/bin/env python3
"""Plot 19 stacked distributions with one RDataFrame event loop per TTree.

Usage:
    python plot_19_stacks.py /path/to/root_files --outdir selection_plots
    python plot_19_stacks.py /path/to/root_files --recursive

Each ROOT file should contain TTree objects named pizero, gamma, KSzero,
pion, and kaon, each with an isSignal branch. File/tree validation is
deferred to data processing for faster startup with many small ROOT files.
"""

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetHatchesLineWidth(1)


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
                 100, -math.pi / 2.0, math.pi / 2.0,
                 "#Delta#phi_{#gamma#gamma} [rad]"),
        PlotSpec(3, "pizero", angle, all_of(mass_window, phi_window, momentum_cut),
                 100, 0.0, math.pi / 2.0,
                 "#Delta#theta_{#gamma#gamma} [rad]"),
        PlotSpec(4, "pizero", "p", all_of(mass_window, phi_window, angle_cut),
                 100, 0.0, 3.0, "p_{#pi^{0}} [GeV/c]"),
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
                 100, -0.5, 20.5, "clusterNHits"),
        PlotSpec(6, "gamma", "theta", all_of(gamma_hits, gamma_energy),
                 100, 0.0, math.pi, "#theta [rad]"),
        PlotSpec(7, "gamma", "E", all_of(gamma_hits, gamma_theta),
                 100, 0.0, 2.0, "E [GeV]"),
        PlotSpec(8, "KSzero", "M", "significanceOfDistance > 50",
                 100, 0.45, 0.55, "Mass of K_{S}^{0} [GeV/c^{2}]"),
        PlotSpec(9, "KSzero", "significanceOfDistance",
                 all_of("M > 0.4876", "M < 0.5076"),
                 100, 0.0, 500.0, "significanceOfDistance"),
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


def book_histograms(dataframes, plots):
    """Book all actions before reading any results (one loop per TTree)."""
    booked = {}
    for spec in plots:
        selected = dataframes[spec.tree].Filter(spec.selection)
        args = (spec.nbins, spec.xmin, spec.xmax)
        h_bkg = selected.Filter("isSignal < 0.5").Histo1D(
            (f"hist_background_{spec.number}", "", *args), spec.variable
        )
        h_sig = selected.Filter("isSignal >= 0.5").Histo1D(
            (f"hist_signal_{spec.number}", "", *args), spec.variable
        )
        booked[spec.number] = (h_bkg, h_sig)
    return booked


def make_plot(spec, booked_histograms, outdir, pdf_filename, is_first, is_last):
    # All actions are already booked. The first GetValue triggers the event
    # loop for this TTree and fills every booked histogram on that dataframe.
    background = booked_histograms[0].GetValue()
    signal = booked_histograms[1].GetValue()
    n_background = int(background.GetEntries())
    n_signal = int(signal.GetEntries())

    # Use a pale solid base for both components so hatching leaves no white gaps.
    # A second, patterned copy of the background will be overlaid after stacking.
    background.SetFillColor(ROOT.TColor.GetColor("#E9C29A"))
    background.SetFillStyle(1001)
    background.SetLineColor(ROOT.TColor.GetColor("#8A5D40"))
    background.SetLineWidth(1)
    signal.SetFillColor(ROOT.TColor.GetColor("#B8D8EA"))
    signal.SetFillStyle(1001)  # Solid pastel blue
    signal.SetLineColor(ROOT.TColor.GetColor("#376586"))
    signal.SetLineWidth(1)

    stack = ROOT.THStack(
        f"stack_{spec.number}", f";{spec.xlabel};Arbitrary units"
    )
    stack.Add(background)
    stack.Add(signal)

    canvas = ROOT.TCanvas(f"canvas_{spec.number}", "", 920, 600)
    canvas.SetLeftMargin(0.13)
    canvas.SetBottomMargin(0.12)
    canvas.SetRightMargin(0.23)  # Keep the legend outside the histogram frame
    canvas.SetTicks(1, 1)
    log_y = (
        (spec.tree == "KSzero" and spec.variable == "significanceOfDistance")
        or (spec.tree in ("pion", "kaon") and
            spec.variable in ("dr", "dz", "pionID", "kaonID"))
    )
    if log_y:
        canvas.SetLogy()
    stack.Draw("HIST")
    stack.SetMinimum(0.5 if log_y else 0.0)

    # Draw hatching only on the *background portion* of the stack.  This is
    # a visual overlay, not another THStack component or counted histogram.
    background_hatch = background.Clone(f"background_hatch_{spec.number}")
    background_hatch.SetDirectory(0)
    background_hatch.SetFillColor(ROOT.TColor.GetColor("#986745"))
    background_hatch.SetFillStyle(3645)  # Slightly denser diagonal hatching
    background_hatch.SetLineColor(ROOT.TColor.GetColor("#8A5D40"))
    background_hatch.SetLineWidth(1)
    background_hatch.Draw("HIST SAME")

    legend = ROOT.TLegend(0.78, 0.77, 0.98, 0.89)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.033)
    particle_names = {
        "pizero": "#pi^{0}",
        "gamma": "#gamma",
        "KSzero": "K_{S}^{0}",
        "pion": "pion",
        "kaon": "kaon",
    }
    particle = particle_names[spec.tree]
    legend.AddEntry(signal, f"true {particle}", "f")
    legend.AddEntry(background_hatch, f"fake {particle}", "f")
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
    files = sorted(p for p in input_path.glob(pattern) if p.is_file())
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
    parser.add_argument("--threads", type=int, default=1,
                        help="RDataFrame worker threads (default: 1)")
    args = parser.parse_args()

    plots = build_plots()
    input_files = find_input_files(args.input_path, args.recursive)
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.threads < 1:
        parser.error("--threads must be at least 1")
    if args.threads > 1:
        ROOT.EnableImplicitMT(args.threads)

    print(f"Reading {len(input_files)} ROOT file(s) from {args.input_path}")
    # Keep each chain alive while its RDataFrame is in use.
    chains = {}
    dataframes = {}
    for tree_name in ("pizero", "gamma", "KSzero", "pion", "kaon"):
        chain = ROOT.TChain(tree_name)
        for filepath in input_files:
            # Default Add() defers file/tree checks until the event loop.
            chain.Add(str(filepath))
        chains[tree_name] = chain
        dataframes[tree_name] = ROOT.RDataFrame(chain)

    # IMPORTANT: book all 38 lazy Histo1D actions before GetValue().
    booked = book_histograms(dataframes, plots)

    pdf_filename = args.outdir / "all_stacks.pdf"
    for i, spec in enumerate(plots):
        make_plot(spec, booked[spec.number], args.outdir, pdf_filename,
                  is_first=(i == 0), is_last=(i == len(plots) - 1))

    print(f"Saved {len(plots)} PNG files and {pdf_filename}")


if __name__ == "__main__":
    main()
