"""CLI entry point for MODIS Nd bias processing."""

import argparse
import glob
import os
from pathlib import Path

from .core import DEFAULT_RESOLUTIONS, process_files


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Process MODIS cloud data to analyze Nd bias."
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default=None,
        help="Directory containing HDF files (can use recursive globbing).",
    )
    parser.add_argument(
        "--input-dirs",
        nargs="+",
        default=None,
        help="One or more directories containing HDF files.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Directory to save processed NetCDF files.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search recursively for HDF files in input directory.",
    )
    parser.add_argument(
        "--plots",
        action="store_true",
        help="Generate plots from processed NetCDF files.",
    )
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Only generate plots (skip processing).",
    )

    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parents[2]
    default_input_dirs = [
        base_dir / "data" / "open",
        base_dir / "data" / "closed",
        base_dir / "data" / "disorganized",
    ]
    output_dir = Path(args.output_dir) if args.output_dir else base_dir / "data" / "processed"

    if args.plots_only:
        try:
            from .plot import generate_all_plots
        except ModuleNotFoundError as exc:
            raise SystemExit(
                "Plotting requires 'matplotlib'. "
                "Install extras with: pip install 'modis-nd-bias[plot]'"
            ) from exc
        output_dir.mkdir(parents=True, exist_ok=True)
        generate_all_plots(output_dir)
        return

    search_pattern = "**/*.hdf" if args.recursive else "*.hdf"

    input_dirs = []
    if args.input_dirs:
        input_dirs = [Path(p) for p in args.input_dirs]
    elif args.input_dir:
        input_dirs = [Path(args.input_dir)]
    else:
        input_dirs = default_input_dirs

    fileset = []
    for d in input_dirs:
        if not d.exists():
            raise SystemExit(f"Input directory not found: {d}")
        fileset.extend(glob.glob(str(d / search_pattern), recursive=args.recursive))

    process_files(fileset, str(output_dir), DEFAULT_RESOLUTIONS)

    if args.plots:
        try:
            from .plot import generate_all_plots
        except ModuleNotFoundError as exc:
            raise SystemExit(
                "Plotting requires 'matplotlib'. "
                "Install extras with: pip install 'modis-nd-bias[plot]'"
            ) from exc
        output_dir.mkdir(parents=True, exist_ok=True)
        generate_all_plots(output_dir)


if __name__ == "__main__":
    main()
