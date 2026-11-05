#!/usr/bin/env python3
"""Create a config copy with generated line-sweep values applied."""

import argparse
import json
from pathlib import Path
import re


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SECTION_SEPARATOR = "=========================================================="
GEOMETRY_PREFIX = "CAGeometryHostESProducer"
ANALYZER_PREFIX = "SimPixelTrackAnalyzer"
CA_HIT_NTUPLET_PREFIX = "CAHitNtuplet"
SUMMARY_GEOMETRY_KEYS = {
    "phiCuts", "minZ", "maxZ", "maxR", "z0cut", "thetaCuts", "dcaCuts",
}


def read_edge_sections(path):
    sections = path.read_text().split(SECTION_SEPARATOR)
    if len(sections) != 2:
        raise ValueError(
            f"{path} must have two sections separated by {SECTION_SEPARATOR!r}"
        )

    values = []
    for index, section in enumerate(sections, start=1):
        try:
            values.append(json.loads("{" + section.strip().rstrip(",") + "}"))
        except json.JSONDecodeError as error:
            raise ValueError(
                f"Invalid JSON-like data in edge file section {index}: {error}"
            ) from error

    geometry, analyzer = values
    geometry_keys = {
        "nLayers", "nPairs", "nModules", "pairGraph", "startingPairs",
        "layerStarts", "phiCuts", "minZ", "maxZ", "maxR", "z0cut",
        "thetaCuts", "dcaCuts",
    }
    analyzer_keys = {
        "nPairs", "nLayers", "pairGraph", "startingPairs", "isBarrel", "ptCuts"
    }
    for name, section, required in (
        ("geometry", geometry, geometry_keys),
        ("analyzer", analyzer, analyzer_keys),
    ):
        missing = required - section.keys()
        if missing:
            raise ValueError(
                f"Edge file {name} section is missing: {', '.join(sorted(missing))}"
            )

    n_pairs = geometry["nPairs"]
    n_layers = geometry["nLayers"]
    if geometry["nModules"] != n_layers:
        raise ValueError("Edge file nModules and nLayers values differ")
    if analyzer["nPairs"] != n_pairs or analyzer["nLayers"] != n_layers:
        raise ValueError("Edge file geometry and analyzer dimensions differ")
    if len(geometry["pairGraph"]) != 2 * n_pairs:
        raise ValueError("Geometry pairGraph length must be twice nPairs")
    if len(analyzer["pairGraph"]) != 2 * n_pairs:
        raise ValueError("Analyzer pairGraph length must be twice nPairs")
    if geometry["pairGraph"] != analyzer["pairGraph"]:
        raise ValueError("Geometry and analyzer pairGraph values differ")
    if geometry["startingPairs"] != analyzer["startingPairs"]:
        raise ValueError("Geometry and analyzer startingPairs values differ")
    if len(geometry["layerStarts"]) != n_layers + 1:
        raise ValueError("Geometry layerStarts length must be nLayers + 1")

    for section, keys, expected in (
        (
            geometry,
            ("startingPairs", "phiCuts", "minZ", "maxZ", "maxR", "z0cut"),
            n_pairs,
        ),
        (geometry, ("thetaCuts", "dcaCuts"), n_layers),
        (analyzer, ("startingPairs", "ptCuts"), n_pairs),
        (analyzer, ("isBarrel",), n_layers),
    ):
        for key in keys:
            if len(section[key]) != expected:
                raise ValueError(
                    f"Edge file {key} has {len(section[key])} values; "
                    f"expected {expected}"
                )

    return geometry, analyzer

def read_sim_pixel_tracks_summary(path, n_pairs, n_layers):
    sections = [
        section.strip()
        for section in re.split(r"(?m)^=+\s*$", path.read_text())
        if section.strip()
    ]
    if len(sections) != 2:
        raise ValueError(
            f"{path} must contain parameter and hardCurvCut sections "
            "separated by lines of '='"
        )

    parsed = []
    for index, section in enumerate(sections, start=1):
        try:
            parsed.append(json.loads("{" + section.rstrip(",") + "}"))
        except json.JSONDecodeError as error:
            raise ValueError(
                f"Invalid summary section {index} in {path}: {error}"
            ) from error

    geometry_values, hit_ntuplet_values = parsed
    missing = SUMMARY_GEOMETRY_KEYS - geometry_values.keys()
    if missing:
        raise ValueError(
            f"Summary is missing geometry values: {', '.join(sorted(missing))}"
        )
    if "hardCurvCut" not in hit_ntuplet_values:
        raise ValueError("Summary is missing hardCurvCut")

    for key in ("phiCuts", "minZ", "maxZ", "maxR", "z0cut"):
        if len(geometry_values[key]) != n_pairs:
            raise ValueError(
                f"Summary {key} has {len(geometry_values[key])} values; "
                f"edge data has {n_pairs} pairs"
            )
    for key in ("thetaCuts", "dcaCuts"):
        if len(geometry_values[key]) != n_layers:
            raise ValueError(
                f"Summary {key} has {len(geometry_values[key])} values; "
                f"edge data has {n_layers} layers"
            )

    return geometry_values, hit_ntuplet_values


def find_module(config, prefix, config_path):
    modules = [
        value
        for name, value in config.items()
        if name.startswith(prefix) and isinstance(value, dict)
    ]
    if len(modules) != 1:
        raise ValueError(
            f"Expected exactly one {prefix} module in {config_path}; "
            f"found {len(modules)}"
        )
    return modules[0]


def compact_array(match):
    try:
        values = json.loads("[" + match.group(1) + "]")
    except json.JSONDecodeError:
        return match.group(0)

    # Keep nested arrays or objects in their existing format.
    if any(isinstance(value, (list, dict)) for value in values):
        return match.group(0)

    return "[" + ", ".join(json.dumps(value) for value in values) + "]"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create a new config from a template and apply the generated "
            "line-sweep values to its geometry and analyzer modules."
        )
    )
    parser.add_argument("-c", "--config", type=Path, help="Template JSON config")
    parser.add_argument(
        "-i", "--input",
        type=Path,
        default= "",
        help="Autograph or SimPixelTracks output (default is left empty; the user has to provide an input file)",
    )
    parser.add_argument(
        "-s", "--simPixelTracks",
        action="store_true",
        help="Interpret --input as a SimPixelTracks summary instead of Autograph output",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        help="Output config (default: <template>_with_edges.json)",
    )
    parser.add_argument(
        "-f", "--force",
        action="store_true",
        help="Overwrite the output file if it already exists",
    )
    args = parser.parse_args()

    if not args.config:
        parser.error("--config is required")
    if not str(args.input):
        parser.error("--input is required")

    output = args.output or args.config.with_name(
        f"{args.config.stem}_with_new_parameters{args.config.suffix}"
    )
    protected_template = (args.config.parent / "templateConfig.json").resolve()
    if output.resolve() == protected_template:
        parser.error("Output must not overwrite templateConfig.json")
    if output.exists() and not args.force:
        parser.error(f"Output already exists: {output} (use --force to overwrite)")
    if not output.parent.is_dir():
        parser.error(f"Output directory does not exist: {output.parent}")

    try:
        config = json.loads(args.config.read_text())
        geometry_module = find_module(config, GEOMETRY_PREFIX, args.config)
        analyzer_module = find_module(config, ANALYZER_PREFIX, args.config)

        if args.simPixelTracks:
            # The template supplies the expected pair/layer counts in this mode.
            summary_geometry, summary_hit_ntuplet = read_sim_pixel_tracks_summary(
                args.input,
                geometry_module["nPairs"],
                geometry_module["nLayers"],
            )
            geometry_module.update(summary_geometry)

            hard_curv_cut = summary_hit_ntuplet["hardCurvCut"]
            analyzer_module["hardCurvCut"] = hard_curv_cut
            find_module(config, CA_HIT_NTUPLET_PREFIX, args.config)[
                "hardCurvCut"
            ] = hard_curv_cut
        else:
            geometry, analyzer = read_edge_sections(args.input)
            geometry_module.update(geometry)
            analyzer_module.update(analyzer)

        serialized_config = json.dumps(config, indent=2)
        serialized_config = re.sub(
            r"\[\s*([^\[\]]*?)\s*\]",
            compact_array,
            serialized_config,
            flags=re.DOTALL,
        )
        serialized_config = re.sub(
            r"(\],)\n(\s*\"[^\"\n]+\":)",
            r"\1\n\n\2",
            serialized_config,
        )
        output.write_text(serialized_config + "\n")
    except (OSError, json.JSONDecodeError, ValueError, KeyError) as error:
        parser.error(str(error))

    print(f"Wrote new config: {output}")


if __name__ == "__main__":
    main()
