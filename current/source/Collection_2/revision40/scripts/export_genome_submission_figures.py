"""Export frozen Stage 39 figure objects to submission EPS and 600-dpi TIFF.

No statistics are refitted and no source/plotting code is edited. Only the source
data writer and export callback are replaced at runtime. TIFF comes directly
from the original figure at 600 dpi, never from enlargement of a PNG.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import re
import sys
import warnings

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
import PIL

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
NAMES = {
    "Fig_G2_method_sensitivity": "Fig4",
    "Fig_G4_cross_study_alignment": "Fig5",
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    def scalar(v):
        if isinstance(v, np.generic):
            return v.item()
        raise TypeError(f"Unsupported JSON type: {type(v).__name__}")
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False, default=scalar) + "\n", encoding="utf-8")


def snapshot(fig):
    """Record artist data, limits, texts and geometry before/after serialization."""
    result = {"size_inches": fig.get_size_inches().tolist(), "axes": [], "texts": []}
    for ax in fig.axes:
        item = {
            "position": list(ax.get_position().bounds),
            "xlim": list(ax.get_xlim()), "ylim": list(ax.get_ylim()),
            "xlabel": ax.get_xlabel(), "ylabel": ax.get_ylabel(),
            "title": ax.get_title(loc="left"),
            "lines": [], "collections": [], "images": [],
        }
        for line in ax.lines:
            item["lines"].append({"x": np.asarray(line.get_xdata(), dtype=float).tolist(),
                                  "y": np.asarray(line.get_ydata(), dtype=float).tolist()})
        for coll in ax.collections:
            a = coll.get_array()
            item["collections"].append({
                "offsets": np.asarray(coll.get_offsets(), dtype=float).tolist(),
                "array": None if a is None else np.asarray(a, dtype=float).tolist(),
                "paths_sha256": hashlib.sha256(b"".join(np.asarray(p.vertices).tobytes() for p in coll.get_paths())).hexdigest(),
            })
        for im in ax.images:
            item["images"].append({"array": np.asarray(im.get_array()).tolist(), "extent": list(im.get_extent())})
        result["axes"].append(item)
    from matplotlib.text import Text
    for text in fig.findobj(Text):
        if text.get_text():
            result["texts"].append({"text": text.get_text(), "position": list(text.get_position()), "fontsize": text.get_fontsize()})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage39", type=Path, default=ROOT.parent / "pfhn_genome39_2026-10-13")
    args = parser.parse_args()
    stage = args.stage39.resolve()
    original_script = stage / "scripts/make_genome_figures.py"
    old_figures = stage / "figures"
    output = ROOT / "figures_new"
    qa = ROOT / "review/genome_submission_exports"
    output.mkdir(parents=True, exist_ok=True)
    qa.mkdir(parents=True, exist_ok=True)
    before_script_hash = sha(original_script)
    report = {
        "utc": datetime.now(timezone.utc).isoformat(),
        "command": sys.argv, "python_executable": sys.executable, "python_version": sys.version,
        "matplotlib": mpl.__version__, "numpy": np.__version__, "pandas": pd.__version__, "pillow": PIL.__version__,
        "original_script": str(original_script), "original_script_sha256": before_script_hash,
        "wrapper_sha256": sha(__file__),
        "runtime_overrides": ["source(): compare generated source tables with frozen originals instead of writing", "exported(): export the unchanged figure object to EPS/TIFF and a QA PNG"],
        "statistical_refit": False, "plot_geometry_modified": False,
        "png_upscaling": False, "source_checks": [], "figures": [],
        "eps_visual_render_status": "Not separately rasterized: no existing EPS renderer was identified; no search beyond initial known-path/PATH check and no installation. EPS structure/font checks and direct 600-dpi TIFF visual inspection are documented separately.",
        "visual_inspection": "Pending actual view_image inspection of TIFF and export comparison images.",
    }

    spec = importlib.util.spec_from_file_location("pfhn_frozen_genome_figure_module", original_script)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mpl.rcParams["ps.fonttype"] == 42 and not mpl.rcParams["ps.useafm"]

    def source_check(frame, name):
        reference = old_figures / "source_data" / name
        generated = frame.to_csv(index=False).encode("utf-8")
        frozen = reference.read_bytes()
        match = generated == frozen
        report["source_checks"].append({
            "name": name, "rows": len(frame), "columns": list(frame.columns),
            "reference_sha256": hashlib.sha256(frozen).hexdigest(),
            "generated_sha256": hashlib.sha256(generated).hexdigest(), "byte_identical": match,
        })
        if not match:
            raise AssertionError("Original source-data export differs: " + name)

    def export(fig, name, axes):
        target = NAMES[name]
        fig.canvas.draw()
        before = snapshot(fig)
        write_json(qa / f"{target}.artist_snapshot.json", before)
        align = mod.require_matplotlib_panel_alignment(
            fig, axes=axes, panel_ids=["a", "b"], row_groups=[{"panels": ["a", "b"]}],
            json_out=qa / f"{target}.alignment.json", overlay_svg=qa / f"{target}.alignment.svg",
            tolerance_pt=1.5, gutter_tolerance_pt=1.5, require_panel_labels=True, strict=True)
        reference_png = old_figures / f"{name}.png"
        preview = qa / f"{target}.direct_figure_300dpi.png"
        fig.savefig(preview, dpi=300)
        with Image.open(preview) as new, Image.open(reference_png) as old:
            same = new.size == old.size and np.array_equal(np.asarray(new), np.asarray(old))
            reference_size = list(old.size)
        if not same:
            raise AssertionError(target + ": direct figure rendering differs from the frozen 300-dpi reference")

        eps = output / f"{target}.eps"
        tiff = output / f"{target}.tiff"
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            fig.savefig(eps, format="eps", dpi=600)
            # This buffer is rendered at 600 dpi by the original figure object.
            # Remove only a verified fully opaque alpha channel; never resize.
            buffer = io.BytesIO()
            fig.savefig(buffer, format="tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
            buffer.seek(0)
            with Image.open(buffer) as raw:
                raw.load()
                if "A" in raw.getbands():
                    assert raw.getchannel("A").getextrema() == (255, 255)
                rgb = raw.convert("RGB")
                direct_pixels = np.asarray(rgb).copy()
                rgb.save(tiff, format="TIFF", compression="tiff_lzw", dpi=(600, 600))
        fig.canvas.draw()
        after = snapshot(fig)
        assert before == after, target + ": figure artist data or geometry changed during export"

        eps_text = eps.read_text(encoding="latin1")
        embedded_fonts = re.findall(r"/FontName\s+/([^\s]+)\s+def", eps_text)
        type42 = len(re.findall(r"/FontType\s+42\s+def", eps_text))
        sfnts = len(re.findall(r"/sfnts\s*\[", eps_text))
        bbox = re.search(r"^%%HiResBoundingBox:\s*(.+)$", eps_text, re.M)
        assert eps_text.startswith("%!PS-Adobe-3.0 EPSF-3.0")
        assert embedded_fonts and type42 > 0 and sfnts == type42
        assert not re.search(r"^%%IncludeResource:\s*font", eps_text, re.M)
        with Image.open(tiff) as rendered:
            rendered.load()
            size = list(rendered.size)
            dpi = tuple(float(v) for v in rendered.info["dpi"])
            assert rendered.mode == "RGB" and dpi == (600.0, 600.0)
            assert np.array_equal(np.asarray(rendered), direct_pixels)
            compression = rendered.info.get("compression")
            # view_image does not accept this TIFF encoding. Decode the final
            # TIFF losslessly for viewing; no resampling or image alteration.
            tiff_preview = qa / f"{target}.decoded_final_tiff_600dpi.png"
            rendered.save(tiff_preview, format="PNG", dpi=(600, 600))
        with Image.open(tiff_preview) as decoded_preview:
            assert np.array_equal(np.asarray(decoded_preview), direct_pixels)
        expected_size = [int(v * 600) for v in fig.get_size_inches()]
        assert size == expected_size
        item = {
            "target": target, "source_figure": name,
            "physical_size_mm": (fig.get_size_inches() * 25.4).tolist(),
            "original_artist_data_and_geometry_unchanged": True,
            "reference_300dpi_pixel_identical": same,
            "reference_300dpi_size_pixels": reference_size,
            "reference_png_sha256": sha(reference_png), "export_comparison_png_sha256": sha(preview),
            "panel_alignment": align.get("verdict"),
            "eps": {"path": str(eps), "sha256": sha(eps), "bytes": eps.stat().st_size,
                    "embedded_fonts": embedded_fonts, "font_type_42_count": type42, "sfnt_program_count": sfnts,
                    "external_font_include_requests": 0, "hires_bounding_box": bbox.group(1) if bbox else None,
                    "separate_visual_rasterization": "not performed; no installed renderer available"},
            "tiff": {"path": str(tiff), "sha256": sha(tiff), "bytes": tiff.stat().st_size,
                     "mode": "RGB", "dpi": list(dpi), "size_pixels": size, "compression": compression,
                     "source": "original figure object rendered directly at 600 dpi", "pixel_values_preserved_after_opaque_alpha_removal": True},
            "tiff_visual_preview": {"path": str(tiff_preview), "sha256": sha(tiff_preview),
                                    "pixel_identical_to_final_tiff": True, "resampled": False},
            "export_warnings": [str(w.message) for w in caught],
        }
        report["figures"].append(item)
        plt.close(fig)
        print(json.dumps({"figure": target, "status": "exported and automatic checks passed", "size_pixels": size}), flush=True)

    mod.source = source_check
    mod.exported = export
    # Call only original plotting functions, never original main/statistical scripts.
    mod.fig_g2()
    mod.fig_g4()
    assert sha(original_script) == before_script_hash
    report["input_hashes"] = mod.INPUTS
    for rel, digest in mod.INPUTS.items():
        assert sha(stage / rel) == digest
    report["original_plotting_code_and_result_inputs_unchanged"] = True
    report["status"] = "Exports and automatic checks passed; actual TIFF inspection pending"
    write_json(qa / "export_execution.json", report)


if __name__ == "__main__":
    main()
