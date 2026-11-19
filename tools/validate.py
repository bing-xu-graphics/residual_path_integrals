import argparse
import json
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET


def main():
    parser = argparse.ArgumentParser(description="Run small renderer regression checks.")
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("output/validation"))
    parser.add_argument("--strict-identity", action="store_true",
                        help="Also require samplewise cancellation for an identity dynamic transform.")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    suffix = ".exe" if (args.bin_dir / "glacier.exe").exists() else ""
    renderer = (args.bin_dir / ("glacier" + suffix)).resolve()
    compare = (args.bin_dir / ("image_compare" + suffix)).resolve()
    scene = root / "scenes/cornell-box/scene.xml"
    commands = []

    def run(command):
        command = [str(part) for part in command]
        result = subprocess.run(command, capture_output=True, text=True)
        commands.append(dict(command=command, returncode=result.returncode,
                             stdout=result.stdout, stderr=result.stderr))
        (output / "commands.json").write_text(json.dumps(commands, indent=2))
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)
        return result.stdout

    def render(name, source=scene, threads=1, mode="residual", frame="new", spp=64):
        destination = output / (name + ".exr")
        run([renderer, source, "--resolution", 32, "--spp", spp,
             "--integrator", mode, "--frame", frame, "-t", threads, "-o", destination])
        return destination

    def metrics(first, second, destination=None):
        command = [compare, first, second]
        if destination:
            command.append(destination)
        text = run(command)
        return {key: float(value) for key, value in
                re.findall(r"(rmse|max_abs|min|max)=([-+0-9.eE]+)", text)}

    def variant(name, remove_dynamic=False):
        tree = ET.parse(scene)
        for entry in tree.findall(".//string[@name='filename']"):
            entry.set("value", str((scene.parent / entry.get("value")).resolve()))
        for shape in tree.findall("./shape"):
            dynamic = shape.find("dynamic")
            if dynamic is not None:
                if remove_dynamic:
                    shape.remove(dynamic)
                else:
                    transform = dynamic.find("transform")
                    for child in list(transform):
                        transform.remove(child)
                    ET.SubElement(transform, "translate", x="0", y="0", z="0")
        destination = output / (name + ".xml")
        tree.write(destination)
        return destination

    residual = render("residual")
    repeat = render("residual-repeat")
    threaded = render("residual-threaded", threads=4)
    unchanged = render("unchanged", source=variant("unchanged"))
    static = render("static", source=variant("static", remove_dynamic=True))
    new = render("pt-new", mode="path", spp=1024)
    old = render("pt-old", mode="path", frame="old", spp=1024)
    difference = output / "pt-difference.exr"
    signed = metrics(new, old, difference)
    zero = output / "zero.exr"
    metrics(new, new, zero)
    summary = dict(repeat=metrics(residual, repeat),
                   threads=metrics(residual, threaded),
                   static=metrics(static, zero),
                   unchanged=metrics(unchanged, static),
                   residual=metrics(residual, zero),
                   pt_difference=signed,
                   residual_vs_pt=metrics(residual, difference))
    (output / "metrics.json").write_text(json.dumps(summary, indent=2))
    if summary["repeat"]["max_abs"] != 0:
        raise RuntimeError("Single-thread repeat differs; inspect commands.json.")
    if summary["threads"]["max_abs"] > 1e-5:
        raise RuntimeError("Thread-count discrepancy exceeds floating-point tolerance.")
    if summary["static"]["max_abs"] != 0:
        raise RuntimeError("Scene without dynamic geometry must have a zero residual.")
    if args.strict_identity and summary["unchanged"]["max_abs"] > 1e-6:
        raise RuntimeError("Unchanged scene residual is not zero.")
    if not signed["min"] < 0 < signed["max"]:
        raise RuntimeError("Movement example must contain both signs.")
    print(json.dumps(summary, indent=2))
    print("Regression checks passed. Identity-transform residual and PT agreement are measurements, not correctness guarantees.")


if __name__ == "__main__":
    main()
