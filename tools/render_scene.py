import argparse
import json
from pathlib import Path
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description="Render the residual for a supplied scene.")
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--scene", choices=["paper-cornell-box", "paper-staircase", "cornell-box"],
                        required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--spp", type=int)
    parser.add_argument("--resolution", type=int)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    for name in ("spp", "resolution", "threads"):
        value = getattr(args, name)
        if value is not None and value <= 0:
            parser.error(name.replace("_", "-") + " must be positive")
    root = Path(__file__).resolve().parents[1]
    scene_name = "scene.xml" if args.scene == "cornell-box" else "residual.xml"
    scene = root / "scenes" / args.scene / scene_name
    suffix = ".exe" if (args.bin_dir / "glacier.exe").exists() else ""
    renderer = (args.bin_dir / ("glacier" + suffix)).resolve()
    if not renderer.is_file():
        parser.error("Executable not found: " + str(renderer))
    output = (args.output or root / "output" / args.scene).resolve()
    output.mkdir(parents=True, exist_ok=True)
    records = []

    def run(command):
        command = [str(part) for part in command]
        print("Running:", subprocess.list2cmdline(command), flush=True)
        start = time.perf_counter()
        result = subprocess.run(command, capture_output=True, text=True)
        records.append(dict(command=command, seconds=time.perf_counter() - start,
                            returncode=result.returncode, stdout=result.stdout,
                            stderr=result.stderr))
        (output / "commands.json").write_text(json.dumps(records, indent=2))
        print(result.stdout, end="")
        print(result.stderr, end="")
        if result.returncode:
            raise SystemExit(result.returncode)

    common = [renderer, scene, "-t", args.threads]
    if args.resolution is not None:
        common.extend(["--resolution", args.resolution])
    residual = output / "residual.exr"
    residual_command = common + ["--integrator", "residual", "-o", residual]
    if args.spp is not None:
        residual_command.extend(["--spp", args.spp])
    run(residual_command)
    print("Outputs and command log:", output)


if __name__ == "__main__":
    main()
