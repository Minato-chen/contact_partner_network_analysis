"""Small file utilities. Each stage owns its calculation and completion record."""
import json
from pathlib import Path
import shutil


class InputUnavailable(ValueError):
    """Required upstream result is absent, incomplete or incompatible."""


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def reset_output(root, name):
    """Clear only the selected simulation. Never follow an old output symlink."""
    path = Path(root) / name
    if path.is_symlink():
        path.unlink()
    elif path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True)
    return path


def completed(path):
    path = Path(path)
    if not (path / "COMPLETE").is_file():
        raise InputUnavailable(f"{path}: incomplete result; run this calculation first")
    return path


def finish(path, simulation, parameters, inputs):
    from datetime import datetime, timezone

    write_json(
        Path(path) / "metadata.json",
        dict(
            simulation_id=simulation["id"],
            parameters=parameters,
            inputs=inputs,
            created_at=datetime.now(timezone.utc).isoformat(),
        ),
    )
    (Path(path) / "COMPLETE").touch()
    print(f"[ok] {path}", flush=True)
