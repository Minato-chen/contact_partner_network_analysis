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


def finish(path, simulation, parameters, inputs, *, scalars=None):
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
    if scalars is not None:
        write_scalars(path, scalars)
    (Path(path) / "COMPLETE").touch()
    print(f"[ok] {path}", flush=True)


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "03_results"


def print_summary(rows, columns):
    """Print batch outcomes only; no persistent state or scheduling."""
    table = [["simulation", *columns], *[[name, *statuses] for name, statuses in rows]]
    widths = [max(len(str(row[i])) for row in table) for i in range(len(columns) + 1)]
    print("\nBatch summary (SUCCESS / FAILED / SKIPPED):", flush=True)
    for row in table:
        print(
            "  ".join(str(value).ljust(width) for value, width in zip(row, widths)), flush=True
        )


def scalar_record(
    metric, value, unit, *, object="", conditions=None, definition=None, sem=None
):
    """One final statistic. Missing values remain null, never zero."""
    import math

    def number(value):
        if value is None:
            return None
        value = float(value)
        if math.isnan(value):
            return None
        if not math.isfinite(value):
            raise ValueError("Scalar summaries cannot contain infinity")
        return value

    error = number(sem)
    if error is not None and error < 0:
        raise ValueError("SEM cannot be negative")
    return dict(
        metric=metric,
        value=number(value),
        unit=unit,
        object=object,
        conditions=conditions or {},
        definition=definition or {},
        sem=error,
    )


def records_from_summary(path, prefix, units, *, object="", conditions=(), definition=None):
    """Select final frame means and their existing block SEM, without recalculation."""
    import pandas as pd

    records = []
    for row in pd.read_csv(path).to_dict("records"):
        for metric, unit in units.items():
            records.append(
                scalar_record(
                    f"{prefix}.{metric}",
                    row[f"{metric}_mean_frames"],
                    unit,
                    object=object,
                    conditions={key: row[key] for key in conditions},
                    definition=definition,
                    sem=row.get(f"{metric}_sem_blocks"),
                )
            )
    return records


def write_scalars(folder, records):
    """Attach a summary to detailed results; the collector also checks COMPLETE."""
    folder = Path(folder)
    metadata = read_json(folder / "metadata.json")
    write_json(
        folder / "scalars.json",
        dict(simulation_id=metadata["simulation_id"], source=metadata, records=records),
    )


def start_scalar(path, simulation, *, parameters=None, inputs=None):
    """For a one-value calculation, invalidate its previous standalone record first."""
    path = Path(path)
    if not path.name.endswith(".scalars.json"):
        raise ValueError("Standalone scalar files must end in .scalars.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(
        path,
        dict(
            simulation_id=simulation["id"],
            complete=False,
            source=dict(parameters=parameters or {}, inputs=inputs or {}),
            records=[],
        ),
    )


def finish_scalar(path, records):
    """Complete the standalone record; no extra folder, CSV or marker file needed."""
    from datetime import datetime, timezone

    value = read_json(path)
    value.update(complete=True, records=records)
    value["source"]["created_at"] = datetime.now(timezone.utc).isoformat()
    write_json(path, value)
