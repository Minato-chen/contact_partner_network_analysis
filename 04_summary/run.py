"""Collect completed scalar records into one wide comparison table per project/object."""
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import sys
from tempfile import TemporaryDirectory

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from common import inputs
from common import results as files


def record_files(folder):
    return sorted({*folder.rglob("scalars.json"), *folder.rglob("*.scalars.json")})


def table_name(category, object_name):
    category = re.sub(r"^\d+_", "", category)
    name = category + ("_" + object_name if object_name else "")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name):
        raise ValueError(f"Invalid summary table name: {name}")
    return name + ".csv"


def load_record(path, simulation_id):
    payload = files.read_json(path)
    if payload["simulation_id"] != simulation_id:
        raise ValueError("Simulation id differs from its result folder")
    if path.name == "scalars.json":
        if not (path.parent / "COMPLETE").is_file():
            return None, "INCOMPLETE"
        if payload["source"] != files.read_json(path.parent / "metadata.json"):
            return None, "STALE"
    elif payload.get("complete") is not True:
        return None, "INCOMPLETE"
    if not isinstance(payload["records"], list) or not payload["records"]:
        raise ValueError("A completed scalar record must contain at least one metric")
    return payload, "SUCCESS"


def validate_record(record):
    for name in ("metric", "unit", "object"):
        if not isinstance(record[name], str):
            raise ValueError(f"{name} must be text")
    if not record["metric"]:
        raise ValueError("Metric name cannot be empty")
    for name in ("conditions", "definition"):
        if not isinstance(record[name], dict):
            raise ValueError(f"{name} must be a dictionary")
    for name in ("value", "sem"):
        value = record[name]
        if value is not None and (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
        ):
            raise ValueError(f"{name} must be finite or null")
    if record["sem"] is not None and record["sem"] < 0:
        raise ValueError("SEM cannot be negative")


def collect(simulations, results_root):
    results_root = Path(results_root)
    # The union reveals missing metrics in simulations that have not finished yet.
    expected = sorted(
        {
            p.relative_to(results_root / sim["id"])
            for sim in simulations
            for p in record_files(results_root / sim["id"])
        }
    )
    if not expected:
        raise ValueError("No scalar records found; run calculations that export scalars first")
    groups, schema, status, sources, names = {}, {}, [], [], {}
    for simulation in simulations:
        sid = simulation["id"]
        for relative in expected:
            path = results_root / sid / relative
            row = dict(simulation_id=sid, source=str(relative))
            if not path.is_file():
                status.append(row | dict(status="MISSING", message="No result record"))
                continue
            try:
                payload, state = load_record(path, sid)
                if payload is None:
                    status.append(
                        row | dict(status=state, message="Result omitted from tables")
                    )
                    continue
                for record in payload["records"]:
                    validate_record(record)
                    category = relative.parts[0]
                    if len(relative.parts) < 2:
                        raise ValueError(
                            "Scalar records must be inside a project/category folder"
                        )
                    group = (category, record["object"])
                    filename = table_name(*group)
                    if filename in names and names[filename] != group:
                        raise ValueError(f"Different projects map to {filename}")
                    names[filename] = group
                    condition = json.dumps(
                        record["conditions"], sort_keys=True, ensure_ascii=False
                    )
                    key = (record["metric"], condition)
                    definition = (
                        record["unit"],
                        json.dumps(record["definition"], sort_keys=True),
                    )
                    schema_key = (group, key)
                    if schema_key in schema and schema[schema_key] != definition:
                        raise ValueError(f"Unit/definition mismatch for {filename}: {key}")
                    schema[schema_key] = definition
                    columns = groups.setdefault(group, {})
                    values = columns.setdefault(key, {})
                    if sid in values:
                        raise ValueError(f"Duplicate metric for {sid}: {filename}/{key}")
                    values[sid] = (record["value"], record["sem"])
                sources.append(row | dict(metadata=payload["source"]))
                status.append(row | dict(status="SUCCESS", message=""))
            except (KeyError, TypeError, ValueError, OSError) as exc:
                # Malformed or conflicting records must not silently become valid comparisons.
                raise ValueError(f"{path}: {exc}") from exc
    tables = {}
    for group, columns in groups.items():
        rows = [
            dict(
                simulation_id=sim["id"],
                variant=sim.get("variant"),
                temperature_K=sim.get("temperature_K"),
            )
            for sim in simulations
        ]
        labels = {"simulation_id", "variant", "temperature_K"}
        for key, values in columns.items():
            metric, condition = key
            unit = schema[(group, key)][0]
            conditions = json.loads(condition)
            suffix = ", ".join(f"{k}={v}" for k, v in conditions.items())
            label = metric + (f" [{suffix}]" if suffix else "") + (f" ({unit})" if unit else "")
            has_sem = any(error is not None for _, error in values.values())
            added = {label, label + " SEM"} if has_sem else {label}
            if labels & added:
                raise ValueError(f"Ambiguous column label: {label}")
            labels.update(added)
            for row in rows:
                value, error = values.get(row["simulation_id"], (None, None))
                row[label] = value
                if has_sem:
                    row[label + " SEM"] = error
        tables[table_name(*group)] = pd.DataFrame(rows)
    return tables, pd.DataFrame(status), sources


def summarize(simulations, results_root=files.RESULTS, output=HERE / "output"):
    # Validate everything before replacing any previous summary.
    tables, status, sources = collect(simulations, results_root)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    previous = (
        files.read_json(output / "sources.json").get("tables", [])
        if (output / "sources.json").exists()
        else []
    )
    with TemporaryDirectory(prefix="summary_", dir=output.parent) as temporary:
        staged = Path(temporary)
        for name, frame in tables.items():
            frame.to_csv(staged / name, index=False)
        status.to_csv(staged / "status.csv", index=False)
        files.write_json(
            staged / "sources.json",
            dict(
                created_at=datetime.now(timezone.utc).isoformat(),
                tables=list(tables),
                sources=sources,
            ),
        )
        output.mkdir(parents=True, exist_ok=True)
        for path in staged.iterdir():
            path.replace(output / path.name)
        for name in set(previous) - set(tables):
            if Path(name).name != name or not name.endswith(".csv"):
                raise ValueError("Invalid previous summary filename")
            (output / name).unlink(missing_ok=True)
    for name, frame in tables.items():
        print(f"[ok] {output / name}: {len(frame)} simulations", flush=True)
    missing = status[status.status.ne("SUCCESS")]
    if len(missing):
        print(missing.to_string(index=False), flush=True)
    print(f"Source status: {status.status.value_counts().to_dict()}", flush=True)
    return tables


def main():
    cli = inputs.parser(__doc__)
    cli.add_argument("--results-dir", type=Path, default=files.RESULTS)
    cli.add_argument("--output-dir", type=Path, default=HERE / "output")
    args = cli.parse_args()
    try:
        summarize(
            inputs.select_simulations(args.config, args.simulation),
            args.results_dir,
            args.output_dir,
        )
    except (ValueError, OSError) as exc:
        print(f"[failed] Summary not updated: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
