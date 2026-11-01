"""Shared file operations and completion markers for result outputs."""
from datetime import datetime, timezone
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class InputUnavailable(ValueError):
    """Required upstream result is absent, incomplete or incompatible."""


def completed(path, hint=None):
    path = Path(path)
    if not (path / "COMPLETE").is_file():
        msg = f"{path}: 缺少上游完成标记 (COMPLETE)"
        if hint:
            msg += f"；提示：{hint}"
        else:
            # 自动推断人类可读的修复命令
            try:
                rel = path.relative_to(ROOT / "03_results")
                parts = rel.parts  # (<sim_id>, <stage>, ...)
                if len(parts) >= 2:
                    stage = parts[1]
                    if stage == "01_basic":
                        msg += f"；请先运行: python 02_basic_calculate/run.py --simulation {parts[0]}"
                    elif stage == "02_contacts":
                        msg += f"；请先运行: python 02_calculate/02_contacts/run.py --simulation {parts[0]}"
            except Exception:
                pass
        raise InputUnavailable(msg)
    return path


def finish(path, simulation, parameters, inputs, *, scalars=None):
    from datetime import datetime, timezone

    write_json(
        path / "metadata.json",
        dict(
            simulation_id=simulation["id"],
            parameters=parameters,
            inputs=inputs,
            created_at=datetime.now(timezone.utc).isoformat(),
        ),
    )
    if scalars is not None:
        write_scalars(path, scalars)
    (path / "COMPLETE").touch()
    print(f"[ok] {path}", flush=True)


def reset_output(parent, name):
    target = Path(parent) / name
    target.mkdir(parents=True, exist_ok=True)
    (target / "COMPLETE").unlink(missing_ok=True)
    return target


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, document):
    Path(path).write_text(json.dumps(document, indent=2, ensure_ascii=False) + "\n")


RESULTS = ROOT / "03_results"


def scalar_record(metric, value, unit, *, sem=None, conditions=None, definition=None, object=""):
    return dict(
        metric=metric,
        value=value,
        sem=sem,
        unit=unit,
        object=object,
        conditions={} if conditions is None else conditions,
        definition={} if definition is None else definition,
    )


def write_scalars(path, records):
    write_json(path / "scalars.json", dict(records=records))


def start_scalar(path, simulation, parameters=None, inputs=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(
        path,
        dict(
            complete=False,
            simulation_id=simulation["id"],
            source=dict(
                parameters={} if parameters is None else parameters,
                inputs={} if inputs is None else inputs,
                created_at=datetime.now(timezone.utc).isoformat(),
            ),
            records=[],
        ),
    )


def finish_scalar(path, records):
    path = Path(path)
    document = read_json(path)
    document["records"] = records
    document["complete"] = True
    write_json(path, document)
    print(f"[ok] {path}", flush=True)


def print_summary(rows, headers):
    import pandas as pd
    col_names = ["simulation"] + list(headers)
    table_data = []
    for sid, statuses in rows:
        table_data.append([sid] + list(statuses))
    df = pd.DataFrame(table_data, columns=col_names)
    print("\nBatch summary (SUCCESS / FAILED / SKIPPED):")
    print(df.to_string(index=False), flush=True)
