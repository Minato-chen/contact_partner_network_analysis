"""Read the raw simulation registry and load stage entry modules."""
import argparse
import importlib.util
import os
from pathlib import Path
import re
import sys
import yaml

ROOT = Path(__file__).resolve().parent


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_config(path):
    path = Path(path).resolve()
    document = yaml.safe_load(path.read_text())
    if not isinstance(document, dict) or not isinstance(document.get('simulations'), list):
        raise ValueError('Configuration requires a simulations list')
    raw_root = os.path.expandvars(str(document.get('data_root', '.')))
    if '$' in raw_root:
        raise ValueError('Unresolved environment variable in data_root')
    root = Path(raw_root).expanduser()
    if not root.is_absolute():
        root = path.parent / root
    simulations = []
    seen = set()
    for record in document['simulations']:
        item = dict(record)
        name = item['id']
        if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name):
            raise ValueError(f'Invalid simulation id: {name!r}')
        if name == 'figures':
            raise ValueError('Simulation id figures is reserved for notebook exports')
        if name in seen:
            raise ValueError(f'Duplicate simulation id: {name}')
        seen.add(name)
        derived = set(item) & {'analysis', 'centering', 'enabled'}
        if derived:
            raise ValueError(f'{name}: {sorted(derived)} belongs in stage entry points, not simulations.yaml')
        for key in ('topology', 'trajectory'):
            raw = os.path.expandvars(str(item[key]))
            if '$' in raw:
                raise ValueError(f'{name}: unresolved environment variable in {key}')
            value = Path(raw).expanduser()
            item[key] = str((value if value.is_absolute() else root / value).resolve())
        simulations.append(item)
    if not simulations:
        raise ValueError('No simulations')
    return simulations


def parser(description):
    result = argparse.ArgumentParser(description=description)
    result.add_argument('--config', type=Path, default=ROOT / 'simulations.yaml')
    result.add_argument('--simulation', action='append', help='Simulation id; repeat to select several')
    return result


def select_simulations(config, selected=None):
    items = load_config(config)
    if selected is not None:
        if not selected:
            raise ValueError('No simulations selected; use None to select all')
        missing = set(selected) - {s['id'] for s in items}
        if missing:
            raise ValueError(f'Unknown simulations: {sorted(missing)}')
        items = [s for s in items if s['id'] in selected]
    return items


def print_summary(rows, columns):
    """Print batch outcomes only; no persistent state or scheduling."""
    table = [['simulation', *columns], *[[name, *statuses] for name, statuses in rows]]
    widths = [max(len(str(row[i])) for row in table) for i in range(len(columns) + 1)]
    print('\nBatch summary (SUCCESS / FAILED / SKIPPED):', flush=True)
    for row in table:
        print('  '.join(str(value).ljust(width) for value, width in zip(row, widths)), flush=True)
