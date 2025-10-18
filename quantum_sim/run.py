from __future__ import annotations

import argparse
import json
import signal
import sys
import time
from pathlib import Path
from typing import Dict

from .config import DEFAULT_COUNTS
from .macros import macro_metrics
from .world import Universe
from .particles import ParticleType


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Quantum-inspired 19-particle simulation")
    p.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    p.add_argument("--steps", type=int, default=0, help="Steps to run; 0 means run forever")
    p.add_argument("--interval", type=int, default=50, help="Metric print interval in steps")
    p.add_argument("--size", type=float, default=None, help="World size (default from constants)")
    p.add_argument("--counts", type=str, default=None,
                   help="JSON mapping of particle type to count to override defaults")
    p.add_argument("--snapshots", type=str, default=None,
                   help="Directory to save JSONL snapshots (optional)")
    return p.parse_args(argv)


def parse_counts(s: str) -> Dict[ParticleType, int]:
    raw = json.loads(s)
    result: Dict[ParticleType, int] = {}
    for k, v in raw.items():
        try:
            pt = ParticleType(k)
        except ValueError:
            # allow keys that match enum values or names
            try:
                pt = ParticleType[k]
            except Exception as e:
                raise SystemExit(f"Unknown particle type key: {k}") from e
        result[pt] = int(v)
    return result


RUNNING = True


def _signal_handler(signum, frame):  # noqa: D401
    global RUNNING
    RUNNING = False


signal.signal(signal.SIGINT, _signal_handler)


def main(argv=None) -> int:
    args = parse_args(argv)
    counts = parse_counts(args.counts) if args.counts else DEFAULT_COUNTS
    u = Universe(seed=args.seed, counts=counts)
    out_dir = None
    if args.snapshots:
        out_dir = Path(args.snapshots)
        out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    step = 0
    target = args.steps if args.steps > 0 else None
    while RUNNING and (target is None or step < target):
        u.step(1)
        step += 1
        if step % args.interval == 0:
            m = macro_metrics(u)
            elapsed = time.time() - t0
            print(f"step={step} elapsed={elapsed:.1f}s metrics=" + json.dumps(m))
        if out_dir and step % args.interval == 0:
            snap = u.snapshot(limit=None)
            with (out_dir / f"snapshot_{step:09d}.jsonl").open("w") as f:
                for item in snap:
                    f.write(json.dumps({
                        "id": item[0],
                        "type": item[1],
                        "x": item[2],
                        "y": item[3],
                        "vx": item[4],
                        "vy": item[5],
                        "color": item[6],
                    }) + "\n")
    print("Simulation finished.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
