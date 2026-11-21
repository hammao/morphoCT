"""Command line: ``morphoct run config.yaml``."""

from __future__ import annotations

import argparse
import json
import sys


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="morphoct", description="microCT analysis of rocks, cores and foraminifera")
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run the pipeline described by a YAML config")
    r.add_argument("config")
    sub.add_parser("version", help="print versions")
    a = p.parse_args(argv)
    if a.cmd == "version":
        from morphoct.pipeline import _versions

        print(json.dumps(_versions()))
        return 0
    from morphoct.pipeline import run

    res = run(a.config)
    out = res["config"].get("output", "")
    print(f"done in {res['total_seconds']} s: porosity {res['morphometrics']['porosity']:.4f}; report in {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
