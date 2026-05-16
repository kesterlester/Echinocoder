#!/usr/bin/env python3
"""
plot_loc_since.py — plot net lines of code added since a base commit.

Usage:
    python scripts/plot_loc_since.py <base-commit>

Example:
    python scripts/plot_loc_since.py 195e7751c2e488b20

For each commit since <base-commit>, computes:
    insertions - deletions  (vs the base commit)
and plots the result against date.

Output is written to /tmp/loc_since_<short-sha>.png and opened automatically.
"""
import subprocess
import sys
import datetime
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
except ImportError:
    sys.exit("matplotlib is required — install it with: pip install matplotlib")


def run(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, text=True).strip()


def main():
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <base-commit>")

    base = sys.argv[1]

    # Verify the commit exists
    try:
        full_sha = run(["git", "rev-parse", base])
    except subprocess.CalledProcessError:
        sys.exit(f"Unknown commit: {base!r}")

    short_sha = full_sha[:12]

    # Collect commits since base
    log_lines = run(
        ["git", "log", "--format=%H %ai", "--reverse", f"{base}..HEAD"]
    ).splitlines()

    if not log_lines:
        sys.exit("No commits found since the base commit.")

    print(f"Base commit : {base} ({short_sha})")
    print(f"Commits     : {len(log_lines)}")
    print("Computing net lines per commit...", flush=True)

    dates: list[datetime.date] = []
    net_lines: list[int] = []

    for i, entry in enumerate(log_lines, 1):
        sha, date_str = entry.split()[0], entry.split()[1]
        stat = run(["git", "diff", "--shortstat", base, sha])
        ins = dels = 0
        for token in stat.split(","):
            token = token.strip()
            if "insertion" in token:
                ins = int(token.split()[0])
            elif "deletion" in token:
                dels = int(token.split()[0])
        dates.append(datetime.datetime.strptime(date_str, "%Y-%m-%d").date())
        net_lines.append(ins - dels)
        print(f"\r  {i}/{len(log_lines)}", end="", flush=True)

    print()  # newline after progress

    # Plot
    date_dts = [datetime.datetime(d.year, d.month, d.day) for d in dates]

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(date_dts, net_lines, marker=".", markersize=4,
            linewidth=1.2, color="#2563eb")
    ax.fill_between(date_dts, net_lines, alpha=0.12, color="#2563eb")

    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
    fig.autofmt_xdate(rotation=45)

    ax.set_xlabel("Date")
    ax.set_ylabel(f"Net lines added since {short_sha}")
    ax.set_title(
        f"Lines of code added since commit {short_sha}\n"
        f"(insertions − deletions vs base)   "
        f"{len(log_lines)} commits  "
        f"{dates[0]} → {dates[-1]}"
    )
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f"{int(x):,}")
    )
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.set_ylim(bottom=0)

    plt.tight_layout()

    out_path = Path(f"/tmp/loc_since_{short_sha}.png")
    plt.savefig(out_path, dpi=150)

    print(f"Output      : {out_path}")
    print(f"Peak        : {max(net_lines):,} lines")
    print(f"Final       : {net_lines[-1]:,} lines")

    subprocess.Popen(["open", str(out_path)])


if __name__ == "__main__":
    main()
