# Usage: python3 plot_rmm.py rmm_bench.csv --save-dir=plots --logx --smooth 3
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt

OPS_ORDER = [
    "rank1",
    "rank0",
    "select1",
    "select0",
    "rank10",
    "select10",
    "excess",
    "fwdsearch",
    "bwdsearch",
    "rmq_pos",
    "rmq_val",
    "rMq_pos",
    "rMq_val",
    "mincount",
    "minselect",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--save-dir", default="plots")
    ap.add_argument("--show", action="store_true")
    ap.add_argument("--logx", action="store_true")
    ap.add_argument("--smooth", type=int, default=0)
    args = ap.parse_args()

    os.makedirs(args.save_dir, exist_ok=True)
    df = pd.read_csv(args.csv)
    df = df.dropna(subset=["ns_per_op"])

    for op in OPS_ORDER:
        d = df[df["op"] == op].copy()
        if d.empty:
            continue

        d = d.groupby("N", as_index=False)["ns_per_op"].median().sort_values("N")

        if args.smooth and args.smooth > 1:
            d["ns_smooth"] = (
                d["ns_per_op"]
                .rolling(window=args.smooth, center=True, min_periods=1)
                .median()
            )
            yplot = d["ns_smooth"]
        else:
            yplot = d["ns_per_op"]

        plt.figure()
        plt.scatter(d["N"], d["ns_per_op"], s=8, alpha=0.3, linewidths=0)
        plt.plot(d["N"], yplot, linewidth=1.5)

        if args.logx:
            plt.xscale("log", base=2)
        plt.xlabel("Размер последовательности N (бит)")
        plt.ylabel("Время на операцию, нс")
        plt.title(f"RmMTree: {op}")
        plt.grid(True, which="both", linestyle="--", alpha=0.4)
        out = os.path.join(args.save_dir, f"{op}.png")
        plt.savefig(out, bbox_inches="tight", dpi=160)
        if not args.show:
            plt.close()
        print(f"[saved] {out}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
