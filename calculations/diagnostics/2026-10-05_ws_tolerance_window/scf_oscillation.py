import re, sys, numpy as np
t = open(sys.argv[1]).read()
its = re.split(r"\n\s+iteration #\s*(\d+)", t)
rows = []
for k in range(1, len(its) - 1, 2):
    n, b = int(its[k]), its[k + 1]
    g = lambda pat: (lambda m: float(m.group(1)) if m else np.nan)(re.search(pat, b))
    tr = re.search(r"Tr\[ns\(\s*5\)\] \(up, down, total\) =\s+([\d.]+)\s+([\d.]+)", b)
    hm = re.search(r"ATOM    5 -.*?Atomic magnetic moment mx, my, mz =\s+([-\d.]+)", b, re.S)
    site = dict((int(a), (float(c), float(m))) for a, c, m in re.findall(
        r"atom number\s+(\d+) relative position.*?\n\s+charge :\s+([\d.]+).*?\n\s+magnetization :\s+([-\d.]+)", b))
    rows.append(dict(it=n, E=g(r"total energy\s+=\s+([-\d.]+) Ry"), acc=g(r"estimated scf accuracy\s+<\s+([\d.Ee-]+)"),
                     EU=g(r"HUBBARD ENERGY =\s+([-\d.]+)"), nup=float(tr.group(1)) if tr else np.nan,
                     ndn=float(tr.group(2)) if tr else np.nan, mU=float(hm.group(1)) if hm else np.nan,
                     q5=site.get(5, (np.nan,))[0], m5=site.get(5, (np.nan, np.nan))[1],
                     q9=site.get(9, (np.nan,))[0], m9=site.get(9, (np.nan, np.nan))[1],
                     q11=site.get(11, (np.nan,))[0], m11=site.get(11, (np.nan, np.nan))[1],
                     q1=site.get(1, (np.nan,))[0], m1=site.get(1, (np.nan, np.nan))[1],
                     nconv=sum(int(x) for x in re.findall(r"c_bands:\s+(\d+) eigenvalues not converged", b))))
keys = ["E", "acc", "EU", "nup", "ndn", "mU", "q5", "m5", "q9", "m9", "q11", "m11", "q1", "m1", "nconv"]
last = [r for r in rows if r["it"] >= max(r["it"] for r in rows) - 24]
print("iter " + " ".join(f"{k:>11s}" for k in keys))
for r in last:
    print(f"{r['it']:4d} " + " ".join(f"{r[k]:11.6f}" if k not in ("acc", "nconv") else (f"{r[k]:11.2e}" if k == "acc" else f"{r[k]:11d}") for k in keys))
print("\nperiod-2 test over the last 20 iterations: |mean(even) - mean(odd)| vs the spread within each parity")
rows = [r for r in rows if not np.isnan(r["E"])]
L = [r for r in rows if r["it"] >= max(r["it"] for r in rows) - 19]
for k in keys:
    v = np.array([r[k] for r in L], float)
    if np.all(np.isnan(v)): continue
    ev, od = v[0::2], v[1::2]
    print(f"  {k:5s}: even-odd difference {abs(ev.mean() - od.mean()):.3e}   spread within parity {max(ev.std(), od.std()):.3e}")
early = [r for r in rows if r["it"] <= 47]
print("\n'eigenvalues not converged' warnings per iteration: iterations 1-47:", sum(r["nconv"] for r in early), "total;",
      "iterations 48+:", sum(r["nconv"] for r in rows if r["it"] > 47), "total over", sum(1 for r in rows if r["it"] > 47), "iterations")
