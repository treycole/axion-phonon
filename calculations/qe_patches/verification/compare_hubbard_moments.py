import re, sys, numpy as np
axis = {5: (-1,-1,-1), 6: (-1,1,1), 7: (1,-1,1), 8: (1,1,-1)}
def load(path):
    blocks, cur, atom = [], None, None
    for ln in open(path).read().splitlines():
        if "HUBBARD OCCUPATIONS" in ln:
            if cur: blocks.append(cur)
            cur = {}
        m = re.search(r"-+ ATOM\s+(\d+) -+", ln)
        if m: atom = int(m.group(1))
        m = re.search(r"Tr\[ns\(\s*(\d+)\)\] \(up, down, total\) =\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", ln)
        if m and cur is not None: cur.setdefault(int(m.group(1)), {})["ud"] = float(m.group(2)) - float(m.group(3))
        m = re.search(r"Atomic magnetic moment mx, my, mz =\s+(\S+)\s+(\S+)\s+(\S+)", ln)
        if m and cur is not None and atom: cur.setdefault(atom, {})["m"] = np.array(list(map(float, m.groups())))
    if cur: blocks.append(cur)
    return [b for b in blocks if all(a in b and "m" in b[a] for a in (5, 6, 7, 8))]
def tilt(a, v):
    n = np.array(axis[a], float); n /= np.linalg.norm(n)
    return np.degrees(np.arccos(min(1, abs(v @ n) / np.linalg.norm(v))))
def report(label, path):
    bl = load(path)
    print(f"\n{label}   ({len(bl)} Hubbard-occupation blocks)")
    print("  block     |m| on Ir1..Ir4                                  spread   tilt from 3-fold axis (deg)")
    for i, b in enumerate(bl):
        mags = [np.linalg.norm(b[a]["m"]) for a in (5,6,7,8)]
        ts = [tilt(a, b[a]["m"]) for a in (5,6,7,8)]
        sp = (max(mags) - min(mags)) / np.mean(mags) * 100
        print(f"  {i+1:4d}   " + " ".join(f"{m:9.5f}" for m in mags) + f"   {sp:6.2f}%   " + " ".join(f"{t:5.2f}" for t in ts))
for lab, p in zip(sys.argv[1::2], sys.argv[2::2]): report(lab, p)
