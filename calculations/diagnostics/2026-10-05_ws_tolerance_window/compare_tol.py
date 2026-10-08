import numpy as np, re
from collections import defaultdict
S=str(__import__('pathlib').Path(__file__).parent) + "/"
r=np.load(S+"base_full_group_tolerances.npz")
t5,t53=r["tol_1e-05"],r["tol_5e-03"]
rows=[l for l in open(S+"base_full_group_widened.log") if l.startswith("  #")]
t3=np.array([[float(x) for x in re.findall(r"([\d.]+)%",l)] for l in rows])/100
g=np.load("/Users/treycole/Repos/axion-phonon/calculations/diagnostics/2026-09-29_centre_symmetrization/logs/point_groups.npz",allow_pickle=True)
ops,pr=g["base_ops"],g["base_priming"]
def cls(R):
    d=np.linalg.det(R); P=R*d; t=round(np.trace(P))
    name={3:"E",-1:"C2",0:"C3",1:"C4"}[t]
    if t==-1:
        w,v=np.linalg.eig(P); ax=np.real(v[:,np.argmin(abs(w-1))]); name="C2<100>" if np.count_nonzero(abs(ax)>1e-6)==1 else "C2<110>"
    return ("" if d>0 else "I*")+name
grp=defaultdict(list)
for i in range(1,48): grp[(cls(ops[i]),str(pr[i]))].append(i)
print("class / median over ops:  total 1e-5 | 1e-3 | 5e-3      external 1e-5 | 1e-3 | 5e-3")
for k,idx in sorted(grp.items(), key=lambda kv: np.median(t53[kv[1],3])):
    f=lambda t,j: 100*np.median(t[idx,j])
    print(f"{k[0]:10s} {k[1]:9s} n={len(idx)}  {f(t5,3):6.2f} {f(t3,3):6.2f} {f(t53,3):6.2f}    {f(t5,2):6.1f} {f(t3,2):6.1f} {f(t53,2):6.1f}")
for j,p in enumerate(["internal","cross","external","total"]):
    print(f"{p:9s} median over 47 ops: {100*np.median(t5[1:,j]):6.2f} {100*np.median(t3[1:,j]):6.2f} {100*np.median(t53[1:,j]):6.2f}  max: {100*t5[1:,j].max():6.2f} {100*t3[1:,j].max():6.2f} {100*t53[1:,j].max():6.2f}")
ratio=t53[1:,:]/t3[1:,:]
print("5e-3/1e-3 per-op ratio, median by part:", np.round(np.median(ratio,0),3), " ops that got worse (total):", int((ratio[:,3]>1).sum()))
