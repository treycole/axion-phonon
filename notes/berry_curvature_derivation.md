# Wannier-interpolated non-Abelian Berry curvature 

---

## 1. Assumptions

An insulator with $M$ occupied bands, gapped from the rest at every $\mathbf k$.
The object of interest is the **gauge-covariant non-Abelian Berry curvature of the
occupied subspace**, an $M\times M$ Hermitian matrix for each pair of directions.

Directions $\mu,\nu$ run over $(k_x,k_y,k_z,\beta)$, with $\beta$ a phonon
amplitude; $i,j,l$ denote the three $k$ directions alone. Everything through §6
is generic in $\mu$; §8 supplies what is special about $\beta$.

Three standing assumptions, used where flagged:

* **(A1)** The occupied bands are separated by a gap from all others at every
  $\mathbf k$. Needed for the energy denominators in §7.
* **(A2)** The $J$ Wannier functions are orthonormal. Used constantly — it is what
  makes the projectors in §3 clean.
* **(A3)** The nominal orbital positions $\boldsymbol\tau_s$ do not depend on
  $\beta$. Needed only in §8.

**Notation.** $s,t,r$ index the $J$ Wannier functions; $v,v'$ the $M$ occupied
(valence) bands; $c$ the $L=J-M$ unoccupied (conduction) bands *inside* the
Wannier space. Blackboard bold ($\mathbb{H}$, $\mathbb{A}$, $\mathbb{\Omega}$)
means a $J\times J$ matrix in the Wannier basis. A superscript block label names
which subspaces an object's two indices live in: $cv$ = rows conduction, columns
valence; $vv$ = valence–valence.

---

## 2. The master identity

Everything below rests on one fact, so it is worth deriving carefully and in full
generality — it will be applied **twice**, to two different subspaces.

**Claim.** Let $\{|w_a\rangle\}_{a=1}^{n}$ be any smooth, orthonormal family
spanning a subspace with projector $\hat P=\sum_a|w_a\rangle\langle w_a|$, and let
$\hat Q=\hat 1-\hat P$. Define the connection and its covariant curvature

$$
\mathcal{A}_{\mu,ab}=i\langle w_a|\partial_\mu w_b\rangle,
\qquad
\mathcal{B}_{\mu\nu}=\partial_\mu\mathcal{A}_\nu-\partial_\nu\mathcal{A}_\mu-i\big[\mathcal{A}_\mu,\mathcal{A}_\nu\big].
\tag{1}
$$

Then

$$
\boxed{\;
\mathcal{B}_{\mu\nu,ab}
=i\langle\partial_\mu w_a|\hat Q|\partial_\nu w_b\rangle
-i\langle\partial_\nu w_a|\hat Q|\partial_\mu w_b\rangle
\;}
\tag{2}
$$

**Derivation.** Differentiate the connection:

$$
\partial_\mu\mathcal{A}_{\nu,ab}
=i\langle\partial_\mu w_a|\partial_\nu w_b\rangle+i\langle w_a|\partial_\mu\partial_\nu w_b\rangle .
$$

The second term is symmetric under $\mu\leftrightarrow\nu$, so it cancels in the
curl:

$$
\partial_\mu\mathcal{A}_{\nu,ab}-\partial_\nu\mathcal{A}_{\mu,ab}
=i\langle\partial_\mu w_a|\partial_\nu w_b\rangle-i\langle\partial_\nu w_a|\partial_\mu w_b\rangle .
\tag{3}
$$

For the commutator, differentiate orthonormality (A2),
$\partial_\mu\langle w_a|w_c\rangle=0$, giving

$$
\langle w_a|\partial_\mu w_c\rangle=-\langle\partial_\mu w_a|w_c\rangle .
\tag{4}
$$

Hence

$$
\sum_c\mathcal{A}_{\mu,ac}\mathcal{A}_{\nu,cb}
=\sum_c\big(i\langle w_a|\partial_\mu w_c\rangle\big)\big(i\langle w_c|\partial_\nu w_b\rangle\big)
\overset{(4)}{=}\sum_c\langle\partial_\mu w_a|w_c\rangle\langle w_c|\partial_\nu w_b\rangle
=\langle\partial_\mu w_a|\hat P|\partial_\nu w_b\rangle,
$$

the two factors of $i$ combining with the sign from (4) to give $+1$. Therefore

$$
-i\big[\mathcal{A}_\mu,\mathcal{A}_\nu\big]_{ab}
=-i\langle\partial_\mu w_a|\hat P|\partial_\nu w_b\rangle
+i\langle\partial_\nu w_a|\hat P|\partial_\mu w_b\rangle .
\tag{5}
$$

Adding (3) and (5) replaces $\hat 1$ by $\hat 1-\hat P=\hat Q$ in both terms,
which is (2). $\;\blacksquare$

**What (2) says.** Only the parts of $|\partial_\mu w\rangle$ that *leave* the
subspace carry curvature. Motion inside the subspace merely re-mixes the basis and
is removed by the commutator — which is exactly why $\mathcal{B}$ is covariant:
under a change of basis $|w_a\rangle\to|w_b\rangle G_{ba}$ with $G$ unitary, the
subspace and $\hat Q$ are untouched, so (2) gives
$\mathcal{B}_{\mu\nu}\to G^\dagger\mathcal{B}_{\mu\nu}G$. Traces and eigenvalues of
$\mathcal{B}$ are therefore basis-independent; individual matrix elements are not.

---

## 3. The Wannier basis and its two complements

### 3.1 Interpolation

Given $J$ Wannier functions $|\mathbf Rs\rangle$ per cell, the Bloch-like
conjugates are

$$
|u_s^W(\mathbf k)\rangle=\sum_{\mathbf R}e^{-i\mathbf k\cdot(\mathbf r-\mathbf R-\boldsymbol\tau_s)}|\mathbf Rs\rangle .
\tag{6}
$$

The $\boldsymbol\tau_s$ in the phase is **convention I**. It is a free *nominal*
position — atomic site or Wannier center — provided one choice is used
consistently. (Setting all $\boldsymbol\tau_s=0$ is convention II.)

For an operator $\hat O$ commuting with lattice translations, insert (6) twice and
translate the summand by $-\mathbf R'$:

$$
\langle u^W_s|\hat O_{\mathbf k}|u^W_t\rangle
=\frac1N\sum_{\mathbf R',\mathbf R}e^{i\mathbf k\cdot(\mathbf R-\mathbf R'+\boldsymbol\tau_t-\boldsymbol\tau_s)}\langle\mathbf R's|\hat O|\mathbf Rt\rangle
=\sum_{\mathbf R}e^{i\mathbf k\cdot\mathbf b_{st}(\mathbf R)}\,O_{st}(\mathbf R),
\tag{7}
$$

$$
\mathbf b_{st}(\mathbf R)=\mathbf R+\boldsymbol\tau_t-\boldsymbol\tau_s,
\qquad
O_{st}(\mathbf R)=\langle\mathbf 0s|\hat O|\mathbf Rt\rangle .
$$

The $1/N$ cancels against the $N$ equivalent choices of $\mathbf R'$. Applying (7)
to $\hat H$ builds $\mathbb{H}_{st}(\mathbf k)$ from `seedname_hr.dat`.

### 3.2 The Wannier-gauge connection

The position operator is **not** lattice periodic, so (7) does not apply directly
to it. Differentiate (6) instead:

$$
i\partial_{k_\mu}|u^W_t\rangle
=\sum_{\mathbf R}(r_\mu-R_\mu-\tau_{t\mu})\,e^{-i\mathbf k\cdot(\mathbf r-\mathbf R-\boldsymbol\tau_t)}|\mathbf Rt\rangle .
$$

Taking the inner product with (6) and translating as in (7),

$$
\mathbb{A}_{\mu,st}(\mathbf k)\equiv\langle u^W_s|i\partial_\mu u^W_t\rangle
=\sum_{\mathbf R}e^{i\mathbf k\cdot\mathbf b_{st}(\mathbf R)}\,
\langle\mathbf 0s|\,\hat r_\mu-R_\mu-\tau_{t\mu}\,|\mathbf Rt\rangle .
$$

The bracket looks asymmetric in $s\leftrightarrow t$, but by orthonormality (A2)
the difference from the symmetric-looking alternative,

$$
\langle\mathbf 0s|(\hat r_\mu-R_\mu-\tau_{t\mu})|\mathbf Rt\rangle
-\langle\mathbf 0s|(\hat r_\mu-\tau_{s\mu})|\mathbf Rt\rangle
=-b_{st,\mu}(\mathbf R)\,\delta_{\mathbf R\mathbf 0}\delta_{st},
$$

is supported only on $\mathbf R=\mathbf 0$, $s=t$, where
$\mathbf b_{ss}(\mathbf 0)=0$. The two forms are therefore identical, and

$$
\boxed{\;
\mathbb{A}_{\mu,st}(\mathbf k)=\sum_{\mathbf R}e^{i\mathbf k\cdot\mathbf b_{st}(\mathbf R)}X_{\mu,st}(\mathbf R),
\qquad
X_{\mu,st}(\mathbf R)=\langle\mathbf 0s|\,\hat r_\mu-\tau_{s\mu}\,|\mathbf Rt\rangle
\;}
\tag{8}
$$

$X_\mu(\mathbf R)$ is read from `seedname_tb.dat`; $\mathbb{A}_\mu$ is Hermitian
because $X_\mu(\mathbf R)=X^\dagger_\mu(-\mathbf R)$. Its plain curl,

$$
\mathbb{\Omega}_{\mu\nu}\equiv\partial_\mu\mathbb{A}_\nu-\partial_\nu\mathbb{A}_\mu
=i\sum_{\mathbf R}e^{i\mathbf k\cdot\mathbf b_{st}(\mathbf R)}
\big[b_{\mu,st}X_{\nu,st}-b_{\nu,st}X_{\mu,st}\big](\mathbf R),
\tag{9}
$$

follows by differentiating (8). It is **not** gauge covariant on its own.

If $X_{\mu}(\mathbf R)=0$ — each orbital a point at its nominal center — then
$\mathbb{A}_\mu=0$ and every external and cross term below vanishes identically.

### 3.3 Hamiltonian gauge and the three-way split

Diagonalizing $\mathbb{H}$ gives $U_{sv}$ ($J\times M$, occupied columns) and
$\bar U_{sc}$ ($J\times L$), with
$U^\dagger U=\mathbb 1_M$, $\bar U^\dagger\bar U=\mathbb 1_L$,
$\bar U^\dagger U=0$, and $UU^\dagger+\bar U\bar U^\dagger=\mathbb 1_J$. Write

$$
|u_v\rangle=\sum_s|u^W_s\rangle U_{sv},
\qquad
|u_c\rangle=\sum_s|u^W_s\rangle\bar U_{sc},
\qquad
p=UU^\dagger,\quad q=\bar U\bar U^\dagger=\mathbb 1_J-p .
\tag{10}
$$

We work in this gauge throughout. By the covariance argument at the end of §2 this
costs nothing: any other basis for the occupied subspace gives
$G^\dagger\mathcal{B}G$, and every invariant is unchanged.

Now split Hilbert space three ways:

$$
\hat 1=\hat P+\hat Q^{\rm in}+\hat Q^{\rm out},
\qquad
\hat Q=\hat 1-\hat P=\hat Q^{\rm in}+\hat Q^{\rm out},
\tag{11}
$$

with $\hat P$ the occupied projector, $\hat Q^{\rm in}=\sum_c|u_c\rangle\langle u_c|$
the unoccupied directions *inside* the Wannier space (represented in the Wannier
basis by $q$), and $\hat Q^{\rm out}$ everything the finite Wannier space omits.
The defining property of the last is

$$
\boxed{\;\hat Q^{\rm out}|u^W_s\rangle=0\quad\text{for every }s.\;}
\tag{12}
$$

**This is the hinge of the whole derivation.** $\hat Q^{\rm out}$ can be reached
only by *moving the Wannier basis*, never by re-mixing its coefficients.

### 3.4 The two mechanisms

Differentiating (10) by the product rule splits the state derivative along exactly
those two mechanisms:

$$
|\partial_\mu u_v\rangle
=\underbrace{\sum_s|u^W_s\rangle\,(\partial_\mu U)_{sv}}_{\text{internal: coefficients move}}
+\underbrace{\sum_s|\partial_\mu u^W_s\rangle\,U_{sv}}_{\text{external: basis moves}} .
\tag{13}
$$

---

## 4. The cross-gap connection

Project (13) onto an unoccupied eigenstate. Using
$\langle u^W_t|\partial_\mu u^W_s\rangle=-i\mathbb{A}_{\mu,ts}$ from (8),

$$
\langle u_c|\partial_\mu u_v\rangle
=\big(\bar U^\dagger\mathbb{A}_\mu U\big)_{cv}\cdot(-i)
+\big(\bar U^\dagger\partial_\mu U\big)_{cv} .
$$

Both terms carry a common factor $-i$ once we name them as connections:

$$
\boxed{\;
\langle u_c|\partial_\mu u_v\rangle=-i\,A^{cv}_{\mu,cv},
\qquad
A^{cv}_\mu=A^{{\rm int},cv}_\mu+A^{{\rm ext},cv}_\mu
\;}
\tag{14}
$$

$$
A^{{\rm int},cv}_\mu=\bar U^\dagger\,i\partial_\mu U ,
\qquad
A^{{\rm ext},cv}_\mu=\bar U^\dagger\,\mathbb{A}_\mu\,U .
\tag{15}
$$

Both are $L\times M$. $A^{cv}_\mu$ is the interband Berry connection between the
occupied manifold and the unoccupied states inside the Wannier space, with one
contribution from each mechanism of (13). The internal block knows only
$\mathbb{H}$; the external block knows only the position matrix.

We also need the valence–valence external block and the rotated curl:

$$
A^{{\rm ext},vv}_\mu=U^\dagger\mathbb{A}_\mu U,
\qquad
\Omega^{{\rm ext},vv}_{\mu\nu}=U^\dagger\mathbb{\Omega}_{\mu\nu}U .
\tag{16}
$$

---

## 5. The two halves of the curvature

Apply the master identity (2) to the occupied subspace and insert the split (11):

$$
\mathcal{B}_{\mu\nu}=\mathcal{B}^{\rm in}_{\mu\nu}+\mathcal{B}^{\rm out}_{\mu\nu}.
\tag{17}
$$

### 5.1 The $\hat Q^{\rm in}$ half

$\hat Q^{\rm in}$ is an explicit sum over the $L$ unoccupied states, so (14) applies
directly. With $\langle\partial_\mu u_v|u_c\rangle=\overline{\langle u_c|\partial_\mu u_v\rangle}=+i\,\overline{A^{cv}_{\mu,cv}}$,

$$
\langle\partial_\mu u_v|\hat Q^{\rm in}|\partial_\nu u_{v'}\rangle
=\sum_c\big(+i\overline{A^{cv}_{\mu,cv}}\big)\big(-iA^{cv}_{\nu,cv'}\big)
=\Big[\big(A^{cv}_\mu\big)^\dagger A^{cv}_\nu\Big]_{vv'},
$$

the factors $(+i)(-i)=1$. Hence

$$
\mathcal{B}^{\rm in}_{\mu\nu}=i\big(A^{cv}_\mu\big)^\dagger A^{cv}_\nu-(\mu\!\leftrightarrow\!\nu).
\tag{18}
$$

Expanding $A^{cv}=A^{{\rm int},cv}+A^{{\rm ext},cv}$ produces **all four** pairings:

$$
\big(A^{cv}_\mu\big)^\dagger A^{cv}_\nu
=\underbrace{\big(A^{{\rm int},cv}_\mu\big)^\dagger A^{{\rm int},cv}_\nu}_{\text{int}\times\text{int}}
+\underbrace{\big(A^{{\rm int},cv}_\mu\big)^\dagger A^{{\rm ext},cv}_\nu
+\big(A^{{\rm ext},cv}_\mu\big)^\dagger A^{{\rm int},cv}_\nu}_{\text{cross}}
+\underbrace{\big(A^{{\rm ext},cv}_\mu\big)^\dagger A^{{\rm ext},cv}_\nu}_{\text{ext}\times\text{ext}} .
\tag{19}
$$

The last pairing is spurious. Keep it and watch it die in §5.3.

### 5.2 The $\hat Q^{\rm out}$ half — the master identity, second use

By (12), $\hat Q^{\rm out}$ annihilates the internal term of (13) outright:

$$
\hat Q^{\rm out}\sum_s|u^W_s\rangle(\partial_\mu U)_{sv}=0
\quad\Longrightarrow\quad
\hat Q^{\rm out}|\partial_\mu u_v\rangle=\sum_s\hat Q^{\rm out}|\partial_\mu u^W_s\rangle U_{sv}.
$$

**The $\hat Q^{\rm out}$ contribution is purely external**, and it is sandwiched
entirely between Wannier states:

$$
\mathcal{B}^{\rm out}_{\mu\nu}=U^\dagger\,\mathbb{K}_{\mu\nu}\,U,
\qquad
\mathbb{K}_{\mu\nu,st}=i\langle\partial_\mu u^W_s|\hat Q^{\rm out}|\partial_\nu u^W_t\rangle-(\mu\!\leftrightarrow\!\nu).
\tag{20}
$$

Now observe that $\{|u^W_s\rangle\}$ is itself a smooth orthonormal family spanning
a subspace, whose orthogonal complement is precisely $\hat Q^{\rm out}$. So
identity (2) applies again, with $n=J$, $\mathcal{A}\to\mathbb{A}_\mu$, and
$\hat Q\to\hat Q^{\rm out}$ — read **right to left** this time, turning the matrix
element into the curl-plus-commutator:

$$
\boxed{\;
\mathbb{K}_{\mu\nu}=\mathbb{\Omega}_{\mu\nu}-i\big[\mathbb{A}_\mu,\mathbb{A}_\nu\big]
\;}
\tag{21}
$$

This is the structural payoff: **curvature leaking into states the Wannier space
does not represent is computable without ever constructing those states.** It is
fixed entirely by $\mathbb{A}_\mu$ and its curl — both available from
`seedname_tb.dat`.

> Two consequences worth recording. First, no second-moment matrix elements
> $\langle\mathbf 0s|r_\mu r_\nu|\mathbf Rt\rangle$ appear anywhere: routing through
> (2) never introduces them. (They *are* needed for the quantum metric, the
> symmetric partner of $\mathcal B$, which is why [VS] carry them.) Second, (21) is
> exact — no truncation, no completeness assumption about the Wannier space.

### 5.3 The cancellation

Resolve the identity *inside* the commutator using $\mathbb 1_J=p+q$ from (10),
and use $U^\dagger\mathbb{A}_\mu\bar U=\big(\bar U^\dagger\mathbb{A}_\mu U\big)^\dagger$
(valid because $\mathbb{A}_\mu$ is Hermitian):

$$
U^\dagger\mathbb{A}_\mu\mathbb{A}_\nu U
=U^\dagger\mathbb{A}_\mu(p+q)\mathbb{A}_\nu U
=A^{{\rm ext},vv}_\mu A^{{\rm ext},vv}_\nu
+\big(A^{{\rm ext},cv}_\mu\big)^\dagger A^{{\rm ext},cv}_\nu .
\tag{22}
$$

Substituting (22) into (21) and then (20):

$$
\mathcal{B}^{\rm out}_{\mu\nu}
=\Omega^{{\rm ext},vv}_{\mu\nu}
-i\big[A^{{\rm ext},vv}_\mu,A^{{\rm ext},vv}_\nu\big]
-i\Big[\big(A^{{\rm ext},cv}_\mu\big)^\dagger A^{{\rm ext},cv}_\nu-(\mu\!\leftrightarrow\!\nu)\Big].
\tag{23}
$$

The last term of (23) is **exactly minus** the ext$\times$ext term of (19). They
cancel.

**What the cancellation means.** A purely external transition
occupied $\to$ unoccupied was counted once explicitly through $\hat Q^{\rm in}$,
and once again inside the $-\mathbb{A}\mathbb{A}$ subtraction that defines "outside
the Wannier space." Only one of them is real.

**Why it matters numerically.** The boundary between $\hat Q^{\rm in}$ and
$\hat Q^{\rm out}$ is an artifact of how many Wannier functions were chosen. This
cancellation is what makes the total independent of that choice; if it is broken by
an implementation bug, the answer drifts with $J$ in a way no window or mesh
convergence can fix.

### 5.4 Result

$$
\boxed{\;\mathcal{B}_{\mu\nu}=\mathcal{B}^{\rm int}_{\mu\nu}+\mathcal{B}^{\rm cross}_{\mu\nu}+\mathcal{B}^{\rm ext}_{\mu\nu}\;}
\tag{24}
$$

$$
\mathcal{B}^{\rm int}_{\mu\nu}
=i\big(A^{{\rm int},cv}_\mu\big)^\dagger A^{{\rm int},cv}_\nu-(\mu\!\leftrightarrow\!\nu),
\tag{25}
$$

$$
\mathcal{B}^{\rm cross}_{\mu\nu}
=i\Big[\big(A^{{\rm int},cv}_\mu\big)^\dagger A^{{\rm ext},cv}_\nu
+\big(A^{{\rm ext},cv}_\mu\big)^\dagger A^{{\rm int},cv}_\nu\Big]-(\mu\!\leftrightarrow\!\nu),
\tag{26}
$$

$$
\mathcal{B}^{\rm ext}_{\mu\nu}
=\Omega^{{\rm ext},vv}_{\mu\nu}-i\big[A^{{\rm ext},vv}_\mu,A^{{\rm ext},vv}_\nu\big].
\tag{27}
$$

Each is $M\times M$, Hermitian in the band indices, antisymmetric in
$(\mu,\nu)$, and covariant. Only the **sum** is independent of $J$: as the Wannier
space grows, directions migrate from $\hat Q^{\rm out}$ into $\hat Q^{\rm in}$ and
weight shifts from (27) into (25)–(26). Two limits bracket the behavior — $J=M$
gives $q=0$ and hence pure (27); $X_\mu(\mathbf R)=0$ gives $\mathbb{A}_\mu=0$ and
hence pure (25).

---

## 6. Computational form

Introduce the two objects an implementation actually holds:

$$
d_\mu\equiv-i\,A^{{\rm int},cv}_\mu,
\qquad
a_\mu\equiv A^{{\rm ext},cv}_\mu .
\tag{28}
$$

Then $A^{cv}_\mu=i d_\mu+a_\mu$, and expanding (25)–(26) gives

$$
\mathcal{B}^{\rm int}_{\mu\nu}=i\big[d^\dagger_\mu d_\nu-d^\dagger_\nu d_\mu\big],
\qquad
\mathcal{B}^{\rm cross}_{\mu\nu}=\big[d^\dagger_\mu a_\nu-a^\dagger_\mu d_\nu\big]-(\mu\!\leftrightarrow\!\nu).
\tag{29}
$$

**The cross term carries no factor of $i$** — the $i$ of (26) is cancelled by the
$i$ hidden in $A^{{\rm int},cv}=id$. Writing out the four terms, the second and
third are the adjoints of the first and fourth, so

$$
\mathcal{B}^{\rm cross}_{\mu\nu}
=\Big(d^\dagger_\mu a_\nu+{\rm h.c.}\Big)-\Big(a^\dagger_\mu d_\nu+{\rm h.c.}\Big),
\tag{30}
$$

which makes Hermiticity manifest and is the cheapest form to evaluate.

---

## 7. The two ingredients

### 7.1 The internal block, from perturbation theory

Differentiate $\mathbb{H}U=UE$ and multiply by $U^\dagger$ on the left, using
$U^\dagger\mathbb{H}=EU^\dagger$. With $W\equiv U^\dagger\partial_\mu U$,

$$
U^\dagger(\partial_\mu\mathbb{H})U+EW=WE+\partial_\mu E .
$$

The $(c,v)$ off-diagonal element ($E_c\neq E_v$ by (A1)) gives
$[U^\dagger(\partial_\mu\mathbb{H})U]_{cv}+E_cW_{cv}=W_{cv}E_v$, so

$$
\boxed{\;
d_{\mu,cv}=\frac{\big[\bar U^\dagger(\partial_\mu\mathbb{H})U\big]_{cv}}{E_v-E_c},
\qquad
A^{{\rm int},cv}_{\mu,cv}=i\,d_{\mu,cv}
\;}
\tag{31}
$$

ordinary first-order perturbation theory: a cross-gap velocity matrix element over
an energy denominator that never vanishes.

**Consistency check.** With no external terms, (25) and (31) give

$$
\mathcal{B}^{\rm int}_{\mu\nu,vv'}
=i\sum_c\frac{\langle u_v|\partial_\mu\mathbb{H}|u_c\rangle\langle u_c|\partial_\nu\mathbb{H}|u_{v'}\rangle}
{(E_v-E_c)(E_{v'}-E_c)}-(\mu\!\leftrightarrow\!\nu),
\tag{32}
$$

the textbook non-Abelian Kubo formula. So the construction reduces correctly.

### 7.2 The external blocks

$a_\mu$, $A^{{\rm ext},vv}_\mu$ and $\Omega^{{\rm ext},vv}_{\mu\nu}$ are just
$\bar U^\dagger(\cdot)U$ and $U^\dagger(\cdot)U$ rotations of the Bloch sums (8)
and (9). No further derivation is required.

---

## 8. The phonon direction $\beta$

Sections 2–7 are generic in $\mu$; nothing about them assumed $\mu$ was a momentum.
What changes for $\beta$ is only how $\partial_\beta\mathbb{H}$ and
$\mathbb{A}_\beta$ are built, because $\partial_\beta$ acts on the Wannier
functions themselves rather than on a Bloch phase.

Assumption **(A3)** enters here: $\partial_\beta\boldsymbol\tau_s=0$. The nominal
positions must be held fixed between the two structures. They need not equal the
Wannier centers — and should not, since the centers *do* move with $\beta$ and that
motion is physical, entering through $\partial_\beta X_\mu$ below.

### 8.1 The parametric connection

Because (A3) freezes $\boldsymbol\tau_t$, $\partial_\beta$ passes through the phase
in (6) and hits only the ket. Repeating the translation argument of (7),

$$
\boxed{\;
\mathbb{A}_{\beta,st}(\mathbf k)=\langle u^W_s|i\partial_\beta u^W_t\rangle
=i\sum_{\mathbf R}e^{i\mathbf k\cdot\mathbf b_{st}(\mathbf R)}\Lambda_{st}(\mathbf R),
\qquad
\Lambda_{st}(\mathbf R)=\langle\mathbf 0s|\partial_\beta|\mathbf Rt\rangle .
\;}
\tag{33}
$$

Note the **explicit factor of $i$**, absent from (8). The asymmetry is real: for a
momentum direction, $i\partial_{k_\mu}$ acting on the Bloch phase manufactures the
position operator and the $i$ is consumed; for $\beta$ there is no phase to act on,
so the $i$ of the definition survives.

**Antisymmetry of $\Lambda$.** Differentiating orthonormality (A2) with respect to
$\beta$, $\langle\partial_\beta\mathbf R's|\mathbf Rt\rangle+\langle\mathbf R's|\partial_\beta\mathbf Rt\rangle=0$.
Translating both to the home cell gives
$\overline{\Lambda_{ts}(\mathbf R'-\mathbf R)}+\Lambda_{st}(\mathbf R-\mathbf R')=0$, i.e.

$$
\Lambda(\mathbf R)=-\Lambda^\dagger(-\mathbf R),
\tag{34}
$$

which is exactly what makes $\mathbb{A}_\beta$ in (33) Hermitian.

### 8.2 The mixed curl

Differentiate (8) with respect to $\beta$ — $\mathbf b_{st}(\mathbf R)$ is
$\beta$-independent by (A3) — and (33) with respect to $k_\mu$:

$$
\partial_\beta\mathbb{A}_\mu=\sum_{\mathbf R}e^{i\mathbf k\cdot\mathbf b}\,\partial_\beta X_{\mu}(\mathbf R),
\qquad
\partial_{k_\mu}\mathbb{A}_\beta=i\cdot i\sum_{\mathbf R}b_\mu e^{i\mathbf k\cdot\mathbf b}\Lambda(\mathbf R)
=-\sum_{\mathbf R}b_\mu e^{i\mathbf k\cdot\mathbf b}\Lambda(\mathbf R).
$$

Hence

$$
\boxed{\;
\mathbb{\Omega}_{\beta\mu,st}=\partial_\beta\mathbb{A}_{\mu,st}-\partial_\mu\mathbb{A}_{\beta,st}
=\sum_{\mathbf R}e^{i\mathbf k\cdot\mathbf b_{st}(\mathbf R)}
\Big[\partial_\beta X_{\mu,st}(\mathbf R)+b_{\mu,st}(\mathbf R)\Lambda_{st}(\mathbf R)\Big].
\;}
\tag{35}
$$

The two terms have distinct physical origins: $\partial_\beta X_\mu$ is the change
of the position matrix elements (centers move, spreads change), while
$b_\mu\Lambda$ is the rotation and deformation of the Wannier functions themselves.

### 8.3 The three new ingredients

| Quantity | Definition | Source |
| --- | --- | --- |
| $\partial_\beta H_{st}(\mathbf R)$ | $\partial_\beta\langle\mathbf 0s\vert\hat H\vert\mathbf Rt\rangle$ | finite difference of two `_hr.dat` |
| $\partial_\beta X_{\mu,st}(\mathbf R)$ | $\partial_\beta\langle\mathbf 0s\vert\hat r_\mu-\tau_{s\mu}\vert\mathbf Rt\rangle$ | finite difference of two `_r.dat` |
| $\Lambda_{st}(\mathbf R)$ | $\langle\mathbf 0s\vert\partial_\beta\vert\mathbf Rt\rangle$ | §8.4 |

### 8.4 Obtaining $\Lambda$ from the AMN and gauge matrices

$\Lambda$ is not written by Wannier90, but it is already determined by data the
workflow produces. Expand both Bloch-like Wannier states in eigenstates of their
own structure:

$$
\Lambda_{st}(\mathbf k)\simeq\frac{\langle w_{s\mathbf k}(0)|w_{t\mathbf k}(\beta)\rangle-\delta_{st}}{\beta}
=\frac{1}{\beta}\Big[\sum_{m,n}\overline{U_{ms}(0)}\,\langle u_{m\mathbf k}(0)|u_{n\mathbf k}(\beta)\rangle\,U_{nt}(\beta)-\delta_{st}\Big].
$$

The first factor $\sum_m\overline{U_{ms}(0)}\langle u_{m\mathbf k}(0)|u_{n\mathbf k}(\beta)\rangle
=\langle w_{s\mathbf k}(0)|u_{n\mathbf k}(\beta)\rangle$ **is an AMN**: the projection
of the *distorted* Bloch states onto the *undistorted* Wannier functions used as
trial orbitals. In Wannier90's storage convention
$A_{ns}=\langle u_{n\mathbf k}(\beta)|w_{s\mathbf k}(0)\rangle$,

$$
\boxed{\;
M(\mathbf k)=A^\dagger(\mathbf k)\,U_{\mathbf k}(\beta),
\qquad
\Lambda^{\rm II}(\mathbf k)=\frac{1}{2\beta}\big[M(\mathbf k)-M^\dagger(\mathbf k)\big].
\;}
\tag{36}
$$

**Take the antihermitian part; do not subtract the identity.** Both remove the
divergent $\delta_{st}/\beta$, since $\mathbb 1$ is Hermitian and drops out of the
antihermitian projection automatically. But the projection is strictly more robust:
at $\beta=0$ the same $A$ generates the gauge, so with the projection gauge
$U=A(A^\dagger A)^{-1/2}$ the polar decomposition $A=UP$ gives
$M(0)=A^\dagger U=PU^\dagger U=P=P^\dagger$. Its antihermitian part vanishes
identically, as $\Lambda(0)$ must. Subtracting $\mathbb 1$ instead leaves
$(P-\mathbb 1)/\beta$, which is zero only in the special case that $A$ is already
semi-unitary. This also gives a free zero test: run (36) at $\beta=0$ and confirm
$\Lambda=0$ to machine precision.

**Convention.** $A$ and $U$ come out of Wannier90 in **convention II**, so (36)
yields $\Lambda^{\rm II}(\mathbf k)=\sum_{\mathbf R}e^{i\mathbf k\cdot\mathbf R}\Lambda(\mathbf R)$
and the transform back is the plain
$\Lambda(\mathbf R)=N_k^{-1}\sum_{\mathbf k}e^{-i\mathbf k\cdot\mathbf R}\Lambda^{\rm II}(\mathbf k)$
— **no bond phases, no Wigner–Seitz weights**. The convention-I phases are applied
downstream when (33) is evaluated; applying them twice is the failure mode to watch
for.

### 8.5 The frozen-Wannier-gauge approximation

Setting $\Lambda=0$, hence $\mathbb{A}_\beta=0$, asserts that the Wannier functions
are transported rigidly with $\beta$. Then $\mathbb{\Omega}_{\mu\beta}=-\partial_\beta\mathbb{A}_\mu$,
evaluated as a finite difference between two Wannierizations. This is legitimate
only if both use the same $J$, orbital ordering, $\boldsymbol\tau_s$, and gauge;
otherwise the difference is a gauge discrepancy rather than a derivative. Note that
"$\Lambda\neq0$" and "the finite difference is contaminated by a gauge change" are
two names for the same thing — which is why (36) is worth computing even if only to
confirm it is small.

---

## 9. The axion response

### 9.1 Why the response is easier than the angle

The Chern–Simons axion angle is

$$
\theta=-\frac{1}{4\pi}\int_{\rm BZ}Q_3,
\qquad
Q_3={\rm tr}\Big(\mathcal{A}\wedge d\mathcal{A}-\tfrac{2i}{3}\mathcal{A}\wedge\mathcal{A}\wedge\mathcal{A}\Big),
\tag{37}
$$

with $\mathcal{A}=\mathcal{A}_\mu dx^\mu$ the **full** occupied connection. $Q_3$ is
not gauge invariant, so evaluating (37) requires an explicitly constructed smooth,
periodic gauge. Its $\beta$-derivative does not.

### 9.2 The prefactor

With $\mathcal{F}=d\mathcal{A}-i\mathcal{A}\wedge\mathcal{A}=\tfrac12\mathcal{B}_{\mu\nu}dx^\mu\wedge dx^\nu$,
the standard descent identity is $dQ_3={\rm tr}(\mathcal{F}\wedge\mathcal{F})$. In
components,

$$
{\rm tr}(\mathcal{F}\wedge\mathcal{F})
=\tfrac14\,\epsilon^{\mu\nu\rho\sigma}\,{\rm tr}\big(\mathcal{B}_{\mu\nu}\mathcal{B}_{\rho\sigma}\big)\,d^4x .
$$

The BZ is a closed 3-torus, so integrating $dQ_3$ over it kills the $d_k$ part and
leaves only the $\beta$ leg:

$$
\boxed{\;
\frac{\partial\theta}{\partial\beta}
=\frac{1}{16\pi}\int_{\rm BZ}d^3k\;\epsilon^{\mu\nu\rho\sigma}\,{\rm tr}\big[\mathcal{B}_{\mu\nu}\mathcal{B}_{\rho\sigma}\big]
\;}
\tag{38}
$$

from $\tfrac{1}{4\pi}\cdot\tfrac14$. **Cross-check:** over a closed $\beta$ loop
this integrates to $\oint d\beta\,\partial_\beta\theta=\tfrac{1}{16\pi}\!\int\! d^4k\,\epsilon\,{\rm tr}(\mathcal{B}\mathcal{B})=2\pi C_2$
with the standard second Chern number
$C_2=\tfrac{1}{32\pi^2}\int d^4k\,\epsilon\,{\rm tr}(\mathcal{B}\mathcal{B})$. Consistent.
(The overall *sign* tracks the sign chosen for $\theta$ in (37) and the orientation
of $\epsilon^{\mu\nu\rho\sigma}$; fix it once against a known case.)

### 9.3 Reduction

$\epsilon^{\mu\nu\rho\sigma}$ needs all four axes distinct, so **every surviving
term carries exactly one $\beta$ leg**. The four placements contribute equally:
moving $\beta$ from slot 4 to slot 3 flips both $\epsilon$ and
$\mathcal{B}_{l\beta}\to\mathcal{B}_{\beta l}$; moving it into the first pair is a
cyclic shift by two indices (even) combined with the cyclicity of the trace.
Therefore, with $\epsilon^{xyz\beta}=+1$,

$$
\boxed{\;
\frac{\partial\theta}{\partial\beta}
=\frac{1}{4\pi}\int_{\rm BZ}d^3k\;\epsilon^{ijl}\,{\rm tr}\big[\mathcal{B}_{ij}\,\mathcal{B}_{l\beta}\big]
=\frac{1}{2\pi}\int_{\rm BZ}d^3k\;\sum_l{\rm tr}\big[\widetilde{\mathcal{B}}_l\,\mathcal{B}_{l\beta}\big]
\;}
\tag{39}
$$

with $\widetilde{\mathcal{B}}_l=\tfrac12\epsilon^{lij}\mathcal{B}_{ij}$. In words:
**the BZ integral of the trace of (ordinary Berry curvature) $\times$ (mixed
$k$–phonon curvature)**. Because it is a matrix trace, the full $M\times M$
curvature must be retained until after the product is formed —
${\rm tr}(\mathcal{B}\mathcal{B})\neq{\rm tr}(\mathcal{B}){\rm tr}(\mathcal{B})$.

---

## 10. Provenance, and the state of the literature

### 10.1 What is derived here vs. assumed

Every step above is derived from (A1)–(A3) except two standard results quoted
without proof: the descent identity $dQ_3={\rm tr}(\mathcal F\wedge\mathcal F)$, and
the expression (37) for $\theta$ itself. Both are textbook.

The derivation is **complete** in the sense that nothing else is needed to go from
Wannier90 output to $\partial_\beta\theta$. It is **not** a derivation of *why*
$\theta$ is the magnetoelectric coupling — that is Qi–Hughes–Zhang / Essin–Moore–Vanderbilt
and is taken as given.

### 10.2 Relation to [VS]

The result (24)–(27) is equivalent to [VS] Eq. (61), reached by a different route.
[VS] work in a general smooth gauge $V$, carry the metric–curvature tensor
$F_{\mu\nu}=\langle\partial_\mu u|\hat Q|\partial_\nu u\rangle$ and the second-moment
matrix $\mathbb{W}_{\mu\nu}$ throughout, and take the antisymmetric part at the end.
The route above works in the Hamiltonian gauge from the start (licensed by the
covariance argument in §2) and antisymmetrizes immediately, so $\mathbb{W}$ never
appears. Same answer, fewer objects. [VS]'s route additionally delivers the quantum
metric, which this one discards.

Correspondence: $\bar A^I_\mu=A^{{\rm int},cv}_\mu$, $\bar A^E_\mu=A^{{\rm ext},cv}_\mu$,
$A^E_\mu=A^{{\rm ext},vv}_\mu$, $\Omega^E_{\mu\nu}=\Omega^{{\rm ext},vv}_{\mu\nu}$.

### 10.3 Errors found in [VS]

1. **Eq. (69), the response prefactor.** [VS] write $1/16$; it should be $1/16\pi$.
   This follows from their own Eq. (68) via §9.2 above, and from the second-Chern
   cross-check. Confirmed independently by both routes.
2. **Eq. (73), the finite-difference formula for $\Lambda$.** Missing the identity
   subtraction — $\Lambda\simeq[\langle\mathbf 0s|\mathbf Rt(\beta)\rangle-\delta_{\mathbf 0s,\mathbf Rt}]/\beta$.
   Reported by DV, to be corrected. (The copy of the notes in this repository
   predates that equation; it ends at Eq. (72).)
3. Minor: their Eq. (68) uses $\tfrac12\mathcal A_\mu\mathcal B_{\nu\sigma}+\tfrac{i}{3}\mathcal A^3$
   where the familiar form is $\mathcal A\partial\mathcal A-\tfrac{2i}{3}\mathcal A^3$.
   These are **identical**, not a discrepancy — substitute
   $\mathcal{B}=\partial\mathcal A-\partial\mathcal A-i[\mathcal A,\mathcal A]$ and use
   $\epsilon^{\mu\nu\sigma}{\rm tr}\,\mathcal A_\mu[\mathcal A_\nu,\mathcal A_\sigma]=2\epsilon^{\mu\nu\sigma}{\rm tr}\,\mathcal A_\mu\mathcal A_\nu\mathcal A_\sigma$.

No error was found in the core chain, Eqs. (38)–(61) of [VS].

### 10.4 What is published

| Result | Status |
| --- | --- |
| External/internal split, **band-summed Abelian** curvature (AHC) | Published: Wang, Yates, Souza, Vanderbilt, PRB **74**, 195118 (2006) |
| Wannier-space formalism with the $\mathbb{W}-\mathbb{A}\mathbb{A}$ structure, band-summed | Published: Lopez, Vanderbilt, Thonhauser, Souza, PRB **85**, 014435 (2012) |
| Extension to $r H r$ operators / spatial dispersion | Published: Urru, Souza, Pozo Ocaña, Tsirkin, Vanderbilt, PRB **112**, 045201 (2025) |
| **Matrix-valued non-Abelian** external curvature, Eq. (24)–(27) | **Apparently unpublished.** [VS] compare only *traces* to Lopez 2012 (their footnotes 5 and 6), and DV's own margin note asks whether it agrees with the literature, with IS replying that he needs to check an unpublished preprint. Treat as new until shown otherwise. |
| Parametric $\beta$ terms, $\Lambda$, Eqs. (33)–(35) | [VS] §3 only; unpublished, and the one equation there with an error |
| $\Lambda$ from AMN $\times$ gauge matrix, Eq. (36) | New (J.-M. Lihm); the antihermitian-projection argument is DV's |
| Smooth-gauge construction needed for $\theta$ itself | Published: Cole & Vanderbilt, PRB **113**, 245106 (2026) |
| $\theta$ in real space / hybrid Wannier | Published: Coh, Vanderbilt, Malashevich, Souza, PRB **83**, 085108 (2011); Olsen *et al.*, PRB **95**, 075137 (2017); Bradlyn & Vanderbilt, PRB **101**, 155130 (2020) |
| Phonon-induced $\partial\theta/\partial\beta$ by Wannier interpolation | No published implementation found |

**Reading of this table.** The Abelian, band-summed machinery is twenty years old
and solid. The non-Abelian matrix generalization — which the axion response
*requires*, since ${\rm tr}(\mathcal{B}\mathcal{B})\neq{\rm tr}\mathcal{B}\,{\rm tr}\mathcal{B}$
— rests on [VS] and on §§2–7 above, and is not in the published literature that I
could locate. The parametric extension rests on [VS] §3 plus §8 above, and is
newer still. That is worth knowing when deciding how much independent numerical
validation the mixed planes deserve: the $k$–$k$ blocks can be checked against
`postw90`, but the $k$–$\beta$ blocks currently have no external reference.
