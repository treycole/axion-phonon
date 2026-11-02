Using Bilbao Crystallographic Server to diagnose axion active phonons in MBT and study its symmetry properties.
## Magnetic/Crystallographic space and point groups

1. Go under 'Magnetic Symmetry and Applications' to [MAGNDATA](https://cryst.ehu.es/magndata/search.php?atoms_srch=Mn+Bi+Te&op_srch=AND&incomm_type_srch=all&submit=Search)
2. Type 'Mn Bi Te' in "Element search", with 'AND' selected
3. There are 6 results, two of which are MnBi$_2$Te$_4$ 
	1. [Entry 1.308](https://cryst.ehu.es/magndata/index.php?this_label=1.308) 
	2. [Entry 1.309](https://cryst.ehu.es/magndata/index.php?this_label=1.309)

| Quantity                            | MBT                                        |
| ----------------------------------- | ------------------------------------------ |
| Parent crystallographic space group | $R\bar{3}m$ (166)                          |
| Parent crystallographic point group | $\bar{3}m = D_{3d}$                        |
| Magnetic space group                | $R_I\bar{3}c$ (167.108)                    |
| Magnetic point group                | $\bar{3}m1^\prime$                         |
| Magnetic propagation vector         | $(0, 0, 3/2)$ parent hexagonal coordinates |
| Magnetic cell transformation        | $(a, b, 2c; 0, 0, 0)$                      |
| Magnetic order irrep                | $mT2^+$                                    |
| Mn1                                 | $(0,0,0)$ multiplicity 6                   |
| Bi1                                 | $(0,0,0.21235)$ multiplicity 12            |
| Te1                                 | $(0,0,0.06620)$ multiplicity 12            |
| Te2                                 | $(0,0,0.14715)$ multiplicity 12            |

## Wyckoff Positions

We start in the parent crystallographic space group $R \bar{3}m$ (166)

1. Go under 'Space-group symmetry' to [WYCKPOS](https://cryst.ehu.es/cryst/get_wp.html)
2. Type the space group number 166 in the search bar

[List of Wyckoff positions for hexagonal axes](https://cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list) 
![[Pasted image 20261027172508.png]]

3. From MAGNDATA, we have that Mn1 is at (0,0,0) in hexagonal coordinates. This corresponds to Wyckoff site **3a** with _site symmetry_ $\bar{3}m$ 

## MSG symmetry operations

We can get all the symmetry operations of the full MSG

1. Go to Magnetic Symmetry and Application -> MGENPOS
2. Enter MSG number 167.108
3. Obtain [full list of 72 operations](https://cryst.ehu.es/cgi-bin/cryst/programs/nph-magtrgen)

At the top of the page, we see the anti-translation $S=(1^\prime| 0, 0, 1/2)_m$. This is proper and primed so $\chi_\theta(S) = -1$. The phonon must also transform as $\chi_Q(S) = -1$ in order for a linear coupling to be possible. Since the phonon is even under time-reversal, it must be odd under half translation. This picks out the $T = (0,0,3/2)_{p, hex}$ point that folds to magnetic $\Gamma$.
### Unitary subgroup of MSG

To get the unitary subgroup we can
- Remove anti-translation subscript $I$ from MSG giving $R\bar{3}c$ (167)
- Enumerate all operations from MSG that exclude time-reversal
	- Rows 1-36 are unprimed
	- Rows 3-36 repeated 1-12 with the $R$-centering translations

![[Pasted image 20261028003442.png]]
![[Pasted image 20261028003502.png]]
## Character table of CPG

1. Under 'Point-group symmetry' go to 'Point'
2. Select MBT CPG (-3m)
3. We get the [character table](https://cryst.ehu.es/cgi-bin/rep/programs/sam/point.py?sg=162&num=20) of -3m 

![[Pasted image 20261027173032.png]]

4. Organize unitary operations 1-12 above into the character table columns and compare axion characters to get phonon irrep (all unitary by construction). The translation phase appears in the character of the phonon since it must be odd under half translation along c. Required spatial character of the phonon follows from matching effect on $\theta$

| CPG column                 | 1    | 3    | 2$_{1-10}$ | -1    | -3    | m$_{1-10}$ |
| -------------------------- | ---- | ---- | ---------- | ----- | ----- | ---------- |
| MGENPOS Ops                | 1    | 2-3  | 4-6        | 7     | 8-9   | 10-12      |
| Translation along c        | 0    | 0    | 1/2        | 0     | 0     | 1/2        |
| Antiunitary                | No + | No + | No +       | No +  | No +  | No +       |
| Improper                   | No + | No + | No +       | Yes - | Yes - | Yes -      |
| Effect on $\theta$         | +    | +    | +          | -     | -     | -          |
| Translation phase          | +    | +    | -          | +     | +     | -          |
| Required spatial character | +    | +    | -          | -     | -     | +          |

5. Looking at the last row, this picks out the $A_{2u}$  irrep. At the $T$ point this is $T2-$.

