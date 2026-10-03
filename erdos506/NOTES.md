# Erdős Problem #506 — the 9-point case

**Problem.** What is the fewest distinct circles determined by *n* points in the plane, not all on one
circle or line? (A collinear triple gives a line, not a circle.)

**Known (erdosproblems.com/506 and its forum, accessed 2026-10-03).**
- Solved for n > 393 (Elliott 1967, corrected by Purdy–Smith): the minimum is
  C(n−1, 2) + 1 − ⌊(n−1)/2⌋.
- Forum post (Aug 2026) gives exact values m(5..8) = 5, 8, 11, 17.
- For **n = 9**, the best known construction gives **25** (8 points on a circle plus one point).
  Open: is m(9) = 25 or smaller? The forum proved any better configuration must have its largest
  circle/line through at most 5 points.

## What this session established

All counts use `circles.py` (exact integer arithmetic) or the float tools in `extend.py`.
These reproduce the forum's values: 8 circles (n=6), 11 (n=7), 17 (n=8), and 18 for the two-squares example.

1. **Cube-based configurations cannot beat 29 (even abstractly).**
   Take the 20-block "cube" structure realised by the forum's 17-circle 8-point set (12 four-point
   circles + 8 triangles), and add a 9th point. Over *all* ways to extend the blocks, choose a
   projection centre, and keep every block to ≤ 5 points, the CP-SAT solver (`cube9.py`) proves
   the minimum number of circles is **29**. So an improvement at n = 9 cannot contain the full
   cube structure on 8 of its points. (Blocks of ≥ 6 points are excluded by the forum's
   proposition, which already shows such configurations have ≥ 25 circles.)

2. **Exhaustive / heuristic searches found nothing below 25.** All searches allow the best
   inversion centre (equivalently, any projection pole on the sphere):
   - all 9-subsets of a 5×5 and 6×6 integer grid: best with largest block ≤ 5 is 32;
     triangular lattice: 29;
   - all 9-subsets of 26 cube/octahedron/cuboctahedron directions on the sphere: 25 (only the
     known 8-on-a-circle type); with largest block ≤ 5: 27;
   - simulated annealing over 74 icosahedral + cubic directions: same, 25 / 27 *(provisional: run
     before a tolerance fix; being re-run)*;
   - ring families (points in layers of regular polygons on the sphere): ≥ 31 with block ≤ 5
     *(provisional, same reason)*;
   - extending the best 6-, 7-, 8-point configurations greedily: ≥ 29;
   - symmetric families (C3, C4, C2, mirror) with incidences forced by solving equations:
     see the run logs; nothing below 25 has been found.

3. **Simple necessary conditions do not prove m(9) = 25.** With exact cover of triples, lines
   through the pole pairwise meeting in ≤ 1 point, Csima–Sawyer at the pole (≥ 5 ordinary lines)
   and at every point (≥ 4 three-point blocks through each point, via inversion), the abstract
   minimum (`lower9.py`) is still ≤ 24 (e.g. B = 33 blocks with 9 lines → 24). Settling n = 9
   needs genuine realizability arguments (like the forum's Fano/radical-axis argument at n = 7).

**Why this is hard.** Up to similarity, 9 planar points have 14 degrees of freedom, but a
≤ 24-circle configuration needs roughly 25+ independent incidences (four-point circles, three-point
lines). So it can only exist if incidence theorems (Miquel-type) make many of those coincidences
automatic — random designs are hopeless; structured ones are the only hope.

Status: **open**. Evidence leans towards m(9) = 25, but nothing here is a proof.
