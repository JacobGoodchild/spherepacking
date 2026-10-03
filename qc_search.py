#!/usr/bin/env python3
"""Search for good quasi-cyclic linear codes over small finite fields.

Codes are one-generator, index-l QC codes with systematic generator
    G = [ I_m | C_2 | ... | C_l ],   n = l*m, k = m,
where each C_j is an m x m circulant defined by its first row (a polynomial).

Minimum distance is computed exactly with a Brouwer-Zimmermann style
enumeration over disjoint information sets (every invertible circulant block
gives one), with early exit as soon as a codeword lighter than the target is
seen. The search perturbs the circulant polynomials with simulated annealing.

A code is only reported as a record if --record-d is given (the published
codetables.de lower bound for [n, k]) and the verified distance exceeds it.
"""
import argparse
import datetime as dt
import itertools
import os
import subprocess
import time

import numpy as np

REPO = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- fields
class Field:
    def __init__(self, q):
        self.q = q
        if q == 4:
            # GF(4) = {0, 1, a, a^2} encoded 0,1,2,3; a^2 = a + 1, addition is XOR.
            self.add = np.array([[a ^ b for b in range(4)] for a in range(4)], np.int8)
            log = {1: 0, 2: 1, 3: 2}
            exp = [1, 2, 3]
            mul = np.zeros((4, 4), np.int8)
            for a in range(1, 4):
                for b in range(1, 4):
                    mul[a, b] = exp[(log[a] + log[b]) % 3]
            self.mul = mul
            self.desc = "GF(4) = {0,1,a,a^2} encoded 0,1,2,3 with a^2 = a + 1"
        elif q in (2, 3, 5, 7, 11, 13):
            r = np.arange(q)
            self.add = ((r[:, None] + r[None, :]) % q).astype(np.int8)
            self.mul = ((r[:, None] * r[None, :]) % q).astype(np.int8)
            self.desc = f"GF({q}) = integers mod {q}"
        else:
            raise ValueError("unsupported field")
        self.neg = np.array([int(np.where(self.add[a] == 0)[0][0]) for a in range(q)], np.int8)
        self.inv = np.zeros(q, np.int8)
        for a in range(1, q):
            self.inv[a] = int(np.where(self.mul[a] == 1)[0][0])

    def matmul(self, A, B):
        """A (r x k) times B (k x n) over the field, via tables."""
        out = np.zeros((A.shape[0], B.shape[1]), np.int8)
        for t in range(A.shape[1]):
            out = self.add[out, self.mul[A[:, t][:, None], B[t][None, :]]]
        return out

    def rref(self, M):
        """Row-reduce M; return (R, pivot_columns)."""
        M = M.copy()
        rows, cols = M.shape
        piv = []
        r = 0
        for c in range(cols):
            if r == rows:
                break
            nz = np.nonzero(M[r:, c])[0]
            if len(nz) == 0:
                continue
            p = r + nz[0]
            M[[r, p]] = M[[p, r]]
            M[r] = self.mul[self.inv[M[r, c]], M[r]]
            for i in range(rows):
                if i != r and M[i, c]:
                    M[i] = self.add[M[i], self.mul[self.neg[M[i, c]], M[r]]]
            piv.append(c)
            r += 1
        return M, piv


def circulant(row):
    m = len(row)
    return np.array([np.roll(row, i) for i in range(m)], np.int8)


def generator(F, polys, m):
    blocks = [np.eye(m, dtype=np.int8)] + [circulant(p) for p in polys]
    return np.hstack(blocks)


# ---------------------------------------------------------- min distance
class MinDist:
    """Exact minimum distance via enumeration over disjoint information sets."""

    def __init__(self, F):
        self.F = F

    def info_set_gens(self, G, m, l):
        """Systematic generators w.r.t. each block that is an information set."""
        F = self.F
        k = G.shape[0]
        gens = []
        for b in range(l):
            cols = list(range(b * m, (b + 1) * m))
            order = cols + [c for c in range(G.shape[1]) if c not in cols]
            R, piv = F.rref(G[:, order])
            if piv == list(range(k)):
                inv_order = np.argsort(order)
                gens.append(R[:, inv_order])
        return gens

    def weight_rows(self, G, coeffs_list, support):
        """Codewords from linear combos of rows in `support` with coefficients."""
        F = self.F
        rows = G[list(support)]
        C = np.array(coeffs_list, np.int8)          # (num, w)
        acc = np.zeros((C.shape[0], G.shape[1]), np.int8)
        for t in range(rows.shape[0]):
            acc = F.add[acc, F.mul[C[:, t][:, None], rows[t][None, :]]]
        return np.count_nonzero(acc, axis=1)

    def compute(self, G, m, l, stop_below=0, max_w=None):
        """Return (d, exact). Stops early with (weight, False) if a codeword of
        weight < stop_below is found."""
        F = self.F
        q = F.q
        k = G.shape[0]
        gens = self.info_set_gens(G, m, l)
        if not gens:
            gens = [G]
        g = len(gens)
        best = G.shape[1]
        nz = list(range(1, q))
        for w in range(1, k + 1):
            if max_w is not None and w > max_w:
                return best, False
            # Coefficient vectors with first entry normalized to 1 (scalar multiples).
            coeffs = [(1,) + c for c in itertools.product(nz, repeat=w - 1)]
            for Gs in gens:
                # QC symmetry: shifting a codeword cyclically within all blocks
                # gives another codeword, so the first support index can be 0.
                for rest in itertools.combinations(range(1, k), w - 1):
                    sup = (0,) + rest
                    wt = self.weight_rows(Gs, coeffs, sup).min()
                    if wt < best:
                        best = int(wt)
                        if best < stop_below:
                            return best, False
            # Any codeword not yet seen has weight >= w+1 on each info set.
            if g * (w + 1) >= best:
                return best, True
        return best, True


def brute_force_d(F, G):
    k, n = G.shape
    best = n
    for msg in itertools.product(range(F.q), repeat=k):
        if any(msg):
            c = F.matmul(np.array([msg], np.int8), G)[0]
            best = min(best, int(np.count_nonzero(c)))
    return best


# ---------------------------------------------------------------- search
NO_PUSH = False


def git_push(files, msg):
    if NO_PUSH:
        return False
    subprocess.run(["git", "-C", REPO, "add", *files], capture_output=True)
    subprocess.run(["git", "-C", REPO, "commit", "-q", "-m", msg +
                    "\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n"
                    "Claude-Session: https://claude.ai/code/session_01VFkLXP7kJjth38Z8secrus"],
                   capture_output=True)
    for delay in (0, 2, 4, 8, 16):
        time.sleep(delay)
        if subprocess.run(["git", "-C", REPO, "push", "-q", "-u", "origin", "main"],
                          capture_output=True).returncode == 0:
            return True
    return False


def write_code(path, F, m, l, polys, d, exact, note):
    G = generator(F, polys, m)
    with open(path, "w") as f:
        f.write(f"# {note}\n# field: {F.desc}\n")
        f.write(f"# [n, k, d] = [{m * l}, {m}, {d}]  (distance {'exact' if exact else 'upper bound'})\n")
        f.write(f"# structure: G = [I_{m} | C_2 .. C_{l}], C_j circulant with first row:\n")
        for j, p in enumerate(polys, 2):
            f.write(f"# C_{j}: {' '.join(map(str, p))}\n")
        for row in G:
            f.write(" ".join(map(str, row)) + "\n")


def search(args):
    F = Field(args.q)
    md = MinDist(F)
    rng = np.random.default_rng(args.seed)
    m, l = args.k, args.n // args.k
    assert m * l == args.n, "n must be a multiple of k for this QC family"
    target = args.target_d
    progress = os.path.join(REPO, "progress.txt")
    end = time.time() + args.minutes * 60
    last_log = 0
    checked = 0
    best_d, best_polys = 0, None

    def score(polys):
        # Cheap screen: weight of low-weight messages bounds d from above.
        d, exact = md.compute(generator(F, polys, m), m, l, stop_below=0, max_w=2)
        return d

    polys = [rng.integers(0, F.q, m).astype(np.int8) for _ in range(l - 1)]
    cur = score(polys)
    temp = 1.0
    try:
        while time.time() < end:
            cand = [p.copy() for p in polys]
            j = rng.integers(len(cand))
            for _ in range(rng.integers(1, 3)):
                cand[j][rng.integers(m)] = rng.integers(F.q)
            s = score(cand)
            checked += 1
            if s >= cur or rng.random() < np.exp((s - cur) / max(temp, 1e-3)):
                polys, cur = cand, s
            temp = max(0.05, temp * 0.9995)
            if checked % 2000 == 0:   # thermal jolt / restart
                polys = [rng.integers(0, F.q, m).astype(np.int8) for _ in range(l - 1)]
                cur, temp = score(polys), 1.0
            if s >= max(best_d, 1):
                d, exact = md.compute(generator(F, cand, m), m, l, stop_below=best_d + 1)
                if exact and d > best_d:
                    best_d, best_polys = d, [p.copy() for p in cand]
                    write_code(os.path.join(REPO, "best_code.txt"), F, m, l, best_polys, d, True,
                               "best QC code found so far")
                    if args.record_d is not None and d > args.record_d:
                        write_code(os.path.join(REPO, "world_record_code.txt"), F, m, l,
                                   best_polys, d, True,
                                   f"beats stated lower bound d={args.record_d}")
                        log(progress, f"{now()}, [{args.n},{args.k}]_GF({args.q}), {checked}, {best_d}  # RECORD")
                        git_push(["progress.txt", "best_code.txt", "world_record_code.txt"],
                                 f"QC search: [{args.n},{args.k},{d}] over GF({args.q})")
                        return
                    if target is not None and d >= target:
                        break
            if time.time() - last_log >= args.log_every:
                last_log = time.time()
                log(progress, f"{now()}, [{args.n},{args.k}]_GF({args.q}), {checked}, {best_d}")
                git_push(["progress.txt", "best_code.txt"], f"QC search progress: best d={best_d}")
    except KeyboardInterrupt:
        pass
    log(progress, f"{now()}, [{args.n},{args.k}]_GF({args.q}), {checked}, {best_d}  # STOP")
    git_push(["progress.txt", "best_code.txt"], f"QC search stopped: best d={best_d}")


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def log(path, line):
    with open(path, "a") as f:
        f.write(line + "\n")
    print(line, flush=True)


def selftest():
    rng = np.random.default_rng(0)
    for q in (2, 3, 4, 5, 7):
        F = Field(q)
        md = MinDist(F)
        for _ in range(15):
            m = int(rng.integers(2, 5 if q > 3 else 7))
            l = int(rng.integers(2, 4))
            polys = [rng.integers(0, q, m).astype(np.int8) for _ in range(l - 1)]
            G = generator(F, polys, m)
            d, exact = md.compute(G, m, l)
            bf = brute_force_d(F, G)
            assert exact and d == bf, (q, m, l, polys, d, bf)
    print("selftest passed: info-set min distance matches brute force over GF(2,3,4,5,7)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--q", type=int, default=4)
    ap.add_argument("--n", type=int)
    ap.add_argument("--k", type=int)
    ap.add_argument("--record-d", type=int, help="published lower bound to beat (from codetables.de)")
    ap.add_argument("--target-d", type=int)
    ap.add_argument("--minutes", type=float, default=90)
    ap.add_argument("--log-every", type=float, default=300)
    ap.add_argument("--seed", type=int)
    ap.add_argument("--no-push", action="store_true", help="for test runs: never commit or push")
    a = ap.parse_args()
    NO_PUSH = a.no_push
    selftest() if a.selftest else search(a)
