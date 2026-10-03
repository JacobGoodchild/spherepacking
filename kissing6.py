#!/usr/bin/env python3
"""Search for 73 unit vectors in R^6 with all pairwise dot products <= 0.5.

Force-directed optimization with simulated annealing and periodic thermal
jolts, run in several independent worker processes. The coordinator logs to
progress.txt and commits/pushes it every LOG_INTERVAL seconds. If a valid
configuration is found it is written to breakthrough_coordinates.txt, pushed,
and the run stops.
"""
import datetime as dt
import multiprocessing as mp
import os
import signal
import subprocess
import sys
import time

import numpy as np

N, D = 73, 6
LIMIT = 0.5
TARGET = 0.5 - 2e-4          # optimize slightly below 0.5 so a hit passes strictly
LOG_INTERVAL = 7 * 60
END_TIME = dt.datetime(2026, 10, 3, 14, 0, 0, tzinfo=dt.timezone.utc)
REPO = os.path.dirname(os.path.abspath(__file__))
PROGRESS = os.path.join(REPO, "progress.txt")
BREAKTHROUGH = os.path.join(REPO, "breakthrough_coordinates.txt")
BEST_FILE = os.path.join(REPO, "best_configuration.txt")
BRANCH = "main"


def normalize(X):
    n = np.linalg.norm(X, axis=1, keepdims=True)
    n[n < 1e-12] = 1.0
    return X / n


def max_dot(X):
    G = X @ X.T
    np.fill_diagonal(G, -np.inf)
    return float(G.max())


def random_config(rng):
    return normalize(rng.standard_normal((N, D)))


def worker(wid, seed, best_val, iters, found, queue):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    rng = np.random.default_rng(seed)
    my_best = np.inf
    eye = np.eye(N, dtype=bool)
    local_iters = 0
    while not found.value:
        X = random_config(rng)
        lr = 0.05
        run_best = np.inf
        stall = 0
        jolts = 0
        step = 0
        while not found.value and jolts < 40:
            G = X @ X.T
            G[eye] = -1.0
            V = G - TARGET
            np.maximum(V, 0.0, out=V)
            # Penalty sum V^2 plus a sharper term focusing on the worst pairs.
            W = 2.0 * V + 40.0 * V ** 3
            grad = W @ X
            # Remove radial component (tangent-space force), then step.
            grad -= np.sum(grad * X, axis=1, keepdims=True) * X
            X = normalize(X - lr * grad)
            step += 1
            local_iters += 1

            if step % 50 == 0:
                m = max_dot(X)
                if not np.isfinite(m):
                    X = random_config(rng)
                    continue
                if m < run_best - 1e-7:
                    run_best = m
                    stall = 0
                else:
                    stall += 50
                if m < my_best:
                    my_best = m
                    with best_val.get_lock():
                        if m < best_val.value:
                            best_val.value = m
                            queue.put((m, X.copy(), wid))
                if m <= LIMIT:
                    found.value = 1
                    queue.put((m, X.copy(), wid))
                    break
                # Cooling schedule
                lr = max(lr * 0.995, 2e-4)
                if stall >= 3000:
                    # Thermal jolt: shake the most crowded points hard,
                    # everyone else gently, then reheat.
                    jolts += 1
                    rowmax = G.max(axis=1)
                    worst = np.argsort(rowmax)[-rng.integers(3, 12):]
                    X = X + 0.02 * rng.standard_normal(X.shape)
                    X[worst] += 0.25 * rng.standard_normal((len(worst), D))
                    X = normalize(X)
                    lr = 0.03
                    stall = 0
                    run_best = np.inf
            if local_iters >= 5000:
                with iters.get_lock():
                    iters.value += local_iters
                local_iters = 0


def git(*args):
    return subprocess.run(["git", "-C", REPO, *args], capture_output=True, text=True)


def push(files, msg):
    git("add", *files)
    git("commit", "-q", "-m", msg + "\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n"
        "Claude-Session: https://claude.ai/code/session_01VFkLXP7kJjth38Z8secrus")
    for delay in (0, 2, 4, 8, 16):
        time.sleep(delay)
        r = git("push", "-u", "origin", BRANCH)
        if r.returncode == 0:
            return True
    print("push failed:", r.stderr, flush=True)
    return False


def log(line):
    with open(PROGRESS, "a") as f:
        f.write(line + "\n")
    print(line, flush=True)


def save_matrix(path, X, header):
    np.savetxt(path, X, fmt="%.17f", header=header)


def main():
    best_val = mp.Value("d", np.inf)
    iters = mp.Value("q", 0)
    found = mp.Value("i", 0)
    queue = mp.Queue()
    nworkers = os.cpu_count() or 4
    seeds = np.random.SeedSequence().spawn(nworkers)
    procs = [mp.Process(target=worker, args=(i, s, best_val, iters, found, queue), daemon=True)
             for i, s in enumerate(seeds)]
    for p in procs:
        p.start()

    best_X, best_m = None, np.inf
    hist = []
    stop_reason = "time limit reached"
    last_log = time.time()
    log(f"# start {dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')} "
        f"N={N} D={D} workers={nworkers} end={END_TIME.isoformat()}")
    log("# timestamp_utc, iterations, best_max_dot, violation(best_max_dot-0.5)")
    push(["progress.txt", "kissing6.py"], "Start 6D kissing search (N=73)")

    try:
        while True:
            try:
                m, X, wid = queue.get(timeout=1.0)
                Xv = normalize(X.astype(np.float64))
                mv = max_dot(Xv)
                if mv < best_m:
                    best_m, best_X = mv, Xv
                    hist.append((time.time(), mv, wid))
                if mv <= LIMIT:
                    stop_reason = "BREAKTHROUGH"
                    break
            except Exception:
                pass
            now = dt.datetime.now(dt.timezone.utc)
            if now >= END_TIME:
                break
            if time.time() - last_log >= LOG_INTERVAL:
                last_log = time.time()
                log(f"{now.isoformat(timespec='seconds')}, {iters.value}, "
                    f"{best_m:.10f}, {best_m - LIMIT:+.10f}")
                if best_X is not None:
                    save_matrix(BEST_FILE, best_X, f"best so far, max dot = {best_m:.12f}")
                push(["progress.txt", "best_configuration.txt"], f"Progress: best max dot {best_m:.6f}")
    except KeyboardInterrupt:
        stop_reason = "interrupted"
    finally:
        found.value = 1
        now = dt.datetime.now(dt.timezone.utc)
        files = ["progress.txt"]
        if best_X is not None:
            save_matrix(BEST_FILE, best_X, f"best so far, max dot = {best_m:.12f}")
            files.append("best_configuration.txt")
        if stop_reason == "BREAKTHROUGH":
            norms = np.linalg.norm(best_X, axis=1)
            save_matrix(BREAKTHROUGH, best_X,
                        f"73 unit vectors in R^6, max pairwise dot = {best_m:.17f}, "
                        f"norm range [{norms.min():.17f}, {norms.max():.17f}]")
            files.append("breakthrough_coordinates.txt")
        log(f"{now.isoformat(timespec='seconds')}, {iters.value}, {best_m:.10f}, "
            f"{best_m - LIMIT:+.10f}  # STOP: {stop_reason}")
        push(files, f"Final: {stop_reason}, best max dot {best_m:.6f}")
        for p in procs:
            p.terminate()
        with open(os.path.join(REPO, "run_summary.txt"), "w") as f:
            f.write(f"stop_reason={stop_reason}\nbest_max_dot={best_m!r}\niters={iters.value}\n")
            for t, mv, wid in hist:
                f.write(f"{dt.datetime.fromtimestamp(t, dt.timezone.utc).isoformat(timespec='seconds')} "
                        f"worker={wid} max_dot={mv:.10f}\n")


if __name__ == "__main__":
    signal.signal(signal.SIGTERM, lambda *a: (_ for _ in ()).throw(KeyboardInterrupt()))
    main()
