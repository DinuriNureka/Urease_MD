#!/usr/bin/env python3
"""Equilibration convergence checks: density (mdout), backbone RMSD drift, Ni–Ni stability (cpptraj), temperature (mdout)."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import math

import numpy as np

# Defaults (overridden by config.env / CLI)
_DEF = {
    "density_min": 0.990,
    "density_max": 1.010,
    "density_std": 0.005,
    "rmsd_drift": 0.5,
    "rmsd_max": 5.0,
    "nini_min": 3.21,
    "nini_max": 3.81,
    "nini_std": 0.20,
    "window": 0.25,
    "dt_ps": 0.002,
    "ntwx": 5000,
}


def _env_float(key: str, default: float) -> float:
    v = os.environ.get(key)
    if v is None or v.strip() == "":
        return default
    try:
        return float(v)
    except ValueError as e:
        raise SystemExit(f"[GATE FAIL] Invalid {key}={v!r}: {e}") from e


def parse_temp_k_lines(mdout_path: Path) -> list[float]:
    """Parse instantaneous temperatures from mdout (TEMP(K) = ... lines)."""
    pat = re.compile(r"TEMP\s*\(\s*K\s*\)\s*=\s*(\S+)", re.I)
    try:
        text = mdout_path.read_text(errors="replace")
    except OSError as e:
        raise SystemExit(f"[GATE FAIL] Cannot read mdout {mdout_path}: {e}") from e
    out: list[float] = []
    for line in text.splitlines():
        m = pat.search(line)
        if not m:
            continue
        tok = m.group(1).strip()
        try:
            v = float(tok)
        except ValueError:
            continue
        if math.isnan(v):
            continue
        out.append(v)
    return out


def parse_density_lines(mdout_path: Path) -> list[float]:
    densities: list[float] = []
    try:
        text = mdout_path.read_text(errors="replace")
    except OSError as e:
        raise SystemExit(f"[GATE FAIL] Cannot read mdout {mdout_path}: {e}") from e
    for line in text.splitlines():
        if "Density" in line and "=" in line:
            try:
                tail = line.split("=", 1)[-1].strip()
                val = float(tail.split()[0])
                densities.append(val)
            except (ValueError, IndexError):
                pass
    return densities


def window_slice(values: list[float], window: float) -> list[float]:
    if not values:
        return []
    n = max(1, int(len(values) * window))
    return values[-n:]


def cpptraj_bin() -> str:
    exe = shutil.which("cpptraj")
    if exe:
        return exe
    amber = os.environ.get("AMBERHOME") or os.environ.get("AMBER_HOME")
    if amber:
        p = Path(amber) / "bin" / "cpptraj"
        if p.is_file():
            return str(p)
    raise RuntimeError(
        "[GATE FAIL] cpptraj not found on PATH and not under AMBERHOME/AMBER_HOME."
    )


def run_cpptraj(inp: str, cwd: Path | None = None) -> None:
    exe = cpptraj_bin()
    try:
        r = subprocess.run(
            [exe, "-i", "-"],
            input=inp,
            text=True,
            capture_output=True,
            cwd=cwd,
            timeout=7200,
        )
    except OSError as e:
        raise RuntimeError(f"[GATE FAIL] cpptraj failed to start: {e}") from e
    if r.returncode != 0:
        raise RuntimeError(
            "[GATE FAIL] cpptraj exited with code "
            f"{r.returncode}\n--- stderr ---\n{r.stderr}\n--- stdout ---\n{r.stdout}"
        )


def read_rmsd_dat(path: Path) -> tuple[np.ndarray, np.ndarray]:
    frames: list[int] = []
    rmsds: list[float] = []
    try:
        for line in path.read_text(errors="replace").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                frames.append(int(float(parts[0])))
                rmsds.append(float(parts[1]))
            except ValueError:
                continue
    except OSError as e:
        raise SystemExit(f"[GATE FAIL] Cannot read RMSD file {path}: {e}") from e
    return np.array(frames, dtype=float), np.array(rmsds, dtype=float)


def read_distance_dat(path: Path) -> tuple[np.ndarray, np.ndarray]:
    frames: list[int] = []
    dists: list[float] = []
    try:
        for line in path.read_text(errors="replace").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                frames.append(int(float(parts[0])))
                dists.append(float(parts[1]))
            except ValueError:
                continue
    except OSError as e:
        raise SystemExit(f"[GATE FAIL] Cannot read distance file {path}: {e}") from e
    return np.array(frames, dtype=float), np.array(dists, dtype=float)


def read_ntwx_dt_from_mdin(path: Path) -> tuple[float, float]:
    """Return (dt_ps, ntwx) from first &cntrl block."""
    try:
        raw = path.read_text(errors="replace")
    except OSError:
        return _DEF["dt_ps"], _DEF["ntwx"]
    chunk = []
    for line in raw.splitlines():
        if "!" in line:
            line = line.split("!", 1)[0]
        chunk.append(line)
    text = "\n".join(chunk)
    m = re.search(r"&cntrl\s*(.*?)\s*/", text, re.S | re.I)
    if not m:
        return _DEF["dt_ps"], _DEF["ntwx"]
    body = m.group(1).replace("\n", " ")
    dt_ps, ntwx = _DEF["dt_ps"], _DEF["ntwx"]
    for tok in body.split(","):
        tok = tok.strip()
        if "=" not in tok:
            continue
        k, v = tok.split("=", 1)
        k = k.strip().lower()
        v = v.strip()
        if k == "dt":
            try:
                dt_ps = float(v)
            except ValueError:
                pass
        if k == "ntwx":
            try:
                ntwx = float(v)
            except ValueError:
                pass
    return dt_ps, ntwx


def main() -> int:
    p = argparse.ArgumentParser(description="Equilibration convergence gate (density, RMSD, Ni–Ni).")
    p.add_argument("--mdout", required=True, type=Path, help="Path to equil_free mdout/out (.out or .mdout)")
    p.add_argument("--prmtop", required=True, type=Path, help="Topology matching trajectory")
    p.add_argument("--traj", required=True, type=Path, help="equil_free NetCDF trajectory")
    p.add_argument("--report", type=Path, help="Write report here (default: systems/<name>/equil_convergence_report.txt)")
    p.add_argument("--equil-mdin", type=Path, help="equil_free.in for dt/ntwx (default: project input_files/equil_free.in)")
    p.add_argument("--ni1", type=int, default=570)
    p.add_argument("--ni2", type=int, default=571)
    p.add_argument("--ohb", type=int, default=572, help="Reserved (OHB residue number; not used in Ni–Ni check)")
    p.add_argument("--window", type=float, default=None, help="Fraction of last frames to use")
    p.add_argument("--density_min", type=float, default=None)
    p.add_argument("--density_max", type=float, default=None)
    p.add_argument("--density_std", type=float, default=None)
    p.add_argument("--rmsd_drift", type=float, default=None)
    p.add_argument("--rmsd_max", type=float, default=None)
    p.add_argument("--nini_min", type=float, default=None)
    p.add_argument("--nini_max", type=float, default=None)
    p.add_argument("--nini_std", type=float, default=None)
    args = p.parse_args()
    _ = args.ohb  # reserved for future OHB-based checks

    dmin = args.density_min if args.density_min is not None else _env_float("EQUIL_DENSITY_MIN", _DEF["density_min"])
    dmax = args.density_max if args.density_max is not None else _env_float("EQUIL_DENSITY_MAX", _DEF["density_max"])
    dstd_thr = args.density_std if args.density_std is not None else _env_float("EQUIL_DENSITY_STD", _DEF["density_std"])
    rmsd_drift_thr = args.rmsd_drift if args.rmsd_drift is not None else _env_float("EQUIL_RMSD_DRIFT", _DEF["rmsd_drift"])
    rmsd_max_thr = args.rmsd_max if args.rmsd_max is not None else _env_float("EQUIL_RMSD_MAX", _DEF["rmsd_max"])
    nini_min = args.nini_min if args.nini_min is not None else _env_float("EQUIL_NINI_MIN", _DEF["nini_min"])
    nini_max = args.nini_max if args.nini_max is not None else _env_float("EQUIL_NINI_MAX", _DEF["nini_max"])
    nini_std_thr = args.nini_std if args.nini_std is not None else _env_float("EQUIL_NINI_STD", _DEF["nini_std"])
    window = args.window if args.window is not None else _env_float("EQUIL_WINDOW", _DEF["window"])
    temp_mean_lo = args.temp_min if args.temp_min is not None else _env_float("EQUIL_TEMP_MIN", _DEF["temp_min"])
    temp_mean_hi = args.temp_max if args.temp_max is not None else _env_float("EQUIL_TEMP_MAX", _DEF["temp_max"])
    temp_obs_ceiling = args.temp_max_obs if args.temp_max_obs is not None else _env_float("EQUIL_TEMP_MAXOBS", _DEF["temp_max_obs"])
    tstd_thr = args.temp_std_max if args.temp_std_max is not None else _env_float("EQUIL_TEMP_STD", _DEF["temp_std"])

    if not args.mdout.is_file():
        raise SystemExit(f"[GATE FAIL] mdout not found: {args.mdout}")
    if not args.traj.is_file():
        raise SystemExit(f"[GATE FAIL] trajectory not found: {args.traj}")
    if not args.prmtop.is_file():
        raise SystemExit(f"[GATE FAIL] prmtop not found: {args.prmtop}")

    try:
        cpptraj_bin()
    except RuntimeError as e:
        print(str(e), file=sys.stderr)
        return 1

    report_path = args.report
    if report_path is None:
        proj = os.environ.get("PROJECT_ROOT", "")
        name = os.environ.get("EQUIL_SYSTEM_NAME", "")
        if proj and name:
            report_path = Path(proj) / "systems" / name / "equil_convergence_report.txt"
        else:
            report_path = Path(args.mdout).resolve().parent / "equil_convergence_report.txt"

    # --- Check 1: density ---
    densities = parse_density_lines(args.mdout)
    dwin = window_slice(densities, window)
    if len(dwin) < 2:
        print(
            "[GATE FAIL] Density not converged: insufficient density samples in mdout "
            f"(parsed {len(densities)} lines).",
            file=sys.stderr,
        )
        dmean = float("nan")
        dstd = float("nan")
        ok_dmean = False
        ok_dstd = False
    else:
        dmean = float(np.mean(dwin))
        dstd = float(np.std(dwin, ddof=1)) if len(dwin) > 1 else 0.0
        ok_dmean = dmin <= dmean <= dmax
        ok_dstd = dstd <= dstd_thr
        if not ok_dmean or not ok_dstd:
            print(
                f"[GATE FAIL] Density not converged: mean={dmean:.4f} std={dstd:.4f}",
                file=sys.stderr,
            )
            print(
                "(target: 0.990-1.010 g/cm³, std < 0.005)",
                file=sys.stderr,
            )

    # --- Check 4: temperature (last window on T>150 K; max over all T>150 K) ---
    temp_series = parse_temp_k_lines(args.mdout)
    temps_use = [t for t in temp_series if t > 150.0 and not math.isnan(t)]
    twin = window_slice(temps_use, window)
    ok_tmean = False
    ok_tstd = False
    ok_tmax = False
    tmean = float("nan")
    tstd = float("nan")
    tmax_observed = float("nan")
    if len(temps_use) < 2:
        print(
            "[GATE FAIL] Temperature check: insufficient TEMP(K) samples (need >150 K, n="
            f"{len(temps_use)}).",
            file=sys.stderr,
        )
    elif len(twin) < 2:
        print(
            "[GATE FAIL] Temperature check: insufficient samples in last window.",
            file=sys.stderr,
        )
    else:
        tmean = float(np.mean(twin))
        tstd = float(np.std(twin, ddof=1)) if len(twin) > 1 else 0.0
        tmax_observed = float(max(temps_use))
        ok_tmean = temp_mean_lo <= tmean <= temp_mean_hi
        ok_tstd = tstd <= tstd_thr
        ok_tmax = tmax_observed < temp_obs_ceiling
        if not ok_tmean or not ok_tstd or not ok_tmax:
            print(
                f"[GATE FAIL] Temperature not stable: mean={tmean:.1f} K std={tstd:.1f} K max={tmax_observed:.1f} K",
                file=sys.stderr,
            )
            print(
                f"(target: {temp_mean_lo:.0f}-{temp_mean_hi:.0f} K, std <= {tstd_thr:.1f} K, no frame > {temp_obs_ceiling:.0f} K)",
                file=sys.stderr,
            )

    # dt / ntwx for RMSD slope
    mdin = args.equil_mdin
    if mdin is None:
        proj = os.environ.get("PROJECT_ROOT")
        if proj:
            cand = Path(proj) / "simulations" / "input_files" / "equil_free.in"
            mdin = cand if cand.is_file() else None
    dt_ps, ntwx = _DEF["dt_ps"], _DEF["ntwx"]
    if mdin is not None and mdin.is_file():
        dt_ps, ntwx = read_ntwx_dt_from_mdin(mdin)
    frames_per_ns = 1000.0 / (ntwx * dt_ps) if ntwx * dt_ps > 0 else 0.0

    drift_ns = float("nan")
    rmean = float("nan")
    ok_drift = False
    ok_rmag = False
    mn = float("nan")
    sn = float("nan")
    ok_nmean = False
    ok_nstd = False

    cpp_rmsd_ok = False
    cpp_nini_ok = False
    with tempfile.TemporaryDirectory(prefix="equil_conv_") as tmp:
        tdir = Path(tmp)
        rmsd_out = tdir / "rmsd.dat"
        nini_out = tdir / "nini.dat"
        prmtop_s = str(args.prmtop.resolve())
        traj_s = str(args.traj.resolve())

        # autoimage: first action after trajin (PBC) before RMSD / distance analysis.
        cpp_rmsd = f"""parm {prmtop_s}
trajin {traj_s}
autoimage
rmsd backbone @CA,C,N,O first out {rmsd_out}
run
"""
        cpp_nini = f"""parm {prmtop_s}
trajin {traj_s}
autoimage
distance ni_ni :{args.ni1}@NI :{args.ni2}@NI out {nini_out}
run
"""
        try:
            run_cpptraj(cpp_rmsd)
            cpp_rmsd_ok = True
        except RuntimeError as e:
            print(str(e), file=sys.stderr)
        try:
            run_cpptraj(cpp_nini)
            cpp_nini_ok = True
        except RuntimeError as e:
            print(str(e), file=sys.stderr)

        if cpp_rmsd_ok:
            _, rmsd = read_rmsd_dat(rmsd_out)
            if len(rmsd) < 2:
                print("[GATE FAIL] cpptraj produced no RMSD data.", file=sys.stderr)
            else:
                n_w = max(1, int(len(rmsd) * window))
                r_win = rmsd[-n_w:]
                idx = np.arange(len(r_win), dtype=float)
                slope_per_frame = float(np.polyfit(idx, r_win, 1)[0])
                drift_ns = slope_per_frame * frames_per_ns
                rmean = float(np.mean(r_win))
                ok_drift = abs(drift_ns) <= rmsd_drift_thr
                ok_rmag = 0.0 <= rmean <= rmsd_max_thr
                if not ok_drift:
                    print(
                        f"[GATE FAIL] RMSD not plateaued: drift={drift_ns:.2f} Å/ns in last "
                        f"{int(window * 100)}% of equil_free",
                        file=sys.stderr,
                    )
                    print("(threshold: <= 0.5 Å/ns)", file=sys.stderr)
                if not ok_rmag:
                    print(
                        f"[GATE FAIL] RMSD too large: mean={rmean:.2f} Å in last "
                        f"{int(window * 100)}% of equil_free",
                        file=sys.stderr,
                    )
                    print(
                        "(threshold: <= 5.0 Å — check for structural instability)",
                        file=sys.stderr,
                    )

        if cpp_nini_ok:
            _, dist = read_distance_dat(nini_out)
            if len(dist) < 2:
                print("[GATE FAIL] cpptraj produced no Ni–Ni distance data.", file=sys.stderr)
            else:
                n_wd = max(1, int(len(dist) * window))
                d_win = dist[-n_wd:]
                mn = float(np.mean(d_win))
                sn = float(np.std(d_win, ddof=1)) if len(d_win) > 1 else 0.0
                ok_nmean = nini_min <= mn <= nini_max
                ok_nstd = sn <= nini_std_thr
                if not ok_nmean:
                    print(
                        f"[GATE FAIL] Ni-Ni distance drifted: mean={mn:.2f} Å",
                        file=sys.stderr,
                    )
                    print(
                        "(crystal reference: 3.51 Å, tolerance: +/-0.30 Å)",
                        file=sys.stderr,
                    )
                if not ok_nstd:
                    print(
                        f"[GATE FAIL] Ni-Ni distance unstable: std={sn:.2f} Å",
                        file=sys.stderr,
                    )
                    print(
                        "(threshold: <= 0.20 Å — check metal center stability)",
                        file=sys.stderr,
                    )

    def sf(x: float, spec: str) -> str:
        return "nan" if x != x else format(x, spec)

    # --- Report table ---
    lines = [
        "\n============================================================",
        "EQUILIBRATION CONVERGENCE REPORT",
        "============================================================",
        "Check                 Value           Threshold    Status",
        "------------------------------------------------------------",
        f"Density (mean)        {sf(dmean, '.4f')} g/cm³    {dmin:.3f}-{dmax:.3f}  {'PASS' if ok_dmean else 'FAIL'}",
        f"Density (std)         {sf(dstd, '.4f')} g/cm³    <= {dstd_thr:.3f}     {'PASS' if ok_dstd else 'FAIL'}",
        f"RMSD drift            {sf(drift_ns, '.2f')} Å/ns       <= {rmsd_drift_thr:.1f}       {'PASS' if ok_drift else 'FAIL'}",
        f"RMSD magnitude        {sf(rmean, '.2f')} Å          <= {rmsd_max_thr:.1f}       {'PASS' if ok_rmag else 'FAIL'}",
        f"Ni-Ni distance (mean) {sf(mn, '.2f')} Å          {nini_min:.2f}-{nini_max:.2f}    {'PASS' if ok_nmean else 'FAIL'}",
        f"Ni-Ni distance (std)  {sf(sn, '.2f')} Å          <= {nini_std_thr:.2f}      {'PASS' if ok_nstd else 'FAIL'}",
        f"Temperature (mean)    {sf(tmean, '.2f')} K          {temp_mean_lo:.0f}-{temp_mean_hi:.0f}      {'PASS' if ok_tmean else 'FAIL'}",
        f"Temperature (std)     {sf(tstd, '.2f')} K          <= {tstd_thr:.1f}       {'PASS' if ok_tstd else 'FAIL'}",
        f"Temperature (max)     {sf(tmax_observed, '.2f')} K          < {temp_obs_ceiling:.0f}         {'PASS' if ok_tmax else 'FAIL'}",
        "------------------------------------------------------------",
    ]
    overall = (
        ok_dmean
        and ok_dstd
        and ok_drift
        and ok_rmag
        and ok_nmean
        and ok_nstd
        and ok_tmean
        and ok_tstd
        and ok_tmax
    )
    lines.append(f"OVERALL: {'PASS' if overall else 'FAIL'}")
    lines.append("============================================================\n")
    report_text = "\n".join(lines)
    print(report_text)
    try:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(report_text, encoding="utf-8")
    except OSError as e:
        print(f"[GATE FAIL] Cannot write report {report_path}: {e}", file=sys.stderr)
        return 1

    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
