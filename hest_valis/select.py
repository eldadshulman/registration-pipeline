"""Per-slide protocol selection: choose micro vs no-micro from the two QC results.

We do not force one setting on every slide. For each slide we keep whichever registration
aligns the nuclei better, with a guard against non-rigid over-fitting:

  - Primary metric: nucleus-coincidence median um. Lower wins.
  - Tie (medians within TIE_UM): take the higher density_r.
  - Over-fit guard: if the lower-um winner's density_r is more than GUARD_R below the loser's
    AND its um advantage is under GUARD_UM, take the other protocol instead.

choose() takes the two metric dicts (from concordance.compute_qc) and returns
(chosen, rule, selected_median_um, selected_density_r).
"""
TIE_UM = 0.15
GUARD_R = 0.10
GUARD_UM = 0.5


def choose(micro, nomicro):
    """micro / nomicro : metric dicts, or None if that protocol was not computed/failed.

    Returns dict: chosen ('micro'|'nomicro'), rule, sel_median_um, sel_density_r.
    """
    if micro is None and nomicro is None:
        return {"chosen": None, "rule": "no_data", "sel_median_um": None, "sel_density_r": None}
    if micro is None:
        return _pack("nomicro", "micro_unavailable", nomicro)
    if nomicro is None:
        return _pack("micro", "nomicro_unavailable", micro)

    mM, mR = _md(micro)
    nM, nR = _md(nomicro)
    if abs(mM - nM) <= TIE_UM:
        chosen, rule = ("micro", "tie->higher_density_r") if mR >= nR else ("nomicro", "tie->higher_density_r")
    else:
        if mM < nM:
            win, wR, lo, lM, lR = "micro", mR, "nomicro", nM, nR
            wM = mM
        else:
            win, wR, lo, lM, lR = "nomicro", nR, "micro", mM, mR
            wM = nM
        if (lR - wR) > GUARD_R and (lM - wM) < GUARD_UM:
            chosen, rule = lo, "overfit_guard_flip"
        else:
            chosen, rule = win, "lower_um"
    return _pack(chosen, rule, micro if chosen == "micro" else nomicro)


def _md(m):
    med = m["nucleus_coincidence"]["median_um"]
    r = m["density_r"]
    # treat NaN density_r as -1 so slides with failed density correlation don't win on that metric
    import math
    if r is None or (isinstance(r, float) and not math.isfinite(r)):
        r = -1.0
    return med, r


def _pack(chosen, rule, m):
    return {"chosen": chosen, "rule": rule,
            "sel_median_um": round(m["nucleus_coincidence"]["median_um"], 3),
            "sel_density_r": round(m["density_r"], 3)}
