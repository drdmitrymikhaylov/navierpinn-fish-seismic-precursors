"""Behavioural state of a swimming fish, measured from bout kinematics.

Larval zebrafish do not swim continuously.  They move in discrete bouts --
a short tail beat, a displacement, a turn -- separated by roughly a second of
stillness.  Every behavioural metric below is therefore defined on bouts
rather than on frames, which is what makes them comparable between recordings
of different length and different frame rates.

Source data: 523 recordings of spontaneous (undisturbed) swimming from 39
larval zebrafish, 16 147 bouts in total, each bout carrying its time, its
position, its heading, its displacement and its turn angle.

The metrics are the ones a continuous-monitoring system would compute, and
they are grouped by what they report:

  activity     bout rate, displacement per bout, distance covered
  turning      dispersion of turn angle, fraction of large turns, and the
               persistence of turn direction from one bout to the next
  space use    distance from the arena centre, which is the standard
               thigmotaxis measure, and how much of the arena is visited
"""

import numpy as np

LARGE_TURN_DEG = 16.5        # the 90th percentile of the pooled |turn| distribution
                             # (14 870 bouts), so "large" is defined by the fish
                             # rather than by a round number


def load(path="data/spontaneous_swim.npz"):
    """Return a list of per-trial dicts, one per recording."""
    z = np.load(path)
    # read every array once: indexing a compressed npz column by column
    # re-inflates the whole array each time and turns a 0.2 s load into 25 s
    T_, X_, Y_, A_, R_, TT_ = (z["TimeBout"], z["x"], z["y"], z["Angle"],
                               z["R"], z["T"])
    fish = z["FishN"].ravel()
    trials = []
    for j in range(T_.shape[1]):
        t = T_[:, j]
        ok = np.isfinite(t)
        n = int(ok.sum())
        if n < 5:
            continue
        r_col, t_col = R_[:, j], TT_[:, j]
        trials.append({
            "fish": int(fish[j]),
            "trial": j,
            "t": t[ok],
            "x": X_[:, j][ok],
            "y": Y_[:, j][ok],
            "angle": A_[:, j][ok],
            # R and T are defined between consecutive bouts, so they are one
            # shorter; take the same leading n-1 entries
            "r": r_col[np.isfinite(r_col)][: max(n - 1, 0)],
            "turn": t_col[np.isfinite(t_col)][: max(n - 1, 0)],
        })
    return trials


def valid(tr, min_bouts=8, min_duration=5.0, arena=(0.0, 1000.0)):
    """Which recordings are usable.

    A recording with a handful of bouts cannot support a rate estimate, and a
    track that leaves the arena bounds is a tracking failure rather than a
    fish.  Both are rejected before anything is measured, because a monitoring
    system that quietly averages over failed tracking will report the failure
    as behaviour.
    """
    if len(tr["t"]) < min_bouts:
        return False
    if tr["t"].max() - tr["t"].min() < min_duration:
        return False
    if not (np.all(tr["x"] >= arena[0]) and np.all(tr["x"] <= arena[1])
            and np.all(tr["y"] >= arena[0]) and np.all(tr["y"] <= arena[1])):
        return False
    if len(tr["r"]) < min_bouts - 1 or len(tr["turn"]) < min_bouts - 1:
        return False
    return True


def circular_dispersion(deg):
    """1 - R, the circular variance of a set of angles.  0 = all identical."""
    a = np.deg2rad(np.asarray(deg, float))
    return float(1.0 - np.abs(np.mean(np.exp(1j * a))))


def state(tr, centre=None):
    """The behavioural state vector of one recording."""
    t, x, y = tr["t"], tr["x"], tr["y"]
    r, turn = np.abs(tr["r"]), tr["turn"]
    dur = float(t.max() - t.min())
    ibi = np.diff(t)
    if centre is None:
        centre = (500.0, 500.0)
    d_centre = np.hypot(x - centre[0], y - centre[1])

    # persistence: does a left turn tend to be followed by another left turn?
    s = np.sign(turn)
    persistence = float(np.mean(s[1:] == s[:-1])) if len(s) > 2 else np.nan

    return {
        "bout_rate": len(t) / dur,
        "ibi_median": float(np.median(ibi)),
        "ibi_cv": float(np.std(ibi) / np.mean(ibi)) if np.mean(ibi) > 0 else np.nan,
        "displacement_median": float(np.median(r)),
        "distance_per_s": float(np.sum(r) / dur),
        "turn_dispersion": circular_dispersion(turn),
        "turn_abs_median": float(np.median(np.abs(turn))),
        "large_turn_fraction": float(np.mean(np.abs(turn) > LARGE_TURN_DEG)),
        "turn_persistence": persistence,
        "centre_distance_median": float(np.median(d_centre)),
        "area_explored": float(np.ptp(x) * np.ptp(y)) / 1e6,
    }


METRICS = ["bout_rate", "ibi_median", "ibi_cv", "displacement_median",
           "distance_per_s", "turn_dispersion", "turn_abs_median",
           "large_turn_fraction", "turn_persistence",
           "centre_distance_median", "area_explored"]

GROUPS = {
    "activity": ["bout_rate", "ibi_median", "ibi_cv", "displacement_median",
                 "distance_per_s"],
    "turning": ["turn_dispersion", "turn_abs_median", "large_turn_fraction",
                "turn_persistence"],
    "space use": ["centre_distance_median", "area_explored"],
}


def state_matrix(trials):
    rows = [state(tr) for tr in trials]
    X = np.array([[r[m] for m in METRICS] for r in rows], float)
    fish = np.array([tr["fish"] for tr in trials])
    return X, fish, METRICS
