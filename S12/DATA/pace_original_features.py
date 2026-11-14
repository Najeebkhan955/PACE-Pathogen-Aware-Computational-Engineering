"""Peptide descriptor engine (PACE-70 feature set).

All scales are literature-standard and hard-coded here so the pipeline is fully
reproducible offline. No values are invented: sources are named in each block.
"""
import numpy as np

AA = list("ACDEFGHIKLMNPQRSTVWY")
AAset = set(AA)

# Kyte & Doolittle (1982) J Mol Biol 157:105-132 hydropathy
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2,
          I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9,
          Y=-1.3, V=4.2)

# Eisenberg et al. (1984) consensus normalised hydrophobicity
EIS = dict(A=0.62, R=-2.53, N=-0.78, D=-0.90, C=0.29, Q=-0.85, E=-0.74, G=0.48,
           H=-0.40, I=1.38, L=1.06, K=-1.50, M=0.64, F=1.19, P=0.12, S=-0.18,
           T=-0.05, W=0.81, Y=0.26, V=1.08)

# Boman (2003) J Intern Med 254:197-215 protein-binding potential (kcal/mol)
BOMAN = dict(L=-4.92, I=-4.92, V=-4.04, F=-2.98, M=-2.35, W=-2.33, A=-1.81,
             C=-1.28, G=-0.94, Y=0.14, T=2.57, S=3.40, H=4.66, Q=5.54, K=5.55,
             N=6.64, E=6.81, D=8.72, R=14.92, P=0.0)

# pKa values (EMBOSS set, as used by Bio.SeqUtils.IsoelectricPoint)
PKA_SIDE = dict(D=3.65, E=4.25, C=8.18, Y=10.07, H=6.00, K=10.53, R=12.48)
PK_NTERM, PK_CTERM = 8.6, 3.6

# Flexibility (Vihinen et al. 1994) average flexibility index
FLEX = dict(A=0.360, C=0.350, D=0.510, E=0.500, F=0.310, G=0.540, H=0.320,
            I=0.460, K=0.470, L=0.370, M=0.300, N=0.460, P=0.510, Q=0.490,
            R=0.530, S=0.510, T=0.440, V=0.390, W=0.310, Y=0.420)

GROUPS = {
    "tiny": set("ACGST"), "small": set("ABCDGNPSTV"), "aliphatic": set("AILV"),
    "aromatic": set("FHWY"), "nonpolar": set("ACFGILMPVWY"), "polar": set("DEHKNQRSTZ"),
    "charged": set("DEHKR"), "basic": set("HKR"), "acidic": set("DE"),
    "cysteine": set("C"), "proline": set("P"), "glycine": set("G"),
}

# CTD attribute groupings (Dubchak et al. 1995; Govindan & Nair 2011)
CTD_ATTR = {
    "hydrophobicity": ("RKEDQN", "GASTPHY", "CLVIMFW"),
    "vdw_volume":     ("GASTPDC", "NVEQIL", "MHKFRYW"),
    "polarity":       ("LIFWCMVY", "PATGS", "HQRKNED"),
    "polarizability": ("GASDT", "CPNVEQIL", "KMHFRYW"),
    "charge":         ("KR", "ANCQGHILMFPSTWYV", "DE"),
    "ss_propensity":  ("EALMQKRH", "VIYCWFT", "GNPSD"),
    "sa_accessibility": ("ALFCGIVW", "RKQEND", "MPSTHY"),
}

GLOBAL_NAMES = [
    "length", "net_charge_pH7", "charge_density", "pI", "gravy_kd",
    "hydrophobicity_eis", "hydrophobic_moment", "amphipathicity_ratio",
    "aliphatic_index", "boman_index", "instability_like", "flexibility",
    "frac_tiny", "frac_small", "frac_aliphatic", "frac_aromatic", "frac_nonpolar",
    "frac_polar", "frac_charged", "frac_basic", "frac_acidic", "frac_cys",
    "frac_pro", "frac_gly", "max_hydrophobic_run", "max_cationic_run",
    "kr_ratio", "hydromoment_x_charge",
]


def clean(seq):
    s = "".join(ch for ch in str(seq).upper() if ch.isalpha())
    return s if s and all(c in AAset for c in s) else None


def net_charge(seq, pH=7.4):
    q = 1.0 / (1.0 + 10 ** (pH - PK_NTERM)) - 1.0 / (1.0 + 10 ** (PK_CTERM - pH))
    for aa, pk in PKA_SIDE.items():
        n = seq.count(aa)
        if not n:
            continue
        if aa in "KRH":
            q += n / (1.0 + 10 ** (pH - pk))
        else:
            q -= n / (1.0 + 10 ** (pk - pH))
    return q


def isoelectric_point(seq):
    """Bio.SeqUtils.IsoelectricPoint (EMBOSS pKa set) — reference implementation."""
    from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
    return float(IP(seq).pi())


def hydrophobic_moment(seq, angle=100, window=11):
    """Eisenberg mean hydrophobic moment, max over sliding windows."""
    n = len(seq)
    w = min(window, n)
    best = 0.0
    rad = np.deg2rad(angle)
    for i in range(n - w + 1):
        sub = seq[i:i + w]
        h = np.array([EIS[c] for c in sub])
        idx = np.arange(w)
        mu = np.sqrt((h * np.sin(idx * rad)).sum() ** 2 + (h * np.cos(idx * rad)).sum() ** 2)
        best = max(best, mu / w)
    return best


def aliphatic_index(seq):
    n = len(seq)
    return 100.0 * (seq.count("A") / n + 2.9 * seq.count("V") / n +
                    3.9 * (seq.count("I") + seq.count("L")) / n)


def instability_like(seq):
    """Dipeptide-instability proxy: mean Boman difference across adjacent pairs.
    Reported as a relative index, not the Guruprasad II value."""
    if len(seq) < 2:
        return 0.0
    return float(np.mean([abs(BOMAN[a] - BOMAN[b]) for a, b in zip(seq, seq[1:])]))


def max_run(seq, members):
    best = cur = 0
    for c in seq:
        cur = cur + 1 if c in members else 0
        best = max(best, cur)
    return best


def global_features(seq):
    n = len(seq)
    q = net_charge(seq)
    mu = hydrophobic_moment(seq)
    h = float(np.mean([EIS[c] for c in seq]))
    f = {g: sum(seq.count(c) for c in mem) / n for g, mem in GROUPS.items()}
    return [
        n, q, q / n, isoelectric_point(seq),
        float(np.mean([KD[c] for c in seq])), h, mu,
        mu / (abs(h) + 1e-6),
        aliphatic_index(seq),
        float(np.mean([BOMAN[c] for c in seq])),
        instability_like(seq),
        float(np.mean([FLEX[c] for c in seq])),
        f["tiny"], f["small"], f["aliphatic"], f["aromatic"], f["nonpolar"],
        f["polar"], f["charged"], f["basic"], f["acidic"], f["cysteine"],
        f["proline"], f["glycine"],
        max_run(seq, set("AILVFWMC")), max_run(seq, set("KR")),
        (seq.count("K") + 1) / (seq.count("R") + 1),
        mu * q,
    ]


def aac(seq):
    n = len(seq)
    return [seq.count(a) / n for a in AA]


def ctd_ct(seq):
    """CTD composition (3) + transition (3) per attribute."""
    out = []
    for attr, groups in CTD_ATTR.items():
        code = []
        for c in seq:
            for gi, g in enumerate(groups):
                if c in g:
                    code.append(gi)
                    break
            else:
                code.append(1)
        code = np.array(code)
        n = len(code)
        out += [float((code == g).sum()) / n for g in range(3)]
        if n > 1:
            pairs = list(zip(code[:-1], code[1:]))
            tot = len(pairs)
            for a, b in [(0, 1), (0, 2), (1, 2)]:
                out.append(sum(1 for x, y in pairs if {x, y} == {a, b}) / tot)
        else:
            out += [0.0, 0.0, 0.0]
    return out


CTD_NAMES = []
for _a in CTD_ATTR:
    CTD_NAMES += [f"{_a}_C{g}" for g in (1, 2, 3)]
    CTD_NAMES += [f"{_a}_T{t}" for t in ("12", "13", "23")]

FEATURE_NAMES = GLOBAL_NAMES + [f"AAC_{a}" for a in AA] + CTD_NAMES


def featurize(seq):
    return global_features(seq) + aac(seq) + ctd_ct(seq)


def featurize_many(seqs):
    return np.asarray([featurize(s) for s in seqs], dtype=float)


if __name__ == "__main__":
    s = "GIGKFLHSAKKFGKAFVGEIMNS"  # magainin-2
    v = featurize(s)
    print("n_features =", len(v), "== len(names)", len(FEATURE_NAMES))
    for k in ["length", "net_charge_pH7", "pI", "gravy_kd", "hydrophobic_moment"]:
        print(f"  {k:22s} {v[FEATURE_NAMES.index(k)]:.3f}")
