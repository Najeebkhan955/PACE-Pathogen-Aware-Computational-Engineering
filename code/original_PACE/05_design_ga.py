"""Stage 5: multi-objective de novo design of ESKAPE-directed peptides.

Fitness combines P(AMP), predicted broad-spectrum potency, predicted safety
(haemolysis-derived selectivity index) and novelty against every known peptide
used anywhere in this study.
"""
import numpy as np, pandas as pd, pickle, json, os, sys, collections, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))
from features import featurize_many, net_charge

W = '/home/claude/work/'
rng = np.random.default_rng(7)

# Cys (disulfide ambiguity) and Met (oxidation liability) excluded by design
ALPHABET = list("ADEFGHIKLNPQRSTVWY")
LMIN, LMAX = 12, 30
POP, GENS, ELITE = 240, 45, 24
TARGETS = ['E. faecium', 'S. aureus', 'K. pneumoniae', 'A. baumannii', 'P. aeruginosa', 'E. cloacae']

cls = pickle.load(open(W + 'model_cls_final.pkl', 'rb'))
CLS = list(cls['members'].values())
REG = {t: pickle.load(open(W + 'reg_' + t.replace('. ', '_').replace(' ', '_') + '.pkl', 'rb'))['model']
       for t in TARGETS}
HEM = pickle.load(open(W + 'reg_Haemolysis_HC50.pkl', 'rb'))['model']

# ---- novelty index over all known peptides ----
pool = set()
for f in ['cls_train.csv', 'cls_test.csv', 'cls_external.csv', 'mic_eskape.csv', 'hemolysis.csv']:
    d = pd.read_csv(W + f)
    if 'label' in d.columns:
        d = d[d.label == 1]
    pool |= set(d['sequence'].astype(str))
pool = sorted(pool)
K = 3
index = collections.defaultdict(list)
pool_k = []
for i, s in enumerate(pool):
    ks = {s[j:j + K] for j in range(len(s) - K + 1)}
    pool_k.append(ks)
    for km in ks:
        index[km].append(i)
print(f'novelty pool: {len(pool)} known peptides', flush=True)


def max_similarity(seq):
    a = {seq[j:j + K] for j in range(len(seq) - K + 1)}
    if not a:
        return 1.0
    cnt = collections.Counter()
    for km in a:
        for i in index.get(km, ()):
            cnt[i] += 1
    best = 0.0
    for i, inter in cnt.items():
        j = inter / (len(a) + len(pool_k[i]) - inter)
        if j > best:
            best = j
    return best


cache = {}


def score(seqs):
    new = [s for s in dict.fromkeys(seqs) if s not in cache]
    if new:
        X = featurize_many(new)
        pa = np.mean([m.predict_proba(X)[:, 1] for m in CLS], axis=0)
        mics = np.column_stack([REG[t].predict(X) for t in TARGETS])
        hc = HEM.predict(X)
        for i, s in enumerate(new):
            sim = max_similarity(s)
            mic_mean = float(mics[i].mean())
            mic_worst = float(mics[i].max())
            ti = float(hc[i] - mic_mean)
            q = net_charge(s)
            fit = (2.0 * pa[i] - 1.0 * mic_mean - 0.5 * mic_worst
                   + 0.6 * float(np.clip(ti, -2, 4))
                   - 2.0 * max(0.0, sim - 0.55)
                   - 0.35 * max(0.0, abs(q - 6.0) - 3.0))
            cache[s] = dict(p_amp=float(pa[i]), mic_mean=mic_mean, mic_worst=mic_worst,
                            hc50=float(hc[i]), ti=ti, sim=sim, charge=float(q), fitness=float(fit),
                            **{'mic_' + t: float(v) for t, v in zip(TARGETS, mics[i])})
    return [cache[s] for s in seqs]


def random_seq():
    return ''.join(rng.choice(ALPHABET, int(rng.integers(LMIN, LMAX + 1))))


def mutate(s):
    s = list(s)
    op = rng.random()
    if op < 0.70 or len(s) <= LMIN:
        for _ in range(int(rng.integers(1, 4))):
            s[int(rng.integers(len(s)))] = str(rng.choice(ALPHABET))
    elif op < 0.85 and len(s) < LMAX:
        s.insert(int(rng.integers(len(s) + 1)), str(rng.choice(ALPHABET)))
    else:
        del s[int(rng.integers(len(s)))]
    return ''.join(s)


def crossover(a, b):
    child = a[:int(rng.integers(1, len(a)))] + b[int(rng.integers(1, len(b))):]
    if len(child) > LMAX:
        child = child[:LMAX]
    while len(child) < LMIN:
        child += str(rng.choice(ALPHABET))
    return child


def tournament(pop, sc, k=3):
    idx = rng.integers(0, len(pop), k)
    return pop[max(idx, key=lambda i: sc[i]['fitness'])]


pop = [random_seq() for _ in range(POP)]
history = []
for g in range(1, GENS + 1):
    sc = score(pop)
    order = np.argsort([-s['fitness'] for s in sc])
    pop = [pop[i] for i in order]
    sc = [sc[i] for i in order]
    history.append(dict(generation=g, best_fitness=sc[0]['fitness'],
                        mean_fitness=float(np.mean([s['fitness'] for s in sc])),
                        best_p_amp=sc[0]['p_amp'],
                        mean_p_amp=float(np.mean([s['p_amp'] for s in sc])),
                        mean_mic=float(np.mean([s['mic_mean'] for s in sc])),
                        mean_ti=float(np.mean([s['ti'] for s in sc])),
                        mean_similarity=float(np.mean([s['sim'] for s in sc]))))
    if g % 5 == 0 or g == 1:
        h = history[-1]
        print(f"gen {g:3d} | best {h['best_fitness']:.3f} | P(AMP) {h['mean_p_amp']:.3f} "
              f"| MIC {h['mean_mic']:.3f} | TI {h['mean_ti']:.2f} | sim {h['mean_similarity']:.3f}", flush=True)
    if g == GENS:
        break
    nxt = pop[:ELITE]
    while len(nxt) < POP:
        c = crossover(tournament(pop, sc), tournament(pop, sc)) if rng.random() < 0.6 else tournament(pop, sc)
        if rng.random() < 0.9:
            c = mutate(c)
        if LMIN <= len(c) <= LMAX:
            nxt.append(c)
    pop = nxt

pd.DataFrame(history).to_csv(W + 'ga_history.csv', index=False)

df = pd.DataFrame([dict(sequence=s, length=len(s), **v) for s, v in cache.items()])
df = df.sort_values('fitness', ascending=False).reset_index(drop=True)
df.to_csv(W + 'ga_all_candidates.csv', index=False)

f1 = df[df.p_amp >= 0.90]
f2 = f1[f1.mic_worst <= np.log10(32)]
f3 = f2[f2.ti >= 1.0]
f4 = f3[f3.sim < 0.55]
lead = f4.head(25).reset_index(drop=True)
lead.insert(0, 'peptide_id', ['PACE-%02d' % i for i in range(1, len(lead) + 1)])
lead.to_csv(W + 'designed_leads.csv', index=False)

funnel = pd.DataFrame([
    dict(step='Candidates evaluated', n=len(df)),
    dict(step='P(AMP) >= 0.90', n=len(f1)),
    dict(step='Predicted MIC <= 32 uM against all six ESKAPE targets', n=len(f2)),
    dict(step='Predicted selectivity index >= 10', n=len(f3)),
    dict(step='Max 3-mer similarity to known peptides < 0.55', n=len(f4)),
    dict(step='Prioritised leads', n=len(lead))])
funnel.to_csv(W + 'design_funnel.csv', index=False)
print('\n' + funnel.to_string(index=False))
print('\nTop 10 leads:')
print(lead[['peptide_id', 'sequence', 'length', 'charge', 'p_amp', 'mic_mean', 'ti', 'sim']].head(10).round(3).to_string(index=False))
