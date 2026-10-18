"""Finite positive S1c/S2c pairs -> fixed-S1 means, sigma and q10/q90."""
from pathlib import Path
import argparse, csv, hashlib, json
import numpy as np

BASE = Path(__file__).resolve().parent

def bounded_interp(x, curve):
    y = np.full(np.asarray(x).shape, np.nan, dtype=float)
    ok = (x >= curve[0, 0]) & (x <= curve[-1, 0])
    y[ok] = np.interp(x[ok], curve[:, 0], curve[:, 1])
    return y

def compute_one(events_path, species, reference_dir, output_csv, min_events=100,
                legacy_er_rounding=False):
    """Erec/true energy are neither sliced nor used as acceptance variables."""
    events_path, reference_dir, output_csv = map(Path, (events_path, reference_dir, output_csv))
    species = species.upper()
    if species not in ('ER', 'NR'):
        raise ValueError('species must be ER or NR')
    with events_path.open(encoding='utf-8-sig', newline='') as f:
        names = [name.strip().lower() for name in next(csv.reader(f))]
    try:
        columns = [names.index('s1c_phd'), names.index('s2c_phd')]
    except ValueError as e:
        raise ValueError('Input needs paired s1c_phd,s2c_phd columns in phd') from e
    pairs = np.loadtxt(events_path, delimiter=',', skiprows=1, usecols=columns, ndmin=2)
    raw_count = len(pairs)
    if species == 'ER' and legacy_er_rounding:
        pairs = np.fromiter((float(format(v, '.6g')) for v in pairs.flat),
                            dtype=float, count=pairs.size).reshape(pairs.shape)
    good = np.isfinite(pairs).all(axis=1) & (pairs[:, 0] > 0) & (pairs[:, 1] > 0)
    x, s2 = pairs[good].T
    y = np.log10(s2)
    full_width = 20. if species == 'ER' else 10.
    centers = np.arange(full_width/2, 800., full_width)
    reference = {k: np.loadtxt(reference_dir/f'pdf_fig2_{species.lower()}_{k}.csv',
                               delimiter=',', skiprows=1) for k in ('mean','lower','upper')}
    for c in reference.values():
        if not np.isfinite(c).all() or not np.all(np.diff(c[:,0]) >= 0):
            raise ValueError('Reference coordinates must be finite and sorted')
    pub = {k: bounded_interp(centers, c) for k,c in reference.items()}
    rows, skipped = [], []
    for i, center in enumerate(centers):
        v = y[(x >= center-full_width/2) & (x < center+full_width/2)]
        n = len(v)
        if n < min_events:
            skipped.append({'s1c_phd': float(center), 'events': n})
            continue
        mu, sd = float(np.mean(v)), float(np.std(v, ddof=0))
        q10, q90 = map(float, np.quantile(v, [0.1,0.9], method='linear'))
        pm, pl, pu = (float(pub[k][i]) for k in ('mean','lower','upper'))
        pw = pu-pl
        rs = 2*sd/pw if np.isfinite(pw) and pw > 0 else float('nan')
        rq = (q90-q10)/pw if np.isfinite(pw) and pw > 0 else float('nan')
        rows.append(dict(recoil=species,s1c_center_phd=float(center),s1c_phd=float(center),
            slice_full_width_phd=full_width,events=n,mean_log10_s2c=mu,
            std_log10_s2c_ddof0=sd,minus1sigma=mu-sd,plus1sigma=mu+sd,q10=q10,q90=q90,
            minus1sigma_log10_s2c=mu-sd,plus1sigma_log10_s2c=mu+sd,
            q10_log10_s2c=q10,q90_log10_s2c=q90,
            published_mean=pm,published_lower=pl,published_upper=pu,published_width=pw,
            sigma1_width_over_published=rs,pm1sigma_width_over_published=rs,
            q10q90_width_over_published=rq))
    if not rows:
        raise ValueError('No S1 slice has the requested minimum event count')
    output_csv.parent.mkdir(parents=True,exist_ok=True)
    with output_csv.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    h=hashlib.sha256()
    with events_path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return {'species':species,'input_filename':events_path.name,'input_sha256':h.hexdigest(),
            'raw_events':raw_count,'finite_positive_pairs':int(good.sum()),'slice_points':len(rows),
            'slice_full_width_phd':full_width,'minimum_events':min_events,'skipped_slices':skipped,
            'legacy_ER_six_significant_digit_serialization':species=='ER' and legacy_er_rounding,
            'definition':'Fixed area-S1c half-open slices; mean/std(ddof0)/linear q10-q90 of log10(S2c/phd).',
            'extra_efficiency_weights':None,'extra_selection':None,'parameter_changes':False}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--er-events',type=Path,required=True)
    p.add_argument('--nr-events',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,default=BASE/'output/statistics')
    p.add_argument('--reference-dir',type=Path,default=BASE/'data/published')
    p.add_argument('--min-events',type=int,default=100)
    p.add_argument('--legacy-er-rounding',action='store_true',
                   help='Optional reproduction of original ER six-significant-digit event export; not a physics scale.')
    a=p.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=True)
    result={'status':'computed_from_events','parameter_changes':False}
    for s,f in (('ER',a.er_events),('NR',a.nr_events)):
        result[s]=compute_one(f,s,a.reference_dir,a.output_dir/f'{s.lower()}_band.csv',
                              a.min_events,a.legacy_er_rounding)
    (a.output_dir/'statistics_metadata.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
