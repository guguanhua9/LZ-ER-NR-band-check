"""Portable redraw of current ER/NR500 comparison. See README_zh.txt."""
from pathlib import Path
import argparse, csv, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--statistics-dir',type=Path,default=BASE/'data/precomputed')
parser.add_argument('--reference-dir',type=Path,default=BASE/'data/published')
parser.add_argument('--output-dir',type=Path,default=BASE/'output')
args=parser.parse_args()
OUT=args.output_dir.resolve(); OUT.mkdir(parents=True,exist_ok=True)
ER=args.statistics_dir/'er_band.csv'; NR=args.statistics_dir/'nr_band.csv'; REF=args.reference_dir
metadata_path=args.statistics_dir/'statistics_metadata.json'
metadata=json.loads(metadata_path.read_text(encoding='utf-8-sig')) if metadata_path.exists() else {}
paths=[ER,NR]+[REF/f'pdf_fig2_{s}_{k}.csv' for s in ('er','nr') for k in ('mean','lower','upper')]
if metadata_path.exists(): paths.append(metadata_path)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
before={str(p.relative_to(BASE)) if p.is_relative_to(BASE) else p.name:sha(p) for p in paths}
def read_table(p):
    return np.genfromtxt(p, delimiter=',', names=True, encoding='utf-8-sig')

curves = {s: {k: np.loadtxt(REF/f'pdf_fig2_{s.lower()}_{k}.csv', delimiter=',', skiprows=1)
              for k in ('mean', 'lower', 'upper')} for s in ('ER', 'NR')}
for spec in curves.values():
    for c in spec.values():
        assert np.isfinite(c).all() and np.all(np.diff(c[:, 0]) >= 0)

def bounded_interp(x, c):
    result = np.full(len(x), np.nan)
    ok = (x >= c[0, 0]) & (x <= c[-1, 0])
    result[ok] = np.interp(x[ok], c[:, 0], c[:, 1])
    return result

data = {}
for species, path in (('ER', ER), ('NR', NR)):
    t = read_table(path)
    if species == 'ER':
        x, mu, sd = t['s1c_center_phd'], t['mean_log10_s2c'], t['std_log10_s2c_ddof0']
        lower, upper = t['minus1sigma'], t['plus1sigma']
        q10, q90 = t['q10'], t['q90']
        saved_sigma_ratio, saved_quantile_ratio = t['sigma1_width_over_published'], t['q10q90_width_over_published']
        assert np.all(t['slice_full_width_phd'] == 20)
    else:
        x, mu, sd = t['s1c_phd'], t['mean_log10_s2c'], t['std_log10_s2c_ddof0']
        lower, upper = t['minus1sigma_log10_s2c'], t['plus1sigma_log10_s2c']
        q10, q90 = t['q10_log10_s2c'], t['q90_log10_s2c']
        saved_sigma_ratio, saved_quantile_ratio = t['pm1sigma_width_over_published'], t['q10q90_width_over_published']
        assert np.all(np.isin(x, np.arange(5., 800., 10.)))
    assert np.allclose(lower, mu-sd, rtol=0, atol=1e-14)
    assert np.allclose(upper, mu+sd, rtol=0, atol=1e-14)
    pl = bounded_interp(x, curves[species]['lower'])
    pu = bounded_interp(x, curves[species]['upper'])
    width = pu-pl
    rs, rq = 2*sd/width, (q90-q10)/width
    assert np.allclose(rs, saved_sigma_ratio, rtol=0, atol=2e-13, equal_nan=True)
    assert np.allclose(rq, saved_quantile_ratio, rtol=0, atol=2e-13, equal_nan=True)
    data[species] = dict(x=x, mean=mu, sigma=sd, lower=lower, upper=upper,
                         q10=q10, q90=q90, sigma_ratio=rs, quantile_ratio=rq)
    with (OUT/f'plotted_{species.lower()}_series.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(('s1c_phd','mc_mean','mc_sigma','mc_mean_minus_sigma','mc_mean_plus_sigma',
                         'mc_q10','mc_q90','published_mean','published_lower','published_upper',
                         'sigma_full_width_over_published','quantile_full_width_over_published'))
        writer.writerows(zip(x,mu,sd,lower,upper,q10,q90,
                             bounded_interp(x,curves[species]['mean']),pl,pu,rs,rq))

plt.rcParams.update({'font.size': 12, 'axes.labelsize': 13, 'axes.titlesize': 13,
                     'legend.fontsize': 9.4, 'svg.fonttype': 'none', 'pdf.fonttype': 42})

def draw_main(ax, species):
    d = data[species]
    color = '#2166ac' if species == 'ER' else '#b2182b'
    x = d['x']
    ax.fill_between(x, d['lower'], d['upper'], color=color, alpha=.09, zorder=1)
    for key in ('lower', 'upper'):
        ax.plot(x, d[key], color=color, lw=2.2, zorder=3)
    ax.plot(x, d['mean'], color=color, lw=1.35, alpha=.92, zorder=5)
    for key in ('q10', 'q90'):
        ax.plot(x, d[key], color='.50', lw=1.55, ls=':', zorder=2)
    for key in ('lower', 'upper'):
        c = curves[species][key]
        ax.plot(c[:,0], c[:,1], color='black', lw=1.55, ls='--', zorder=4)
    c = curves[species]['mean']
    ax.plot(c[:,0], c[:,1], color='black', lw=.9, alpha=.9, zorder=4)
    handles = [Line2D([],[],color=color,lw=1.35,label='This work: slice mean'),
               Line2D([],[],color='black',lw=.9,label='Published mean'),
               Line2D([],[],color=color,lw=2.2,label=r'This work: $\mu\pm\sigma$'),
               Line2D([],[],color='.50',lw=1.55,ls=':',label='This work: q10 / q90'),
               Line2D([],[],color='black',lw=1.55,ls='--',label='Published q10 / q90')]
    ax.legend(handles=handles, loc='lower right', framealpha=.96)
    ax.set(xlim=(0,800), ylim=(3.25,5.4) if species=='ER' else (2.75,4.2),
           ylabel=r'$\log_{10}(S2_c/\mathrm{phd})$')
    ax.set_title('ER: 20-phd S1c slices' if species=='ER'
                 else 'NR: flat-spectrum endpoint 500 keVnr; 10-phd S1c slices', pad=11)
    ax.set_xticks(np.arange(0,801,100))
    ax.grid(alpha=.18)
    for mark in (100,600):
        ax.axvline(mark,color='.63',lw=.75,ls='--',zorder=0)

def draw_ratio(ax, species):
    d=data[species]
    color='#2166ac' if species=='ER' else '#b2182b'
    ax.plot(d['x'],d['sigma_ratio'],color=color,lw=1.9,label=r'$2\sigma$ / published width')
    ax.plot(d['x'],d['quantile_ratio'],color='.50',lw=1.55,ls=':',label='(q90 - q10) / published width')
    ax.axhline(1,color='black',lw=.8,ls='--')
    ax.set(xlim=(0,800),ylim=(.65,1.85),xlabel='S1c [phd]',ylabel='Full-width ratio')
    ax.set_xticks(np.arange(0,801,100))
    ax.legend(loc='upper right',framealpha=.96)
    ax.grid(alpha=.18)
    for mark in (100,600):
        ax.axvline(mark,color='.63',lw=.75,ls='--',zorder=0)

fig, axes=plt.subplots(2,2,figsize=(14.5,9.0),sharex='col',
                      gridspec_kw={'height_ratios':[2.25,1]})
for col,species in enumerate(('ER','NR')):
    draw_main(axes[0,col],species)
    draw_ratio(axes[1,col],species)
fig.suptitle(r'ER and updated NR: mean, $\mu\pm\sigma$, and 10%–90% bands',fontsize=17,y=.987)
fig.text(.5,.019,
         r'No empirical scaling. Same MC for $\mu\pm\sigma$ and q10/q90; '
         r'$\sigma=\mathrm{std}[\log_{10}(S2_c/\mathrm{phd})]$ (ddof=0). '
         'Published boundaries retain their q10/q90 labels.',ha='center',fontsize=10)
fig.subplots_adjust(left=.075,right=.985,top=.895,bottom=.105,wspace=.19,hspace=.075)
for ext in ('png','pdf','svg'):
    fig.savefig(OUT/f'ER_NR500_combined.{ext}',dpi=240,facecolor='white')
plt.close(fig)

# A compact two-panel copy is convenient for messages or a slide.
fig,axes=plt.subplots(1,2,figsize=(14.5,6.35))
for ax,species in zip(axes,('ER','NR')):
    draw_main(ax,species)
    ax.set_xlabel('S1c [phd]')
fig.suptitle(r'ER and updated NR: mean, $\mu\pm\sigma$, and 10%–90% bands',fontsize=17,y=.985)
fig.text(.5,.025,'Fixed response parameters; no empirical scaling. '
         'NR flat-spectrum endpoint: 500 keVnr. Published boundaries retain q10/q90 labels.',
         ha='center',fontsize=10)
fig.subplots_adjust(left=.075,right=.985,top=.855,bottom=.16,wspace=.19)
for ext in ('png','pdf'):
    fig.savefig(OUT/f'ER_NR500_main_panels.{ext}',dpi=240,facecolor='white')
plt.close(fig)

after={str(p.relative_to(BASE)) if p.is_relative_to(BASE) else p.name:sha(p) for p in paths}
assert before==after
record={'status':'complete','scope':'Data-native redraw, not raster editing; no new MC or response changes.',
        'ER_input':'ER flat0.1-350keVee preset; input event counts are in statistics_metadata.',
        'NR_input':'NR flat0.1-500keVnr preset; input event counts are in statistics_metadata.',
        'cut_interpretation':'The previously adjusted quantity is the input flat-spectrum energy endpoint, not a new event-selection cut.',
        'mean_definition':'Arithmetic mean of log10(S2c/phd) within each S1c slice.',
        'sigma_definition':'Population standard deviation of log10(S2c/phd), ddof=0.',
        'empirical_scale':1,'published_definition_changed':False,'new_MC':False,'parameters_changed':False,
        'statistics_metadata':metadata,'input_hashes':before,'all_inputs_unchanged':before==after,
        'published_curve_support':{s:{k:[float(c[0,0]),float(c[-1,0])]
                                      for k,c in cc.items()} for s,cc in curves.items()},
        'support_policy':'Original vector paths drawn as given; width ratios only interpolate within both boundaries support; no endpoint extrapolation.',
        'point_counts':{s:len(d['x']) for s,d in data.items()},
        'outputs':{p.name:sha(p) for p in sorted(OUT.iterdir()) if p.suffix in ('.png','.pdf','.svg','.csv')}}
(OUT/'plot_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(json.dumps({'point_counts':record['point_counts'],'all_inputs_unchanged':True,
                  'outputs':list(record['outputs'])},indent=2))
