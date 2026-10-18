"""Redraw the supplied bands, or build and regenerate the frozen ER/NR chain."""
from pathlib import Path
import argparse, json, os, subprocess, sys

BASE=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--simulate',action='store_true',help='Build bundled C++ sources and generate events before statistics/plotting.')
    p.add_argument('--quick',action='store_true',help='10000 events/species, min 5/slice; functional test only, not reference-width precision.')
    p.add_argument('--events-er',type=int,default=1000000)
    p.add_argument('--events-nr',type=int,default=3000000)
    p.add_argument('--wsl',metavar='DISTRIBUTION',help='Optional Windows compiler/run environment, e.g. Ubuntu-22.04. Python plotting still uses this interpreter.')
    p.add_argument('--legacy-er-rounding',action='store_true',help='Match the original ER event serialization (six significant digits). No response change.')
    p.add_argument('--output-dir',type=Path,default=BASE/'output')
    a=p.parse_args()
    out=a.output_dir.resolve(); out.mkdir(parents=True,exist_ok=True)
    commands=[]
    def run(cmd,label):
        commands.append({'step':label,'argv':[str(v) for v in cmd]})
        (out/'execution_commands.json').write_text(json.dumps(commands,indent=2),encoding='utf-8')
        log=out/'logs'/f'{label}.log'; log.parent.mkdir(exist_ok=True)
        with log.open('wb') as f:
            result=subprocess.run([str(v) for v in cmd],stdout=f,stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f'{label} failed; inspect {log}')
        print(f'{label}: complete',flush=True)
    def linux_path(path):
        path=Path(path).resolve()
        if not path.drive or len(path.drive)!=2:
            raise ValueError('WSL mode expects a Windows drive path')
        return '/mnt/'+path.drive[0].lower()+path.as_posix()[2:]
    def native(cmd):
        if not a.wsl: return cmd
        converted=[linux_path(c) if isinstance(c,Path) else c for c in cmd]
        return ['wsl','-d',a.wsl,'--exec',*converted]
    if a.simulate:
        build=out/'build_response'; events=out/'events'; events.mkdir(exist_ok=True)
        run(native(['cmake','-S',BASE/'response','-B',build]),'configure_response')
        run(native(['cmake','--build',build,'--parallel','2']),'build_response')
        runner=build/'lz_response'
        if not a.wsl and os.name=='nt':
            choices=[build/'lz_response.exe',build/'Release/lz_response.exe',build/'Debug/lz_response.exe']
            runner=next((q for q in choices if q.exists()),choices[0])
        ne,nn=(10000,10000) if a.quick else (a.events_er,a.events_nr)
        for species,n,spectrum in (('er',ne,'flat_er_0p1_350.csv'),('nr',nn,'flat_nr_0p1_500.csv')):
            run(native([runner,species,str(n),BASE/'response/inputs'/spectrum,events/f'{species}_events.csv']),f'generate_{species}')
        stats=out/'statistics'
        cmd=[sys.executable,BASE/'compute_bands.py','--er-events',events/'er_events.csv',
             '--nr-events',events/'nr_events.csv','--output-dir',stats,
             '--min-events','5' if a.quick else '100']
        if a.legacy_er_rounding: cmd.append('--legacy-er-rounding')
        run(cmd,'compute_fixed_s1_bands')
        run([sys.executable,BASE/'plot_bands.py','--statistics-dir',stats,'--output-dir',out/'figures'],'plot_regenerated_bands')
        print('Regeneration complete. Quick mode is a functional check, not a precision reproduction.' if a.quick else 'Full frozen-response regeneration complete.')
    else:
        run([sys.executable,BASE/'plot_bands.py','--output-dir',out/'figures'],'redraw_precomputed_bands')
        print('Redraw complete; no event generation or parameter changes.')

if __name__=='__main__': main()
