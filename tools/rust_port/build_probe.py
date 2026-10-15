"""Build validation probes linked to the already-built exact reference library."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
from reference import ROOT


def build(compiler):
    compiler = Path(compiler).resolve()
    out = ROOT / '.reference'
    env = dict(os.environ)
    env['PATH'] = str(compiler.parent) + os.pathsep + env['PATH']
    source = ROOT / 'tools/rust_port/reference_probe.f90'
    library = out / ('thermopack.dll' if os.name == 'nt' else 'libthermopack.so')
    target = out / ('reference_probe.dll' if os.name == 'nt' else 'libreference_probe.so')
    assert library.exists(), 'Build the exact-source reference first'
    flags = ['-O3', '-fdefault-real-8', '-fdefault-double-8', '-frecursive']
    if os.name != 'nt':
        flags += ['-fPIC', '-Wl,-rpath,$ORIGIN']
    subprocess.run([str(compiler), *flags, '-shared', '-I', str(out), str(source),
                    str(library), '-o', str(target)], cwd=out, env=env, check=True)
    metadata = {'source_sha256': hashlib.sha256(source.read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
                'library_sha256': hashlib.sha256(library.read_bytes()).hexdigest(),
                'probe_sha256': hashlib.sha256(target.read_bytes()).hexdigest(),
                'compiler': subprocess.check_output([str(compiler), '--version'], text=True).splitlines()[0],
                'flags': flags}
    (out / 'probe.metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(f'Built {target}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--compiler', required=True)
    build(parser.parse_args().compiler)
