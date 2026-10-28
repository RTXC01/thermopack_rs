"""Single bounded exact-source PC-SAFT F_n oracle process."""
import ctypes, json, os, sys
from pathlib import Path
source = Path(sys.argv[1]); compiler = Path(sys.argv[2]); request = json.loads(sys.stdin.readline())
sys.path.insert(0, str(source / 'tools/rust_port'))
import reference
reference.load(compiler)
from thermopack.pcsaft import pcsaft
model = pcsaft(request['components'], parameter_reference=request['reference'], simplified=request['simplified'])
model.activate()
probe = ctypes.CDLL(str(Path(__file__).resolve().parents[2] / '.reference/pcsaft-mole/pcsaft_mole_probe.dll'))
f = probe.pcsaft_mole_gradient
f.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int)]
f.restype = None
n = request['moles']; amounts = (ctypes.c_double * 16)(*n, *([0.0] * (16-len(n))))
values = (ctypes.c_double * 16)(); status = ctypes.c_int()
f(request['temperature'], request['volume'], amounts, values, ctypes.byref(status))
if status.value: raise RuntimeError(f'unsupported active eos type: {status.value}')
print('RESULT='+json.dumps(list(values)[:len(n)], allow_nan=False), flush=True)
