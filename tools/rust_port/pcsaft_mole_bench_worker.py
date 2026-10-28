"""Isolated worker: initialize the original model before its native timed loop."""
import ctypes
import json
import sys
from pathlib import Path

source, compiler, helper = map(Path, sys.argv[1:4])
request = json.loads(sys.stdin.readline())
sys.path.insert(0, str(source / 'tools/rust_port'))
import reference
reference.load(compiler)
from thermopack.pcsaft import pcsaft
model = pcsaft(request['components'], parameter_reference='Default', simplified=request['simplified'])
model.activate()
dll = ctypes.CDLL(str(helper))
f = dll.pcsaft_mole_bench
D, I = ctypes.c_double, ctypes.c_int
f.argtypes = [D, D, ctypes.POINTER(D), I, I, ctypes.POINTER(D), ctypes.POINTER(D), ctypes.POINTER(D), ctypes.POINTER(I)]
f.restype = None
moles = request['moles']
amounts = (D * 16)(*moles, *([0.] * (16 - len(moles))))
vector = (D * 16)()
checksum, seconds, status = D(), D(), I()
f(request['temperature'], request['volume'], amounts, request['warmup'], request['iterations'], vector, ctypes.byref(checksum), ctypes.byref(seconds), ctypes.byref(status))
if status.value:
    raise RuntimeError('unsupported active eos type: %d' % status.value)
print('RESULT=' + json.dumps(dict(vector=list(vector)[:len(moles)], checksum=checksum.value, seconds=seconds.value), allow_nan=False), flush=True)
