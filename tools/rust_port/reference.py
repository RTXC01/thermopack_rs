"""Load the checked-in wrapper against the exact-source local build."""
import os
from pathlib import Path
import sys
import types
import json
import subprocess

ROOT = Path(__file__).resolve().parents[2]
_DLL_HANDLES = []


def load(compiler_directory=None):
    if compiler_directory and os.name == "nt":
        _DLL_HANDLES.append(os.add_dll_directory(str(Path(compiler_directory).resolve())))
    library = ROOT / ".reference" / ("thermopack.dll" if os.name == "nt" else "libthermopack.so")
    if not library.exists():
        raise RuntimeError("Build the exact-source reference first")
    settings = types.ModuleType("thermopack.platform_specifics")
    settings.DIFFERENTIAL_RETURN_MODE = "v2"
    settings.get_platform_specifics = lambda: {
        "os_id": "win" if os.name == "nt" else "linux", "prefix": "__",
        "module": "_MOD_", "postfix": "", "postfix_no_module": "_",
        "dyn_lib": str(library), "diff_return_mode": "v2"}
    sys.modules[settings.__name__] = settings
    sys.path.insert(0, str(ROOT / "addon/pycThermopack"))
    from thermopack.cubic import cubic
    return cubic


class ReferenceProcess:
    """Persistent oracle isolated from fatal Fortran STOP/exit calls.

    A terminal worker produces a recorded reference failure. The following
    independent case gets a new worker; failed cases are never silently retried.
    """
    def __init__(self, compiler_directory=None):
        self.compiler_directory=compiler_directory
        self.process=None

    def call(self, request):
        if self.process is None or self.process.poll() is not None:
            command=[sys.executable,"-u",str(ROOT/"tools/rust_port/oracle_worker.py")]
            if self.compiler_directory:
                command.extend(["--compiler-directory",self.compiler_directory])
            self.process=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                                          stderr=subprocess.STDOUT,text=True,encoding="utf-8",errors="replace")
        process=self.process
        process.stdin.write(json.dumps(request)+'\n')
        process.stdin.flush()
        output=[]
        while True:
            line=process.stdout.readline()
            if not line:
                code=process.wait()
                if any('Traceback (most recent call last)' in message for message in output):
                    raise RuntimeError('\n'.join(output))
                return {"reference_error":{"exit_code":code,"output":output}}
            if line.startswith("THERMOPACK_RESULT="):
                return json.loads(line.split('=',1)[1])
            output.append(line.rstrip())

    def close(self):
        if self.process is not None and self.process.poll() is None:
            self.process.stdin.close()
            self.process.wait()

    def __enter__(self):
        return self

    def __exit__(self,*args):
        self.close()
