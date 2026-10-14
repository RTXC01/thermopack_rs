"""Load the checked-in wrapper against the exact-source local build."""
import os
from pathlib import Path
import sys
import types

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
