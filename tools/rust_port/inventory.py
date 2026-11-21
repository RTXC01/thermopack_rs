"""Inventory the checked-in reference without importing or running it."""
import ast
import hashlib
import json
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]


def _path_sort_key(path):
    return (path.name.lower(), path.name)


def inventory():
    files = []
    for path in sorted((ROOT / "src").glob("*.f90"), key=_path_sort_key):
        source = path.read_bytes().decode("utf-8", errors="replace")
        files.append({"path": path.relative_to(ROOT).as_posix(),
                      "lines": len(source.splitlines()),
                      "sha256": hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest(),
                      "modules": re.findall(r"^\s*module\s+(\w+)", source, re.M | re.I)})
    interfaces = []
    for path in sorted((ROOT / "addon/pycThermopack/thermopack").glob("*.py"), key=_path_sort_key):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in tree.body:
            if isinstance(node, ast.ClassDef):
                methods = [{"name": method.name, "line": method.lineno,
                            "arguments": ast.unparse(method.args), "status": "pending"}
                           for method in node.body
                           if isinstance(method, ast.FunctionDef)]
                interfaces.append({"path": path.relative_to(ROOT).as_posix(),
                                   "class": node.name, "methods": methods})
    return {"baseline": "d68c794", "fortran": files, "python": interfaces,
            "fluid_files": [p.name for p in sorted((ROOT / "fluids").glob("*.json"), key=_path_sort_key)],
            "binary_files": [p.name for p in sorted((ROOT / "binaries").glob("*.json"), key=_path_sort_key)],
            "unit_tests": [p.name for p in sorted((ROOT / "unittests").glob("*.pf"), key=_path_sort_key)]}


if __name__ == "__main__":
    result = inventory()
    target = ROOT / "docs/rust/inventory.json"
    target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"{len(result['fortran'])} Fortran files, "
          f"{sum(f['lines'] for f in result['fortran'])} lines, "
          f"{sum(len(c['methods']) for c in result['python'])} Python methods, "
          f"{len(result['fluid_files'])} fluid files")
