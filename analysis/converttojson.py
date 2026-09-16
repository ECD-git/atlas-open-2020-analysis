import json
import importlib.util
import pathlib

spec = importlib.util.spec_from_file_location("mcinfofile", str(pathlib.Path(__file__).parent.resolve())+"/mcinfofile.py")
mcinfofile = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mcinfofile)

with open(str(pathlib.Path(__file__).parent.resolve())+"/mcinfofile.json", "w") as f:
    json.dump(mcinfofile.infos, f, indent=2)

print(f"Wrote {len(mcinfofile.infos)} entries to mc_info.json")