"""Run the standalone three-body tests against maintained notebook definitions."""
import json,subprocess,tempfile,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
n=json.loads((ROOT/'三体问题数值求解.ipynb').read_text())
with tempfile.TemporaryDirectory() as td:
    p=Path(td)/'runtests.jl'
    source='using Test, LinearAlgebra, Printf, DelimitedFiles\n'
    source+='\n'.join(''.join(n['cells'][i]['source']) for i in [4,6])
    source+='\n'+(ROOT/'src/threebody/validation.jl').read_text()
    source+='\n'+(ROOT/'tests/threebody_validation.jl').read_text()
    p.write_text(source)
    subprocess.run([shutil.which('julia') or 'julia','--startup-file=no',str(p)],cwd=ROOT,check=True)
