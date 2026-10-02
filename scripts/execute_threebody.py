"""Execute and save the lesson with the installed Julia and this project's environment."""
import json,os,shutil,tempfile
from pathlib import Path
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager
ROOT=Path(__file__).resolve().parents[1]
os.environ['GKSwstype']='100'
path=ROOT/'三体问题数值求解.ipynb'
nb=nbformat.read(path,as_version=4)
with tempfile.TemporaryDirectory() as td:
    spec=Path(td)/'threebody';spec.mkdir()
    (spec/'kernel.json').write_text(json.dumps({'argv':[shutil.which('julia'),'-i','--startup-file=no','--project='+str(ROOT),'-e','import IJulia; IJulia.run_kernel()','{connection_file}'],'display_name':'Julia project','language':'julia'}))
    km=KernelManager(kernel_name='threebody',kernel_spec_manager=KernelSpecManager(kernel_dirs=[td]))
    def done(cell,cell_index,**kwargs):
        print('PASS cell',cell_index,flush=True)
        for output in cell.get('outputs',[]):
            if output.output_type=='stream':print(output.text,flush=True)
    client=NotebookClient(nb,km=km,timeout=300,resources={'metadata':{'path':str(ROOT)}},on_cell_executed=done)
    try:
        client.execute()
        nbformat.write(nb,path)
    finally:
        if km.has_kernel:km.shutdown_kernel(now=True)
print('Saved executed notebook:',path,flush=True)
