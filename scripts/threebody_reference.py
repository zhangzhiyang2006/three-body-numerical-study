"""Independent pair-force/DOP853 references. Run: uv run --with scipy python scripts/threebody_reference.py"""
import hashlib
import json
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/threebody_reference'
CASES={
 'single':([1.],[0.,2.,1.,0.],50.,500),
 'two':([.3,.03],[2.,2.,0.,0.,.2,-.2,-.01,.01],20.,200),
 'three':([.3,.03,.03],[2.,2.,0.,0.,-2.,-2.,.2,-.2,0.,0.,-.2,.2],30.,300),
 'figure8':([1.,1.,1.],[-.97000436,.24308753,0.,0.,.97000436,-.24308753,.466203685,.432365730,-.932407370,-.864731460,.466203685,.432365730],20*6.32591398,12000),
}
def force(masses):
 m=np.array(masses);n=len(m)
 def rhs(t,z):
  q=z[:2*n].reshape(n,2);a=np.zeros_like(q)
  if n==1:a[0]=-q[0]/np.linalg.norm(q[0])**3 # central mass M=1, G=1
  else:
   for i in range(n):
    for j in range(i+1,n):
     d=q[j]-q[i];f=d/np.linalg.norm(d)**3
     a[i]+=m[j]*f;a[j]-=m[i]*f
  return np.r_[z[2*n:],a.ravel()]
 return rhs

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 meta={'method':'DOP853','scipy':scipy.__version__,'numpy':np.__version__,'G':1.,'rtol':3e-14,'atol':3e-16,'comparison_rtol':1e-12,'comparison_atol':1e-14,'cases':{},'files':{}}
 for name,(m,z,T,n) in CASES.items():
  t=np.linspace(0,T,n+1);f=force(m)
  fine=solve_ivp(f,(0,T),z,method='DOP853',rtol=3e-14,atol=3e-16,t_eval=t)
  coarse=solve_ivp(f,(0,T),z,method='DOP853',rtol=1e-12,atol=1e-14,t_eval=t)
  assert fine.success and coarse.success and np.isfinite(fine.y).all()
  path=OUT/f'{name}.csv';np.savetxt(path,np.c_[t,fine.y.T],delimiter=',')
  meta['cases'][name]={'masses':m,'initial_state':z,'end_time':T,'samples':n+1,'reference_tolerance_difference':float(np.max(np.abs(fine.y-coarse.y))),'nfev':fine.nfev}
  if name=='single':
   D=2*np.sinh(np.arcsinh(3*t/8)/3);exact=np.array([4*D,2-2*D**2,1/(1+D**2),-D/(1+D**2)])
   meta['cases'][name]['analytic_error']=float(np.max(np.abs(fine.y-exact)))
  if name in ('two','three'):
   pert=np.array(z);pert[7 if name=='two' else 10]+=1e-6
   p=solve_ivp(f,(0,T),pert,method='DOP853',rtol=3e-14,atol=3e-16,t_eval=t)
   assert p.success
   np.savetxt(OUT/f'{name}_pert.csv',np.c_[t,p.y.T],delimiter=',')
  print(name,meta['cases'][name]['reference_tolerance_difference'],flush=True)
 for path in OUT.glob('*.csv'):meta['files'][path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
 (OUT/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
 (OUT/'tolerance_checks.csv').write_text(''.join(f"{name},{case['reference_tolerance_difference']:.17e}\n" for name,case in meta['cases'].items()))
if __name__=='__main__':main()
