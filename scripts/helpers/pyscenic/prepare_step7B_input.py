import sys
from pathlib import Path
import numpy as np
import pandas as pd
import loompy
from scipy import sparse
from scipy.io import mmread
from ctxcore.rnkdb import FeatherRankingDatabase
mtx,gf,bf,mf,tf,ff,loom,smoke,smoketf,audit = map(Path, sys.argv[1:11])
force = sys.argv[11].lower() == 'true'
eg=int(sys.argv[12]); ec=int(sys.argv[13])
def onecol(p):
    return [x.strip().split('\t')[0] for x in p.open(encoding='utf-8-sig') if x.strip()]
genes=onecol(gf); cells=onecol(bf); tfs=onecol(tf); md=pd.read_csv(mf, low_memory=False)
req=['cell','step5CR_state','sample_id','timepoint','retained_step5BR_pseudotime','Y29_signed_progress']
miss=[x for x in req if x not in md.columns]
assert not miss, f'Missing metadata: {miss}'
assert len(genes)==eg and len(cells)==ec
assert len(set(genes))==len(genes) and len(set(cells))==len(cells)
assert md.shape[0]==len(cells) and md['cell'].astype(str).tolist()==cells
states=md['step5CR_state'].astype(str).value_counts().to_dict()
assert states=={'Shared_prebranch':116,'Arm_1':2289,'Arm_2':1798}, states
x=mmread(str(mtx)); x=sparse.csr_matrix(x, dtype=np.float32)
assert x.shape==(len(genes),len(cells)) and x.nnz>0
assert np.all(np.isfinite(x.data)) and np.all(x.data>=0)
tfset=set(tfs); overlap=[g for g in genes if g in tfset]
db=FeatherRankingDatabase(fname=str(ff),name='mm10_10kb'); dbg=set(db.genes); dbo=[g for g in genes if g in dbg]
rebuild=force or not loom.exists()
if not rebuild:
    try:
        with loompy.connect(str(loom),'r') as ds:
            rebuild=tuple(ds.shape)!=tuple(x.shape) or list(ds.ra.Gene.astype(str))!=genes or list(ds.ca.CellID.astype(str))!=cells
    except Exception:
        rebuild=True
if rebuild:
    if loom.exists(): loom.unlink()
    ca={'CellID':np.asarray(cells,dtype=object),'BranchState':np.asarray(md['step5CR_state'].astype(str),dtype=object),
        'SampleID':np.asarray(md['sample_id'].astype(str),dtype=object),'Timepoint':np.asarray(md['timepoint'].astype(str),dtype=object),
        'Pseudotime':pd.to_numeric(md['retained_step5BR_pseudotime'],errors='coerce').to_numpy(dtype=float),
        'Y29SignedProgress':pd.to_numeric(md['Y29_signed_progress'],errors='coerce').to_numpy(dtype=float)}
    loompy.create(str(loom),x,{'Gene':np.asarray(genes,dtype=object)},ca)
# smoke subset: all prebranch + 100 cells from each arm; 1000 variable genes + 60 variable TFs
idx=[]
for s,n in [('Shared_prebranch',116),('Arm_1',100),('Arm_2',100)]: idx += np.flatnonzero(md['step5CR_state'].astype(str).to_numpy()==s)[:n].tolist()
idx=np.asarray(sorted(set(idx)),dtype=int)
mu=np.asarray(x.mean(axis=1)).ravel(); va=np.maximum(np.asarray(x.power(2).mean(axis=1)).ravel()-mu**2,0)
top=np.argsort(va)[-1000:]; tfi=np.asarray([i for i,g in enumerate(genes) if g in tfset],dtype=int); stfi=tfi[np.argsort(va[tfi])[-60:]]
gi=np.asarray(sorted(set(top.tolist()+stfi.tolist())),dtype=int); sx=x[gi,:][:,idx].tocsr(); sg=[genes[i] for i in gi]; sc=[cells[i] for i in idx]; st=[g for g in sg if g in tfset]
assert len(st)>=20
if smoke.exists(): smoke.unlink()
loompy.create(str(smoke),sx,{'Gene':np.asarray(sg,dtype=object)},{'CellID':np.asarray(sc,dtype=object)})
smoketf.write_text('\n'.join(st)+'\n',encoding='utf-8')
metrics=[('matrix_genes',len(genes)),('matrix_cells',len(cells)),('matrix_nnz',x.nnz),('TF_overlap',len(overlap)),
('feather_overlap',len(dbo)),('feather_overlap_fraction',len(dbo)/len(genes)),('loom_rebuilt',rebuild),
('smoke_genes',len(sg)),('smoke_cells',len(sc)),('smoke_TFs',len(st)),
('state_Shared_prebranch',states['Shared_prebranch']),('state_Arm_1',states['Arm_1']),('state_Arm_2',states['Arm_2'])]
with audit.open('w',encoding='utf-8') as h:
    h.write('metric\tvalue\n')
    [h.write(f'{k}\t{v}\n') for k,v in metrics]
print('STEP7B_PREFLIGHT_PASS'); [print(k,v) for k,v in metrics]
