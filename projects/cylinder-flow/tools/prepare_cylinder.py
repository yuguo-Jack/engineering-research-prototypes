"""Lossless-input preprocessing for the supplied 2D case. Final solution is never read."""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import time
import h5py
import numpy as np


def prepare(case, data, output):
    started = time.perf_counter()
    with h5py.File(case, 'r') as f:
        mesh = f['meshes/1']
        points = mesh['nodes/coords/2'][()].astype('<f8')
        c0 = mesh['faces/c0/1'][()].astype('<i4')-1
        c1 = np.full(len(c0),-1,dtype='<i4')
        dataset = mesh['faces/c1/1']
        start = int(dataset.attrs['minId'][0])-1
        c1[start:start+len(dataset)] = dataset[()].astype('<i4')-1
        assert np.all(mesh['faces/nodes/1/nnodes'][()] == 2)
        nodes = mesh['faces/nodes/1/nodes'][()].reshape(-1,2).astype('<i4')-1
        topology = mesh['faces/zoneTopology']
        types = np.zeros(len(c0),dtype='<i4')
        for lo, hi, typ in zip(topology['minId'][()],topology['maxId'][()],topology['zoneType'][()]):
            types[int(lo)-1:int(hi)] = typ
        nc = int(mesh['cells/zoneTopology/maxId'][0])
    a,b = points[nodes[:,0]],points[nodes[:,1]]
    centers = (a+b)*0.5
    areas = np.column_stack((b[:,1]-a[:,1],a[:,0]-b[:,0]))
    cross = a[:,0]*b[:,1]-b[:,0]*a[:,1]
    volume = np.bincount(c0,weights=cross*.5,minlength=nc)
    moment = np.column_stack([np.bincount(c0,weights=(a[:,i]+b[:,i])*cross/6,minlength=nc) for i in range(2)])
    interior = c1>=0
    volume += np.bincount(c1[interior],weights=-cross[interior]*.5,minlength=nc)
    for i in range(2):
        moment[:,i] += np.bincount(c1[interior],weights=-(a[interior,i]+b[interior,i])*cross[interior]/6,minlength=nc)
    assert np.all(volume>0), 'face winding or cell volume error'
    cellcenters = moment/volume[:,None]
    closure = np.column_stack([np.bincount(c0,weights=areas[:,i],minlength=nc)-np.bincount(c1[interior],weights=areas[interior,i],minlength=nc) for i in range(2)])
    assert np.max(np.abs(closure))<1e-12
    orientation = np.einsum('ij,ij->i',areas[interior],cellcenters[c1[interior]]-cellcenters[c0[interior]])
    assert np.all(orientation>0)
    with h5py.File(data,'r') as f:
        cells = f['results/1/phase-1/cells']
        rho,u,v,p,t,h = [cells['SV_'+key+'/1'][()] for key in ['DENSITY','U','V','P','T','H']]
    R = float(np.mean(p/(rho*t)))
    cp = float(np.mean(h/t))
    gamma = cp/(cp-R)
    a0 = (gamma*R*float(t[0]))**0.5
    rho0 = float(rho[0]);p0=rho0*a0*a0
    initial = np.column_stack((rho/rho0,u/a0,v/a0,p/p0))
    output.mkdir(parents=True,exist_ok=True)
    with (output/'cylinder.bin').open('wb') as stream:
        stream.write(b'CFDCYL01')
        stream.write(struct.pack('<II',nc,len(c0)))
        stream.write(struct.pack('<7d',gamma,R,cp,rho0,a0,float(p[0]),float(t[0])))
        for array in [volume,cellcenters,c0,c1,types,centers,areas,initial]:
            stream.write(np.ascontiguousarray(array).tobytes())
    info={'case':str(case),'case_sha256':hashlib.sha256(case.read_bytes()).hexdigest(),'data':str(data),'data_sha256':hashlib.sha256(data.read_bytes()).hexdigest(),'cells':nc,'faces':len(c0),'nodes':len(points),'R':R,'cp':cp,'gamma':gamma,'density_scale':rho0,'velocity_scale':a0,'pressure_scale':p0,'initial_mach':float(u[0]/a0),'volume_total':float(volume.sum()),'minimum_volume':float(volume.min()),'geometry_closure_Linf':float(np.max(np.abs(closure))),'interior_normals_outward':bool(np.all(orientation>0)),'boundary_types':{str(int(k)):int(np.sum(types==k)) for k in np.unique(types)},'preprocessing_seconds':time.perf_counter()-started,'input_contract':'little-endian fixed-width integers/IEEE binary64; primitive rho,u,v,p scaled by initial rho and sound speed; original mesh/initial field only'}
    (output/'input-contract.json').write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(info,ensure_ascii=False,indent=2))


if __name__=='__main__':
    root=Path(__file__).resolve().parents[2]
    base=root/'2026-命题一_赛题资料/cylinder/case'
    parser=argparse.ArgumentParser()
    parser.add_argument('--case',type=Path,default=base/'cylinder_initial.cas.h5')
    parser.add_argument('--data',type=Path,default=base/'cylinder_initial.dat.h5')
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'data/cylinder')
    args=parser.parse_args();prepare(args.case,args.data,args.output)
