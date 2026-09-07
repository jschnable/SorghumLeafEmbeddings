#!/usr/bin/env python3
"""Generate regional GWAS, all-panel LD and representative-transcript gene tracks.

Replaces six copied locus scripts. Regions/features are frozen in locus_regions.json.
Kinship is recomputed for the exact phenotype sample set instead of trusting the
historical cache based only on its number of samples. Existing figure data are
untouched unless this workflow is explicitly run for that locus.
"""
from __future__ import annotations
import argparse
import json
import re
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from figure_data_io import save_region_gwas
from marker_utils import DEFAULT_GENOTYPE, DEFAULT_COVARIATE_FILE
ROOT=Path(__file__).resolve().parents[1]

def gene_tracks(gff,chrom,start,end):
    lines=[]
    with open(gff) as f:
        for line in f:
            if line.startswith('#'):continue
            p=line.rstrip('\n').split('\t')
            if len(p)==9 and p[0]==f'Chr{int(chrom):02d}':lines.append(p)
    def attr(s,k):
        m=re.search(rf'{k}=([^;]+)',s);return m.group(1) if m else None
    genes=[];ids={};rep={};seen=set()
    for p in lines:
        s,e=int(p[3]),int(p[4])
        if p[2]=='gene' and not(e<start or s>end):
            name=attr(p[8],'Name');genes.append((name,s,e,p[6]));ids[attr(p[8],'ID')]=name
    for p in lines:
        par=attr(p[8],'Parent')
        if p[2]=='mRNA' and par in ids and par not in seen:
            seen.add(par);rep[attr(p[8],'ID')]=ids[par]
    exons=[]
    for p in lines:
        par=attr(p[8],'Parent')
        if p[2] in ['CDS','five_prime_UTR','three_prime_UTR'] and par in rep:
            exons.append((rep[par],p[6],int(p[3]),int(p[4]),'CDS' if p[2]=='CDS' else 'UTR'))
    return (pd.DataFrame(genes,columns=['gene_id','start','end','strand']),
        pd.DataFrame(exons,columns=['gene_id','strand','seg_start','seg_end','kind']))

def ld_track(vcf,chrom,start,end,lead):
    gm={'0/0':0,'0|0':0,'0/1':1,'0|1':1,'1|0':1,'1/1':2,'1|1':2}
    q=subprocess.run(['bcftools','query','-r',f'{chrom}:{start}-{end}','-f','%POS[\t%GT]\n',str(vcf)],capture_output=True,text=True,check=True)
    pos=[];mat=[]
    for line in q.stdout.splitlines():
        f=line.split('\t');d=np.array([gm.get(x,np.nan) for x in f[1:]],float)
        if np.isfinite(d).sum()<50:continue
        pos.append(int(f[0]));mat.append(d)
    if lead not in pos:raise ValueError(f'Lead {lead} missing from LD population')
    m=np.vstack(mat);m=np.where(np.isnan(m),np.nanmean(m,axis=1,keepdims=True),m)
    z=m-m.mean(axis=1,keepdims=True);sd=z.std(axis=1,keepdims=True);sd[sd==0]=np.nan;z/=sd
    r2=((z@z[pos.index(lead)])/m.shape[1])**2
    return pd.DataFrame({'POS':pos,'r2':r2}).dropna()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,default=ROOT/'data/provided/locus_regions.json')
    p.add_argument('--locus',required=True)
    p.add_argument('--genotype',type=Path,default=DEFAULT_GENOTYPE)
    p.add_argument('--gff',type=Path,required=True,help='BTx623 v5.1 gene GFF3')
    p.add_argument('--blues',type=Path,default=ROOT/'data/generatable/blues/nebraska_sam3_embeddings_2016crop/blues_Nebraska2025.csv')
    p.add_argument('--covariates',type=Path,default=DEFAULT_COVARIATE_FILE)
    p.add_argument('--out-dir',type=Path)
    p.add_argument('--tracks-only',action='store_true',help='Prepare gene and LD tracks without regional association fitting')
    p.add_argument('--cpu',type=int,default=1)
    a=p.parse_args();cfg=json.loads(a.config.read_text())[a.locus]
    out=a.out_dir or ROOT/cfg['output'];out.mkdir(parents=True,exist_ok=True)
    genes,exons=gene_tracks(a.gff,cfg['chrom'],cfg['start'],cfg['end'])
    genes.to_csv(out/'gene_models.csv',index=False);exons.to_csv(out/'gene_exons.csv',index=False)
    ld_track(a.genotype,cfg['chrom'],cfg['start'],cfg['end'],cfg['lead']).to_csv(out/'ld_track.csv',index=False)
    if not a.tracks_only:
        from panicle.data.loaders import load_genotype_file
        from panicle.matrix.pca import PANICLE_PCA
        from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
        from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI
        geno,ids,gmap=load_genotype_file(str(a.genotype),file_format='vcf',precompute_alleles=False)
        mdf=gmap.to_dataframe();mdf.CHROM=mdf.CHROM.astype(str)
        region=np.flatnonzero((mdf.CHROM==cfg['chrom']) & mdf.POS.between(cfg['start'],cfg['end']))
        blues=pd.read_csv(a.blues).set_index('genotype');cov=pd.read_csv(a.covariates).set_index('genotype')
        selected=np.array([i for i,g in enumerate(ids) if g in blues.index and g in cov.index and np.isfinite(cov.loc[g].values).all()])
        samples=[ids[i] for i in selected];sub=geno.subset_individuals(selected)
        pcs=PANICLE_PCA(M=sub,pcs_keep=5,verbose=False)
        c=cov.loc[samples,['mask_pixels_blue','days_to_flower_blue']].to_numpy(float)
        c=(c-c.mean(axis=0))/c.std(axis=0)
        loco=PANICLE_K_VanRaden_LOCO(sub,gmap,maxLine=5000,cpu=a.cpu,verbose=False)
        res=PANICLE_MLM_LOCO_MULTI(phe=blues.loc[samples,cfg['traits']].to_numpy(float),geno=sub.subset_markers(region),
            map_data=gmap.subset_markers(region),trait_names=cfg['traits'],loco_kinship=loco,CV=np.column_stack([pcs,c]),
            maxLine=5000,cpu=a.cpu,lrt_refinement=True,verbose=False)
        rows=[(t,int(pos),float(pv)) for t in cfg['traits'] for pos,pv in zip(mdf.POS.iloc[region],np.asarray(res[t].pvalues,float))]
        save_region_gwas(out/'region_gwas.npz',pd.DataFrame(rows,columns=['trait','POS','p_value']))
    me=cfg['effective_tests']
    meta=dict(cfg.get("metadata", {}))
    meta.update(dict(Me=me,bonferroni_threshold=.05/me,neglog10_threshold=float(np.log10(me/.05)),
        region_chrom=cfg['chrom'],region_lo=cfg['start'],region_hi=cfg['end'],peak_marker=cfg['lead']))
    (out/'meta.json').write_text(json.dumps(meta,indent=2)+'\n')

if __name__=='__main__':main()
