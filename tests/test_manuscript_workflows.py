"""Numerical/sample-selection regression checks for the consolidation."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import run_embedding_correlations as correlations
import run_embedding_replication as replication
import prepare_locus_data as locus




def test_regional_gene_track_uses_first_transcript_and_inclusive_overlap(tmp_path):
    gff=tmp_path/'gene.gff3'
    gff.write_text('Chr04\t.\tgene\t10\t30\t.\t+\t.\tID=g;Name=G\n'
        'Chr04\t.\tmRNA\t10\t30\t.\t+\t.\tID=t1;Parent=g\n'
        'Chr04\t.\tmRNA\t10\t30\t.\t+\t.\tID=t2;Parent=g\n'
        'Chr04\t.\tCDS\t12\t20\t.\t+\t0\tParent=t1\n'
        'Chr04\t.\tCDS\t15\t25\t.\t+\t0\tParent=t2\n')
    genes,exons=locus.gene_tracks(gff,'4',30,40)
    assert genes.gene_id.tolist()==['G']
    assert exons.seg_start.tolist()==[12]

def test_replication_uses_each_embeddings_own_strongest_marker(tmp_path):
    peaks=pd.DataFrame({'chrom':[4],'peak_start_bp':[100],'peak_end_bp':[199],
        'peak_window_Mb':[.0001],'max_sam3_embeddings':[2],'max_dino_embeddings':[0],
        'top_marker_pos':[110],'top_marker_p':[1e-12]})
    hits=pd.DataFrame({'CHROM':[4,4,4],'POS':[110,120,190],
        'trait':['embedding_mean_0','embedding_mean_0','embedding_std_1'],
        'p_value':[1e-9,1e-12,1e-11],'MARKER':['4:110:A:T','4:120:A:T','4:190:G:C'],'REF':['A','A','G'],'ALT':['T','T','C'],'effect':[.1,.2,.3]})
    a,b=tmp_path/'peaks.csv',tmp_path/'hits.csv';peaks.to_csv(a,index=False);hits.to_csv(b,index=False)
    result=replication.select_all_hotspot_embedding_pairs(a,b)
    assert len(result)==2
    # The column name is the historical representative_marker, but selection is per embedding.
    assert set(result.lead_marker)=={'4:120:A:T','4:190:G:C'}

@pytest.mark.parametrize('scope',['within','cross'])
def test_pairwise_outputs_match_reference_tables(scope,tmp_path,monkeypatch):
    """Check independently captured reference results and complete-case populations."""
    rng=np.random.default_rng(19);ids=[f'G{i}' for i in range(18)]
    hs=pd.DataFrame({'chrom':[4,9],'peak_start_bp':[100,100],'peak_end_bp':[200,200],
        'peak_marker':['4:150:A:T','9:150:G:C']})
    hs.to_csv(tmp_path/'hotspots.csv',index=False)
    pd.DataFrame({'CHROM':[4,4,9],'POS':[150,155,150],
        'trait':['embedding_mean_0','embedding_std_1','embedding_std_2']}).to_csv(tmp_path/'hits.csv',index=False)
    pd.DataFrame(columns=['CHROM','POS','trait']).to_csv(tmp_path/'empty_hits.csv',index=False)
    frame=pd.DataFrame(rng.normal(size=(18,3)),columns=['embedding_mean_0','embedding_std_1','embedding_std_2'])
    frame.loc[0,'embedding_std_2']=np.nan  # within chr4 retains this sample; cross drops it
    frame.insert(0,'genotype',ids);frame.to_csv(tmp_path/'blue.csv',index=False)
    pd.DataFrame({'genotype':ids}).to_csv(tmp_path/'empty_blue.csv',index=False)
    cov=pd.DataFrame(rng.normal(size=(18,7)),columns=[*[f'PC{i}' for i in range(1,6)],'human_score','ExG_P20_disease_pct'])
    cov.insert(0,'genotype',ids);cov.to_csv(tmp_path/'cov.csv',index=False)
    settings=dict(HOTSPOT_MASTER=tmp_path/'hotspots.csv',SIGNIFICANT_MARKERS={'sam3':tmp_path/'hits.csv','dino2':tmp_path/'empty_hits.csv'},
        EMBEDDING_BLUES={'sam3':tmp_path/'blue.csv','dino2':tmp_path/'empty_blue.csv'},
        PC_FILE=tmp_path/'cov.csv',HUMAN_SCORE_FILE=tmp_path/'cov.csv',EXG_LOGIT_FILE=tmp_path/'cov.csv')
    dose=pd.DataFrame({m:np.arange(18)%3 for m in hs.peak_marker},index=ids)
    for module in [correlations]:
        for k,v in settings.items():monkeypatch.setattr(module,k,v)
        if scope=='within':monkeypatch.setattr(module,'load_peak_marker_dosages',lambda _:dose)
    monkeypatch.setattr(correlations,'OUT_CSV',tmp_path/'new.csv')
    (correlations.run_within if scope=='within' else correlations.run_cross)()
    pd.testing.assert_frame_equal(pd.read_csv(ROOT/'tests/fixtures'/f'partial_correlations_{scope}.csv'),pd.read_csv(tmp_path/'new.csv'),check_exact=True)



def test_hotspots_count_distinct_traits_merge_bins_and_resolve_ties(tmp_path):
    # Two equally dense bins 200 kb apart form one peak; a 300-kb gap starts another.
    rows = []
    for start, pvalue in [(100_000, 1e-10), (300_000, 1e-12), (600_000, 1e-11)]:
        rows.extend(dict(CHROM=4, POS=start+10, trait=f'e{i}', p_value=pvalue) for i in range(10))
    rows.append(rows[0].copy())  # repeated association must not inflate feature count
    sam3 = tmp_path/'sam3.csv';dino2 = tmp_path/'dino2.csv'
    pd.DataFrame(rows).to_csv(sam3,index=False)
    pd.DataFrame(columns=['CHROM','POS','trait','p_value']).to_csv(dino2,index=False)
    replication.generate_hotspot_tables(sam3,dino2,tmp_path/'out')
    result = pd.read_csv(tmp_path/'out/sam3_peaks_ge10_embeddings.csv')
    assert result.peak_start_bp.tolist() == [100_000,600_000]
    assert result.peak_end_bp.tolist() == [399_999,699_999]
    assert result.peak_window_Mb.tolist() == [.3,.6]
    assert result.max_sam3_embeddings.tolist() == [10,10]
    assert result.max_dino_embeddings.tolist() == [0,0]
    assert result.top_marker_pos.tolist() == [300_010,600_010]
    assert pd.read_csv(tmp_path/'out/dino2_peaks_ge10_embeddings.csv').empty


def test_common_cohort_uses_all_three_fitted_populations(tmp_path):
    for env, ids in zip(replication.ENVIRONMENTS, [['SC 1','G2','G3'],['SC1','G3'],['SC1','G2']]):
        pd.DataFrame({'genotype':ids}).to_csv(tmp_path/f'blues_{env}.csv',index=False)
    output = tmp_path/'out/common.csv'
    replication.generate_common_genotypes(tmp_path,output)
    assert pd.read_csv(output).genotype.tolist() == ['SC1']
