import pandas as pd
import scanpy as sc
import drugreflector as dr
from utils import compute_vscores_adata

adata = sc.datasets.pbmc3k()   
annots = sc.datasets.pbmc3k_processed().obs

adata.obs = pd.merge(adata.obs, annots, how='left', left_index=True, right_index=True)

adata

vscores = compute_vscores_adata(
    adata, 
    group_col='louvain',
    group1_value='CD14+ Monocytes',    # Classical monocytes
    group2_value='FCGR3A+ Monocytes', 
    layer = None   # Non-classical monocytes
)
print(f"V-score comparison: {vscores.name}")
print(f"Computed v-scores for {len(vscores)} genes")
print(f"Top upregulated genes in FCGR3A+ vs CD14+ monocytes:")
print(vscores.nlargest(10))


model_paths = [
    'checkpoints/model_fold_0.pt',
    'checkpoints/model_fold_1.pt', 
    'checkpoints/model_fold_2.pt'
]
model = dr.DrugReflector(checkpoint_paths=model_paths)


predictions = model.predict(vscores, n_top=50)
print(predictions)

print(f"Prediction results shape: {predictions.shape}")
print(f"Columns: {predictions.columns.names}")

Drugs = pd.read_csv('../compoundinfo_beta.txt', sep= '\t')

predictions = predictions.reset_index()
predictions.columns = ["pert_id", "rank", "logit", "prob"]
print(predictions)
merged = pd.merge(predictions, Drugs, 
                  on= 'pert_id')
print(merged)


NK_data = sc.read_h5ad("../NK.h5ad")
NK_data
vscores_NK = compute_vscores_adata(
    NK_data, 
    group_col='group',
    group1_value='NDR',    
    group2_value='DR', 
    layer = None   
)

predictions = model.predict(vscores_NK, n_top=50)
predictions

merged.to_csv('../drug.csv', sep='\t', index= False)