## Batch Effect Removal from RNA-Seq data
This repo contains methods for removing batch-effects from RNA-Seq data.
At this point, the ComBat algorithm is the default pre-processing step for most downstream analyisis such as drug response prediction and tissue type classification.

The PCA plots demonstrate the effect of ComBat method applied to RNA-Seq profiles. 
The RNA-Seq comes from multiple data sources (i.e., drug sensitivity studies).
PCA of raw expression data shows that the source carries the strognest component in the combined dataset.
ComBat normalization significantly suppresses the signal that's coming from source component.
This often improves predictive algorithms that aim to extract important biological information.

<!---
Data Before Dec 2019
<p float="left">
  <img src="README/Dec2019/pca_raw.png" width="300" height="300">
  <img src="README/Dec2019/pca_src_scale.png" width="300" height="300">
  <img src="README/Dec2019/pca_combat.png" width="300" height="300">
</p>
--->

<!--- Data Since Dec 2020 --->
<p float="left">
  <img src="README/July2020/pca_raw.png" width="300" height="250">
  <img src="README/July2020/pca_src_scale.png" width="300" height="250">
  <img src="README/July2020/pca_combat.png" width="300" height="250">
</p>

## Example
Get the required files from nciftp:
- combined_rnaseq_data_lincs1000
- cl_mapping
```shell
bash scripts/run_Dec2019.bash
```

## Combat
We adapted the Python implementation of combat from this repo https://github.com/brentp/combat.py. Thanks!
