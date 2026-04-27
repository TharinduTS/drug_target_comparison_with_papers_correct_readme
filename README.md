# drug_target_comparison_with_papers_correct_readme

Here I am trying to do my drug target identification prioritization in a similar way to previous papers so I can easily compare these

Following Ryaboshapkina et al

With my previous tests I realized using enrichment values within tissue types is the better measure to use for drug target identification.
There fore I am starting from the place I calculated enrichment values for enrichment values within tissue types - 

https://github.com/TharinduTS/cell_type_enrichment_v2#11-i-tissue-expression-layer

Here I have the enrichment values in a table that looks like

```txt
Gene    Gene name       Tissue:CellType avg_nCPM        weight_sum      clusters_used   specificity_tau Enrichment score        log2_enrichment Enrichment score (tau penalized)      log2_enrichment_penalized       single_cell_type_gene
ENSG00000167531 LALBA   breast:breast lactating cells   483428.6        420     1       0.9999997849673347      4606146.319293296       22.135128809974248      4606145.328821376     22.135128499747655      False
ENSG00000135222 CSN2    breast:breast lactating cells   333043.0        420     1       0.9999991703984089      1201050.8623900507      20.195865817257094      1201049.8659963442    20.1958646203945        False
```

I started by extracting only the essential columns

```bash
less enrichment_values_for_tissue_types.tsv | cut -f 1,2,3,11 >selected_data_for_comparison.tsv
```

and that gave me 
```txt
Gene    Gene name       Tissue:CellType log2_enrichment_penalized
ENSG00000167531 LALBA   breast:breast lactating cells   22.135128499747655
ENSG00000135222 CSN2    breast:breast lactating cells   20.1958646203945
```

Then I selected the highest and runners up values for each gene and calculated the ratio between them with python

When the runners up has a negative value, the ratio becomes 'not_enough_secondary_expression'

cal_ratio_all_genes.py
```
import pandas as pd
import numpy as np

df = pd.read_csv("selected_data_for_comparison.tsv", sep="\t")

df_sorted = df.sort_values(
    ["Gene", "Gene name", "log2_enrichment_penalized"],
    ascending=[True, True, False]
)

top2 = (
    df_sorted
    .groupby(["Gene", "Gene name"])
    .head(2)
    .copy()
)

top2["rank"] = (
    top2
    .groupby(["Gene", "Gene name"])
    .cumcount() + 1
)

highest = top2[top2["rank"] == 1].set_index(["Gene", "Gene name"])
runner_up = top2[top2["rank"] == 2].set_index(["Gene", "Gene name"])

result = pd.concat(
    [
        highest["Tissue:CellType"]
            .rename("highest_Tissue:CellType"),
        highest["log2_enrichment_penalized"]
            .rename("highest_log2_enrichment_penalized"),
        runner_up["Tissue:CellType"]
            .rename("runners_up_Tissue:CellType"),
        runner_up["log2_enrichment_penalized"]
            .rename("runners_up_log2_enrichment_penalized"),
    ],
    axis=1
)

# ✅ Conditional ratio logic
result["ratio"] = np.where(
    result["runners_up_log2_enrichment_penalized"] < 0,
    "not_enough_secondary_expression",
    result["highest_log2_enrichment_penalized"] /
    result["runners_up_log2_enrichment_penalized"]
)

result = result.reset_index()

result.to_csv("gene_top2_enrichment.tsv", sep="\t", index=False)
```

The output looks like following 

```
Gene    Gene name       highest_Tissue:CellType highest_log2_enrichment_penalized       runners_up_Tissue:CellType      runners_up_log2_enrichment_penalized    ratio
ENSG00000000003 TSPAN6  testis:late spermatids  7.093081007217921       esophagus:esophageal apical cells       6.157569660424798       1.1519286664031965
ENSG00000000005 TNMD    breast:vascular endothelial cells       6.877755591875127       adipose tissue:adipocytes       4.802532933968559       1.4321100316102808
```
