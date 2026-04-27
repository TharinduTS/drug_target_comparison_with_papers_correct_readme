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
