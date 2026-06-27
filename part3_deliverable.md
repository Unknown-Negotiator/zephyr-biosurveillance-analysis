# Part 3B - Embedding-Based Clustering

For the optional coding track, I chose embedding-based clustering. I used a lightweight nucleotide-composition approach rather than downloading a large sequence model. Each reference sequence and each sample-virus read set was represented as a normalized 5-mer frequency vector. I then centered the matrix, projected it into six dimensions with SVD, and clustered the embedded points with k-means using seven clusters. Sample-virus read sets were included when the Part 2 classifier assigned at least 25 reads to that virus in that sample.

This analysis embedded 38 respiratory-virus reference records and 22 sample-virus read sets. The reference panel included common respiratory-virus genomes or segments across Coronaviridae, Picornaviridae, Orthomyxoviridae, Pneumoviridae, Paramyxoviridae, Adenoviridae, and Parvoviridae. The sample read sets came from the Part 2 calls, mainly rhinoviruses plus one Influenza C and one HCoV OC43 high-confidence signal.

## Results

The clustering agreed well with the broad taxonomy profile. The weighted family purity across clusters was 0.917. For sample read sets, the nearest-reference family match rate was 1.000: every sample read set had its nearest reference in the expected viral family. The exact-virus match rate was lower, 0.455, mostly because many rhinovirus A/B/C sample read sets were nearest to another rhinovirus reference rather than the exact assigned rhinovirus species.

The main sample-containing clusters were:

| Cluster | Records | References | Sample read sets | Main family | Family purity | Interpretation |
|---:|---:|---:|---:|---|---:|---|
| 2 | 24 | 4 | 20 | Picornaviridae | 95.8% | Contains all rhinovirus sample read sets and the rhinovirus references; this supports the Part 2 rhinovirus calls at family/genus level. |
| 4 | 9 | 8 | 1 | Orthomyxoviridae | 100.0% | Contains the Influenza C sample read set from 260608-Copl-NAS-P1 near Influenza C references. |
| 6 | 6 | 5 | 1 | Coronaviridae | 100.0% | Contains the HCoV OC43 sample read set from 260609-Cent-NAS-P2 near coronavirus references. |

The remaining clusters contained references only. Influenza references split across several clusters, which is expected because the reference panel contains multiple segmented influenza records rather than one genome-length record per virus. Pneumoviridae and Paramyxoviridae grouped together in one reference-only cluster, suggesting that this simple 5-mer/SVD approach captures broad sequence-composition structure but does not perfectly separate every respiratory-virus family.

## Comparison to Taxonomic Profile

The embedding results support the main Part 2 conclusions. The high-confidence Influenza C call in 260608-Copl-NAS-P1 clustered with Orthomyxoviridae and had Influenza C as its nearest reference. The high-confidence HCoV OC43 call in 260609-Cent-NAS-P2 clustered with Coronaviridae and had HCoV OC43 as its nearest reference. The rhinovirus calls all clustered in the Picornaviridae/rhinovirus region, consistent with the strong read-count and coverage evidence in Part 2.

Exact rhinovirus species agreement was weaker. Many sample read sets assigned as Rhinovirus A had Rhinovirus B or C as the nearest reference in 5-mer space, and some Rhinovirus B calls were nearest to Rhinovirus C. I interpret this as a limitation of the simple composition embedding, not as evidence against rhinovirus presence. The Part 2 mappings already showed that many rhinovirus calls had broad genome coverage but relatively low identity to prototype references, so family/genus-level agreement is the more appropriate readout here.

## Limitations

This is a deliberately lightweight clustering analysis. Normalized 5-mer composition ignores alignment, gene order, and evolutionary substitution models, and it can be affected by read length, coverage bias, and ONT error profiles. The reference panel is also small and uses prototype references, which limits exact species resolution for diverse rhinoviruses. A stronger follow-up would compare assembled contigs or higher-confidence ORFs against a broader reference set and add ANI or phylogenetic analysis. For this work task, the embedding provides an independent sanity check that the sample read sets occupy the expected viral-family neighborhoods.

The full outputs are in `results/part3_combined_kmer_embeddings.csv`, `results/part3_cluster_summary.csv`, and `results/part3_sample_reference_comparison.csv`. The main figure is `figures/part3_sample_reference_kmer_clusters.png`.
