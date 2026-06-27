# Pivotal Biodefense Work Task: Metagenomic Biosurveillance

**Time spent:** [fill in actual time]  
**Code:** [paste GitHub repository link]

## Part 1 - Evaluation Design for AI-Generated Nucleic Acids

To evaluate a classifier's ability to detect AI-generated nucleic acids, I would use a balanced held-out benchmark containing both authentic biological sequences and red-teamed synthetic sequences. Natural sequences should span viruses, bacteria, and eukaryotes, and should be stratified by length, GC content, codon usage, taxonomy, and genomic context. Synthetic sequences should come from multiple generative tools, including adversarial "synthetic homologs" designed to preserve function while reducing detectable similarity to known sequences [1]. I would report accuracy, precision, recall, false-negative rate, and ROC-AUC, with particular emphasis on minimizing false negatives. Stress tests should include obfuscations such as decoy flanking regions, fragmented hazards, synonymous recoding, and sequences designed to evade alignment-based screening [2]. Performance should be compared with simple baselines and published alignment-free methods such as Synsor [3]. Finally, I would bootstrap confidence intervals and report subgroup performance to identify where detection fails.

## Methods

For Part 2, I analyzed 12 Zephyr respiratory-read pools from June 5-13, 2026. I mapped ONT respiratory viral reads against a curated panel of common respiratory-virus references using minimap2, kept the best viral alignment per read, and summarized read support, median alignment identity, genome coverage breadth, and mean depth. I assigned confidence labels conservatively: high confidence required at least 50 reads, at least 20% reference coverage, and median identity at least 85%; medium required at least 10 reads, at least 2% coverage, and median identity at least 75%; low required at least 3 reads and median identity at least 70%.

For Part 3B, I used a lightweight embedding-based clustering method. Each reference sequence and each sample-virus read set with at least 25 assigned reads was represented by normalized 5-mer composition. I projected these vectors into six dimensions with SVD and clustered them with k-means into seven clusters. This avoided downloading a large DNA language model while still providing an independent check on whether sample-derived reads grouped near the expected viral families.

## Results

The dominant signal across the pools was rhinovirus. Rhinovirus A and B were detected in all 12 pools, often with broad genome coverage. Rhinovirus C was detected in 6 pools with weaker support. I also found two clear non-rhinovirus respiratory-virus signals: Human coronavirus OC43 in sample 260609-Cent-NAS-P2 and Influenza C in sample 260608-Copl-NAS-P1. No SARS-CoV-2, influenza A/B, RSV A/B, metapneumovirus, parainfluenza viruses, adenovirus C, or bocavirus 1 were called under this reference-panel approach.

| Virus | Positive pools | Total assigned reads | Max reads in one pool | Median identity | Max coverage breadth |
|---|---:|---:|---:|---:|---:|
| Rhinovirus A | 12 | 170,550 | 89,662 | 72.5% | 100.0% |
| Rhinovirus B | 12 | 62,681 | 30,881 | 74.8% | 100.0% |
| HCoV OC43 | 6 | 20,394 | 20,386 | 97.2% | 88.6% |
| Rhinovirus C | 6 | 677 | 329 | 73.2% | 63.5% |
| Influenza C | 1 | 413 | 413 | 96.0% | 85.8% |

**Figure 1.** Genome coverage for the high-confidence HCoV OC43 call in sample 260609-Cent-NAS-P2. Reads covered 27,247 / 30,741 reference bases (88.6% breadth), with 20,386 assigned reads and high median identity, supporting a robust OC43 detection in this pool.  
**Insert:** `figures/coverage/260609-Cent-NAS-P2__Human_coronavirus_OC43.png`

The rhinovirus calls should be interpreted carefully. Many had thousands of reads and near-complete genome coverage, but only about 72-77% identity to prototype references. I therefore treat them as strong rhinovirus-family/genus signals, but more cautious exact species assignments.

For Part 3B, the k-mer embedding clustered sample read sets with the expected viral families. The analysis embedded 38 reference records and 22 sample-virus read sets. Weighted family purity across clusters was 0.917. For sample read sets, the nearest-reference family match rate was 1.000, while the exact-virus match rate was 0.455. The lower exact-virus agreement mostly came from rhinovirus A/B/C calls being nearest to another rhinovirus reference in 5-mer space. The high-confidence Influenza C sample clustered with Orthomyxoviridae, and the high-confidence HCoV OC43 sample clustered with Coronaviridae.

**Figure 2.** 5-mer/SVD embedding of respiratory-virus references and sample-derived read sets. Sample read sets cluster with the expected viral families: rhinovirus calls fall in the Picornaviridae cluster, the Influenza C sample falls with Orthomyxoviridae, and the HCoV OC43 sample falls with Coronaviridae. This supports the Part 2 calls at family level, while exact rhinovirus species separation remains weaker.  
**Insert:** `figures/part3_sample_reference_kmer_clusters.png`

## Limitations

This is a deliberately simple analysis. The taxonomic classifier uses a small curated reference panel and best-read alignment rather than a full metagenomic classifier, so it can miss viruses outside the panel and may under-resolve divergent rhinoviruses. The coverage estimates are reference-based and can be affected by ONT error, read length, and uneven sampling. The clustering analysis uses 5-mer composition rather than ANI, phylogenetics, or a pretrained sequence model, so it is best viewed as a sanity check for broad family-level agreement rather than exact species assignment.

## References

[1] Microsoft Research. *Toward AI-resilient screening of nucleic acid synthesis orders: process, results, and recommendations.* https://www.microsoft.com/en-us/research/publication/toward-ai-resilient-screening-of-nucleic-acid-synthesis-orders-process-results-and-recommendations/  

[2] IBBIS. *Common Mechanism performance: resilience to AI.* https://ibbis.bio/our-work/common-mechanism/common-mechanism-performance/  

[3] Synsor paper. *Detection of engineered DNA sequences using alignment-free methods.* https://pmc.ncbi.nlm.nih.gov/articles/PMC11272466/
