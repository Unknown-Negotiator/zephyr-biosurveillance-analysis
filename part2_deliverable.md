# Part 2 - Taxonomic Classification

I analyzed 12 Zephyr respiratory-read swab pools from June 5-13, 2026, covering Boston Park/Tremont St, Copley Square, Central Square, MBTA Davis Square, and MBTA Harvard Station. I mapped the ONT respiratory viral reads against a curated panel of common respiratory-virus references using minimap2, kept the best viral alignment per read, and summarized read support, median alignment identity, genome coverage breadth, and mean depth per virus.

Confidence labels were assigned from read support, coverage, and median identity: high confidence required at least 50 reads, at least 20% reference coverage, and median identity at least 85%; medium confidence required at least 10 reads, at least 2% coverage, and median identity at least 75%; low confidence required at least 3 reads and median identity at least 70%; smaller calls were treated as very low confidence. These thresholds are conservative for rhinoviruses: many rhinovirus calls have thousands of reads and near-complete genome coverage, but only ~72-77% identity to the prototype reference, so I interpret them as strong evidence for rhinovirus-like viruses while being cautious about exact species assignment.

## Main Findings

The dominant signal across the pools was rhinovirus. Rhinovirus A and B were detected in all 12 pools, often with broad genome coverage. Rhinovirus C was detected in 6 pools, but generally with less read support and lower coverage. I also detected two non-rhinovirus respiratory-virus signals with high confidence: Influenza C in 260608-Copl-NAS-P1 and Human coronavirus OC43 in 260609-Cent-NAS-P2. Human coronavirus NL63 appeared as a single-read, very-low-confidence signal and should not be treated as a robust positive.

No SARS-CoV-2, influenza A/B, RSV A/B, metapneumovirus, parainfluenza viruses, adenovirus C, or bocavirus 1 were called in these pools under this reference-panel mapping approach.

## Virus Summary Across Pools

| Virus | Positive pools | Total assigned reads | Max reads in one pool | Median identity | Max coverage breadth | Median coverage breadth |
|---|---:|---:|---:|---:|---:|---:|
| Rhinovirus A | 12 | 170,550 | 89,662 | 72.5% | 100.0% | 97.3% |
| Rhinovirus B | 12 | 62,681 | 30,881 | 74.8% | 100.0% | 88.8% |
| Human coronavirus OC43 | 6 | 20,394 | 20,386 | 97.2% | 88.6% | 16.7% |
| Rhinovirus C | 6 | 677 | 329 | 73.2% | 63.5% | 25.7% |
| Influenza C | 1 | 413 | 413 | 96.0% | 85.8% | 85.8% |
| Human coronavirus NL63 | 1 | 1 | 1 | 98.1% | 4.6% | 4.6% |

## Per-Pool Calls

| Sample | Location/date | Main calls |
|---|---|---|
| 260605-BC-NAS-P2 | Boston Park/Tremont, Jun 5 | Rhinovirus A: 5,116 reads, 100.0% coverage, low exact-reference confidence; Rhinovirus B: 5 reads, 84.7% coverage, low; OC43 and Rhinovirus C: single-read very-low calls |
| 260605-BC-NAS-P3 | Boston Park/Tremont, Jun 5 | Rhinovirus A: 6,561 reads, 100.0% coverage, low exact-reference confidence; Rhinovirus B: 2 reads, very low |
| 260608-Copl-NAS-P1 | Copley Square, Jun 8 | Rhinovirus A: 89,662 reads, 100.0% coverage, low exact-reference confidence; Rhinovirus B: 4,999 reads, 100.0% coverage, medium; Influenza C: 413 reads, 85.8% aggregate coverage, high; Rhinovirus C: 329 reads, 63.5% coverage, low |
| 260608-Copl-NAS-P2 | Copley Square, Jun 8 | Rhinovirus A: 83 reads, 88.3% coverage, low; Rhinovirus B: 7 reads, 41.5% coverage, low |
| 260609-Cent-NAS-P1 | Central Square, Jun 9 | Rhinovirus B: 5,465 reads, 96.8% coverage, medium; Rhinovirus A: 2,814 reads, 87.1% coverage, low exact-reference confidence; HCoV NL63: single-read very-low call |
| 260609-Cent-NAS-P2 | Central Square, Jun 9 | Human coronavirus OC43: 20,386 reads, 88.6% coverage, high; Rhinovirus B: 1,554 reads, 98.6% coverage, low exact-reference confidence; Rhinovirus A: 256 reads, 98.5% coverage, low |
| 260610-MBTA_Da-NAS-P1 | MBTA Davis, Jun 10 | Rhinovirus A: 30,366 reads, 100.0% coverage, low exact-reference confidence; Rhinovirus B: 764 reads, 92.9% coverage, low; Rhinovirus C: 55 reads, 63.5% coverage, low; OC43: single-read very-low call |
| 260610-MBTA_Da-NAS-P2 | MBTA Davis, Jun 10 | Only very-low calls: OC43, Rhinovirus A, and Rhinovirus B, with 2 reads each |
| 260611-BC-NAS-P1 | Boston Park/Tremont, Jun 11 | Rhinovirus A: 1,050 reads, 97.6% coverage, medium; Rhinovirus C: 257 reads, 42.3% coverage, low; Rhinovirus B: 21 reads, 33.4% coverage, medium; OC43: 3 reads, low |
| 260611-BC-NAS-P2 | Boston Park/Tremont, Jun 11 | Rhinovirus B: 30,881 reads, 98.7% coverage, low exact-reference confidence; Rhinovirus A: 20,907 reads, 95.2% coverage, low; Rhinovirus C: 19 reads, 9.1% coverage, low |
| 260611-BC-NAS-P3 | Boston Park/Tremont, Jun 11 | Rhinovirus A: 12,611 reads, 97.0% coverage, low exact-reference confidence; Rhinovirus B: 17 reads, 55.6% coverage, low; Rhinovirus C: 16 reads, 5.7% coverage, medium |
| 260613-MBTA_Ha-NAS-P1 | MBTA Harvard, Jun 13 | Rhinovirus B: 18,964 reads, 96.4% coverage, medium; Rhinovirus A: 1,122 reads, 72.5% coverage, low; OC43: single-read very-low call |

## Genome Coverage

Rhinovirus A showed broad coverage in nearly every pool, reaching complete reference coverage in four samples and near-complete coverage in many others. Rhinovirus B also showed broad coverage, including complete coverage in 260608-Copl-NAS-P1 and 98-99% coverage in 260609-Cent-NAS-P1, 260609-Cent-NAS-P2, 260611-BC-NAS-P2, and 260613-MBTA_Ha-NAS-P1. Rhinovirus C was weaker: the best-supported pools reached 63.5% breadth, but most calls were much lower.

The clearest non-rhinovirus calls were Influenza C and HCoV OC43. Influenza C in 260608-Copl-NAS-P1 had high median identity and broad segment coverage: PB2 86.8%, PB1 87.9%, P3 90.7%, HEF 97.7%, NP 91.3%, and MP 41.3%. HCoV OC43 in 260609-Cent-NAS-P2 covered 27,247 of 30,741 reference bases, or 88.6%, with high depth. Other OC43 calls were supported by only 1-3 reads and are best treated as very-low or low-confidence background signals rather than robust positives.

The full result tables are in `results/virus_calls_by_sample.csv`, `results/virus_summary_across_samples.csv`, and `results/coverage_by_reference_segment.csv`; coverage plots are listed in `results/coverage_figures.csv` and stored under `figures/coverage/`.
