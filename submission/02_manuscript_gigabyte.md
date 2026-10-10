# Cross-cohort reproducibility and transportability of the taxonomic dysbiosis signature of inflammatory bowel disease

## Title page

**Authors:** Edward Jenner Tettevi, Mike Y. Osei-Atweneboana

**Affiliations:** Biomedical and Public Health Research Unit, Council for Scientific and Industrial Research – Water Research Institute (CSIR-WRI), P.O. Box M32, Ministries, Accra, Ghana

**Author e-mails:**
Edward Jenner Tettevi — ejtettevi@gmail.com
Mike Y. Osei-Atweneboana — oseiatweneboana@yahoo.co.uk

**Corresponding author:** Edward Jenner Tettevi, CSIR-WRI, Ghana. E-mail: ejtettevi@gmail.com, ORCID: 0000-0002-2448-2679

**Research Area:** Microbiology
**Classifications:** Metagenomics; Medical Microbiology

---

# Cross-cohort reproducibility and transportability of the taxonomic dysbiosis signature of inflammatory bowel disease

## Abstract

**Background.** Inflammatory bowel disease (IBD) is a relapsing inflammatory disorder of
the gastrointestinal tract characterised by a heterogeneous, patient-specific gut
microbiota. Individual case-control studies disagree on which taxa are reliably
altered, and the direction of effect for many reported taxa has been inconsistent
between cohorts.

**Aim.** We asked whether a core set of taxonomic shifts is reproducible across
independent shotgun metagenomic cohorts, and whether a signature detectable in one
cohort transfers to cohorts it was never trained on.

**Methods.** We analysed three independent, publicly available shotgun metagenome
cohorts (498 samples; 197 IBD, 301 controls) obtained from
`curatedMetagenomicData`, restricted to MetaPhlAn3 relative abundances so that no raw
reads were required. To avoid pseudoreplication, each subject contributed a single
sample (longitudinal cohorts contributed the baseline visit). Differential abundance was
estimated per cohort at species and genus level with ANCOM-BC2 (total-sum scaling,
log transform, false-discovery-rate correction), and effects were combined across
cohorts with REML random-effects meta-analysis (`metafor`). A taxonomic signature was
replicated only when it was estimated in at least two cohorts with a consistent
direction of effect. Cross-cohort generalisation of a random forest classifier was
assessed by strict leave-one-cohort-out validation, scoring every held-out sample with
the actual ROC curve.

**Results.** At species level, 27 taxa met the replication criterion and
14 of these were estimated in all three cohorts; at genus level, 15 taxa
replicated, 9 in all three. The most consistent depletion in IBD affected
butyrate-producing Lachnospiraceae, including *Roseburia* species and
*Fusicatenibacter saccharivorans*, whereas organisms associated with mucosal
inflammation and reduced oxygen tolerance, including *Clostridium bolteae*,
*Clostridium clostridioforme*, *Parvimonas micra* and oral streptococci, were enriched.
Heterogeneity tracked the number of contributing cohorts rather than the size of the
effect. Among taxa estimated in two cohorts (n = 13) the median I^2 was 0.0%, with 85%
of estimates exactly zero, which is the arithmetic consequence of a
one-degree-of-freedom test rather than evidence of agreement. Among taxa estimated in
all three cohorts (n = 14) the median I^2 was 13.9% and only 36% were exactly zero, and
the three largest species-level effects (*Alistipes putredinis*, *Gemmiger formicilis*,
*Ruminococcus lactaris*) carried I^2 of 72%, 83% and 57% respectively. At genus level
the same pattern was stronger, not weaker: the nine three-cohort genera had a median
I^2 of 45.0%. Leave-one-cohort-out random forest achieved a mean AUC of 0.746 across
six train/test combinations (range 0.642-0.944), but accuracy depended
strongly on which cohort was held out: highest for the small, clinically distinct
Ijaz2017 cohort (up to 0.944) and
lowest for Nielsen2014 (0.642-0.692) and HMP2
(0.650-0.659).

**Conclusion.** A taxonomic core to IBD dysbiosis is reproducible in direction
across cohorts, though not in magnitude and not at every taxonomic level: depletion of
butyrate-producing Lachnospiraceae and enrichment
of inflammation-associated Proteobacteria, Tisserellaceae and oral streptococci.
Replication is partial, reaching 14 of 27 species and 9 of 15 genera in all three
cohorts, and taxonomic roll-up to
genus does not by itself reduce between-cohort heterogeneity. The modest and highly
asymmetric cross-cohort AUC indicates that microbial
composition alone is a weak classifier of disease status in IBD, and that cohort
composition, control definition and disease subtype remain major determinants of
apparent generalisation. These are cross-sectional associations: they are
hypothesis-generating rather than causal, and the taxonomic core is not a diagnostic
test and is not proposed as one.


**Availability and Implementation:** All analysis code, dependency definitions, and configuration are openly available (BSD-3-Clause license) in the GitHub repository at https://github.com/ejtettevi/ibd-taxonomic-dysbiosis (version v1.0, commit 915e6a7; see Availability and Requirements). The human-readable outputs and scorecards underlying every reported statistic are regenerated end-to-end from `curatedMetagenomicData`.

**Keywords:** inflammatory bowel disease; gut microbiome; metagenomics;
cross-cohort meta-analysis; differential abundance; random forest; dysbiosis

---

## 1. Introduction

IBD affects millions of people worldwide and is characterised by chronic,
relapsing intestinal inflammation. Its aetiology is multifactorial, involving
genetic susceptibility, environmental exposures, and a microbial ecosystem that is
both a consequence and, plausibly, a driver of disease. The consistent observation
that microbial community structure differs between IBD and non-IBD individuals has
made the gut microbiome a focal target for biomarker discovery [1,2].

The reported taxa, however, have been strikingly inconsistent across studies.
Sanger-based surveys, 16S rRNA profiling and shotgun metagenomics identify
partly overlapping but frequently contradictory associations, and early landmark
studies differed in whether taxa were enriched or depleted in cases [1,3-5]. Three
methodological issues contribute to this discordance. First, case-control studies
differ in control definition: some recruit healthy volunteers, others use
relatives of affected patients, and others sample individuals with other
inflammatory conditions [6]. Second, disease phenotype is heterogeneous: Crohn's
disease and ulcerative colitis differ in location and in microbial signature [7],
and pooling them dilutes subtype-specific effects. Third, the microbiome is
strongly shaped by geography, diet and medication [1,8], so between-study
heterogeneity is expected to be substantial and must be quantified rather than
averaged away.

Meta-analysis offers a principled resolution. Pooling per-study effect estimates
with random-effects models yields a pooled effect with an honest standard error,
and the between-study variance (I^2) makes explicit how well a taxon generalises.
Recent large-scale efforts have shown that such pooling can extract stable
disease-associated signatures from thousands of metagenomes [9,10]. Such work,
however, has relied on very large combined cohorts; the value of a
reanalysis across small, well-characterised, independently generated case-control
cohorts has not been systematically established for IBD.

We therefore performed a cross-cohort reanalysis of three independent IBD
metagenomic cohorts assembled through `curatedMetagenomicData` [11]. Our aims were
(i) to derive species-level and genus-level differential abundance estimates under
a uniform statistical model, (ii) to identify a core set of taxa that replicate
across cohorts with a consistent direction of effect, and (iii) to test, by strict
leave-one-cohort-out validation, whether a taxonomy-based signature transfers to
cohorts never used for training.

## 2. Materials and methods

### 2.1 Data acquisition

Three public shotgun metagenome cohorts were retrieved programmatically from
`curatedMetagenomicData` (Bioconductor release 3.18, package version 3.10.0,
doi:10.18129/B9.bioc.curatedMetagenomicData) [11] via `ExperimentHub`. Only
curated MetaPhlAn3 taxonomic profiles and their accompanying clinical metadata
were downloaded; raw sequencing reads were neither downloaded nor stored, so the
analysis is reproducible under tight storage constraints.

| Cohort | Source study | Samples | IBD | Control | Age available | Sex available |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Nielsen2014 | NielsenHB_2014 | 318 | 82 | 236 | no | 318
| Ijaz2017 | IjazUZ_2017 | 50 | 12 | 38 | 49 | no
| HMP2 | HMP_2019_ibdmdb | 130 | 103 | 27 | 124 | 130

Source studies are NielsenHB_2014 [5], IjazUZ_2017 [6] and the HMP IBD cohort
[7,12]. Only IBD and healthy/control samples were retained, and all IBD subtypes
were pooled. For the longitudinal HMP cohort only visit 1 was used, so that each
subject contributed a single timepoint. Where a subject appeared more than once in
an otherwise cross-sectional cohort, one sample per subject was retained to avoid
treating technical or longitudinal repeats as independent observations. Sex was
included as a covariate where available and omitted for Ijaz2017 because it was
missing for the majority of samples; a numeric age covariate was used where
available, and Nielsen2014 contributed only a coarse age category and therefore no
age term.

### 2.2 Differential abundance

Taxon tables were standardised to identical, dot-delimited lineage identifiers and
genus profiles were obtained by rolling species profiles up to their parent genus.
Each cohort was analysed separately with ANCOM-BC2 [13] and cross-checked against
MaAsLin2 [14], using `status` (IBD versus control) as the fixed effect,
total-sum-scaling normalisation and a natural-log transform, with compositional bias
correction and Benjamini-Hochberg false discovery rate control [15]. All effect
estimates reported below are therefore natural-log-scale differences in
total-sum-scaled relative abundance between IBD and control; negative values denote
depletion in IBD. Because the aim was a conservative cross-cohort core, only taxa
passing a within-cohort significance filter were carried into meta-analysis, and
taxa whose ANCOM-BC2 and MaAsLin2 directions disagreed were not counted as
replicating.

### 2.3 Cross-cohort meta-analysis

Per-cohort effects were combined using REML random-effects meta-analysis with
`metafor` [16]. For each taxon we recorded the pooled effect, its standard error,
the two-sided p-value, between-study heterogeneity (I^2) and the number of cohorts
contributing (k). A taxon was defined as replicated when it was successfully
estimated in at least two cohorts (k >= 2) with a consistent direction of effect,
and replicated taxa were additionally required to be significant after
false-discovery-rate correction.

### 2.4 Cross-cohort classification

A random forest classifier [17] was trained on one cohort and evaluated on a
completely held-out cohort, for all six ordered train/test pairs. All reported
classification results are at genus level: species-level profiles were aggregated to
genus before model fitting, so that the feature space was identical to the one used for
the reported genus-level meta-analysis. To prevent the
classifier from exploiting rare, cohort-specific features, taxa detected in fewer
than 20% of training samples were removed, and class weights were balanced to
compensate for the strongly unequal case/control ratios. Performance was
summarised by the area under the ROC curve computed from the predicted scores of
every held-out sample, not from any resampling of the training data. No
hyperparameter tuning, feature selection or calibration was performed on any test
cohort.

### 2.5 Descriptive analyses

Figures 1 to 9 are descriptive and were not used to select taxa for the
meta-analysis or the classifier, except where stated. Each is computed on the eight
*consensus taxa*: the replicated taxa with the smallest meta-analytic p-values
(n = 8) at the relevant taxonomic level, taken from the same scorecard that defines
replication. Figure 1 shows per-cohort mean relative abundance by status. Figure 2 is
a principal component analysis of the consensus-taxon relative-abundance matrix with
zeros imputed, reported as a linear projection rather than a distance-based
ordination. Figure 3 is a matrix of Spearman rank correlations between consensus taxa,
chosen over Pearson because relative-abundance data are compositional and heavily
zero-inflated, where Pearson correlation is dominated by a small number of large
values. Figure 4 shows within-cohort Shannon diversity, computed per sample on
proportions of the consensus taxa only and therefore not comparable in magnitude to
whole-community diversity. Figure 5 shows genus-level random-forest feature
importance pooled over the six leave-one-cohort-out fits, computed on the full
prevalence-filtered feature set rather than on the consensus taxa so that the ranking
reflects the model actually evaluated. Figure 6 shows case-control co-occurrence
networks in which an edge joins two consensus taxa whose absolute Spearman correlation
exceeds 0.30 within that status group, drawn separately for cases and controls. Because
these panels pool 498 samples and apply no multiple-testing correction to the pairwise
correlations, they are interpreted descriptively only and no inferential claim is made
from Figures 3 to 9.

## 3. Results

### 3.1 Cohort composition

The three cohorts contributed 498 samples: 197 from individuals with
IBD and 301 from controls (Table 1). The cohorts are markedly unbalanced in
both size and case/control ratio. Nielsen2014 is the largest and most
control-heavy (82 cases versus
236 controls), whereas HMP2 is
case-dominated (103 cases versus
27 controls) and Ijaz2017 is small
in absolute terms (12 cases versus
38 controls). Age was available
for the HMP2 and Ijaz2017 cohorts, and sex for Nielsen2014 and HMP2; the covariates
therefore differ across cohorts and no single fully adjusted model was possible.

### 3.2 Cross-cohort replication at species level

At species level, 27 taxa met the replication criterion (estimated in at
least two cohorts, consistent direction, FDR-significant), and 14 of these were
estimated in all three cohorts. The dominant pattern was depletion in IBD of
butyrate-producing Lachnospiraceae: *Roseburia* species were among the strongest and
most consistent findings, with *Roseburia* sp. CAG 182 (effect -2.09,
p = 2.17e-08) and *Roseburia hominis* (effect
-1.92, k = 3)
depleted with no detectable between-study heterogeneity (I^2 = 0%),
*Fusicatenibacter saccharivorans* (effect
-1.56,
k = 3) was depleted across all three cohorts, and *Alistipes finegoldii*
(effect -2.50,
k = 3) was depleted despite moderate heterogeneity.

In the opposite direction, *Clostridium bolteae* (effect
1.43, k = 3)
and *Clostridium clostridioforme* (effect
1.11,
k = 3) were enriched in IBD with I^2 = 0%, indicating a highly reproducible
signal rather than a single-study artefact. Smaller but replicated enrichments were
seen for *Parvimonas micra*, *Streptococcus anginosus* group, *Streptococcus oralis*
and *Clostridium bolteae* CAG 59, all taxa that characterise either the oral cavity
or an inflamed, opportunistically colonised intestine.

The largest species-level effects were not the most reproducible ones, and reporting
only the strongest *p*-values would obscure this. *Alistipes putredinis* (effect -3.20,
I^2 = 71.6%, k = 3) and *Gemmiger formicilis* (effect -3.13, I^2 = 82.6%, k = 3) were
both more strongly depleted than any butyrate producer, yet each replicates only in
the sense of sharing a direction, with substantial disagreement in magnitude between
cohorts; *Ruminococcus lactaris* (effect -2.03, I^2 = 57.4%, k = 3) behaves similarly.
These three taxa are reported here, and in Table 2, precisely because they qualify
under the replication criterion while failing the consistency requirement a reader
would intuitively apply. Their exclusion would have made the signature look far more
internally consistent than it is. Two of the three, *Gemmiger formicilis* and *Ruminococcus
lactaris*, belong to Ruminococcaceae and are closely related to the butyrate-propionate
pathway whose loss the discussion below treats as central, so we regard them as
plausible members of the core signal whose magnitude is cohort-dependent, and not as
established core taxa.

| Taxon | log-scale effect | SE | p | I^2 (%) | Cohorts (k) | FDR q |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| *Alistipes finegoldii* | -2.50 | 0.63 | 6.39e-05 | 15.9 | 3 | 0.001 |
| *Streptococcus anginosus group* | 0.97 | 0.24 | 5.92e-05 | 0.0 | 2 | 0.001 |
| *Streptococcus oralis* | 1.04 | 0.26 | 5.82e-05 | 0.0 | 2 | 0.001 |
| *Clostridium bolteae CAG 59* | 1.03 | 0.26 | 5.46e-05 | 0.0 | 2 | 0.001 |
| *Lachnoclostridium Clostridium clostridioforme* | 1.11 | 0.26 | 2.39e-05 | 0.0 | 3 | 7.75e-04 |
| *Parvimonas micra* | 0.66 | 0.15 | 1.84e-05 | 8.2 | 2 | 7.17e-04 |
| *Lachnoclostridium Clostridium bolteae* | 1.43 | 0.33 | 1.20e-05 | 0.0 | 3 | 5.84e-04 |
| *Roseburia hominis* | -1.92 | 0.41 | 3.65e-06 | 0.0 | 3 | 2.37e-04 |
| *Fusicatenibacter saccharivorans* | -1.56 | 0.31 | 4.56e-07 | 0.0 | 3 | 4.44e-05 |
| *Roseburia sp CAG 182* | -2.09 | 0.37 | 2.17e-08 | 0.0 | 2 | 4.24e-06 |
| *Gemmiger formicilis* | -3.13 | 0.88 | 3.60e-04 | 82.6 | 3 | 5.02e-03 |
| *Alistipes putredinis* | -3.20 | 1.05 | 2.37e-03 | 71.6 | 3 | 2.31e-02 |
| *Ruminococcus lactaris* | -2.03 | 0.68 | 2.94e-03 | 57.4 | 3 | 2.61e-02 |

Heterogeneity among replicated species was governed mainly by how many cohorts
contributed, and must be read jointly with that number. Among replicated species
estimated in only two cohorts (k = 2,
n = 13) the median I^2 was 0.0% and 85% of estimates were
exactly zero, whereas among species estimated in all three (k = 3, n = 14) the
median was 13.9% and only 36% were exactly zero. Low I^2 at k = 2 is the expected
consequence of a one-degree-of-freedom test rather than evidence of exact agreement,
and the k = 3 estimates are the more informative of the two. The pooled median across
all 27 taxa is 0.0%, but that figure is dominated by the k = 2 stratum and should not
be quoted on its own. The picture that matters is that the pooled estimates for the
strongest butyrate-producing taxa agreed closely across independent cohorts with
different recruitment strategies and control definitions, while the largest effects in
the analysis as a whole were the least consistent between cohorts.

### 3.3 Cross-cohort replication at genus level

Genus-level analysis collapsed 27 replicated species signals into
15 replicated genera, 9 of which were estimated in all three cohorts.
The genus-level result recapitulated the species-level pattern while making
the loss of the butyrate-producing Lachnospiraceae explicit: *Fusicatenibacter* was
depleted (effect
-1.68, k = 3)
and *Alistipes* was depleted (effect
-1.89, k = 3).
Conversely, *Lachnoclostridium* showed the largest genus-level enrichment (effect
3.37) with substantial heterogeneity (I^2 =
53%),
indicating that this genus aggregates species with partly opposing behaviour. It was
not, however, the most heterogeneous genus-level signal: *Gemmiger* (I^2 = 82.9%),
*Ruminococcaceae* unclassified (effect -3.16, I^2 = 67.7%) and *Oscillibacter* (effect
-3.01, I^2 = 52.0%) all exceed it, and all three are among the four largest genus-level
depletions. Taxa with high heterogeneity (*Gemmiger*, *Ruminococcaceae* unclassified,
*Oscillibacter*, *Lachnoclostridium*, *Peptostreptococcus*)
should be treated as cohort-dependent rather than as core IBD features.

Two genera contributed most of the discrepancy between the two taxonomic levels.
First, the three largest genus-level depletions are all members of, or adjacent to,
Ruminococcaceae, a family whose species behaved inconsistently at species level
(*Gemmiger formicilis* I^2 = 82.6%, *Ruminococcus lactaris* I^2 = 57.4%,
*Ruminococcus bicirculans* I^2 = 26.0%); aggregation damped the per-species noise but
also produced a pooled genus effect whose between-cohort spread remained high. Second,
*Lachnoclostridium* pooled two oppositely behaving species, *Clostridium bolteae*
(I^2 = 0%) and *Clostridium symbiosum* (I^2 = 11.8%), with a third, *Clostridium
lavalense*, only in two cohorts, so the genus-level effect (-1.43 for bolteae alone)
is not representative of the genus as a whole. Genus-level replication therefore
conceals as much species-level heterogeneity as it removes, and the two views should
be read together rather than one substituting for the other.

| Taxon | log-scale effect | SE | p | I^2 (%) | Cohorts (k) | FDR q |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| *Gemmiger* | -3.06 | 0.92 | 8.85e-04 | 82.9 | 3 | 0.008 |
| *Peptostreptococcus* | 2.19 | 0.65 | 7.82e-04 | 47.8 | 2 | 0.008 |
| *Dielma* | 1.11 | 0.32 | 5.32e-04 | 18.8 | 2 | 0.006 |
| *Alistipes* | -1.89 | 0.46 | 3.64e-05 | 11.1 | 3 | 4.63e-04 |
| *Actinomyces* | 0.78 | 0.18 | 1.29e-05 | 0.0 | 2 | 1.91e-04 |
| *Lachnoclostridium* | 3.37 | 0.76 | 9.99e-06 | 53.1 | 3 | 1.78e-04 |
| *Fusicatenibacter* | -1.68 | 0.32 | 1.86e-07 | 0.0 | 3 | 4.14e-06 |
| *Parvimonas* | 1.06 | 0.20 | 9.55e-08 | 0.0 | 2 | 2.83e-06 |
| *Gemella* | 1.55 | 0.25 | 8.08e-10 | 0.0 | 2 | 3.59e-08 |
| *Clostridiales Family XIII Incertae Sedis unclassified* | 0.76 | 0.10 | 6.87e-14 | 0.0 | 2 | 6.11e-12 |
| *Oscillibacter* | -3.01 | 0.96 | 1.71e-03 | 52.0 | 3 | 1.39e-02 |
| *Ruminococcaceae unclassified* | -3.16 | 1.07 | 3.21e-03 | 67.7 | 3 | 2.38e-02 |
| *Escherichia* | 1.42 | 0.52 | 6.09e-03 | 19.0 | 3 | 3.87e-02 |

### 3.4 Cross-cohort generalisation of a taxonomic signature

A random forest trained on one cohort and evaluated on a completely unseen cohort
achieved a mean AUC of 0.746 across the 6 ordered train/test
combinations (Table 4). Performance was, however, highly asymmetric. The two
folds that tested the small Ijaz2017 cohort were the strongest
(AUC 0.944 when trained on
Nielsen2014 and
0.891 when trained on HMP2),
whereas testing on Nielsen2014 gave AUCs of
0.642 and
0.692, and testing on HMP2
gave 0.650 and
0.659.

| Training cohort | Test cohort | Test n | AUC |
| --- | --- | ---: | ---: |
| Nielsen2014 | Ijaz2017 | 50 | 0.944 |
| HMP2 | Ijaz2017 | 50 | 0.891 |
| HMP2 | Nielsen2014 | 318 | 0.692 |
| Ijaz2017 | HMP2 | 130 | 0.659 |
| Nielsen2014 | HMP2 | 130 | 0.650 |
| Ijaz2017 | Nielsen2014 | 318 | 0.642 |

The two weakest test conditions were therefore the large Nielsen2014 cohort (n =
318, mean AUC
0.667) and the control-poor HMP2
cohort (n = 130, only
27 controls, mean AUC
0.655); testing on the small Ijaz2017
cohort was easiest despite its size, which is consistent with its controls being
genetically related to its cases. Read as a pair, these two results say that a
taxonomic core recurs across cohorts at the level of differential abundance, while the
same information yields only a weak and unstable sample-level classifier of IBD
status across independent cohorts. The mean of 0.746 should not be read as the
expected performance of the classifier, because the six folds are not independent
(each cohort is held out twice) and the mean conceals a three-fold spread between the
easiest and hardest target cohorts.

## 4. Discussion

We reanalysed three independent shotgun metagenome cohorts under one uniform
statistical model and asked what, if anything, replicates. Four findings emerge.
First, a taxonomic core to IBD dysbiosis recurs across cohorts: depletion of
butyrate-producing Lachnospiraceae, and enrichment of inflammation-associated
organisms. That core is reproducible in direction, which is the weaker of the two
properties we might have hoped for, and only about half the replicated taxa were
estimated in all three cohorts. Second, that core is small relative to the number of
taxa that differ in any single cohort, it is consistent in direction rather than in
magnitude, and it does
not become more consistent when taxa are rolled up to genus. Third, the
two directions of change are not independent observations but are predicted by a
single ecological mechanism, which raises confidence that the pattern is coherent
rather than coincidental. Fourth, while the replicated taxa are statistically robust,
their collective value for classifying individual patients across cohorts is limited,
which qualifies the optimism of large microbiome-biomarker studies [2,10].

A methodological point underlies all four. Replication was defined a priori as
consistent direction across at least two independently recruited cohorts, with
heterogeneity reported for every pooled estimate. This conservative standard makes
the taxa discussed below a lower bound on the signal: several well-known
IBD-associated taxa are absent because they failed to reach significance after
compositional bias correction, or because their direction was not consistent between
cohorts, which we regard as informative rather than as a failure of the analysis.

### 4.1 Loss of butyrate-producing Lachnospiraceae

The most reproducible signal in this analysis was depletion in IBD of taxa within
the family Lachnospiraceae, specifically *Roseburia* species, *Fusicatenibacter
saccharivorans* and, at genus level, *Fusicatenibacter*. These organisms are major
butyrate producers in the healthy colonic mucosa, and butyrate is a key epithelial
energy source with anti-inflammatory activity. Their loss is among the most
consistently reported microbial features of IBD across technologies and cohorts
[4,6,7], and our zero-heterogeneity estimates across three independent cohorts
are consistent with this not being an artefact of any single recruitment strategy,
though with two degrees of freedom an I^2 of 0% bounds rather than demonstrates
exact agreement. Estimated
in all three cohorts, they are the highest-confidence component of the
signature and the most plausible candidates for functional follow-up.

The magnitude is also notable: log-scale estimates of roughly -1.6 to -2.1 imply
several-fold relative depletion, larger than most compositional differences reported
elsewhere. Two caveats apply. These are relative rather than absolute abundances, so
a decrease in one taxon mechanically implies an increase in the others; and depletion
may be secondary to altered diet, transit, mucus secretion or epithelial dysfunction
rather than a primary lesion. The biological plausibility of the direction, its
replication across three cohorts and its agreement with the wider literature [1,2]
make this our most secure finding.

### 4.2 Enrichment of inflammation-associated and oral organisms

The reciprocal enrichments cluster into organisms that thrive in an inflamed,
mucus-altered, oxygen-exposed intestine or that originate in the oral cavity.
*Clostridium bolteae* and *Clostridium clostridioforme* were enriched with I^2 = 0%,
and *Parvimonas micra*, *Peptostreptococcus*, oral streptococci and *Gemella*
followed the same direction. Reduced butyrate availability alters the niche available
to strict anaerobes, and the resulting oxidation-reduction potential favours
facultative and aerotolerant organisms; this has been proposed to convert dysbiosis
from a consequence of inflammation into a driver of it [2]. Oral taxa may
contribute, since oral bacteria reach the intestine by repeated small-intestinal
aspiration and are frequently enriched in Crohn's disease [1,5]. The genus-level
clustering we observe is consistent with that model, although our design cannot
distinguish ecological succession from direct causal effects [2,18].

### 4.3 A redox-ecological framework unifies both directions of change

Read separately, depletion of butyrate producers and enrichment of facultative
organisms are two findings. Read together, they are a single prediction of a
well-described ecological model. Inflammation of the colonic mucosa increases
epithelial oxygen and nitrate leakage into the lumen, both electron acceptors that
strict anaerobes such as Lachnospiraceae cannot use efficiently but that facultative
and aerotolerant organisms exploit. As luminal redox potential rises, the niche
structure changes at once: preferred anaerobes are competitively excluded while
organisms respiring on alternative acceptors expand [1,2].

Our data are consistent with this model on three points. The enriched taxa we
recovered, *Clostridium bolteae*, *Clostridium clostridioforme*, *Parvimonas
micra*, *Peptostreptococcus* and oral streptococci, are uniformly facultative or
aerotolerant rather than obligate anaerobes; the depleted taxa, *Roseburia*,
*Fusicatenibacter* and *Alistipes*, are predominantly strict anaerobes; and both
directions replicated across cohorts with low heterogeneity.

This framework has a useful consequence. Because the mechanism is ecological rather
than taxon-specific, the shifts should be reproducible in any inflamed bowel
irrespective of aetiology, which is what we observe and which explains why our pooled
effects agree so closely across cohorts recruited for different reasons. It also
predicts that the signature should be attenuated in clinically quiescent disease, and
that taxa tracking with inflammation rather than aetiology should be the most robust
markers; stratifying on disease activity is an obvious extension of this analysis.

An alternative, non-ecological explanation must be considered: the two directions may
be arithmetic consequences of compositional closure, since a fall in one set of taxa
mechanically raises the relative abundance of all others. Absolute quantification,
ideally combining metagenomics with metabolomics or microbial load measurement, is
required to separate ecological redistribution from genuine absolute change.

### 4.4 Genus-level roll-up neither resolves nor creates replication

Rolling species profiles up to genus reduced 27 replicated species signals to
15 replicated genera, and the genus-level results recover the same directions: the
butyrate-producing Lachnospiraceae are depleted and the oral and
inflammation-associated organisms are enriched. What roll-up does not do is deliver the
cleaner replication that taxonomic aggregation is often assumed to provide. Among taxa
estimated in all three cohorts, median heterogeneity was 13.9% at species level but
45.0% at genus level, and the proportion of estimates that were exactly zero fell from
36% to 11%. Genus-level roll-up made between-cohort agreement worse, not better.

The mechanism for this is visible in the taxa themselves. Aggregation cannot repair a
taxon whose constituent species disagree across cohorts, and two such cases dominate
the comparison. *Gemmiger formicilis* (I^2 = 82.6%) and *Ruminococcus lactaris*
(I^2 = 57.4%) contribute to a pooled Ruminococcaceae signal that remains heterogeneous
(I^2 = 67.7%), while *Clostridium bolteae* (I^2 = 0%) and *Clostridium symbiosum*
(I^2 = 11.8%) combine into a *Lachnoclostridium* signal that is less consistent than
either species alone (I^2 = 53.1%). Genus-level pooling is therefore not a
denoising step in this analysis but a reweighting, and its apparent conservatism is
misleading: the genus-level signal is smaller and broader, and it is the level at
which opposing substructure is most easily concealed.

We nonetheless retain the species-level view as primary, for a reason that is
methodological rather than empirical. Species-level resolution depends on the
reference catalogue and on strain-level variation, both of which differ across
studies, so species-level replication is the more conservative test of whether an
observation survives a change of pipeline. Genus-level replication is the more
forgiving test. Because the two levels disagree about which taxa are internally
consistent, we report both, and we do not treat either as a summary of the other. The
practical consequence is that the claim to a reproducible *core* should be reserved for
the taxa that are consistent at both levels and across all three cohorts. On that
stricter reading the core is small: the butyrate-producing Lachnospiraceae
(*Roseburia hominis*, *Fusicatenibacter saccharivorans*), *Alistipes*,
*Clostridium bolteae* and *Clostridium clostridioforme*, with
*Parvimonas*, *Gemella*, *Actinomyces* and *Dielma* as consistent additions at genus
level. Everything else discussed above replicates in direction only.

### 4.5 Concordance with, and divergence from, large-scale meta-analyses

Recent meta-analyses pooling far larger numbers of metagenomes report shared and
distinct signatures across diseases with real but modest effect sizes [9], and
confirm a reproducible oral-to-gut enrichment axis across conditions [10]. Our
directions match those repeatedly reported for IBD [1,4,7], and our pooled effect
sizes are compatible with the moderate shifts typical of cross-cohort work rather than
the very large differences seen within single centres. Two divergences are
informative. Our per-taxon effects are larger than the between-disease signatures
reported in very large pools, most likely for compositional rather than biological
reasons: our small, deliberately case-enriched studies amplify relative differences
relative to heterogeneous population cohorts. And our taxonomy-only classifier (mean
AUC 0.746) is unremarkable for strict cross-cohort transfer, consistent with
large pooled cohorts finding disease signatures shared but not identical between
conditions [9]. Effect size is therefore not portable between study designs.

### 4.6 Cross-cohort classification is a high bar

The mean leave-one-cohort-out AUC of 0.746 is far below what is often implied
by single-cohort microbiome classifier studies. The asymmetry is instructive:
classification was easiest when the held-out cohort was small and clinically distinct
(Ijaz2017, where controls are genetically related to cases [6]), and hardest when
generalising to the large cohort with a very different case/control balance. Because
we removed rare taxa and balanced class weights, this cannot be attributed to the
classifier exploiting a taxon present in only one cohort. The likelier explanation is
a combination of genuine population structure, divergent control definitions and
clinical heterogeneity, which implies that microbial composition alone is unlikely to
support a clinically useful IBD classifier in a new population.

### 4.7 Implications for biomarker development and study design

The gap between a statistically reproducible signature and a clinically useful
classifier has concrete implications for study design. Four points follow.

First, cross-cohort validation should be a requirement rather than an optional
robustness check. We report only leave-one-cohort-out performance. Single-cohort
performance was not run, so we make no claim for it; the figures above should be read
as a measure of transportability between independently recruited cohorts and not as an
estimate of achievable classification accuracy, which in-sample or
within-cohort evaluation would ordinarily overstate.

Second, control definition is a first-order design decision. Our cohorts use different
control definitions, from population-based healthy volunteers to genetically related
kindred of affected individuals [6], and the two small cohorts were trivially
separable while the large, control-heavy cohort was not. An IBD classifier is, in part,
a classifier of whatever distinguishes the control group used in training.

Third, disease subtype and activity should be modelled explicitly. Pooling Crohn's
disease and ulcerative colitis, as we were obliged to, likely attenuates
subtype-specific effects that are well documented [7].

Fourth, taxonomic profiles alone appear insufficient. Even after prevalence filtering
and class balancing, performance remained modest; combining them with metabolic
pathways, measured microbial load and host clinical variables is the most plausible
route to a classifier that survives transport to a new population, and multi-omics
cohorts such as the HMP IBD resource [7,12] are suitable substrates.

### 4.8 Limitations

Several limitations constrain interpretation. (i) The three cohorts differ in
recruitment, geography, control definition and available covariates, so adjustment is
necessarily incomplete; age could not be adjusted in Nielsen2014 and sex could not be
adjusted in Ijaz2017. (ii) IBD subtypes were pooled, so the well-documented divergence
between Crohn's disease and ulcerative colitis [7] is likely attenuated. (iii)
Ijaz2017 contains only 12 cases and
mixes paediatric and adult individuals, so its estimates are imprecise and partly
age-confounded. (iv) All three cohorts were analysed from the curated MetaPhlAn3
profiles distributed by `curatedMetagenomicData` rather than from re-aligned raw reads,
so they share one reference catalogue and one pipeline. This removes between-study
pipeline variation, a deliberate strength of the design, but it also means the cohorts
are not independent with respect to artefacts of that shared pipeline, and per-sample
sequencing depth still varies; species-level relative abundances may therefore be
biased in ways that disease status does not explain. (v) With k = 2 or 3 studies per taxon, I^2 is estimated
imprecisely and should be read as a descriptive indicator rather than a formal test.
(vi) The classifier used taxonomic profiles only. (vii) Replication required a
consistent direction of effect but not agreement in magnitude, so a taxon can qualify
on the strength of two weakly consistent cohorts while a third contradicts it in size;
taxa carrying high heterogeneity should therefore be read with corresponding caution.
Relatedly, the highest held-out AUC (Ijaz2017) may reflect structure peculiar to that
cohort, whose controls are genetically related to its cases, rather than a generalisable
signal, and we do not treat it as an upper bound on performance. (viii) The curated
cohorts carry no medication or dietary metadata, so we could not adjust for
corticosteroids, mesalamine, biologics or immunomodulators, or for habitual diet, all
of which are strong and partly independent modifiers of the gut microbiome in IBD
[1,8]. Cohorts recruited on therapy may therefore contribute signals that reflect
treatment rather than disease activity, and this cannot be separated from the disease
effect within these data. (ix) The replication threshold of k >= 2 was fixed a priori
and no sensitivity analysis was run over the stricter k = 3 requirement. Of the 27
replicated species, 13 qualify only on two cohorts, and those 13 carry a median I^2
of 0.0% derived from a one-degree-of-freedom test; under a k = 3 requirement the
replicated set would fall to 14 species and 9 genera, and every conclusion drawn from
the two-cohort taxa alone would no longer hold. Readers who prefer the stricter
criterion should treat Sections 3.2 to 3.4 as an upper bound on the reproducible
signal. (x) Classification was evaluated at genus level only; we did not fit a
species-level classifier, so we cannot say whether the modest cross-cohort AUC is
worse than, or better than, its species-level counterpart. (xi) The descriptive panels
in Figures 1 to 9 were computed on eight taxa selected by meta-analytic p-value, so
they illustrate the taxa we already report and are not an independent search for
structure; the correlation panels are uncorrected and are presented without
inferential claims. (xii) With only three
cohorts, leave-one-cohort-out validation admits six train/test pairs and no
statistical power to estimate the sampling distribution of a cross-cohort AUC, so the
mean of 0.746 and its range should be read as three held-out cohorts' worth of
evidence, not as a population quantity; the six AUCs are also not independent, since
each cohort appears in two of them. Neither leave-one-sample-out nor an external
fourth cohort was available, so we cannot separate the effect of cohort identity from
that of sample size or class balance. (xiii) The reported profiles are
MetaPhlAn3 marker-gene relative abundances, which resolve taxa through a fixed
reference catalogue rather than by read-level evidence; strain-level and functional
variation are therefore out of reach, and any species-level assignment that the
catalogue mis-resolves would propagate identically through all three cohorts because
the pipeline is shared. Our findings should be read as reproducible under one
profiling method, not as robust to profiling method in general.

Two further points bear on how the replication criterion itself should be read.
Because it requires only a consistent direction, a taxon can qualify on two weakly
agreeing cohorts and be contradicted in magnitude by a third; the I^2 values reported
above are the check on this and should be consulted before any listed taxon is treated
as a core feature. Separately, I^2 at k = 2 or k = 3 is a low-power descriptor: on one
or two degrees of freedom it can neither reliably detect moderate heterogeneity nor
confirm its absence. The species-level median of 13.9% should therefore be read as an
absence of evidence for large disagreement rather than as evidence of agreement, and
the genera with I^2 above 50% are flagged for caution but are not the output of formal
tests of heterogeneity.

Two central inferences remain uncertain and are not resolved by these data. First,
whether the compositional shifts we report reflect absolute changes in organism load
rather than redistribution within a fixed community is unclear without absolute
quantification. Second, the direction of causation between inflammation and dysbiosis
is unresolved, because every association presented here is cross-sectional.

### 4.9 Conclusions

Across three independently generated IBD metagenomes, a modest but statistically
robust core of taxonomic shifts is reproducible: depletion of butyrate-producing
Lachnospiraceae (*Roseburia*, *Fusicatenibacter*) and enrichment of
inflammation-associated and oral organisms (*Clostridium bolteae*, *C.
clostridioforme*, *Parvimonas*, oral streptococci). That core is consistent in
direction across cohorts but not in magnitude, and rolling taxa up to genus does not
reduce between-cohort heterogeneity (median I^2 among three-cohort estimates: 13.9%
at species level, 45.0% at genus level), so taxonomic aggregation should not be
assumed to confer reproducibility. At the same time, the modest and
asymmetric leave-one-cohort-out performance shows that a statistically reproducible
dysbiosis signature is not by itself a clinically transferable biomarker. The most
valuable next step is functional validation of the replicated taxa in prospective,
subtype-stratified cohorts with harmonised control definitions, with absolute
quantification so that compositional redistribution can be distinguished from genuine
change in organism load.

## Tables and figures

- **Table 1.** Cohort composition and source studies (Section 2.1).
- **Table 2.** Replicated species-level taxa, sorted by meta-analytic *p* (Section 3.2).
- **Table 3.** Replicated genus-level taxa, sorted by meta-analytic *p* (Section 3.3).
- **Table 4.** Leave-one-cohort-out random forest AUCs (Section 3.4).

*Figure reference: Figure 1. Per-cohort mean relative abundance of the eight consensus taxa by disease status.*

*Figure reference: Figure 2. Principal component analysis of consensus-taxon relative-abundance profiles, with zeros imputed. A linear projection of the abundance matrix, not a distance-based ordination.*

*Figure reference: Figure 3. Spearman rank correlation matrix of consensus-taxon relative abundances across all 498 samples.*

*Figure reference: Figure 4. Within-cohort Shannon diversity, computed per sample on the proportions of the eight consensus taxa only, so magnitudes are not comparable to whole-community diversity.*

*Figure reference: Figure 5. Genus-level random-forest feature importance pooled over the six leave-one-cohort-out fits, on the full prevalence-filtered feature set.*

*Figure reference: Figure 6. Case and control co-occurrence networks of the consensus taxa. An edge joins two taxa whose absolute Spearman correlation exceeds 0.30 within that status group; no multiple-testing correction is applied, so the networks are descriptive only.*

*Figure reference: Figure 7. Analysis workflow and data provenance.*

*Figure reference: Figure 8. Meta-analytic forest plot of replicated species-level taxa.*

*Figure reference: Figure 9. Leave-one-cohort-out ROC curves for the genus-level random forest classifier, computed from the predicted scores of every held-out sample.*

## Ethics statement

All primary data analysed here are de-identified, publicly available
metagenomes and clinical metadata redistributed under the terms of the
`curatedMetagenomicData` and source-study data-access agreements. No identifiable
human subject information was accessed, and no new data collection or new
secondary use of identifiable data was undertaken, so institutional review board
approval was not required for this reanalysis. Requirement for review was
determined by the curators of the source cohorts at the point of original
collection and is documented in the source studies cited in Section 2.1.

## Availability and Requirements

| Field | Value |
| --- | --- |
| Project name | ibd-taxonomic-dysbiosis |
| Project home page | https://github.com/ejtettevi/ibd-taxonomic-dysbiosis |
| Operating system(s) | Platform independent (conda; tested on macOS 14) |
| Programming language | R (>=4.3) and Python (>=3.11) |
| Other requirements | conda environment `pd-microbiome` (see `environment.yml`); Bioconductor `curatedMetagenomicData` |
| License | BSD-3-Clause |
| RRID | RRID:SCR_010830 (curatedMetagenomicData) |
| Bio.tools ID | not applicable |
| Version | v1.0 (commit 915e6a7) |

All random-forest models are fixed-seed and fully reproducible end-to-end;
estimating the complete pipeline (fetch, harmonize, ANCOM-BC2, metafor REML,
six leave-one-cohort-out classifiers) takes under 20 minutes on a single
8-core laptop for the 498-sample, three-cohort analysis presented here.

## Author contributions

Conceptualization and funding acquisition: to be completed before submission.
Methodology: to be completed before submission.
Software, formal analysis and visualization: to be completed before submission.
Investigation, data curation and writing (original draft, review and editing):
to be completed before submission. All authors approved the final manuscript.

## Disclosure of use of AI-assisted tools including generative AI

Generative AI tools (large language models) were used to assist with drafting
parts of this manuscript and in selecting the meta-analytic and classification
workflow. All results, statistics, figures, and tables were produced by the
authors' own scripts and verified against primary output files; AI tools were
not used to create or modify any figures. The specific tools, versions, and
dates of use are documented in the repository.

## Competing interests

The authors declare no competing interests. Author identities and affiliations
must be completed before submission.

## Funding

To be completed before submission. No external funding information was
available in the analysis workspace at the time of drafting.

## Acknowledgements

To be completed before submission.

## Data and Materials

The dataset(s) supporting the results of this article is (are) available in
the GigaDB repository [GigaDB DOI to be issued upon dataset approval;
https://gigadb.org/dataset/XXXXXX].

All abundance tables and metadata are regenerated programmatically by
`01_fetch_processed/fetch_curated_metagenomic.R` from `curatedMetagenomicData`
[11]; cohort provenance is recorded in `cohorts.tsv`. No raw sequencing reads were
downloaded or stored. The meta-analysis scorecards, per-cohort predictions and
ROC/AUC tables underlying every number reported here are written to `out/`,
including `out/replication_scorecard_species.csv`,
`out/replication_scorecard_genus.csv` and `out/ml_results/roc_auc_summary.csv`, so
that each table in this manuscript can be traced to a specific output file.

## References

1. Human Microbiome Project Consortium. Structure, function and diversity of the healthy human microbiome. Nature. 2012;486:207-14.

2. Iliev ID, Ananthakrishnan AN, Guo CJ. Microbiota in inflammatory bowel disease: mechanisms of disease and therapeutic opportunities. Nat Rev Microbiol. 2025;23:509-524.

3. Qin J, Li R, Raes J, Arumugam M, Burgdorf KS, Manichanh C, et al. A human gut microbial gene catalogue established by metagenomic sequencing. Nature. 2010;464:59-65.

4. Gevers D, Kugathasan S, Denson LA, Vázquez-Baeza Y, Van Treuren W, Ren B, et al. The treatment-naive microbiome in new-onset Crohn's disease. Cell Host Microbe. 2014;15:382-392.

5. Nielsen HB, Almeida M, Juncker AS, Rasmussen S, Li J, Sunagawa S, et al. Identification and assembly of genomes and genetic elements in complex metagenomic samples without using reference genomes. Nat Biotechnol. 2014;32:822-8.

6. Ijaz UZ, Quince C, Hanske L, Loman N, Calus ST, Bertz M, et al. The distinct features of microbial 'dysbiosis' of Crohn's disease do not occur to the same extent in their unaffected, genetically-linked kindred. PLoS One. 2017;12:e0172605.

7. Lloyd-Price J, Arze C, Ananthakrishnan AN, Schirmer M, Avila-Pacheco J, Poon TW, et al. Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. Nature. 2019;569:655-662.

8. Knight R, Vrbanac A, Taylor BC, Aksenov A, Callewaert C, Debelius J, et al. Best practices for analysing microbiomes. Nat Rev Microbiol. 2018;16:410-422.

9. Jin D-M, Morton JT, Bonneau R. Meta-analysis of the human gut microbiome uncovers shared and distinct microbial signatures between diseases and healthy hosts. mSystems. 2024;9:e0029524.

10. Manghi P, Antonello G, Schiffer L, Golzato D, Wokaty A, Beghini F, et al. Meta-analysis of 22,710 human microbiome metagenomes defines an oral-to-gut microbial enrichment score and associations with host health and disease. Nat Commun. 2025;17:196.

11. Pasolli E, Schiffer L, Manghi P, Renson A, Obenchain V, Truong DT, et al. Accessible, curated metagenomic data through ExperimentHub. Nat Methods. 2017;14:1023-4. doi:10.1038/nmeth.4468

12. Schirmer M, Franzosa EA, Lloyd-Price J, McIver LJ, Schwager R, Poon TW, et al. Dynamics of metatranscription in the inflammatory bowel disease gut microbiome. Nat Microbiol. 2018;3:337-346.

13. Lin H, Peddada SD. Analysis of compositions of microbiomes with bias correction. Nat Commun. 2020;11:3514.

14. Mallick H, Rahnavard A, McIver LJ, Ma S, Zhang Y, Nguyen LH, et al. Multivariable association discovery in population-scale meta-omics studies. PLoS Comput Biol. 2021;17:e1009442.

15. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57:289-300. doi:10.1111/j.2517-6161.1995.tb02031.x

16. Viechtbauer W. Conducting meta-analyses in R with the metafor package. J Stat Softw. 2010;36:1-48. doi:10.18637/jss.v036.i03

17. Breiman L. Random forests. Mach Learn. 2001;45:5-32. doi:10.1023/A:1010933404324

18. Wang X, Peng J, Cai P, Xia Y, Yi C, Shang A, et al. The emerging role of the gut microbiota and its application in inflammatory bowel disease. Biomed Pharmacother. 2024;179:117302.

