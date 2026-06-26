# Background: why this dataset

Short scientific context for the K-CHOPORE anac017 study. For the full plain-language
version see the lab notes; this is the repo-level summary.

## The system

mRNA carries chemical marks that tune its stability, splicing and translation. The set of
those marks is the epitranscriptome, and the most abundant internal mark on mRNA is
**m6A** (N6-methyladenosine), written by MTA/MTB/FIP37/VIR/HAKAI, read in an **RRACH**
context. Whether and how stress signalling reshapes this layer in plants is largely open
(Barbieri & Kouzarides 2022; Shen et al. 2023).

Mitochondria sense stress. When the electron transport chain stalls, electrons leak as
reactive oxygen species (ROS), and that signal travels from the mitochondrion back to the
nucleus to change gene expression. This is mitochondrial retrograde signalling
(Mittler et al. 2022; Van Aken 2021). The transcription factor **ANAC017** is the master
regulator of the pathway: it sits inactive at the ER membrane and, on the redox signal, is
released to the nucleus to drive the stress response. AOX1a (AT3G22370) is its canonical
target, an alternative oxidase that relieves the stalled chain.

## The question and the design

Does mitochondrial redox stress remodel the m6A epitranscriptome, and does that remodelling
depend on ANAC017?

To answer both at once, the experiment is a **2×2 factorial**: genotype (WT vs *anac017-1*)
crossed with treatment (Control vs **Antimycin A**), three biological replicates, 12 Direct
RNA-seq libraries. Antimycin A blocks complex III, which triggers the retrograde response on
demand (Shapiguzov et al. 2019). The *anac017-1* knockout is the test: changes that vanish
without ANAC017 are the ANAC017-dependent ones.

Oxford Nanopore **Direct RNA-seq** reads native RNA, so modifications are detected directly
in the signal, with no antibody and no chemical conversion (Furlan et al. 2021).

## What the data show

Mitochondrial stress changes thousands of m6A sites in WT; most of that change collapses in
*anac017-1*. The changed sites are genuine m6A (they overlap the MTA/METTL3-dependent set and
carry RRACH), and the responding genes are enriched in hypoxia, oxygen and energy processes,
consistent with complex-III inhibition. The working model: **Antimycin A → ROS → ANAC017 →
m6A remodelling**. How a transcription factor reaches the m6A machinery is the open question.

This repo reproduces and extends that analysis; see the README for the headline numbers and
`MANUAL.md` for the methods. m6A machinery reference for plants: Wong et al. 2023,
*Plant Physiology*.

## References

- Barbieri I, Kouzarides T (2022). Role of RNA modifications in cancer. *Nat Rev Cancer*.
- Shen L et al. (2023). Messenger RNA modifications in plants. *Trends Plant Sci* / review.
- Mittler R et al. (2022). Reactive oxygen species signalling. *Nat Rev Mol Cell Biol*.
- Van Aken O (2021). Mitochondrial redox systems as central hubs in plant stress signalling.
- Shapiguzov A et al. (2019). Arabidopsis RCD1 coordinates chloroplast and mitochondrial
  functions through interaction with ANAC transcription factors. *eLife*.
- Furlan M et al. (2021). Computational methods for RNA modification detection from nanopore
  direct RNA sequencing. *RNA Biology*.
- Wong CE et al. (2023). m6A in the Arabidopsis response. *Plant Physiology* 191(3):2045-2063.
