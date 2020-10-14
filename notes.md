- Az is heterozygous
- mouse showing masses of inflamation, which is supposed to be the later stages of the disease

bits
- 5 different transgenic mouse brain models, found to be discordant
- Only showing concordance for up-regulated genes


transcriptome Knock-in model
- simulating early-onset Az
- 4 diff genes that cause Az in humans: mutations in all 3 genes
- All 3 affect oxphos
- energy production is affected, but this probably isn't specific to Az
- interestingly, ALS risk factors are opposite to Az
- "astrocytes"
- neurons are extemely energy-hungry
- iron important in mitochondria, MoA lead-in?

Hypothesis
- iron important for controlling choice of glycolysis and oxphos - iron deficiency hypothesis --> hypoxia --> shift to glycolysis
- at gene regulation level, iron and hypoxia are very similar
- brains accumulates (too much) iron over age
- endolysosome system problem (common in Az): accumulation of ferric iron (less reactive, so safer to store). underlying ferrous iron deficiency, along with ferric iron excess.
- iron deficiency typically drives oxidative stress in mito_c.
- haemochromatosis allele in humans (gene HFE). Iron accumulation nullifies Az risk from FOE4 gene allele

Tm McC
- Extracellular metabolite measurements
- Lactate production
- Max consistency between fluxes to correspond with metab data (eflux2, standalone)
- 1. do some model-refinement - tissue-specific fine-tuning. See Nathan E Lewis (search "Corda", "tropo" - should be in matlab)
- 0. transcriptomics model reduction?
- 2. Thermodynamic metabolic flux analysis (PyTFA). Set metab concentration +/- 1 s.d. This resolves reaction directions for all in the network. This gives info about "pathways turning on and off".

More notes from mlardelli
* oxphos quite a central pathway
* neuron cells have huge energy requirements (see glucose)
* astrocytes
* when neurons get too active, they die (excito-toxicity). Astrocytes are there to "soak this up"
* astrocytes and neurons have different metabolisms
* transcriptome/metabolome data is based on the whole brain; it's a weighted average of astrocytes and neurons
* NB.: cultured cells are normally in a stressed state, and wouldn't behave as they would naturally. mlardelli is skeptical of a lot of *in vitro* work.

See table 1:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5303712/
See fig 2:
https://tbiomed.biomedcentral.com/articles/10.1186/1742-4682-4-48

Alz always enriched for oxphos

APOE2 - protective against alzheimer's
APOE3 - most people are this
APOE4 - increased risk of alz

Nhi: not much diff between
