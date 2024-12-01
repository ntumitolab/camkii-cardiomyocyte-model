# CaMKII cardiomyocyte model

- https://github.com/drgrandilab/Morotti-et-al-2014-Mouse-Ventricular-Model
- https://github.com/ntumitolab/Hopkins-CaMKII (private repo)

---

- CaMKII model: "Mechanisms of Ca2+/calmodulin-dependent kinase II activation in single dendritic spines" (Chang et al., 2019, Nature Communications); https://doi.org/10.1038/s41467-019-10694-z
- ROS activation model: Oxidized Calmodulin Kinase II Regulates Conduction Following Myocardial Infarction: A Computational Analysis (Christensen et al. 2009); https://doi.org/10.1371/journal.pcbi.1000583
- Isoproterenol and CaMKII effects: A novel computational model of mouse myocyte electrophysiology to assess the synergy between Na+ loading and CaMKII. (Moritti et al. 2014); https://doi.org/10.1113/jphysiol.2013.266676
- Neonatal rat ventricular myocyte model: Model of Excitation-Contraction Coupling of Rat Neonatal Ventricular Myocytes (Korhonen et al. 2009); https://pmc.ncbi.nlm.nih.gov/articles/PMC2716686/

## Adjusted from the original models

- Caffeine activation of RyR: increasing sensitivity to Sub-SR Ca instead of constant opening.
- Corrected compartment volumes for SR Ca and Sub-SR Ca ODEs.
- Added a simplified beta-adrenergic system with a steady-state algebraic representation.
- Fast sodium channel recovery rate increased by 3 times to accomodate 2Hz and 3Hz pacing.

### Parameters adjusted

- kCaM4_on: 30 -> 15Hz / μM
- kCaM2C_on: 0.92 -> 0.5Hz / μM
- k_phosCaM: 30 -> 10Hz
- CmdnTotal: 50 -> 30μM
