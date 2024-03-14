# Scripts description and how to run

How to run these macros:
` ` ` 
root -l /home/espedica/macros_fairmu/name.cpp
` ` ` 


## nth_new_skim.cpp or nth_new_skim_eff97.cpp

These macros do Giovanni's skimming algorithm on minimum bias MC data. Eff97 simulates also the efficiencies of single modules (97%). The final version is the higest one, in the fifth version there is also the limitation on 
the position of  stubs in station0 (3x3cm^2) in order to select a fiducial region.

##  nth_new_skim_eff97_1hitallowed.cpp
 Same as before but allowing just 1 hit shared for classe of single clean mu interactions.

## data_singleclean_afterskim.cpp

Study the real data sample after skimming of the class single clean mu interaxtions. N of stubs, N of reco tracks etc.

## dati_noReco.cpp and notReco.cpp

'dati_noReco.cpp' is used to study real data (with a selection on the signal) and compare the number of events, in certain angular region, with the theoretical MC expectation, given N number of golden muons. This comparison is made in macro
` ` `
/home/espedica/fair_install/instFairRoot/share/MUonE/macros/divide
` ` `
which compare the results coming from 'notReco.cpp' (mesmer events). The same selection is applyied also on MC MESMER events.
