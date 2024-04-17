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

## nth_mesmer_skim.cpp or nth_mesmer_skim_eff97.cpp

These macros do Giovanni's skimming algorithm on MESMER MC data. Eff97 simulates also the efficiencies of single modules (97%). The final version is the higest one. It was used to evaluate efficiency of skimming code, and also
to study the problem of the events with small number of stubs (probably radoative events).


## data_singleclean_afterskim.cpp

Study the real data sample after skimming of the class single clean mu interaxtions. N of stubs, N of reco tracks etc.

## dati_noReco.cpp and notReco.cpp

'dati_noReco.cpp' is used to study real data (with a selection on the signal) and compare the number of events, in certain angular region, with the theoretical MC expectation, given N number of golden muons. 
This comparison is made in macro
` ` `
/home/espedica/fair_install/instFairRoot/share/MUonE/macros/divide
` ` `
which compare the results coming from 'notReco.cpp' (mesmer events). The same selection is applyied also on MC MESMER events.

## n_pileup.cpp
Test the pile-up code, counting the number of generated muons.

## old_effMC.cpp
Evaluate the reconstruction efficiencies as a function of leptons scattering angle using 6 samples in 6 regions of theta e. Because of that, you need to normalize using the weighted cross section. If you want to use
a single data sample, just remember to put ```wnorm=1``` or to take the value from mesmer ouptu in the ```*.root_mesmer_data/s>txt``` file. You need to pay attention if you use samples with different versions of fairmu,
because in older version you had PDG code of incoming muon = 13, while the outgoing was correctly = -13.

## gen_eff.cpp
It is used to compare LO and NNLO cross section, using ` ` `/home/espedica/fair_install/instFairRoot/share/MUonE/macros/MC.cpp` ` `. It doesn't look to reconstruction, because to compare NLO and LO, or other orders, it is
enoough to just use the generated MCTracks.

## ev_draw.cpp
It is an event drawing that for the moment works on reconstructed events and stubs in fairmu.
