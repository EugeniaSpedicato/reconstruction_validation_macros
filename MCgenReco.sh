FILE=$1
NAME=$2
NHITS=$3

cat << EOF > ${FILE}

outputFile: /mnt/raid10/DATA/espedica/fairmu/TB2025/reco/${NAME}_${NHITS}hit.root #used as output filename
saveParametersFile: false #required for event display to draw detector geometry, not used otherwise
saveGeoTracks: false #required for event display to draw tracks, off by default as the current default Geant physics list produces accurate cascades in the calo and thus A LOT of tracks
#numberOfEvents: 1000000 #required for generation, you can also set it for jobs without generation if you want to process only the first N events

#detectorConfiguration: TR2025_geometry_MC.yaml #no need to give extension (added automatically), configs in run directory are prioritized, if not found, they are looked for in the default directory (common/geom$
detectorConfiguration: 3stations_minbias_patrick.yaml #no need to give extension (added automatically), configs in run directory are prioritized, if not found, they are looked for in the default directory (commo$

#inputDirectory: /mnt/raid10/DATA/espedica/fairmu/TB2025/run8/single_muon_interaction_1/
#inputFiles: ["*"]

inputFiles: [
/eos/user/e/espedica/minbias_3stations/${NAME}.root,
]


runDataPreprocessor: true           # Unify different data input (real data/digitization) in a single format ready to be used by reconstruction, required to be true if runReconstruction: true
dataPreprocessorConfiguration:
  dataFormat: MC                  # MC, TR23, TR25
  saveRawDataEvent: false


runReconstruction: true # if true,  then require runDataPreProcessor: true 
reconstructionConfiguration:

  verbose: false

  runAdaptiveVertexFitter: false

  savedVerticesChi2PerNdfThreshold: 100 #only relevant to the vertices saved to the container with all reconstructed vertices, the event is still considered to be reconstructed if the best vertex has a higher chi2/ndf

  weightLinkedTracksByDepositedCharge: false
  allowTrivialIncomingTracks: false #allows for tracks without stereo hits; implemented for testrun

  addMSCorrectionFromTargetInKinematicFit: true
  addBaselineMSCorrectionToAllTracks: true

  maxNumberOfSharedHits: ${NHITS}
  restrictHitSharingToFirstNModulesInSector: 6

  reassignHitsAfterKinematicFit: false

  refitMuonTrackWithMSCorrection: false
  refitIncomingTrackWithMSCorrection: false

  useFittedZPositionInKinematicFit: true
  restrictVertexPositionToTarget: false

  outputFormat: analysis

  useSeedingSensorOnly: false

  xyHitAssignmentWindow: 0.2 #cm
  uvHitAssignmentWindow: 0.3 #cm


  useMuonFilterTagForPid: false # use results from muon filter reconstruction to tag muons and electron in vertex when possible

  muonTagging: # MUON FILTER RECONSTRUCTION
    # match track to muon filter hits:
    muonMatchingNSigmaTrackToHit: 5.0     # distance max in sigma between track and hit to be associated, default 5.0
    minMuonFilterHits: 3                  # minimum number of hit in muon filter a track is required to match in order to identify a muon, default 3
    sharedHitsAllowedMuonFilter: 1        # nb of hit shared allowed in muon filter, default 1
    # match track to muon track on the next station:
    muonMatchingNSigmaTrackToTrack: 5.0   # distance max in sigma between two tracks from different sector to be associated, default 5.0
    cutAmbiguousMuonIdentification: 1.2   # [> 1] if two incoming track match the same outgoing track this cut is the region below the one we do not try Chi2_ndf_track1/Chi2_ndf_track2 < cut (same for flipped numerator and denominator) , default 1.2


#  alignmentFile: starting_alignment_idealGeom_it10
#  alignmentFile: alignment_run8_MF_it3
  alignmentFile: tracker_MF_it6_run${RUN}


runEventFilter: false # not working for data will be fixed soon 
eventFilterConfiguration:

  generatedEvents:


  digitizedEvents:


  reconstructedEvents:

    saveNotReconstructedEvents: false #do not save events with isReconstructed == false
    maxBestVertexChi2PerNdf: 150 #do not save events with the chi2/ndf of best vertex > 150


EOF
