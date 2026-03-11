FILE=$1
RUN=$2
TYPE=$3
NAME=$4
NHITS=$5

# Imposta TAR in base a TYPE
if [ "$TYPE" = "single_muon_interaction_0" ]; then
    TAR=0
elif [ "$TYPE" = "single_muon_interaction_1" ]; then
    TAR=1
else
    TAR=-1
fi

cat << EOF > ${FILE}

outputFile: /mnt/raid10/DATA/espedica/fairmu/TB2025/run${RUN}/${TYPE}/${NAME}_${NHITS}hit_WiP11x.root #used as output filename
saveParametersFile: false #required for event display to draw detector geometry, not used otherwise
saveGeoTracks: false #required for event display to draw tracks, off by default as the current default Geant physics list produces accurate cascades in the calo and thus A LOT of tracks
#numberOfEvents: 1000000 #required for generation, you can also set it for jobs without generation if you want to process only the first N events
detectorConfiguration: TR2025_geometry.yaml #no need to give extension (added automatically), configs in run directory are prioritized, if not found, they are looked for in the default directory (common/geometry)

#inputFiles: [/scratch/cdevanne/data/run10/muedaq03-1750368426.root,] #in case you're not running generation
#inputDirectory: /mnt/raid10/DATA/espedica/fairmu/TB2025/run8/single_muon_interaction_1/
#inputFiles: ["*"]

inputDirectory: /eos/experiment/mu-e/staging/daq/2025/decoded/${RUN}/${TYPE}/
inputFiles: [${NAME}.root]


#############################
##### DATA PREPROCESSOR #####
#############################

runDataPreprocessor: true         # Unify different data input (real data/digitization) in a single format ready to be used by reconstruction, required to be true if runRec$
dataPreprocessorConfiguration:
  verbose: false                  # if true: check trigger consistency between OnLine bits and Offline results and print warnings in case of mismatch
  dataFormat: TR25                # MC, TR23, TR25
  triggerVersion: 3               # 2 (for 2025 runs with RunNumber<24), 3 (for 2025 runs with RunNumber>=24
  saveRawDataEvent: false          # save the unified data format in output tree as 'RawDataEvent'

  serenityToCrystal: [[5,  4,  25, 26, 31],   # Map of calorimeter crystal for real data
                      [7,  6,  30, 23, 24],
                      [2,  1,  0,  28, 29],
                      [20, 12, 3,  11, 27],
                      [8,  22, 21, 9,  10]]

##########################
##### RECONSTRUCTION #####
##########################

runReconstruction: true
reconstructionConfiguration:

  general:
    outputFormat: analysis                             # minimal | analysis | full
#    alignmentFile: alignment_2025_run${RUN}
    alignmentFile: Alignment_Files_TR2025_Data/alignment_run26_muedaq03-1753551112_tchain_startIdealGeom_it10_MFit3.yaml

# alignment_2025_run8 alignment_2025_run11 alignment_2025_run12 alignment_2025_run18 alignment_2025_run32 alignment_2025_run48


  tracking:
    maxNumberOfSharedHits: ${NHITS}
    restrictHitSharingToFirstNModulesInSector: 6         # Apply the above only to early modules (-1 means all)
    disableNonDiagonalCovarianceMatrixTerms: false

  muonFilter:

  vertexing:
    useMuonFilterPid: true                                # Use muon filter tag when for PID in vertex
    reassignHitsAfterKinematicFit: false
    addMSCorrectionFromTargetInKinematicFit: true
    refitMuonTrackWithMSCorrection: false            # Re‑fit outgoing muon with PID specific MS
    refitElectronTrackWithMSCorrection: false        # Re‑fit outgoing electron with PID specific MS
    addBaselineMSCorrectionToAllTracks: true
    restrictVertexPositionToTarget: false
    useFittedZPositionInKinematicFit: true
    forceVertexReconstructionInTarget: ${TAR}
    #notUseMuonFilterPidAboveAngleThreshold: 0.005 # rad   # if both track angle are bigger this value, then muon is identified using scattered angles even if useMuonFilterTagForPid: true 
    #removeVertexMuonAngleAboveThreshold: 0.006 #rad       # if muon angles is bigger than this value,  reject the vertex 

  calorimeter:
    energyToUse: 1  # for function that use calo energy, will use [0: peak Value, 1:cmsFit, 2: crystalBallFit]

    baselineNoiseReco: [[3.51046, 3.64642, 3.23283, 3.67022, 4.04895],
                        [3.25993, 3.07851, 3.64882, 3.54215, 3.13898],
                        [3.33528, 3.28898, 3.79472, 3.687,   3.71632],
                        [2.9034,  3.0851,  3.97006, 3.63701, 3.97955],
                        [3.15612, 3.56306, 3.27315, 3.89946, 3.82752]]

    ### for run <= 8

    # gainFactors:  [[0.965538, 0.957609, 0.954613, 0.951408, 0.963937],
    #               [0.935792, 0.943421, 0.93187,  0.914109, 0.941811],
    #               [0.94212,  0.935432, 0.933649, 0.932904, 0.940076],
    #               [0.945161, 0.950441, 0.927815, 0.929331, 0.958035],
    #               [0.960411, 0.947162, 0.959919, 0.941569, 0.961076]]

    # hitPositions: [6,7]

    ## for run > 8
    gainFactors:   [[1.11267, 1.09588, 1.09524, 1.03781, 1.02751],
                    [1.08095, 1.11344, 1.15435, 1.04072, 0.936157],
                    [0.960796, 0.883295, 1.03181, 0.901111, 0.894644],
                    [0.864763, 0.893827, 0.903561, 0.887046, 0.892432],
                    [0.90795, 0.87826, 0.843291, 0.859894, 0.920323]]

    hitPositions: [5,6]




########################
##### EVENT FILTER #####
########################

runEventFilter: false # I haven't merged yet the fix of eventFilter from master


EOF
