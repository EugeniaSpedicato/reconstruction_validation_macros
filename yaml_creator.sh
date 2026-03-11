FILE=$1
RANGE=$2
NHITS=$3
INFO=$4
NEVENTS=$5

cat << EOF > ${FILE}


outputFile: /mnt/raid10/DATA/espedica/fairmu/TB2025/reco/snakemake/range_${RANGE}_${NHITS}hit_${INFO}.root
numberOfEvents: ${NEVENTS} #required for generation, you can also set it for jobs without generation if you want to process only the first N events
saveParametersFile: false #required for event display to draw detector geometry, not used otherwise
saveGeoTracks: false
detectorConfiguration: TR2025_geometry_ideal.yaml

inputFiles: [/eos/user/r/rpilato/shared_with_Eugenia/MCSamples_validation_WiP_1.1.x/idealGeometry/generated_${RANGE}.root]

#############################
##### DATA PREPROCESSOR #####
#############################

runDataPreprocessor: true         # Unify different data input (real data/digitization) in a single format ready to be used by re$
dataPreprocessorConfiguration:
  dataFormat: MC     # MC, TR23, TR25
  saveRawDataEvent: true          # save the unified data format in output tree as 'RawDataEvent'

##########################
##### RECONSTRUCTION #####
##########################

runReconstruction: true
reconstructionConfiguration:

  general:
    outputFormat: full                             # minimal | analysis | full

  tracking:
    maxNumberOfSharedHits: ${NHITS}                            # Remove tracks that share >N hits with a better one
    #disableNonDiagonalCovarianceMatrixTerms: true

  muonFilter:
    minMuonFilterHits: 3                  # Min hits in the filter required to tag a muon

  vertexing:
    restrictVertexPositionToTarget: false
    useFittedZPositionInKinematicFit: true
    useMuonFilterPid: true                               # Use muon filter tag when for PID in vertx
    forceVertexReconstructionInTarget: 0
    #addMSCorrectionFromTargetInKinematicFit: false
    #refitIncomingTrackWithMSCorrection: false
    #refitMuonTrackWithMSCorrection: false
    #refitElectronTrackWithMSCorrection: false

  calorimeter:
    energyToUse: 1  # for function that use calo energy, will use [0: peak Value, 1:cmsFit, 2: crystalBallFit]

    baselineNoiseReco: [[3.51046, 3.64642, 3.23283, 3.67022, 4.04895], # in adc, 1 Gev = 72 adc * gainFactors
                        [3.25993, 3.07851, 3.64882, 3.54215, 3.13898],
                        [3.33528, 3.28898, 3.79472, 3.687,   3.71632],
                        [2.9034,  3.0851,  3.97006, 3.63701, 3.97955],
                        [3.15612, 3.56306, 3.27315, 3.89946, 3.82752]]

    gainFactors:   [[1.11267, 1.09588, 1.09524, 1.03781, 1.02751],
                    [1.08095, 1.11344, 1.15435, 1.04072, 0.936157],
                    [0.960796, 0.883295, 1.03181, 0.901111, 0.894644],
                    [0.864763, 0.893827, 0.903561, 0.887046, 0.892432],
                    [0.90795, 0.87826, 0.843291, 0.859894, 0.920323]]

    hitPositions: [5,6]


########################
##### EVENT FILTER #####
########################

runEventFilter: false
eventFilterConfiguration:

  generatedEvents:

  digitizedEvents:

  reconstructedEvents:


EOF
