#include "MUonEReconstruction.h"

#include "MUonEAlignmentContainer.h"

#include <iostream>
#include <algorithm>

#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairBaseParSet.h"
#include <fairlogger/Logger.h>

#include "TClonesArray.h"
#include "TVector3.h"

MUonEReconstruction::MUonEReconstruction()
	: FairTask("MUonEReconstruction", 0),
		m_output(new MUonERecoOutput())
	{}


MUonEReconstruction::~MUonEReconstruction()
{}


Bool_t MUonEReconstruction::setConfiguration(MUonEDetectorConfiguration const& detectorConfig, MUonEReconstructionConfiguration const& reconstructionConfig) {

	resetDefaultConfiguration();

	m_detectorConfiguration = detectorConfig;
	m_configuredCorrectly = true;
	m_configuredCorrectly = detectorConfig.configuredCorrectly() & reconstructionConfig.configuredCorrectly();

	LOG(info) << "";
	LOG(info) << "Creating Reconstruction object.";

	m_reconstructionConfiguration = reconstructionConfig;

    if(m_configuredCorrectly)
        LOG(info) << "Reconstruction object created successfully.";
    else
        LOG(info) << "Digitization configuration contains errors. Please fix them and try again.";
	
	return m_configuredCorrectly;
}

void MUonEReconstruction::logCurrentConfiguration() {

	m_reconstructionConfiguration.logCurrentConfiguration();
}

void MUonEReconstruction::resetDefaultConfiguration() {

	m_configuredCorrectly = false;

	m_detectorConfiguration.resetDefaultConfiguration();
	m_reconstructionConfiguration.resetDefaultConfiguration();

}



InitStatus MUonEReconstruction::ReInit()
{
	return kSUCCESS;
}

InitStatus MUonEReconstruction::Init()
{
	LOG(info) << "Initializing reconstruction.";
	FairRootManager* ioman = FairRootManager::Instance();

	if(!ioman) 
	{ 
		LOG(fatal) << "No FairRootManager"; 
		return kERROR;
	} 


	MUonEAlignmentContainer alignment;

	if(m_reconstructionConfiguration.isMC()) {

		m_mcTracks = static_cast<TClonesArray*>(ioman->GetObject("MCTrack"));
		m_signalTracks = static_cast<TClonesArray*>(ioman->GetObject("SignalTracks"));

		if(m_detectorConfiguration.hasCalorimeter())
			m_calorimeterDeposit = ioman->InitObjectAs<const MUonECalorimeterDigiDeposit*>("CalorimeterDigiDeposit");

		if(m_detectorConfiguration.hasModules())
			m_stubs = static_cast<TClonesArray*>(ioman->GetObject("TrackerStubs"));

		if(m_reconstructionConfiguration.alignmentFileSet()) {
			
			if(!alignment.readConfiguration(m_reconstructionConfiguration.alignmentFile()))
				return kERROR;
		}

	} else {

		m_dataBx = ioman->InitObjectAs<const std::vector<unsigned short>*>("Bx");
		m_dataLink = ioman->InitObjectAs<const std::vector<unsigned short>*>("Link");
		m_dataStrip = ioman->InitObjectAs<const std::vector<float>*>("LocalX");
		m_dataCicID = ioman->InitObjectAs<const std::vector<float>*>("LocalY");
		m_dataBend = ioman->InitObjectAs<const std::vector<float>*>("Bend");
		m_dataSuperID = ioman->InitObjectAs<const std::vector<unsigned int>*>("SuperID");

		//prepare linkID to internal IDs map
		auto const& stations = m_detectorConfiguration.stations();

		m_linkIDstationIDmoduleIDMap.reserve(stations.size() * 10);

		for(Int_t station_index = 0; station_index < stations.size(); ++station_index) {

			auto const& station = stations[station_index];
			auto const& station_modules = station.modules();

			for(Int_t module_index = 0; module_index < station_modules.size(); ++module_index) {

				m_linkIDstationIDmoduleIDMap.emplace_back(station_index, module_index);
			}
		}

		if(m_reconstructionConfiguration.alignmentFileSet()) {
			
			if(!alignment.readConfiguration(m_reconstructionConfiguration.alignmentFile()))
				return kERROR;
		}
	}


	Register();

	//targets and modules have to be sorted by their Z position
	//this is guaranted when MUonEDetectorConfiguration is used
	//as long as defined stations do not overlap

	//this requirement could be removed with a lot of index tracking involved, especially since the targets wouldn't be sorted
	//since we need to assume some convention for loading data anyway, it's convenient to assume everything is sorted

	//we will use absolute positions of every target and module during reconstruction.

	//load targets (we need to know their positions before we start loading modules)

	//if last station is to be ignored (to be used as muon detection system), do not add the target - this will stop vertexing in it

	auto const& stations = m_detectorConfiguration.stations();

	Int_t N_stations = stations.size();
	if(m_reconstructionConfiguration.ingoreLastStationInVertexing())
		N_stations -= 1;


	//N_stations is one less than actual number of stations if last one is to be ignored
	for(Int_t station_index = 0; station_index < N_stations; ++station_index) {

		auto const& station = stations[station_index];
		auto const& target_conf = station.targetConfiguration();

		if(station.hasTarget())
			m_targets.emplace_back(m_targets.size(), station_index, station.targetPosition(), target_conf);
	}

 	//load modules
	//modules are split into sectors between targets (normally that's equivalent to stations, but 
	//since we allow stations without targets, we need to account for this; we're also doing this separately
	//from targets loading for safety in case there are modules before the target defined in a station;
	//this is just an initialization so performance is irrelevant)
	//sectors are as follows (T - target)
	// 0 T0 1 T1 2 T2 3 T3 4 ... Tn number_of_targets
	//so the index of sector before target is the same as index of the target, and the one after is n+1	

	//note that order in m_modules is fixed and doesn't depend on the sector
	//it is used to map hits to modules when event is processed
	//this is also a case for requiring sorting, as otherwise detector configurations with modules or stations
	//defined in different order wouldn't be compatible

	//we'll also fill a simple map to make assigning hits to modules easier
	m_hitModuleMap.resize(stations.size());

	//if last station is to be ignored, it will be added separately as a new sector
	//this will stop the modules from being used with the previous station

	if(m_reconstructionConfiguration.ingoreLastStationInVertexing()) {

		m_modulesPerSector.resize(m_targets.size() + 2);

	} else {

		m_modulesPerSector.resize(m_targets.size() + 1);
	}

	//N_stations is one less than actual number of stations if last one is to be ignored
	for(Int_t station_index = 0; station_index < N_stations; ++station_index) {

		auto const& station = stations[station_index];
		auto const& station_modules = station.modules();

		m_hitModuleMap[station_index].resize(station_modules.size());

		for(Int_t module_index = 0; module_index < station_modules.size(); ++module_index) {

			auto const& module = station_modules[module_index];

			//assume module is behind the last target
			Int_t module_sector = m_targets.size();
			Double_t module_position = station.modulePosition(module);

			//find the first target such that module is before it
			//since we assume targets (stations) are sorted by their Z position ascending, this approach is valid
			for(Int_t target_index = 0; target_index < m_targets.size(); ++target_index) {

				if(module_position < m_targets[target_index].z()) {

					module_sector = target_index;
					break;					
				}
			}

			//add module
			m_hitModuleMap[station_index][module_index] = std::make_pair(module_sector, m_modulesPerSector[module_sector].size());
			m_modulesPerSector[module_sector].emplace_back(m_modulesPerSector[module_sector].size(), module_sector, station_index, module_index, module_position, m_reconstructionConfiguration.isMC(), module);			
		
			if(m_reconstructionConfiguration.alignmentFileSet()) {

				auto const& module_alignment = alignment.moduleAlignment(station_index, module_index);
				m_modulesPerSector[module_sector].back().setAlignmentParameters(module_alignment.xOffset(), module_alignment.yOffset(), module_alignment.zOffset(), module_alignment.angleOffset(), module_alignment.tiltOffset());
			}

		}
	}

	//repeat the above for the last station and add its entirety as a new sector
	if(m_reconstructionConfiguration.ingoreLastStationInVertexing()) {

		Int_t station_index = stations.size() - 1;
		auto const& station = stations[station_index];
		auto const& station_modules = station.modules();

		m_hitModuleMap[station_index].resize(station_modules.size());

		//note that m_targets.size() is the sector behind last target, since ordering starts from 0 (before first target)
		Int_t module_sector = m_targets.size() + 1;		

		for(Int_t module_index = 0; module_index < station_modules.size(); ++module_index) {

			auto const& module = station_modules[module_index];

			Double_t module_position = station.modulePosition(module);

			//add module
			m_hitModuleMap[station_index][module_index] = std::make_pair(module_sector, m_modulesPerSector[module_sector].size());
			m_modulesPerSector[module_sector].emplace_back(m_modulesPerSector[module_sector].size(), module_sector, station_index, module_index, module_position, m_reconstructionConfiguration.isMC(), module);			
		
			if(m_reconstructionConfiguration.alignmentFileSet()) {

				auto const& module_alignment = alignment.moduleAlignment(station_index, module_index);
				m_modulesPerSector[module_sector].back().setAlignmentParameters(module_alignment.xOffset(), module_alignment.yOffset(), module_alignment.zOffset(), module_alignment.angleOffset(), module_alignment.tiltOffset());
			}

		}
	}


	//we'll also count the number of x/y modules per sector, because if there are only 2, we don't have to look for additional hits
	m_numberOfXModulesInSector.resize(m_modulesPerSector.size());
	std::fill(m_numberOfXModulesInSector.begin(), m_numberOfXModulesInSector.end(), 0);

	m_numberOfYModulesInSector.resize(m_modulesPerSector.size());
	std::fill(m_numberOfYModulesInSector.begin(), m_numberOfYModulesInSector.end(), 0);

	//and stereo modules, so we can avoid unnecessary looping
	m_numberOfStereoModulesInSector.resize(m_modulesPerSector.size());
	std::fill(m_numberOfStereoModulesInSector.begin(), m_numberOfStereoModulesInSector.end(), 0);

	for(Int_t sector = 0; sector < m_modulesPerSector.size(); ++sector) {
		for(auto const& mod : m_modulesPerSector[sector]) {

			if('x' == mod.projection())
				++m_numberOfXModulesInSector[sector];
			else if('y' == mod.projection())
				++m_numberOfYModulesInSector[sector];
			else //any stereo
				++m_numberOfStereoModulesInSector[sector];
		}
	}

	//m_targets, m_modules and the map are not modified in any way during reconstruction
    
	return kSUCCESS;
}

void MUonEReconstruction::SetParContainers()
{}


void MUonEReconstruction::Exec(Option_t* option)
{

	Bool_t event_reconstructed = true;
	Bool_t event_skipped = false;

	//fill these with either digitized MC data or real data
	std::vector<std::vector<std::vector<MUonERecoHit>>> hits_per_module_per_sector; //the indices of the outer vector correspond to the m_modulesPerSector vector

	hits_per_module_per_sector.resize(m_modulesPerSector.size());
	for(Int_t sector = 0; sector < m_modulesPerSector.size(); ++sector) {

		hits_per_module_per_sector[sector].resize(m_modulesPerSector[sector].size());

		for(auto& hits_per_module : hits_per_module_per_sector[sector])
			hits_per_module.reserve(10);
	}

	Double_t event_energy{0};
	Int_t number_of_hits{0};

	if(m_reconstructionConfiguration.isMC())
		number_of_hits = fillMCData(hits_per_module_per_sector, event_energy);
	else
		number_of_hits = fillData(hits_per_module_per_sector, event_energy);



	//indices correspond to the sector's index
	std::vector<std::vector<MUonERecoTrack3D>> reconstructed_tracks_per_sector;
	reconstructed_tracks_per_sector.resize(m_modulesPerSector.size());


	//indices correspond to m_targets vector
	std::vector<MUonERecoVertex> reconstructed_vertices;
	reconstructed_vertices.reserve(10);	

	std::vector<MUonERecoAdaptiveFitterVertex> reconstruced_AF_vertices;
	reconstruced_AF_vertices.reserve(10);

	std::vector<MUonERecoGenericVertex> reconstrucedGenericVertices;
	reconstruced_AF_vertices.reserve(10);


	if(m_reconstructionConfiguration.skipEventsWithMoreThanNHits() > 0 && number_of_hits > m_reconstructionConfiguration.skipEventsWithMoreThanNHits()) {

		event_skipped = true;
		event_reconstructed = false;
		goto end_of_reconstruction;
	}


	if(!reconstruct3DTracks(hits_per_module_per_sector, reconstructed_tracks_per_sector)) {

		if(m_reconstructionConfiguration.verbose())
			std::cout << "Reconstruction failed!" << std::endl;

		event_reconstructed = false;
		goto end_of_reconstruction;
	}

	//construct vertices from a single track before and pairs of tracks after target
	if(!reconstructVertices(hits_per_module_per_sector, reconstructed_tracks_per_sector, reconstructed_vertices, reconstruced_AF_vertices, reconstrucedGenericVertices)) {

		if(m_reconstructionConfiguration.verbose())
			std::cout << "Reconstruction failed!" << std::endl;

		event_reconstructed = false;
		goto end_of_reconstruction;
	}

	end_of_reconstruction:


	if(event_reconstructed) {

		*m_output = MUonERecoOutput(m_eventNumber, event_energy, hits_per_module_per_sector, reconstructed_tracks_per_sector, reconstructed_vertices, reconstruced_AF_vertices, reconstrucedGenericVertices, m_reconstructionConfiguration);

	} else {

		*m_output = MUonERecoOutput(m_eventNumber, event_energy, hits_per_module_per_sector, reconstructed_tracks_per_sector, reconstruced_AF_vertices, m_reconstructionConfiguration, event_skipped);
	}

	++m_eventNumber;
}

Int_t MUonEReconstruction::fillMCData(std::vector<std::vector<std::vector<MUonERecoHit>>>& hitsPerModulePerSector, Double_t& eventEnergy) const {

	eventEnergy = m_detectorConfiguration.hasCalorimeter() ? m_calorimeterDeposit->totalSmearedEnergy() : 0;

	if(!m_detectorConfiguration.hasModules())
		return 0;

	auto number_of_stubs = m_stubs->GetEntriesFast();
	for(Int_t i = 0; i < number_of_stubs; ++i) {

		auto temp_stub = *static_cast<MUonETrackerStub*>(m_stubs->At(i));
		auto module_indices = m_hitModuleMap[temp_stub.stationID()][temp_stub.moduleID()];

		hitsPerModulePerSector[module_indices.first][module_indices.second].emplace_back(temp_stub, module_indices.first, module_indices.second, &m_modulesPerSector, i, m_reconstructionConfiguration.useSeedingSensorOnly());
	}

	return number_of_stubs;
}

Int_t MUonEReconstruction::fillData(std::vector<std::vector<std::vector<MUonERecoHit>>>& hitsPerModulePerSector, Double_t& eventEnergy) const {

	//for now
	eventEnergy = 0;

	if(!m_detectorConfiguration.hasModules())
		return 0;

	for(Int_t stub_index = 0; stub_index < m_dataLink->size(); ++stub_index) {

		auto indices = m_linkIDstationIDmoduleIDMap[(*m_dataLink)[stub_index]];
		auto module_indices = m_hitModuleMap[indices.first][indices.second];
		hitsPerModulePerSector[module_indices.first][module_indices.second].emplace_back((*m_dataStrip)[stub_index], (*m_dataCicID)[stub_index], (*m_dataBend)[stub_index], module_indices.first, module_indices.second, &m_modulesPerSector, stub_index, m_reconstructionConfiguration.useSeedingSensorOnly(), (*m_dataBx)[stub_index], (*m_dataSuperID)[stub_index]);
	}

	return m_dataLink->size();
}


Bool_t MUonEReconstruction::reconstruct2DTracks(std::vector<std::vector<std::vector<MUonERecoHit>>> const& hitsPerModulePerSector, Int_t sector, std::vector<MUonERecoTrack2D>& reconstructed2DTracksInXProjection, std::vector<MUonERecoTrack2D>& reconstructed2DTracksInYProjection) const {

	reconstructed2DTracksInXProjection.reserve(50);
	reconstructed2DTracksInYProjection.reserve(50);

	//note we're only considering modules in the same sector, so it's unnecessary to check for this
	auto const hits_per_module = hitsPerModulePerSector[sector];
	auto const modules = m_modulesPerSector[sector]; 

	if(modules.empty()) return false;

	for(Int_t first_module_index = 0; first_module_index < hits_per_module.size() - 1; ++first_module_index) {

		auto const& first_module = modules[first_module_index];
		auto const& first_module_hits = hits_per_module[first_module_index];


		//ignore U and V modules
		if(first_module.isStereo())
			continue;


		for(Int_t second_module_index = first_module_index + 1; second_module_index < hits_per_module.size(); ++second_module_index) {


			auto const& second_module = modules[second_module_index];
			auto const& second_module_hits = hits_per_module[second_module_index];


			//only modules of the same type
			if(!second_module.sameProjectionAs(first_module))
				continue;

			//create 2D track and assign hits
			//both hits are of the same type and in the same sector
			for(auto const& first_hit : first_module_hits) {
				for(auto const& second_hit : second_module_hits) {

					MUonERecoTrack2D temp_line(first_hit, second_hit);

					//otherwise we can't find more than the two hits already assigned
					//note also that no fit is needed, as the slope and x0 are calculated on creation
					if('x' == first_module.projection() && m_numberOfXModulesInSector[sector] > 2
					|| 'y' == first_module.projection() && m_numberOfYModulesInSector[sector] > 2) {

						//assign hits
						for(Int_t module_index = 0; module_index < hits_per_module.size(); ++module_index) {


							auto const& module = modules[module_index];

							//only compatible modules
							if(!module.sameProjectionAs(first_module))
								continue;						

							//skip modules from which the first two hits are
							if(module.sameModuleAs(first_module) || module.sameModuleAs(second_module))
								continue;

							temp_line.addClosestHit(hits_per_module[module_index], m_reconstructionConfiguration.xyHitAssignmentWindow());

						}//! end of hit assigning


						//can't do that since we now have only 2 x and 2 y modules per station
						//reject tracks with no additional hits assigned
						//if(2 == temp_line.numberOfHits())
						//	continue;


						//refit and check if hits from other modules can be assigned to a fitted line
						Bool_t new_hit_assigned{false};
						Int_t number_of_additional_hits{0};

						do{
							new_hit_assigned = false;
							m_fitter.addFitInfo(temp_line);

							//assign hits
							for(Int_t module_index = 0; module_index < hits_per_module.size(); ++module_index) {

								auto const& module = modules[module_index];

								//only compatible modules
								if(!module.sameProjectionAs(first_module))
									continue;						

								//skip modules from which the first two hits are
								if(module.sameModuleAs(first_module) || module.sameModuleAs(second_module))
									continue;

								//skip modules from which the hits are already assigned to this track
								auto const& track_hits = temp_line.hits();
								if(std::find_if(track_hits.begin(), track_hits.end(), [&module_index](MUonERecoHit const& hit){return module_index == hit.moduleIndex();}) != track_hits.end())
									continue;

								if(temp_line.addClosestHit(hits_per_module[module_index], m_reconstructionConfiguration.xyHitAssignmentWindow())){
									new_hit_assigned = true;
									++number_of_additional_hits;
								}
							}//! end of hit assigning
						} while(new_hit_assigned);						
					}//!end of number of x/y modules check


					auto lineAlreadyInCollection = [](MUonERecoTrack2D const& new_line, std::vector<MUonERecoTrack2D> const& lines) {

						for(auto const& line : lines)
							if(new_line == line) //same hits, see MUonERecoTrack2D.h and MUonERecoHit.h
								return true;

						return false;
					};

					//save track, both modules are of the same type (either X or Y) and in the same sector

					if('x' == first_module.projection() && !lineAlreadyInCollection(temp_line, reconstructed2DTracksInXProjection)) 
						reconstructed2DTracksInXProjection.emplace_back(std::move(temp_line));
						
					else if('y' == first_module.projection() && !lineAlreadyInCollection(temp_line, reconstructed2DTracksInYProjection)) 
						reconstructed2DTracksInYProjection.emplace_back(std::move(temp_line));
						
				} //!loop on hits from second module
			} //! loop on hits from first module
		} //!second module loop
	} //!first module loop



	if(m_reconstructionConfiguration.verbose()) {

		std::cout << "Found " << reconstructed2DTracksInXProjection.size() << " 2D x track(s) in sector " << sector <<"." << std::endl;
		std::cout << "Found " << reconstructed2DTracksInYProjection.size() << " 2D y track(s) in sector " << sector <<"." << std::endl;
	}

	if(reconstructed2DTracksInXProjection.empty() || reconstructed2DTracksInYProjection.empty())
		return false;

	return true;
}

Bool_t MUonEReconstruction::reconstruct3DTracks(std::vector<std::vector<std::vector<MUonERecoHit>>> const& hitsPerModulePerSector, std::vector<std::vector<MUonERecoTrack3D>>& reconstructedTracksPerSector) const {


	for(Int_t sector = 0; sector < reconstructedTracksPerSector.size(); ++sector) {

		std::vector<MUonERecoTrack2D> reconstructed_2DTracks_in_x_projection;
		std::vector<MUonERecoTrack2D> reconstructed_2DTracks_in_y_projection;

		if(!reconstruct2DTracks(hitsPerModulePerSector, sector, reconstructed_2DTracks_in_x_projection, reconstructed_2DTracks_in_y_projection))
			continue;


		//note we're only considering modules in the same sector, so it's unnecessary to check for this
		auto const hits_per_module = hitsPerModulePerSector[sector];
		auto const modules = m_modulesPerSector[sector]; 

		reconstructedTracksPerSector[sector].reserve(reconstructed_2DTracks_in_x_projection.size() * reconstructed_2DTracks_in_y_projection.size());

		
		//create 3D track and assign stereo hits
		for(auto const& x_track : reconstructed_2DTracks_in_x_projection) {
			for(auto const& y_track : reconstructed_2DTracks_in_y_projection) {

				MUonERecoTrack3D track(x_track, y_track, sector, reconstructedTracksPerSector[sector].size());

				if(m_numberOfStereoModulesInSector[sector] > 0) {

					//assign stereo hits
					for(Int_t module_index = 0; module_index < hits_per_module.size(); ++module_index) {

						auto const& module = modules[module_index];

						//only compatible stereo modules
						if(module.isStereo())
							track.addClosestStereoHit(hits_per_module[module_index], m_reconstructionConfiguration.uvHitAssignmentWindow());

					}//! end of hit assigning					
				}

				if(!(m_reconstructionConfiguration.allowTrivialIncomingTracks() && 0 == sector)) {

					//reject trivial tracks
					if(0 == track.numberOfStereoHits()) {

						if(track.xTrack().numberOfHits() < 3 || track.yTrack().numberOfHits() < 3)
							continue;
					}					
				}


				if(m_reconstructionConfiguration.addBaselineMSCorrectionToAllTracks()) {

					m_fitter.addMS(track, m_reconstructionConfiguration.MSCorrectionAssumedBeamEnergy(), m_modulesPerSector);
				}


				//during fitting some outlier hits may be removed and the track may no longer have enough hits
				if(!m_fitter.addFitInfo(track, 0 == sector ? m_reconstructionConfiguration.allowTrivialIncomingTracks() : false, m_reconstructionConfiguration.maxOutlierChi2()))
					reconstructedTracksPerSector[sector].emplace_back(std::move(track));


			} //!ytrack loop
		} //!xtrack loop	


		if(m_reconstructionConfiguration.verbose()) 
			std::cout << "Found " << reconstructedTracksPerSector[sector].size() << " 3D track(s) in sector " << sector << ". ";


		//filter tracks

		//sort the tracks from best to worst, see operator definitions in MUonERecoTrack3D.h
		std::sort(reconstructedTracksPerSector[sector].begin(), reconstructedTracksPerSector[sector].end(), std::greater<MUonERecoTrack3D>());	


		if(m_reconstructionConfiguration.maxNumberOfSharedHits() < 0) {
                if(m_reconstructionConfiguration.verbose())
			std::cout << "-1 HIT SHARED" << std::endl;
			//for each of the track pairs, remove shared hits from the worse track
			//and remove those tracks if at least one of their 2D tracks has 0 or 1 hit left
			//ex. start with track 1 and remove shared hits from tracks 2,3,4,..., go to track 2 and remove hits from tracks 3,4,5... and so on

			auto begin = reconstructedTracksPerSector[sector].begin();
			auto it = reconstructedTracksPerSector[sector].begin();
			auto end = reconstructedTracksPerSector[sector].end();

//  reconstructedTracksPerSector[sector].erase(std::remove_if(it + 1, end, [&](MUonERecoTrack3D& track) {return track.removeHitsSharedWith(*it);}), end);

			while(it != end) {

			auto& tr = *it;
			auto track_stereo_hits = tr.stereoHits();



if(m_reconstructionConfiguration.verbose()) std::cout << "------------REFERENCE TRACK------------" << std::endl;
 if(m_reconstructionConfiguration.verbose()){
                        for ( auto& hit : tr.xTrack().hits())
                                {
                                        std::cout << hit.moduleIndex() << " X module position " << hit.position() << std::endl;
                                }
                        for ( auto& hit : tr.yTrack().hits())
                                {
                                        std::cout << hit.moduleIndex() << " Y module position " << hit.position() << std::endl;

                                }
                        for ( auto& hit : tr.stereoHits())
                                {
                                        std::cout << hit.moduleIndex() << " Stereo module position " << hit.position() << std::endl;

                                }
 }

//remove XY if shared
reconstructedTracksPerSector[sector].erase(std::remove_if(it + 1, end, [&](MUonERecoTrack3D& track) {return track.removeXYHitsSharedWith(*it);}), end);

if(it==begin){


if(m_reconstructionConfiguration.verbose())
std::cout << "N tracks after removal " <<  reconstructedTracksPerSector[sector].size() << " tracks, reference track has " << tr.stereoHits().size() << " stereo" << std::endl;

std::vector<MUonERecoTrack3D> matches;
matches.reserve(reconstructedTracksPerSector[sector].size()-1);

std::copy_if(it + 1, reconstructedTracksPerSector[sector].end(), std::back_inserter(matches), [&track_stereo_hits](MUonERecoTrack3D& track){

                if(std::find_if(track_stereo_hits.begin(),track_stereo_hits.end(), [&](MUonERecoHit const& hit){

                        auto const& stereo_hits = track.stereoHits();
                        for(auto const& shits : stereo_hits)
                                if(hit == shits)
                                        return true;
                        return false;

                        })!=track_stereo_hits.end()) return true;
                else return false;

                });

if(matches.size()==0){

if(m_reconstructionConfiguration.verbose())
std::cout << "no stereo hit shared" << std::endl;
                                ++it;
                                end = reconstructedTracksPerSector[sector].end();
}
else{

                                for ( auto sh_tr : matches){
//auto& sh_tr=*sh;

                        auto hits_per_module_removed = hitsPerModulePerSector[sector];
                        Int_t nstereoU_window=0;
                        Int_t nstereoV_window=0;

                for(Int_t h=0; h<hits_per_module_removed.size(); h++){
                        for(Int_t hit=0; hit < hits_per_module_removed.at(h).size(); hit++){

                        if(modules[h].isStereo() and sh_tr.distanceToHit(hits_per_module_removed.at(h).at(hit))<=m_reconstructionConfiguration.uvHitAssignmentWindow() and hits_per_module_removed.at(h).at(hit).moduleIndex()==2) nstereoU_window++;
                        if(modules[h].isStereo() and sh_tr.distanceToHit(hits_per_module_removed.at(h).at(hit))<=m_reconstructionConfiguration.uvHitAssignmentWindow() and hits_per_module_removed.at(h).at(hit).moduleIndex()==3) nstereoV_window++;
                        }
                }
                if(m_reconstructionConfiguration.verbose()) std::cout << "nstereoU_window and nstereoV_window are " << nstereoU_window << ", " << nstereoV_window << std::endl;

                        if(nstereoU_window>1 or nstereoV_window>1){


if(m_reconstructionConfiguration.verbose()) std::cout << "stereo hits shared track " << sh_tr.stereoHits().size() << std::endl;


                        std::vector<MUonERecoHit> removed_hit;
                        removed_hit.reserve(2);
                        auto sh_track_stereo_hits = sh_tr.stereoHits();


               if(m_reconstructionConfiguration.verbose()){ std::cout << " before removal "
<<" N x " << sh_tr.xTrack().numberOfHits()
<<" N y " << sh_tr.yTrack().numberOfHits()
<<" N stereo " << sh_tr.numberOfStereoHits()
<< std::endl;}



                        for(Int_t h=0; h<sh_track_stereo_hits.size(); h++){

            if(m_reconstructionConfiguration.verbose()) std::cout << " inedex " << sh_track_stereo_hits.at(h).module().index() << std::endl;

                         if(std::find(track_stereo_hits.begin(),track_stereo_hits.end(), sh_track_stereo_hits.at(h))!=track_stereo_hits.end()){
                                    if(m_reconstructionConfiguration.verbose()) std::cout << "found"<<std::endl;
                                        if(sh_track_stereo_hits.at(h).module().index()==2 and nstereoU_window>1){removed_hit.push_back(sh_track_stereo_hits.at(h));

//              sh_tr.stereoHits().erase(std::find(sh_tr.stereoHits().begin(), sh_tr.stereoHits().end(), sh_track_stereo_hits.at(h)), sh_tr.stereoHits().end());
                sh_tr.removeStereoHitsIfEqual(sh_track_stereo_hits.at(h));

                if(m_reconstructionConfiguration.verbose()) std::cout << "index removed hit " << track_stereo_hits.at(h).module().index() << std::endl;}
                                        else if(sh_track_stereo_hits.at(h).module().index()==3 and nstereoV_window>1){removed_hit.push_back(sh_track_stereo_hits.at(h));

//                sh_tr.stereoHits().erase(std::find(sh_tr.stereoHits().begin(), sh_tr.stereoHits().end(), sh_track_stereo_hits.at(h)), sh_tr.stereoHits().end());

                sh_tr.removeStereoHitsIfEqual(sh_track_stereo_hits.at(h));

                                            if(m_reconstructionConfiguration.verbose()) std::cout << "index removed hit " << track_stereo_hits.at(h).module().index() << std::endl;}
                                }
                        }


               if(m_reconstructionConfiguration.verbose()){ std::cout << " after removal "
<<" N x " << sh_tr.xTrack().numberOfHits()
<<" N y " << sh_tr.yTrack().numberOfHits()
<<" N stereo " << sh_tr.numberOfStereoHits()              
<< std::endl;}

for(int r=0; r< removed_hit.size(); r++){
if(m_reconstructionConfiguration.verbose())
{std::cout << r << "th removed hit at module " << removed_hit.at(r).module().index() << " over " << removed_hit.size() << std::endl;
std::cout << "hitsPerModulePerSector " << hits_per_module_removed[removed_hit.at(r).module().index()].size() << std::endl;}

MUonERecoHit rhit=removed_hit.at(r);

hits_per_module_removed[removed_hit.at(r).module().index()].erase(std::remove_if(hits_per_module_removed[removed_hit.at(r).module().index()].begin(), 
 hits_per_module_removed[removed_hit.at(r).module().index()].end(), [&rhit](MUonERecoHit const& hit){

                                if(hit == rhit)
                                        return true;

                                else return false;
                        }), hits_per_module_removed[rhit.module().index()].end());


if(m_reconstructionConfiguration.verbose())
std::cout << "hitsPerModulePerSector after removal " << hits_per_module_removed[removed_hit.at(r).module().index()].size() << std::endl;
                        Int_t nstereo=hits_per_module_removed[removed_hit.at(r).module().index()].size();

if(m_reconstructionConfiguration.verbose())
std::cout << "how many stereo in module "<<removed_hit.at(r).module().index()<< "!? " << nstereo << std::endl;

                                if(nstereo > 0) {

                if(m_reconstructionConfiguration.verbose()) std::cout << "for loop " << std::endl;

                                        //assign stereo hits
                                        for(Int_t stubmodule_index = 0; stubmodule_index < hits_per_module_removed.size(); ++stubmodule_index) {

                if(m_reconstructionConfiguration.verbose()) std::cout << "for loop " << stubmodule_index << std::endl;


                                                auto const& module = modules[stubmodule_index];

                                                //only compatible stereo modules
                                                if(module.isStereo() and module.index()==removed_hit.at(r).module().index())
                                                        sh_tr.addClosestStereoHit(hits_per_module_removed[stubmodule_index], m_reconstructionConfiguration.uvHitAssignmentWindow());

                                        }//! end of hit assigning
                }//nstereo
        }// end for loop on removed hits


                                        if(!m_fitter.addFitInfo(sh_tr, 0 == sector ? m_reconstructionConfiguration.allowTrivialIncomingTracks() : false, m_reconstructionConfiguration.maxOutlierChi2()))
                                        {
               if(m_reconstructionConfiguration.verbose()) std::cout << "i am not trivial" << std::endl;
                                //reconstructedTracksPerSector[sector].emplace_back(std::move(sh_tr));



        auto element=std::find_if(it+1, reconstructedTracksPerSector[sector].end(), [&sh_tr](MUonERecoTrack3D& track){
                        if(track.index()==sh_tr.index())
                        return true;

                        return false;} );

                reconstructedTracksPerSector[sector].erase(element);
                reconstructedTracksPerSector[sector].insert(element,sh_tr);


if(m_reconstructionConfiguration.verbose()) std::cout << "------------MODIFIED TRACK------------" << std::endl;

                        for ( auto& hit : sh_tr.xTrack().hits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " X module position " << hit.position() << std::endl;
                                }
                        for ( auto& hit : sh_tr.yTrack().hits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " Y module position " << hit.position() << std::endl;

                                }
                        for ( auto& hit : sh_tr.stereoHits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " Stereo module position " << hit.position() << std::endl;

                                }


               if(m_reconstructionConfiguration.verbose()){ std::cout
<<" N x " << sh_tr.xTrack().numberOfHits()
<<" N y " << sh_tr.yTrack().numberOfHits()
<<" N stereo " << sh_tr.numberOfStereoHits()
<< std::endl;}

                                }else{

               if(m_reconstructionConfiguration.verbose()) std::cout << "i am trivial" << std::endl;


        auto element=std::find_if(it+1, reconstructedTracksPerSector[sector].end(), [&sh_tr](MUonERecoTrack3D& track){
                        if(track.index()==sh_tr.index())
                        return true;

                        return false;} );

                reconstructedTracksPerSector[sector].erase(element);
                reconstructedTracksPerSector[sector].insert(element,sh_tr);

                                        }
                        }//end nstereo_window
                else{

        if(m_reconstructionConfiguration.verbose()) std::cout << "not enough stereo hits " << std::endl;

        auto element=std::find_if(it+1, reconstructedTracksPerSector[sector].end(), [&sh_tr](MUonERecoTrack3D& track){
                        if(track.index()==sh_tr.index())
                        return true;

                        return false;} );

                reconstructedTracksPerSector[sector].erase(element);
                reconstructedTracksPerSector[sector].insert(element,sh_tr);

if(m_reconstructionConfiguration.verbose()) std::cout << "------------MODIFIED TRACK------------" << std::endl;

                        for ( auto& hit : sh_tr.xTrack().hits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " X module position " << hit.position() << std::endl;
                                }
                        for ( auto& hit : sh_tr.yTrack().hits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " Y module position " << hit.position() << std::endl;

                                }
                        for ( auto& hit : sh_tr.stereoHits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " Stereo module position " << hit.position() << std::endl;

                                }

               if(m_reconstructionConfiguration.verbose()){ std::cout
<<" N x " << sh_tr.xTrack().numberOfHits()
<<" N y " << sh_tr.yTrack().numberOfHits()
<<" N stereo " << sh_tr.numberOfStereoHits()
<< std::endl;}

        if(m_reconstructionConfiguration.verbose()) std::cout << "removed?"<<std::endl;
                }
                                                }//end for(sh_tr)
         if(m_reconstructionConfiguration.verbose()) std::cout << "end for(sh_tr)"<<std::endl;

                                ++it;
                                end = reconstructedTracksPerSector[sector].end();

                                        }//end if(sh!=end)
                                   }//if it==begin
else{
                 if(m_reconstructionConfiguration.verbose()) std::cout << "checking other tracks" << std::endl;


		    reconstructedTracksPerSector[sector].erase(std::remove_if(it + 1, reconstructedTracksPerSector[sector].end(), [&](MUonERecoTrack3D& track) {return track.StereoHitsSharedWith(*it);}), reconstructedTracksPerSector[sector].end());

                                ++it;
                                end = reconstructedTracksPerSector[sector].end();

					}
			}
		} else { //number of hit shared >0

			//for each of the track pairs, remove tracks which share more than N hits with better track
			auto begin = reconstructedTracksPerSector[sector].begin();
			auto it = reconstructedTracksPerSector[sector].begin();
			auto end = reconstructedTracksPerSector[sector].end();

			if(m_reconstructionConfiguration.restrictHitSharingToFirstNModulesInSector() >= 0) {

				while(it != end) {

if(m_reconstructionConfiguration.verbose()) std::cout << "------------REFERENCE TRACK------------" << std::endl;

                        auto& tr = *it;
                        auto track_stereo_hits = tr.stereoHits();


                        for ( auto& hit : tr.xTrack().hits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " X module position " << hit.position() << std::endl;
                                }
                        for ( auto& hit : tr.yTrack().hits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " Y module position " << hit.position() << std::endl;

                                }
                        for ( auto& hit : tr.stereoHits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " Stereo module position " << hit.position() << std::endl;

                                }

if(it==begin){

if(m_reconstructionConfiguration.verbose()) std::cout << "N tracks before removal allowing XY " <<  reconstructedTracksPerSector[sector].size() << " tracks"<< std::endl;
					reconstructedTracksPerSector[sector].erase(std::remove_if(it + 1, end, [&](MUonERecoTrack3D& track) {auto count = track.numberOfXYHitsSharedWith(*it, m_reconstructionConfiguration.restrictHitSharingToFirstNModulesInSector()); return (count.first > m_reconstructionConfiguration.maxNumberOfSharedHits() || count.second > 0);}), end);


if(m_reconstructionConfiguration.verbose())
std::cout << "N tracks after removal " <<  reconstructedTracksPerSector[sector].size() << " tracks, reference track has " << tr.stereoHits().size() << " stereo" << std::endl;

std::vector<MUonERecoTrack3D> matches;
matches.reserve(reconstructedTracksPerSector[sector].size()-1);

std::copy_if(it + 1, reconstructedTracksPerSector[sector].end(), std::back_inserter(matches), [&track_stereo_hits](MUonERecoTrack3D& track){

                if(std::find_if(track_stereo_hits.begin(),track_stereo_hits.end(), [&](MUonERecoHit const& hit){

                        auto const& stereo_hits = track.stereoHits();
                        for(auto const& shits : stereo_hits)
                                if(hit == shits)
                                        return true;
                        return false;

                        })!=track_stereo_hits.end()) return true;
                else return false;

                });


if(matches.size()==0){

if(m_reconstructionConfiguration.verbose())
std::cout << "no stereo hit shared" << std::endl;
                                ++it;
                                end = reconstructedTracksPerSector[sector].end();
}
else{

				for ( auto sh_tr : matches){
//auto& sh_tr=*sh;

                        auto hits_per_module_removed = hitsPerModulePerSector[sector];
			Int_t nstereoU_window=0;
			Int_t nstereoV_window=0;

		for(Int_t h=0; h<hits_per_module_removed.size(); h++){
			for(Int_t hit=0; hit < hits_per_module_removed.at(h).size(); hit++){

/*		if(m_reconstructionConfiguration.verbose()){
std::cout << "module " << h << " index " << hits_per_module_removed.at(h).at(hit).moduleIndex() << " distance " << sh_tr.distanceToHit(hits_per_module_removed.at(h).at(hit)) << std::endl;
std::cout << "UV window " << m_reconstructionConfiguration.uvHitAssignmentWindow() << std::endl;
}*/
			if(modules[h].isStereo() and sh_tr.distanceToHit(hits_per_module_removed.at(h).at(hit))<=m_reconstructionConfiguration.uvHitAssignmentWindow() and hits_per_module_removed.at(h).at(hit).moduleIndex()==2) nstereoU_window++;
			if(modules[h].isStereo() and sh_tr.distanceToHit(hits_per_module_removed.at(h).at(hit))<=m_reconstructionConfiguration.uvHitAssignmentWindow() and hits_per_module_removed.at(h).at(hit).moduleIndex()==3) nstereoV_window++;
			}
		}
                if(m_reconstructionConfiguration.verbose()) std::cout << "nstereoU_window and nstereoV_window are " << nstereoU_window << ", " << nstereoV_window << std::endl;

			if(nstereoU_window>1 or nstereoV_window>1){


if(m_reconstructionConfiguration.verbose()) std::cout << "stereo hits shared track " << sh_tr.stereoHits().size() << std::endl;


                        std::vector<MUonERecoHit> removed_hit;
                        removed_hit.reserve(2);
                        auto sh_track_stereo_hits = sh_tr.stereoHits();


               if(m_reconstructionConfiguration.verbose()){ std::cout << " before removal "
<<" N x " << sh_tr.xTrack().numberOfHits()
<<" N y " << sh_tr.yTrack().numberOfHits()
<<" N stereo " << sh_tr.numberOfStereoHits()
<< std::endl;}



                        for(Int_t h=0; h<sh_track_stereo_hits.size(); h++){

            if(m_reconstructionConfiguration.verbose()) std::cout << " inedex " << sh_track_stereo_hits.at(h).module().index() << std::endl;

                         if(std::find(track_stereo_hits.begin(),track_stereo_hits.end(), sh_track_stereo_hits.at(h))!=track_stereo_hits.end()){
			            if(m_reconstructionConfiguration.verbose()) std::cout << "found"<<std::endl;
					if(sh_track_stereo_hits.at(h).module().index()==2 and nstereoU_window>1){removed_hit.push_back(sh_track_stereo_hits.at(h));

//		sh_tr.stereoHits().erase(std::find(sh_tr.stereoHits().begin(), sh_tr.stereoHits().end(), sh_track_stereo_hits.at(h)), sh_tr.stereoHits().end());
                sh_tr.removeStereoHitsIfEqual(sh_track_stereo_hits.at(h));

                                                                                                                                        if(m_reconstructionConfiguration.verbose()) std::cout << "index removed hit " << track_stereo_hits.at(h).module().index() << std::endl;}

					else if(sh_track_stereo_hits.at(h).module().index()==3 and nstereoV_window>1){removed_hit.push_back(sh_track_stereo_hits.at(h));

//                sh_tr.stereoHits().erase(std::find(sh_tr.stereoHits().begin(), sh_tr.stereoHits().end(), sh_track_stereo_hits.at(h)), sh_tr.stereoHits().end());

		sh_tr.removeStereoHitsIfEqual(sh_track_stereo_hits.at(h));

                                                                                                                                        if(m_reconstructionConfiguration.verbose()) std::cout << "index removed hit " << track_stereo_hits.at(h).module().index() << std::endl;}
				}

                        }


               if(m_reconstructionConfiguration.verbose()){ std::cout << " after removal "
<<" N x " << sh_tr.xTrack().numberOfHits()
<<" N y " << sh_tr.yTrack().numberOfHits()
<<" N stereo " << sh_tr.numberOfStereoHits()              
<< std::endl;}

for(int r=0; r< removed_hit.size(); r++){
if(m_reconstructionConfiguration.verbose())
{std::cout << r << "th removed hit at module " << removed_hit.at(r).module().index() << " over " << removed_hit.size() << std::endl;
std::cout << "hitsPerModulePerSector " << hits_per_module_removed[removed_hit.at(r).module().index()].size() << std::endl;}

MUonERecoHit rhit=removed_hit.at(r);

hits_per_module_removed[removed_hit.at(r).module().index()].erase(std::remove_if(hits_per_module_removed[removed_hit.at(r).module().index()].begin(), 
 hits_per_module_removed[removed_hit.at(r).module().index()].end(), [&rhit](MUonERecoHit const& hit){

                                if(hit == rhit)
                                        return true;

                                else return false;
                        }), hits_per_module_removed[rhit.module().index()].end());


if(m_reconstructionConfiguration.verbose())
std::cout << "hitsPerModulePerSector after removal " << hits_per_module_removed[removed_hit.at(r).module().index()].size() << std::endl;

                        Int_t nstereo=hits_per_module_removed[removed_hit.at(r).module().index()].size();

if(m_reconstructionConfiguration.verbose())
std::cout << "how many stereo in module "<<removed_hit.at(r).module().index()<< "!? " << nstereo << std::endl;

                                if(nstereo > 0) {

                if(m_reconstructionConfiguration.verbose()) std::cout << "for loop " << std::endl;

                                        //assign stereo hits
                                        for(Int_t stubmodule_index = 0; stubmodule_index < hits_per_module_removed.size(); ++stubmodule_index) {

                if(m_reconstructionConfiguration.verbose()) std::cout << "for loop " << stubmodule_index << std::endl;


                                                auto const& module = modules[stubmodule_index];

                                                //only compatible stereo modules
                                                if(module.isStereo() and module.index()==removed_hit.at(r).module().index())
                                                        sh_tr.addClosestStereoHit(hits_per_module_removed[stubmodule_index], m_reconstructionConfiguration.uvHitAssignmentWindow());

                                        }//! end of hit assigning
                }//nstereo
        }// end for loop on removed hits


                                        if(!m_fitter.addFitInfo(sh_tr, 0 == sector ? m_reconstructionConfiguration.allowTrivialIncomingTracks() : false, m_reconstructionConfiguration.maxOutlierChi2()))
                                        {
               if(m_reconstructionConfiguration.verbose()) std::cout << "i am not trivial" << std::endl;
                                //reconstructedTracksPerSector[sector].emplace_back(std::move(sh_tr));

auto count = sh_tr.numberOfHitsSharedWith(*it, m_reconstructionConfiguration.restrictHitSharingToFirstNModulesInSector());

        if(m_reconstructionConfiguration.verbose()) std::cout << "numberOfHitsSharedWith " << count.first << std::endl;

	auto element=std::find_if(it+1, reconstructedTracksPerSector[sector].end(), [&sh_tr](MUonERecoTrack3D& track){
                        if(track.index()==sh_tr.index())
                        return true;

                        return false;} );

if(count.first > m_reconstructionConfiguration.maxNumberOfSharedHits() || count.second > 0){
		reconstructedTracksPerSector[sector].erase(element);}
else{
                reconstructedTracksPerSector[sector].erase(element);
		reconstructedTracksPerSector[sector].insert(element,sh_tr);

	}
if(m_reconstructionConfiguration.verbose()) std::cout << "------------MODIFIED TRACK------------" << std::endl;

                        for ( auto& hit : sh_tr.xTrack().hits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " X module position " << hit.position() << std::endl;
                                }
                        for ( auto& hit : sh_tr.yTrack().hits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " Y module position " << hit.position() << std::endl;

                                }
                        for ( auto& hit : sh_tr.stereoHits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " Stereo module position " << hit.position() << std::endl;

                                }


               if(m_reconstructionConfiguration.verbose()){ std::cout
<<" N x " << sh_tr.xTrack().numberOfHits()
<<" N y " << sh_tr.yTrack().numberOfHits()
<<" N stereo " << sh_tr.numberOfStereoHits()
<< std::endl;}
                                }else{
               if(m_reconstructionConfiguration.verbose()) std::cout << "i am trivial" << std::endl;
auto count = sh_tr.numberOfHitsSharedWith(*it, m_reconstructionConfiguration.restrictHitSharingToFirstNModulesInSector());


        if(m_reconstructionConfiguration.verbose()) std::cout << "numberOfHitsSharedWith " << count.first << std::endl;

        auto element=std::find_if(it+1, reconstructedTracksPerSector[sector].end(), [&sh_tr](MUonERecoTrack3D& track){
                        if(track.index()==sh_tr.index())
                        return true;

                        return false;} );

if(count.first > m_reconstructionConfiguration.maxNumberOfSharedHits() || count.second > 0){
                reconstructedTracksPerSector[sector].erase(element);}
else{
                reconstructedTracksPerSector[sector].erase(element);
                reconstructedTracksPerSector[sector].insert(element,sh_tr);    

        }
if(m_reconstructionConfiguration.verbose()) std::cout << "------------MODIFIED TRACK------------" << std::endl;

                        for ( auto& hit : sh_tr.xTrack().hits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " X module position " << hit.position() << std::endl;
                                }
                        for ( auto& hit : sh_tr.yTrack().hits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " Y module position " << hit.position() << std::endl;

                                }
                        for ( auto& hit : sh_tr.stereoHits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " Stereo module position " << hit.position() << std::endl;

                                }

               if(m_reconstructionConfiguration.verbose()){ std::cout
<<" N x " << sh_tr.xTrack().numberOfHits()
<<" N y " << sh_tr.yTrack().numberOfHits()
<<" N stereo " << sh_tr.numberOfStereoHits()
<< std::endl;}
                                }

			}//end nstereo_window
		else{
	if(m_reconstructionConfiguration.verbose()) std::cout << "not enough stereo hits " << std::endl;

auto count = sh_tr.numberOfHitsSharedWith(*it, m_reconstructionConfiguration.restrictHitSharingToFirstNModulesInSector());

        if(m_reconstructionConfiguration.verbose()) std::cout << "numberOfHitsSharedWith " << count.first << std::endl;

        auto element=std::find_if(it+1, reconstructedTracksPerSector[sector].end(), [&sh_tr](MUonERecoTrack3D& track){
                        if(track.index()==sh_tr.index())
                        return true;

                        return false;} );

if(count.first > m_reconstructionConfiguration.maxNumberOfSharedHits() || count.second > 0){
                reconstructedTracksPerSector[sector].erase(element);}
else{
                reconstructedTracksPerSector[sector].erase(element);
                reconstructedTracksPerSector[sector].insert(element,sh_tr);    

        }
if(m_reconstructionConfiguration.verbose()) std::cout << "------------MODIFIED TRACK------------" << std::endl;

                        for ( auto& hit : sh_tr.xTrack().hits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " X module position " << hit.position() << std::endl;
                                }
                        for ( auto& hit : sh_tr.yTrack().hits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " Y module position " << hit.position() << std::endl;

                                }
                        for ( auto& hit : sh_tr.stereoHits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " Stereo module position " << hit.position() << std::endl;

                                }

               if(m_reconstructionConfiguration.verbose()){ std::cout
<<" N x " << sh_tr.xTrack().numberOfHits()
<<" N y " << sh_tr.yTrack().numberOfHits()
<<" N stereo " << sh_tr.numberOfStereoHits()
<< std::endl;}

        if(m_reconstructionConfiguration.verbose()) std::cout << "removed?"<<std::endl;
		}
						}//end for(sh_tr)
         if(m_reconstructionConfiguration.verbose()) std::cout << "end for(sh_tr)"<<std::endl;

                                ++it;
                                end = reconstructedTracksPerSector[sector].end();

					}//end if(sh!=end)
				   }//if it==begin
else{
	         if(m_reconstructionConfiguration.verbose()) std::cout << "checking other tracks" << std::endl;

		reconstructedTracksPerSector[sector].erase(std::remove_if(it + 1, end, [&](MUonERecoTrack3D& track) {auto count = track.numberOfHitsSharedWith(*it, m_reconstructionConfiguration.restrictHitSharingToFirstNModulesInSector()); return (count.first > m_reconstructionConfiguration.maxNumberOfSharedHits() || count.second > 0);}), end);
                                ++it;
                                end = reconstructedTracksPerSector[sector].end();

	}
				}//end while it!=end

			} else {

				while(it != end) {
					reconstructedTracksPerSector[sector].erase(std::remove_if(it + 1, end, [&](MUonERecoTrack3D& track) {return (track.numberOfHitsSharedWith(*it)).second > m_reconstructionConfiguration.maxNumberOfSharedHits();}), end);
					++it;
					end = reconstructedTracksPerSector[sector].end();
				}
			}
		}


		//loop on remaining tracks and perform kalman fit and linking
		for(auto& tr : reconstructedTracksPerSector[sector]) {
if(m_reconstructionConfiguration.verbose()) std::cout << "------------FINAL TRACKS------------" << std::endl;

			for ( auto& hit : tr.xTrack().hits())
				{
					if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " X module position " << hit.position() << std::endl;
				}
                        for ( auto& hit : tr.yTrack().hits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " Y module position " << hit.position() << std::endl;

                                }
                        for ( auto& hit : tr.stereoHits())
                                {
                                        if(m_reconstructionConfiguration.verbose()) std::cout << hit.moduleIndex() << " Stereo module position " << hit.position() << std::endl;

                                }

			m_fitter.addFirstStateInfo(tr);

			if(m_reconstructionConfiguration.isMC())
				tr.linkToMC(m_mcTracks, m_reconstructionConfiguration.weightLinkedTracksByDepositedCharge(), m_reconstructionConfiguration.useSeedingSensorOnly());
		}


		if(m_reconstructionConfiguration.verbose()) 
			std::cout << reconstructedTracksPerSector[sector].size() << " track(s) left after filtering." << std::endl;

	} //end of reconstructing 3D tracks in all sectors


	if(std::all_of(reconstructedTracksPerSector.begin(), reconstructedTracksPerSector.end(), [](std::vector<MUonERecoTrack3D> const& vec){return vec.empty();}))
		return false;

		
	return true;
}


Bool_t MUonEReconstruction::reconstructVertices(std::vector<std::vector<std::vector<MUonERecoHit>>> const& hitsPerModulePerSector, std::vector<std::vector<MUonERecoTrack3D>> const& reconstructedTracksPerSector, std::vector<MUonERecoVertex>& reconstrucedVertices, std::vector<MUonERecoAdaptiveFitterVertex>& reconstrucedAFVertices, std::vector<MUonERecoGenericVertex>& reconstrucedGenericVertices) const {

	for(Int_t target_index = 0; target_index < m_targets.size(); ++target_index) {

		auto const& target = m_targets[target_index];
		auto const& tracks_before_target = reconstructedTracksPerSector[target_index];
		auto const& tracks_after_target = reconstructedTracksPerSector[target_index + 1];

		if(tracks_after_target.size() < 2)
			continue;

		//run adaptive fitter to reconstruct vertices with variable number of outgoing tracks
		//allowed to run even if there's no incoming track
		if(m_reconstructionConfiguration.runAdaptiveVertexFitter())
			m_fitter.runAdaptiveFitter(target, tracks_after_target, reconstrucedAFVertices, m_reconstructionConfiguration, !m_reconstructionConfiguration.disableAdaptiveFitterSeeding());		

		if(tracks_before_target.empty())
			continue;

		//reconstruct signal vertices using linear kinematic fit
		//single track going in, two going out
		for(auto const& incoming_track : tracks_before_target) {

			for(Int_t first_track_index = 0; first_track_index < tracks_after_target.size() - 1; ++first_track_index) {

				auto const& first_track = tracks_after_target[first_track_index];

				for(Int_t second_track_index = first_track_index + 1; second_track_index < tracks_after_target.size(); ++second_track_index) {
					
					auto const& second_track = tracks_after_target[second_track_index];

					auto vertex = MUonERecoVertex(incoming_track, first_track, second_track, target);
					Bool_t add_alternative_PID = m_reconstructionConfiguration.vertexSplittingAngleThreshold() > 0 && vertex.firstTheta() < m_reconstructionConfiguration.vertexSplittingAngleThreshold() && vertex.secondTheta() < m_reconstructionConfiguration.vertexSplittingAngleThreshold();
					
					auto fit_status = m_fitter.addFitInfo(target, m_modulesPerSector, vertex, m_reconstructionConfiguration);
					if(!fit_status && m_reconstructionConfiguration.reassignHitsAfterKinematicFit()) {

						if(vertex.reassignHits(hitsPerModulePerSector[target_index + 1], m_reconstructionConfiguration.maxNumberOfSharedHits(), m_reconstructionConfiguration.xyHitAssignmentWindow(), m_reconstructionConfiguration.uvHitAssignmentWindow()))
							fit_status = m_fitter.addFitInfo(target, m_modulesPerSector, vertex, m_reconstructionConfiguration);
					} 
					
					if(!fit_status) {

						if(m_reconstructionConfiguration.maxVertexChi2PerDegreeOfFreedom() > 0) {

							if(vertex.chi2PerDegreeOfFreedom() < m_reconstructionConfiguration.maxVertexChi2PerDegreeOfFreedom()) {

								reconstrucedVertices.emplace_back(std::move(vertex));
							}

						} else {

							reconstrucedVertices.emplace_back(std::move(vertex));
						}

					}


					if(add_alternative_PID) {

						auto vertex = MUonERecoVertex(incoming_track, first_track, second_track, target, true);
						
						auto fit_status = m_fitter.addFitInfo(target, m_modulesPerSector, vertex, m_reconstructionConfiguration);
						if(!fit_status && m_reconstructionConfiguration.reassignHitsAfterKinematicFit()) {

							if(vertex.reassignHits(hitsPerModulePerSector[target_index + 1], m_reconstructionConfiguration.maxNumberOfSharedHits(), m_reconstructionConfiguration.xyHitAssignmentWindow(), m_reconstructionConfiguration.uvHitAssignmentWindow()))
								fit_status = m_fitter.addFitInfo(target, m_modulesPerSector, vertex, m_reconstructionConfiguration);
						} 
						
						if(!fit_status) {

							if(m_reconstructionConfiguration.maxVertexChi2PerDegreeOfFreedom() > 0) {

								if(vertex.chi2PerDegreeOfFreedom() < m_reconstructionConfiguration.maxVertexChi2PerDegreeOfFreedom()) {

									reconstrucedVertices.emplace_back(std::move(vertex));
								}

							} else {

								reconstrucedVertices.emplace_back(std::move(vertex));
							}
							
						}						
					}

				}
			}
		} //!loop on incoming tracks



		//generic vertices
		if(m_reconstructionConfiguration.maxGenericVertexOutgoingTrackMultiplicity() > 2) {


			//max 10 multiplicity supported with SMatrices for now, if needed a larger one or specific multiplicities can be added in the future 
			Int_t N = tracks_after_target.size();
			std::vector<MUonERecoTrack3D> out_tracks;

			auto generateGenericVertices = [&](Int_t offset, Int_t k, Int_t current_track, auto&& ggv) {

				if(0 == k) {

					//a set of outgoing tracks found
					for(auto const& incoming_track : tracks_before_target) {

						reconstrucedGenericVertices.emplace_back(incoming_track, out_tracks, target);
					}

					return;
				}

				for(Int_t i = offset; i <= N - k; ++i) {

					out_tracks[current_track] = tracks_after_target[i];

					ggv(i+1, k-1, current_track+1, ggv);
				}
			};


			Int_t maxMult = std::min(std::min(m_reconstructionConfiguration.maxGenericVertexOutgoingTrackMultiplicity(), 10), N);

			Int_t mult_start = 0;

			switch(maxMult) {

			case 10:

				out_tracks.clear();
				out_tracks.resize(10);
				generateGenericVertices(0, 10, 0, generateGenericVertices);
				for(Int_t i = mult_start; i < reconstrucedGenericVertices.size(); ++i)
					m_fitter.addFitInfo<10>(target, m_modulesPerSector, reconstrucedGenericVertices[i], m_reconstructionConfiguration);
				mult_start = reconstrucedGenericVertices.size();

			case 9:

				out_tracks.clear();
				out_tracks.resize(9);
				generateGenericVertices(0, 9, 0, generateGenericVertices);
				for(Int_t i = mult_start; i < reconstrucedGenericVertices.size(); ++i)
					m_fitter.addFitInfo<9>(target, m_modulesPerSector, reconstrucedGenericVertices[i], m_reconstructionConfiguration);
				mult_start = reconstrucedGenericVertices.size();

			case 8:

				out_tracks.clear();
				out_tracks.resize(8);
				generateGenericVertices(0, 8, 0, generateGenericVertices);
				for(Int_t i = mult_start; i < reconstrucedGenericVertices.size(); ++i)
					m_fitter.addFitInfo<8>(target, m_modulesPerSector, reconstrucedGenericVertices[i], m_reconstructionConfiguration);
				mult_start = reconstrucedGenericVertices.size();

			case 7:

				out_tracks.clear();
				out_tracks.resize(7);
				generateGenericVertices(0, 7, 0, generateGenericVertices);
				for(Int_t i = mult_start; i < reconstrucedGenericVertices.size(); ++i)
					m_fitter.addFitInfo<7>(target, m_modulesPerSector, reconstrucedGenericVertices[i], m_reconstructionConfiguration);
				mult_start = reconstrucedGenericVertices.size();

			case 6:

				out_tracks.clear();
				out_tracks.resize(6);
				generateGenericVertices(0, 6, 0, generateGenericVertices);
				for(Int_t i = mult_start; i < reconstrucedGenericVertices.size(); ++i)
					m_fitter.addFitInfo<6>(target, m_modulesPerSector, reconstrucedGenericVertices[i], m_reconstructionConfiguration);
				mult_start = reconstrucedGenericVertices.size();

			case 5:

				out_tracks.clear();
				out_tracks.resize(5);
				generateGenericVertices(0, 5, 0, generateGenericVertices);
				for(Int_t i = mult_start; i < reconstrucedGenericVertices.size(); ++i)
					m_fitter.addFitInfo<5>(target, m_modulesPerSector, reconstrucedGenericVertices[i], m_reconstructionConfiguration);
				mult_start = reconstrucedGenericVertices.size();

			case 4:

				out_tracks.clear();
				out_tracks.resize(4);
				generateGenericVertices(0, 4, 0, generateGenericVertices);
				for(Int_t i = mult_start; i < reconstrucedGenericVertices.size(); ++i)
					m_fitter.addFitInfo<4>(target, m_modulesPerSector, reconstrucedGenericVertices[i], m_reconstructionConfiguration);
				mult_start = reconstrucedGenericVertices.size();

			case 3:

				out_tracks.clear();
				out_tracks.resize(3);
				generateGenericVertices(0, 3, 0, generateGenericVertices);
				for(Int_t i = mult_start; i < reconstrucedGenericVertices.size(); ++i)
					m_fitter.addFitInfo<3>(target, m_modulesPerSector, reconstrucedGenericVertices[i], m_reconstructionConfiguration);
				mult_start = reconstrucedGenericVertices.size();

				break;

			default:
				break;
			}

			if(m_reconstructionConfiguration.maxGenericVertexChi2PerDegreeOfFreedom() > 0) {

				reconstrucedGenericVertices.erase(std::remove_if(reconstrucedGenericVertices.begin(), reconstrucedGenericVertices.end(), [&](MUonERecoGenericVertex& vtx) {return vtx.chi2PerDegreeOfFreedom() > m_reconstructionConfiguration.maxGenericVertexChi2PerDegreeOfFreedom();}), reconstrucedGenericVertices.end());
			}
			
		}

	}

	if(m_reconstructionConfiguration.verbose()) {
		std::cout << "Found " << reconstrucedVertices.size() << " vertices." << std::endl;

		if(m_reconstructionConfiguration.runAdaptiveVertexFitter())
			std::cout << "Found " << reconstrucedAFVertices.size() << " vertices from adaptive fitter." << std::endl;
	}

	if(!reconstrucedAFVertices.empty())
		std::sort(reconstrucedAFVertices.begin(), reconstrucedAFVertices.end(), std::greater<MUonERecoAdaptiveFitterVertex>());

	if(!reconstrucedGenericVertices.empty())
		std::sort(reconstrucedGenericVertices.begin(), reconstrucedGenericVertices.end(), std::greater<MUonERecoGenericVertex>());




	if(reconstrucedVertices.empty())
		return false;	

	//sort vertices from best to worst, see operator definitions in MUonERecoVertex.h
	std::sort(reconstrucedVertices.begin(), reconstrucedVertices.end(), std::greater<MUonERecoVertex>());	

	return true;

}


void MUonEReconstruction::Finish()
{
	if(m_output)
		delete m_output;
}

void MUonEReconstruction::Register()
{
	FairRootManager::Instance()->RegisterAny("ReconstructionOutput", m_output, true);
}


ClassImp(MUonEReconstruction)


