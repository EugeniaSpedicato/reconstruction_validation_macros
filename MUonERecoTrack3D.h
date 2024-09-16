#ifndef MUONERECOTRACK3D_H
#define MUONERECOTRACK3D_H

#include "Rtypes.h"
#include "TMath.h"
#include "TVector3.h"
#include "TClonesArray.h"

#include "Math/SMatrix.h"

#include <vector>
#include <utility>

#include "MUonETrack.h"
#include "MUonERecoHit.h"
#include "MUonERecoTrack2D.h"

/*
	Used to store 3D tracks (x- and y-projection 2D tracks combined) with their assigned hits and hits from stereo modules.
*/

class MUonERecoTrack3D {

public:

	MUonERecoTrack3D() = default;


	//creates 3D track from a pair of x- and y-projection 2D tracks
	MUonERecoTrack3D(MUonERecoTrack2D const& xTrack, MUonERecoTrack2D const& yTrack, Int_t sector, Int_t index)
		: m_xTrack(xTrack), m_yTrack(yTrack), m_sector(sector), m_index(index),
		m_chi2(xTrack.chi2() + yTrack.chi2())
		{}

	Int_t sector() const {return m_sector;}

	Int_t index() const {return m_index;}

	MUonERecoTrack2D const& xTrack() const {return m_xTrack;}
	MUonERecoTrack2D xTrackCopy() const {return m_xTrack;}

	MUonERecoTrack2D const& yTrack() const {return m_yTrack;}
	MUonERecoTrack2D yTrackCopy() const {return m_yTrack;}


	//returns a copy of all track hits
	std::vector<MUonERecoHit> hitsCopy() const {

		std::vector<MUonERecoHit> hits;
		hits.reserve(m_xTrack.numberOfHits() + m_yTrack.numberOfHits() + numberOfStereoHits());

		auto const& x_hits = m_xTrack.hits();
		auto const& y_hits = m_yTrack.hits();
		auto const& stereo_hits = m_stereoHits;

		hits.insert(hits.end(), x_hits.begin(), x_hits.end());
		hits.insert(hits.end(), y_hits.begin(), y_hits.end());
		hits.insert(hits.end(), stereo_hits.begin(), stereo_hits.end());

		return hits;	
	}

	//for access to the rest of hits see xTrack and yTrack and their appropriate methods
	std::vector<MUonERecoHit> const& stereoHits() const {return m_stereoHits;};
	std::vector<MUonERecoHit> stereoHitsCopy() const {return m_stereoHits;};


	Int_t numberOfStereoHits() const {return m_stereoHits.size();}
	Int_t numberOfHits() const {return m_xTrack.numberOfHits() + m_yTrack.numberOfHits() + numberOfStereoHits();}		

	Double_t chi2() const {return m_chi2;};

	Int_t degreesOfFreedom() const {return m_xTrack.degreesOfFreedom() + m_yTrack.degreesOfFreedom() + numberOfStereoHits();} //no new parameters added going from 2D to 3D track
	Double_t chi2PerDegreeOfFreedom() const {return m_chi2/degreesOfFreedom();}

	Bool_t kalmanFitSuccessful() const {return m_kalmanFitSuccessful;}

	Double_t firstStateZ() const {return m_firstStateZ;}
	Double_t firstStateX() const {return m_firstStateX;}
	Double_t firstStateY() const {return m_firstStateY;}

	Double_t firstStateXSlope() const {return m_firstStateXSlope;}
	Double_t firstStateYSlope() const {return m_firstStateYSlope;}

	ROOT::Math::SMatrix<Double_t, 4> const& firstStateCovarianceMatrix() const {return m_firstStateCovarianceMatrix;}	
	
	ROOT::Math::SMatrix<Double_t, 4> const& linearFitDerivativeMatrix() const {return m_linearFitDerivativeMatrix;}		
	ROOT::Math::SMatrix<Double_t, 4> const& linearFitCovarianceMatrix() const {return m_linearFitCovarianceMatrix;}	

	//defined for convenience
	Double_t x(Double_t z) const {return m_xTrack.valueAtZ(z);}
	Double_t xErrorAtZ(Double_t z) const {return m_xTrack.errorAtZ(z);}

	Double_t y(Double_t z) const {return m_yTrack.valueAtZ(z);}
	Double_t yErrorAtZ(Double_t z) const {return m_yTrack.errorAtZ(z);}


	TVector3 directionVector() const {return TVector3(m_xTrack.slope(), m_yTrack.slope(), 1.0).Unit();}
	TVector3 directionVectorFirstState() const {return TVector3(m_firstStateXSlope, m_firstStateYSlope, 1.0).Unit();}

	Bool_t hitsReassigned() const {return m_hitsReassigned;}
	

	//distance from track to hit in plane perpendicular to the Z axis (see MUonERecoHit.h for definition of perpendicular position)
	Double_t distanceToHit(const MUonERecoHit& hit) const {return fabs(distanceToHitSigned(hit));}
	Double_t distanceToHitSigned(const MUonERecoHit& hit) const {return hit.module().localCoordinatePerpendicular(x(hit.z()), y(hit.z())) - hit.positionPerpendicular();}


	//methods for track creation and filtering


	//add closest hit within distance threshold from given stereo module
	//returns true if hit was found and added
	Int_t findClosestStereoHit(const std::vector<MUonERecoHit>& hits_from_module, Double_t threshold) {

		Int_t closest_hit_index = -1;
		Double_t closest_hit_distance = 9999999999999;

		for(Int_t hit_index = 0; hit_index < hits_from_module.size(); ++hit_index) {

			auto const& hit = hits_from_module[hit_index];
			Double_t distance_from_track = distanceToHit(hit);

			if(distance_from_track < threshold && distance_from_track < closest_hit_distance) {

				closest_hit_distance = distance_from_track;
				closest_hit_index = hit_index;
			}
		}

		return closest_hit_index;
	}


	bool addClosestStereoHit(const std::vector<MUonERecoHit>& hits_from_module, Double_t threshold) {

		//find closest hit within threshold
		Int_t closest_hit_index = findClosestStereoHit(hits_from_module, threshold);


		if(-1 != closest_hit_index) {

			m_stereoHits.emplace_back(hits_from_module[closest_hit_index]);
			return true;
		}

		return false;
	}

	//used for hit reassignment after the kinematic fit
	Bool_t reassignClosestStereoHit(const std::vector<MUonERecoHit>& hitsInModule, Int_t hitIndex, Double_t threshold) {

		Int_t closest_hit_index = findClosestStereoHit(hitsInModule, threshold);

		if(-1 != closest_hit_index && (m_stereoHits[hitIndex] != hitsInModule[closest_hit_index])) {

			m_stereoHits[hitIndex] = hitsInModule[closest_hit_index];
			return true;
		}

		return false;
	}

	Bool_t reassignHit(std::vector<MUonERecoHit> const& hitsInModule, MUonERecoModule const& mod, Int_t hitIndex, Double_t xyThreshold, Double_t uvThreshold) {

		m_hitsReassigned = false;

		if (mod.isStereo()) {

			m_hitsReassigned = m_hitsReassigned | reassignClosestStereoHit(hitsInModule, hitIndex, uvThreshold);
		}
		else if ('x' == mod.projection()) {

			m_hitsReassigned = m_hitsReassigned | m_xTrack.reassignClosestHit(hitsInModule, hitIndex, xyThreshold);
		}
		else if ('y' == mod.projection()) {

			m_hitsReassigned = m_hitsReassigned | m_yTrack.reassignClosestHit(hitsInModule, hitIndex, xyThreshold);
		}

		return m_hitsReassigned;
	}


	//remove hits shared with given track from m_hits collection
	//returns true if any of the projection tracks is left with only one hit, so that the track can be removed
	//see MUonERecoTrack2D.h for details on the methods below
	Bool_t removeHitsSharedWith(MUonERecoTrack3D const& track) {

		if(m_xTrack.removeHitsSharedWith(track.xTrack()))
			return true;

		if(m_yTrack.removeHitsSharedWith(track.yTrack()))
			return true;

		/*auto const& track_hits = track.stereoHits();
		m_stereoHits.erase(std::remove_if(m_stereoHits.begin(), m_stereoHits.end(), [&track_hits](MUonERecoHit const& hit){

			for(auto const& track_hit : track_hits)
				if(hit == track_hit)
					return true;

			return false;
		}), m_stereoHits.end());*/

		return false;
	}

        Bool_t StereoHitsSharedWith(MUonERecoTrack3D const& track) {

                auto const& track_hits = track.stereoHits();

                m_stereoHits.erase(std::remove_if(m_stereoHits.begin(), m_stereoHits.end(), [&track_hits](MUonERecoHit const& hit){


                    for(auto const& track_hit : track_hits)
                                if(hit == track_hit)
                                        return true;

                        return false;
                }), m_stereoHits.end());

                return false;
        }

        Bool_t removeStereoHitsIfEqual(MUonERecoHit const& st_hit) {

                m_stereoHits.erase(std::find(m_stereoHits.begin(), m_stereoHits.end(), st_hit)), m_stereoHits.end();

		return true;
        }




	std::pair<Int_t, Int_t> numberOfHitsSharedWith(MUonERecoTrack3D const& track, Int_t restrictHitSharingToFirstNModulesInSector = -1) {

		auto x_count = m_xTrack.numberOfHitsSharedWith(track.xTrack(), restrictHitSharingToFirstNModulesInSector);
		auto y_count = m_yTrack.numberOfHitsSharedWith(track.yTrack(), restrictHitSharingToFirstNModulesInSector);


		//count stereo hits
		Int_t below_count = 0;
		Int_t above_count = 0;

		auto const& track_hits = track.stereoHits();

		for(auto const& hit : m_stereoHits) {

			for(auto const& track_hit : track_hits) {

				if(hit == track_hit) {

					hit.moduleIndex() < restrictHitSharingToFirstNModulesInSector ? ++below_count : ++above_count;
					break;
				}
			}
		}

		return std::make_pair(x_count.first + y_count.first + below_count, x_count.second + y_count.second + above_count); 
	}

        std::pair<Int_t, Int_t> numberOfXYHitsSharedWith(MUonERecoTrack3D const& track, Int_t restrictHitSharingToFirstNModulesInSector = -1) {

                auto x_count = m_xTrack.numberOfHitsSharedWith(track.xTrack(), restrictHitSharingToFirstNModulesInSector);
                auto y_count = m_yTrack.numberOfHitsSharedWith(track.yTrack(), restrictHitSharingToFirstNModulesInSector);



                return std::make_pair(x_count.first + y_count.first, x_count.second + y_count.second); 
        }


	std::pair<Int_t, Int_t> numberOfStereoHitsSharedWith(MUonERecoTrack3D const& track, Int_t restrictHitSharingToFirstNModulesInSector = -1) {


                //count stereo hits
                Int_t below_count = 0;
                Int_t above_count = 0;

                auto const& track_hits = track.stereoHits();

                for(auto const& hit : m_stereoHits) {

                        for(auto const& track_hit : track_hits) {

                                if(hit == track_hit) {

                                        hit.moduleIndex() < restrictHitSharingToFirstNModulesInSector ? ++below_count : ++above_count;
                                        break;
                                }
                        }
                }

                return std::make_pair(below_count,above_count); 
        }



	//methods for fitter

	//call this after new track parameters have been set or hits vectors have been modified
	Double_t recalculateChi2() {

		//xy hits have to be recalculated to account for alignment
		//m_chi2 = m_xTrack.recalculateChi2() + m_yTrack.recalculateChi2();

		m_chi2 = 0;

		for(auto const& hit : m_xTrack.hits()){
			Double_t chi2_sqrt = distanceToHit(hit) / hit.positionError();
			m_chi2 += chi2_sqrt * chi2_sqrt;
		}

		for(auto const& hit : m_yTrack.hits()){
			Double_t chi2_sqrt = distanceToHit(hit) / hit.positionError();
			m_chi2 += chi2_sqrt * chi2_sqrt;
		}

		for(auto const& hit : m_stereoHits){
			Double_t chi2_sqrt = distanceToHit(hit) / hit.positionError();
			m_chi2 += chi2_sqrt * chi2_sqrt;
		}

		return m_chi2;
	}


	void setXTrackParameters(Double_t slope, Double_t slope_error, Double_t x0, Double_t x0_error) {

		m_xTrack.setSlope(slope, slope_error);
		m_xTrack.setX0(x0, x0_error);
	}

	void setYTrackParameters(Double_t slope, Double_t slope_error, Double_t y0, Double_t y0_error) {

		m_yTrack.setSlope(slope, slope_error);
		m_yTrack.setY0(y0, y0_error);
	}

	void setLinearFitDerivativeMatrix(ROOT::Math::SMatrix<Double_t, 4> const& m) {

		 m_linearFitDerivativeMatrix = m;
	}

	void setLinearFitCovarianceMatrix(ROOT::Math::SMatrix<Double_t, 4> const& m) {

		 m_linearFitCovarianceMatrix = m;
	}

	void setFirstState(Bool_t successful, Double_t z, Double_t x, Double_t y, Double_t xSlope, Double_t ySlope, ROOT::Math::SMatrix<Double_t, 4> const& cov) {

		m_kalmanFitSuccessful = successful;

		m_firstStateZ = z;
		m_firstStateX = x;
		m_firstStateY = y;

		m_firstStateXSlope = xSlope;
		m_firstStateYSlope = ySlope;

		m_firstStateCovarianceMatrix = cov;
	}

	void setZ0(Double_t z0) {

		m_xTrack.setZ0(z0);
		m_yTrack.setZ0(z0);
	}	


	//used by fitter if any outliers were removed
	void setHits(std::vector<MUonERecoHit>& hits) {

		//split into x and y hits for appropriate 2D tracks
		std::vector<MUonERecoHit> x_hits;
		std::vector<MUonERecoHit> y_hits;
		
		x_hits.reserve(m_xTrack.numberOfHits());
		y_hits.reserve(m_yTrack.numberOfHits());
		
		//stereo hits can be added directly
		m_stereoHits.clear();

		for(auto& hit : hits) {

			if('x' == hit.module().projection())
				x_hits.emplace_back(hit);

			else if('y' == hit.module().projection())
				y_hits.emplace_back(hit);

			else
				m_stereoHits.emplace_back(hit);			
		}

		m_xTrack.setHits(x_hits);
		m_yTrack.setHits(y_hits);
	}


    Int_t linkedTrackID() const {return m_linkedTrackID;}
    Double_t fractionOfHitsSharedWithLinkedTrack() const {return m_fractionOfHitsSharedWithLinkedTrack;}
    Int_t processIDofLinkedTrack() const {return m_processIDofLinkedTrack;}


	void linkToMC(const TClonesArray* mcTracks, Bool_t weightLinkedTracks, Bool_t useSeedingSensorOnly) {

        //linked tracks are sorted by their charge contribution to both stub clusters, so weights are applied
        //and track with highest sum is chosen

        std::unordered_map<Int_t, Int_t> id_map; //key - track ID, value - index in vectors below
        std::vector<Double_t> track_sums; //sum of weights for all hits for a given track ID
        track_sums.reserve(20);
            
		std::vector<Int_t> track_hit_count; //how many assigned hits had this track linked (track ids don't repeat within a single hit) 
        track_hit_count.reserve(20);	

		Int_t sensorsCount = 0;

        auto processHits = [&id_map, &track_sums, &track_hit_count, &weightLinkedTracks, &useSeedingSensorOnly, &sensorsCount](const std::vector<MUonERecoHit>& hits) {

            for (auto const& hit : hits) { //loop on all hits assigned to track

                auto processDepositsFromSensor = [&id_map, &track_sums, &track_hit_count, &weightLinkedTracks](const std::vector<Int_t> ids) {
                        
                    for(Int_t track_index = 0; track_index < ids.size(); ++track_index) { //loop on all track ids linked to hit, in one sensor

                        Int_t track_id = ids[track_index];
                        Double_t track_weight = weightLinkedTracks ? 1/(track_index + 1) : 1; //first track - 1, second - 0.5 etc. 

                        //try to map this track ID to last index in vectors
                        std::pair<std::unordered_map<Int_t, Int_t>::iterator, bool> tmp_ret = id_map.insert({track_id, track_sums.size()});

                        if(!tmp_ret.second) { //track ID already present in the map

                            Int_t vector_index = tmp_ret.first->second;

                            track_sums[vector_index] += track_weight;
                            ++track_hit_count[vector_index];

                        } else { //track ID not yet considered


                            track_sums.push_back(track_weight);
                            track_hit_count.push_back(1);
                        }
                    }
                };

                processDepositsFromSensor(hit.seedingClusterLinkedTrackIds());
				++sensorsCount;

                if(hit.module().twoSensorsPerModule() && !useSeedingSensorOnly) {

                    processDepositsFromSensor(hit.correlationClusterLinkedTrackIds());
					++sensorsCount;
				}
            }
        };	

        processHits(xTrack().hits());
        processHits(yTrack().hits());
        processHits(stereoHits());	

        //find track ID with highest weight
        Double_t max_weight = -9999;
        Int_t max_weight_track_id = -9999;
        Int_t max_weight_track_hit_count = -9999;

        for(auto id_pair : id_map) {

            auto const& track_id = id_pair.first;
            auto const& vector_index = id_pair.second;
            auto const& current_track_weight = track_sums[vector_index];
            auto const& current_track_hit_count = track_hit_count[vector_index];

            if(current_track_weight > max_weight) {

                max_weight = current_track_weight;
                max_weight_track_id = track_id;
                max_weight_track_hit_count = current_track_hit_count;
            }
        }
           
        m_linkedTrackID = max_weight_track_id;
        m_fractionOfHitsSharedWithLinkedTrack = static_cast<Double_t>(max_weight_track_hit_count) / static_cast<Double_t>(sensorsCount);

        m_processIDofLinkedTrack = (static_cast<MUonETrack*>(mcTracks->At(m_linkedTrackID)))->interactionID();	
	}



private:

	Int_t m_sector{-1};

	Int_t m_index{-1};

	MUonERecoTrack2D m_xTrack;
	MUonERecoTrack2D m_yTrack;	

	std::vector<MUonERecoHit> m_stereoHits;

	Double_t m_chi2{0};

	Bool_t m_kalmanFitSuccessful{false};
	Double_t m_firstStateZ{0};
	Double_t m_firstStateX{0};
	Double_t m_firstStateY{0};

	Double_t m_firstStateXSlope{0};
	Double_t m_firstStateYSlope{0};

	ROOT::Math::SMatrix<Double_t, 4> m_firstStateCovarianceMatrix;

	ROOT::Math::SMatrix<Double_t, 4> m_linearFitDerivativeMatrix;
	ROOT::Math::SMatrix<Double_t, 4> m_linearFitCovarianceMatrix;

    Int_t m_linkedTrackID{-1};
    Double_t m_fractionOfHitsSharedWithLinkedTrack{0};
    Int_t m_processIDofLinkedTrack{-1};

	Bool_t m_hitsReassigned{false};

public:


	//if sorting is required, prefer tracks with more hits, more stereo hits and lower chi2
	//note that after the first if, both tracks have the same number of hits, so it's enough to compare chi2 and not chi2 per degree of freedom
	//worse track is considered to be of lower value, that is worse < better is true
    friend inline bool operator< (MUonERecoTrack3D const& lhs, MUonERecoTrack3D const& rhs) {
        
		if(lhs.numberOfHits() != rhs.numberOfHits())
			return lhs.numberOfHits() < rhs.numberOfHits();

		else if(lhs.numberOfStereoHits() != rhs.numberOfStereoHits())
			return lhs.numberOfStereoHits() < rhs.numberOfStereoHits();

		else
        	return lhs.chi2() > rhs.chi2();
    }
    friend inline bool operator> (MUonERecoTrack3D const& lhs, MUonERecoTrack3D const& rhs) {return  operator< (rhs,lhs);}
    friend inline bool operator<=(MUonERecoTrack3D const& lhs, MUonERecoTrack3D const& rhs) {return !operator> (lhs,rhs);}
    friend inline bool operator>=(MUonERecoTrack3D const& lhs, MUonERecoTrack3D const& rhs) {return !operator< (lhs,rhs);}


    ClassDef(MUonERecoTrack3D,2)

};

#endif//MUONERECOTRACK3D_H

