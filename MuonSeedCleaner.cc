/**
 *  See header file for a description of this class.
 *  
 *  \author Shih-Chuan Kao, Dominique Fortin - UCR
 */

#include <RecoMuon/MuonSeedGenerator/src/MuonSeedCleaner.h>

// Data Formats 
#include <DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h>
#include <DataFormats/TrajectorySeed/interface/TrajectorySeed.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>
#include <DataFormats/CSCRecHit/interface/CSCSegment.h>
#include <DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h>
#include <DataFormats/DTRecHit/interface/DTRecSegment4D.h>

// Geometry
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include <TrackingTools/DetLayers/interface/DetLayer.h>
#include <RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h>
#include <RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h>
#include <RecoMuon/Records/interface/MuonRecoGeometryRecord.h>

// muon service
#include <TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h>
#include <DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h>

// Framework
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include "FWCore/Framework/interface/Event.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 
#include <DataFormats/Common/interface/Handle.h>

// C++
#include <vector>
#include <deque>
#include <utility>

//typedef std::pair<double, TrajectorySeed> seedpr ;
//static bool ptDecreasing(const seedpr s1, const seedpr s2) { return ( s1.first > s2.first ); }
static bool lengthSorting(const TrajectorySeed s1, const TrajectorySeed s2) { return ( s1.nHits() > s2.nHits() ); }

/*
 * Constructor
 */
MuonSeedCleaner::MuonSeedCleaner(const edm::ParameterSet& pset){

  // Local Debug flag
  debug                = pset.getParameter<bool>("DebugMuonSeed");

  // muon service
  edm::ParameterSet serviceParameters = pset.getParameter<edm::ParameterSet>("ServiceParameters");
  theService        = new MuonServiceProxy(serviceParameters);

}

/*
 * Destructor
 */
MuonSeedCleaner::~MuonSeedCleaner(){

  if (theService) delete theService;
}

/* 
 * Cleaner
 *
 */

std::vector<TrajectorySeed> MuonSeedCleaner::seedCleaner(const edm::EventSetup& eventSetup, std::vector<TrajectorySeed>& seeds ) {

  theService->update(eventSetup);
  
  std::vector<TrajectorySeed> FinalSeeds;
  
  // group the seeds
  std::vector<SeedContainer> theCollection = GroupSeeds(seeds);

  // ckeck each group and pick the good one
  for (size_t i=0; i< theCollection.size(); i++ ) {

      // pick the seed w/ more than 1 segments and w/ 1st layer segment information
      SeedContainer goodSeeds   = SeedCandidates( theCollection[i], true );
      SeedContainer otherSeeds  = SeedCandidates( theCollection[i], false );
      if ( MomentumFilter( goodSeeds ) )  {
         std::vector<TrajectorySeed> betterSeeds = LengthFilter( goodSeeds ); 
	 TrajectorySeed bestSeed   = BetterChi2( betterSeeds );
	 //TrajectorySeed bestSeed   = BetterDirection( betterSeeds );
	 FinalSeeds.push_back( bestSeed );
      } 
      else if( MomentumFilter( otherSeeds ) ) {
         //std::cout<<" == type2 "<<std::endl;
         SeedContainer betterSeeds = LengthFilter( otherSeeds ); 
	 TrajectorySeed bestSeed   = BetterChi2( betterSeeds );
	 //TrajectorySeed bestSeed   = BetterDirection( betterSeeds );
	 FinalSeeds.push_back( bestSeed );
      } 
      else {
         //std::cout<<" == type3 "<<std::endl;
         SeedContainer betterSeeds = LengthFilter( theCollection[i] ); 
	 TrajectorySeed bestSeed   = BetterChi2( betterSeeds );
	 //TrajectorySeed bestSeed   = BetterDirection( betterSeeds );
	 FinalSeeds.push_back( bestSeed );
      }
  }  
  return FinalSeeds ; 

}


TrajectorySeed MuonSeedCleaner::BetterDirection(std::vector<TrajectorySeed>& seeds ) {

  float bestProjErr  = 9999.; 
  int winner = 0 ;
  AlgebraicSymMatrix mat(5,0) ; 
  for ( size_t i = 0; i < seeds.size(); i++ ) {

      edm::OwnVector<TrackingRecHit>::const_iterator r1 = seeds[i].recHits().first ;
      mat = r1->parametersError().similarityT( r1->projectionMatrix() );
      float ddx = mat[1][1]; 
      float ddy = mat[2][2]; 
      float dxx = mat[3][3]; 
      float dyy = mat[4][4];
      float projectErr = sqrt( (ddx*10000.) + (ddy*10000.) + dxx + dyy ) ;

      if ( projectErr > bestProjErr ) continue;
      winner = static_cast<int>(i) ;
      bestProjErr = projectErr ;
  }
  return seeds[winner];

}

TrajectorySeed MuonSeedCleaner::BetterChi2(std::vector<TrajectorySeed>& seeds ) {

  if ( seeds.size() == 1 ) return seeds[0];

  int winner = 0 ;
  std::vector<int>    moreHits(4,0);  
  std::vector<double> bestChi2(4,99999.);  
  for ( size_t i = 0; i < seeds.size(); i++ ) {

      // 1. fill out the Nchi2 of segments of the seed 
      //GlobalVector mom = SeedMomentum( seeds[i] ); // temporary use for debugging
      //double pt = sqrt( (mom.x()*mom.x()) + (mom.y()*mom.y()) );
      //std::cout<<" > SEED"<<i<<"  pt:"<<pt<< std::endl;

      int it = -1;
      std::vector<double> theChi2(4,99999.);  
      std::vector<int>    theHits(4,0);  
      for (edm::OwnVector<TrackingRecHit>::const_iterator r1 = seeds[i].recHits().first; r1 != seeds[i].recHits().second; r1++){
          it++;
          //std::cout<<"    segmet : "<<it <<std::endl; 
          theHits[it] = NRecHitsFromSegment( *r1 );
          theChi2[it] = NChi2OfSegment( *r1 );
      }
      //std::cout<<" -----  "<<std::endl;

      // 2. longer segment
      for (int j =0; j<4; j++) {

          if ( theHits[j] <  moreHits[j] ) break;

          if ( theHits[j] == moreHits[j] ) { 
            ///  compare the chi2 list
            bool decisionMade = false ;
            for (int k =0; k<4; k++) {
               if ( theChi2[k] >  bestChi2[k] ) break;
               if ( theChi2[k] == bestChi2[k] ) continue;
               if ( theChi2[k] <  bestChi2[k] ) {
                  bestChi2 = theChi2 ;
                  winner = static_cast<int>(i) ;
                  decisionMade = true;
               }
               break;
            }
            if ( decisionMade) break;
            if (!decisionMade) continue;
          }

          if ( theHits[j] >  moreHits[j] ) {
             moreHits = theHits ;
             winner = static_cast<int>(i) ;
          }
          break;
      }
  }
  //std::cout<<" Winner is "<< winner <<std::endl;
  return seeds[winner];
}

SeedContainer MuonSeedCleaner::LengthFilter(std::vector<TrajectorySeed>& seeds ) {
 
  SeedContainer longSeeds; 
  int NSegs = 0;
  for (size_t i = 0; i< seeds.size(); i++) {
    
      int theLength = static_cast<int>( seeds[i].nHits());
      if ( theLength > NSegs ) {
         NSegs = theLength ;
         longSeeds.clear();
         longSeeds.push_back( seeds[i] );
      } 
      else if ( theLength == NSegs ) {
         longSeeds.push_back( seeds[i] );
      } else {
         continue;
      } 
  }
  //std::cout<<" final Length :"<<NSegs<<std::endl;

  return longSeeds ; 
  
}

bool MuonSeedCleaner::MomentumFilter(std::vector<TrajectorySeed>& seeds ) {

  bool findgoodMomentum = false;
  SeedContainer goodMomentumSeeds = seeds;
  seeds.clear();
  for ( size_t i = 0; i < goodMomentumSeeds.size(); i++ ) {
       GlobalVector mom = SeedMomentum( goodMomentumSeeds[i] );
       double pt = sqrt( (mom.x()*mom.x()) + (mom.y()*mom.y()) );
       //if ( pt < 6.  || pt > 3000. ) continue;
       if ( pt < 6. ) continue;
       //std::cout<<" passed momentum :"<< pt <<std::endl;
       seeds.push_back( goodMomentumSeeds[i] );
       findgoodMomentum = true;  
  }
  if ( seeds.size() == 0 ) seeds = goodMomentumSeeds;

  return findgoodMomentum;
}

SeedContainer MuonSeedCleaner::SeedCandidates( std::vector<TrajectorySeed>& seeds, bool good ) {

  SeedContainer theCandidate;
  theCandidate.clear();

  bool longSeed = false;  
  bool withFirstLayer = false ;

  //std::cout<<"***** Seed Classification *****"<< seeds.size() <<std::endl;
  for ( size_t i = 0; i < seeds.size(); i++ ) {

      if (seeds[i].nHits() > 1 ) longSeed = true ;
      //std::cout<<"  Seed: "<<i<< std::endl;
      int idx = 0;
      for (edm::OwnVector<TrackingRecHit>::const_iterator r1 = seeds[i].recHits().first; r1 != seeds[i].recHits().second; r1++){

         idx++;
         const GeomDet* gdet = theService->trackingGeometry()->idToDet( (*r1).geographicalId() );
         DetId geoId = gdet->geographicalId();

         if ( geoId.subdetId() == MuonSubdetId::DT ) {
            DTChamberId DT_Id( (*r1).geographicalId() );
            //std::cout<<" ID:"<<DT_Id <<" pos:"<< r1->localPosition()  <<std::endl;
            if (DT_Id.station() != 1)  continue;
            withFirstLayer = true;
         }
         if ( geoId.subdetId() == MuonSubdetId::CSC ) {
            idx++;
            CSCDetId CSC_Id = CSCDetId( (*r1).geographicalId() );
            //std::cout<<" ID:"<<CSC_Id <<" pos:"<< r1->localPosition()  <<std::endl;
            if (CSC_Id.station() != 1)  continue;
            withFirstLayer = true;
         }
      }
      bool goodseed = (longSeed && withFirstLayer) ? true : false ;
  
      if ( goodseed == good )  theCandidate.push_back( seeds[i] );
  }
  return theCandidate;

}

std::vector<SeedContainer> MuonSeedCleaner::GroupSeeds( std::vector<TrajectorySeed>& seeds) {

  std::vector<SeedContainer> seedCollection;
  seedCollection.clear();
  std::vector<TrajectorySeed> theGroup ;
  std::vector<bool> usedSeed(seeds.size(),false);

  // categorize seeds by comparing overlapping segments or a certian eta-phi cone 
  for (unsigned int i= 0; i<seeds.size(); i++){
    
    if (usedSeed[i]) continue;
    theGroup.push_back( seeds[i] );
    usedSeed[i] = true ;

    GlobalPoint pos1 = SeedPosition( seeds[i]);

    for (unsigned int j= i+1; j<seeds.size(); j++){
 
       // seeds with overlaaping segments will be grouped together
       unsigned int overlapping = OverlapSegments(seeds[i], seeds[j]) ;
       if ( !usedSeed[j] && overlapping > 0 ) {
          // reject the identical seeds
          if ( seeds[i].nHits() == overlapping && seeds[j].nHits() == overlapping ) {
             usedSeed[j] = true ;
             continue;
          }
          theGroup.push_back( seeds[j] ); 
          usedSeed[j] = true ;
       }
       if (usedSeed[j]) continue;

       // seeds in a certain cone are grouped together  
       GlobalPoint pos2 = SeedPosition( seeds[j]);
       double dh = pos1.eta() - pos2.eta() ;
       double df = pos1.phi() - pos2.phi() ;
       double dR = sqrt( (dh*dh) + (df*df) );

       if ( dR > 0.3 && seeds[j].nHits() == 1) continue;
       if ( dR > 0.2 && seeds[j].nHits() >  1) continue;
       theGroup.push_back( seeds[j] ); 
       usedSeed[j] = true ;
    }
    sort(theGroup.begin(), theGroup.end(), lengthSorting ) ;
    seedCollection.push_back(theGroup);
    theGroup.clear(); 
  }
  return seedCollection;

}

unsigned int MuonSeedCleaner::OverlapSegments( TrajectorySeed seed1, TrajectorySeed seed2 ) {

  unsigned int overlapping = 0;
  for (edm::OwnVector<TrackingRecHit>::const_iterator r1 = seed1.recHits().first; r1 != seed1.recHits().second; r1++){
      
      DetId id1 = (*r1).geographicalId();
      const GeomDet* gdet1 = theService->trackingGeometry()->idToDet( id1 );
      GlobalPoint gp1 = gdet1->toGlobal( (*r1).localPosition() );

      for (edm::OwnVector<TrackingRecHit>::const_iterator r2 = seed2.recHits().first; r2 != seed2.recHits().second; r2++){

          DetId id2 = (*r2).geographicalId();
          if (id1 != id2 ) continue;

	  const GeomDet* gdet2 = theService->trackingGeometry()->idToDet( id2 );
	  GlobalPoint gp2 = gdet2->toGlobal( (*r2).localPosition() );

	  double dx = gp1.x() - gp2.x() ;
	  double dy = gp1.y() - gp2.y() ;
	  double dz = gp1.z() - gp2.z() ;
	  double dL = sqrt( dx*dx + dy*dy + dz*dz);

          if ( dL < 1. ) overlapping ++;
      
      }
  }
  return overlapping ;

}

GlobalPoint MuonSeedCleaner::SeedPosition( TrajectorySeed seed ) {

  TrajectoryStateTransform tsTransform;

  PTrajectoryStateOnDet pTSOD = seed.startingState();
  DetId SeedDetId(pTSOD.detId());
  const GeomDet* geoDet = theService->trackingGeometry()->idToDet( SeedDetId );
  TrajectoryStateOnSurface SeedTSOS = tsTransform.transientState(pTSOD, &(geoDet->surface()), &*theService->magneticField());
  GlobalPoint  pos  = SeedTSOS.globalPosition();

  return pos ;

}

GlobalVector MuonSeedCleaner::SeedMomentum( TrajectorySeed seed ) {

  TrajectoryStateTransform tsTransform;

  PTrajectoryStateOnDet pTSOD = seed.startingState();
  DetId SeedDetId(pTSOD.detId());
  const GeomDet* geoDet = theService->trackingGeometry()->idToDet( SeedDetId );
  TrajectoryStateOnSurface SeedTSOS = tsTransform.transientState(pTSOD, &(geoDet->surface()), &*theService->magneticField());
  GlobalVector  mom  = SeedTSOS.globalMomentum();

  return mom ;

}

int MuonSeedCleaner::NRecHitsFromSegment( const TrackingRecHit& rhit ) {

      int NRechits = 0 ; 
      const GeomDet* gdet = theService->trackingGeometry()->idToDet( rhit.geographicalId() );
      MuonTransientTrackingRecHit::MuonRecHitPointer theSeg = 
      MuonTransientTrackingRecHit::specificBuild(gdet, rhit.clone() );

      DetId geoId = gdet->geographicalId();
      if ( geoId.subdetId() == MuonSubdetId::DT ) {
         DTChamberId DT_Id( rhit.geographicalId() );
	 std::vector<TrackingRecHit*> DThits = theSeg->recHits();
	 int dt1DHits = 0;
	 for (size_t j=0; j< DThits.size(); j++) {
             dt1DHits += (DThits[j]->recHits()).size();
         }
         NRechits = dt1DHits ;
      }

      if ( geoId.subdetId() == MuonSubdetId::CSC ) {
         NRechits = (theSeg->recHits()).size() ;
      }
      return NRechits ;
}

int MuonSeedCleaner::NRecHitsFromSegment( MuonTransientTrackingRecHit *rhit ) {

    int NRechits = 0 ; 
    DetId geoId = rhit->geographicalId();
    if ( geoId.subdetId() == MuonSubdetId::DT ) {
       DTChamberId DT_Id( geoId );
       std::vector<TrackingRecHit*> DThits = rhit->recHits();
       int dt1DHits = 0;
       for (size_t j=0; j< DThits.size(); j++) {
           dt1DHits += (DThits[j]->recHits()).size();
       }
       NRechits = dt1DHits ;
       //std::cout<<" D_rh("<< dt1DHits  <<") " ;
    }
    if ( geoId.subdetId() == MuonSubdetId::CSC ) {
       NRechits = (rhit->recHits()).size() ;
       //std::cout<<" C_rh("<<(rhit->recHits()).size() <<") " ;
    }
    return NRechits;

}

double MuonSeedCleaner::NChi2OfSegment( const TrackingRecHit& rhit ) {

      double NChi2 = 999999. ; 
      const GeomDet* gdet = theService->trackingGeometry()->idToDet( rhit.geographicalId() );
      MuonTransientTrackingRecHit::MuonRecHitPointer theSeg = 
      MuonTransientTrackingRecHit::specificBuild(gdet, rhit.clone() );

      double dof = static_cast<double>( theSeg->degreesOfFreedom() );
      NChi2 = theSeg->chi2() / dof ;
      //std::cout<<" Chi2 = "<< NChi2  <<" |" ;

      return NChi2 ;
}


