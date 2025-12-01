#ifndef SVTXCALOEVAL_H
#define SVTXCALOEVAL_H

#include <fun4all/SubsysReco.h>

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <string>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase/TrkrDefs.h>
#include <g4eval/SvtxTruthEval.h>
#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/CaloRawClusterEval.h>
#include <calobase/RawCluster.h>

#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Particle.h>

#include <TDatabasePDG.h>

class PHCompositeNode;

struct TrackData {
        double charge;
        double energy;
        Acts::Vector3 momentum;
        Acts::Vector3 vertex;
        int pId;
    };

class SvtxCaloEval : public SubsysReco
{
 public:

  SvtxCaloEval(PHCompositeNode *topNode, const std::string& svtxName, const std::string& caloName);

  ~SvtxCaloEval() override;

  void next_event(PHCompositeNode* topNode);

  bool isSvtxCaloTruth( SvtxTrack* track, RawCluster* cluster);

  bool isSvtxTruth( SvtxTrack* track);

  std::set<int> get_SvtxCalo_trkId();
  
  std::set<int> get_Svtx_trkId();

 private:

  std::string _caloname;
  std::string _svtxname;

  std::set<int> SvtxCalo_trackId_set;
  std::set<int> Svtx_trackId_set;

  std::map<int, TrackData> track_map;

  SvtxEvalStack *svtxeval;
  SvtxTruthEval *svtxtrutheval;
  SvtxClusterEval *svtxcluseval;
  CaloRawClusterEval *caloeval;

  void AddTrackMap(PHG4Particle* p);

  std::set<PHG4Particle*> ReturnTruthParticle(RawCluster* cluster);

  std::set<PHG4Particle*> ReturnTruthParticle(SvtxTrack* track);
};

#endif // SVTXCALOEVAL_H
