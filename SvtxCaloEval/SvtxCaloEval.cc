#include "SvtxCaloEval.h"

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

#include <g4main/PHG4Particle.h>

//____________________________________________________________________________..
SvtxCaloEval::SvtxCaloEval(PHCompositeNode *topNode, const std::string& svtxName, const std::string& caloName):
  _caloname(caloName),
  _svtxname(svtxName)
{
  std::cout << "SvtxCaloEval::SvtxCaloEval(const std::string &name) Calling ctor" << std::endl;
  svtxeval = new SvtxEvalStack(topNode);
  svtxeval->set_track_nodename(svtxName);
  caloeval = new CaloRawClusterEval(topNode, caloName);
  svtxeval->set_strict(false);
  caloeval->set_strict(false);
  svtxcluseval = svtxeval->get_cluster_eval();
  svtxtrutheval = svtxeval->get_truth_eval();
}
//____________________________________________________________________________..
SvtxCaloEval::~SvtxCaloEval()
{
  std::cout << "SvtxCaloEval::~SvtxCaloEval() Calling dtor" << std::endl;
}
//____________________________________________________________________________..
void SvtxCaloEval::next_event(PHCompositeNode* topNode)
{
  caloeval->next_event(topNode);
  svtxeval->next_event(topNode);
}
//____________________________________________________________________________..
bool SvtxCaloEval::isSvtxCaloTruth( SvtxTrack* track,  RawCluster* cluster)
{
  // 1. Track → truth particles
  std::set<PHG4Particle*> truth_particles = ReturnTruthParticle(track);
  if (truth_particles.empty()){
    return false;
  }

  std::set<RawCluster*> truth_clusters;

  for (auto* p : truth_particles){
    std::set<RawCluster*> clus_from_p = caloeval->all_clusters_from(p);

    truth_clusters.insert(clus_from_p.begin(), clus_from_p.end());
  }
  if (truth_clusters.empty()){
    return false;
  }
  if (truth_clusters.find(cluster) != truth_clusters.end()){
    for (auto* p : truth_particles){
      int id = p->get_track_id();
      SvtxCalo_trackId_set.insert(id);
      AddTrackMap(p);
    }
    return true;
  }
  return false;
}

//____________________________________________________________________________..
bool SvtxCaloEval::isSvtxTruth( SvtxTrack* track)
{
  std::set<PHG4Particle*> ps = ReturnTruthParticle(track);
  if(!ps.empty()){
    return true;
  }
  for (auto* particle : ps) {
    if (particle) {
      int id = particle->get_track_id(); 
      Svtx_trackId_set.insert(id);
      AddTrackMap(particle);
    }
  }
  
  return false;
}
//____________________________________________________________________________..
std::set<PHG4Particle*> SvtxCaloEval::ReturnTruthParticle(SvtxTrack* track)
{
  std::map<PHG4Particle*, int> particle_count;
  std::vector<TrkrDefs::cluskey> out;
  for (const auto &seed : {track->get_silicon_seed(), track->get_tpc_seed()}){
    if (seed){
      std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
    }
  }
  std::vector<TrkrDefs::cluskey> ckeys = out;
  for (const auto &cluskey : ckeys) {
    std::set<PHG4Particle*> truth_particles = svtxcluseval->all_truth_particles(cluskey);
    for (auto* particle : truth_particles) {
      if (particle) {
	particle_count[particle] += 1;
      }
    }
  }
  if (particle_count.empty()) {
    return {}; 
  }

  auto max_iter = std::max_element(
				   particle_count.begin(),
				   particle_count.end(),
				   [](auto& a, auto& b) { return a.second < b.second; });

  int max_count = max_iter->second;
  std::set<PHG4Particle*> dominant_particles;
  for (auto& kv : particle_count) {
    if (kv.second == max_count) {
      dominant_particles.insert(kv.first);
    }
  }
  int nClusters = track->size_cluster_keys();
  if (max_count != nClusters) {
    return {};
  }
  return dominant_particles;
}
//____________________________________________________________________________..
void SvtxCaloEval::AddTrackMap(PHG4Particle* p)
{
  if (!p) return;  
  int trkId = p->get_track_id();
  auto it = track_map.find(trkId);
  if (it != track_map.end()) {
    return;
  }

  int pdg = p->get_pid();
  auto pdginfo = TDatabasePDG::Instance()->GetParticle(pdg);
  double charge = (pdginfo ? pdginfo->Charge()/3.0 : 0.0);

  TrackData data;
  data.charge   = charge;
  data.energy   = p->get_e();      
  data.pId      = p->get_pid();


  PHG4VtxPoint* vtx = svtxtrutheval->get_vertex(p);

  data.momentum = Acts::Vector3(p->get_px(), p->get_py(), p->get_pz());
  data.vertex   = Acts::Vector3(vtx->get_x(), vtx->get_y(), vtx->get_z());

  track_map.insert({trkId, data});
}
//____________________________________________________________________________..
std::set<int> SvtxCaloEval::get_SvtxCalo_trkId()
{
  return SvtxCalo_trackId_set;
}
//____________________________________________________________________________..
std::set<int> SvtxCaloEval::get_Svtx_trkId()
{
  return Svtx_trackId_set;
}
