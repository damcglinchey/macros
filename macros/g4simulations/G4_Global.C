
void Global_Reco(int verbosity = 0) {
  
  //---------------
  // Load libraries
  //---------------

  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4vertex.so");

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();

  GlobalVertexReco* gblvertex = new GlobalVertexReco();
  se->registerSubsystem(gblvertex);

  return;
}  

void Global_FastSim(int verbosity = 0) {
  
  //---------------
  // Load libraries
  //---------------

  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4vertex.so");

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();

  GlobalVertexFastSimReco* gblvertex = new GlobalVertexFastSimReco();
  gblvertex->set_x_smearing(0.0100); // 100 um
  gblvertex->set_y_smearing(0.0100); // 100 um
  gblvertex->set_z_smearing(0.0150); // 150 um
  gblvertex->set_t_smearing(0.002);  // 20 ps
  se->registerSubsystem(gblvertex);

  return;
}  
