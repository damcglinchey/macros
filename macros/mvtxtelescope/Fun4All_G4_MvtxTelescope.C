int Fun4All_G4_MvtxTelescope(
  const int nEvents = 1,
  const char * inputFile = "/sphenix/data/data02/review_2017-08-02/single_particle/spacal2d/fieldmap/G4Hits_sPHENIX_e-_eta0_8GeV-0002.root",
  const char * outputFile = "G4MvtxTelescope.root"
)
{
  //===============
  // Input options
  //===============
  int verbosity = 2;

  int n_maps_layer = 4;
  double maps_layer_radius[4] = {23.4, 31.5, 39.3, 47.1};

  // type 1 = inner barrel stave, 2 = middle barrel stave, 3 = outer barrel stave
  // we use only type 0 here
  int stave_type[4] = {0, 0, 0, 0};
  int staves_in_layer[4] = {1, 1, 1, 1};   // Number of staves per layer in sPHENIX MVTX
  // double phi_tilt[3] = {0.304, 0.304, 0.304};  // radians, from the gdml file, 0.304 radians is 17.4 degrees
  double phi_tilt[4] = {0.0, 0.0, 0.0, 0.0};  // radians, from the gdml file, 0.304 radians is 17.4 degrees

  bool maps_overlapcheck = false; // set to true if you want to check for overlaps

  // Either:
  // read previously generated g4-hits files, in this case it opens a DST and skips
  // the simulations step completely. The G4Setup macro is only loaded to get information
  // about the number of layers used for the cell reco code
  //
  // read files in HepMC format (typically output from event generators like hijing or pythia)
  const bool readhepmc = false; // read HepMC files
  // Or:
  // Use pythia
  const bool runpythia8 = false;

  // Besides the above flags. One can further choose to further put in following particles in Geant4 simulation
  // Use multi-particle generator (PHG4SimpleEventGenerator), see the code block below to choose particle species and kinematics
  const bool particles = false;
  // or gun/ very simple single particle gun generator
  const bool usegun = true;


  //======================
  // What to run
  //======================

  bool do_svtx = true;
  bool do_svtx_cell = do_svtx && true;
  bool do_svtx_track = do_svtx_cell && true;
  bool do_svtx_eval = do_svtx_track && true;

  bool do_dst_compress = false;

  //Option to convert DST to human command readable TTree for quick poke around the outputs
  bool do_DSTReader = false;


  //---------------
  // Load libraries
  //---------------

  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libphhepmc.so");
  gSystem->Load("libg4testbench.so");
  gSystem->Load("libg4hough.so");
  gSystem->Load("libcemc.so");
  gSystem->Load("libg4eval.so");

  int absorberactive = 1; // set to 1 to make all absorbers active volumes
  //  const string magfield = "1.5"; // if like float -> solenoidal field in T, if string use as fieldmap name (including path)
  const string magfield = "/phenix/upgrades/decadal/fieldmaps/sPHENIX.2d.root"; // if like float -> solenoidal field in T, if string use as fieldmap name (including path)
  const float magfield_rescale = 1.4 / 1.5; // scale the map to a 1.4 T field

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);
  // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();
  // By default every random number generator uses
  // PHRandomSeed() which reads /dev/urandom to get its seed
  // if the RANDOMSEED flag is set its value is taken as seed
  // You ca neither set this to a random value using PHRandomSeed()
  // which will make all seeds identical (not sure what the point of
  // this would be:
  //  rc->set_IntFlag("RANDOMSEED",PHRandomSeed());
  // or set it to a fixed value so you can debug your code
  //  rc->set_IntFlag("RANDOMSEED", 12345);

  //-----------------
  // Event generation
  //-----------------

  if (readhepmc)
  {
    // this module is needed to read the HepMC records into our G4 sims
    // but only if you read HepMC input files
    HepMCNodeReader *hr = new HepMCNodeReader();
    se->registerSubsystem(hr);
  }
  else if (runpythia8)
  {
    gSystem->Load("libPHPythia8.so");

    PHPythia8* pythia8 = new PHPythia8();
    // see coresoftware/generators/PHPythia8 for example config
    pythia8->set_config_file("phpythia8.cfg");
    se->registerSubsystem(pythia8);

    HepMCNodeReader *hr = new HepMCNodeReader();
    se->registerSubsystem(hr);
  }

  // If "readhepMC" is also set, the particles will be embedded in Hijing events
  if (particles)
  {
    // toss low multiplicity dummy events
    PHG4SimpleEventGenerator *gen = new PHG4SimpleEventGenerator();
    gen->add_particles("pi-", 2); // mu+,e+,proton,pi+,Upsilon
    //gen->add_particles("pi+",100); // 100 pion option
    if (readhepmc || do_embedding)
    {
      gen->set_reuse_existing_vertex(true);
      gen->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
    }
    else
    {
      gen->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                            PHG4SimpleEventGenerator::Uniform,
                                            PHG4SimpleEventGenerator::Uniform);
      gen->set_vertex_distribution_mean(0.0, 0.0, 0.0);
      gen->set_vertex_distribution_width(0.0, 0.0, 5.0);
    }
    gen->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
    gen->set_vertex_size_parameters(0.0, 0.0);
    gen->set_eta_range(-1.0, 1.0);
    gen->set_phi_range(-1.0 * TMath::Pi(), 1.0 * TMath::Pi());
    //gen->set_pt_range(0.1, 50.0);
    gen->set_pt_range(0.1, 20.0);
    gen->Embed(1);
    gen->Verbosity(0);

    se->registerSubsystem(gen);
  }
  if (usegun)
  {
    PHG4ParticleGun *gun = new PHG4ParticleGun();
    //  gun->set_name("anti_proton");
    gun->set_name("geantino");
    gun->set_vtx(0, 0, 0);
    gun->set_mom(10, 0, 0.01);
    // gun->AddParticle("geantino",1.7776,-0.4335,0.);
    // gun->AddParticle("geantino",1.7709,-0.4598,0.);
    // gun->AddParticle("geantino",2.5621,0.60964,0.);
    // gun->AddParticle("geantino",1.8121,0.253,0.);
    //    se->registerSubsystem(gun);
    PHG4ParticleGenerator *pgen = new PHG4ParticleGenerator();
    pgen->set_name("geantino");
    pgen->set_z_range(0, 0);
    pgen->set_eta_range(0.01, 0.01);
    pgen->set_mom_range(10, 10);
    pgen->set_phi_range(5.3 / 180.*TMath::Pi(), 5.7 / 180.*TMath::Pi());
    se->registerSubsystem(pgen);
  }

  //---------------------
  // Detector description
  //---------------------

  PHG4Reco* g4Reco = new PHG4Reco();
  g4Reco->set_rapidity_coverage(1.1); // according to drawings
// uncomment to set QGSP_BERT_HP physics list for productions
// (default is QGSP_BERT for speed)
  //  g4Reco->SetPhysicsList("QGSP_BERT_HP");

  g4Reco->set_field(0); // use const soleniodal field


  // MVTX Ladders
  for (int ilayer = 0; ilayer < n_maps_layer; ilayer++)
  {
    if (verbosity)
      cout << "Create Maps layer " << ilayer  << " with radius " << maps_layer_radius[ilayer] << " mm, stave type " << stave_type[ilayer]
           << " pixel size 30 x 30 microns " << " active pixel thickness 0.0018 microns" << endl;

    PHG4MapsTelescopeSubsystem  *lyr = new PHG4MapsTelescopeSubsystem("MAPS", ilayer, stave_type[ilayer]);
    lyr->Verbosity(verbosity);

    lyr->set_double_param("layer_nominal_radius", maps_layer_radius[ilayer]); // thickness in cm
    lyr->set_int_param("N_staves", staves_in_layer[ilayer]);    // uses fixed number of staves regardless of radius, if set. Otherwise, calculates optimum number of staves

    // The cell size is used only during pixilization of sensor hits, but it is convemient to set it now because the geometry object needs it
    lyr->set_double_param("pixel_x", 0.0030); // pitch in cm
    lyr->set_double_param("pixel_z", 0.0030); // length in cm
    lyr->set_double_param("pixel_thickness", 0.0018); // thickness in cm
    lyr->set_double_param("phitilt", phi_tilt[ilayer]);

    lyr->set_int_param("active", 1);
    lyr->OverlapCheck(maps_overlapcheck);

    lyr->set_string_param("stave_geometry_file",
                          string(getenv("CALIBRATIONROOT")) + string("/Tracking/geometry/ALICE_ITS_tgeo.gdml"));

    g4Reco->registerSubsystem( lyr );

  }

  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);
  se->registerSubsystem( g4Reco );

  //------------------
  // Detector Division
  //------------------

  if (do_svtx_cell)
  {
    // MAPS cells
    PHG4MapsCellReco *maps_cells = new PHG4MapsCellReco("MAPS");
    maps_cells->Verbosity(verbosity);
    for (int ilayer = 0; ilayer < n_maps_layer; ilayer++)
    {
      maps_cells->set_timing_window(ilayer, -2000, 2000);
    }
    se->registerSubsystem(maps_cells);

  }

  //--------------
  // SVTX tracking
  //--------------

  if (do_svtx_track)
  {
    //----------------------------------
    // Digitize the cell energy into ADC
    //----------------------------------
    PHG4SvtxDigitizer* digi = new PHG4SvtxDigitizer();
    digi->Verbosity(0);
    for (int i = 0; i < n_maps_layer; ++i)
    {
      digi->set_adc_scale(i, 255, 0.4e-6);  // reduced by a factor of 2.5 when going from maps thickess of 50 microns to 18 microns
    }

    se->registerSubsystem( digi );

    //-------------------------------------
    // Apply Live Area Inefficiency to Hits
    //-------------------------------------
    // defaults to 1.0 (fully active)

    PHG4SvtxDeadArea* deadarea = new PHG4SvtxDeadArea();

    for (int i = 0; i < n_maps_layer; i++)
    {
      deadarea->Verbosity(verbosity);
      //deadarea->set_hit_efficiency(i,0.99);
      deadarea->set_hit_efficiency(i, 1.0);
    }
    se->registerSubsystem( deadarea );

    //-----------------------------
    // Apply MIP thresholds to Hits
    //-----------------------------

    PHG4SvtxThresholds* thresholds = new PHG4SvtxThresholds();
    thresholds->Verbosity(verbosity);

    // maps
    for (int i = 0; i < n_maps_layer; i++)
    {
      // reduced by x2.5 when going from cylinder maps with 50 microns thickness to actual maps with 18 microns thickness
      // Note the non-use of set_using_thickness here, this is so that the shortest dimension of the cell sets the mip energy loss
      thresholds->set_threshold(i, 0.1);
    }

    se->registerSubsystem( thresholds );

    //-------------
    // Cluster Hits
    //-------------

    PHG4SvtxClusterizer* clusterizer = new PHG4SvtxClusterizer("PHG4SvtxClusterizer", 0, n_maps_layer - 1);
    clusterizer->Verbosity(verbosity);
    // Reduced by 2 relative to the cylinder cell maps macro. I found this necessary to get full efficiency
    // Many hits in the present simulation are single cell hits, so it is not clear why the cluster threshold should be higher than the cell threshold
    clusterizer->set_threshold(0.1);   // fraction of a mip

    se->registerSubsystem( clusterizer );


  }

  //----------------------
  // Simulation evaluation
  //----------------------

  // if (do_svtx_eval) Svtx_Eval(string(outputFile) + "_g4svtx_eval.root");


  //----------------------
  // Visualization
  //----------------------
  gROOT->LoadMacro("DisplayOn.C");
  PHG4Reco *g4 = DisplayOn();
  //g4->ApplyCommand("/vis/scene/add/axes 0 0 0 50cm");
  //g4->ApplyCommand("/vis/viewer/viewpointThetaPhi 0 0");


  //--------------
  // IO management
  //--------------

  if (readhepmc)
  {
    Fun4AllInputManager *in = new Fun4AllHepMCInputManager( "DSTIN");
    se->registerInputManager( in );
    se->fileopen( in->Name().c_str(), inputFile );
  }
  else
  {
    // for single particle generators we just need something which drives
    // the event loop, the Dummy Input Mgr does just that
    Fun4AllInputManager *in = new Fun4AllDummyInputManager( "JADE");
    se->registerInputManager( in );
  }

  // if (do_DSTReader)
  // {
  //   //Convert DST to human command readable TTree for quick poke around the outputs
  //   gROOT->LoadMacro("G4_DSTReader.C");

  //   G4DSTreader( outputFile, //
  //                /*int*/ absorberactive ,
  //                /*bool*/ do_svtx ,
  //                /*bool*/ do_preshower ,
  //                bool do_cemc ,
  //                /*bool*/ do_hcalin ,
  //                /*bool*/ do_magnet ,
  //                /*bool*/ do_hcalout ,
  //                /*bool*/ do_cemc_twr ,
  //                /*bool*/ do_hcalin_twr ,
  //                /*bool*/ do_magnet  ,
  //                /*bool*/ do_hcalout_twr
  //              );
  // }

  Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", outputFile);
  // if (do_dst_compress) DstCompress(out);
  se->registerOutputManager(out);

  //-----------------
  // Event processing
  //-----------------
  if (nEvents < 0)
  {
    return;
  }
  // if we run the particle generator and use 0 it'll run forever
  if (nEvents == 0 && !readhepmc)
  {
    cout << "using 0 for number of events is a bad idea when using particle generators" << endl;
    cout << "it will run forever, so I just return without running anything" << endl;
    return;
  }

  se->run(nEvents);

  //-----
  // Exit
  //-----

  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
}
