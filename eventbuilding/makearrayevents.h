
Int_t           ntel=3;
Int_t           array_nteltrig;
Int_t           array_event_num;
Int_t           array_tel_event[3];
Int_t           array_npix[3];
UShort_t        array_scope_id[3];
Double_t        array_pcap_time[3];
Double_t        array_pcap_time_since_start[3];
UChar_t         array_acq_mode[3];
UInt_t          array_TAI[3];
UInt_t          array_nanosec[3];
Double_t        array_WR_time[3];
Short_t        array_pix_data[3][32][32];
Bool_t         array_triggered[3];
TH2D * array_2D_hist[3];

TTree * define_arraydata()
{
  TTree *arraydata=new TTree("arraydata","Tree with data for a full PANOSETI array");

  arraydata->Branch("ntel", &ntel,"ntel/I");
  arraydata->Branch("array_event_num", &array_event_num,"array_event_num/I");
  arraydata->Branch("array_nteltrig", &array_nteltrig,"array_nteltrig/I");

  arraydata->Branch("array_tel_event", &array_tel_event,"array_tel_event[ntel]/I");
  arraydata->Branch("array_npix", &array_npix,"array_npix[ntel]/I");
  arraydata->Branch("array_scope_id", &array_scope_id,"array_scope_id[ntel]/s");
  arraydata->Branch("array_pcap_time",&array_pcap_time,"array_pcap_time[ntel]/D");
  arraydata->Branch("array_pcap_time_since_start",&array_pcap_time_since_start,"array_pcap_time_since_start[ntel]/D");
  arraydata->Branch("array_acq_mode",&array_acq_mode,"array_acq_mode[ntel]/b");
  arraydata->Branch("array_TAI",&array_TAI,"array_TAI[ntel]/i");
  arraydata->Branch("array_nanosec",&array_nanosec,"array_nanosec[ntel]/i");
  arraydata->Branch("array_WR_time",&array_WR_time,"array_WR_time[ntel]/D");
  arraydata->Branch("array_pix_data",&array_pix_data,"array_pix_data[ntel][32][32]/S");
  arraydata->Branch("array_triggered",&array_triggered,"array_triggered[ntel]/O");
  //arraydata->Branch("array_2D_hist",&array_2D_hist,16000,0);

  

  for (int i=0;i<ntel;i++)
    {
      array_npix[i]=256*4;
      //char arrayname[200];
      //snprintf(arrayname,200,"array_2Dhist%d",i);	
      //array_2D_hist[i]=new TH2D(arrayname,"camera ",32,0,32,32,0,32);
    }
  return arraydata;
}

Int_t           cam_npix=256*4;
UShort_t        cam_scope_id;
Double_t        cam_pcap_time;
Double_t        cam_pcap_time_since_start;
UChar_t         cam_acq_mode;
UInt_t          cam_TAI;
UInt_t          cam_nanosec;
Double_t        cam_WR_time;
Short_t        cam_pix_data[32][32];
TH2D * cam_2D_hist=new TH2D("cam_2D_hist","camera",32,0,32,32,0,32);

TTree * load_camdata(const char *infile)
{
  cout << "reading from " << infile << endl;
  TFile *root_infile = TFile::Open(infile);
  TTree *camdata=(TTree*)root_infile->Get("camdata"); //variable declarations are in the .h file

  camdata->SetBranchAddress("cam_npix", &cam_npix);
  camdata->SetBranchAddress("cam_scope_id", &cam_scope_id);
  camdata->SetBranchAddress("cam_pcap_time",&cam_pcap_time);
  camdata->SetBranchAddress("cam_pcap_time_since_start",&cam_pcap_time_since_start);
  camdata->SetBranchAddress("cam_acq_mode",&cam_acq_mode);
  camdata->SetBranchAddress("cam_TAI",&cam_TAI);
  camdata->SetBranchAddress("cam_nanosec",&cam_nanosec);
  camdata->SetBranchAddress("cam_WR_time",&cam_WR_time);
  camdata->SetBranchAddress("cam_pix_data",cam_pix_data);
  camdata->SetBranchAddress("cam_2D_hist",&cam_2D_hist);
  return camdata;
}
