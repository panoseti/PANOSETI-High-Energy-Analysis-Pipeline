
Int_t           cam_npix;
UShort_t        cam_scope_id;
Double_t        cam_pcap_time;
Double_t        cam_pcap_time_since_start;
UChar_t         cam_acq_mode;
UInt_t          cam_TAI;
UInt_t          cam_nanosec;
Double_t        cam_WR_time;
Short_t        cam_pix_data[32][32];
TH2D * cam_2D_hist;

double ped[32][32];
double pedvar[32][32];
TTree *camdata;

TTree * loadcamdata(const char *infile)
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
