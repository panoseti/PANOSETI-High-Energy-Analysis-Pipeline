
Int_t           npix=256;
Double_t        pcap_time;
Double_t        pcap_time_since_start;
UChar_t         acq_mode;
UChar_t         packer_ver;
UShort_t        packet_no;
UShort_t        boardloc;
UInt_t          TAI;
UInt_t          nanosec;
Double_t        WR_time;
UShort_t        dummy;
UShort_t        pix_data_unsigned[256];  
Short_t        pix_data_signed[256];   
Short_t        pix_data[256];


TTree * loadpdata(const char *infile)
{
   cout << "reading from " << infile << endl;
  TFile *root_infile = TFile::Open(infile);
  TTree *pdata=(TTree*)root_infile->Get("pdata"); //variable declarations are in the .h file

  pdata->SetBranchAddress("npix", &npix);
  pdata->SetBranchAddress("pcap_time", &pcap_time);
  pdata->SetBranchAddress("pcap_time_since_start", &pcap_time_since_start);
  pdata->SetBranchAddress("acq_mode", &acq_mode);
  pdata->SetBranchAddress("packer_ver", &packer_ver);
  pdata->SetBranchAddress("packet_no", &packet_no);
  pdata->SetBranchAddress("boardloc", &boardloc);
  pdata->SetBranchAddress("TAI", &TAI);
  pdata->SetBranchAddress("nanosec", &nanosec);
  pdata->SetBranchAddress("WR_time", &WR_time);
  pdata->SetBranchAddress("dummy", &dummy);
  pdata->SetBranchAddress("pix_data_unsigned", pix_data_unsigned);
  pdata->SetBranchAddress("pix_data_signed", pix_data_signed);

  return pdata;
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

TTree * define_camdata()
{
  TTree *camdata=new TTree("camdata","Tree with data for a full PANOSETI camera");

  camdata->Branch("cam_npix", &cam_npix,"cam_npix/I");
  camdata->Branch("cam_scope_id", &cam_scope_id,"cam_scope_id/s");
  camdata->Branch("cam_pcap_time",&cam_pcap_time,"cam_pcap_time/D");
  camdata->Branch("cam_pcap_time_since_start",&cam_pcap_time_since_start,"cam_pcap_time_since_start/D");
  camdata->Branch("cam_acq_mode",&cam_acq_mode,"cam_acq_mode/b");
  camdata->Branch("cam_TAI",&cam_TAI,"cam_TAI/i");
  camdata->Branch("cam_nanosec",&cam_nanosec,"cam_nanosec/i");
  camdata->Branch("cam_WR_time",&cam_WR_time,"cam_WR_time/D");
  camdata->Branch("cam_pix_data",&cam_pix_data,"cam_pix_data[32][32]/S");
  camdata->Branch("cam_2D_hist",&cam_2D_hist,16000,0);
  return camdata;
}

UShort_t getCameraPosition(UShort_t i_boardloc)
{
  // The two least significant bits of the boardloc variable correspond to the QUABO board location in the camera 
  // So I do some bit-wise logic to find them (because I am smart, and have access to stackoverflow)

  UShort_t mask=0b11; // a mask to AND with the least significant two bits
  UShort_t cam_pos; 
  
  // loop over the four possible positions
  for (UShort_t i=0;i<4;i++) { 
    if (((i_boardloc^i) & mask) == 0) cam_pos=i;
  }
  return cam_pos;
}

UShort_t getScopeID(UShort_t i_boardloc)
{
  UShort_t scope_id=i_boardloc-getCameraPosition(i_boardloc);
  return scope_id;
}
#include "TH2D.h"


TH2D* transpose(TH2D *hin)
{
  // 1. Create a clone with swapped X and Y bin specifications
  // Cloning the original directly creates the wrong axes layout for non-square histograms
  TString name = TString::Format("%s_transposed", hin->GetName());
  TH2D *hout = new TH2D(name, hin->GetTitle(), 
                        hin->GetNbinsY(), hin->GetYaxis()->GetXmin(), hin->GetYaxis()->GetXmax(),
                        hin->GetNbinsX(), hin->GetXaxis()->GetXmin(), hin->GetXaxis()->GetXmax());
  
  hout->SetDirectory(nullptr); // <-- FIX: Disables automatic ROOT tracking and removes warnings
  int nBinsX = hin->GetNbinsX();
  int nBinsY = hin->GetNbinsY();
  
  // 2. Loop using global bin indices rather than 2D coordinate calculations
  for (int ix = 0; ix <= nBinsX + 1; ix++)
    {
      for (int iy = 0; iy <= nBinsY + 1; iy++)
        {
          int inIdx = hin->GetBin(ix, iy);
          int outIdx = hout->GetBin(iy, ix);
          
          hout->SetBinContent(outIdx, hin->GetBinContent(inIdx));
          hout->SetBinError(outIdx, hin->GetBinError(inIdx)); // Keep errors intact
        }
    }
  return hout;
}

TH2D* flipX(TH2D *hin)
{
  TH2D *hout = (TH2D*)hin->Clone("hout_flipX");
  
  int nBinsX = hin->GetNbinsX();
  int nBinsY = hin->GetNbinsY();
  
  // Grab direct pointers to memory arrays to bypass Get/Set method overhead completely
  double* inArray  = hin->GetArray();
  double* outArray = hout->GetArray();
  double* inErrArray  = hin->GetSumw2()->GetArray();
  double* outErrArray = hout->GetSumw2()->GetArray();
  bool hasErrors = (hin->GetSumw2N() > 0);

  for (int ix = 1; ix <= nBinsX; ix++)
    {
      int srcX = nBinsX + 1 - ix;
      for (int iy = 1; iy <= nBinsY; iy++)
        {
          int destIdx = hin->GetBin(ix, iy);
          int srcIdx  = hin->GetBin(srcX, iy);
          
          outArray[destIdx] = inArray[srcIdx];
          if (hasErrors) outErrArray[destIdx] = inErrArray[srcIdx];
        }
    }
  return hout;
}

TH2D* flipY(TH2D *hin)
{
  TH2D *hout = (TH2D*)hin->Clone("hout_flipY");
  
  int nBinsX = hin->GetNbinsX();
  int nBinsY = hin->GetNbinsY();
  
  double* inArray  = hin->GetArray();
  double* outArray = hout->GetArray();
  double* inErrArray  = hin->GetSumw2()->GetArray();
  double* outErrArray = hout->GetSumw2()->GetArray();
  bool hasErrors = (hin->GetSumw2N() > 0);

  for (int ix = 1; ix <= nBinsX; ix++)
    {
      for (int iy = 1; iy <= nBinsY; iy++)
        {
          int srcY = nBinsY + 1 - iy; // FIXED: Changed from GetNbinsX() to GetNbinsY()
          int destIdx = hin->GetBin(ix, iy);
          int srcIdx  = hin->GetBin(ix, srcY);
          
          outArray[destIdx] = inArray[srcIdx];
          if (hasErrors) outErrArray[destIdx] = inErrArray[srcIdx];
        }
    }
  return hout;
}


TH2D* old_transpose(TH2D *hin)
{
  TH2D *hout = new TH2D(*hin);
  hout->SetName("hout");
  hout->Reset();
  for (int ix=1;ix < hin->GetNbinsX()+1;ix++)
    {
      for (int iy=1;iy < hin->GetNbinsY()+1;iy++)
	{
	  hout->SetBinContent(ix,iy,hin->GetBinContent(iy,ix));
	}
    }
  return hout;
}

/*
TH2D* flipX(TH2D *hin)
{
  TH2D *hout = new TH2D(*hin);
  hout->Reset();
  hout->SetName("hout");
  for (int ix=1;ix < hin->GetNbinsX()+1;ix++)
    {
      for (int iy=1;iy < hin->GetNbinsY()+1;iy++)
	{
	  hout->SetBinContent(ix,iy,hin->GetBinContent(hin->GetNbinsX()+1-ix,iy));
	}
    }
  return hout;
}

TH2D* flipY(TH2D *hin)
{
  TH2D *hout = new TH2D(*hin);
  hout->Reset();
  hout->SetName("hout");
  for (int ix=1;ix < hin->GetNbinsX()+1;ix++)
    {
      for (int iy=1;iy < hin->GetNbinsY()+1;iy++)
	{
	  hout->SetBinContent(ix,iy,hin->GetBinContent(ix,hin->GetNbinsX()+1-iy));
	}
    }
  return hout;
}
*/

#include <vector>
// Convert a clean 1D vector back to a TH2D (Row-Major: y * nBinsX + x)
TH2D* convertVectorToTH2(const std::vector<double>& vin, 
                         int nBinsX, double xmin, double xmax,
                         int nBinsY, double ymin, double ymax,
                         const std::string& name = "h2_out", 
                         const std::string& title = "Transformed Histogram")
{
    // Safety check: Ensure the vector matches the specified matrix dimensions
    if (vin.size() != static_cast<size_t>(nBinsX * nBinsY)) {
        return nullptr; 
    }

    // 1. Create the new empty histogram
    TH2D* h2 = new TH2D(name.c_str(), title.c_str(), 
                        nBinsX, xmin, xmax, 
                        nBinsY, ymin, ymax);

    // 2. Populate the bins (skipping ROOT's underflow/overflow bins)
    for (int iy = 0; iy < nBinsY; ++iy) {
        for (int ix = 0; ix < nBinsX; ++ix) {
            // Map 2D loops to flat vector row-major index
            int vectorIdx = iy * nBinsX + ix; 
            
            // Map to ROOT's 1-indexed internal structure
            int rootBinX = ix + 1;
            int rootBinY = iy + 1;

            h2->SetBinContent(rootBinX, rootBinY, vin[vectorIdx]);
        }
    }

    return h2;
}

std::vector<double> convertTH2ToVector(TH2D* h2) 
{
    if (!h2) return std::vector<double>();

    int nBinsX = h2->GetNbinsX();
    int nBinsY = h2->GetNbinsY();
    
    // Allocate exact memory space upfront (excluding overflows)
    std::vector<double> vout;
    vout.reserve(nBinsX * nBinsY);

    // Loop through active bins only (skipping index 0 and NBins+1)
    for (int iy = 1; iy <= nBinsY; ++iy) {
        for (int ix = 1; ix <= nBinsX; ++ix) {
            vout.push_back(h2->GetBinContent(ix, iy));
        }
    }
    return vout;
}


// Helper to convert 2D coordinates to a 1D flat index
inline int getIndex(int x, int y, int nBinsX) {
    return y * nBinsX + x;
}

// 1. Transpose Matrix
std::vector<double> vtranspose(const std::vector<double>& vin, int nBinsX, int nBinsY)
{
    // The output dimensions are swapped (nBinsY becomes the new X/columns)
    std::vector<double> vout(vin.size(), 0.0);

    for (int ix = 0; ix < nBinsX; ++ix) {
        for (int iy = 0; iy < nBinsY; ++iy) {
            int inputIdx  = getIndex(ix, iy, nBinsX);
            // Swapped geometry: ix becomes the row, iy becomes the column
            int outputIdx = getIndex(iy, ix, nBinsY); 
            
            vout[outputIdx] = vin[inputIdx];
        }
    }
    return vout;
}

// 2. Flip Horizontally (X-axis)
std::vector<double> vflipX(const std::vector<double>& vin, int nBinsX, int nBinsY)
{
    std::vector<double> vout(vin.size(), 0.0);

    for (int ix = 0; ix < nBinsX; ++ix) {
        int flippedX = nBinsX - 1 - ix; // 0-indexed flip
        for (int iy = 0; iy < nBinsY; ++iy) {
            int inputIdx  = getIndex(flippedX, iy, nBinsX);
            int outputIdx = getIndex(ix, iy, nBinsX);
            
            vout[outputIdx] = vin[inputIdx];
        }
    }
    return vout;
}

// 3. Flip Vertically (Y-axis)
std::vector<double> vflipY(const std::vector<double>& vin, int nBinsX, int nBinsY)
{
    std::vector<double> vout(vin.size(), 0.0);

    for (int ix = 0; ix < nBinsX; ++ix) {
        for (int iy = 0; iy < nBinsY; ++iy) {
            int flippedY = nBinsY - 1 - iy; // FIX: Use nBinsY for vertical flip
            int inputIdx  = getIndex(ix, flippedY, nBinsX);
            int outputIdx = getIndex(ix, iy, nBinsX);
            
            vout[outputIdx] = vin[inputIdx];
        }
    }
    return vout;
}


/*
  // useful code for testing 2D histogram rotations
  TH2D *htest= new TH2D("htest","htest",16,0,16,16,0,16);;
  htest->Reset();
  int count=0;
  for (int ix=1;ix < htest->GetNbinsX()+1;ix++)
    {
      for (int iy=1;iy < htest->GetNbinsY()+1;iy++)
	{
	  count+=1;
	  htest->SetBinContent(ix,iy,count);
	}
    }

 
  htest->SetStats(0);
  TCanvas *mycanv=new TCanvas();
  mycanv->Divide(4,1);

  mycanv->cd(1);
  htest->Draw("TEXT");

  TH2D *htest_transposed=transpose(htest);
  //htest_transposed->Draw("TEXT");

  TH2D *htest_rotated90=flipY(htest_transposed);
  mycanv->cd(2);
  htest_rotated90->Draw("TEXT");

  TH2D *htest_rotated270=flipX(htest_transposed);
  mycanv->cd(3);
  htest_rotated270->Draw("TEXT");

  TH2D *htest_flipX=flipX(htest);
  TH2D *htest_rotated180=flipY(htest_flipX);
  mycanv->cd(4);
  htest_rotated180->Draw("TEXT");
*/

TH2D * fixquabo763(TH2D *hq)
{
  //quabo 763 was using the wrong firmware during the March 2024 observations
  // this means that the data from two of the four SiPM arrays making up that QUABO need to be rotated by 180 degrees.
  
  TH2D *hll=new TH2D("hll","hll",8,0,8,8,0,8);
  TH2D *hur=new TH2D("hur","hur",8,0,8,8,0,8);
  
  for (int i=0;i<8;i++)
    {
      for (int j=0;j<8;j++)
	{
	  hll->SetBinContent(i+1,j+1,hq->GetBinContent(i+1,j+1));
	  hur->SetBinContent(i+1,j+1,hq->GetBinContent(i+1+8,j+1+8));
	}
    }
  
  TH2D *hll_flippedX=flipX(hll);
  hll_flippedX->SetName("hll_flippedX");
  TH2D *hll_rotated180=flipY(hll_flippedX);
  hll_rotated180->SetName("hll_rotated180");
  
  TH2D *hur_flippedX=flipX(hur);
  hur_flippedX->SetName("hur_flippedX");
  TH2D *hur_rotated180=flipY(hur_flippedX);
  hur_rotated180->SetName("hur_rotated180");
  
  for (int i=0;i<8;i++)
    {
      for (int j=0;j<8;j++)
	{
	  hq->SetBinContent(i+1,j+1,hll_rotated180->GetBinContent(i+1,j+1));
	  hq->SetBinContent(i+1+8,j+1+8,hur_rotated180->GetBinContent(i+1,j+1));
	}
    } 
  return hq;
}

TH2D * makequabohist(TTree * pdata, bool rotate, int zmin, int zmax)
{
  char quabo_id[20];
  snprintf(quabo_id, 20, "QUABO %d", boardloc);
  
  // 1. Create the initial histogram and immediately unregister from the global directory.
  // This safely eliminates memory warnings without globally setting AddDirectory(kFALSE).
  TH2D *hq = new TH2D("hq", quabo_id, 16, 0, 16, 16, 0, 16);
  hq->SetDirectory(nullptr);
  hq->SetStats(0);
  
  // 2. Optimization: Use raw pointers for direct array filling instead of SetBinContent loops
  double* hq_array = hq->GetArray();
  
  for (int i = 0; i < 256; i++) {
    int xbin = 16 - (i / 16); 
    int ybin = 1 + (i % 16);
    
    // Choose the data source based on acquisition mode
    if (acq_mode == 1)      pix_data[i] = pix_data_signed[i];
    else if (acq_mode >= 2) pix_data[i] = pix_data_unsigned[i]; // Merged 2 and 3 into single logic check
    
    // Convert 2D bin coordinates directly to internal ROOT global bin index
    int global_bin = hq->GetBin(xbin, ybin);
    hq_array[global_bin] = pix_data[i];
  }
  
  // Create a pointer to keep track of our final result histogram
  TH2D *hq_transformed = nullptr;
  
  if (rotate) {
    UShort_t cam_pos = getCameraPosition(boardloc); 
    
    // 0b00: top right (180 deg) -> flip X then flip Y
    if (cam_pos == 0b00) {
      TH2D *hq_flippedX = flipX(hq);
      hq_flippedX->SetDirectory(nullptr);
      
      hq_transformed = flipY(hq_flippedX);
      hq_transformed->SetDirectory(nullptr);
      
      delete hq_flippedX; // Clean up intermediate allocation immediately
    }
    
    // 0b01: top left (90 deg) -> transpose then flip Y
    else if (cam_pos == 0b01) {
      TH2D *hq_transposed = transpose(hq);
      hq_transposed->SetDirectory(nullptr);
      
      hq_transformed = flipY(hq_transposed);
      hq_transformed->SetDirectory(nullptr);
      
      delete hq_transposed;
    }

    // 0b10: bottom left (no rotation) -> use hq directly
    else if (cam_pos == 0b10) {
      hq_transformed = hq;
      hq = nullptr; // Clear pointer so we don't accidentally delete it later
    }

    // 0b11: bottom right (270 deg) -> transpose then flip X
    else if (cam_pos == 0b11) {
      TH2D *hq_transposed = transpose(hq);
      hq_transposed->SetDirectory(nullptr);
      
      hq_transformed = flipX(hq_transposed);
      hq_transformed->SetDirectory(nullptr);
      
      delete hq_transposed;
      //if (boardloc==763) hq_transformed = fixquabo763(hq_transformed);
    }
    
  } else {
    // No rotation -> use hq directly
    hq_transformed = hq;
    hq = nullptr; 
  }
  
  // 3. Clean up the original hq histogram if it wasn't adopted as the final object
  if (hq != nullptr) {
    delete hq;
  }
  
  // 4. Safely set final name string and return the object without cloning it again!
  hq_transformed->SetName("hq_transformed_return");
  return hq_transformed;
}

TH2D * old_makequabohist(TTree * pdata, bool rotate, int zmin, int zmax)
{
  char quabo_id[20];
  snprintf(quabo_id,20,"QUABO %d",boardloc);
  //TH1::AddDirectory(kFALSE); //this is dangerous. don't do it.
  TH2D *hq=new TH2D("hq",quabo_id,16,0,16,16,0,16);
  
  //first we get the data from the QUABO into a 2D histogram
  for (int i=0;i<256;i++){
    int xbin=16-(int)i/16; // the "16-" flips the x-axis - seems to be required based on Jerome's star images
    int ybin=1+i%16;
    if (acq_mode==1) pix_data[i]=pix_data_signed[i]; //PH mode. Signed integers
    if (acq_mode==2) pix_data[i]=pix_data_unsigned[i]; // Imaging mode. Unsigned integers
    if (acq_mode==3) pix_data[i]=pix_data_unsigned[i]; // Imaging packet in Dual mode. Unsigned integers
    hq->SetBinContent(xbin,ybin,pix_data[i]);
  }
  hq->SetStats(0);

  

  TH2D *hq_transformed=new TH2D(*hq);
  hq_transformed->SetName("hq_transformed");
  hq_transformed->Reset();
  UShort_t cam_pos=getCameraPosition(boardloc); 
  
  if (rotate){
    // Now we need to rotate the QUABO pixel data, based on the QUABO's location in the camera  
    // positions in the camera are:
    // 0b00: top right;    180 degrees clockwise: flip X and flip Y
    // 0b01: top left;     90 degrees clockwise: transpose and flip Y
    // 0b10: bottom left;  no rotation
    // 0b11: bottom right; 270 degrees clockwise: transpose and flip X

      
    TH2D *hq_transposed = transpose(hq);
    hq_transposed->SetName("hq_transposed");
    if (cam_pos==0b00) {
      //cout << "board " << boardloc << " is top right" << endl;
      //cout << "rotating 180 degrees clockwise" << endl;
      TH2D *hq_flippedX=flipX(hq);
      hq_flippedX->SetName("hq_flippedX");
      TH2D *hq_rotated180=flipY(hq_flippedX);
      hq_rotated180->SetName("hq_rotated180");
      hq_transformed=hq_rotated180;
    }
    
    if (cam_pos==0b01) {
      //cout << "board " << boardloc << " is top left" << endl;
      //cout << "rotating 90 degrees clockwise" << endl;
      TH2D *hq_rotated90=flipY(hq_transposed);
      hq_rotated90->SetName("hq_rotated90");
      hq_transformed=hq_rotated90;
    }

    if (cam_pos==0b10) {
      //cout << "board " << boardloc << " is bottom left" << endl;
      //cout << "no rotation" << endl;
      ////hq_transformed=(TH2D*)hq->Clone("hq_transformed");
      hq_transformed=hq;
    }

    if (cam_pos==0b11) {
      //cout << "board " << boardloc << " is bottom right" << endl;
      //cout << "rotating 270 degrees clockwise" << endl;
      TH2D *hq_rotated270=flipX(hq_transposed);
      hq_rotated270->SetName("hq_rotated270");
      hq_transformed=hq_rotated270;
      //if (boardloc==763) hq_transformed = fixquabo763(hq_transformed);
    }
    
  } else {
    ////hq_transformed=(TH2D*)hq->Clone("hq_transformed");    
    hq_transformed=hq;
  }
  
  hq_transformed->SetName("hq_transformed");
  TH2D *hq_transformed_return = (TH2D*)hq_transformed->Clone("hq_transformed_return"); //Hack to allow deleting hq
  delete hq;  //Now delete to avoid a memory leak
  return hq_transformed_return;
  //return hq_transformed;  
 
  
}

void initializeCamera()
{
  for (int x=0;x<32;x++)
    {
      for (int y=0;y<32;y++)
	{
	  cam_pix_data[x][y]=0;
	  int xbin=x+1;
	  int ybin=y+1;
	  cam_2D_hist->SetBinContent(xbin,ybin,0);
	}
    }
  cam_2D_hist->SetStats(0);
  //cam_2D_hist->SetMinimum(-50);
  //cam_2D_hist->SetMaximum(100);
}

void fillCamera(UShort_t i_boardloc, TH2D *hquabo)
{
  //cout << "filling camera. board:" << i_boardloc << endl;
  int xoff=0;
  int yoff=0;
  if (i_boardloc<0 || i_boardloc>3) {
    cout << "invalid quabo board location in fillCamera:" << i_boardloc << endl;
    return;
  }
  
  if (i_boardloc==0b00) {
    xoff=16;
    yoff=16;
  }
  if (i_boardloc==0b01) {
    xoff=0;
    yoff=16;
  }
  if (i_boardloc==0b10) {
    xoff=0;
    yoff=0;
  }
  if (i_boardloc==0b11) {
    xoff=16;
    yoff=0;
  }
  
  int xbin=0;
  int ybin=0;
  if (hquabo->GetNbinsX()!=16) cout << "WHy doesn't the quabo have 16 X bins?" << endl;
    for (int x=0;x<16;x++)
    {
      for (int y=0;y<16;y++)
	{
	  cam_pix_data[x+xoff][y+yoff]=hquabo->GetBinContent(x+1,y+1);
	  int xbin=x+1+xoff;
	  int ybin=y+1+yoff;
	  cam_2D_hist->SetBinContent(xbin,ybin,hquabo->GetBinContent(x+1,y+1));
	}
    }
}
 

//TH2D *cam_2D_hist;
//camdata->SetBranchAddress("cam_2D_hist",&cam_2D_hist);
//camdata->GetEntry(13)
//cam_2D_hist->Draw("COLZ");
