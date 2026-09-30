#include "./processcamevents.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "TH2F.h"
#include "TRandom3.h"


struct HillasParams {
    double size;
    double cogX, cogY;
    double length, width;
    double psi;   // Orientation of the ellipse major axis
    double alpha; // Angle between major axis and vector to origin
};

HillasParams calculateHillas(double image[32][32], double threshold) {
    HillasParams p = {0, 0, 0, 0, 0, 0, 0};
    double sumQ = 0, sumX = 0, sumY = 0;
    double sumX2 = 0, sumY2 = 0, sumXY = 0;

    // Center of the camera (assuming origin is center of 32x32 array)
    double centerX = 15.5; 
    double centerY = 15.5;

    // 1. First Moments (COG)
    for (int i = 0; i < 32; ++i) {
        for (int j = 0; j < 32; ++j) {
            if (image[i][j] > threshold) {
                sumQ += image[i][j];
                sumX += image[i][j] * i;
                sumY += image[i][j] * j;
            }
        }
    }

    if (sumQ <= 0) return p;

    p.size = sumQ;
    p.cogX = sumX / sumQ;
    p.cogY = sumY / sumQ;

    // 2. Second Moments
    for (int i = 0; i < 32; ++i) {
        for (int j = 0; j < 32; ++j) {
            if (image[i][j] > threshold) {
                sumX2 += image[i][j] * pow(i - p.cogX, 2);
                sumY2 += image[i][j] * pow(j - p.cogY, 2);
                sumXY += image[i][j] * (i - p.cogX) * (j - p.cogY);
            }
        }
    }

    double varX = sumX2 / sumQ;
    double varY = sumY2 / sumQ;
    double covXY = sumXY / sumQ;

    // Orientation (Psi)
    double d = varX - varY;
    double z = sqrt(d * d + 4 * pow(covXY, 2));
    
    // Psi in radians
    double psi_rad = 0.5 * atan2(2.0 * covXY, d);
    p.psi = psi_rad * (180.0 / M_PI);

    // Shape (Length/Width)
    p.length = sqrt(0.5 * (varX + varY + z));
    p.width  = sqrt(0.5 * (varX + varY - z));

    // 3. Alpha Calculation
    // Vector from camera center to COG
    double dx = p.cogX - centerX;
    double dy = p.cogY - centerY;
    
    if (dx != 0 || dy != 0) {
        // Angle of the vector to COG
        double phi = atan2(dy, dx);
        // Alpha is the difference between psi and phi
        double alpha_rad = asin(fabs(sin(psi_rad - phi)));
        p.alpha = alpha_rad * (180.0 / M_PI);
    }

    return p;
}

void calcSkygains_old(const char *infile,int ii=0,int jj=0)
{
  camdata=loadcamdata(infile);
  camdata->GetEntry(0);
  cout << "Acquisition mode " << (int)cam_acq_mode << endl;
  TH2D *skygains_2D_hist=new TH2D("h","h",32,0,32,32,0,32);
  skygains_2D_hist->SetName("skygains_2D_hist");
  skygains_2D_hist->SetTitle("Sky Gains");
  
  TH2D *gaincorrected_pedvars=new TH2D("h","h",32,0,32,32,0,32);
  gaincorrected_pedvars->SetName("gaincorrected_pedvars");
  gaincorrected_pedvars->SetTitle("Gain-corrected pedvars");

  double hmin,hmax,hbins;
  hmin=-100;
  hmax=1000;
  hbins=1100;
  TH1D *phd=new TH1D("phd","Pulse Height Distribution",hbins,hmin,hmax);
  TH1D *iphd=new TH1D("iphd","Integral Pulse Height Distribution",hbins,hmin,hmax);

  TCanvas *mycanv=new TCanvas();
  mycanv->SetLogx();
  mycanv->SetLogy();
  mycanv->Divide(1,2);
  mycanv->Draw();
  mycanv->cd(1);
  double sumgain=0;
  int ngain=0;
  for (int i=0;i<32;i++)
    {
      for (int j=0;j<32;j++)
	{
	  //int i=ii;
	  //int j=jj;
	  phd->Reset();
	  iphd->Reset();
	  char variable[200];
	  char cut[200];
	  snprintf(variable,200,"cam_pix_data[%d][%d]>>phd",i,j);
	  snprintf(cut,200,"cam_pix_data[%d][%d]<%f",i,j,hmax);
	  camdata->Draw(variable,cut);
	  int nevents=phd->GetEntries();
	  double int_ph=0;
	  double thresh1=0;
	  double thresh2=0;
	  bool thresh1flag=0;
	  bool thresh2flag=0;
	  for (int k=0;k<phd->GetNbinsX();k++)
	    {
	      int this_bin=phd->GetNbinsX()-k;
	      
	      int_ph+=phd->GetBinContent(this_bin);
	      iphd->SetBinContent(this_bin,int_ph);
	      //cout << this_bin << "\t" << int_ph << endl;
	      
	      if (int_ph>0.03*nevents && !thresh1flag)
		{
		  //cout << thresh1 << "\t" << thresh1flag << "\t" << this_bin << "\t" << int_ph << endl;
		  thresh1=this_bin;
		  thresh1flag=1;
		 
		}
	      if (int_ph>0.3*nevents && !thresh2flag)
		{
		  //cout << thresh2 << "\t" << thresh2flag << "\t" << this_bin << "\t" << int_ph << endl;
		  thresh2=this_bin;
		  thresh2flag=1;
		}
	      
	    }
	  
	  cout << i << "\t" << j << "\t" << (thresh1 - thresh2) << endl;
	  
	  iphd->GetXaxis()->SetRangeUser(10,150);
	  iphd->SetXTitle("pulse height (digital counts)");
	  phd->SetXTitle("pulse height (digital counts)");
	  phd->GetXaxis()->SetRangeUser(10,150);
	  iphd->SetStats(0);
	  iphd->SetMaximum(nevents);
	  iphd->Draw();
	  TLine *l1=new TLine(thresh1+hmin,0,thresh1+hmin,nevents);
	  TLine *l2=new TLine(thresh2+hmin,0,thresh2+hmin,nevents);
	  l2->SetLineColor(2);
	  l1->Draw("SAME");
	  l2->Draw("SAME");
	  skygains_2D_hist->SetBinContent(i+1,j+1,thresh1 - thresh2);
	  sumgain+=thresh1 - thresh2;
	  ngain+=1;
	  mycanv->cd(2);
	  phd->Draw();
	  //std::cin.get();
	}
    }
  double meangain=sumgain/ngain;
  //normalize gains to 1;
  for (int i=0;i<32;i++)
    {
      for (int j=0;j<32;j++)
	{
	  double this_gain=skygains_2D_hist->GetBinContent(i+1,j+1)/meangain;
	  skygains_2D_hist->SetBinContent(i+1,j+1,this_gain);
	  if (this_gain>0)
	    {
	      gain[i][j]=this_gain;
	    }
	  else
	    {
	      gain[i][j]=1;
	    }
	  gaincorrected_pedvars->SetBinContent(i+1,j+1,pedvar[i][j]/gain[i][j]);
	}
    }

  skygains_2D_hist->Draw("COLZ");
}

void calcSkygains(const char *infile)
{
  // Assume loadcamdata returns a pointer containing the TTree
  TTree* camdata = loadcamdata(infile);
  if (!camdata) return;

  // 1. Turn off all branches, then turn on only the pixel data branch to save I/O time
  camdata->SetBranchStatus("*", 0);
  camdata->SetBranchStatus("cam_pix_data", 1);

  // 2. Bind the tree branch directly to a local array buffer
  // Adjust the data type (e.g., float, int, double) to match your tree definition exactly
  float cam_pix_data[32][32]; 
  camdata->SetBranchAddress("cam_pix_data", &cam_pix_data);

  // 3. Instantiate the histograms
  TH2D *skygains_2D_hist = new TH2D("skygains_2D_hist","Sky Gains",32,0,32,32,0,32);
  TH2D *gaincorrected_pedvars = new TH2D("gaincorrected_pedvars","Gain-corrected pedvars",32,0,32,32,0,32);

  double hmin = -100;
  double hmax = 1000;
  int hbins = 1100;

  // Create an array of histograms so we can fill them in a single file-pass
  TH1D* phd[32][32];
  TH1D* iphd[32][32];
  
  for (int i=0; i<32; i++) {
    for (int j=0; j<32; j++) {
      TString hname = Form("phd_%d_%d", i, j);
      phd[i][j] = new TH1D(hname, "Pulse Height Distribution", hbins, hmin, hmax);
      
      TString iphd_name = Form("iphd_%d_%d", i, j);
      iphd[i][j] = new TH1D(iphd_name, "Integral Pulse Height Distribution", hbins, hmin, hmax);
    }
  }

  // 4. CRITICAL: Single loop over data entries
  Long64_t nentries = camdata->GetEntries();
  for (Long64_t entry = 0; entry < nentries; ++entry) {
    camdata->GetEntry(entry); // Reads only the cam_pix_data branch for this entry
    
    for (int i=0; i<32; i++) {
      for (int j=0; j<32; j++) {
        float val = cam_pix_data[i][j];
        if (val < hmax) { // Applies your original cut string rule
          phd[i][j]->Fill(val);
        }
      }
    }
  }

  // 5. Calculate thresholds and gains sequentially in memory
  TCanvas *mycanv = new TCanvas("mycanv", "Canvas", 800, 800);
  mycanv->SetLogx();
  mycanv->SetLogy();
  mycanv->Divide(1,2);

  double sumgain = 0;
  int ngain = 0;

  for (int i=0; i<32; i++) {
    for (int j=0; j<32; j++) {
      int nevents = phd[i][j]->GetEntries();
      if (nevents == 0) continue;

      double int_ph = 0;
      double thresh1 = 0;
      double thresh2 = 0;
      bool thresh1flag = false;
      bool thresh2flag = false;
      int nbinsX = phd[i][j]->GetNbinsX();

      for (int k=0; k<nbinsX; k++) {
        int this_bin = nbinsX - k;
        int_ph += phd[i][j]->GetBinContent(this_bin);
        iphd[i][j]->SetBinContent(this_bin, int_ph);
        
        if (int_ph > 0.005 * nevents && !thresh1flag) {
          thresh1 = this_bin;
          thresh1flag = true;
        }
        if (int_ph > 0.01 * nevents && !thresh2flag) {
          thresh2 = this_bin;
          thresh2flag = true;
        }
      }
      
      double delta_thresh = thresh1 - thresh2;
      skygains_2D_hist->SetBinContent(i+1, j+1, delta_thresh);
      sumgain += delta_thresh;
      ngain += 1;

      // Optional: Draw only the very last pixel's metrics on screen to avoid Canvas freezing
      if (i == 31 && j == 31) {
        mycanv->cd(1);
        iphd[i][j]->GetXaxis()->SetRangeUser(10,150);
        iphd[i][j]->SetXTitle("pulse height (digital counts)");
        iphd[i][j]->SetStats(0);
        iphd[i][j]->SetMaximum(nevents);
        iphd[i][j]->Draw();

        TLine *l1 = new TLine(thresh1+hmin, 0, thresh1+hmin, nevents);
        TLine *l2 = new TLine(thresh2+hmin, 0, thresh2+hmin, nevents);
        l2->SetLineColor(2);
        l1->Draw("SAME");
        l2->Draw("SAME");

        mycanv->cd(2);
        phd[i][j]->GetXaxis()->SetRangeUser(10,150);
        phd[i][j]->SetXTitle("pulse height (digital counts)");
        phd[i][j]->Draw();
      }
    }
  }
  
  // 6. Normalization pass
  double meangain = (ngain > 0) ? (sumgain / ngain) : 1.0;
  for (int i=0; i<32; i++) {
    for (int j=0; j<32; j++) {
      double this_gain = skygains_2D_hist->GetBinContent(i+1, j+1) / meangain;
      skygains_2D_hist->SetBinContent(i+1, j+1, this_gain);
      
      gain[i][j] = (this_gain > 0) ? this_gain : 1.0;
      gaincorrected_pedvars->SetBinContent(i+1, j+1, pedvar[i][j] / gain[i][j]);
    }
  }

  // Draw final 2D summary mapping
  TCanvas *c2 = new TCanvas("c2","Final Gains Mapping", 600, 600);
  skygains_2D_hist->Draw("COLZ");
}


void calcPedestals(const char *infile, bool do_fit = true) {
    camdata = loadcamdata(infile);
    long nEntries = camdata->GetEntries();
    
    // 1. Create a 1D histogram for EVERY pixel once
    // Using a 1D array of pointers to manage 1024 histograms
    TH1D* h_pixels[32][32];
    for(int i=0; i<32; ++i) {
        for(int j=0; j<32; ++j) {
            h_pixels[i][j] = new TH1D(Form("h_%d_%d", i, j), "", 1000, -500, 500);
            h_pixels[i][j]->SetDirectory(0); // Faster, keeps them out of global ROOT list
        }
    }

    // 2. Single loop over the TTree (The "Speed-Up")
    cout << "Processing " << nEntries << " entries for pedestals..." << endl;
    for (long entry = 0; entry < nEntries; ++entry) {
        camdata->GetEntry(entry);
        for (int i = 0; i < 32; ++i) {
            for (int j = 0; j < 32; ++j) {
                h_pixels[i][j]->Fill(cam_pix_data[i][j]);
            }
        }
    }

    // 3. Extract results (Mean/RMS or Fit)
    TH2D *peds_hist = new TH2D("peds_2D_hist", "Pedestals", 32, 0, 32, 32, 0, 32);
    TH2D *vars_hist = new TH2D("pedvars_2D_hist", "Pedestal Variances", 32, 0, 32, 32, 0, 32);

    for (int i = 0; i < 32; ++i) {
        for (int j = 0; j < 32; ++j) {
            double thisPed, thisVar;
            if (do_fit && h_pixels[i][j]->GetEntries() > 10) {
                h_pixels[i][j]->Fit("gaus", "Q0"); // "Q" for quiet, "0" to not draw
                TF1 *g = h_pixels[i][j]->GetFunction("gaus");
                thisPed = g ? g->GetParameter(1) : h_pixels[i][j]->GetMean();
                thisVar = g ? g->GetParameter(2) : h_pixels[i][j]->GetRMS();
            } else {
                thisPed = h_pixels[i][j]->GetMean();
                thisVar = h_pixels[i][j]->GetRMS();
            }
            
            ped[i][j] = thisPed;
            pedvar[i][j] = thisVar;
            peds_hist->SetBinContent(i+1, j+1, thisPed);
            vars_hist->SetBinContent(i+1, j+1, thisVar);
            
            delete h_pixels[i][j]; // Clean up memory immediately
        }
    }

    // 4. Save results
    TFile *fOut = TFile::Open(Form("%s.pedvars", infile), "RECREATE");
    peds_hist->Write();
    vars_hist->Write();

    TCanvas *c1 = new TCanvas("c1", "Sequential Image Display", 800, 400);
    c1->Divide(2,1);
    c1->cd(1);
    peds_hist->Draw("COLZ");
    c1->cd(2);
    vars_hist->Draw("COLZ");

    fOut->Close();
}

void calcPedestals_old(const char *infile, bool do_fit=1)
{
  // this was the initial pedestal calculation routine. It is very slow. 
  camdata=loadcamdata(infile);
  camdata->GetEntry(0);
  cout << "Acquisition mode " << (int)cam_acq_mode << endl;
  TH2D *pedvars_2D_hist=new TH2D("h","h",32,0,32,32,0,32);
  pedvars_2D_hist->SetName("pedvars_2D_hist");
  pedvars_2D_hist->SetTitle("Pedestal Variances");
  
  TH2D *peds_2D_hist=new TH2D("h","h",32,0,32,32,0,32);
  peds_2D_hist->SetName("peds_2D_hist");
  peds_2D_hist->SetTitle("Pedestals");
  
  double hmin,hmax,hbins;
  if ((int)cam_acq_mode>1){
    //imaging data
    hmin=0;
    hmax=20000;
    hbins=20000;
  } else {
    //pulse height data
    hmin=-500;
    hmax=500;
    hbins=1000;
  }
  TH1D *h1=new TH1D("h1","h1",hbins,hmin,hmax);

  for (int i=0;i<32;i++)
    {
      for (int j=0;j<32;j++)
	{
	  h1->Reset();
	  char variable[200];
	  char cut[200];
	  snprintf(variable,200,"cam_pix_data[%d][%d]>>h1",i,j);
	  snprintf(cut,200,"cam_pix_data[%d][%d]<%f",i,j,hmax);
	  camdata->Draw(variable,cut);
	  double firstmean=0;
	  double firstRMS=0;
	  if (h1->GetEffectiveEntries()>10)
	    {
	      firstmean=h1->GetMean();
	      firstRMS=h1->GetRMS();
	    }
	  //cout << "first mean: " << firstmean;
	  //cout << " first RMS: " << firstRMS << endl;
	  if ((int)cam_acq_mode==1)  // only need to iterate for pulse height mode
	    {
	      h1->Reset();
	      snprintf(cut,200,"cam_pix_data[%d][%d]< (%f+5*%f)",i,j,firstmean,firstRMS);
	      camdata->Draw(variable,cut);
	      //cout << " second mean: " << h1->GetMean();
	      //cout << " second RMS: " << h1->GetRMS() << endl;
	    }
	  double thisped;
	  double thispedvar;
	  if (do_fit)
	    {
	      cout << i << " " << j << " "<< h1->GetEntries() << endl;
	      if (h1->GetEffectiveEntries()>10)
		{
		  h1->Fit("gaus");
		  TF1 *g = (TF1*)h1->GetListOfFunctions()->FindObject("gaus");
		  thisped=g->GetParameter(1);
		  thispedvar=g->GetParameter(2);
		}
	    }
	  else
	    {
	      thisped=h1->GetMean();
	      thispedvar=h1->GetRMS();
	    }
	  ped[i][j]= thisped;
	  pedvar[i][j]=thispedvar;
	  peds_2D_hist->SetBinContent(i+1,j+1,thisped);
	  pedvars_2D_hist->SetBinContent(i+1,j+1,thispedvar);
	}
    }
  pedvars_2D_hist->SetStats(0);
  pedvars_2D_hist->Draw("COLZ");
  char outfile[200];
  snprintf(outfile,200,"%s.pedvars",infile);
  TFile *root_outfile = TFile::Open(outfile,"RECREATE");
  pedvars_2D_hist->Write();
  
  peds_2D_hist->SetStats(0);
  peds_2D_hist->Write();
  //peds_2D_hist->Draw("COLZ");
}

void draw_cam(int event)
{
  //TTree *camdata=loadcamdata(infile);
  //TH2D *cam_2D_hist;
  //camdata->SetBranchAddress("cam_2D_hist",&cam_2D_hist);
  camdata->GetEntry(event);
  cam_2D_hist->SetStats(0);
  cam_2D_hist->Draw("COLZ");
}

void read_pedestals(const char *infile)
{
  
  TFile *pedfile = TFile::Open(infile,"READ");
  TH2D *peds_2D_hist=(TH2D*)pedfile->Get("peds_2D_hist");
  TH2D *pedvars_2D_hist=(TH2D*)pedfile->Get("pedvars_2D_hist");
  for (int i=0;i<32;i++)
    {
      for (int j=0;j<32;j++)
	{
	  ped[i][j]=peds_2D_hist->GetBinContent(i+1,j+1);
	  pedvar[i][j]=pedvars_2D_hist->GetBinContent(i+1,j+1);
	}
    }
  

}
void read_gains(const char *infile)
{
  
  TFile *gainfile = TFile::Open(infile,"READ");
  TH2D *relgain_2D_hist=(TH2D*)gainfile->Get("relgain_2D_hist");
  for (int i=0;i<32;i++)
    {
      for (int j=0;j<32;j++)
	{
	  gain[i][j]=relgain_2D_hist->GetBinContent(i+1,j+1)*relgain_2D_hist->GetBinContent(i+1,j+1);
	}
    }
  

}

void prefilter(const char *infile, const char *pedfile, double image_thresh=4.25, double border_thresh=2.25)
{
  //apply a low level cleaning to the camera files to remove obvious noise events. Makes the following stages of the analysis much more manageable when the rates are high.
  gSystem->mkdir("./prefiltered/", kTRUE);
  
  camdata = loadcamdata(infile);
  long nEntries = camdata->GetEntries();

  char outfile_name[200];
  snprintf(outfile_name,200,"prefiltered/%s",infile);
  TFile *newfile = new TFile(outfile_name,"recreate");
  TTree *camdata_prefilter = camdata->CloneTree(0);

  read_pedestals(pedfile);
  for (int i=0; i<nEntries; i++)
    {
      camdata->GetEntry(i);
      int nimage=0;
      int nborder=0;
      int nborderorimage=0;
      for (int i=0;i<32;i++)
	{
	  for (int j=0;j<32;j++)
	    {
	      double SN=(cam_pix_data[i][j]-ped[i][j])/pedvar[i][j];
	      if (SN>image_thresh) nimage+=1;
	      if (SN>border_thresh && SN<image_thresh) nborder+=1;
	      if (SN>border_thresh) nborderorimage+=1;
	    }
	}
      if (nimage>2 && nborderorimage>4) camdata_prefilter->Fill();
    }
  camdata_prefilter->AutoSave();
  cout << "Started with " << nEntries << endl;
  cout << "Filtered down to  " << camdata_prefilter->GetEntries() << endl;
  cout << "Written to: "<< outfile_name << endl;
  newfile->Close();
}

void eventdisplay_fast(int start=0, int end=-1, bool display=false) {
    double image_thresh=4.25;
    double border_thresh=2.25;
    
    if (!display) gROOT->SetBatch(kTRUE);
    
    // 1. Move canvas creation inside the conditional block
    TCanvas *c1 = nullptr;
    if (display) {
        c1 = new TCanvas("c1", "sequential image display", 1200, 400);
        c1->Divide(3,1);
    }
    
    long long max_entries = camdata->GetEntries();
    if (end > max_entries || end < start) end = max_entries;
    
    std::ofstream outfile("output.txt");
    //camdata->GetEntry(1);
    TH2D *cam_ped_subtracted = cam_ped_subtracted=new TH2D(*cam_2D_hist);
    TH2D *cam_gain_corrected = cam_gain_corrected=new TH2D(*cam_2D_hist);

    // 2. Reuse histogram objects out of the loop to avoid heavy allocations
    //TH2D *cam_ped_subtracted = (TH2D*)cam_2d_hist->Clone("cam_ped_subtracted");
    //TH2D *cam_gain_corrected = (TH2D*)cam_2d_hist->Clone("cam_gain_corrected");
    
    cam_ped_subtracted->SetStats(0);
    cam_gain_corrected->SetStats(0);
    
    char mytitle[128]; // Smaller stack allocation
    
    // 3. Main Loop
    for (int i=start; i<end; i++) {
        camdata->GetEntry(i);
        
        // Update titles without recreating objects
        snprintf(mytitle, sizeof(mytitle), "event %i pedestal subtracted camera", i);
        cam_ped_subtracted->SetTitle(mytitle);
        
        snprintf(mytitle, sizeof(mytitle), "event %i pedestal subtracted and gain corrected camera", i);
        cam_gain_corrected->SetTitle(mytitle);
        
        // 4. Cache multi-dimensional array data to continuous pointers if possible, 
        // or ensure fast 1D array indexing into ROOT histograms
        for (int x=0; x<32; x++) {
            for (int y=0; y<32; y++) {
                // Pedestal subtraction
                double ped_sub_val = cam_pix_data[x][y] - ped[x][y];
                cam_ped_subtracted->SetBinContent(x+1, y+1, ped_sub_val);
                
                // Gain correction logic (safeguard division by zero)
                double current_gain = gain[x][y];
                if (current_gain <= 0.0) current_gain = 1.0;
                
                double gain_corr_val = ped_sub_val * current_gain; 
                cam_gain_corrected->SetBinContent(x+1, y+1, gain_corr_val);
                
                // (Your cleaning logic using image_thresh / border_thresh goes here)
            }
        }
        
        if (display && c1) {
            c1->cd(1);
            cam_ped_subtracted->Draw("colz");
            c1->Update(); // Force canvas updates if displaying interactively
        }
    }
    
    // 5. Clean up ROOT objects from memory safely when finished
    if (!display) {
        delete cam_ped_subtracted;
        delete cam_gain_corrected;
    }
}
 
void eventdisplay(int start=0, int end=-1, bool display=false)
{
  //If you calculated pedestals in a previous ROOT session, you have to do this first:
  //read_pedestals("../Fern/rawdata/Mrk421_preflip.root.pedvars")
  //read_gains("../Fern/rawdata/Mrk421_preflip.root.gain")
  //camdata=loadcamdata("../Fern/rawdata/Mrk421_preflip.root")
  
  double image_thresh=4.25;
  double border_thresh=2.25;

  if (!display) gROOT->SetBatch(kTRUE);
  
  TCanvas *c1 = new TCanvas("c1", "Sequential Image Display", 1200, 400);
  
  if (display) c1->Divide(3,1);
    
  if (end>camdata->GetEntries()) end=camdata->GetEntries();
  if (end<start) end=camdata->GetEntries();
  std::ofstream outFile("output.txt");
  TH2D *cam_ped_subtracted;
     for (int i=start;i<end;i++)
    {
      camdata->GetEntry(i);
      cam_ped_subtracted=new TH2D(*cam_2D_hist);
      cam_ped_subtracted->SetName("cam_ped_subtracted");
      char mytitle[200];
      snprintf(mytitle,200,"Event %i Pedestal Subtracted Camera",i);
      cam_ped_subtracted->SetTitle(mytitle);
      for (int i=0;i<32;i++)
	{
	  for (int j=0;j<32;j++)
	    {
	      cam_ped_subtracted->SetBinContent(i+1,j+1,(cam_pix_data[i][j]-ped[i][j]));
	    }
	}
      cam_ped_subtracted->SetStats(0);
      if (display) c1->cd(1);
      //cam_ped_subtracted->SetMinimum(1300);
      //cam_ped_subtracted->SetMaximum(2000);
      if (display) cam_ped_subtracted->Draw("COLZ");

      TH2D *cam_gain_corrected=new TH2D(*cam_2D_hist);
      cam_gain_corrected->SetName("cam_gain_corrected");
      snprintf(mytitle,200,"Event %i Pedestal Subtracted and Gain Corrected Camera",i);
      cam_gain_corrected->SetTitle(mytitle);
      double pix_value[32][32];
      double calibrated_pix_value[32][32];
      double cleaned_pix_value[32][32];
      //double* cleaned_data = cam_cleaned->GetArray();
      //double* gain_data = cam_gain_corrected->GetArray();

      for (int i=0;i<32;i++)
	{
	  for (int j=0;j<32;j++)
	    {
	      if (gain[i][j]<=0) gain[i][j]=1.0;
	      calibrated_pix_value[i][j]=(cam_pix_data[i][j]-ped[i][j])/(gain[i][j]);
	      cam_gain_corrected->SetBinContent(i+1,j+1,calibrated_pix_value[i][j]);
	    }
	}
      cam_gain_corrected->SetStats(0);
      if (display) c1->cd(2);
      //cam_gain_corrected->SetMinimum(1300);
      //cam_gain_corrected->SetMaximum(2000);
      if (display) cam_gain_corrected->Draw("COLZ");


      TH2D *cam_cleaned=new TH2D(*cam_2D_hist);
      cam_cleaned->SetName("cam_cleaned");
      snprintf(mytitle,200,"Event %i Cleaned Camera",i);
      cam_cleaned->SetTitle(mytitle);
      double Size=0;
      double SN[32][32];
      bool border[32][32];
      bool image[32][32];
      for (int i=0;i<32;i++)
	{
	  for (int j=0;j<32;j++)
	    {
	      pix_value[i][j]=cam_pix_data[i][j]-ped[i][j];
	      if (pedvar[i][j]!=0)
		{
		  SN[i][j]=pix_value[i][j]/pedvar[i][j];
		} else {
		SN[i][j]=0;
	      }
	      border[i][j]=false;
	      image[i][j]=false;
	      if (SN[i][j]>border_thresh && SN[i][j]<image_thresh) border[i][j]=true;
	      if (SN[i][j]>image_thresh) image[i][j]=true;
	    }	  
	}
      for (int i=0;i<32;i++)
	{
	  for (int j=0;j<32;j++)
	    {
	      int count_adjacent_image=0;
	      int count_adjacent_border=0;
	      if (!(border[i][j] || image[i][j]))
		{
		  pix_value[i][j]=0;
		}
	      else
		{
		  //set bounds to search neighbours
		  int minx=max(i-1,0);
		  int maxx=min(i+2,32);
		  int miny=max(j-1,0);
		  int maxy=min(j+2,32);
	      

		  for (int ii=minx;ii<maxx;ii++)
		    {
		      for (int jj=miny;jj<maxy;jj++)
			{
			  if (ii==i && jj==j) continue; 
			  if (image[ii][jj]) count_adjacent_image+=1;
			  if (border[ii][jj]) count_adjacent_border+=1;			  
			}
		    }
		}

	      if (border[i][j] && count_adjacent_image==0) pix_value[i][j]=0;
	      if (image[i][j] && (count_adjacent_border==0 && count_adjacent_image==0)) pix_value[i][j]=0;
	      if (pix_value[i][j]>0)
		{
		  cleaned_pix_value[i][j]=calibrated_pix_value[i][j];
		  cam_cleaned->SetBinContent(i+1,j+1,cam_gain_corrected->GetBinContent(i+1,j+1));
		  Size+=cam_gain_corrected->GetBinContent(i+1,j+1);
		} else {
		cam_cleaned->SetBinContent(i+1,j+1,0);
		cleaned_pix_value[i][j]=0;
	      }
	    }
	}
      
      cam_cleaned->SetStats(0);
      if (display) c1->cd(3);
      //cam_cleaned->SetMinimum(0);
      //cam_cleaned->SetMaximum(500);
      if (display) cam_cleaned->Draw("COLZ");
      
      
      
      
      HillasParams res = calculateHillas(cleaned_pix_value, 0.00000001);

      printf("Hillas Results:\n");
      printf("Size: %.2f | COG: (%.2f, %.2f)\n", res.size, res.cogX, res.cogY);
      printf("Length: %.2f | Width: %.2f\n", res.length, res.width);
      printf("Psi: %.2f deg | Alpha: %.2f deg\n\n", res.psi, res.alpha);

      outFile << res.size << " ";
      outFile << res.cogX << " ";
      outFile << res.cogY << " ";
      outFile << res.width << " ";
      outFile << res.length << " ";
      outFile << res.alpha << " ";
      outFile << endl;
      if (display)
      {
      TMarker *M=new TMarker(res.cogX+0.5, res.cogY+0.5, 29);
      M->SetMarkerColor(1);
      M->Draw("SAME");
      TEllipse *ellipse=new TEllipse(res.cogX+0.5, res.cogY+0.5, res.width,res.length,0,360,res.psi-90);
      ellipse->SetLineColor(1);
      ellipse->SetLineWidth(2);
      ellipse->SetFillStyle(0);
      ellipse->Draw("SAME");
      }

      if (display) c1->Modified();
      if (display) c1->Update();
      gSystem->ProcessEvents(); // Handle window interactions
      gSystem->Sleep(200);      // 100ms delay between frames	
      char a[200];
      if (display) (read(1,a,200));
    }
}



void draw_cam_ped_subtracted(int event)
{
  //TTree *camdata=loadcamdata(infile);
  camdata->GetEntry(event);
  TH2D *cam_ped_subtracted=new TH2D(*cam_2D_hist);
  cam_ped_subtracted->SetName("cam_ped_subtracted");
  cam_ped_subtracted->SetTitle("Pedestal Subtracted Camera");
  for (int i=0;i<32;i++)
    {
      for (int j=0;j<32;j++)
	{
	  cam_ped_subtracted->SetBinContent(i+1,j+1,(cam_pix_data[i][j]-ped[i][j]));
	}
    }
  cam_ped_subtracted->SetStats(0);
  cam_ped_subtracted->Draw("COLZ");
}

void draw_cam_sigma(int event)
{
  //TTree *camdata=loadcamdata(infile);
  camdata->GetEntry(event);
  TH2D *cam_sigma=new TH2D(*cam_2D_hist);
  cam_sigma->SetName("cam_sigma");
  cam_sigma->SetTitle("sigma Camera");
  for (int i=0;i<32;i++)
    {
      for (int j=0;j<32;j++)
	{
	  cam_sigma->SetBinContent(i+1,j+1,(cam_pix_data[i][j]-ped[i][j])/pedvar[i][j]);
	}
    }
  cam_sigma->SetStats(0);
  cam_sigma->Draw("COLZ");
}

