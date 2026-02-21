#include "./processcamevents.h"

void calcPedestals(const char *infile, bool do_fit=0)
{
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
	  double firstmean=h1->GetMean();
	  double firstRMS=h1->GetRMS();
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
	      h1->Fit("gaus");
	      TF1 *g = (TF1*)h1->GetListOfFunctions()->FindObject("gaus");
	      thisped=g->GetParameter(1);
	      thispedvar=g->GetParameter(2);
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
 
void eventdisplay(int start=0, int end=-1, bool wait=false)
{
  //If you calculated pedestals in a previous ROOT session, you have to do this first:
  //read_pedestals("../Fern/rawdata/Mrk421_preflip.root.pedvars")
  //camdata=loadcamdata("../Fern/rawdata/Mrk421_preflip.root")
  double image_thresh=4;

  TCanvas *c1 = new TCanvas("c1", "Sequential Image Display", 1200, 600);
  c1->Divide(2,1);
    
  if (end>camdata->GetEntries()) end=camdata->GetEntries();
  if (end<start) end=camdata->GetEntries();
  for (int i=start;i<end;i++)
    {
      camdata->GetEntry(i);

      TH2D *cam_ped_subtracted=new TH2D(*cam_2D_hist);
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
      c1->cd(1);
      cam_ped_subtracted->Draw("COLZ");

      TH2D *cam_cleaned=new TH2D(*cam_2D_hist);
      cam_cleaned->SetName("cam_cleaned");
      snprintf(mytitle,200,"Event %i Cleaned Camera",i);
      cam_cleaned->SetTitle(mytitle);
      for (int i=0;i<32;i++)
	{
	  for (int j=0;j<32;j++)
	    {
	      double pix_value=0;
	      pix_value=cam_pix_data[i][j]-ped[i][j];
	      double SN=0;
	      if (pedvar[i][j]!=0) SN=pix_value/pedvar[i][j];
	      if (SN<image_thresh) pix_value=0;
	      cam_cleaned->SetBinContent(i+1,j+1,pix_value);
	    }
	}
      
      cam_cleaned->SetStats(0);
      c1->cd(2);
      cam_cleaned->Draw("COLZ");
      
      
      c1->Modified();
      c1->Update();
      gSystem->ProcessEvents(); // Handle window interactions
      gSystem->Sleep(200);      // 100ms delay between frames	
      char a[200];
      if (wait) (read(1,a,200));
      
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

