void gain_from_pedvars(const char *prefile, const char *postfile, const char *gainfile)
// This script calculates pixel-pixel relative gains using the pedestal variances.
// It minimizes the impact of bright stars by using two star fields (e.g. Markarian 421, pre- and post-meridian flip).
// If a pixel is bright in one field, but not the other, the lower value is used.
// Results are saved to a single 2D relative histogram in a .root.gain file
{
  TFile *fpreflip = TFile::Open(prefile);
  TH2D *pedvars_2D_hist_pre = (TH2D*)gDirectory->Get("pedvars_2D_hist");
  double pedvars_preflip[32][32];
  for (int i=0;i<32;i++)
    {
      for (int j=0;j<32;j++)
	{
	  pedvars_preflip[i][j]=pedvars_2D_hist_pre->GetBinContent(i+1,j+1);
	}
    }
  fpreflip->Close();
  
  TFile *fpostflip = TFile::Open(postfile);
  TH2D *pedvars_2D_hist_post = (TH2D*)gDirectory->Get("pedvars_2D_hist");
  double pedvars_postflip[32][32];
  for (int i=0;i<32;i++)
    {
      for (int j=0;j<32;j++)
	{
	  pedvars_postflip[i][j]=pedvars_2D_hist_post->GetBinContent(i+1,j+1);
	}
    }
  fpostflip->Close();
  double vratio[32][32];
  TGraph *g=new TGraph(1024);
  TH1D *hratio=new TH1D("hratio","Pedestal variance ratio (postflip/preflip)",40,0.5,1.5);
    for (int i=0;i<32;i++)
    {
      for (int j=0;j<32;j++)
	{
	  double x=pedvars_preflip[i][j];
	  double y=pedvars_postflip[i][j];
	  double ratio=1;
	  if (y>0) ratio=y/x;
	  vratio[i][j]=ratio;
	  hratio->Fill(ratio);
	  g->SetPoint(i*32+j,x,y);
	}
    }
    g->SetTitle("Pedestal Variance comparison for PTI-Heli. Mrk 421 Jan 15th 2026");
    g->GetXaxis()->SetTitle("PedVar before meridian flip");
    g->GetYaxis()->SetTitle("PedVar after meridian flip");
    g->SetMarkerStyle(7);
    hratio->SetLineWidth(3);
    hratio->SetXTitle("Pedestal variance ratio (postflip/preflip)");
    hratio->Draw();
    double brightpix_upper=hratio->GetMean()+3*hratio->GetRMS();
    double brightpix_lower=hratio->GetMean()-3*hratio->GetRMS();

    TH1D *hrelgain=new TH1D("hrelgain","Relative gain",50,0.0,2.0);

    double relgain[32][32];
    double mean=0;
    double count=0;
    for (int i=0;i<32;i++)
      {
	for (int j=0;j<32;j++)
	  {
	    relgain[i][j]=(pedvars_preflip[i][j]+pedvars_postflip[i][j])/2.0; // no stars - use the average
	    if (vratio[i][j]>brightpix_upper) relgain[i][j]=pedvars_preflip[i][j]; // star in the postflip
	    if (vratio[i][j]<brightpix_lower) relgain[i][j]=pedvars_postflip[i][j]; // star in the preflip
	    mean+=relgain[i][j];
	    count+=1;
	  }
      }
    mean=mean/count;
    cout << "Mean: " << mean << endl;
    for (int i=0;i<32;i++)
      {
	for (int j=0;j<32;j++)
	  {
	    relgain[i][j]=relgain[i][j]/mean;
	    hrelgain->Fill(relgain[i][j]);
	  }
      }
    hrelgain->SetLineWidth(3);
    hrelgain->Draw();

    TFile *fgain = TFile::Open(gainfile,"RECREATE");

    TH2D *relgain_2D_hist=new TH2D("relgain_2D_hist","Relative Gain",32,0,32,32,0,32);
    for (int i=0;i<32;i++)
      {
	for (int j=0;j<32;j++)
	  {
	    relgain_2D_hist->SetBinContent(i+1,j+1,relgain[i][j]);
	  }
      }
    relgain_2D_hist->SetStats(0);
    relgain_2D_hist->SetMinimum(0.5);
    relgain_2D_hist->SetMaximum(1.6);
    
    relgain_2D_hist->Draw("COLZ");
    relgain_2D_hist->Write();
    //fgain->Close();
}
