#include "./makeevents.h"

TH2D *hdiff;
TH2D *hdiffcleaned;

void starfields(const char *infile1, const char *infile2)
{
  cout << "comparing before: " << infile1 << endl;
  cout << "with after:      " << infile2 << endl;
  cout << "difference is after - before      " << endl;

  TFile *_file0 = TFile::Open(infile1);
  TH2D * before=(TH2D*)_file0->Get("pedvars_2D_hist");
  //TH2D * before=(TH2D*)_file0->Get("peds_2D_hist");
  TH2D mybefore=*before;
  mybefore.SetTitle("Before");

  cout << "here OK" << endl;

  TFile *_file1 = TFile::Open(infile2);
  TH2D * after=(TH2D*)_file1->Get("pedvars_2D_hist");
  //TH2D * after=(TH2D*)_file1->Get("peds_2D_hist");
  TH2D myafter=*after;
  myafter.SetTitle("After");

  TCanvas *cjh=new TCanvas("cjh","cjh",1200,600);
  cjh->Draw();
  cjh->Divide(2,1);
  //cjh->cd(1);
  //before->Draw("COLZ");
  //cjh->cd(2);
  //after->Draw("COLZ");

  hdiff=new TH2D("hdiff","Difference",32,0,32,32,0,32);
  TH2D *hratio=new TH2D("hratio","Ratio",32,0,32,32,0,32);
  for (int i=0; i<32;i++)
    {
      for (int j=0; j<32;j++)
	{
	  double difference=myafter.GetBinContent(i+1,j+1) - mybefore.GetBinContent(i+1,j+1);

	  if (difference>0) difference=difference*difference;
	  if (difference<0) difference=-1*difference*difference;
	  double myratio=1;
	  if (mybefore.GetBinContent(i+1,j+1)!=0)
	    myratio=myafter.GetBinContent(i+1,j+1)/mybefore.GetBinContent(i+1,j+1);
	  hdiff->SetBinContent(i+1,j+1,difference);
	  hratio->SetBinContent(i+1,j+1,myratio);
	}
    }
  hdiff->SetStats(0);
  hdiffcleaned = (TH2D*)hdiff->Clone("hdiffcleaned");
  //hdiff->SetMinimum(-5);
  //hdiff->SetMaximum(5);
  hratio->SetStats(0);
  hratio->SetMinimum(0.9);
  hratio->SetMaximum(1.8);
  cjh->cd(1);
  hdiff->Draw("COLZ");
  cjh->cd(2);
  hratio->Draw("COLZ");
  //  cjh->cd();
  //  TLatex *t =new TLatex(10,0,infile1);
  //  t->Draw();
  cjh->SetTitle(infile1);
  cjh->Print("./test.png");
  /*
  TH2D *hdiff_flipX=flipX(hdiff);
  TH2D *hdiff_rotated180=flipY(hdiff_flipX);

  TH2D *hdiff2=new TH2D("hdiff2","Difference2",32,0,32,32,0,32);
  hdiff2->SetStats(0);
  cjh->cd(3);
  hdiff_rotated180->SetTitle("Difference rotated 180 deg");
  hdiff_rotated180->SetName("hdiff_rotated180");
  hdiff_rotated180->Draw("COLZ");
  //This works for T760
  //int xoffset=-2;
  //int yoffset=+4;
  //This works for T12
  //int xoffset=-3;
  //int yoffset=+1;
  //This works for T1016 unrotated
  int xoffset=-1;
  int yoffset=0;
  double difference=0;
  for (int i=0; i<32;i++)
    {
      for (int j=0; j<32;j++)
	{
	  if ((i+1+xoffset)<33 && (j+1+yoffset)<33 && (i+1+xoffset)>0 && (j+1+yoffset)>0)
	    {
	      //difference=hdiff->GetBinContent(i+1,j+1) - hdiff_rotated180->GetBinContent(i+1+xoffset,j+1+yoffset);
	      difference=hdiff->GetBinContent(i+1,j+1) - hdiff->GetBinContent(i+1+xoffset,j+1+yoffset);
	    }
	  else
	    {
	      difference = -100;
	    }
	  hdiff2->SetBinContent(i+1,j+1,difference);
	}
    }
  cjh->cd(4);
  hdiff2->SetMinimum(-30);
  hdiff2->SetMaximum(30);
  hdiff2->Draw("COLZ");
  */
}


int findandremovestar()
// this finds the brightest pixel in the map, and returns its bin number. It also zeros that bin and adjacent pixels in the cleaned map
{
  
  int maxbin=hdiffcleaned->GetMaximumBin();
  double maxbincontent=hdiffcleaned->GetBinContent(maxbin);
  int maxxbin,maxybin,maxzbin;
  hdiffcleaned->GetBinXYZ(maxbin,maxxbin,maxybin,maxzbin);
  
    
  cout << "max bin " << maxbin << " contains " << maxbincontent << " counts at X,Y bins " <<  maxxbin << "," << maxybin << endl;
  for (int i=0;i<3;i++)
    {
      for (int j=0;j<3;j++)
	{
	  
	  hdiffcleaned->SetBinContent(maxxbin-1+i,maxybin-1+j,0);
	}
    }
  
  return maxbin;
  
}

void fitstar(double xstar=-100,double ystar=-100)
{
  int starbin=findandremovestar();
  int starxbin,starybin,starzbin;
  hdiff->GetBinXYZ(starbin,starxbin,starybin,starzbin);
  double maxbincontent=hdiff->GetBinContent(starbin);
  double range=2.5;
  if (maxbincontent<5) range=1.5;
  cout << maxbincontent << " " << range << endl;
  TF2 *f2 = new TF2("f2", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[2])", starxbin-0.5-range, starxbin-0.5+range, starybin-0.5-range, starybin-0.5+range);
  cout << "here: " <<  starxbin-0.5 << " " << starybin-0.5 << endl ;
  f2->SetParameters(maxbincontent, starxbin-0.5, 0.3, starybin-0.5);
  f2->SetParNames("Constant", "MeanX", "Sigma", "MeanY");

// 4. Fit the Histogram
  hdiff->Fit(f2,"R");
  cout << f2->GetParameter(1) << " " <<  f2->GetParameter(3) << endl;

  cout << maxbincontent << " " << range <<  " " << f2->GetParameter(2) << " " << f2->GetParError(2) << endl;
  
  

}
