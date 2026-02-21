#include "./makeevents.h"

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

  TCanvas *cjh=new TCanvas("cjh","cjh",800,800);
  cjh->Draw();
  cjh->Divide(2,2);
  //cjh->cd(1);
  //before->Draw("COLZ");
  //cjh->cd(2);
  //after->Draw("COLZ");

  TH2D *hdiff=new TH2D("hdiff","Difference",32,0,32,32,0,32);
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
  //hdiff->SetMinimum(-5);
  //hdiff->SetMaximum(5);
  hratio->SetStats(0);
  hratio->SetMinimum(0.7);
  hratio->SetMaximum(1.7);
  cjh->cd(1);
  hdiff->Draw("COLZ");
  cjh->cd(2);
  hratio->Draw("COLZ");
  //  cjh->cd();
  //  TLatex *t =new TLatex(10,0,infile1);
  //  t->Draw();
  cjh->SetTitle(infile1);
  cjh->Print("./test.png");

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
  
}

