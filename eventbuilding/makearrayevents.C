#include "./makearrayevents.h"

void findtimingoffsets(const char *infile,int refscope=2, int testscope=0)
{
  //Find the timing offsets between telescope clocks using matched events
  // The output is analyzed using analyzetiming()
  
  TCanvas *mycanv=new TCanvas();
  mycanv->Divide(2);
  
  TTree *arraydata=define_arraydata();

  char cam_infile[200];
  char scope_name[ntel][100];
  snprintf(scope_name[0],100,"PTI");
  snprintf(scope_name[1],100,"Fern");
  snprintf(scope_name[2],100,"Winter");
  int scope_id[ntel];
  scope_id[0]=1000; // PTI
  scope_id[1]=1008; // Fern
  scope_id[2]=1012; // Winter
  TTree * camdata[ntel];

  char outfile_name[200];
  snprintf(outfile_name,200,"%s_T%d-T%d.timing",infile,scope_id[testscope],scope_id[refscope]);
  cout << "writing to: " << outfile_name << endl;
  
  double toffset[ntel];
  double match_tolerence=0.025;

  for (int i=0;i<ntel;i++)
    {
      //      snprintf(cam_infile,200,"%s/rawdata/%s.T%d",scope_name[i],infile,scope_id[i]);
      snprintf(cam_infile,200,"%s/rawdata/%s",scope_name[i],infile);
      cout << "here  " << cam_infile << endl;
      camdata[i]=load_camdata(cam_infile);
    }

  // find the start and end times of the data 
  double tstart;
  double tend;
  camdata[refscope]->GetEntry(0);
  tstart=cam_pcap_time;
  camdata[refscope]->GetEntry(camdata[refscope]->GetEntries()-1);
  tend=cam_pcap_time;
  int ntimebins=(tend-tstart)/60.;

  // placing the times for the reference telescope into a vector makes everything run alot quicker
  vector<double> refscopetime;
  for (int i=0;i<camdata[refscope]->GetEntries();i++)
    {
      camdata[refscope]->GetEntry(i);
      refscopetime.push_back(cam_pcap_time);
    }
  
  TH1D *htdiff=new TH1D("tdiff","tdiff",500,-1*match_tolerence,match_tolerence);
  TH2D *htdiff_t=new TH2D("tdiff_t","Time Offset vs Time",ntimebins,tstart, tend, 1250,-1*match_tolerence,match_tolerence);
  for (int i=0;i<camdata[testscope]->GetEntries();i++)
    {
      camdata[testscope]->GetEntry(i);
      double testscopetime=cam_pcap_time;
      double tdiff=0;
      for (int j=0;j<refscopetime.size();j++)
	{
	  tdiff=testscopetime-refscopetime[j];	  
	  if(fabs(tdiff)<match_tolerence)
	    {
	      htdiff->Fill(tdiff);
	      htdiff_t->Fill(refscopetime[j],tdiff);
	      
	      cout << "matched! " << i << " " << j <<  " t=" << cam_pcap_time-1.76847e9  << " tdiff=" << tdiff << endl; 
	    }
	}
    }
  TFile *outfile = TFile::Open(outfile_name,"RECREATE");
  htdiff->Draw();
  htdiff_t->SetXTitle("Time (s)");
  char ytitle[200];
  snprintf(ytitle,200,"Time offset %s(T%d)-%s(T%d)",scope_name[testscope],scope_id[testscope], scope_name[refscope],scope_id[refscope]);
  htdiff_t->SetYTitle(ytitle);
  cout <<  endl;
  cout << htdiff->GetMaximumBin() << endl;
  cout << htdiff->GetBinCenter(htdiff->GetMaximumBin()) << endl;
  htdiff_t->GetYaxis()->SetTitleOffset(1.4);
  htdiff_t->SetStats(0);
  htdiff_t->Draw();

  htdiff->Write();
  htdiff_t->Write();
  outfile->Close();
  return;
}

void makefakearrayevents(const char *infile, int real_tel=1)
{
  TTree *arraydata=define_arraydata();

  char outfile_name[200];
  snprintf(outfile_name,200,"%s.array",infile);
  cout << "writing to: " << outfile_name << endl;
 
  char cam_infile[200];
  char scope_name[ntel][100];
  snprintf(scope_name[0],100,"PTI");
  snprintf(scope_name[1],100,"Fern");
  snprintf(scope_name[2],100,"Winter");
  int scope_id[ntel];
  scope_id[0]=1000; // PTI
  scope_id[1]=1008; // Fern
  scope_id[2]=1012; // Winter
  TTree * camdata[ntel];
  
  double toffset[ntel];
  double match_tolerence;

  //I'm just going to use the same telescope data for all telescopes 
  for (int i=0;i<ntel;i++)
    {
      //snprintf(cam_infile,200,"%s/rawdata/%s",scope_name[real_tel],infile);
       camdata[i]=load_camdata(infile);
    }

  int max_nevents=camdata[real_tel]->GetEntries();
  
  TFile *outfile = TFile::Open(outfile_name,"RECREATE");
  
  for (int i=0;i<max_nevents;i++)
    {
      for (int jtel=0;jtel<ntel;jtel++)
	{
	  if (real_tel==jtel)
	    {
	      camdata[jtel]->GetEntry(i);
	      array_npix[jtel]=cam_npix;
	      array_tel_event[jtel]=real_tel;
	      array_scope_id[jtel]=cam_scope_id;
	      array_pcap_time[jtel]=cam_pcap_time;
	      array_pcap_time_since_start[jtel]=cam_pcap_time_since_start;
	      array_acq_mode[jtel]=cam_acq_mode;
	      array_TAI[jtel]=cam_TAI;
	      array_nanosec[jtel]=cam_nanosec;
	      array_WR_time[jtel]=cam_WR_time;
	      
	      for (int ip=0;ip<32;ip++)
		{
		  for (int jp=0;jp<32;jp++)
		    {
		      array_pix_data[jtel][ip][jp]=cam_pix_data[ip][jp];
		    }
		}
	    }
	  else
	    {
	      array_npix[jtel]=-999;
	      array_tel_event[jtel]=-999;
	      array_scope_id[jtel]=0;
	      array_pcap_time[jtel]=-999;
	      array_pcap_time_since_start[jtel]=-999;
	      array_acq_mode[jtel]=0;
	      array_TAI[jtel]=-999;
	      array_nanosec[jtel]=-999;
	      array_WR_time[jtel]=-999;
	      for (int ip=0;ip<32;ip++)
		{
		  for (int jp=0;jp<32;jp++)
		    {
		      array_pix_data[jtel][ip][jp]=-999;
		    }
		}	    
	    }
	}
      arraydata->Fill();	          
    }
  cout << " writing to " << outfile_name << endl;  
  arraydata->Write();
  outfile->Close();
}



void analyzetiming(const char *infile)
{
  // this analyzes the time offset vs time histograms produced by findtimingoffsets, and produces a smoothed graph of timing offset as a function of time, with points evenly spaced in time. 
  TFile *_file0 = TFile::Open(infile);
  TH2D *mytdiff_t=_file0->Get<TH2D>("tdiff_t");
  int npoints=0;
  vector<double> xp;
  vector<double> yp;
  
  for (int i=0;i<mytdiff_t->GetNbinsX();i++)
    {
      for (int j=0;j<mytdiff_t->GetNbinsY();j++)
	{
	  int neventsinbin=mytdiff_t->GetBinContent(i+1,j+1);
	  double binx=mytdiff_t->GetXaxis()->GetBinCenter(i+1);	  
	  double biny=mytdiff_t->GetYaxis()->GetBinCenter(j+1);
	  if (neventsinbin>=3)
	    {
	      xp.push_back(binx);
	      yp.push_back(biny);
	      npoints+=1;
	    }
	  
	}
    }
  cout << npoints << endl;
  TGraph *g=new TGraph(npoints);
  for (int i=0;i<npoints;i++)
    {
      cout << i << " " << xp[i] << " " << yp[i] << endl;
      g->SetPoint(i,xp[i],yp[i]);
    }

  //g->Draw("AW*");

      TGraphSmooth *gs=new TGraphSmooth("normal");
      TGraph *smoothedGraph = gs->SmoothLowess(g);

      double tmin=mytdiff_t->GetXaxis()->GetXmin();
      double tmax=mytdiff_t->GetXaxis()->GetXmax();
      int nbins_smoothed=(int) (tmax-tmin)/60.;
      gsmoothtime = gs->Approx(g,"linear",nbins_smoothed);
      gsmoothtime->SetMarkerStyle(21);
      gsmoothtime->SetMarkerSize(0.4);
      gsmoothtime->SetMarkerColor(2);
      gsmoothtime->Draw("AWP");
      g->Draw("*SAME");
      gsmoothtime->GetXaxis()->SetTitle("Time(s)");
      gsmoothtime->GetYaxis()->SetTitle("Time difference (s)");
      gsmoothtime->GetYaxis()->SetTitleOffset(1.3);
      gsmoothtime->SetTitle();

}
  
  
void correcttiming(const char * cam_infile)
{
// this uses the smoothed timing offset graph produced by analyzetiming to correct the times in the camera event trees.
// it produces a new file with a .corr extension, with timing corrections applied.
// for the telescope which is used as a reference, you should simply copy the camera event file and give it a .corr extension

  if (gsmoothtime!=NULL){
    gsmoothtime->Draw("AWP");
  } else {
    cout << "Failed. You must make the timing offset graph first using analyzetiming()" << endl;
    return;
  }
  TFile* originalFile = TFile::Open(cam_infile, "READ");
  TTree* camdata = (TTree*)originalFile->Get("camdata");
  char newfilename[200];
  snprintf(newfilename, 200,"%s.corr",cam_infile);
  TFile* newFile = TFile::Open(newfilename, "RECREATE");

  
  TTree *camdata_tcorr=(TTree*)camdata->CloneTree(0);
  
  double cam_pcap_time;
  double cam_pcap_time_since_start;
  camdata->SetBranchAddress("cam_pcap_time",&cam_pcap_time);
  camdata->SetBranchAddress("cam_pcap_time_since_start",&cam_pcap_time_since_start);
  //camdata_tcorr->SetBranchAddress("cam_pcap_time",&pcap_time_corr);
  //camdata_tcorr->SetBranchAddress("cam_pcap_time_since_start",&pcap_time_since_start_corr);
  vector<double> tcorrs;
  vector<double> t_tcorrs;
  for (int i=0;i<gsmoothtime->GetN();i++)
    {
      tcorrs.push_back(gsmoothtime->GetPointY(i));
      t_tcorrs.push_back(gsmoothtime->GetPointX(i));
    }

  int nentries = camdata->GetEntries();  
  for (int i = 0; i < nentries; i++) {
    camdata->GetEntry(i);
    double test_closest_time=9.e99;
    int closest_time_element;
    for (int j=0;j<tcorrs.size();j++)
      {
	if (fabs(cam_pcap_time-t_tcorrs[j])<test_closest_time)
	  {
	    test_closest_time=fabs(cam_pcap_time-t_tcorrs[j]);
	    closest_time_element=j;
	  }
      }
    cout << "correction:" << t_tcorrs[closest_time_element] << " " << tcorrs[closest_time_element] << endl;
    cam_pcap_time=cam_pcap_time-tcorrs[closest_time_element];
    cam_pcap_time_since_start=cam_pcap_time_since_start-tcorrs[closest_time_element];
    camdata_tcorr->Fill();
  }
    
  newFile->Write();
  newFile->Close();
  originalFile->Close();
}
