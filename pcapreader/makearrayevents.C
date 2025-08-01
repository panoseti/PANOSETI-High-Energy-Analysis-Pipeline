#include "./makearrayevents.h"

void findtimingoffsets(const char *infile)
{
   TTree *arraydata=define_arraydata();

  char outfile_name[200];
  snprintf(outfile_name,200,"%s.array",infile);
  cout << "writing to: " << outfile_name << endl;
 
  char cam_infile[200];
  int scope_id[ntel];
  scope_id[0]=12;
  scope_id[1]=760;
  scope_id[2]=1016;
  TTree * camdata[ntel];

  
  double toffset[ntel];
  double match_tolerence=0.1;

  for (int i=0;i<ntel;i++)
    {
      snprintf(cam_infile,200,"%s.T%d",infile,scope_id[i]);
      camdata[i]=load_camdata(cam_infile);
    }
    //first find the timing offsets between telescope clocks using matched events
  TH1D *htdiff=new TH1D("tdiff","tdiff",5000,-1*match_tolerence,match_tolerence);
  for (int i=0;i<camdata[0]->GetEntries();i++)
    {
      camdata[0]->GetEntry(i);
      double t1=cam_pcap_time_since_start;
      double tdiff=0;
      for (int j=0;j<camdata[1]->GetEntries();j++)
	{
	  camdata[1]->GetEntry(j);
	  //cam_pcap_time_since_start+= 0.00037;
	  tdiff=t1-cam_pcap_time_since_start;	  
	  if(fabs(tdiff)<match_tolerence)
	    {
	      htdiff->Fill(tdiff);
	      cout << "matched! " << i << " " << j <<  " t=" << cam_pcap_time_since_start  << " tdiff=" << tdiff << endl; 
	    }
	}
      
    }
  htdiff->Draw();
  cout <<  endl;
  cout << htdiff->GetMaximumBin() << endl;
  cout << htdiff->GetBinCenter(htdiff->GetMaximumBin()) << endl;
  return;
}

void makearrayevents(const char *infile)
{

  TTree *arraydata=define_arraydata();

  char outfile_name[200];
  snprintf(outfile_name,200,"%s.array",infile);
  cout << "writing to: " << outfile_name << endl;
 
  char cam_infile[200];
  int scope_id[ntel];
  scope_id[0]=12;
  scope_id[1]=760;
  scope_id[2]=1016;
  TTree * camdata[ntel];
  
  double toffset[ntel];
  double match_tolerence=0.1;

  
  for (int i=0;i<ntel;i++)
    {
      snprintf(cam_infile,200,"%s.T%d",infile,scope_id[i]);
      camdata[i]=load_camdata(cam_infile);
    }
  
  //add these offsets to the cam_pcap_time_since_start 
  // the numbers you put here are the numbers spit out by the findtimingoffsets function above - don't invert them
  //toffset[0]=0;
  //toffset[1]=-0.00149;
  //toffset[2]=0.00071;
  toffset[0]=0;
  toffset[1]=-0.00111;
  toffset[2]=+0.00037;
  match_tolerence=5.e-4;
  int max_nevents=0;
  double maxtime=0;
  int max_arrayevents=0;
  //first find the maximum number of events in the camera files
  for (int i=0;i<ntel;i++)
    {
      if (camdata[i]->GetEntries()>max_nevents) max_nevents=camdata[i]->GetEntries();

      camdata[i]->GetEntry(camdata[i]->GetEntries()-1);
      if ((cam_pcap_time_since_start+toffset[i])>maxtime) maxtime=cam_pcap_time_since_start+toffset[i];

      max_arrayevents+=camdata[i]->GetEntries();
    }
  double times[ntel][max_nevents];

  // now make an array of the offset-corrected event times
  for (int i=0;i<ntel;i++)
    {
      for (int j=0;j<max_nevents;j++)
	{
	  if (j<camdata[i]->GetEntries())
	    {
	      camdata[i]->GetEntry(j);
	      times[i][j]=cam_pcap_time_since_start+toffset[i];
	    }
	  else
	    {
	      times[i][j]=-100;
	    }
	}
    }

  // now implement an array trigger
  // this is slow and clunky, but I think it will work...

  int array_event[ntel][max_arrayevents];
  double trig_width=2*match_tolerence;
  double trig_step=trig_width/3.;

  double tmin=0;
  double tmax=trig_width;
  bool tend=false;
  
  int tsteps = 1+(int) (maxtime/trig_step);
  cout << "here " << ntel << " " << tsteps << endl;
  vector<vector<int> >   vevent_triggered;

  for (int it=0;it<tsteps;it++)
    {
      vector<int> vtelhits(ntel);
      bool triggered=0;
      for (int jtel=0;jtel<ntel;jtel++)
	{
	  vtelhits[jtel]=-100;
	  for (int k=0;k<max_nevents;k++)
	    {
	      if (times[jtel][k]>tmin && times[jtel][k]<tmax)
		{
		  vtelhits[jtel]=k;
		  triggered=1;
		}
	    }
	}
      if (triggered) vevent_triggered.push_back(vtelhits);
      tmin+=trig_step;
      tmax+=trig_step;
    }
  
  //now we just need to identify _unique_ array events
  TFile *outfile = TFile::Open(outfile_name,"RECREATE");

  int narrayevents=0;
  int n3folds=0;
  int n2folds=0;
  int n1folds=0;
  int n0folds=0;
  int prev_event[ntel];
  for (int jtel=0;jtel<ntel;jtel++) prev_event[jtel]=-100;
  int nmatch_prev=0;
  
  cout << "size vev " << vevent_triggered.size() << endl; 
  for (int i=0;i<vevent_triggered.size();i++)
    {
      int nteltrig=0;
      bool array_event_triggered=false;
      for (int jtel=0;jtel<ntel;jtel++)
	{
	  if (vevent_triggered[i][jtel]>0)
	    {
	      array_event_triggered=true;
	      nteltrig+=1;
	    }
	}
      if (array_event_triggered) {
	bool match_prev=false;
	if (vevent_triggered[i][0]==prev_event[0] && vevent_triggered[i][1]==prev_event[1] && vevent_triggered[i][2]==prev_event[2] ) {
	  match_prev=true;
	  nmatch_prev+=1;
	}
	cout << i << " " << vevent_triggered[i][0];
	cout << " " << vevent_triggered[i][1];
	cout << " " << vevent_triggered[i][2];	
	cout << " " << match_prev;
	cout << " p" << prev_event[0]; 
	cout << " p" << prev_event[1]; 
	cout << " p" << prev_event[2]; 
	cout << " " << nteltrig << endl;
	
	if (!match_prev)
	  {
	    
	    array_nteltrig=nteltrig;
	    for (int jtel=0;jtel<ntel;jtel++)
	      {
		if (vevent_triggered[i][jtel]==-100)
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
		else
		  {
		    camdata[jtel]->GetEntry(vevent_triggered[i][jtel]);
		    array_npix[jtel]=cam_npix;
		    array_tel_event[jtel]=vevent_triggered[i][jtel];
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
	      }
	    array_event_num=narrayevents;
	    arraydata->Fill();	    
	    narrayevents+=1;
	    if (array_nteltrig==0) n0folds+=1;
	    if (array_nteltrig==1) n1folds+=1;
	    if (array_nteltrig==2) n2folds+=1;
	    if (array_nteltrig==3) n3folds+=1;
	  }
	for (int jtel=0;jtel<ntel;jtel++)
	  prev_event[jtel]=vevent_triggered[i][jtel];
      }
    }
  cout << narrayevents << " unique array events found" << endl;
  cout << n0folds << " 0 telescope array events found" << endl;
  cout << n1folds << " 1 telescope array events found" << endl;
  cout << n2folds << " 2 telescope array events found" << endl;
  cout << n3folds << " 3 telescope array events found" << endl;
  cout << nmatch_prev << " match previous" << endl;
  cout << " writing to " << outfile_name << endl;
  
  arraydata->Write();
  outfile->Close();
}
