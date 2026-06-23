#include "./makeevents.h"



void plotaquabo(const char *infile, int packet, bool rotate)
// reads the pdata tree and plots a single quadrant board for a given data packet
{
  TTree *pdata=loadpdata(infile);
  pdata->GetEntry(packet);
  cout << "boardloc: " << boardloc << "\t time: " << pcap_time_since_start<< endl;
  TH2D *hq = makequabohist(pdata,rotate,0,8500);
  hq->Draw("COLZ");
}

TH2D *  averageaquabo(const char *infile, int quabo_id, bool rotate)
// calculates the average quabo output over many packets
// mostly useful for imaging mode
{
  char full_infile[200];
  snprintf(full_infile,200,"%s",infile);
  TTree *pdata=loadpdata(full_infile);
  pdata->GetEntry(0);
  TH2D *hq = makequabohist(pdata,rotate,0,8500);
  TH2D *hq_ave=(TH2D*)hq->Clone("hq_ave");
  hq_ave->Reset();
  char hq_ave_name[200];
  snprintf(hq_ave_name, 200,"QUABO %d %s",quabo_id, infile);
  hq_ave->SetName(hq_ave_name);
  hq_ave->SetTitle(hq_ave_name);
  
  int npackets=0;
  for (int i=10000;i<11000;i++)
    //for (int i=0;i<pdata->GetEntries();i++)
    {
      pdata->GetEntry(i);
      if (boardloc!=quabo_id) continue;
      npackets+=1;
      TH2D *hq = makequabohist(pdata,rotate,0,8500);
      for (int j=1;j<hq->GetNbinsX()+1;j++)
	{
	  for (int k=1;k<hq->GetNbinsY()+1;k++)
	    hq_ave->SetBinContent(j,k,hq_ave->GetBinContent(j,k)+hq->GetBinContent(j,k));
	}
    }
  cout << "npackets found=" << npackets << endl;
  
  for (int j=1;j<hq->GetNbinsX()+1;j++)
    {
      for (int k=1;k<hq->GetNbinsY()+1;k++)
	{
	  if (npackets>0) hq_ave->SetBinContent(j,k,hq_ave->GetBinContent(j,k)/(npackets*1.0));
	  //hq_ave->SetBinContent(j,k,hq_ave->GetBinContent(j,k)-hbg->GetBinContent(j,k));
	}
    }
  return hq_ave;
}

void makeavecameraimage(const char *infile, int scope)
{
  
  //char bg_filename[200];
  //snprintf(bg_filename,200,"d20241101/T%d_mid.root",scope);
  //cout << "opening " << bg_filename << endl;
  //TFile *_file0 = TFile::Open(bg_filename);
  //TCanvas *mycanv= (TCanvas*)_file0->Get("mycanv2");
  //TH2D *hbg=(TH2D*)mycanv->GetPrimitive("cam_2D_hist");
  //hbg->SetName("hbg");
  //hbg->SetTitle("hbg");
  //mycanv->Draw();
  //hbg->Draw("COLZ");
  
  
  TH2D *hq_ave;
  for (int i=0;i<4;i++)
    {
      int quabo_id=scope+i;
      hq_ave=averageaquabo(infile, quabo_id, 1);
      fillCamera(getCameraPosition(quabo_id),hq_ave);
    }
  for (int j=1;j<cam_2D_hist->GetNbinsX()+1;j++)
    {
      for (int k=1;k<cam_2D_hist->GetNbinsY()+1;k++)
	{
	  int bin_content=cam_2D_hist->GetBinContent(j,k);
	  //int bin_content=cam_2D_hist->GetBinContent(j,k)-hbg->GetBinContent(j,k);
	  //do this to fill with sqrt
	  //if (bin_content>0) bin_content=sqrt(bin_content);
	  //if (bin_content<0) bin_content=-1*sqrt(-1*bin_content);	  
	  cam_2D_hist->SetBinContent(j,k,bin_content);
	}
    }
  
  TCanvas *mycanv2= new TCanvas("mycanv2",infile,500,500);
  mycanv2->cd();
  cam_2D_hist->SetStats(0);
  cam_2D_hist->Draw("COLZ");
  //char full_outfile[200];
  //snprintf(full_outfile,200,"d20241031/ima/images_nobg/T%d/%s.png",scope,infile);
  //mycanv2->SaveAs(full_outfile);
}

void makeevents_softwaretrigger(const char *infile, int this_scope_id)
{
  //Let's try a simple brute force approach to QUABO matching, which is a software coincidence trigger between exactly four modules within a very short coincidence window. Everything else gets discarded. Since the coincidence time is so short (<20ns), I am going to allow duplicates (i.e. a QUABO packet can be used twice). But this should very rarely happen.
  //WR_time is TAI+nanosec. This seems to clock over at 1024 (10 bits). But so long as it's counting, and the run is less than 1024 seconds long (17 minutes) that should be OK for this application. 

  TTree *pdata=loadpdata(infile);    
  int npdata=pdata->GetEntries();
  pdata->GetEntry(npdata-1);
  double duration=pcap_time_since_start+1;
  cout << "Run duration = " << duration << endl;					   
  if (duration>1023) cout << "CAUTION: run duration is longer than 1023 seconds. This QUABO matching algorithm may get confused (TAI time only uses 10 bits = 1024s)" << endl; 

  double WR_time;
  pdata->SetBranchAddress("WR_time", &WR_time);
   // 3. Extract and fill individual time vectors
    std::vector<double> ts1, ts2, ts3, ts4;    
    for (Long64_t i = 0; i < npdata; ++i) {
        pdata->GetEntry(i);
        if (boardloc-this_scope_id==0) ts1.push_back(WR_time);
        if (boardloc-this_scope_id==1) ts2.push_back(WR_time);
        if (boardloc-this_scope_id==2) ts3.push_back(WR_time);
        if (boardloc-this_scope_id==3) ts4.push_back(WR_time);
    }
    //inFile->Close(); // Close input file since data is in memory
    std::sort(ts1.begin(), ts1.end());
    std::sort(ts2.begin(), ts2.end());
    std::sort(ts3.begin(), ts3.end());
    std::sort(ts4.begin(), ts4.end());

    std::vector<double> mt1,mt2,mt3,mt4;
      
    // 6. Sliding window coincidence search
    const Double_t window = 20e-9; // 20 ns window
    size_t i2 = 0, i3 = 0, i4 = 0;
    long long matchCount = 0;
  
    // I used Gemini for this bit. So it's confusing, but fast... 
    for (size_t i1 = 0; i1 < ts1.size(); ++i1) {
        Double_t t1 = ts1[i1];

        // Advance baseline pointers to minimize search iterations
        while (i2 < ts2.size() && ts2[i2] < t1 - window) i2++;
        while (i3 < ts3.size() && ts3[i3] < t1 - window) i3++;
        while (i4 < ts4.size() && ts4[i4] < t1 - window) i4++;

        // Evaluate nearby candidate events
        for (size_t k2 = i2; k2 < ts2.size() && ts2[k2] <= t1 + window; ++k2) {
            for (size_t k3 = i3; k3 < ts3.size() && ts3[k3] <= t1 + window; ++k3) {
                for (size_t k4 = i4; k4 < ts4.size() && ts4[k4] <= t1 + window; ++k4) {
                    
                    Double_t min_t = std::min({t1, ts2[k2], ts3[k3], ts4[k4]});
                    Double_t max_t = std::max({t1, ts2[k2], ts3[k3], ts4[k4]});

                    if ((max_t - min_t) <= window) {
		      mt1.push_back(t1);
		      mt2.push_back(ts2[k2]);
		      mt3.push_back(ts3[k3]);
		      mt4.push_back(ts4[k4]);
		      //out_t3 = ts3[k3];
		      //out_t4 = ts4[k4];
		      //out_dt = max_t - min_t;
		      
		      //outTree->Fill();
		      matchCount++;
                    }
                }
            }
        }
    }
    std::cout << "Found " << matchCount << " coincidences." << std::endl;

    //OK - I've matched up the QUABO packets into camera events.
    //Now create those camera events
    char outfile_name[200];
    snprintf(outfile_name,200,"%s.T%d",infile,this_scope_id);
    cout << "Done. Writing to: " << outfile_name << endl;
    TFile *root_outfile = TFile::Open(outfile_name,"RECREATE");
    TTree *camdata=define_camdata();
    TH2D *hquabo;

    // Clever ROOT trick: Maps "WR_time" branch directly to entry number - allows quick lookup
    pdata->BuildIndex("WR_time");

    for (size_t i = 0; i < mt1.size(); ++i)    
      {
	//I just want to loop over QUABO packet entries close in time to this event  
	int this_entry = pdata->GetEntryNumberWithBestIndex(mt1[i]);
	int start_entry=this_entry-10;
	if (start_entry<0) start_entry=0;
	int end_entry=this_entry+10;
	if (end_entry>npdata) end_entry=npdata;      
	//cout << i << " " << start_entry << " " << end_entry << endl;
							       
	for (int j=start_entry;j<end_entry;j++)
	  {	    
	    pdata->GetEntry(j);
	    bool found1=0,found2=0,found3=0,found4=0;
	    if (WR_time==mt1[i] && boardloc-this_scope_id==0)
	      {
		UShort_t scope_id=getScopeID(boardloc);		
		if (scope_id!=this_scope_id) continue;
		// assign event time stamps etc based on the first quabo
		cam_pcap_time=pcap_time;
		cam_pcap_time_since_start=pcap_time_since_start;
		cam_acq_mode=acq_mode;
		cam_TAI=TAI;
		cam_nanosec=nanosec;
		cam_WR_time=WR_time;
		cam_scope_id=scope_id;
		hquabo = makequabohist(pdata,1,-30,100); //puts the quabo packet into a 2D hist and rotates it based on its position in the camera
		fillCamera(getCameraPosition(boardloc),hquabo);
		found1=1;
	      }
	    if (WR_time==mt2[i] && boardloc-this_scope_id==1)
	      {
		hquabo = makequabohist(pdata,1,-30,100); //puts the quabo packet into a 2D hist and rotates it based on its position in the camera
		fillCamera(getCameraPosition(boardloc),hquabo);
		found2=1;
	      }
	    if (WR_time==mt3[i] && boardloc-this_scope_id==2)
	      {
		hquabo = makequabohist(pdata,1,-30,100); //puts the quabo packet into a 2D hist and rotates it based on its position in the camera
		fillCamera(getCameraPosition(boardloc),hquabo);
		found3=1;
	      }
	    if (WR_time==mt4[i] && boardloc-this_scope_id==3)
	      {
		hquabo = makequabohist(pdata,1,-30,100); //puts the quabo packet into a 2D hist and rotates it based on its position in the camera
		fillCamera(getCameraPosition(boardloc),hquabo);
		found4=1;
	      }
	    if (found1 && found2 && found3 && found4) break;
	
	  }
	camdata->Fill();
	//cout << "found event " << i << endl;
				       
      }
    camdata->Write();    
    root_outfile->Close();
    cout << endl;
}




// There are three scopes in the triangle data. IDs are: 760 (dorm), 12 (kron) and 1016 (crocker)
// This method matches quabo packets together using a time histogram approach. Simple and robust, but might struggle at high rates.
void makeevents_timehistomatch(const char *infile, int this_scope_id)
{

  TTree *pdata=loadpdata(infile);    
  int npdata=pdata->GetEntries();
  pdata->GetEntry(npdata-1);
  double duration=pcap_time_since_start;
  //How many histogram bins do we need such that there are a maximum of 4 packets per bin?
    int nbins=(int)(duration*100);
    TCanvas *mycanv=new TCanvas();
    mycanv->cd();
    TH1D *h1;
    for (int i=0;i<10;i++)
      {
	h1=new TH1D("h1","h1",nbins,0,duration);
	h1->SetXTitle("Time (s)");
	h1->SetYTitle("Number of packets");
	h1->SetStats(0);
	char cut[200];
	snprintf(cut,200,"boardloc>=%d && boardloc<%d+4",this_scope_id,this_scope_id);
	pdata->Draw("pcap_time_since_start>>h1",cut);
	if (h1->GetMaximum()<=4) break;
	nbins*=2;
      }
    //nbins*=10;
    cout << "using " << nbins << " bins in the rate histogram " << h1->GetNbinsX() << endl;

    // these vectors bracket the start and end times of the camera event, there should be the same number of elements for each  
    vector<double> event_tstart;
    vector<double> event_tend;
    vector<int> event_nquabos;

    //there are lots of ways the data can be corrupted. For example:
    //An event may not have all four quabos read-out
    //A single quabo might read out multiple times over a very short time interval (i.e. a burst)
    //Quabo packets may not arrive in sequential order
    for (int i=0;i<nbins;i++)
      {
	int nquabos=h1->GetBinContent(i);
	if (nquabos>0) cout << "in here " << i << " " << nquabos<< endl; 
	if (nquabos>0)
	  {
	    //cout << "initial number of quabos= "<< nquabos << " at t=" << h1->GetBinLowEdge(i) << endl;
	    event_tstart.push_back(h1->GetBinLowEdge(i));
	    if (nquabos==4) {
	      //cout << "A nice clean event found with 4 quabos" << endl;
	      event_tend.push_back(h1->GetBinLowEdge(i)+h1->GetBinWidth(i));
	      event_nquabos.push_back(nquabos);
	      continue;
	    }
	    // if I'm here, there must be less than 4 quabos in this bin - search over following 10 bins for additional quabos	    
	    cout << "initial number of quabos= "<< nquabos << " at t=" << h1->GetBinLowEdge(i) << endl;
	    for (int j=0;j<10;j++)
	      {
		int next_bin=i+j+1;
		if (next_bin<nbins)
		  nquabos+=h1->GetBinContent(next_bin);
		if (nquabos==4) {
		  // OK - I found enough quabos nearby to make a full camera. Let's assume these are part of the same event
		  cout << " adding " << h1->GetBinContent(next_bin)  << " gives event found with " << "4 quabos ending at " << h1->GetBinLowEdge(next_bin)+h1->GetBinWidth(next_bin) << endl;
		  event_tend.push_back(h1->GetBinLowEdge(next_bin)+h1->GetBinWidth(next_bin));
		  event_nquabos.push_back(nquabos);
		  i=next_bin;
		  break;
		}
		if (nquabos>4){
		  //Additional quabos cannot be associated with this event - there are too many (i.e. this event must have less than 4 quabos)
		  cout << "event found with " << nquabos << " quabos. Too many! reverting to an incomplete event ending at " << h1->GetBinLowEdge(i)+h1->GetBinWidth(i) << endl;
		  //event_tend.push_back(h1->GetBinLowEdge(i)+h1->GetBinWidth(i));
		  nquabos-=h1->GetBinContent(next_bin);
		  //event_nquabos.push_back(nquabos);
		  break;
		}
	      }
	    if (nquabos<4){
	      cout << "Not enough additional quabos found. Incomplete event" << endl;
	      event_tend.push_back(h1->GetBinLowEdge(i)+h1->GetBinWidth(i));
	      event_nquabos.push_back(nquabos);
	    }
	  }
      }
    
    cout  << "nquabos size "<< event_nquabos.size() << " tstart size " << event_tstart.size() << " tend size " << event_tend.size() << endl;
    for (int i=0;i<event_tstart.size(); i++)
      {
	if (event_nquabos[i] !=4)
	  // at this stage, events with fewer than 4 quabos seem to be either genuine events missing a quabo (OK), or duplicate reads of the same quabo (Bad)
	  // let's look at the board_ids to determine which.
	  cout << i << " nquabos " << event_nquabos[i] << " tstart " << event_tstart[i] << " tend " << event_tend[i] << " diff=" << event_tend[i] - event_tstart[i] << endl;
	if (event_nquabos[i]>4) cout << "more than 4 quabos. This should never happen." << endl;
      }
    
    //OK - I've matched up the QUABO packets into camera events.
    //Now create those events
    char outfile_name[200];
    snprintf(outfile_name,200,"%s.T%d",infile,this_scope_id);
    cout << "writing to: " << outfile_name << endl;
    TFile *root_outfile = TFile::Open(outfile_name,"RECREATE");
    TTree *camdata=define_camdata();
    TH2D *hquabo;
    int start_bin=0;

    
    for (int i=0;i<event_tstart.size(); i++)
      {
	if (event_nquabos[i]<3) continue; // not worth saving camera events with less than 3 quabos;
	int nquabos_found=0;
	vector<int> boards; // list of quabo board numbers in this event  
	for (int j=0;j<pdata->GetEntries();j++)
	  {
	    pdata->GetEntry(j);
	    UShort_t scope_id=getScopeID(boardloc);
	    if (scope_id!=this_scope_id) continue;
	    if (pcap_time_since_start<event_tstart[i]) continue;
	    if (pcap_time_since_start>event_tend[i]) break;
	    nquabos_found+=1;
	    if (nquabos_found==1)
	      {
		// assign event time stamps etc based on the first quabo
		cam_pcap_time=pcap_time;
		cam_pcap_time_since_start=pcap_time_since_start;
		cam_acq_mode=acq_mode;
		cam_TAI=TAI;
		cam_nanosec=nanosec;
		cam_WR_time=WR_time;
		cam_scope_id=scope_id;
	      }
	    hquabo = makequabohist(pdata,1,-30,100); //puts the quabo packet into a 2D hist and rotates it based on its position in the camera
	    fillCamera(getCameraPosition(boardloc),hquabo);
	    boards.push_back(boardloc);
	  }
	//this identifies events with duplicate quabos
	std::sort(boards.begin(), boards.end());
	bool duplicate=0;
	for (int j=1;j<boards.size();j++) {
	  if (boards[j]==boards[j-1]) duplicate=1;
	}

	if (!duplicate) camdata->Fill();
      }
    camdata->Write();    
    root_outfile->Close();
    cout << endl;
    delete mycanv;
    
}	    

    
// there are three scopes in the triangle data. IDs are: 760 (dorm), 12 (kron) and 1016 (crocker)
void makeevents_packetmatch(const char *infile, int this_scope_id)
{
  TTree *pdata=loadpdata(infile);
    
  int npdata=pdata->GetEntries();
  cout << "Number of QUABO board packets: " << npdata << endl;

  int camera_events[npdata][4];
  double tcamera_events[npdata][4];
  int ncamera_events=0;
  
  int packet_numbers_start[4];
  double tstart[4];
  bool packet_numbers_start_read[4];
  int packet_numbers_end[4];
  double t_last[4];
  double pn_last[4];
  int duplicates[4];
  for (int i=0;i<4;i++)
    {
      packet_numbers_start_read[i]=0;
      duplicates[i]=0;
    }
    
  // searching through all of the packets, looking for packets just from this scope
  for (int i=0;i<npdata;i++)
    {
      pdata->GetEntry(i);
      UShort_t scope_id=getScopeID(boardloc);
      UShort_t board_id=boardloc-scope_id;
      if (scope_id!=this_scope_id) continue;

      //record some starting values
      if (!packet_numbers_start_read[board_id])
	{
	  packet_numbers_start[board_id]=packet_no;
	  packet_numbers_start_read[board_id]=true;
	  tstart[board_id]=pcap_time_since_start;
	  camera_events[0][board_id]=i;
	}
      else
	{
	  cout << "board_id" << board_id << "\t packet no. - start" << packet_no - packet_numbers_start[board_id] << "\t tdiff" << pcap_time_since_start-t_last[board_id] << "\t packet no. diff" << packet_no - pn_last[board_id] << endl;
	  if  (fabs(pcap_time_since_start-t_last[board_id])<5.e-4)
	    {
	      cout << "found a duplicate packet at event " << i << " ; for board " << board_id << " ; at packet number " <<  packet_no << " which is " << packet_no - packet_numbers_start[board_id] << " from start" << endl;
	      duplicates[board_id]+=1;
	    }
	}
      t_last[board_id]=pcap_time_since_start;
      pn_last[board_id]=packet_no;
      
      packet_numbers_end[board_id]=packet_no;
      //cout << board_id << "\t" << packet_no - packet_numbers_start[board_id] << "\t" << pcap_time_since_start - tstart[board_id] << "\t" << pix_data_signed[10] << "\t" << pix_data_signed[100]<< endl;
      int this_camera_event=packet_no - packet_numbers_start[board_id]-duplicates[board_id];
      camera_events[this_camera_event][board_id]=i;      
      tcamera_events[this_camera_event][board_id]=pcap_time_since_start;      
      ncamera_events=this_camera_event;
    }
  cout << endl;

  double sum=0;
  for (int i=0;i<4;i++)
    {
      cout << i << "\t" << packet_numbers_start[i] << endl;
      cout << i << "\t" << packet_numbers_end[i] << endl;
      cout << i << "\t" << packet_numbers_end[i]-packet_numbers_start[i] << endl;
      sum+=packet_numbers_end[i]-packet_numbers_start[i];
      cout << i << "duplicates: " << duplicates[i] << endl;
    }
  cout << "sum=: " << sum << endl;
  
  for (int i=0;i<ncamera_events;i++)
    {
      
      
      
      double timesum=0;
      for (int j=0;j<4;j++)
	{
	  timesum+=tcamera_events[i][j];
	}

      double timemean=timesum/4;
      double timevar=0;
      for (int j=0;j<4;j++)
	{
	  timevar+=(tcamera_events[i][j]-timemean)*(tcamera_events[i][j]-timemean);
	}
      timevar/=3.0;
      if (timevar>1.e-5) {
	cout << "WARNING: large time variance: " << timevar << endl;
	      cout << i << " " << camera_events[i][0];
	      cout << " " << camera_events[i][1];
	      cout << " " << camera_events[i][2];
	      cout << " " << camera_events[i][3] << endl;	
	      
	      cout << i << " " << tcamera_events[i][0];
	      cout << " " << tcamera_events[i][1];
	      cout << " " << tcamera_events[i][2];
	      cout << " " << tcamera_events[i][3] << endl;	
	      cout << endl;
      }
      
    }
  if ((sum/4.0-(int)(sum/4.0))!=0) {
    cout << "QUABO packet numbers do not match. sum remainder=" << (sum/4.0-(int)(sum/4.0)) << endl;    
  }

  //OK - I've matched up the QUABO packets into camera events.
  //Now create those events
  char outfile_name[200];
  snprintf(outfile_name,200,"%s.T%d",infile,this_scope_id);
  cout << "writing to: " << outfile_name << endl;
  TFile *outfile = TFile::Open(outfile_name,"RECREATE");
  TTree *camdata=define_camdata();
  TH2D *hquabo;
  for (int i=0;i<ncamera_events;i++)
    {
      for (int j=0;j<4;j++)
	{
	  pdata->GetEntry(camera_events[i][j]);
	  UShort_t scope_id=getScopeID(boardloc);
	  if (scope_id!=this_scope_id) continue;
      
	  //cout << camera_events[i][j] << " this packet is telescope ID " << scope_id << endl;

	  cam_pcap_time=pcap_time;
	  cam_pcap_time_since_start=pcap_time_since_start;
	  cam_acq_mode=acq_mode;
	  cam_TAI=TAI;
	  cam_nanosec=nanosec;
	  cam_WR_time=WR_time;
	  cam_scope_id=scope_id;
      
	  //cout << "boardloc: " << boardloc << "\t time: " << pcap_time_since_start<< endl;
	  hquabo = makequabohist(pdata,1,-30,100); //puts the quabo packet into a 2D hist and rotates it based on its position in the camera
	  fillCamera(getCameraPosition(boardloc),hquabo);


	}
      camdata->Fill();
    }
  camdata->Write();
  gROOT->Reset();
}

void makeevents_timematch(const char *infile, int this_scope_id)
{
  TTree *pdata=loadpdata(infile);
    
  int npdata=pdata->GetEntries();
  cout << "Number of QUABO board packets: " << npdata << endl;
  TFile *outfile = TFile::Open("test.root","RECREATE");

  TTree *camdata=define_camdata();
  TH2D *hquabo;
  double tstart=-100;
  UShort_t pos_in_cam=99;

  bool used[npdata+1000];
  for (int i=0;i<npdata+1000;i++) used[i]=0;

  int max_search_range_packets=100; // how far forward from this packet am I prepared to search (in packets);
  double max_search_range_time=1; // how far forward from this packet am I prepared to search (in seconds);

  //for (int packet=0;packet<npdata; packet++)
  for (int packet=0;packet<npdata; packet++)
    {
      if (used[packet]) {
	cout << "This packet has already been used: " << packet << endl;
	continue;
      }
      cout << "now working on packet " << packet << endl;
      pdata->GetEntry(packet);
      UShort_t scope_id=getScopeID(boardloc);
      if (scope_id!=this_scope_id) continue;
      
      cout << "this is telescope ID " << scope_id << endl;

      cam_pcap_time=pcap_time;
      cam_pcap_time_since_start=pcap_time_since_start;
      cam_acq_mode=acq_mode;
      cam_TAI=TAI;
      cam_nanosec=nanosec;
      cam_WR_time=WR_time;
      
      cout << "boardloc: " << boardloc << "\t time: " << pcap_time_since_start<< endl;
      hquabo = makequabohist(pdata,1,-30,100); //puts the quabo packet into a 2D histogram and rotates it based on its position in the camera

      //hquabo->Draw("COLZ");

      //Now we identify which QUABO packets correspond to the same event
      // Approach: first make individual telescope events. Run through the whole file for each scope separately. Deal with array events later. Slow, but robust and less confusing. Worry about speed later.  
      // Algorithm: get the time of a QUABO packet. Then search through packets until a packet from the same telescope arrives with a time difference greater than some small value. Or a packet from the same quabo arrives. Or the end of file is reached. (could also check for packets from different telescopes with a time difference greater than some larger value - leave that for now).
 
      int upper_bound_packets=packet+max_search_range_packets;

      double max_tdiff=0.004; // this is the maximum time difference between QUABO packets within the same camera for an event
      if (acq_mode==3) max_tdiff=0.01; // might need to make this time difference acquisition-mode-dependent
      initializeCamera();
      bool foundfirst=false;
      if (getScopeID(boardloc)==this_scope_id)
	{
	  fillCamera(getCameraPosition(boardloc),hquabo);
	  pos_in_cam=getCameraPosition(boardloc); 
	  tstart=pcap_time;
	  cam_scope_id=scope_id;
	  cout << "the first QUABO found is board " << boardloc << " at " << tstart << endl;
	  foundfirst=true;
	}

      if (upper_bound_packets>npdata) upper_bound_packets=npdata;
      for (int i=packet+1;i<upper_bound_packets;i++)
	{
	  pdata->GetEntry(i);
	  hquabo = makequabohist(pdata,1,-30,100); //puts the quabo packet into a 2D histogram and rotates it based on its position in the camera
	  if (!foundfirst) {
	    cout << "first quabo not found yet" << endl;
	    break;
	  }
	  if (getScopeID(boardloc)==this_scope_id)
	    {
	      cout<< "this packet is from the requested telescope" << endl;
	      double tdiff=fabs(pcap_time-tstart);
	      if (getCameraPosition(boardloc) == pos_in_cam) {
		cout << "stopping the search: found another packet from the same QUABO " << getCameraPosition(boardloc) <<" at tdiff=" << tdiff<<endl;
		break;
	      }
	      if (tdiff>max_tdiff){
		cout << "stopping the search: time difference is too great: tdiff=" << tdiff << endl;
		break;
	      }
	      cout << "found QUABO " << getCameraPosition(boardloc) << " at tdiff=" << tdiff<< endl;
	      fillCamera(getCameraPosition(boardloc),hquabo);	  
	      used[i]=true;
	      cout << "i=" << i << endl;
	    }
	}
      camdata->Fill();
      //char htitle[20];
      //snprintf(htitle,20,"TELESCOPE %d",cam_scope_id);
      //cam_2D_hist->SetTitle(htitle);
      //cam_2D_hist->Draw("COLZ");
      //cin.get();
    }
  camdata->Write();
  //delete outfile;

}

void ProcessFiles(const char* listFileName = "filelist.dat",int scope_id=0) {
    std::ifstream inputFile(listFileName);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Cannot open list file " << listFileName << std::endl;
        return;
    }

    std::string line;
    // Loop through each line (filename) in the text file
    while (std::getline(inputFile, line)) {
        // Remove potential leading/trailing whitespace or newline characters
        // (basic example; robust handling might need TString methods)
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        if (!line.empty() && line[0] != '#') { // Ignore empty lines and comments
	  //makeevents_timehistomatch(line.c_str(),scope_id);
	  makeevents_softwaretrigger(line.c_str(),scope_id);
        }
    }

    inputFile.close();
}
