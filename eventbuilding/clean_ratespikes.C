// script to clean rate spikes from PANOSETI packet data
// all packets occurring at a rate greater than maxrate are removed from the output file

#include "./makeevents.h"

//maxrate is the maximum number of events per bin (i.e. it is not in Hz)



void clean_ratespikes(const char *infile, double maxrate)
{
  double bin_duration=0.1;
  TCanvas *mycanv=new TCanvas();
  mycanv->Divide(1,2);
  TTree *pdata=loadpdata(infile);
  char outfile_name[200];
  snprintf(outfile_name,200,"%s.cleaned",infile);
  TFile *newfile = new TFile(outfile_name,"recreate");
  TTree *pdataclean = pdata->CloneTree(0);
  pdata->GetEntry(0);
  double tstart=0; // note that pcap_time_since_start of first entry may not be zero - so this is safer.
  pdata->GetEntry(pdata->GetEntries()-1);
  double tend=pcap_time_since_start;
  int nbins=(int)((tend-tstart)/bin_duration); 
  cout << pdata->GetEntries() << " " << tstart << " " << tend << " " << nbins << endl;
  
  TH1D *hrate=new TH1D("hrate","hrate",nbins,tstart,tend);										 
  TH1D *hrateclean=new TH1D("hrateclean","hrateclean",nbins,tstart,tend);

  mycanv->cd(1);	 
  pdata->Draw("pcap_time_since_start>>hrate");
  hrate->Draw();
  cout << "maxrate " << maxrate << endl;

    
  for (int i=0;i<pdata->GetEntries();i++)
    {
      pdata->GetEntry(i);
      double t=pcap_time_since_start;
      //int this_bin=hrate->GetBinWithContent(t);
      int this_bin=1+(int)(t/bin_duration);
      
      // write events to new tree only if they are not in a time bin with a high event rate
      if (hrate->GetBinContent(this_bin)>maxrate) continue;

      // this also checks and removes events in adjacent time bins      
      if (this_bin>1) {
	if (hrate->GetBinContent(this_bin-1)>maxrate) continue;
      }
      if (this_bin<nbins) {
	if (hrate->GetBinContent(this_bin+1)>maxrate) continue;
      }
      
      hrateclean->Fill(t);
      pdataclean->Fill();    
    }
  mycanv->cd(2);
  hrateclean->Draw();
  mycanv->Modified();
  mycanv->Update();
  gSystem->ProcessEvents(); // Handle window interactions
  gSystem->Sleep(200);      // 100ms delay between frames	
  char a[200];
  bool wait=0;
  if (wait) (read(1,a,200));
      
  pdataclean->AutoSave();
  hrate->Reset();
  hrateclean->Reset();
}

void ProcessFiles(const char* listFileName = "filelist.dat") {
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
	  clean_ratespikes(line.c_str(),15);
        }
    }

    inputFile.close();
}
