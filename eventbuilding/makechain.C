TChain ch;

void MakeChainFromList(const char* listFileName = "filelist.dat") {
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
	// use an empty line as a break point
        if (line.empty()) break;
        if (!line.empty() && line[0] != '#') { // Ignore empty lines and comments
	  char treename[200];
	  snprintf(treename,200,"%s/camdata",line.c_str());
	  ch.Add(treename);
        }
    }

    inputFile.close();
}

void make_chain(const char *infile = "filelist.dat", const char *outfile = "all.root")
{

  MakeChainFromList(infile);
  
  double pt=0;
  ch.SetBranchAddress("cam_pcap_time",&pt);

  ch.GetEntry(1);
  double mintime=pt;
  ch.GetEntry(ch.GetEntries()-1);
  double maxtime=pt;
  
  cout << mintime << " " << maxtime << " " << (maxtime-mintime)/60/60. << endl;
  double tstart=-500;
  double tend=maxtime-mintime+500;

  TCanvas *mycanv=new TCanvas("mycanv","mycanv");
  mycanv->Draw();
  mycanv->cd();
  TH1D *hrate=new TH1D("hrate","Event rate",(tend-tstart)/100.,tstart,tend);
 hrate->SetLineWidth(2);
 hrate->SetStats(0);
 //hrate->SetMaximum(80);
 hrate->SetXTitle("Time(s)");
 hrate->SetYTitle("Events per 100 secs");

 for (int i=0;i<ch.GetEntries();i++)
   {
     ch.GetEntry(i);
     hrate->Fill(pt-mintime);
   }

   // ch.Draw("cam_pcap_time-tstart>>hrate");
 hrate->Draw();
 
 TFile *outputFile = TFile::Open(outfile, "RECREATE");
 TTree *newTree = ch.CloneTree(-1, "fast");

 outputFile->Write();
 outputFile->Close();
 ch.Reset();
 delete outputFile;
 
}

