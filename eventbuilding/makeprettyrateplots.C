{
  static const int ntel=3;
  TFile *_file0 = TFile::Open("all.root.corr.array");
  TTree *arraydata=(TTree*)_file0->Get("arraydata");
  double array_pcap_time[ntel];
  arraydata->SetBranchAddress("array_pcap_time",array_pcap_time);

  arraydata->GetEntry(0);
  double mintime=9e99;
  for (int i=0;i<ntel;i++)
    {
      if (array_pcap_time[i] >-999 && array_pcap_time[i]<mintime) mintime=array_pcap_time[i];
    }

  arraydata->GetEntry(arraydata->GetEntries()-1);
  double maxtime=0;
  for (int i=0;i<ntel;i++)
    {
      if (array_pcap_time[i] >-999 && array_pcap_time[i]>maxtime) maxtime=array_pcap_time[i];
    }

  cout << mintime << " " << maxtime << " " << maxtime-mintime << endl;
  double tstart=mintime-500;
  double tend=maxtime+500;
  int nbins100=ceil((tend-tstart)/100);
  int nbins500=ceil((tend-tstart)/500);

  char scope_name[ntel][100];
  snprintf(scope_name[0],100,"PTI");
  snprintf(scope_name[1],100,"Fern");
  snprintf(scope_name[2],100,"Winter");
  int scope_id[ntel];
  scope_id[0]=1000; // PTI
  scope_id[1]=1008; // Fern
  scope_id[2]=1012; // Winter

  int ymax=100;

  TCanvas *singles=new TCanvas("singles","singles",1350,350);
  singles->Divide(ntel,1);
  
  TH1D *hrate[ntel];
  for (int i=0;i<ntel;i++)
    {
      char hname[200];
      snprintf(hname,200,"hrate%d",i);	
      char htitle[200];
      snprintf(htitle,200,"Rate for %s (T%d)",scope_name[i], scope_id[i]);	
      hrate[i]=new TH1D(hname,htitle,nbins100,tstart,tend);
      hrate[i]->SetLineWidth(3);
      hrate[i]->SetStats(0);
      hrate[i]->SetXTitle("Time(s)");
      hrate[i]->SetYTitle("Events/100s");
    }
for (int i=0;i<ntel;i++)
    {
      singles->cd(i+1);
      char var[200];
      snprintf(var,200,"array_pcap_time[%d]>>hrate%d",i,i);
      char cut[200];
      snprintf(cut,200,"array_pcap_time[%d]>0",i);
      arraydata->Draw(var,cut);
      hrate[i]->SetMaximum(ymax);
      hrate[i]->Draw();
    }

  TCanvas *esingles=new TCanvas("esingles","exclusive singles",1350,350);
  esingles->Divide(ntel,1);
  
  TH1D *ehrate[ntel];
  for (int i=0;i<ntel;i++)
    {
      char hname[200];
      snprintf(hname,200,"ehrate%d",i);	
      char htitle[200];
      snprintf(htitle,200,"Rate for %s (T%d)",scope_name[i], scope_id[i]);	
      ehrate[i]=new TH1D(hname,htitle,nbins100,tstart,tend);
      ehrate[i]->SetLineWidth(3);
      ehrate[i]->SetStats(0);
      ehrate[i]->SetXTitle("Time(s)");
      ehrate[i]->SetYTitle("Events/100s");
    }
for (int i=0;i<ntel;i++)
    {
      esingles->cd(i+1);
      char var[200];
      snprintf(var,200,"array_pcap_time[%d]>>ehrate%d",i,i);
      char cut[200];
      snprintf(cut,200,"array_pcap_time[%d]>0 && array_nteltrig==1",i);
      arraydata->Draw(var,cut);
      ehrate[i]->SetMaximum(ymax);
      ehrate[i]->Draw();
    }


  TCanvas *doubles=new TCanvas("doubles","doubles",1350,350);
  doubles->Divide(ntel,1);
  
  TH1D *dhrate01;
  TH1D *dhrate02;
  TH1D *dhrate12;
  TH1D *hrate012;

  dhrate01=new TH1D("dhrate01","Rate for PTI+Fern",nbins100,tstart,tend);
  dhrate02=new TH1D("dhrate02","Rate for PTI+Winter",nbins100,tstart,tend);
  dhrate12=new TH1D("dhrate12","Rate for Fern+Winter",nbins100,tstart,tend);
  hrate012=new TH1D("hrate012","Rate for Fern+Winter+PTI",nbins100,tstart,tend);

  dhrate01->SetLineWidth(3);
  dhrate01->SetLineWidth(3);
  dhrate01->SetStats(0);
  dhrate01->SetXTitle("Time(s)");
  dhrate01->SetYTitle("Events/400s");

  dhrate02->SetLineWidth(3);
  dhrate02->SetStats(0);
  dhrate02->SetXTitle("Time(s)");
  dhrate02->SetYTitle("Events/400s");

  dhrate12->SetLineWidth(3);
  dhrate12->SetStats(0);
  dhrate12->SetXTitle("Time(s)");
  dhrate12->SetYTitle("Events/400s");

  hrate012->SetLineWidth(3);
  hrate012->SetStats(0);
  hrate012->SetXTitle("Time(s)");
  hrate012->SetYTitle("Events/400s");

  ymax=60;
  doubles->cd(1);
  arraydata->Draw("array_pcap_time[0]>>dhrate01","array_pcap_time[0]>0 && array_pcap_time[1]>0 && array_nteltrig==2");
  dhrate01->SetMaximum(ymax);
  dhrate01->Rebin(4);
  doubles->cd(2);
  arraydata->Draw("array_pcap_time[0]>>dhrate02","array_pcap_time[0]>0 && array_pcap_time[2]>0 && array_nteltrig==2");
  dhrate02->SetMaximum(ymax);
  dhrate02->Rebin(4);
  doubles->cd(3);
  arraydata->Draw("array_pcap_time[1]>>dhrate12","array_pcap_time[1]>0 && array_pcap_time[2]>0 && array_nteltrig==2");
  dhrate12->SetMaximum(ymax);
  dhrate12->Rebin(4);
  
  TCanvas *triples=new TCanvas("triples","triples",450,350);
  triples->cd();
  arraydata->Draw("array_pcap_time[1]>>hrate012","array_pcap_time[0]>0 && array_pcap_time[1]>0 && array_pcap_time[2]>0 && array_nteltrig==3");
  hrate012->SetMaximum(ymax);
  hrate012->Rebin(4);
  
 
}
