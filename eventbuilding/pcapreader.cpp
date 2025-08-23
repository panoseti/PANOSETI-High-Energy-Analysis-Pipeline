/* 
written by Jamie Holder
Use pcap libraries to read PANOSETI data and convert to ROOT format.
Should compile with: g++ -Wall -o pcapreader pcapreader.cpp -lpcap `root-config --cflags --glibs`

This works for root v6-30-06

For root v6.22  I had to remove -lfreetype from the root libraries returned by root-config so I compile like so instead:

g++ -Wall -o pcapreader pcapreader.cpp -lpcap -stdlib=libc++ -D_REENTRANT -std=c++11 -m64 -I/Applications/root_v6.22.02/include -L/Applications/root_v6.22.02/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame  -stdlib=libc++ -lpthread -lm -ldl

Each panoseti pcap packet consists of:
* a pcap header which holds the packet length and a timestamp (https://www.winpcap.org/docs/docs_412/html/structpcap__pkthdr.html
* a 42-byte ethernet/UDP header
* the PANOSETI data payload described here: 
* https://github.com/panoseti/panoseti/wiki/Quabo-packet-interface#science-packets
* Wireshark is a very useful tool for examining the raw bytes in a pcap file. https://www.wireshark.org/
*/

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <pcap.h>
#include <TTree.h>
#include <TFile.h>
#include "pcapreader.h"
using namespace std;

#define IFSZ 16
#define FLTRSZ 120
#define MAXHOSTSZ 256
#define ETHER_ADDR_LEN	6

int packets = 0;   /* running count of packets read in */
long start_secs = 0; /* start time (secs) in pcap */
int start_usecs = 0; /* start time (usec) in pcap */

TTree* pdata;

//struct panodata pdata;


int usage(char *progname)
{
        printf("Usage: pcapreader  [<file name>]\n");
	
        return(0);
}


void print_panodata(struct panodata i_panodata, int npix)
{
  printf("Acquisition mode: %X\n",i_panodata.acq_mode);
  printf("Packer version  : %X\n",i_panodata.packer_ver);  
  printf("Packet Number   : %u\n",i_panodata.packet_no);
  printf("Board Location  : %u\n",i_panodata.boardloc);
  printf("TAI WR time     : %u\n",i_panodata.TAI);
  printf("nanosecs        : %u\n",i_panodata.nanosec);
  printf("unused          : %u\n",i_panodata.dummy);
  for (int i=0;i<npix;i++)
    {
      printf("pix_data unsigned %d : %u\n",i,i_panodata.pix_data_unsigned[i]);
      printf("pix_data signed %d : %u\n",i,i_panodata.pix_data_signed[i]);
    }
}

void parse_data(u_char *user, const struct pcap_pkthdr *hdr, const u_char *packet)
{
  //cout << "This is the user variable, in case it's needed: " << user << endl;
  const struct panodata_bytes *read_panodata;
  struct panodata store_panodata;
  int npix=256;
  double pcap_start_time=0;
  double pcap_time=0;
  double pcap_time_since_start=0;
  double WR_time=0;

  
  //printf("packet no.: %u caplen: %u \n",packets, hdr->caplen);
  if (packets==0)
    {
      start_secs=hdr->ts.tv_sec;
      start_usecs=hdr->ts.tv_usec;
      //initialize the root Tree;
      pdata = new TTree("pdata","Panoseti data");
      pdata->Branch("npix",&npix,"npix/I");
      pdata->Branch("pcap_time",&pcap_time,"pcap_time/D");
      pdata->Branch("pcap_time_since_start",&pcap_time_since_start,"pcap_time_since_start/D");
      pdata->Branch("acq_mode",&store_panodata.acq_mode,"acq_mode/b");
      pdata->Branch("packer_ver",&store_panodata.packer_ver,"packer_ver/b");
      pdata->Branch("packet_no",&store_panodata.packet_no,"packet_no/s");
      pdata->Branch("boardloc",&store_panodata.boardloc,"boardloc/s");
      pdata->Branch("TAI",&store_panodata.TAI,"TAI/i");
      pdata->Branch("nanosec",&store_panodata.nanosec,"nanosec/i");
      pdata->Branch("WR_time",&WR_time,"WR_time/D");
      pdata->Branch("dummy",&store_panodata.dummy,"dummy/s");
      pdata->Branch("pix_data_unsigned",&store_panodata.pix_data_unsigned,"pix_data_unsigned[npix]/s");
      pdata->Branch("pix_data_signed",&store_panodata.pix_data_signed,"pix_data_signed[npix]/S");

    }
  pcap_start_time=start_secs+start_usecs/1.e6;
  pcap_time=hdr->ts.tv_sec+hdr->ts.tv_usec/1.e6;
  pcap_time_since_start=pcap_time - pcap_start_time;
  //printf("pcap time since first packet: %lf \n",pcap_time_since_start);

  // Just selecting science data on packet length for the moment, which is not very robust. 
  //if (hdr->caplen!=570) printf("not a science packet\n");  
  if (hdr->caplen==570)
    {
      read_panodata = (struct panodata_bytes*)(packet+42);
      
      // JH - not sure if we really need all this memcpy stuff,
      // but it seems to unpack the bytes correctly.
      
      u_short tempshort;
      short temp_signed_short;
      u_int tempint;

      store_panodata.acq_mode=read_panodata->acq_mode;
      store_panodata.packer_ver=read_panodata->packer_ver;
    
      memcpy(&tempshort, read_panodata->packet_no, sizeof(u_short));    
      store_panodata.packet_no=tempshort;
      
      memcpy(&tempshort, read_panodata->boardloc, sizeof(u_short));    
      store_panodata.boardloc=tempshort;
    
      memcpy(&tempint, read_panodata->TAI, sizeof(u_int));    
      store_panodata.TAI=tempint;

      memcpy(&tempint, read_panodata->nanosec, sizeof(u_int));    
      store_panodata.nanosec=tempint;

      memcpy(&tempshort, read_panodata->dummy, sizeof(u_short));    
      store_panodata.dummy=tempshort;
    
      for (int i=0;i<256;i++)
	{
	  memcpy(&tempshort, read_panodata->pix_data[i], sizeof(u_short));
	  store_panodata.pix_data_unsigned[i]=tempshort;
	  
	  memcpy(&temp_signed_short, read_panodata->pix_data[i], sizeof(short));
	  store_panodata.pix_data_signed[i]=temp_signed_short;

	}    

      //print_panodata(store_panodata, 0);
      //pdata=store_panodata;
      WR_time=store_panodata.TAI+store_panodata.nanosec/1.e9;

      pdata->Fill();
    }
  packets++; /* keep a running total of number of packets read in */
}

  
int main(int argc, char **argv)
{
        pcap_t *p;               /* packet capture descriptor */
        char filename[200];       /* name of file to read packet data from */
        char errbuf[PCAP_ERRBUF_SIZE];  /* buffer to hold error text */
        char prestr[80];         /* prefix string for errors from pcap_perror */

	if (argc != 2)
	  {
	    usage(argv[0]);
	    return(0);
	  }
        if (argc == 2)
                strcpy(filename,argv[1]);

	printf("filename %s \n",filename);

        /*
         * Open a file containing packet capture data. 
         */
        if (!(p = pcap_open_offline(filename, errbuf))) {
                fprintf(stderr,
                        "Error in opening file, %s, for reading: %s\n",
                        filename, errbuf);
                return(2);
        }

	char outfile[200];
	snprintf(outfile,200,"%s.root",filename);
	printf("%s",outfile);

        /*
         * Call pcap_dispatch() with a count of 0 which will cause
         * pcap_dispatch() to read and process packets until an error or EOF
         * occurs. For each packet read from the file, the output routine,
         * parse_data(), will be called to print the data.
	 * JH: This is a little annoying because it means that everything has to 
	 * happen inside of the output routine (parse_data). Not sure how to avoid that. 
         */
        if (pcap_dispatch(p, 0, &parse_data, (u_char *)"test") < 0) {
                /*
                 * Otherwise print out appropriate text, followed by the error message
                 * generated by the packet capture library.
                 */
	  snprintf(prestr,80,"Error reading packets:");
                pcap_perror(p,prestr);
                return(4);
        }
        printf("\nNumber of Raw packets read in: %d\n", packets);

        pcap_close(p);
	cout << "Number of PANOSETI QUABO Board entries in ROOT output: " << pdata->GetEntries() << endl;
	TFile *f = new TFile(outfile,"recreate");
	pdata->Write();
	f->Close();
}
