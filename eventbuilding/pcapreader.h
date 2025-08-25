struct panodata_bytes {
  // This holds the panoseti data as individual bytes and arrays of bytes
  u_char acq_mode;
  u_char packer_ver;
  u_char packet_no[2];
  u_char boardloc[2];
  u_char TAI[4];
  u_char nanosec[4];
  u_char dummy[2];
  u_char pix_data[255][2];
};

struct panodata {
  // This holds the panoseti data as data types of the correct size
  u_char acq_mode;
  u_char packer_ver;
  u_short packet_no;
  u_short boardloc;
  u_int TAI;
  u_int nanosec;
  u_short dummy;
  u_short pix_data_unsigned[255];
  short pix_data_signed[255];
};
