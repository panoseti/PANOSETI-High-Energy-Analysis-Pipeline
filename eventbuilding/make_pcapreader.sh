g++ -Wall -o pcapreader pcapreader.cpp -lpcap `root-config --cflags --glibs`

# For root v6.22  I had to remove -lfreetype from the root libraries returned by root-config so I compile like so instead:
#g++ -Wall -o pcapreader pcapreader.cpp -lpcap -stdlib=libc++ -D_REENTRANT -std=c++11 -m64 -I/Applications/root_v6.22.02/include -L/Applications/root_v6.22.02/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame  -stdlib=libc++ -lpthread -lm -ldl
