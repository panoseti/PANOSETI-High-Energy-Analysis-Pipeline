// Simple runner for panodisplay to view example CORSIKA simulation
// Usage: root -l run_example.C

void run_example(int eventNumber = 0) {
    std::cout << "Reading example.root" << std::endl;
    std::cout << "Displaying event " << eventNumber << std::endl;
    gROOT->LoadMacro("panodisplay.C");

    // Read the file and display event - check both possible locations
    const char* rootFile = "example.root";
    if (gSystem->AccessPathName(rootFile)) {
        // Check if it's in ../example/
        std::string altPath = "../example/example.root";
        if (!gSystem->AccessPathName(altPath.c_str())) {
            rootFile = altPath.c_str();
        }
    }

    readFile(rootFile);

    // Set visibility to true to open display window
    gROOT->SetBatch(kFALSE);

    panodisplay(eventNumber);
}