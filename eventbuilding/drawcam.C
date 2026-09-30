#include <TCanvas.h>
#include <TBox.h>
#include <TLine.h>
#include <TText.h>
#include <TLatex.h>
#include <TColor.h>

void draw_panoseti_camera() {
    // 1. Create the main layout canvas
    TCanvas *c1 = new TCanvas("c1", "PANOSETI SiPM Camera Layout", 900, 900);
    c1->Range(-5, -5, 37, 37); // Buffer margins surrounding the 32x32 pixel area
    c1->SetFillColor(kWhite);

    // Color definitions for visual hierarchy
    Int_t pixelColor   = TColor::GetColor("#2c3e50"); // Dark slate for individual pixels
    Int_t arrayColor   = TColor::GetColor("#16a085"); // Teal border for 8x8 arrays
    Int_t quadColor    = TColor::GetColor("#c0392b"); // Solid dark red for quadrants
    Int_t segmentColor = TColor::GetColor("#7f8c8d"); // Gray lines separating pixels

    // 2. Draw 32x32 individual pixels
    for (int x = 0; x < 32; x++) {
        for (int y = 0; y < 32; y++) {
            TBox *pixel = new TBox(x, y, x + 0.9, y + 0.9);
            pixel->SetFillColor(pixelColor);
            pixel->SetFillStyle(3001); // Soft solid fill style
            pixel->SetLineColor(segmentColor);
            pixel->SetLineWidth(1);
            pixel->Draw();
        }
    }

    // 3. Overlay the 16 separate Hamamatsu 8x8 SiPM matrices
    for (int arrX = 0; arrX < 4; arrX++) {
        for (int arrY = 0; arrY < 4; arrY++) {
            // Calculate boundaries for each 8x8 tile
            double xLow  = arrX * 8 - 0.1;
            double yLow  = arrY * 8 - 0.1;
            double xHigh = (arrX + 1) * 8;
            double yHigh = (arrY + 1) * 8;

            TBox *arrayBorder = new TBox(xLow, yLow, xHigh, yHigh);
            arrayBorder->SetFillStyle(0); // Transparent interior
            arrayBorder->SetLineColor(arrayColor);
            arrayBorder->SetLineWidth(3);
            arrayBorder->Draw();

            // Label each individual 8x8 sub-array matrix
            TText *tArr = new TText(arrX * 8 + 4, arrY * 8 + 4, Form("8x8\nArray %d", (arrY * 4) + arrX + 1));
            tArr->SetTextAlign(22); // Center alignment
            tArr->SetTextColor(kWhite);
            tArr->SetTextSize(0.018);
            tArr->Draw();
        }
    }

    // 4. Trace the distinct 4 Quadrant Boards (Split across x=16 and y=16)
    TLine *vQuadSplit = new TLine(16, -0.5, 16, 32.5);
    vQuadSplit->SetLineColor(quadColor);
    vQuadSplit->SetLineWidth(6);
    vQuadSplit->Draw();

    TLine *hQuadSplit = new TLine(-0.5, 16, 32.5, 16);
    hQuadSplit->SetLineColor(quadColor);
    hQuadSplit->SetLineWidth(6);
    hQuadSplit->Draw();

    // 5. Annotate Quadrant Hardware Identifiers
    TText *q1 = new TText(8, 24, "QUADRANT 1");  q1->SetTextAlign(22); q1->SetTextColor(kRed+2); q1->SetTextSize(0.035); q1->Draw();
    TText *q2 = new TText(24, 24, "QUADRANT 2"); q2->SetTextAlign(22); q2->SetTextColor(kRed+2); q2->SetTextSize(0.035); q2->Draw();
    TText *q3 = new TText(8, 8, "QUADRANT 3");   q3->SetTextAlign(22); q3->SetTextColor(kRed+2); q3->SetTextSize(0.035); q3->Draw();
    TText *q4 = new TText(24, 8, "QUADRANT 4");  q4->SetTextAlign(22); q4->SetTextColor(kRed+2); q4->SetTextSize(0.035); q4->Draw();

    // 6. Descriptive metadata header blocks
    TLatex *title = new TLatex(16, 35, "PANOSETI SiPM Camera Focal Plane Schematic");
    title->SetTextAlign(22); title->SetTextSize(0.04); title->SetTextFont(62);
    title->Draw();

    TLatex *sub = new TLatex(16, 33.5, "Total Scale: 32 #times 32 Pixels (1,024 Channels) | Field of View: 9.9^{#circ} #times 9.9^{#circ}");
    sub->SetTextAlign(22); sub->SetTextSize(0.022); sub->SetTextFont(52);
    sub->Draw();

    // Technical specifications legend text box
    TLatex *leg = new TLatex(-3, -3, "#bf{Legend Layout:}   #color[7934]{Pixel Matrix (1,024 Total)}   |   #color[1424]{16 #times Hamamatsu S13361 (8#times8 Arrays)}   |   #color[1179]{4 #times Hardware Quadrant Boards}");
    leg->SetTextAlign(12); leg->SetTextSize(0.02);
    leg->Draw();

    // Render configuration pipeline update
    c1->Update();
    c1->SaveAs("panoseti_sipm_camera_schematic.png");
}
