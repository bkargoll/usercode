
//======================================================================
/**
 * @file
 *
 * @author      RWTH Aachen IIIb Top
 *
 * @version     0.0
 *
 * @brief
 *              Common plot macro functions.
 *
 */

//======================================================================
/**
 * setup ROOT's gStyle
 */

#include "TROOT.h"
#include "TStyle.h"

void
InitgStyle()
{
// {{{
    // reset everything
    gROOT->Reset();

    // canvas
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasDefH(600);
    gStyle->SetCanvasDefW(800);
    gStyle->SetCanvasDefX(0);
    gStyle->SetCanvasDefY(0);

    // pad
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);

    // frame
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetFrameLineColor(1);
    gStyle->SetFrameLineStyle(1);
    gStyle->SetFrameLineWidth(1);

    // fit and fit-function
    gStyle->SetOptFit(1111);
    gStyle->SetFitFormat("5.4g");

    // statistics box
    //gStyle->SetOptStat(1001100);
    gStyle->SetOptStat("nmri");
    gStyle->SetStatColor(0);
    gStyle->SetStatFont(42);
    gStyle->SetStatFontSize(0.025);
    gStyle->SetStatTextColor(1);
    gStyle->SetStatFormat("6.4g");
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatX(0.95);
    gStyle->SetStatY(0.95);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.15);

    // margins
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetPadRightMargin(0.04);

    // global title
    gStyle->SetOptTitle(0);

    // axis titles
    gStyle->SetTitleColor(1, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetTitleSize(0.06, "XYZ");
    gStyle->SetTitleXOffset(1.1);
    gStyle->SetTitleYOffset(1.3);

    // axis labels
    gStyle->SetLabelColor(1, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLabelOffset(0.007, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");

    // axis
    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(kFALSE);
    gStyle->SetTickLength(0.03, "XYZ");
    gStyle->SetNdivisions(510, "XYZ");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // markers and error markers
    gStyle->SetErrorX(0);
    gStyle->SetEndErrorSize(2);
//    gStyle->SetErrorMarker(20); alte root Version
    gStyle->SetMarkerStyle(20);

    // force this style on loaded objects, esp. histograms
    gROOT->ForceStyle();
};



//======================================================================
/**
 * FIXME: CINT templating not working for simultaneous use of TH1F and TH2F
 *
 * Loads a histogram from a ROOT-file, scales it, and either returns it as \c
 * returnedHistogram or adds it to \c returnedHistogram . The first if
 * \c returnedHistogram is \c NULL otherwise the latter.
 *
 * @param histogramName
 *                              name of the histogram.
 * @param inputFilename
 *                              name of the ROOT-file.
 * @param scale
 *                              scale factor for histogram.
 * @param returnedHistogram
 *                              the returned histogram (input, pass-by-reference does not work).
 * @param debug
 *                              switch on debug output, defaults to \c false
 *
 * @return
 *                              the returned histogram (output, pass-by-reference does not work).
 */
// template<class histogramTyp>
// void
// LoadHistogram(const TString& histogramName, const TString& inputFilename, double scale, histogramTyp*& returnedHistogram, bool debug = false)
// {
// // {{{
//     TFile inputFile(inputFilename);
//     if (!inputFile.IsOpen())
//     {
//         cerr << "Could not open '" << inputFilename << "' for reading." << endl;
//     }
//     else
//     {
//         histogramTyp* histogram = dynamic_cast<histogramTyp*>( inputFile.Get(histogramName) );
//         if (!histogram)
//         {
//             cerr << "No histogram named '" << histogramName << "' in file '" << inputFilename << "'" << endl;
//         }
//         else
//         {
//             if (debug) cerr << inputFilename << " " << histogramName << " entries=" << histogram->GetEntries() << " integral=" << histogram->Integral() << " scale*integral=" << scale*histogram->Integral() << endl;
//             histogram->Scale(scale);
//             if (!returnedHistogram)
//             {
//                 returnedHistogram = new histogramTyp(*histogram);
//                 returnedHistogram->SetDirectory(0); // otherwise "TFile inputFile" owns this returnedHistogram and deletes it when "TFile inputFile" goes out of scope
//             }
//             else
//             {
//                 returnedHistogram->Add(histogram);
//             }
//         }
//         inputFile.Close();
//     }
// // }}}
// };



// //======================================================================
// /**
//  * Loads a histogram from a ROOT-file, scales it, and either returns it as \c
//  * returnedHistogram or adds it to \c returnedHistogram . The first if
//  * \c returnedHistogram is \c NULL otherwise the latter.
//  *
//  * @param histogramName
//  *                              name of the histogram.
//  * @param inputFilename
//  *                              name of the ROOT-file.
//  * @param scale
//  *                              scale factor for histogram.
//  * @param returnedHistogram
//  *                              the returned histogram (input, pass-by-reference does not work).
//  * @param debug
//  *                              switch on debug output, defaults to \c false
//  *
//  * @return
//  *                              the returned histogram (output, pass-by-reference does not work).
//  */
// void
// LoadHistogramTH1F(const TString& histogramName, const TString& inputFilename, double scale, TH1F*& returnedHistogram, bool debug = false)
// {
// // {{{
//     TFile inputFile(inputFilename);
//     if (!inputFile.IsOpen())
//     {
//         cerr << "Could not open '" << inputFilename << "' for reading." << endl;
//     }
//     else
//     {
//         TH1F* histogram = dynamic_cast<TH1F*>( inputFile.Get(histogramName) );
//         if (!histogram)
//         {
//             cerr << "No histogram named '" << histogramName << "' in file '" << inputFilename << "'" << endl;
//         }
//         else
//         {
//             if (debug) cerr << inputFilename << " " << histogramName << " entries=" << histogram->GetEntries() << " integral=" << histogram->Integral() << " scale*integral=" << scale*histogram->Integral() << endl;
//             histogram->Scale(scale);
//             if (!returnedHistogram)
//             {
//                 returnedHistogram = new TH1F(*histogram);
//                 returnedHistogram->SetDirectory(0); // otherwise "TFile inputFile" owns this returnedHistogram and deletes it when "TFile inputFile" goes out of scope
//             }
//             else
//             {
//                 returnedHistogram->Add(histogram);
//             }
//         }
//         inputFile.Close();
//     }
// // }}}
// };



// //======================================================================
// /**
//  * Loads a histogram from a ROOT-file, scales it, and either returns it as \c
//  * returnedHistogram or adds it to \c returnedHistogram . The first if
//  * \c returnedHistogram is \c NULL otherwise the latter.
//  *
//  * @param histogramName
//  *                              name of the histogram.
//  * @param inputFilename
//  *                              name of the ROOT-file.
//  * @param scale
//  *                              scale factor for histogram.
//  * @param returnedHistogram
//  *                              the returned histogram (input, pass-by-reference does not work).
//  * @param debug
//  *                              switch on debug output, defaults to \c false
//  *
//  * @return
//  *                              the returned histogram (output, pass-by-reference does not work).
//  */
// void
// LoadHistogramTH2F(const TString& histogramName, const TString& inputFilename, double scale, TH2F*& returnedHistogram, bool debug = false)
// {
// // {{{
//     TFile inputFile(inputFilename);
//     if (!inputFile.IsOpen())
//     {
//         cerr << "Could not open '" << inputFilename << "' for reading." << endl;
//     }
//     else
//     {
//         TH2F* histogram = dynamic_cast<TH2F*>( inputFile.Get(histogramName) );
//         if (!histogram)
//         {
//             cerr << "No histogram named '" << histogramName << "' in file '" << inputFilename << "'" << endl;
//         }
//         else
//         {
//             if (debug) cerr << inputFilename << " " << histogramName << " entries=" << histogram->GetEntries() << " integral=" << histogram->Integral() << " scale*integral=" << scale*histogram->Integral() << endl;
//             histogram->Scale(scale);
//             if (!returnedHistogram)
//             {
//                 returnedHistogram = new TH2F(*histogram);
//                 returnedHistogram->SetDirectory(0); // otherwise "TFile inputFile" owns this returnedHistogram and deletes it when "TFile inputFile" goes out of scope
//             }
//             else
//             {
//                 returnedHistogram->Add(histogram);
//             }
//         }
//         inputFile.Close();
//     }
// // }}}
// };



// //======================================================================
// /**
//  * Returns the number of events in a file as defined by convention through the
//  * content of the histogram "numberOfEvents".
//  *
//  * @param inputFilename
//  *                      name of the input root file.
//  * @return
//  *                      the number of events in histogram "numberOfEvents"
//  */
// int
// numberOfEvents(const TString& inputFilename)
// {
// // {{{
//     TH1F *h = 0;
//     LoadHistogramTH1F("numberOfEvents", inputFilename, 1.0, h);
//     int nEvents = 1;
//     if (h)
//     {
//         nEvents = h->GetEntries();
//         delete h;
//     }
//     return nEvents;
// // }}}
// };



// //======================================================================
// /**
//  * Returns a suitable epsFilename build from a file or directory name and
//  * histogram name.
//  *
//  * @param fileOrDirectoryName
//  *                            name of the file or directory
//  * @param histogramNamw
//  *                            name of the histogram
//  * @return
//  *                            a suitable epsFilename
//  */
// TString
// epsFilename(const TString& fileOrDirectoryName, const TString& histogramName)
// {
// // {{{
//     TString s = fileOrDirectoryName;
//     if (s.EndsWith("/")) s.Remove(s.Length()-1);
//     s = gSystem->BaseName(s);
//     s.ReplaceAll(".root", "");
//     s.Append(".").Append(histogramName);
//     s.ReplaceAll(".", "_");
//     return s;
// // }}}
// };



// //
// // Local Variables:
// // End:
// //


