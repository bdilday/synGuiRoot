#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMinuit.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TCanvas.h"

#include "TGComboBox.h"
#include "TGButton.h"
#include "TGButtonGroup.h"
#include "TRootEmbeddedCanvas.h"
#include "TGLayout.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGTextEntry.h"
#include "TGSlider.h"
#include "TGLabel.h"

#include <stdlib.h>
#include <stdio.h>

#define MAXWAVELENGTHS 50000
#define MAXION 100
#define FSLIDERID 100
#define FFRAMEID  200
#define FTEXTID   300
#define FCHECKID  400
#define FBUFFERID 500
#define FXTEXTMIN 600
#define FXTEXTMAX 601
#define FRADIOTPHOTID 700
#define FRADIOETEMPID 800
#define FCOMBOBOXID 900
#define NETEMP 3
#define NTPHOT 3
#define NLOGTAU 8

static const double SLIDER_WEIGHT_SCALE = 0.0001;
static const double SLIDER_Z1_SCALE = 0.0001;

static const int IBB = 0;
static const int IZ1 = (MAXION-1);
static const int ISNSCALE = (MAXION-2);


enum ETestCommandIdentifiers {
   HId1,
   HId2,
   HId3,
   HCId1,
   HCId2,

   HSId1
};

class TSynGui : public TGMainFrame {
private:
  ;
public:
  TRootEmbeddedCanvas *fCanvas;
  TGLayoutHints       *fLcan;
  TF1                 *fFitFcn;
  TGHorizontalFrame   *fHframe0, *fHframe1, *fHframe2;
  TGLayoutHints       *fBly, *fBfly1, *fBfly2, *fBfly3;
  TGHSlider     *fHslider1;
  TGTextEntry         *fTeh1, *fTeh2, *fTeh3;
  TGTextBuffer        *fTbh1, *fTbh2, *fTbh3;
  TGCheckButton       *fCheck1, *fCheck2;

  TGraph *errorGraphUpper, *errorGraphLower;
  TGraph *snGraph;
  TGraph *modelGraph;
  TGraph *ionGraphs[MAXION];

  char activeIonFile[1024];
  char refsFile[1024];
  char snFile[1024];
  
  TGHSlider *fSliders[MAXION];
  TGHorizontalFrame *fSliderFrames[MAXION], *fSliderFrameRight, *fSliderFrameRadio1;
  TGLayoutHints     *fSliderBlys[MAXION];
  TGLayoutHints     *fSliderBflys[MAXION];
  TGTextEntry *fSliderTextEntries[MAXION];
  TGCheckButton *fSliderCheckButtons[MAXION];
  TGTextBuffer        *fSliderTextBuffers[MAXION];
  TGLabel *fLabelChi2;

  TGButtonGroup *fRadioGroupTphot;
  TGButtonGroup *fRadioGroupEtemp;

  TGComboBox *fComboBox;

  TGRadioButton *fRadioTphot[NTPHOT];
  TGRadioButton *fRadioEtemp[NETEMP];

  double tphotValues[NTPHOT];
  double etempValues[NETEMP];
  double logtauValues[NLOGTAU];

  double currentTphotValue;
  double currentEtempValue;
  double currentLogtauValue;

  TGTextEntry *fXminText;
  TGTextEntry *fXmaxText;
  TGTextBuffer *fXminTextBuffer;
  TGTextBuffer *fXmaxTextBuffer;

  double fSlidersCurrentPos[MAXION];

  int NION;
  int ionList[MAXION];
  char ionNames[MAXION][32];
  int ionIsActive[MAXION];

  int NSNWAVELENGTHS;
  double snWavelengths[MAXWAVELENGTHS];
  double snFlux[MAXWAVELENGTHS];
  double snFluxError[MAXWAVELENGTHS];

  int NIONWAVELENGTHS;
  double ionWavelengths[MAXION][MAXWAVELENGTHS];
  double ionFlux[MAXION][MAXWAVELENGTHS];
  double ionFluxError[MAXION][MAXWAVELENGTHS];

  double modelWavelengths[MAXWAVELENGTHS];
  double modelFlux[MAXWAVELENGTHS];

  double fChi2;

  double t_phot;
  double etemp;
  double log_tau;

/**********************/

  TSynGui();
  virtual ~TSynGui();

  void LinearBBFit();
  void ClearAll();
  void GetChi2();
  void DoDraw();
  void ReadSNDATA();
  void ReadAPPDATA();  
  void CloseWindow();
  void DoText(const char *text);
  void DoRadio(const Int_t iradio);
  void DoLogtauCombo();
  void DoSlider();
  void HandleButtons();
  void SetActiveIonFile(const char* s);
  void SetSNFile(const char* s);
  void ReadActiveIonFile();
  void SetRefsFile(const char* s);
  void ScaleToMean(int n, double *a, double target);
  double GetMean(int n, double *a);
  void MultiplyArray(int n, double *oldA, double *newA, double multiplier);
  void DoXFitText(const char *text);

  void GetModelFlux(int nW, double *a);
  void DoFit();
  void SaveFigure();
  void SaveFile();

  double minXForFit;
  double maxXForFit;

  ClassDef(TSynGui, 0)
};

/** globals **/
TSynGui *Gsyngui;
void GGetChi2();
void Gmy_mnf_lk(int &npar, double *gin, double &f, double *x, int iflag);
double GfChi2;
double currentWeightVals[MAXION];

//______________________________________________________________________________
void TSynGui::SaveFigure() {
  fCanvas->GetCanvas()->Print("blah.ps");
  /**
  TCanvas *tc;
  tc = ((TObject *) fCanvas->GetCanvas())->Copy(tc);
  tc->SetCanvasSize(10, 8);
  tc->Print("blah.ps");
  **/
}

//______________________________________________________________________________
void TSynGui::SaveFile() {
  TFile *tf;
  tf = new TFile("blah.root", "RECREATE");
  Gsyngui->Write("myGsyngui");
  tf->Close();

}


//______________________________________________________________________________
void TSynGui::LinearBBFit() {
  char buf[32];
  double xx[MAXWAVELENGTHS];
  double yy[MAXWAVELENGTHS];
  double yye[MAXWAVELENGTHS];
  double ans, g2, gy, s2, t8, z1, target;

  t8 = SLIDER_WEIGHT_SCALE*Gsyngui->fSliders[ISNSCALE]->GetPosition();
  z1 = 1.0+SLIDER_Z1_SCALE*Gsyngui->fSliders[IZ1]->GetPosition();
  
  Gsyngui->MultiplyArray(Gsyngui->NSNWAVELENGTHS, Gsyngui->snWavelengths, xx, 1.0/z1);
  Gsyngui->MultiplyArray(Gsyngui->NSNWAVELENGTHS, Gsyngui->snFlux, yy, t8);
  Gsyngui->MultiplyArray(Gsyngui->NSNWAVELENGTHS, Gsyngui->snFluxError, yye, t8);

  //  target = currentWeightVals[IBB];
  target = 1.0;

  Gsyngui->GetModelFlux(Gsyngui->NIONWAVELENGTHS, Gsyngui->modelFlux);
  Gsyngui->ScaleToMean(Gsyngui->NIONWAVELENGTHS, Gsyngui->modelFlux, target);
  TGraph *modelGraph = new TGraph(Gsyngui->NIONWAVELENGTHS, 
				  Gsyngui->ionWavelengths[IBB], 
				  Gsyngui->modelFlux);  

  gy = 0.0;
  g2 = 0.0;
  double g;
  for (int iW=0;iW<Gsyngui->NSNWAVELENGTHS;iW++) {
    if(xx[iW]>=Gsyngui->minXForFit && xx[iW]<=Gsyngui->maxXForFit) {    
      s2 = yye[iW]*yye[iW];
      g = modelGraph->Eval(xx[iW]);
      g2 += g*g/s2;
      gy += g*yy[iW]/s2;
    }
  }

  ans = gy/g2;
  printf("LinearFit n= %10d gy= %.12e g2= %.12e ans= %.12e \n", NSNWAVELENGTHS, gy, g2, ans);

  fSliders[IBB]->SetPosition((int)(ans/SLIDER_WEIGHT_SCALE));
  currentWeightVals[IBB] = ans;
  sprintf(buf, "%.3f", fSliders[IBB]->GetPosition()*SLIDER_WEIGHT_SCALE);
  fSliderTextBuffers[IBB]->Clear();
  fSliderTextBuffers[IBB]->AddText(0, buf);
  gClient->NeedRedraw(fSliderTextEntries[IBB]);
  DoDraw();
}


//______________________________________________________________________________
void TSynGui::DoFit() {
  int jerr;
  int &ierr = jerr;
  int iIon; 
  double scale, t8, fittedPars[MAXION], fittedParErrors[MAXION];
  //  char parName[64];
  char buf[32], thisStr[32];
  TMinuit *tm = new TMinuit(NION);
  Double_t dumval, dumerr;
  Double_t &rval=dumval, &rerr=dumerr;

  //  tm->SetFCN(&TSynGui::mdouble ionParCurrentValues[MAXION];y_mnf_lk);
  //  tm->SetFCN(this->my_mnf_lk);
  tm->SetPrintLevel(0);
  tm->SetFCN(Gmy_mnf_lk);
  //  tm->SetFCN(GetChi2);

  for (iIon=0; iIon<NION; iIon++) {
    scale = SLIDER_WEIGHT_SCALE;
    t8 = scale*fSliders[iIon]->GetPosition();

    //    tm->mnparm(iIon, ionNames[iIon], t8, fabs(0.1*t8), 0, 10, ierr);
    tm->mnparm(iIon, ionNames[iIon], t8, 0.1, -1e-4, 10, ierr);

    if (ionIsActive[iIon]==0) {
      tm->FixParameter(iIon);
    }

    currentWeightVals[iIon] = t8;
  }
  
  ierr = 1;
  int ncall = 0;
  while(ierr && ncall<1) {
    tm->mnexcm("MIGRAD",0,0,ierr);
    //    tm->mnimpr();
    ncall++;
  }

  for (iIon=0;iIon<NION;iIon++) {
    tm->GetParameter(iIon,rval,rerr);
    *(fittedPars+iIon) = rval;
    *(fittedParErrors+iIon) = rerr;

    if (ionIsActive[iIon]) {
      fSliders[iIon]->SetPosition((int)(rval/scale));
      sprintf(buf, "%.3f", fSliders[iIon]->GetPosition()*scale);
      fSliderTextBuffers[iIon]->Clear();
      fSliderTextBuffers[iIon]->AddText(0, buf);
      gClient->NeedRedraw(fSliderTextEntries[iIon]);
      sprintf(thisStr, "%04d", ionList[iIon]);
      printf("fittedVal iIon= %5s name= %6s --> %+.6e +- %+.6e | chi2= %.6e \n", thisStr, ionNames[iIon], rval, rerr, GfChi2);      
    }
  }
  printf("Minuit iflag= %d \n", ierr); 
  DoDraw();
  
}

//______________________________________________________________________________
void TSynGui::SetActiveIonFile(const char* s) {
  sprintf(activeIonFile, "%s", s);
}

//______________________________________________________________________________
void TSynGui::GetModelFlux(int nW, double *a) {
  int iIon, iW;
  double tmp;

  for (iW=0; iW<nW; iW++) {
    tmp = 0.0;
    for (iIon=1; iIon<NION; iIon++) {
      if (ionIsActive[iIon]) {
	//	printf("ION %d is ACTIVE!\n", iIon);
	//	tmp += SLIDER_WEIGHT_SCALE*fSliders[iIon]->GetPosition()*ionFlux[iIon][iW];
	//	tmp += currentWeightVals[iIon]*TMath::Log(ionFlux[iIon][iW]);
	tmp += currentWeightVals[iIon]*ionFlux[iIon][iW];
      }
    }
    //    printf("a1= %.6e \n", *(a+iW));
    tmp = exp(tmp);
    //    printf("a2= %.6e \n", *(a+iW));

    //    tmp *= SLIDER_WEIGHT_SCALE*fSliders[IBB]->GetPosition();
    //    tmp = SLIDER_WEIGHT_SCALE*fSliders[IBB]->GetPosition()*ionFlux[IBB][iW];
    //    tmp *= SLIDER_WEIGHT_SCALE*fSliders[IBB]->GetPosition()*ionFlux[IBB][iW];
    tmp *= currentWeightVals[IBB]*ionFlux[IBB][iW];
    *(a+iW) = tmp;
    //    printf("a3= %.6e \n", *(a+iW));
  }

}

//______________________________________________________________________________
void TSynGui::SetSNFile(const char* s) {
  sprintf(snFile, "%s", s);
}

//______________________________________________________________________________
void TSynGui::MultiplyArray(int n, double *oldA, double *newA, double multiplier) {

  for(int i=0; i<n; i++) {
    *(newA+i) = *(oldA+i)*multiplier;
  }
  
}


//______________________________________________________________________________
void TSynGui::ReadActiveIonFile() {
  FILE *ifp;
  char line[256], str[32];
  int iIon, i;

  ifp = fopen(activeIonFile, "r");
  if (ifp == NULL ) {
    fprintf(stderr, "ERROR: CANT OPEN FILE %s! \n", activeIonFile);
    //    exit(0);
  }

  i = 0;
  while(fgets(line,256,ifp)!=NULL) {
    //    printf("line %s \n", line);
    if (line[0]=='#') {
      continue;
    }
    sscanf(line,"%d %s", &iIon, str);
    printf("%6s %5d \n", str, iIon);
    ionList[i] = iIon;
    sprintf(ionNames[i] ,"%s", str);
    ionIsActive[i] = 0;
    i++;
  }

  ionIsActive[IBB] = 1;
  NION = i;
  printf("READ %d IONS.\n", NION);

}

//______________________________________________________________________________
void TSynGui::ScaleToMean(int n, double *a, double target) {
  double mean;

  mean = GetMean(n, a);
  if (mean<=0) {
    printf("ERROR: MEAN <=0 (%.6e) \n", mean);
    //    exit(0);
  }
  
  for (int i=0; i<n;i++ ) {
    *(a+i) = a[i]*target/mean;
  }

}

//______________________________________________________________________________
void TSynGui::SetRefsFile(const char* s) {
  sprintf(refsFile, "%s", s);
}


//______________________________________________________________________________
TSynGui::TSynGui() : TGMainFrame(gClient->GetRoot(), 10, 10)
{

  double fMean1, fMean2, eMean1, eMean2;
  char buf[32];

  tphotValues[0] = 5.0;
  tphotValues[1] = 8.0;
  tphotValues[2] = 11.0;

  etempValues[0] = 5.0;
  etempValues[1] = 8.0;
  etempValues[2] = 11.0;

  logtauValues[0] = -2.0;
  logtauValues[1] = -1.0;
  logtauValues[2] = +0.0;
  logtauValues[3] = +1.0;
  logtauValues[4] = +2.0;
  logtauValues[5] = +3.0;
  logtauValues[6] = +4.0;
  logtauValues[7] = +5.0;

  currentTphotValue = tphotValues[0];
  currentEtempValue = etempValues[0];
  currentLogtauValue = logtauValues[2];

  minXForFit = 3500.0;
  maxXForFit = 8000.0;
  //  SetActiveIonFile("./APPDATA/activeIonList.txt");
  SetActiveIonFile("./APPDATA/allIons.dat");

  SetSNFile("./APPDATA/1991bg_19911214_3640_9282_00.dat");
  //  SetSNFile("./APPDATA/2004dt_20040908_3432_7704_00.dat");
  //  SetSNFile("/Users/bdilday/data/otherSNe/PTF12bho/12bho_20120315_Keck1_v1.ascii");
  //  SetSNFile("./myAPPDATA/11kx_20110326_Keck1_v1.ascii");
  ReadActiveIonFile();

  ReadSNDATA();
  fMean1 = GetMean(NSNWAVELENGTHS, snFlux);
  eMean1 = GetMean(NSNWAVELENGTHS, snFluxError);
  printf("MEAN IS %.12e %.12e\n", fMean1, eMean1);
  ScaleToMean(NSNWAVELENGTHS, snFlux, 1.0);
  fMean2 = GetMean(NSNWAVELENGTHS, snFlux);
  MultiplyArray(NSNWAVELENGTHS, snFluxError, snFluxError, fMean2/fMean1);
  eMean2 = GetMean(NSNWAVELENGTHS, snFluxError);
  printf("MEAN IS %.12e %.12e\n", fMean2, eMean2);


  ReadAPPDATA();
  printf("MEAN IS %.12e\n", GetMean(NIONWAVELENGTHS, ionFlux[IBB]));
  ScaleToMean(NIONWAVELENGTHS, ionFlux[IBB], 1.0);
  printf("MEAN IS %.12e\n", GetMean(NIONWAVELENGTHS, ionFlux[IBB]));

  SetRefsFile("/Users/bdilday/local/share/es/refs.dat");

  SetCleanup(kDeepCleanup);
  // Create an embedded canvas and add to the main frame, centered in x and y
  // and with 30 pixel margins all around
  fCanvas = new TRootEmbeddedCanvas("Canvas", this, 240, 180);
  fLcan = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 10);
  //  fLcan = new TGLayoutHints(kLHintsExpandY, 10, 10, 10, 10);
  AddFrame(fCanvas, fLcan);
  fCanvas->GetCanvas()->SetFillColor(33);
  fCanvas->GetCanvas()->SetFrameFillColor(41);
  fCanvas->GetCanvas()->SetBorderMode(0);
  fCanvas->GetCanvas()->SetGrid();
  //  fCanvas->GetCanvas()->SetLogy();
  
  int irow, icol;
  int SLIDER_WIDTH = 100;

  /**
  **/
  fSliderFrameRight = new TGHorizontalFrame(this, 0, 0, 0);
  fSliderFrameRight->Resize(150, 100);

  
  
  /**************************/
  fSliders[IZ1] = new TGHSlider(fSliderFrameRight, 
				SLIDER_WIDTH, kScaleBoth, FSLIDERID+IZ1,
				kHorizontalFrame,
				GetDefaultFrameBackground());
  
  fSliders[IZ1]->Connect("PositionChanged(Int_t )", "TSynGui", 
			 this, "DoSlider()");
  
  fSliders[IZ1]->SetRange(0,(int)(1.0/SLIDER_Z1_SCALE));
  
  fSliderTextEntries[IZ1] = new TGTextEntry(fSliderFrameRight, 
					    fSliderTextBuffers[IZ1]=new TGTextBuffer(5), FTEXTID+IZ1);

  fSliderTextBuffers[IZ1]->AddText(0, "0.0");
  fSliders[IZ1]->SetPosition(0);
  fSliderTextEntries[IZ1]->Connect("TextChanged(char*)", "TSynGui", this,
			      "DoText(char*)");

  /**************************/

  fSliders[ISNSCALE] = new TGHSlider(fSliderFrameRight, 
				SLIDER_WIDTH, kScaleBoth, FSLIDERID+ISNSCALE,
				kHorizontalFrame,
				GetDefaultFrameBackground());
  
  fSliders[ISNSCALE]->Connect("PositionChanged(Int_t )", "TSynGui", 
			      this, "DoSlider()");
  
  fSliders[ISNSCALE]->SetRange(0,(int)(10.0/SLIDER_WEIGHT_SCALE));
  //fSliders[ISNSCALE]->SetRange(0,(int)(20.0/SLIDER_WEIGHT_SCALE));
  
  fSliderTextEntries[ISNSCALE] = new TGTextEntry(fSliderFrameRight, 
						 fSliderTextBuffers[ISNSCALE]=new TGTextBuffer(5), FTEXTID+ISNSCALE);

  fSliderTextBuffers[ISNSCALE]->AddText(0, "1.0");
  fSliderTextEntries[ISNSCALE]->Connect("TextChanged(char*)", "TSynGui", this,
					"DoText(char*)");
  fSliders[ISNSCALE]->SetPosition((int)(1.0/SLIDER_WEIGHT_SCALE));

  /*****************/
  fSliderFrameRadio1 = new TGHorizontalFrame(this, 0, 0, 0);
  fSliderFrameRadio1->Resize(150, 100);

  //  TGHorizontalFrame *thisFrame = fSliderFrameRight;
  TGHorizontalFrame *thisFrame = fSliderFrameRadio1;

  TGTextButton *fTextButtonChi2 = new TGTextButton(thisFrame, "GetChi2");
  //  fTextButtonChi2->Connect("Clicked()", "TSynGui", this, "GGetChi2()");
  fTextButtonChi2->Connect("Clicked()", "TSynGui", this, "GetChi2()");

  fLabelChi2 = new TGLabel(thisFrame, "----------------------");

  /*****************/
  TGTextButton *fTextButtonLinear = new TGTextButton(thisFrame, "LinearBBFit");
  //  fTextButtonChi2->Connect("Clicked()", "TSynGui", this, "GGetChi2()");
  fTextButtonLinear->Connect("Clicked()", "TSynGui", this, "LinearBBFit()");

  /*****************/
  TGTextButton *fTextButtonDoFit = new TGTextButton(thisFrame, "DoFit");
  fTextButtonDoFit->Connect("Clicked()", "TSynGui", this, "DoFit()");

  /*****************/
  TGTextButton *fTextButtonClearAll = new TGTextButton(thisFrame, "ClearAll");
  fTextButtonClearAll->Connect("Clicked()", "TSynGui", this, "ClearAll()");

  /*****************/
  TGTextButton *fTextButtonSaveFigure = new TGTextButton(thisFrame, "SaveFigure");
  fTextButtonSaveFigure->Connect("Clicked()", "TSynGui", this, "SaveFigure()");

  /*****************/
  //  TGTextButton *fTextButtonSaveFile = new TGTextButton(thisFrame, "SaveFile");
  //  fTextButtonSaveFile->Connect("Clicked()", "TSynGui", this, "SaveFile()");


  /*****************/
  fXminText = new TGTextEntry(fSliderFrameRight, 
			      fXminTextBuffer=new TGTextBuffer(6), FXTEXTMIN);
  fXmaxText = new TGTextEntry(fSliderFrameRight, 
			      fXmaxTextBuffer=new TGTextBuffer(6), FXTEXTMAX);

  sprintf(buf, "%.1f", minXForFit);
  fXminTextBuffer->AddText(0, buf);
  sprintf(buf, "%.1f", maxXForFit);
  fXmaxTextBuffer->AddText(0, buf);

  fXminText->Connect("TextChanged(char*)", "TSynGui", this,
		     "DoXFitText(char*)");
  fXmaxText->Connect("TextChanged(char*)", "TSynGui", this,
		     "DoXFitText(char*)");

  /**
  fSliderFrameRight->AddFrame(fSliders[IZ1]);
  fSliderFrameRight->AddFrame(fSliderTextEntries[IZ1]);
  fSliderFrameRight->AddFrame(fSliders[ISNSCALE]);
  fSliderFrameRight->AddFrame(fSliderTextEntries[ISNSCALE]);
  **/

  TGLayoutHints *tgl1 = new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1);
  fSliderFrameRight->AddFrame(fSliders[IZ1], tgl1);
  fSliderFrameRight->AddFrame(fSliderTextEntries[IZ1], tgl1);
  fSliderFrameRight->AddFrame(fSliders[ISNSCALE], tgl1);
  fSliderFrameRight->AddFrame(fSliderTextEntries[ISNSCALE], tgl1);

  /**
  fSliderFrameRight->AddFrame(fTextButtonDoFit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
  fSliderFrameRight->AddFrame(fTextButtonClearAll, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
  fSliderFrameRight->AddFrame(fTextButtonLinear, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
  fSliderFrameRight->AddFrame(fTextButtonSaveFigure, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
  **/

  fSliderFrameRight->AddFrame(fXminText, tgl1);
  fSliderFrameRight->AddFrame(fXmaxText, tgl1);

  AddFrame(fSliderFrameRight);  



  fSliderFrameRadio1->AddFrame(fTextButtonChi2, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
  fSliderFrameRadio1->AddFrame(fLabelChi2);
  fSliderFrameRadio1->AddFrame(fTextButtonDoFit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
  fSliderFrameRadio1->AddFrame(fTextButtonClearAll, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
  fSliderFrameRadio1->AddFrame(fTextButtonLinear, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
  fSliderFrameRadio1->AddFrame(fTextButtonSaveFigure, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
  //  fSliderFrameRadio1->AddFrame(fTextButtonSaveFile, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));



  AddFrame(fSliderFrameRadio1);  

  /**************************/
  /** radio buttons to control t_phot and electron temperture **/

  /**
  fRadioGroupTphot = new TGButtonGroup(fSliderFrameRadio1, "t_phot", kHorizontalFrame);
  fRadioGroupEtemp = new TGButtonGroup(fSliderFrameRadio1, "e_temp", kHorizontalFrame);
  **/

  fRadioGroupTphot = new TGButtonGroup(fSliderFrameRight, "t_phot", kHorizontalFrame);
  fRadioGroupEtemp = new TGButtonGroup(fSliderFrameRight, "e_temp", kHorizontalFrame);


  //  fComboBox = new TGComboBox(fSliderFrameRight, "logtau", FCOMBOBOXID);
  sprintf(buf, "%+.2f", currentLogtauValue);
  fComboBox = new TGComboBox(fSliderFrameRadio1, buf, FCOMBOBOXID);

  for (int i=0;i<NLOGTAU;i++) {
    sprintf(buf, "%+6.2f", logtauValues[i]);
    fComboBox->AddEntry(buf, FCOMBOBOXID+i+1);
  }

  fComboBox->Select(FCOMBOBOXID+2);
  fComboBox->Resize(80, 20);
  fComboBox->Connect("Selected(Int_t)", "TSynGui", this, "DoLogtauCombo()");
 
  for (int i=0;i<NTPHOT;i++) {
    sprintf(buf, "%.1f", tphotValues[i]);
    fRadioTphot[i] = new TGRadioButton(fRadioGroupTphot, buf, FRADIOTPHOTID+i);
  } 
  fRadioTphot[0]->SetState(kButtonDown);

  for (int i=0;i<NETEMP;i++) {
    sprintf(buf, "%.1f", etempValues[i]);
    fRadioEtemp[i] = new TGRadioButton(fRadioGroupEtemp, buf, FRADIOETEMPID+i); 
  }
  fRadioEtemp[0]->SetState(kButtonDown);

  //  fSliderFrameRadio1->AddFrame(fRadioTphot, tgl1);
  /**
  fRadioGroupTphot->AddFrame(fRadioTphot[0], tgl1);
  fRadioGroupTphot->AddFrame(fRadioTphot[1], tgl1);
  fRadioGroupTphot->AddFrame(fRadioTphot[2], tgl1);
  **/

  fRadioGroupTphot->Connect("Pressed(Int_t)", "TSynGui", this, 
			    "DoRadio(Int_t)");

  fRadioGroupEtemp->Connect("Pressed(Int_t)", "TSynGui", this, 
			    "DoRadio(Int_t)");

  /**
  **/
  fSliderFrameRight->AddFrame(fRadioGroupTphot, tgl1);
  fSliderFrameRight->AddFrame(fRadioGroupEtemp, tgl1);

  //  fSliderFrameRight->AddFrame(fComboBox, tgl1);
  fSliderFrameRadio1->AddFrame(fComboBox, tgl1);


  //  AddFrame(fSliderFrameRadio1);  

  fRadioGroupTphot->SetLayoutHints(tgl1);
  fRadioGroupEtemp->SetLayoutHints(tgl1);
  fRadioGroupTphot->Show();
  fRadioGroupEtemp->Show();



  /**************************/
  //  int NROW = 8;
  int NCOL = 5;
  char thisStr[64];

  /**
     TGHSlider *fSliders[MAXION];
     TGHorizontalFrame *fSliderFrames[MAXION];
     TGLayoutHints     *fSliderBlys[MAXION];
     TGTextEntry *fSliderTextEntries[MAXION];
     TGCheckButton *fSliderCheckButtons[MAXION];
  **/
  

  for (Int_t i=0; i<NION; i++ ) {
    //    fSliders[i] = 
    irow = i/NCOL;
    icol = i % NCOL;
    if (i % NCOL == 0) {
      fSliderFrames[irow] = new TGHorizontalFrame(this, 0, 0, 0);
      fSliderFrames[irow]->Resize(100, 25);
    }

    fSliders[i] = new TGHSlider(fSliderFrames[irow], 
				SLIDER_WIDTH, kScaleBoth, FSLIDERID+i,
				kHorizontalFrame,
				GetDefaultFrameBackground());
    fSliderBlys[i] = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5);
    fSliderBflys[i] = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 5, 5, 5, 5);
    
    //    sprintf(thisStr, "%6s", ionNames[i]);
    sprintf(thisStr, "%04d", ionList[i]);
    fSliderCheckButtons[i] = new TGCheckButton(fSliderFrames[irow], 
					       thisStr, FCHECKID+i);
    fSliderCheckButtons[i]->SetToolTipText(ionNames[i]);

    fSliders[i]->Connect("PositionChanged(Int_t )", "TSynGui", 
			 this, "DoSlider()");

    //    fSliders[i]->SetRange(0,(int)(5.0/SLIDER_WEIGHT_SCALE));
    fSliders[i]->SetRange(0,(int)(20.0/SLIDER_WEIGHT_SCALE));

    fSliderTextEntries[i] = new TGTextEntry(fSliderFrames[irow], 
					    fSliderTextBuffers[i]=new TGTextBuffer(5), FTEXTID+i);
    fSliderTextEntries[i]->SetToolTipText(ionNames[i]);

    //    fSliderTextBuffers[i]->AddText(0, "0.0");
    fSliderTextEntries[i]->Connect("TextChanged(char*)", "TSynGui", this,
				   "DoText(char*)");
    
    fSliderCheckButtons[i]->Connect("Clicked()", "TSynGui", this,
				    "HandleButtons()");
  }


  for (Int_t i=0; i<NION; i++ ) {
    //    fMain->AddFrame(fSliderFrames[i]);
    irow = i/NCOL;
    icol = i % NCOL;
    fSliderFrames[irow]->AddFrame(fSliderCheckButtons[i]);
    fSliderFrames[irow]->AddFrame(fSliderTextEntries[i]);
    fSliderFrames[irow]->AddFrame(fSliders[i]);
  }

  for (Int_t i=0; i<NION; i++ ) {
    irow = i/NCOL;
    icol = i % NCOL;
    if (icol==0) {
      AddFrame(fSliderFrames[irow], fSliderBlys[i]);
    }
  }

  // Set main frame name, map sub windows (buttons), initialize layout
  // algorithm via Resize() and map main frame
  SetWindowName("Interactive SYN++ with ROOT");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();
  
  //   fFitFcn = new TF1("fFitFcn", "TMath::LogNormal(x, [0], [1], [2])", 0, 5);

  
  fCanvas->GetCanvas()->SetGrid();  
    
  //  modelGraph->Draw();

  Int_t i;
  for (i=1; i<NION; i++ ) {  
    fSliders[i]->SetPosition(0);
    currentWeightVals[i] = 0.0;
    sprintf(buf, "%.3f", SLIDER_WEIGHT_SCALE*fSliders[i]->GetPosition());
    fSliderTextBuffers[i]->Clear();
    fSliderTextBuffers[i]->AddText(0, buf);
  }

  i = 0;
  currentWeightVals[i] = 1.0;
  fSliders[i]->SetPosition((int)(1.0/SLIDER_WEIGHT_SCALE));
  fSliderTextBuffers[i]->Clear();
  fSliderTextBuffers[i]->AddText(0, "1.000");

  DoDraw();

}

//______________________________________________________________________________
double TSynGui::GetMean(int n, double *a) {
  double sum = 0.0;
  double mean;
  for (int i=0; i<n; i++ ) {
    sum += a[i];
  }

  mean = sum/(1.0*n);
  return (mean);
}

//______________________________________________________________________________
TSynGui::~TSynGui()
{
   // Clean up

   Cleanup();
}

//______________________________________________________________________________
void TSynGui::CloseWindow()
{
   // Called when window is closed via the window manager.

   delete this;
}

//______________________________________________________________________________
void TSynGui::ReadSNDATA() {
  FILE *ifp;
  char line[256];
  double lam, flux, fluxError;
  double errfloor = 5e-16;

  NSNWAVELENGTHS = 0;
  ifp = fopen(snFile, "r");
  if (ifp == NULL) {
    printf("WARNING: CANNOT SN FILE %s. SKIP. \n", snFile);
    exit(0);
  } else {
    while(fgets(line,256,ifp)!=NULL) {
      sscanf(line,"%lf %lf %lf ", &lam, &flux, &fluxError);
      //      printf("%lf %.12e \n", lam, flux);
      snWavelengths[NSNWAVELENGTHS] = lam;
      snFlux[NSNWAVELENGTHS] = flux;
      if ( fluxError>0 ) {
	snFluxError[NSNWAVELENGTHS] = fluxError;
      } else {
	if (lam<5000) {
	  snFluxError[NSNWAVELENGTHS] = sqrt( (0.2*flux)*(0.2*flux) + errfloor*errfloor);
	} else {
	  snFluxError[NSNWAVELENGTHS] = sqrt( (0.1*flux)*(0.1*flux) + errfloor*errfloor);
	}
      }

      NSNWAVELENGTHS++;
    }
    fclose(ifp);
  }

  printf("MEAN IS %.12e\n", GetMean(NSNWAVELENGTHS, snFlux));

}

//______________________________________________________________________________
void TSynGui::ReadAPPDATA() {
  FILE *ifp;
  char line[256];
  char ifile[1024];
  double lam, flux;

  t_phot = currentTphotValue;
  etemp  = currentEtempValue;
  log_tau = currentLogtauValue;

  for (int i=0; i<NION; i++ ){
    sprintf(ifile, "./APPDATA/ion%04d_%05.2f_%05.2f_%+06.3f.dat", 
	    ionList[i], t_phot, etemp, log_tau);
    ifp = fopen(ifile, "r");
    if (ifp == NULL) {
      printf("WARNING: CANNOT FILE %s. SKIP. \n", ifile);
    } else {
      NIONWAVELENGTHS = 0;
      while(fgets(line,256,ifp)!=NULL) {
	sscanf(line,"%lf %lf ", &lam, &flux);
	//	printf("%lf %.12e \n", lam, flux);
	ionWavelengths[i][NIONWAVELENGTHS] = lam;
	if (i==0) {
	  ionFlux[i][NIONWAVELENGTHS] = flux;
	} else {
	  ionFlux[i][NIONWAVELENGTHS] = TMath::Log(flux);
	}
	NIONWAVELENGTHS++;
      }
      fclose(ifp);
    }
  }

}


//______________________________________________________________________________
void TSynGui::DoRadio(const Int_t iradio) {
  printf("DoRadio called with argument %d \n", iradio);
  int FID = (iradio/100)*100;
  int index = iradio-FID;

  if (FID==FRADIOTPHOTID) {
    currentTphotValue = tphotValues[index];
  } else if (FID==FRADIOETEMPID) {
    currentEtempValue = etempValues[index];
  }

  ReadAPPDATA();
  //  DoDraw();
  GetChi2();

}

//______________________________________________________________________________
void TSynGui::DoLogtauCombo() {
  int index = fComboBox->GetSelected()-1-FCOMBOBOXID;
  currentLogtauValue = logtauValues[index];
  ReadAPPDATA();
  GetChi2();

}


//______________________________________________________________________________
void TSynGui::DoText(const char * /*text*/)
{
  // Handle text entry widgets.
  
  TGTextEntry *te = (TGTextEntry *) gTQSender;
  Int_t textId = te->WidgetId();
  Int_t index;
  double t8;

  index = textId - FTEXTID;
  printf("textid= %d index= %d \n", textId, index);

  if (index==IZ1) {
    t8 = atof(fSliderTextBuffers[index]->GetString());
    fSliders[index]->SetPosition((int)(t8/SLIDER_Z1_SCALE));
  } else {
    t8 = atof(fSliderTextBuffers[index]->GetString());
    fSliders[index]->SetPosition((int)(t8/SLIDER_WEIGHT_SCALE));
    currentWeightVals[index] = t8;
  }


  DoDraw();

}

//______________________________________________________________________________
void TSynGui::DoXFitText(const char *)
{
  // Handle text entry widgets.
  
  //  TGTextEntry *te = (TGTextEntry *) gTQSender;
  //  Int_t textId = te->WidgetId();
  
  minXForFit = atof(fXminTextBuffer->GetString());
  maxXForFit = atof(fXmaxTextBuffer->GetString());
  printf("new Xmin= %.1f Xmax = %.1f \n", minXForFit, maxXForFit);

}

//______________________________________________________________________________
void TSynGui::DoSlider()
{
  // Handle slider widgets.
  char buf[32];
  TGSlider *te = (TGSlider *) gTQSender;
  Int_t sliderId = te->WidgetId();
  Int_t index;

  double scale, t8;
  
  index = sliderId - FSLIDERID;

  if (index==IZ1) {
    scale = SLIDER_Z1_SCALE;  
    sprintf(buf, "%.3f", scale*fSliders[index]->GetPosition());
  } else {
    scale = SLIDER_WEIGHT_SCALE;  
    t8 = scale*fSliders[index]->GetPosition();
    sprintf(buf, "%.3f", t8);
    currentWeightVals[index] = t8;
  }

  fSliderTextBuffers[index]->Clear();
  fSliderTextBuffers[index]->AddText(0, buf);
  fSliderTextEntries[index]->Deselect();
  gClient->NeedRedraw(fSliderTextEntries[index]);

  DoDraw();

}


//______________________________________________________________________________
void TSynGui::DoDraw() {
  int iW;
  double t8, z1, target;
  double xx[MAXWAVELENGTHS];
  double yy[MAXWAVELENGTHS];
  double yye[MAXWAVELENGTHS];
  double yyeUpper[MAXWAVELENGTHS];
  double yyeLower[MAXWAVELENGTHS];

  t8 = SLIDER_WEIGHT_SCALE*fSliders[ISNSCALE]->GetPosition();
  z1 = 1.0+SLIDER_Z1_SCALE*fSliders[IZ1]->GetPosition();

  MultiplyArray(NSNWAVELENGTHS, snWavelengths, xx, 1.0/z1);
  MultiplyArray(NSNWAVELENGTHS, snFlux, yy, t8);
  MultiplyArray(NSNWAVELENGTHS, snFluxError, yye, t8);

  snGraph = new TGraph(NSNWAVELENGTHS, xx, yy);
  snGraph->Clear();
  snGraph->Draw("AL");
  snGraph->GetXaxis()->SetRangeUser(2500, 12000);
  snGraph->GetYaxis()->SetRangeUser(0,5);

  for (iW=0;iW<NSNWAVELENGTHS;iW++) {
    yyeUpper[iW] = yy[iW]+yye[iW];
    yyeLower[iW] = yy[iW]-yye[iW];
  }

  errorGraphUpper = new TGraph(NSNWAVELENGTHS, xx, yyeUpper);
  errorGraphLower = new TGraph(NSNWAVELENGTHS, xx, yyeLower);
  errorGraphUpper->SetLineStyle(3);
  errorGraphLower->SetLineStyle(3);
  errorGraphUpper->Draw("L");
  errorGraphLower->Draw("L");

  target = SLIDER_WEIGHT_SCALE*fSliders[IBB]->GetPosition();
  GetModelFlux(NIONWAVELENGTHS, modelFlux);
  ScaleToMean(NIONWAVELENGTHS, modelFlux, target);
  modelGraph = new TGraph(NIONWAVELENGTHS, ionWavelengths[IBB], modelFlux);  
  modelGraph->SetLineColor(kRed);
  //  modelGraph->SetLineWidth(2);
  //  modelGraph->Draw("L SAME");
  modelGraph->Draw("L");


  fCanvas->GetCanvas()->Update();

}


//______________________________________________________________________________
void TSynGui::HandleButtons()
{
   // Handle different buttons.
  int index;
  TGButton *btn = (TGButton *) gTQSender;
  Int_t id = btn->WidgetId();

  index = id-FCHECKID;
  
  //  printf("button %d state %d \n", index, (int)btn->GetState());

  ionIsActive[index] = (int)btn->GetState();

  DoDraw();

}


//______________________________________________________________________________
void TSynGui::ClearAll() {
  int iIon;
  char buf[32];

  fSliders[0]->SetPosition((int)(1.0/SLIDER_WEIGHT_SCALE));
  fSliderCheckButtons[0]->SetState(kButtonDown);
  ionIsActive[0] = 1;

  for (iIon=1;iIon<NION;iIon++) {
    fSliders[iIon]->SetPosition(0);
    fSliderCheckButtons[iIon]->SetState(kButtonUp);
    ionIsActive[iIon] = 0;
  }
    
  for (iIon=0;iIon<NION;iIon++) {
    sprintf(buf, "%.3f", fSliders[iIon]->GetPosition()*SLIDER_WEIGHT_SCALE);
    fSliderTextBuffers[iIon]->Clear();
    fSliderTextBuffers[iIon]->AddText(0, buf);
    gClient->NeedRedraw(fSliderTextEntries[iIon]);
  }

  DoDraw();

}


//______________________________________________________________________________
void GGetChi2() {
  double t8, z1, target;
  double xx[MAXWAVELENGTHS];
  double yy[MAXWAVELENGTHS];
  double yye[MAXWAVELENGTHS];
  char sChi2[64];

  t8 = SLIDER_WEIGHT_SCALE*Gsyngui->fSliders[ISNSCALE]->GetPosition();
  z1 = 1.0+SLIDER_Z1_SCALE*Gsyngui->fSliders[IZ1]->GetPosition();
  
  Gsyngui->MultiplyArray(Gsyngui->NSNWAVELENGTHS, Gsyngui->snWavelengths, xx, 1.0/z1);
  Gsyngui->MultiplyArray(Gsyngui->NSNWAVELENGTHS, Gsyngui->snFlux, yy, t8);
  Gsyngui->MultiplyArray(Gsyngui->NSNWAVELENGTHS, Gsyngui->snFluxError, yye, t8);

  //  TGraph *snGraph = new TGraph(Gsyngui->NSNWAVELENGTHS, xx, yy);
 
  //  target = SLIDER_WEIGHT_SCALE*Gsyngui->fSliders[IBB]->GetPosition();
  target = currentWeightVals[IBB];

  printf("GGetChi2 target = %.6e \n", target);

  Gsyngui->GetModelFlux(Gsyngui->NIONWAVELENGTHS, Gsyngui->modelFlux);
  Gsyngui->ScaleToMean(Gsyngui->NIONWAVELENGTHS, Gsyngui->modelFlux, target);
  TGraph *modelGraph = new TGraph(Gsyngui->NIONWAVELENGTHS, 
				  Gsyngui->ionWavelengths[IBB], 
				  Gsyngui->modelFlux);  

  GfChi2 = 0.0;
  double chi;
  double mm;
  for (int iW=0;iW<Gsyngui->NSNWAVELENGTHS;iW++) {
    if(xx[iW]>=Gsyngui->minXForFit && xx[iW]<=Gsyngui->maxXForFit) {    
      mm = modelGraph->Eval(xx[iW]);
      chi = (mm-yy[iW])/(yye[iW]);    
      //    chi = (mm-yy[iW]);    
      GfChi2 += chi*chi;
      if (fabs(xx[iW]-4192.00)<1e-2) {
	printf("w= %.3f model= %.6e data= %.6e +- %.6e chi= %.6e chi2= %.6e \n", xx[iW], mm, yy[iW], yye[iW], chi, GfChi2);
      }
    }
  }

  //  printf("Chi2= %.12e \n", GfChi2);
  sprintf(sChi2,"%.9e", GfChi2);
  Gsyngui->fLabelChi2->SetText(sChi2);

  Gsyngui->DoDraw();

}

//______________________________________________________________________________
void TSynGui::GetChi2() {
  char sChi2[64];

  GGetChi2();
  printf("Chi2= %.12e \n", GfChi2);
  sprintf(sChi2,"%.9e", GfChi2);
  fLabelChi2->SetText(sChi2);

}

//______________________________________________________________________________
void Gmy_mnf_lk(int &npar, double *gin, double &f, double *x, int iflag) {

  int iIon;
  double val8;
  int ival;

  for (iIon=0;iIon<Gsyngui->NION;iIon++) {
    //    val8 = SLIDER_WEIGHT_SCALE*Gsyngui->fSliders[iIon]->GetPosition();
    currentWeightVals[iIon] = x[iIon];
    val8 = currentWeightVals[iIon];
    ival = ((int)(x[iIon]/SLIDER_WEIGHT_SCALE));
    if (Gsyngui->ionIsActive[iIon]>=1) {
      printf("Gmy_mnf_lk iIon= %4d ionName= %6s val8= %+.6e xIon= %.6e \n", iIon, Gsyngui->ionNames[iIon], val8, x[iIon]);
    }
    Gsyngui->fSliders[iIon]->SetPosition(ival);
  }

  GGetChi2();
  //  printf("x[0]= %.12e chi2= %.12e \n", x[0], GfChi2);
  f = GfChi2;

}

//______________________________________________________________________________
void synGuiRoot() 
{
  //  TSynGui *tsg = new TSynGui();
  Gsyngui = new TSynGui();
  GGetChi2();

  //  tsg->ReadAPPDATA();

}

