#include "LTimeAnalysis.h"

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMath.h"
#include "TString.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "MRText.h"

#include "RDK2Constants.h"

using namespace std;

void analyzeProtonBackscattering(TString idString)
{
	//Read in results file and tree
	TString filePath=TString(MRKGEANTDATA_DIR)+"Results_"+idString+"_Set0.root";
	cout << "Analyzing proton backscattering from flux tree of : " << filePath << endl;

	TFile resultsFile(filePath,"READ");

	TTree* sbdFluxTree=(TTree*) resultsFile.Get("SBDFluxTree");
	TTree* resultsTree=(TTree*) resultsFile.Get("g");
	int numEvents=resultsTree->GetEntries();

	Int_t event;
	Char_t hitCode,particle;
	Double_t time, posX,posY,posZ,momX,momY,momZ;
	Double_t KE=0;

	sbdFluxTree->SetBranchAddress("Event",&event);
	sbdFluxTree->SetBranchAddress("HitCode",&hitCode);
	sbdFluxTree->SetBranchAddress("Particle",&particle);
	sbdFluxTree->SetBranchAddress("KE",&KE);
	sbdFluxTree->SetBranchAddress("Time",&time);
	sbdFluxTree->SetBranchAddress("PosX",&posX);
	sbdFluxTree->SetBranchAddress("PosY",&posY);
	sbdFluxTree->SetBranchAddress("PosZ",&posZ);
	sbdFluxTree->SetBranchAddress("MomX",&momX);
	sbdFluxTree->SetBranchAddress("MomY",&momY);
	sbdFluxTree->SetBranchAddress("MomZ",&momZ);

	//Make reflect file to collect results of this analysis
	TString reflectFilePath=TString(MRKGEANTDATA_DIR)+"Results_"+idString+"_Reflect_Set0.root";
	TFile reflectFile(reflectFilePath,"RECREATE");
	reflectFile.cd();
	TTree* reflectTree=new TTree("reflectTree","Records reflections off SBD");

	Int_t refEvent=-1;
	Int_t NumExits=0;
	Int_t NumEntries=0;
	Double_t KEEntry1=0;
	Double_t KEExit1=0;
	Double_t KEEntry2=0;
	Double_t PosXEntry1=0;
	Double_t PosYEntry1=0;
	Double_t PosZEntry1=0;
	Double_t PosXEntry2=0;
	Double_t PosYEntry2=0;
	Double_t PosZEntry2=0;
	Double_t MomXEntry1=0;
	Double_t MomYEntry1=0;
	Double_t MomZEntry1=0;
	Double_t MomXEntry2=0;
	Double_t MomYEntry2=0;
	Double_t MomZEntry2=0;

	reflectTree->Branch("Event",&refEvent,"Event/I");
	reflectTree->Branch("NumEntries",&NumEntries,"NumEntries/I");
	reflectTree->Branch("NumExits",&NumExits,"NumExits/I");
	reflectTree->Branch("KEEntry1",&KEEntry1,"KEEntry1/D");
	reflectTree->Branch("KEExit1",&KEExit1,"KEExit1/D");
	reflectTree->Branch("KEEntry2",&KEEntry2,"KEEntry2/D");
	reflectTree->Branch("PosXEntry1",&PosXEntry1,"PosXEntry1/D");
	reflectTree->Branch("PosYEntry1",&PosYEntry1,"PosYEntry1/D");
	reflectTree->Branch("PosZEntry1",&PosZEntry1,"PosZEntry1/D");
	reflectTree->Branch("PosXEntry2",&PosXEntry2,"PosXEntry2/D");
	reflectTree->Branch("PosYEntry2",&PosYEntry2,"PosYEntry2/D");
	reflectTree->Branch("PosZEntry2",&PosZEntry2,"PosZEntry2/D");
	reflectTree->Branch("MomXEntry1",&MomXEntry1,"MomXEntry1/D");
	reflectTree->Branch("MomYEntry1",&MomYEntry1,"MomYEntry1/D");
	reflectTree->Branch("MomZEntry1",&MomZEntry1,"MomZEntry1/D");
	reflectTree->Branch("MomXEntry2",&MomXEntry2,"MomXEntry2/D");
	reflectTree->Branch("MomYEntry2",&MomYEntry2,"MomYEntry2/D");
	reflectTree->Branch("MomZEntry2",&MomZEntry2,"MomZEntry2/D");

	TVector3 aVec;
	int lastEvent=-1;
	for(int i=0;i < sbdFluxTree->GetEntries();i++)
	{
		sbdFluxTree->GetEntry(i);


		if(refEvent==event)
		{
			if(hitCode==1) // Entry
			{
				if(NumEntries==1) //Second entry
				{
					KEEntry2=KE;
					aVec.SetXYZ(posX,posY,posZ);
					aVec.RotateY(-9.5*TMath::Pi()/180.); //Rotate to bend coordinate system
					PosXEntry2=aVec.X();
					PosYEntry2=aVec.Y();
					PosZEntry2=aVec.Z();

					aVec.SetXYZ(momX,momY,momZ);
					aVec.RotateY(-9.5*TMath::Pi()/180.); //Rotate to bend coordinate system
					MomXEntry2=aVec.X();
					MomYEntry2=aVec.Y();
					MomZEntry2=aVec.Z();
				}
				NumEntries++;
			}
			else //Exit
			{
				if(NumExits==0) //First exit
				{
					KEExit1=KE;
				}
				NumExits++;
			}

		}
		else //New event
		{
			lastEvent=refEvent;
			//Fill previous event
			if(refEvent != -1)
			{
				reflectTree->Fill();
			}
			NumExits=0;
			NumEntries=0;
			KEEntry1=0;
			KEExit1=0;
			KEEntry2=0;
			PosXEntry1=0;
			PosYEntry1=0;
			PosZEntry1=0;
			PosXEntry2=0;
			PosYEntry2=0;
			PosZEntry2=0;
			MomXEntry1=0;
			MomYEntry1=0;
			MomZEntry1=0;
			MomXEntry2=0;
			MomYEntry2=0;
			MomZEntry2=0;

			//Any skipped events set to defaults
			for(int i=0; i < event-lastEvent-1; i++)
			{
				refEvent=lastEvent+i+1;
				reflectTree->Fill();
			}

			//1st entry
			refEvent=event;
			NumEntries=1;
			KEEntry1=KE;
			aVec.SetXYZ(posX,posY,posZ);
			aVec.RotateY(-9.5*TMath::Pi()/180.); //Rotate to bend coordinate system
			PosXEntry1=aVec.X();
			PosYEntry1=aVec.Y();
			PosZEntry1=aVec.Z();

			aVec.SetXYZ(momX,momY,momZ);
			aVec.RotateY(-9.5*TMath::Pi()/180.); //Rotate to bend coordinate system
			MomXEntry1=aVec.X();
			MomYEntry1=aVec.Y();
			MomZEntry1=aVec.Z();
		}

	}

	//Fill last detected event and skip remaining
	//Fill L event
	if(refEvent != -1)
	{
		reflectTree->Fill();
	}
	//This and any skipped events set to defaults

	lastEvent=refEvent;
	NumExits=0;
	NumEntries=0;
	KEEntry1=0;
	KEExit1=0;
	KEEntry2=0;
	PosXEntry1=0;
	PosYEntry1=0;
	PosZEntry1=0;
	PosXEntry2=0;
	PosYEntry2=0;
	PosZEntry2=0;
	MomXEntry1=0;
	MomYEntry1=0;
	MomZEntry1=0;
	MomXEntry2=0;
	MomYEntry2=0;
	MomZEntry2=0;
	for(int i=0; i < numEvents-lastEvent-1; i++)
	{
		refEvent=lastEvent+i+1;
		reflectTree->Fill();
	}



	reflectFile.cd();
	reflectTree->Write("",TObject::kOverwrite);

	resultsFile.Close();
	reflectFile.Close();


}

//reflectTree->Draw("PosXEntry1/1000.:PosYEntry1/1000.","KEEntry1 > 0")

void makePosXYPlotFromReflect(TString idString)
{
	//Get reflectTree
	TString reflectFilePath=TString(MRKGEANTDATA_DIR)+"Results_"+idString+"_Reflect_Set0.root";
	TFile reflectFile(reflectFilePath,"READ");
	TTree* reflectTree=(TTree*) reflectFile.Get("reflectTree");

	//Make reflect file to collect results of this analysis
	TString eventFilePath=TString(MRKGEANTDATA_DIR)+"Events_STD_3B_130429_Set0.root";
	TFile eventFile(eventFilePath,"READ");
	TTree* eventTree=(TTree*) eventFile.Get("t");
	reflectTree->AddFriend(eventTree,"e");

	TCanvas theCanvas("theCanvas", "theCanvas", 1200, 1200);

	int numFirstEntries = reflectTree->Draw("PosYEntry1*1000.:PosXEntry1*1000.","KEEntry1 > 0 && e.z0 > 0 && e.z0 < 0.32","goff");
	TGraph* firstEntryGraph = new TGraph(numFirstEntries,reflectTree->GetV1(),reflectTree->GetV2());
	firstEntryGraph->SetNameTitle("firstEntryGraph","Protons detected entering SBD Active Layer (Geant4.9.6.p02 with deadlayer/gold RDK2 Geo);x (mm);y (mm)");
	firstEntryGraph->Draw("AP");

	int numSecondEntries = reflectTree->Draw("PosYEntry2*1000.:PosXEntry2*1000.","KEEntry2 > 0 && e.z0 > 0 && e.z0 < 0.32","goff");
	TGraph* secondEntryGraph = new TGraph(numSecondEntries,reflectTree->GetV1(),reflectTree->GetV2());
	secondEntryGraph->SetNameTitle("firstEntryGraph","Protons detected entering SBD");
	secondEntryGraph->SetMarkerColor(kRed);
	secondEntryGraph->SetMarkerSize(1.5);
	secondEntryGraph->SetMarkerStyle(8);
	secondEntryGraph->Draw("P same");

	int numExits = reflectTree->Draw("KEExit1","KEExit1 > 0 && e.z0 > 0 && e.z0 < 0.32","goff");

	double percentExitsEntries=100.*numExits/(double) numFirstEntries;
	double percentSecondEntries=100.*numSecondEntries/(double) numFirstEntries;

	TLegend theLegend(0.65, .75, 0.90, 0.90);  //Top Right
	theLegend.AddEntry(firstEntryGraph,"First Entry Location","P");
	theLegend.AddEntry(secondEntryGraph,"Second Entry Location","P");
	theLegend.AddEntry(firstEntryGraph,"Percent exited: "+d2str(percentExitsEntries,2,1)+"%","");
	theLegend.AddEntry(firstEntryGraph,"Percent returned: "+d2str(percentSecondEntries,2,1)+"%","");
	theLegend.Draw("same");

	TString imagePath=TString(GRAPHS_DIR)+"Plot_LifeTimeBackscattering_PosX_v_PosY_1.png";
	theCanvas.SaveAs(imagePath);
}

