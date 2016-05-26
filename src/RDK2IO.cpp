#include "RDK2IO.h"

#include <iostream>
#include <sstream>

#include "TString.h"
#include "TChain.h"

#include "MRText.h"

#include "RDK2Constants.h"


using namespace std;

TString getNextEventFileName(TString inpFileName)
{

	//Example: MJB-Evn-hmg-n1.B3-147-20100421.root
	TString outString = fileNameFromFullPath(inpFileName);
	int fileNumber;
	stringstream theStream;
	int digits;

	for (int i = 0; i < outString.Length() - 3; i++)
	{
		if(outString[i] == '.')
		{
			if(outString[i + 1] == 'B')
			{
				if(outString[i + 5] == '-') //one digit number
				{
					fileNumber = (outString[i + 4] - '0');
					digits = 1;
				}
				else if(outString[i + 6] == '-') // two digit number
				{
					fileNumber = (outString[i + 5] - '0') + (outString[i + 4] - '0') * 10;
					digits = 2;
				}
				else
				{
					fileNumber = (outString[i + 6] - '0') + (outString[i + 5] - '0') * 10 + (outString[i + 4] - '0') * 100;
					digits = 3;
				}
				fileNumber++;

				theStream << fileNumber;

				outString.Replace(i + 4, digits, theStream.str());

				return filePathFromFullPath(inpFileName) + outString;

			}
		}
	}

	return (outString = "");
}

TString getNextResultFileName(TString inpFileName)
{
	TString outString = fileNameFromFullPath(inpFileName);
	int fileNumber;
	int digits;
	stringstream theStream;

	for (int i = 0; i < outString.Length() - 3; i++)
	{
		if(outString[i] == 'S' && outString[i + 1] == 'e' && outString[i + 2] == 't')
		{
			if(outString[i + 4] == '_' || outString[i + 4] == '.')
			{
				fileNumber = (outString[i + 3] - '0');
				fileNumber++;
			}
			else if(outString[i + 5] == '_' || outString[i + 5] == '.')
			{
				fileNumber = (outString[i + 4] - '0') + (outString[i + 3] - '0') * 10;
				fileNumber++;
			}
			else
			{
				fileNumber = (outString[i + 5] - '0') + (outString[i + 4] - '0') * 10 + (outString[i + 3] - '0') * 100;
				fileNumber++;
			}

			digits = 1;
			for (int fileCheck = 10;;)
			{
				if(fileNumber <= fileCheck) break;
				digits++;
				fileCheck *= 10;
			}

			theStream << fileNumber;

			outString.Replace(i + 3, digits, theStream.str());

			return filePathFromFullPath(inpFileName) + outString;

		}
	}

	return (outString = "");
}

TChain* makeEventChain(TString eventString,TString treeName,unsigned int numEventFiles, bool isEventGenerator)
{

    TChain* outChain=new TChain(treeName);
    unsigned int totalLoaded=numEventFiles;



    for(unsigned int i=0;i<numEventFiles;i++)
    {
        TString nextFileName=eventString;
        if(isEventGenerator)
        {
            nextFileName=string(EVENTS_DIR)+ "Events_"+eventString+"_Set"+int2str(i)+".root";
        }
        else if(i >=1)
        {
             getNextEventFileName(nextFileName);
        }

        if(!FileExists(nextFileName))
        {
            cout << "Could not load event file: " << nextFileName << endl;
            totalLoaded=i;
            break;

        }
        cout << "Adding to chain: " << nextFileName << endl;
        if(outChain->Add(nextFileName) <= 0)
        {
            totalLoaded=i;
            break;

        }

    }
    cout << "Total files loaded from file chain: " <<totalLoaded << " of " << numEventFiles << endl;

    if(totalLoaded < numEventFiles)
        typeAnythingToContinue("Less files found than requested.");

    return outChain;
}

TChain* makeLilBChain(TString eventString,TString treeName,unsigned int numEventFiles)
{

    TChain* outChain=new TChain(treeName);
    unsigned int totalLoaded=numEventFiles;



    for(unsigned int i=0;i<numEventFiles;i++)
    {
        TString lilBFilePath=string(LITTLEB_DIR)+ "LilB_"+eventString+"_Set"+int2str(i)+".root";


        if(!FileExists(lilBFilePath))
        {
            cout << "Could not load little b file: " << lilBFilePath << endl;
            totalLoaded=i;
            break;

        }
        cout << "Adding to little b chain: " << lilBFilePath << endl;
        if(outChain->Add(lilBFilePath) <= 0)
        {
            totalLoaded=i;
            break;

        }

    }
    cout << "Total files loaded from file chain: " <<totalLoaded << " of " << numEventFiles << endl;

    if(totalLoaded < numEventFiles)
        typeAnythingToContinue("Less files found than requested.");

    return outChain;
}


TChain* makeResultChain(TString firstFileName,TString treeName,unsigned int numEventFiles)
{
    TString nextFileName=firstFileName;
    TChain* outChain=new TChain(treeName);
    unsigned int totalLoaded=numEventFiles;

    for(unsigned int i=0;i<numEventFiles;i++)
    {
        if(!FileExists(nextFileName))
        {
            cout << "Could not load file: " << nextFileName << endl;
            totalLoaded=i;
            break;

        }
        cout << "Adding to chain: " << nextFileName << endl;
        if(outChain->Add(nextFileName) <= 0)
        {
            totalLoaded=i;
            break;

        }

        nextFileName=getNextResultFileName(nextFileName);

    }
    cout << "Total files loaded from file chain: " <<totalLoaded << " of " << numEventFiles << endl;
    if(totalLoaded < numEventFiles)
        typeAnythingToContinue("Less files found than requested.");

    return outChain;
}

TChain* makeExpChain(TString PIDString, const int numFiles,const int* seriesIDArray)
{

    TString expPath=EXPANALYSIS_DIR;
    TString runFilePath;

    TChain* outChain=new TChain("expTree");
    int totalLoaded=numFiles;


    for(int i=0;i<numFiles;i++)
    {

        runFilePath=expPath+"Processed_"+PIDString+"_"+"O_v8d10_"+int2str(seriesIDArray[i])+"++.itx.root";

        if(!FileExists(runFilePath))
        {
            cout << "Couldn't Load: " << runFilePath << endl;
            totalLoaded=i;
            break;

        }
        cout << "Adding to chain: " << runFilePath << endl;
        if(outChain->Add(runFilePath) <= 0)
        {
            totalLoaded=i;
            break;

        }


    }
    cout << "Total files loaded from file chain: " <<totalLoaded << " of " << numFiles << endl;

    if(totalLoaded < numFiles)
        typeAnythingToContinue("Less files found than requested.");

    return outChain;
}

void makeGeantMacroFilesForClusterEventGenerator(int startSet,int numSets,TString outResultTypeString,
                         ParType partype,TString eventString,TString addLine)
{
    TString partypeString,partypeLetter;
    TString timeString;
    int seedOffset=1;
    if(partype == PROTON)
    {
        partypeString="proton";
        partypeLetter="P";
        timeString="50:00:00";
//        timeString="70:00:00";
    }
    else if(partype == ELECTRON)
    {
        partypeString="e-";
        partypeLetter="E";
        seedOffset=1000;
        timeString="100:00:00";
    }
    else if(partype == PHOTON)
    {
        partypeString="gamma";
        partypeLetter="G";
        seedOffset=2000;
        timeString="8:00:00";
    }




    for(int set=startSet;set < numSets+startSet; set++)
    {

            ofstream outputFile;
            TString outputFileName;
            outputFileName=MACROS_DIR+outResultTypeString+"_"+int2str(set)+".mac";
            outputFile.open(outputFileName);

            outputFile << "/run/verbose 0" << endl;
            outputFile << "/event/verbose 0" << endl;
            outputFile << "/tracking/verbose 0" << endl;
            outputFile << "#" << endl;

            if(addLine != "")
                outputFile << addLine << endl;
            outputFile << "/process/eLoss/StepFunction .2 .1 mm" << endl;
            outputFile << "/process/eLoss/minKinEnergy 250 eV" << endl;
            outputFile << "/process/msc/RangeFactor .01" << endl;
            outputFile << "/process/msc/Skin 3" << endl;
            outputFile << "/MRK/phys/setCuts 1 nm" << endl;
            outputFile << "/cuts/setMaxCutEnergy 99 eV" << endl;
            outputFile << "/MRK/file/genMode 3" << endl;
            outputFile << "/MRK/field/update" << endl;
            outputFile << "/run/initialize" << endl;
            outputFile << "/gun/particle "<< partypeString << endl;
            outputFile << "/MRK/file/fluxmap Fluxmap_RDK2_Matt1.txt" << endl;
            outputFile << "/MRK/file/eventSetNumber "<<int2str(set)<<endl;
            outputFile << "/MRK/file/inpfile EventSettings_"+eventString+".root"<<endl;
            outputFile << "/MRK/file/rootout Results_"<<outResultTypeString<<"_"<< partypeLetter <<"_"+eventString+"_Set"<< set <<".root" << endl;
            outputFile << "/run/beamOn 1000000" << endl;
            outputFile.close();
    }

    ofstream subFile;
    TString subName=MACROS_DIR;
    subName+=outResultTypeString+".sub";

    subFile.open(subName);

    TString clusterIDString="g4_"+ outResultTypeString;
    TString macroFileString=outResultTypeString + "_${PBS_ARRAYID}.mac";

    subFile << "#!/bin/bash -l" << endl;
    subFile << "#PBS -A chupp_fluxod" << endl;
    subFile << "#PBS -l qos=flux" << endl;
    subFile << "#PBS -q fluxod" << endl;
    subFile << "#PBS -l qos=flux" << endl;
    subFile << "#PBS -N " << clusterIDString << endl;
    subFile << "#PBS -S /bin/bash" << endl;
    subFile << "#PBS -l procs=1,walltime="+timeString+",pmem=2000mb" << endl;
    subFile << "#PBS -t "<< startSet <<"-" <<startSet+numSets-1 << endl;
    subFile << "#PBS -M beam.that.slew@gmail.com" << endl;
    subFile << "#PBS -V" << endl;
    subFile << "#PBS -m abe" << endl;
    subFile << "#PBS -joe" << endl;
    subFile << "mkdir /tmp/${PBS_JOBID}" << endl;
    subFile << "cd /tmp/${PBS_JOBID}" << endl;
    subFile << "cp /home2/mjbales/geant/geant4/work/mrkgeant/build/4.9.6.p02/mrkgeant mrkgeant" << endl;
    subFile << "cp /home2/mjbales/exe/filesToAlwaysCopy/* ." << endl;
    subFile << "cp /home2/mjbales/macros/" << macroFileString << " ." << endl;
    subFile << "source /home2/mjbales/geant/geant4/geant4.9.6.p02-install/share/Geant4-9.6.2/geant4make/geant4make.sh" << endl;
    subFile << "module add lsa" << endl;
    subFile << "echo Running mrkgeant" << endl;
    subFile << "./mrkgeant " << macroFileString <<" $((${PBS_ARRAYID}+" << seedOffset <<"))" << endl;
    subFile << "#Cleanup" << endl;
    subFile << "cp Results* /home2/mjbales/exe" << endl;
    subFile << "cd" << endl;
    subFile << "#Retrieve output" << endl;
    subFile << "/bin/rm -rf /tmp/${PBS_JOBID}" << endl;

    subFile.close();

}


void makeGeantMacroFilesForClusterEventGeneratorMakeEventsOnly(int startSet,int numSets,TString eventString,TString addLine)
{
    TString timeString;

	timeString="00:10:00";


    for(int set=startSet;set < numSets+startSet; set++)
    {

            ofstream outputFile;
            TString outputFileName;
            outputFileName=TString(MACROS_DIR)+"makeEvents_"+int2str(set)+".mac";
            outputFile.open(outputFileName);

            outputFile << "/run/verbose 0" << endl;
            outputFile << "/event/verbose 0" << endl;
            outputFile << "/tracking/verbose 0" << endl;
            outputFile << "#" << endl;

            if(addLine != "")
                outputFile << addLine << endl;
            outputFile << "/MRK/file/genMode 3" << endl;
            outputFile << "/run/initialize" << endl;
            outputFile << "/MRK/file/fluxmap Fluxmap_RDK2_Matt1.txt" << endl;
            outputFile << "/MRK/file/eventSetNumber "<<int2str(set)<<endl;
            outputFile << "/MRK/file/inpfile EventSettings_"+eventString+".root"<<endl;
            outputFile.close();
    }

    ofstream subFile;
    TString subName=MACROS_DIR;
    subName+="makeEvents.sub";

    subFile.open(subName);

    TString clusterIDString="g4_makeEvents";
    TString macroFileString="makeEvents_${PBS_ARRAYID}.mac";

    subFile << "#!/bin/bash -l" << endl;
    subFile << "#PBS -A chupp_fluxod" << endl;
    subFile << "#PBS -l qos=flux" << endl;
    subFile << "#PBS -q fluxod" << endl;
    subFile << "#PBS -l qos=flux" << endl;
    subFile << "#PBS -N " << clusterIDString << endl;
    subFile << "#PBS -S /bin/bash" << endl;
    subFile << "#PBS -l procs=1,walltime="+timeString+",pmem=2000mb" << endl;
    subFile << "#PBS -t "<< startSet <<"-" <<startSet+numSets-1 << endl;
    subFile << "#PBS -M beam.that.slew@gmail.com" << endl;
    subFile << "#PBS -V" << endl;
    subFile << "#PBS -m abe" << endl;
    subFile << "#PBS -joe" << endl;
    subFile << "mkdir /tmp/${PBS_JOBID}" << endl;
    subFile << "cd /tmp/${PBS_JOBID}" << endl;
    subFile << "cp /home2/mjbales/geant/geant4/work/mrkgeant/build/4.9.6.p02/mrkgeant mrkgeant" << endl;
    subFile << "cp /home2/mjbales/exe/filesToAlwaysCopy/* ." << endl;
    subFile << "cp /home2/mjbales/macros/" << macroFileString << " ." << endl;
    subFile << "source /home2/mjbales/geant/geant4/geant4.9.6.p02-install/share/Geant4-9.6.2/geant4make/geant4make.sh" << endl;
    subFile << "module add lsa" << endl;
    subFile << "echo Running mrkgeant" << endl;
    subFile << "./mrkgeant " << macroFileString <<" $((${PBS_ARRAYID}+1))" << endl;
    subFile << "#Cleanup" << endl;
    subFile << "cp Events* /home2/mjbales/exe/events" << endl;
    subFile << "cd" << endl;
    subFile << "#Retrieve output" << endl;
    subFile << "/bin/rm -rf /tmp/${PBS_JOBID}" << endl;

    subFile.close();

}
