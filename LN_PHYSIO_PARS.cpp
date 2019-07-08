// PhysioParse.cpp : Program to parse Siemens Physiolog files files.
// The data is written into a tab seperated text file with a the name provided on the command line
// 
// This is larely taken from the idea discussion boars. Thus, I belive the fist version is from Peter Kochunov
// 
//  







//#include "stdafx.h" // Enable for Microsoft compilers, disable for linux and mac
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

enum PhysioMethod
{
	METHOD_NONE = 0x01,
	METHOD_TRIGGERING = 0x02,
	METHOD_GATING = 0x04,
	METHOD_RETROGATING = 0x08,
	METHOD_SOPE = 0x10,
	METHOD_ALL = 0x1E
};

enum ArrhythmiaDetection
{
	AD_NONE = 0x01,
	AD_TIMEBASED = 0x02,
	AD_PATTERNBASED = 0x04
};

enum PhysioSignal
{
	SIGNAL_NONE = 0x01,
	SIGNAL_EKG = 0x02,
	SIGNAL_PULSE = 0x04,
	SIGNAL_EXT = 0x08,
	SIGNAL_CARDIAC = 0x0E, /* the sequence usually takes this */
	SIGNAL_RESPIRATION = 0x10,
	SIGNAL_ALL = 0x1E,
};

#define DELTAT 0.0025f


int show_help( void )
{
   printf(
      "LN_PHYSIO_PARS: short program that pases SIEMENS physio logs\n"
      "\n"
      "    This program takes SIEMENS physio files ,\n"
      "    ECGlog_*.ecg\n"
      "    EXTlog_*.ext\n"
      "    Pulslog_*.puls\n"
      "    Resplog_*.resp\n"
      "    and parses them into txt files that \n"
      "    can be used in RETRIOCOR \n"
      "\n"
      "    basic usage: LN_PHYSIO_PARS input.puls output.txt  \n"
      "\n"
      "\n");
   return 0;
}



int main(int argc, char* argv[])
{
	if ( argc != 3 )
	{
		//cerr << "No file on command line\n";
		//cerr << "Provide the name of the input file (Siemens log file) and the output file (Tab separted log file)\n";
		
		
		show_help();

		return 1;
	}

	// get file from command line
	ifstream pfile(argv[1]);
	if (pfile.bad())
    {
        // Dump the contents of the file to cout.
        cout << "Cannot open file\n";
        pfile.close();
		return 1;
    }

	// read header
	int pMethod;
	pfile >> pMethod;
	cerr<<"Method ="<<pMethod<<endl;
	int ArrDect;
	pfile >> ArrDect;
	cerr<<"Detection of Arrytmia ="<<ArrDect<<endl;
	int SigSource;
	pfile >> SigSource;

	int GateOpen,GateClose;
	pfile >> GateOpen;
	pfile >> GateClose;
	cerr<<" Opening file "<<argv[1]<<endl; 
	// find out how long the file is. Parse to ECG in file at end of data
	long cnt=0;
	int input;
	do
	{
		pfile >> input;
	//	cerr<<cnt<<","<<input<<",";
		++cnt;
	}
	
	while ( input != 5003 || pfile.eof() );

	pfile.close();
	cerr<<" The file is "<<cnt<<" lines long\n";
	// now we know how long it is reopen and get data
	pfile.open(argv[1],ios::in);

	// skip header
	for ( int i=0 ; i<5 ; ++i )
		pfile >> input;

	// allocate for two interleved sets of data
	int n=cnt/2;
	int *pWform1,*pWform2;
	int *pTrigOn,*pTrigOff;
	int *pTrigCnt;
	float *pFreq;

	pWform1 = new int[n];
	pWform2 = new int[n];
	pTrigOn = new int[n];
	pTrigOff = new int[n];
	pTrigCnt = new int[n];
	pFreq = new float[n];

	// fill up the array
	int cnt2=0;
	int ntrig=0;
	int inval;
	for ( int i=0 ; i<n ; ++i )
	{
		//  zero the trigger on/off signals
		pTrigOn[i] = pTrigOff[i] = 0;
		// get values from file

		// get first data channel
		pfile >> inval;
		if ( inval == 5000 ) // trigger on
		{
			pTrigOn[i] = 500;
			pTrigCnt[ntrig++] = cnt2;
			pfile >> inval;	// get next data value
		}
		else if ( inval == 6000 ) // trigger off
		{
			pTrigOff[i] = 600;
			pfile >> inval;	// get next data value
		}
		else if ( inval == 5003 ) // end of data
		{
			break;
		}
		pWform1[i] = inval;

		// get second data channel
		pfile >> inval;
		if ( inval == 5000 ) // trigger on
		{
			pTrigOn[i] = 500;
			pTrigCnt[ntrig++] = cnt2;
			pfile >> inval;	// get next data value
		}
		else if ( inval == 6000 ) // trigger off
		{
			pTrigOff[i] = 600;
			pfile >> inval;	// get next data value
		}
		else if ( inval == 5003 ) // end of data
		{
			break;
		}
		pWform2[i] = inval;

		// subtract offsets
		pWform1[i] -= 10240; // it seems that the big values come first
		pWform2[i] -= 2048;  // if random then it could be checked for automatically

		++cnt2;
	}

	pfile.close();

	// calculate frequency from triggers

	// zero up to first trigger
	for ( int j=0 ; j<pTrigCnt[0] ; ++j )
		pFreq[j] = 0;

	// calculate frequency between trigger pulses
	for ( int i=0 ; i<ntrig-1 ; ++i )
	{
		float freq = 1.0f/((pTrigCnt[i+1] - pTrigCnt[i])*DELTAT);
		for ( int j=pTrigCnt[i] ; j<pTrigCnt[i+1] ; ++j )
			pFreq[j] = freq;
	}

	// frequency before first trigger same as first interval
	for ( int j=0 ; j<pTrigCnt[0] ; ++j )
		pFreq[j] = pFreq[pTrigCnt[0]];

	// frequency after last trigger same as last interval
	for ( int j=pTrigCnt[ntrig-1] ; j<cnt2 ; ++j )
		pFreq[j] = pFreq[pTrigCnt[ntrig-1]-1];

	// create an output file name from input file name
	string outname=argv[2];

	// write data to tab seperated file
	ofstream tfile(outname.c_str(),ios::out);

	// header
	tfile << "Time\tChan1\tChan2\tTrigOn\tTrigOff\tFreq\n";
	for ( int i=0 ; i<cnt2 ; ++i )
	{
		float t = i*DELTAT;
		tfile << t << '\t' << pWform1[i] << '\t' << pWform2[i] << '\t' << pTrigOn[i] <<
			'\t' << pTrigOff[i] << '\t' << pFreq[i] << endl;
	}

	tfile.close();

	return 0;
}

