/*
*    This file is part of Partanalyzer.
*
*    Partanalyzer is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Partanalyzer is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with Partanalyzer.  If not, see <http://www.gnu.org/licenses/>.
*
*/
/** Partanalyzer
Copyright (C) Miguel A. Santos, HSC, Toronto, 2008-2011.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )
*/


#ifndef _PARTANALYZER_MAIN
#define  _PARTANALYZER_MAIN 1
#endif

#include "partanalyzer.h"

//extern void printCommandLineError(const string label);
//extern void printCommandLineError(char* lastSeeOption);
//extern void exitWithMsg(const string);
//Reading from standard input
inline bool checkIfStandardInput(char* file1, char* file2="nofile", char* file3="nofile"){
			if(VERBOSE)cout<<"#Checking for standard input...";
			if(strcmp(file1,"-")==0&&strcmp(file2,"-")==0)printCommandLineError("Standard input can be assigned to only 1 input file");
			if(strcmp(file1,"-")==0&&strcmp(file3,"-")==0)printCommandLineError("Standard input can be assigned to only 1 input file");
			if(strcmp(file2,"-")==0&&strcmp(file3,"-")==0)printCommandLineError("Standard input can be assigned to only 1 input file");
			if(strcmp(file1,"-")==0){
				strcpy(file1,"/dev/stdin");
				if(!QUIET)cout<<"#Reading file1 from standard input "<<file1<<endl;
				return true;
			}
			else if(strcmp(file2,"-")==0){
				strcpy(file2,"/dev/stdin");
				if(!QUIET)cout<<"#Reading file2 from standard input "<<file2<<endl;
				return true;
			}
			else if(strcmp(file3,"-")==0){
				strcpy(file3,"/dev/stdin");
				if(!QUIET)cout<<"#Reading file3 from standard input "<<file3<<endl;
				return true;
			}
			if(VERBOSE)cout<<" NO"<<endl;
			return false;
}

//void readListInputFiles(ifstream is, vector<Charr> infilenames){
inline void readListInputFiles(char* argv0, vector<Charr>& infilenames){
	const string errmsg="readListInputFilesI: Standard input can be assigned to only 1 input file";
	bool usestdin=false;
	size_t inif=infilenames.size();
        ifstream is(argv0);
        if(!is){
                cout<<"ERROR: Cannot open file "<<argv0<<endl;
                exit(1);
        }
        string fn;
        Charr f;
        if(!QUIET)cout<<"#Getting list of partitions from "<<argv0<<endl;
        while(is>>fn){
                if(fn.compare(0,1,"#")==0){ // It's a comment line...
                        getline(is,fn);
                        continue;               // ignore the whole line and get to the next one.
                }
                f.car = new char[fn.size()+1];
                strcpy(f.car,fn.c_str());
	                if(checkIfStandardInput(f.car)){
				if(usestdin)printCommandLineError(errmsg);
				usestdin=true;
			}
                infilenames.push_back(f);
        }
        if(!QUIET)cout<<"#Partitions list file contains "<<infilenames.size()<<" entries.\n#First entry seen "<<infilenames[inif].car<<"\n#Last ("<<infilenames.size()<<") entry seen "<<infilenames[infilenames.size()-1].car<<endl;
}

inline void readListInputFiles(int& argc, char* argv[], vector<Charr>& infilenames, int nfiles=1000){
	const string errmsg="readListInputFilesI: Standard input can be assigned to only 1 input file";
	bool usestdin=false;
	size_t inif=infilenames.size();
	if(strcmp(*argv,"-f")==0){
		if(argc<2)printCommandLineError();
		argc--;argv++;
	        ifstream is(*argv);
	        if(!is){
	                cout<<"ERROR: Cannot open file "<<*argv<<endl;
	                exit(1);
	        }
	        string fn;
	        Charr f;
	        if(!QUIET)cout<<"#Getting list of partitions from "<<*argv<<endl;
	        while(is>>fn){
	                if(fn.compare(0,1,"#")==0){ // It's a comment line...
	                        getline(is,fn);
	                        continue;               // ignore the whole line and get to the next one.
	                }
	                f.car = new char[fn.size()+1];
	                strcpy(f.car,fn.c_str());
	                if(checkIfStandardInput(f.car)){
				if(usestdin)printCommandLineError(errmsg);
				usestdin=true;
			}
	                infilenames.push_back(f);
	        }
	}else {
		if(nfiles>argc)exitWithMsg("readListInputFilesI: Expect more files than arguments provided. Check code!");
		for(int i=0;i<nfiles;i++){
			Charr f={argv[0]};
			argc--;argv--;
	                if(checkIfStandardInput(f.car)){
				if(usestdin)printCommandLineError(errmsg);
				usestdin=true;
			}
			infilenames.push_back( f );
		}
	}
        if(!QUIET)cout<<"#Partitions list file contains "<<infilenames.size()<<" entries.\n#First entry seen "<<infilenames[inif].car<<"\n#Last ("<<infilenames.size()<<") entry seen "<<infilenames[infilenames.size()-1].car<<endl;
}


///Current name (full path) of Partanalyzer
char* program;
double beta;
double mu;

int main(int argc, char* argv[]) {
	program=argv[0];
	beta=BETA_UNSET;
	mu=MU_UNSET;
        if(argc<2) exitWithHelp();
	if(!QUIET&& !( strcmp(argv[1],"-q")==0 || strcmp(argv[1],"--help")==0 || strcmp(argv[1],"-help")==0 || strcmp(argv[1],"-h")==0 )){
		printCopyright();
		systemDate();
	}
	vector<Charr > infilenames;
	vector<Charr > infilenames2; //2nd list of input files for evaluating distances only between those in infilenames and infilenames2
	bool SETUPCONSENSUSP=false;
	bool MCLTABF=false;
	bool DIST_SUBSPROJECT=false;
	char* msaf=NULL;
	char* msafb=NULL;
	char* partitionf1=NULL;
	char* partitionf2=NULL;
	char* mcltabfile=NULL;
	char* mxofval=NULL;
	char* mxofvalb=NULL;
	string ClusterPrefix="C";
	double doublea;
	double doubleb;
	double extensivity=EXTENSIVITY_DEFAULT;
        double threshold=-1; //Default threshold
        double Id_threshold=50.0; //Default Seq. Id threshold
        double graph_pruning_thr; //Threshold for pruning graph edges
	int cluster1_offset=2; //Default cluster offset
	int cluster2_offset=cluster1_offset;
	int clstat_normalization_ofs=0;
	int Nsamples=1;
	int Nneighbors=1;
	int Nseq=1;
	int seed=-1;
	vector<int> positions;
	svect namelist;
	pmetricv metric=shannon;
	//int analysis=2;
	prganalysis analysis=prgVIPP;
	///Default partition input format: partanalyzer's own format.
	partFileFormat piformat=partFmtPART;
	///Default partition output format: partanalyzer's own format.
	partFileFormat poformat=partFmtPART;
	///Null MSA format
	MSAformat msafmt=msaFmtNULL;
	/* TESTING CUSTOM GRAPH CLASS
	graph_base<string, double> mygraph;
	string a,b;
	a="a";b="b"; 
	//mygraph.insert(pair<string,string>(a,b),0.0);
	mygraph[pair<string,string>(a,b)]=0.5;
	mygraph.insert(pair<string,string>(a,a),1.0);
	mygraph.insert(pair<string,string>(b,b),1.0);
	cout<<"GRAPH"<<endl;
	cout<<mygraph[pair<string,string>(a,a)]<<"\t"<<mygraph[pair<string,string>(a,b)]<<endl;
	cout<<mygraph[pair<string,string>(a,b)]<<"\t"<<mygraph[pair<string,string>(b,b)]<<endl;
	exit(0);
	*/
	if(argc>1){
		argc--;
		argv++;
///General Options
//We want these to be specified in any order
		int argcold=argc+1;
		while(argcold>argc){
			argcold=argc;
			if (strcmp(*argv,"--help")==0 || strcmp(*argv,"-help")==0){
				cout<<"\t\t\t"<<_programb_<<" (Partition Analyzer)\n"<<endl;
				printHelp();
				printHelpLong();
				cout<<"\nLicense"<<endl;
				printCopyright();
				exit(0);
			}
			if (strcmp(*argv,"-h")==0){
				printHelp();
				exit(0);
			}
			if(strcmp(*argv,"--version")==0){
				printCopyright();
				exit(0);
			}
			if(strcmp(*argv,"-q")==0||strcmp(*argv,"--quiet")==0){
				argc--;
				argv++;
				QUIET=true;
			}
			if(strcmp(*argv,"--verbose")==0){
				argc--;
				argv++;
				VERBOSE=true;
				cout<<"#VERBOSITY on"<<endl;
			}
			if(strcmp(*argv,"--debug")==0){
				argc--;
				argv++;
				DEBUG=true;
				VERBOSE=true;
				QUIET=false;
				cout<<"#EXTENSIVE VERBOSITY mode on"<<endl;
			}
			if(strcmp(*argv,"-z")==0||strcmp(*argv,"--pid-normalization")==0){
				argc--;argv++;
				if(argc<2)printCommandLineError();
				if(strcmp(*argv,"shorter-sequence")==0||strcmp(*argv,"s")==0){
					PIDNORMALIZATION=shorterSequence;
				}
				else if(strcmp(*argv,"aligned-positions")==0||strcmp(*argv,"p")==0){
					PIDNORMALIZATION=numberOfAlignedPositions;
				}
				else if(strcmp(*argv,"aligned-residues")==0||strcmp(*argv,"r")==0){
					PIDNORMALIZATION=numberOfAlignedResiduePairs;
				}
				else if(strcmp(*argv,"average-length")==0||strcmp(*argv,"l")==0){
					PIDNORMALIZATION=arithmeticMeanSequenceLenth;
				}
				else{
					cout<<"ERROR: Invalid PIDNORMALIZATION value : "<<*argv<<endl;
					exit(1);
				}
				argc--;argv++;
			}
			if(strcmp(*argv,"-t")==0||strcmp(*argv,"--format")==0||strcmp(*argv,"--fmt")==0){
				argc--;
				argv++;
				if(argc<2)printCommandLineError();
				if(strcmp(*argv,"MCL")==0||strcmp(*argv,"mcl")==0)
					piformat=partFmtMCL;
				else if(strcmp(*argv,"FREE")==0||strcmp(*argv,"free")==0)
					piformat=partFmtFREE;
				else if(strcmp(*argv,"PART")==0||strcmp(*argv,"part")==0)
					piformat=partFmtPART;
				else
					piformat=partFmtPART;
				argc--;
				argv++;
			}
			if(strcmp(*argv,"--oformat")==0||strcmp(*argv,"--ofmt")==0){
				argc--;
				argv++;
				if(argc<2)printCommandLineError();
				if(strcmp(*argv,"MCL")==0||strcmp(*argv,"mcl")==0)
					poformat=partFmtMCL;
				else if(strcmp(*argv,"FREE")==0||strcmp(*argv,"free")==0)
					poformat=partFmtFREE;
				else if(strcmp(*argv,"PART")==0||strcmp(*argv,"part")==0)
					poformat=partFmtPART;
				else
					poformat=partFmtPART;
				argc--;
				argv++;
			}
			if(strcmp(*argv,"--tab")==0){
				argc--;argv++;
				mcltabfile=argv[0]; //Its a tab file
				MCLTABF=true;
				argc--;argv++;
			}
			if(strcmp(*argv,"--DIST_SUBSPROJECT")==0){
				argc--;
				argv++;
				DIST_SUBSPROJECT=true;
			}
			if(strcmp(*argv,"--beta")==0){
				if(argc<5)printCommandLineError();
				argc--;
				argv++;
				beta=atof(argv[0]);
				argc--;
				argv++;
			}
			if(strcmp(*argv,"--mu")==0){
				if(argc<5)printCommandLineError();
				argc--;
				argv++;
				mu=atof(argv[0]);
				argc--;
				argv++;
			}
			if(strcmp(*argv,"--fuzzy")==0){
				argc--;
				argv++;
				FUZZYPARTITION=true;
				if(!QUIET)cout<<"#Partition type: Fuzzy"<<endl;
			}
		}
///(For analyzing partitions)
		///Calculate intra-cluster and inter-cluster distribution of weights
		if(strcmp(*argv,"-d")==0||strcmp(*argv,"--intra-inter-edge-dist")==0){
			analysis=prgCDIS;
			if(argc<3)printCommandLineError();
			mxofval=argv[1];
			partitionf1=argv[2];
			checkIfStandardInput(mxofval,partitionf1);
        		if(argc>3) cluster1_offset=atoi(argv[4]);
		}
		///Check consistency of partition according to weights in maxofval
		else if(strcmp(*argv,"-c")==0||strcmp(*argv,"--check-consistency-of-partition")==0||strcmp(*argv,"--ccop")==0){
			argc--;argv++;
			analysis=prgCCOP;
			if(argc<2)printCommandLineError();
			piformat=partFmtPART; //Enforcing partanalyzer's own input format. We don't want to deal with graphs in other formats yet.
			if(strcmp(*argv,"-tab")==0){
				if(argc<4)printCommandLineError();
				argc--;argv++;
				//partitionf2=argv[0]; //Its a tab file, but we use the available char* variable partitionf2.
				mcltabfile=argv[0]; //Its a tab file, but we use the available char* variable partitionf2.
				MCLTABF=true;
				argc--;argv++;
			}
			mxofval=argv[0];
			partitionf1=argv[1];
			checkIfStandardInput(mxofval,partitionf1);
			argc--;argv++;
			argc--;argv++;
        		if(argc>0) threshold=atof(argv[0]);
        		if(argc>1) cluster1_offset=atoi(argv[1]);
			
		}
		else if (strcmp(*argv,"-v")==0||strcmp(*argv,"-e")==0||strcmp(*argv,"-p")==0\
			||strcmp(*argv,"-i")==0||strcmp(*argv,"--intersection")==0\
			||strcmp(*argv,"-m")==0||strcmp(*argv,"--meet")==0\
			||strcmp(*argv,"-u")==0||strcmp(*argv,"--union")==0\
			||strcmp(*argv,"-j")==0||strcmp(*argv,"--join")==0\
			||strcmp(*argv,"--vi-distance")==0||strcmp(*argv,"--edit-distance")==0\
			||strcmp(*argv,"--purity-scores")==0){
			if(strcmp(*argv,"-v")==0||strcmp(*argv,"--vi-distance")==0)///analysis="vipp";
				analysis=prgVIPP;
			else if(strcmp(*argv,"-e")==0||strcmp(*argv,"--edit-distance")==0)
				analysis=prgEDSC;
			else if(strcmp(*argv,"-p")==0||strcmp(*argv,"--purity-scores")==0)
				analysis=prgPSPP;
			else if(strcmp(*argv,"-i")==0||strcmp(*argv,"--intersection")==0||strcmp(*argv,"-m")==0||strcmp(*argv,"--meet")==0)
				analysis=prgMEET;
			else if(strcmp(*argv,"-u")==0||strcmp(*argv,"--union")==0||strcmp(*argv,"-j")==0||strcmp(*argv,"--join")==0)
				analysis=prgJOIN;
			else {
				cout<<"ERROR: invalid binary option"<<endl;
				exit(1);
			}
			if(argc<3)printCommandLineError();
			partitionf1=argv[1];
			partitionf2=argv[2];
			checkIfStandardInput(partitionf1,partitionf2);
			/*if(strcmp(partitionf1,"-")==0&&strcmp(partitionf2,"-")==0)printCommandLineError("Standard input can be assigned to only 1 input file");
			if(strcmp(partitionf1,"-")==0){
				partitionf1="/dev/stdin";
			}
			else if(strcmp(partitionf2,"-")==0){
				partitionf2="/dev/stdin";
			}*/
			if(argc>3) cluster1_offset=atoi(argv[3]);
			if(argc>4) cluster2_offset=atoi(argv[4]);
		}
		else if (strcmp(*argv,"-n")==0||strcmp(*argv,"--ipot")==0||strcmp(*argv,"--iPot")==0\
			||strcmp(*argv,"--cpot")==0||strcmp(*argv,"--conditional-potential")==0\
			||strcmp(*argv,"--jpot")==0||strcmp(*argv,"--joint-potential")==0\
			||strcmp(*argv,"--mpot")==0||strcmp(*argv,"--mutual-potential")==0||strcmp(*argv,"--SA")==0||strcmp(*argv,"--subadditivity")==0\
			||strcmp(*argv,"--cmpot")==0||strcmp(*argv,"--conditional-mutual-potential")==0||strcmp(*argv,"--SSA")==0||strcmp(*argv,"--strong-subadditivity")==0\
			||strcmp(*argv,"--SSSA")==0||strcmp(*argv,"--soft-strong-subadditivity")==0||strcmp(*argv,"--SSA")==0\
			||strcmp(*argv,"--v-measure-a")==0||strcmp(*argv,"--v-measure-arithmetic")==0\
			||strcmp(*argv,"--v-measure-g")==0||strcmp(*argv,"--v-measure-geometric")==0\
			||strcmp(*argv,"--v-measure-h")==0||strcmp(*argv,"--v-measure-harmonic")==0){
			if(strcmp(*argv,"-n")==0||strcmp(*argv,"--ipot")==0||strcmp(*argv,"--iPot")==0)
				analysis=prgIPOT;
			else if(strcmp(*argv,"--cpot")==0||strcmp(*argv,"--conditional-potential")==0)
				analysis=prgCPOT;
			else if(strcmp(*argv,"--jpot")==0||strcmp(*argv,"--joint-potential")==0)
				analysis=prgJPOT;
			else if(strcmp(*argv,"--mpot")==0||strcmp(*argv,"--mutual-potential")==0||strcmp(*argv,"--SA")==0||strcmp(*argv,"--subadditivity")==0)
				analysis=prgMPOT;
			else if(strcmp(*argv,"--cmpot")==0||strcmp(*argv,"--conditional-mutual-potential")==0||strcmp(*argv,"--SSA")==0||strcmp(*argv,"--strong-subadditivity")==0)
				analysis=prgCMPO;
			else if(strcmp(*argv,"--SSSA")==0||strcmp(*argv,"--soft-strong-subadditivity")==0||strcmp(*argv,"--SSA")==0)
				analysis=prgSSSA;
			else if(strcmp(*argv,"--v-measure-a")==0||strcmp(*argv,"--v-measure-arithmetic")==0)
				analysis=prgVMAM;
			else if(strcmp(*argv,"--v-measure-g")==0||strcmp(*argv,"--v-measure-geometric")==0)
				analysis=prgVMGM;
			else if(strcmp(*argv,"--v-measure-h")==0||strcmp(*argv,"--v-measure-harmonic")==0){
				analysis=prgVMHM;
				if(beta==BETA_UNSET)beta=1.0;
				if(!QUIET)cout<<"#beta= "<<beta<<endl;
			}
			else {
				cout<<"ERROR: main : Sorry, I forgot what you said."<<endl;
				exit(1);
			}
			if(argc<3)printCommandLineError();
			argc--;argv++;
			if(strcmp(*argv,"v")==0||strcmp(*argv,"s")==0||strcmp(*argv,"vonneumann")==0||strcmp(*argv,"shannon")==0)
				metric=shannon;
			else if(strcmp(*argv,"b")==0||strcmp(*argv,"boltzmann")==0)
				metric=boltzmann;
			else if(strcmp(*argv,"e")==0||strcmp(*argv,"c")==0||strcmp(*argv,"cardinality")==0)
				metric=cardinality;
			else if(strcmp(*argv,"r")==0||strcmp(*argv,"renyi")==0)
				metric=renyi;
			else if(strcmp(*argv,"t")==0||strcmp(*argv,"tsallis")==0)
				metric=tsallis;
			else if(strcmp(*argv,"q")==0||strcmp(*argv,"tarantola")==0||strcmp(*argv,"jeffrey")==0||strcmp(*argv,"tjqn")==0){ //TARANTOLA-JEFFREY QNORM
				metric=jeffreyQnorm;
				extensivity=1.0;//DEFAULT VALUE FOR TARANTOLA QNORM
			}
			else {
				cout<<"ERROR: main : --xPot : unknown metric. Did you specified it?"<<endl;
				exit(1);
			}
			argc--;argv++;
			if(strcmp(*argv,"-e")==0||strcmp(*argv,"-ext")==0||strcmp(*argv,"--extensivity")==0){
				if(argc<3)printCommandLineError();
				argc--;argv++;
				extensivity=atof(argv[0]);
				if(argc<2)printCommandLineError();
				argc--;argv++;
			}
			if(strcmp(*argv,"-g")==0||strcmp(*argv,"--graph")==0||strcmp(*argv,"--matrix")==0){
				argc--;argv++;
				if(argc<2)printCommandLineError();
				mxofval=argv[0];
				argc--;argv++;
				if(argc<1)printCommandLineError();
				switch(analysis){
					case prgIPOT: analysis=prgIPOG; break;
					case prgCPOT: analysis=prgCPOG; break;
					case prgJPOT: analysis=prgJPOG; break;
					case prgMPOT: analysis=prgMPOG; break;
					case prgCMPO: analysis=prgCMPG; break;
					case prgSSSA: analysis=prgSSSG; break;
					default: 
						printCommandLineError("Weighted potential not available for chosen measure, only for: iPot, jPot, cPot, SA, SSA and SSSA"); 
						break;
				}
			}
			bool usestdin=false;
			const string errmsg="readListInputFilesI: Standard input can be assigned to only 1 input file";
			if(strcmp(*argv,"-ref")==0){
				switch(analysis){
					case prgIPOT: analysis=prgIPOR; break;
					case prgCPOT: analysis=prgCPOR; break;
					case prgJPOT: analysis=prgJPOR; break;
					case prgMPOT: analysis=prgMPOR; break;
					case prgCMPO: analysis=prgCMPR; break;
					case prgSSSA: analysis=prgSSSR; break;
					case prgVMAM: analysis=prgVAMR; break;
					case prgVMGM: analysis=prgVGMR; break;
					case prgVMHM: analysis=prgVHMR; break;
					case prgPSYM: analysis=prgPSYR; break;
					default: 
						printCommandLineError("Using -ref for weighted potentials not yet implemented"); 
						break;
				}
				if(argc<3)printCommandLineError();
				argc--;argv++;
				if(strcmp(*argv,"-f")!=0){
					Charr f={argv[0]};
			        	        if(checkIfStandardInput(f.car)){
							if(usestdin)printCommandLineError(errmsg);
							usestdin=true;
						}
					infilenames.push_back( f );
					argc--;argv++;
				}
			}
			if(strcmp(*argv,"-f")==0){
				if(argc<2)printCommandLineError();
				argc--;argv++;
				readListInputFiles(argv[0],infilenames);
			}
			else{
				for(int i=0;i<argc;i++){
					Charr f={argv[i]};
			                if(checkIfStandardInput(f.car)){
						if(usestdin)printCommandLineError(errmsg);
						usestdin=true;
					}
					infilenames.push_back( f );
				}
			}
		}
		else if (strcmp(*argv,"-J")==0||strcmp(*argv,"--Join")==0||strcmp(*argv,"-U")==0||strcmp(*argv,"--Union")==0||strcmp(*argv,"-M")==0||strcmp(*argv,"--Meet")==0||strcmp(*argv,"-I")==0||strcmp(*argv,"--Intersection")==0){
			if(strcmp(*argv,"-J")==0||strcmp(*argv,"--Join")==0||strcmp(*argv,"-U")==0||strcmp(*argv,"--Union")==0)
				analysis=prgJOST;
			else if(strcmp(*argv,"-M")==0||strcmp(*argv,"--Meet")==0||strcmp(*argv,"-I")==0||strcmp(*argv,"--Intersection")==0)
				analysis=prgMEST;
			bool usestdin=false;
			const string errmsg="readListInputFilesI: Standard input can be assigned to only 1 input file";
			if(argc<3)printCommandLineError();
			argc--;argv++;
			if(strcmp(*argv,"-f")==0){
				if(argc<2)printCommandLineError();
				argc--;argv++;
				if(!QUIET)cout<<"#Detected input file list"<<endl;
				readListInputFiles(argv[0],infilenames);
				argc--;argv++;
				if(argc>0&&strcmp(*argv,"-f2")==0){
					if(argc!=2)printCommandLineError();
					argc--;argv++;
					if(!QUIET)cout<<"#Detected second input file list"<<endl;
					readListInputFiles(argv[0],infilenames2);
					argc--;argv++;
				}
			}
			else{
				for(int i=0;i<argc;i++){
					Charr f={argv[i]};
			                if(checkIfStandardInput(f.car)){
						if(usestdin)printCommandLineError(errmsg);
						usestdin=true;
					}
					infilenames.push_back( f );
				}
			}
		}
		else if (strcmp(*argv,"-V")==0||strcmp(*argv,"-E")==0||strcmp(*argv,"-B")==0||strcmp(*argv,"-T")==0||strcmp(*argv,"-R")==0||strcmp(*argv,"-Q")==0\
			||strcmp(*argv,"--vstat")==0||strcmp(*argv,"--estat")==0||strcmp(*argv,"--bstat")==0||strcmp(*argv,"--tstat")==0||strcmp(*argv,"--rstat")==0||strcmp(*argv,"--qstat")==0){
			if(strcmp(*argv,"--vstat")==0||strcmp(*argv,"-V")==0)///analysis="vipp";
				analysis=prgVIST;
			else if(strcmp(*argv,"--bstat")==0||strcmp(*argv,"-B")==0)
				analysis=prgBDST;
			else if(strcmp(*argv,"--tstat")==0||strcmp(*argv,"-T")==0)
				analysis=prgTDST;
			else if(strcmp(*argv,"--rstat")==0||strcmp(*argv,"-R")==0)
				analysis=prgRDST;
			else if(strcmp(*argv,"--qstat")==0||strcmp(*argv,"-Q")==0){
				analysis=prgQDST;
				extensivity=1.0;//SET 1 AS DEFAULT FOR TARANTOLA-JEFFREY QNORM DISTANCE
			}
			else
				analysis=prgEDST;
			bool usestdin=false;
			const string errmsg="readListInputFilesI: Standard input can be assigned to only 1 input file";
			if(argc<3)printCommandLineError();
			argc--;argv++;
			if(strcmp(*argv,"-e")==0||strcmp(*argv,"-ext")==0||strcmp(*argv,"--extensivity")==0){
				if(argc<4)printCommandLineError();
				argc--;argv++;
				extensivity=atof(argv[0]);
				if(argc<3)printCommandLineError();
				argc--;argv++;
			}
			if(strcmp(*argv,"-ref")==0){
				switch(analysis){
					case prgVIST: analysis=prgVISR; break;
					case prgEDST: analysis=prgEDSR; break;
					case prgBDST: analysis=prgBDSR; break;
					case prgTDST: analysis=prgTDSR; break;
					case prgRDST: analysis=prgRDSR; break;
					case prgQDST: analysis=prgQDSR; break;
				}
				if(argc<3)printCommandLineError();
				argc--;argv++;
				if(strcmp(*argv,"-f")!=0){
					Charr f={argv[0]};
			        	        if(checkIfStandardInput(f.car)){
							if(usestdin)printCommandLineError(errmsg);
							usestdin=true;
						}
					infilenames.push_back( f );
					argc--;argv++;
				}
			}
			if(strcmp(*argv,"-f")==0){
				if(argc<2)printCommandLineError();
				argc--;argv++;
				if(!QUIET)cout<<"#Detected input file list"<<endl;
				readListInputFiles(argv[0],infilenames);
				argc--;argv++;
				if(argc>0&&strcmp(*argv,"-f2")==0){
					if(argc!=2)printCommandLineError();
					argc--;argv++;
					if(!QUIET)cout<<"#Detected second input file list"<<endl;
					readListInputFiles(argv[0],infilenames2);
					argc--;argv++;
				}
			}
			else{
				for(int i=0;i<argc;i++){
					Charr f={argv[i]};
			                if(checkIfStandardInput(f.car)){
						if(usestdin)printCommandLineError(errmsg);
						usestdin=true;
					}
					infilenames.push_back( f );
				}
			}
		}
		else if (strcmp(*argv,"-P")==0||strcmp(*argv,"--pstat")==0\
			||strcmp(*argv,"--pstat-sym")==0||strcmp(*argv,"--pstat-symmetric")==0){
			analysis=prgPSST;
			if(strcmp(*argv,"--pstat-sym")==0||strcmp(*argv,"--pstat-symmetric")==0)
				analysis=prgPSYM;
			if(argc<3)printCommandLineError();
			argc--;argv++;
			if(strcmp(*argv,"-ref")==0){
				if(analysis==prgPSYM) analysis=prgPSYR;
				else analysis=prgPSSR;
				if(argc<3)printCommandLineError();
				argc--;argv++;
				Charr f={argv[0]};
				infilenames.push_back( f );
				argc--;argv++;
			}
			if(strcmp(*argv,"-target")==0){
				if(analysis==prgPSYM) analysis=prgPSYR;
				else analysis=prgPSTG;
				if(argc<3)printCommandLineError();
				argc--;argv++;
				Charr f={argv[0]};
				infilenames.push_back( f );
				argc--;argv++;
			}
			if(strcmp(*argv,"-ofs")==0){
				if(argc<3)printCommandLineError();
				argc--;argv++;
				cluster1_offset=atoi(argv[0]);
				argc--;argv++;
			}
			if(strcmp(*argv,"-f")==0){
				if(argc<2)printCommandLineError();
				argc--;argv++;
				readListInputFiles(argv[0],infilenames);
				if(argc>0&&strcmp(*argv,"-f2")==0){
					if(argc!=2)printCommandLineError();
					argc--;argv++;
					readListInputFiles(argv[0],infilenames2);
				}
			}
			else
				for(int i=0;i<argc;i++){
					Charr f={argv[i]};
					infilenames.push_back( f );
				}
		}
		else if (strcmp(*argv,"-S")==0||strcmp(*argv,"--split-merge-analysis")==0||strcmp(*argv,"--splitstat")==0){
			analysis=prgSPSO;
			if(argc<3)printCommandLineError("--splitstat");
			argc--;argv++;
			if(strcmp(*argv,"-sim")==0){
				if(argc<3)printCommandLineError("--splitstat -sim");
				argc--;argv++;
				analysis=prgSPSS;
			}
			if(strcmp(*argv,"-over")==0||strcmp(*argv,"--overlap")==0||strcmp(*argv,"--fraction")==0){
				if(argc<3)printCommandLineError("--splitstat -over");
				argc--;argv++;
				analysis=prgSPSO;
			}
			if(strcmp(*argv,"-split")==0){
				if(argc<3)printCommandLineError("--splitstat -sim");
				argc--;argv++;
				analysis=prgSPST;
			}
			if(strcmp(*argv,"-ofs")==0){
				if(argc<3)printCommandLineError();
				argc--;argv++;
				cluster1_offset=atoi(argv[0]);
				argc--;argv++;
			}
			if(strcmp(*argv,"-f")==0){
				if(argc<2)printCommandLineError();
				argc--;argv++;
				readListInputFiles(argv[0],infilenames);
			}
			else
				for(int i=0;i<argc;i++){
					Charr f={argv[i]};
					infilenames.push_back( f );
				}
		}
		else if (strcmp(*argv,"-A")==0||strcmp(*argv,"--adjacency-stat")==0||strcmp(*argv,"--adjstat")==0){
			analysis=prgADST;
			if(argc<2)printCommandLineError();
			argc--;argv++;
			if(strcmp(*argv,"--print-fuzzy-partition")==0||strcmp(*argv,"-fuzzy")==0){
				if(argc<2)printCommandLineError();
				argc--;argv++;
				analysis=prgASFP;
			}
			if(strcmp(*argv,"-ofs")==0){
				if(argc<3)printCommandLineError();
				argc--;argv++;
				cluster1_offset=atoi(argv[0]);
				argc--;argv++;
			}
			if(strcmp(*argv,"-f")==0){
				if(argc<2)printCommandLineError();
				argc--;argv++;
				readListInputFiles(argv[0],infilenames);
			}
			else
				for(int i=0;i<argc;i++){
					Charr f={argv[i]};
					infilenames.push_back( f );
				}
		}
		else if (strcmp(*argv,"-C")==0||strcmp(*argv,"--cluster-stat")==0||strcmp(*argv,"--clstat")==0){
			analysis=prgCLST;
			if(argc<2)printCommandLineError();
			argc--;argv++;
			if(strcmp(*argv,"-ofs")==0){
				if(argc<3)printCommandLineError();
				argc--;argv++;
				cluster1_offset=atoi(argv[0]);
				argc--;argv++;
			}
			if(strcmp(*argv,"-norm")==0){
				if(argc<3)printCommandLineError();
				argc--;argv++;
				clstat_normalization_ofs=atoi(argv[0]);
				argc--;argv++;
			}
			if(strcmp(*argv,"-cons")==0||strcmp(*argv,"-consensus")==0){
				if(argc<2)printCommandLineError();
				argc--;argv++;
				SETUPCONSENSUSP=true;
			}
			if(strcmp(*argv,"-f")==0){
				if(argc<2)printCommandLineError();
				argc--;argv++;
				readListInputFiles(argv[0],infilenames);
				//if(VERBOSE)cout<<"#Option -f read "<<infilenames.size()<<" partitions"<<endl;
			}
			else
				for(int i=0;i<argc;i++){
					Charr f={argv[i]};
					infilenames.push_back( f );
				}
		}
		else if (strcmp(*argv,"--Info")==0||strcmp(*argv,"--info")==0||strcmp(*argv,"--isPart")==0||strcmp(*argv,"--is-partition")==0||strcmp(*argv,"--isaPart")==0){
			analysis=prgIPAR;
			INFO=true;
			if(argc<2)printCommandLineError();
			argc--;argv++;
			if(strcmp(*argv,"-ofs")==0){
				if(argc<3)printCommandLineError();
				argc--;argv++;
				cluster1_offset=atoi(argv[0]);
				argc--;argv++;
			}
			if(strcmp(*argv,"-f")==0){
				if(argc<2)printCommandLineError();
				argc--;argv++;
				readListInputFiles(argv[0],infilenames);
			}
			else
				for(int i=0;i<argc;i++){
					Charr f={argv[i]};
					infilenames.push_back( f );
				}
		}
		else if (strcmp(*argv,"-H")==0||strcmp(*argv,"-hasse")==0){
			analysis=prgHASS;
			if(argc<2)printCommandLineError();
			argc--;argv++;
			if(strcmp(*argv,"-ofs")==0){
				if(argc<3)printCommandLineError();
				argc--;argv++;
				cluster1_offset=atoi(argv[0]);
				argc--;argv++;
			}
			if(strcmp(*argv,"-f")==0){
				if(!QUIET)cout<<"#Seen option -f"<<endl;
				if(argc<2)printCommandLineError();
				argc--;argv++;
				readListInputFiles(argv[0],infilenames);
			}
			else
				for(int i=0;i<argc;i++){
					Charr f={argv[i]};
					infilenames.push_back( f );
				}
		}
///(For creating partitions)
		else if(strcmp(*argv,"--cluster-robust")==0||strcmp(*argv,"--robust-cluster")==0||strcmp(*argv,"--cluster-robust-self-consistent")==0||strcmp(*argv,"--RDC")==0){
			analysis=prgRPAR;
			if(strcmp(*argv,"--cluster-robust-self-consistent")==0||strcmp(*argv,"--RDC")==0){
				analysis=prgRPSC;
			}
			if(argc<2)printCommandLineError("Missing graph");
			argc--;argv++;
			mxofval=argv[0];
			argc--;argv++;
			Nsamples=RDC_DEFAULT_TOP_BEST_PARTITIONS;
			Nneighbors=RDC_DEFAULT_NUMBER_NEIGHBORS;
			while(argc>0){
				//if(strcmp(*argv,"below")==0||strcmp(*argv,"above")==0){
				if(strcmp(*argv,"-below")==0){
					if(analysis==prgRPSC){
						analysis=prgRDCB;
					}else{
						analysis=prgRPAB;
					}
					argc--;argv++;
				}
				else if(strcmp(*argv,"-above")==0){
					if(analysis==prgRPSC){
						analysis=prgRDCA;
					}else{
						analysis=prgRPAA;
					}
					argc--;argv++;
				}
				else if(strcmp(*argv,"-s")==0){
					if(argc<2)printCommandLineError("Missing NumberOfSample");
					argc--;argv++;
					Nsamples=atoi(argv[0]);
					argc--;argv++;
				}
				else if(strcmp(*argv,"-n")==0){
					if(argc<2)printCommandLineError("Missing NumberOfNeighbors");
					argc--;argv++;
					Nneighbors=atoi(argv[0]);
					argc--;argv++;
				}
				else if(strcmp(*argv,"-ext")==0||strcmp(*argv,"--extensivity")==0){
					if(argc<2)printCommandLineError("Missing NumberOfSample");
					argc--;argv++;
					extensivity=atof(argv[0]);
					argc--;argv++;
				}
				else if(strcmp(*argv,"-v")==0||strcmp(*argv,"-V")==0){
					metric=shannon;
					argc--;argv++;
				}
				else if(strcmp(*argv,"-e")==0||strcmp(*argv,"-E")==0){
					metric=cardinality;
					argc--;argv++;
				}
				else if(strcmp(*argv,"-r")==0||strcmp(*argv,"-R")==0){
					metric=renyi;
					argc--;argv++;
				}
				else if(strcmp(*argv,"-t")==0||strcmp(*argv,"-T")==0){
					metric=tsallis;
					argc--;argv++;
				}
				else if(strcmp(*argv,"-q")==0||strcmp(*argv,"-Q")==0){
					metric=jeffreyQnorm;
					argc--;argv++;
				}
			}
		}
		else if(strcmp(*argv,"--cluster")==0){
			analysis=prgPART;
			if(argc<2)printCommandLineError("Missing graph");
			argc--;argv++;
			if(strcmp(*argv,"-below")==0||strcmp(*argv,"-above")==0)printCommandLineError("File name must be first argument");
			mxofval=argv[0];
			argc--;argv++;
			if(argc>0){
				if(strcmp(*argv,"-above")==0){
					analysis=prgPARA;
					if(argc<2)printCommandLineError("Missing threshold");
					argc--;argv++;
				}else{
					analysis=prgPARB;
					if(strcmp(*argv,"-below")==0){
						if(argc<2)printCommandLineError("Missing threshold");
						argc--;argv++;
					}
				}
				graph_pruning_thr=atof(argv[0]);
				argc--;argv++;
			}
		}
///(For editing partitions)
		else if (strcmp(*argv,"--part-swap-names")==0||strcmp(*argv,"--part-swap-labels")==0){
			analysis=prgPACN;
			if(!MCLTABF)printCommandLineError("Use option --tab tabfile before command to specify mapping for new names");
			if(argc<2)printCommandLineError("Missing partition");
			argc--;argv++;
			partitionf1=argv[0];
			argc--;argv++;
		}
		else if (strcmp(*argv,"--part-sort")==0||strcmp(*argv,"--part-sort-rename")==0){
			analysis=prgPASO;
			if(strcmp(*argv,"--part-sort-rename")==0) analysis=prgPASR;
			if(argc<2)printCommandLineError();
			argc--;argv++;
			partitionf1=argv[0];
			argc--;argv++;
			if(argc>0 && analysis==prgPASR){
				ClusterPrefix=argv[0];
				argc--;argv++;
			}
		}
		else if (strcmp(*argv,"--part-extract-elements")==0||strcmp(*argv,"--extract-elements")==0){
			analysis=prgPAXE;
			if(argc<3)printCommandLineError();
			argc--;argv++;
			if(!QUIET)cout<<"#Extracting elements from partition"<<endl;
			readListFromFile(argv[0],namelist);
			argc--;argv++;
			if(argc>0&&strcmp(*argv,"-tab")==0){
				if(argc<2)printCommandLineError();
				argc--;argv++;
				//partitionf2=argv[0]; //Its a tab file, but we use the available char* variable partitionf2.
				mcltabfile=argv[0]; 
				MCLTABF=true;
				argc--;argv++;
			}
			partitionf1=argv[0];
			argc--;argv++;
        		if(argc>0) cluster1_offset=atoi(argv[0]);
		}
///(For converting between different partition formats)
		else if (strcmp(*argv,"--toMCL")==0||strcmp(*argv,"--toFREE")==0||strcmp(*argv,"--MCLtoPART")==0 ||\
			strcmp(*argv,"--MCL2PART")==0||\
			strcmp(*argv,"--MSAtoPART")==0 ||strcmp(*argv,"--MSA2PART")==0 ||\
			 strcmp(*argv,"--toPART")==0||strcmp(*argv,"--FREEtoPART")==0){
			if(strcmp(*argv,"--toMCL")==0){
				analysis=prg2MCL;
				poformat=partFmtMCL;
			}
			else if(strcmp(*argv,"--toFREE")==0){
				analysis=prg2FRE;
				poformat=partFmtFREE;
			}
			else if(strcmp(*argv,"--toPART")==0){
				analysis=prg2FRE;
				poformat=partFmtPART;
			}
			else if(strcmp(*argv,"--MCLtoPART")==0||strcmp(*argv,"--MCL2PART")==0){
				analysis=prgM2PA;
				piformat=partFmtMCL;
				poformat=partFmtPART;
			}
			else if(strcmp(*argv,"--FREEtoPART")==0){
				analysis=prgM2PA;
				piformat=partFmtFREE;
				poformat=partFmtPART;
			}
			else if(strcmp(*argv,"--MSAtoPART")==0||strcmp(*argv,"--MSA2PART")==0){
				analysis=prgA2PA;
				piformat=partFmtMSA;
				poformat=partFmtPART;
			}
			else {
				cout<<"ERROR: unknown conversion option"<<endl;
				exit(1);
			}
			argc--;argv++;
			if(argc<1)printCommandLineError();
			if( (analysis==prgM2PA || analysis==prg2MCL )&& strcmp(*argv,"-tab")==0){
				argc--;argv++;
				if(argc<2)printCommandLineError();
				//partitionf2=argv[0]; //Its a MCL tab file, but we use the available char* variable partitionf2.
				mcltabfile=argv[0]; //Its a tab file, but we use the available char* variable partitionf2.
				MCLTABF=true;
				argc--;argv++;
			}
			partitionf1=argv[0];
		}
///(For dealing with (fasta) sequence files
                else if (strcmp(*argv,"--drop-clone-sequences")==0||strcmp(*argv,"--seq-noclone-sequences")==0||strcmp(*argv,"--msa-noclone-sequences")==0){
                        analysis=prgSNOC;
                        if(argc<2)printCommandLineError();
                        argc--;argv++;
                        msaf=argv[0];
                        argc--;argv++;
			if(argc>0){
				msafb=argv[0];
                        	argc--;argv++;
			}
                }
///(For analyzing Multiple Sequence Alignments)
                else if (strcmp(*argv,"--msa-print")==0||strcmp(*argv,"--print-msa")==0){
                        analysis=prgPMSA;
                        if(argc<2)printCommandLineError();
                        argc--;argv++;
			while(argc>0){
				if(strcmp(*argv,"-nosort")==0){
                        		argc--;argv++;
                        		if(argc<1)printCommandLineError("Missing MSA filename");
				}
				else if(strcmp(*argv,"-sort")==0){
                        		argc--;argv++;
                        		if(argc<1)printCommandLineError("Missing MSA filename");
                        		analysis=prgPMAS;
				}
				else {
                        		msaf=argv[0];
                        		argc--;argv++;
				}
			}
                        //msaf=argv[0];
                        //argc--;argv++;
                }
                else if (strcmp(*argv,"--msa-map-partition")==0){
                        analysis=prgMSMP;
                        if(argc<3)printCommandLineError();
                        argc--;argv++;
			partitionf1=argv[0];
                        argc--;argv++;
                        msaf=argv[0];
                        argc--;argv++;
			while(argc>0){
				if(strcmp(*argv,"FASTA")==0||strcmp(*argv,"fasta")==0){
					if(!QUIET)cout<<"#MSA output format: FASTA"<<endl;
					msafmt=FASTA;
				}
				else if(strcmp(*argv,"FASTA2")==0||strcmp(*argv,"fasta2")==0){
					if(!QUIET)cout<<"#MSA output format: FASTA2"<<endl;
					msafmt=FASTA2;
				}
				else if(strcmp(*argv,"FASTA3")==0||strcmp(*argv,"fasta3")==0){
					if(!QUIET)cout<<"#MSA output format: FASTA3"<<endl;
					msafmt=FASTA3;
				}
				else if(strcmp(*argv,"GDE")==0||strcmp(*argv,"gde")==0){
					if(!QUIET)cout<<"#MSA output format: GDE"<<endl;
					msafmt=GDE;
				}
				else if(strcmp(*argv,"GDE2")==0||strcmp(*argv,"gde2")==0){
					if(!QUIET)cout<<"#MSA output format: GDE2"<<endl;
					msafmt=GDE2;
				}
				else if(strcmp(*argv,"GDE3")==0||strcmp(*argv,"gde3")==0){
					if(!QUIET)cout<<"#MSA output format: GDE3"<<endl;
					msafmt=GDE3;
				}
				else if(strcmp(*argv,"SPEER")==0||strcmp(*argv,"speer")==0){
					if(!QUIET)cout<<"#MSA output format: SPEER"<<endl;
					msafmt=SPEER;
				}
				else if(strcmp(*argv,"SPEER2")==0||strcmp(*argv,"speer2")==0){
					if(!QUIET)cout<<"#MSA output format: SPEER2"<<endl;
					msafmt=SPEER2;
				}
				else if(strcmp(*argv,"SPEER3")==0||strcmp(*argv,"speer3")==0){
					if(!QUIET)cout<<"#MSA output format: SPEER3"<<endl;
					msafmt=SPEER3;
				}
				else if(strcmp(*argv,"GSIM")==0||strcmp(*argv,"gsim")==0){
					if(!QUIET)cout<<"#MSA output format: GSIM"<<endl;
					msafmt=GSIM;
				}
				else if(strcmp(*argv,"GSIM2")==0||strcmp(*argv,"gsim2")==0){
					if(!QUIET)cout<<"#MSA output format: GSIM2"<<endl;
					msafmt=GSIM2;
				}
				else if(strcmp(*argv,"GSIM3")==0||strcmp(*argv,"gsim3")==0){
					if(!QUIET)cout<<"#MSA output format: GSIM3"<<endl;
					msafmt=GSIM3;
				}
				else{
					if(!QUIET)cout<<"#WARNING : Main : --msa-map-partition : Unknown MSA format : Ignoring option"<<endl;
					msafmt=msaFmtNULL;
				}
                        	argc--;argv++;
			}
                }
                else if (strcmp(*argv,"--msa-seqid-stat")==0){
                        analysis=prgMSPI;
                        if(argc<2)printCommandLineError();
                        argc--;argv++;
                        if(strcmp(*argv,"--positions")==0){
                                if(argc<3)printCommandLineError();
                                argc--;argv++;
                                ifstream is(argv[0]);
                                if(!is){
                                        cout<<"ERROR: Cannot open file "<<argv[0]<<endl;
                                        exit(1);
                                }
                                argc--;argv++;
                                int pos;
                                if(!QUIET)cout<<"#Using only specific MSA columns (positions):\n#";
                                while(is>>pos){
                                        if(!QUIET)cout<<"\t"<<pos;
                                        positions.push_back(pos);
                                }
                                if(!QUIET)cout<<endl;
                        }
                        msaf=argv[0];
			argc--;argv++;
                        if(argc>0){
				msafb=argv[0];
				analysis=prgMMPI;
			}
                }
                else if (strcmp(*argv,"--msa-seqid-avg")==0){
                        analysis=prgMAPI;
                        if(argc<2)printCommandLineError();
			argc--;argv++;
			if(strcmp(*argv,"-thr")==0){
                        	if(argc<3)printCommandLineError();
				argc--;argv++;
				Id_threshold=atof(argv[0]);
				argc--;argv++;
			}
                        msaf=argv[0];
			argc--;argv++;
                        if(argc>0){
				msafb=argv[0];
				analysis=prgMMAI;
			}
                }
                else if (strcmp(*argv,"--msa-extract-positions")==0){
                        analysis=prgMSXP;
                        if(argc<3)printCommandLineError();
                        argc--;argv++;
                        ifstream is(argv[0]);
                        if(!is){
                                cout<<"ERROR: Cannot open file "<<argv[0]<<endl;
                                exit(1);
                        }
                        argc--;argv++;
                        int pos;
                        if(!QUIET)cout<<"#Using only specific MSA columns (positions):\n#";
                        while(is>>pos){
                                if(!QUIET)cout<<"\t"<<pos;
				//User gives positions starting at 1; Sequence starts them at 0
                                positions.push_back(--pos);
                        }
                        if(!QUIET)cout<<endl;
                        if(argc<1)printCommandLineError("Missing MSA file");
                        msaf=argv[0];
		}
		else if(strcmp(*argv,"--msa-extract-sequences-by-topid")==0||strcmp(*argv,"--msa-drop-sequences-by-topid")==0){
                        analysis=prgMSXT;
			string Upp="Extract";
			// DROP not implemented yet
			if(strcmp(*argv,"--msa-drop-sequences-by-topid")==0){
				if(!QUIET)cout<<"#WARNING: Not implemented yet : --msa-drop-sequences-by-topid. Proceeding with --msa-extract-sequences-by-topid"<<endl;
				/*
				analysis=prgMSDT;
				Upp="Drop";
				*/
			}
                        if(argc<3)printCommandLineError("--msa(extract/drop)-sequences-by-topid");
                        argc--;argv++;
			msaf=argv[0];
                        argc--;argv++;
			msafb=argv[0];	//We want to extract/drop sequences such that their Id's against msafb
                        argc--;argv++;
                        Nsamples=0;
			if(argc>0){
				Nsamples=atoi(argv[0]);
                        	argc--;argv++;
				if(argc>0){ 
                        		ifstream is(argv[0]);
                        		if(!is){
                        		        cout<<"ERROR: Cannot open file "<<argv[0]<<endl;
                        		        exit(1);
                        		}
                        		argc--;argv++;
                        		int pos;
                        		if(!QUIET)cout<<"#Using only specific MSA columns (positions):\n#";
                        		while(is>>pos){
                        		        if(!QUIET)cout<<"\t"<<pos;
						//User gives positions starting at 1; Sequence starts them at 0
                        		        positions.push_back(--pos);
                        		}
                        		if(!QUIET)cout<<endl;
				}
			}
		}
                else if (strcmp(*argv,"--msa-extract-sequences-by-id")==0||strcmp(*argv,"--msa-drop-sequences-by-id")==0){
                        analysis=prgMSXI;
			string Upp="Extract";
			if(strcmp(*argv,"--msa-drop-sequences-by-id")==0){
				analysis=prgMSDI;
				Upp="Drop";
			}
                        if(argc<5)printCommandLineError("--msa(extract/drop)-sequences-by-id");
                        argc--;argv++;
			msaf=argv[0];
                        argc--;argv++;
			msafb=argv[0];	//We want to extract/drop sequences such that their Id's against msafb
                        argc--;argv++;
			doublea=30; //B
			doubleb=100;
			if(argc>0&&argc<2)printCommandLineError("--msa(extract/drop)-sequences-by-id");
			doublea=atof(argv[0]); // is doublea < Id <= doubleb
                        argc--;argv++;
			doublea=atof(argv[0]);
                        argc--;argv++;
			if(argc>0){
                        	ifstream is(argv[0]);
                        	if(!is){
                        	        cout<<"ERROR: Cannot open file "<<argv[0]<<endl;
                        	        exit(1);
                        	}
                        	argc--;argv++;
                        	int pos;
                        	if(!QUIET)cout<<"#Using only specific MSA columns (positions):\n#";
                        	while(is>>pos){
                        	        if(!QUIET)cout<<"\t"<<pos;
					//User gives positions starting at 1; Sequence starts them at 0
                        	        positions.push_back(--pos);
                        	}
                        	if(!QUIET)cout<<endl;
			}
		}
                else if (strcmp(*argv,"--msa-extract-sequences")==0||strcmp(*argv,"--msa-drop-sequences")==0){
                        analysis=prgMSXS;
			string Upp="Extract";
			if(strcmp(*argv,"--msa-drop-sequences")==0){
				analysis=prgMSDS;
				Upp="Drop";
			}
                        if(argc<3)printCommandLineError();
                        argc--;argv++;
                        ifstream is(argv[0]);
                        if(!is){
                                cout<<"ERROR: Cannot open file "<<argv[0]<<endl;
                                exit(1);
                        }
                        argc--;argv++;
			string st,stb;
			while(is>>st){ //Sequence names expected in a file
				namelist.push_back(st);
			}
			if(!QUIET)cout<<"#"<<Upp<<"ing Sequences from multiple sequence alignment\n#Number of sequences to "<<Upp<<" "<<namelist.size()<<"\n#Last name found : "<<namelist[namelist.size()-1]<<endl;
			msaf=argv[0];
                        argc--;argv++;
			if(argc>0&&argc!=3)printCommandLineError("--msa(extract/drop)-sequences");
			msafb=argv[0];	//We want to extract/drop sequences such that their Id's against msafb
                        argc--;argv++;
			doublea=atof(argv[0]); // is doublea < Id <= doubleb
                        argc--;argv++;
			doublea=atof(argv[0]);
                        argc--;argv++;
		}
		else if (strcmp(*argv,"--msa-redundant")==0){
			analysis=prgMRED;
			if(argc<2)printCommandLineError();
			argc--;argv++;
			if(strcmp(*argv,"-nsam")==0){
				if(argc<3)printCommandLineError();
				argc--;argv++;
				Nsamples=atoi(argv[0]);
				argc--;argv++;
			}
			if(strcmp(*argv,"-nseq")==0){
				if(argc<3)printCommandLineError();
				argc--;argv++;
				Nseq=atoi(argv[0]);
				argc--;argv++;
			}
			if(strcmp(*argv,"-seed")==0){
				if(argc<3)printCommandLineError();
				argc--;argv++;
				seed=atoi(argv[0]);
				if(seed<0)seed*=-1;
				argc--;argv++;
			}
                        msaf=argv[0];
			argc--;argv++;
		}
///(For dealing with -interaction- matrices) 
		else if (strcmp(*argv,"--edge-dist")==0){
			analysis=prgEDMX;
			if(argc<2)printCommandLineError();
			mxofval=argv[1];
			if(argc>2){
				analysis=prgEDMP;
				partitionf1=argv[2];
			}
			if(argc>3)cluster1_offset=atoi(argv[4]);
		}
		else if (strcmp(*argv,"-m")==0||strcmp(*argv,"--merge-graphs")==0){
			///analysis="merge matrices";
			analysis=prgMGMX;
			if(argc<3)printCommandLineError();
			mxofval=argv[1];
			mxofvalb=argv[2];
		}
		else if (strcmp(*argv,"-r")==0||strcmp(*argv,"--merge-graphs-color")==0){
			///analysis="merge matrices with color";
			analysis=prgMMXC;
			//piformat=partFmtPART; //Enforcing partanalyzer's own input format. We don't want to deal with graphs in other formats yet.
			if(argc<4)printCommandLineError();
			mxofval=argv[1];
			mxofvalb=argv[2];
			partitionf1=argv[3];
			if(argc>4) cluster1_offset=atoi(argv[4]);
		}
		else if (strcmp(*argv,"-l")==0 || strcmp(*argv,"--cull-edges")==0){
			analysis=prgCEMX;
			if(argc<3)printCommandLineError();
			mxofval=argv[1];
			mxofvalb=argv[2];
		}
		else if(strcmp(*argv,"--prune-edges-above")==0){
			analysis=prgGPRA;
			if(argc<3)printCommandLineError("Missing arguments");
			argc--;argv++;
			graph_pruning_thr=atof(argv[0]);
			argc--;argv++;
			mxofval=argv[0];
			argc--;argv++;
		}
		else if(strcmp(*argv,"--prune-edges-below")==0){
			analysis=prgGPRB;
			if(argc<3)printCommandLineError("Missing arguments");
			argc--;argv++;
			graph_pruning_thr=atof(argv[0]);
			argc--;argv++;
			mxofval=argv[0];
			argc--;argv++;
		}
		else if(strcmp(*argv,"--print-matrix")==0||strcmp(*argv,"--graph-print")==0||strcmp(*argv,"--matrix-print")==0){
			argc--;argv++;
			int col=EDGES_DEFAULT_COLUMN;
			//if(argc<1)printCommandLineError("Missing graph file");
			mxofval="/dev/stdin";
			while(argc>0){
				if(strcmp(*argv,"-c")==0||strcmp(*argv,"-col")==0||strcmp(*argv,"-column")==0){
					argc--;argv++;
					//if(argc<2)printCommandLineError("Missing graph file");
					if(argc<1)printCommandLineError("Missing column number");
					col=atoi(*argv);
					argc--;argv++;
				}else{
					mxofval=argv[0];
					argc--;argv++;
				}
			}
        		MatrixOfValues MX(mxofval,col);
			MX.printMatrix();
			exit(0);
		}
		else if(strcmp(*argv,"--graph-nodes")==0||strcmp(*argv,"--matrix-nodes")==0){
			argc--;argv++;
			mxofval="/dev/stdin";
			while(argc>0){
				mxofval=argv[0];
				argc--;argv++;
			}
        		MatrixOfValues MX(mxofval);
			MX.printNodes();
			exit(0);
		}
		else {
			printCommandLineError();
			exit(1);
		}
	}
#ifdef DEBUG
	QUIET=false;
	cout<<"#DEBUG MODE\tRUNNING IN DEBUGGING MODE\tQUIET=false (ignoring -q option). Additional debugging info with '-V' option"<<endl;
	cout<<"#DEBUG MODE: program task : analysis="<<analysis<<endl;
#endif
	switch(analysis){
///(For analyzing partitions)
		case prgCDIS:
		case prgCCOP:{
        		MatrixOfValues MX(mxofval);
			if(beta==BETA_UNSET)beta=0.01;
			if(mu==MU_UNSET)mu=0.0;
		        if(!QUIET) systemDate();
		        Partition partition(partitionf1,piformat,cluster1_offset);
			//if(MCLTABF)partition.mclTabFile(partitionf2);
			if(MCLTABF)partition.swapLabels(mcltabfile);
			ccop do_ccop(&MX,&partition,threshold);
			switch(analysis){
				case prgCCOP:	
					do_ccop.checkConsistency();
					break;
				case prgCDIS:	
					do_ccop.distribution();
					break;
			}
			break;
			}
		case prgVIPP:
		case prgEDSC:
		case prgMEET:
		case prgJOIN:
		{ 
			Partition partition1(partitionf1,piformat,cluster1_offset);
			Partition partition2(partitionf2,piformat,cluster2_offset);
			if(!QUIET) systemDate();
			if(analysis==prgJOIN){
				Partition join_1_2=partition1+partition2;
				join_1_2.printPartition(piformat);
				exit(0);
			}
			Partition intersection_1_2=partition1*partition2;
			if(analysis==prgMEET){
				if(!QUIET)cout<<"#Printing Intersection"<<endl;
				intersection_1_2.printPartition(piformat);
				exit(0);
			}
			if(!QUIET)cout<<"#dist-inters-p1"<<endl;
			if(analysis==prgVIPP)
				partition1.vipp(&intersection_1_2);
			else
				partition1.edsc(&intersection_1_2);
			if(!QUIET)cout<<"#dist-inters-p2"<<endl;
			if(analysis==prgVIPP)
				partition2.vipp(&intersection_1_2);
			else
				partition2.edsc(&intersection_1_2);
			if(!QUIET)cout<<"#dist-p1-p2"<<endl;
			if(analysis==prgVIPP)
		        	partition1.vipp(&partition2);
			else
		        	partition1.edsc(&partition2);
			partition1.missing();
			break;
			}
		case prgPSPP:{
			//Reference partition
			Partition partition1(partitionf1,piformat,cluster1_offset);
			//Target partition
			Partition partition2(partitionf2,piformat,cluster2_offset);
			//target againts reference
			vector<double> pscores=partition2.purityScore(&partition1);
			//cout<<"purityLax= "<<pscores[1]<<" purityStrict= "<<pscores[0]<<endl;
			cout<<" purityStrict= "<<pscores[0]<<"\tpurityLax= "<<pscores[1]<<endl;
			break;
			}
		case prgIPOT:
		case prgCPOT:
		case prgJPOT:
		case prgMPOT:
		case prgCMPO:
		case prgSSSA:
		case prgPSYM:
		case prgVMAM:
		case prgVMGM:
		case prgVMHM:
		case prgCPOR:
		case prgJPOR:
		case prgMPOR:
		case prgCMPR:
		case prgSSSR:
		case prgVAMR:
		case prgVGMR:
		case prgVHMR:
		case prgPSYR:
		case prgSEXV: //Scan extensivity values
		case prgIPOG:
		case prgCPOG:
		case prgJPOG:
		case prgMPOG:
		case prgCMPG:
		case prgSSSG:
		{
			PartitionStats partstats(infilenames, piformat, extensivity, cluster1_offset, clstat_normalization_ofs);
			pmeasure measure;
			switch(analysis){
				case prgCPOT:
				case prgCPOR:
				case prgCPOG:
					measure=conditionalEntropy;
					break;
				case prgJPOT:
				case prgJPOR:
				case prgJPOG:
					measure=jointEntropy;
					break;
				//case prgMPOT: //Testing for subadditivity
				case prgCKSA: //Testing for subadditivity
				case prgMPOR: //Testing for subadditivity
				case prgMPOG: //Testing for subadditivity
					measure=mutualInformation;
					break;
				//case prgCMPO: //Check if Strong Subadditivity holds
				case prgCSSA: //Check if Strong Subadditivity holds
				case prgCMPR: //Check if Strong Subadditivity holds
				case prgCMPG: //Check if Strong Subadditivity holds
					measure=conditionalMutualInformation;
					break;
				case prgSSSA: //Check if Weak-Strong Subadditivity holds
				case prgSSSR: //Check if Weak-Strong Subadditivity holds
				case prgSSSG: //Check if Weak-Strong Subadditivity holds
					measure=SSSA;
					break;
				case prgVMAM:
				case prgVAMR:
					measure=vmeasureArithmetic;
					break;
				case prgVMGM:
				case prgVGMR:
					measure=vmeasureGeometric;
					break;
				case prgVMHM:
				case prgVHMR:
					measure=vmeasureHarmonic;
					break;
				case prgPSYM:
				case prgPSYR:
					measure=symmetricPurity;
					break;
			}
			switch(analysis){
				case prgIPOT: //==prgIPOTR
					partstats.iPotential(metric);
					break;
				case prgCPOT:
				case prgJPOT:
				case prgCKSA: //Testing for subadditivity
				case prgCSSA: //Check if Strong Subadditivity holds
				case prgSSSA: //Check if Weak-Strong Subadditivity holds
				case prgVMAM:
				case prgVMGM:
				case prgVMHM:
					partstats.pmeasures(measure,metric);
					break;
				case prgPSYM:
					partstats.pmeasures(measure);
					break;
				case prgCPOR:
				case prgJPOR:
				case prgMPOR: //Testing for subadditivity
				case prgKSSR: //Check if Strong Subadditivity holds
				case prgSSSR: //Check if Weak-Strong Subadditivity holds
				case prgVAMR:
				case prgVGMR:
				case prgVHMR:
					partstats.pmeasuresRef(measure,metric);
					break;
				case prgPSYR:
					partstats.pmeasuresRef(measure);
					break;
				case prgIPOG:
				case prgCPOG:
				case prgJPOG:
				case prgKSAG: //Testing for subadditivity
				case prgKSSG: //Check if Strong Subadditivity holds
				case prgSSSG: //Check if Weak-Strong Subadditivity holds
					printCommandLineError("Not yet implemented");
					break;
			}
			break;
			}
		case prgJOST:
		case prgMEST:
		case prgVIST:
		case prgEDST:
		case prgBDST:
		case prgTDST:
		case prgRDST:
		case prgQDST:
		case prgVISR:
		case prgEDSR:
		case prgBDSR:
		case prgTDSR:
		case prgRDSR:
		case prgQDSR:
		case prgPSSR:
		case prgPSTG:
		case prgPSST:{
			if(VERBOSE){
				cout<<"#Read "<<infilenames.size()<<" partitions: ";
				for(vector<Charr>::iterator fn=infilenames.begin();fn!=infilenames.end();fn++)
					cout<<fn->car<<"\t";
				cout<<endl;
				if(infilenames2.size()>0){
					cout<<"#Read 2nd list of "<<infilenames2.size()<<" partitions: ";
					for(vector<Charr>::iterator fn=infilenames2.begin();fn!=infilenames2.end();fn++)
						cout<<fn->car<<"\t";
					cout<<endl;
				}
			}
			PartitionStats partstats(infilenames,mcltabfile , extensivity);
			//PartitionStats partstats(infilenames, piformat, extensivity, cluster1_offset, clstat_normalization_ofs,mcltabfile);
			if(infilenames2.size()>0){
				partstats.defineListOfPartitions(infilenames2,mcltabfile,piformat);
			}
			//pmetricv metric=shannon;
			switch(analysis){
				case prgJOST:
					partstats.getJoin().print();
					break;
				case prgMEST:
					partstats.getMeet().print();
					break;
				case prgPSSR:
					partstats.getPurityRef();
					break;
				case prgPSTG:
					partstats.getPurityTarget();
					break;
				case prgPSST:
					partstats.getPurity();
					break;
				case prgVIST:
				case prgEDST:
				case prgBDST:
				case prgTDST:
				case prgRDST:
				case prgQDST:
					switch(analysis){
						case prgVIST: metric=shannon; break;
						case prgEDST: metric=cardinality ; break;
						case prgBDST: metric=boltzmann; break;
						case prgTDST: metric=tsallis; break;
						case prgRDST: metric=renyi; break;
						case prgQDST: metric=jeffreyQnorm; break;
						default: cout<<"ERROR: unknown metric option"; exit(1); break;
					}
					if(!DIST_SUBSPROJECT) partstats.distances(metric);
					else	partstats.distances_Subsprojection(metric);
					break;
				case prgVISR:
				case prgEDSR:
				case prgBDSR:
				case prgTDSR:
				case prgRDSR:
				case prgQDSR:
					switch(analysis){
						case prgVISR: metric=shannon; break;
						case prgEDSR: metric=cardinality; break;
						case prgBDSR: metric=boltzmann; break;
						case prgTDSR: metric=tsallis; break;
						case prgRDSR: metric=renyi; break;
						case prgQDSR: metric=jeffreyQnorm; break;
					}
					if(!DIST_SUBSPROJECT) partstats.distancesRef(metric);
					else	partstats.distancesRef_Subsprojection(metric);
					break;
			}
			break;
			}
		case prgSPST:
		case prgSPSO:
		case prgSPSS:{
			if(VERBOSE){
				cout<<"#Reading "<<infilenames.size()<<" partitions: ";
				for(vector<Charr>::iterator fn=infilenames.begin();fn!=infilenames.end();fn++)
					cout<<fn->car<<"\t";
				cout<<endl;
			}
			PartitionStats partstats(infilenames, piformat, extensivity, cluster1_offset, clstat_normalization_ofs);
			splitmethod method=overlap;
			switch(analysis){
				case prgSPSS:
					method=cosine;
					break;
				case prgSPSO:
					method=overlap;
					break;
				case prgSPST:
					method=split;
					break;
			}
			partstats.getSplitsRef(method); //Without arguments=>method=split
			break;
			}
		case prgADST:
		case prgASFP:{
			PartitionStats partstats(infilenames, piformat, extensivity, cluster1_offset, clstat_normalization_ofs);
			partstats.get_Ad();
			partstats.pgm_Ad();
			if(analysis==prgASFP) partstats.get_FuzzyConsensusPartition();
			break;
			}
		case prgCLST:{
			if(VERBOSE){
				cout<<"Reading "<<infilenames.size()<<" partitions: ";
				for(vector<Charr>::iterator fn=infilenames.begin();fn!=infilenames.end();fn++)
					cout<<fn->car<<"\t";
				cout<<endl;
			}
			if(!QUIET)cout<<"#Analyzing neighborhoods..."<<endl;
			PartitionStats partstats(infilenames, piformat, extensivity, cluster1_offset, clstat_normalization_ofs);
			partstats.getCover();
			if(!QUIET)cout<<"#Printing neighborhoods..."<<endl;
			partstats.printCover(SETUPCONSENSUSP);
			if(SETUPCONSENSUSP){
				if(!QUIET)cout<<"#Printing consensus partition..."<<endl;
				partstats.printConsensusPart();
			}
			Partition ConsensusPart = partstats.getConsensusPartition();
			if(!ConsensusPart.isaPartition()){
				cout<<"ERROR: consensus Partition is not a sound partition"<<endl;
				exit(1);
			}else{ if(!QUIET)cout<<"#Consensus Partition is a sound partition\n#clusters/elements: "<<ConsensusPart.n_clusters()<<" "<<ConsensusPart.n_items()<<endl;}
			break;
			}
		case prgIPAR:{
			if(VERBOSE){
				cout<<"Reading "<<infilenames.size()<<" partitions: ";
				for(vector<Charr>::iterator fn=infilenames.begin();fn!=infilenames.end();fn++)
					cout<<fn->car<<"\t";
				cout<<endl;
			}
			if(!QUIET)cout<<"#Analyzing if partitions are sound..."<<endl;
			PartitionStats partstats(infilenames, piformat, extensivity, cluster1_offset, clstat_normalization_ofs);
			if(!QUIET){systemDate(); cout<<"#Finished instantiating all partitions. Starting analysis..."<<endl;}
			int nopart= partstats.arePartitions();
			if(!QUIET)if(nopart>0){
					cout<<"#WRONG partitions found:\n"<<nopart<<endl;
					exit(1);
				}else
					cout<<"#OK All partitions are sound partitions"<<endl;
			break;
			}
		case prgHASS:{
			if(VERBOSE){
				cout<<"Reading "<<infilenames.size()<<" partitions: ";
				for(vector<Charr>::iterator fn=infilenames.begin();fn!=infilenames.end();fn++)
					cout<<fn->car<<"\t";
				cout<<endl;
			}
			PartitionStats partstats(infilenames, piformat, extensivity, cluster1_offset, clstat_normalization_ofs);
			partstats.printHasseNodes();
			partstats.printHasseDiagram();
			break;
			}
///(For creating partitions)
		case prgRPAR:
		case prgRPSC:
		case prgRDCB:
		case prgRDCA:
		case prgRPAA:
		case prgRPAB:{
			bool pruneBelow=true;
			bool RDCselfconsistently=false;
			switch(analysis){
				case prgRPSC:
				case prgRDCB:
					RDCselfconsistently=true;
					break;
				case prgRDCA:
					RDCselfconsistently=true;
				case prgRPAA:{
					pruneBelow=false;
					break;
					}
			}
			RobustDivisiveClustering pruneCluster(mxofval,pruneBelow,Nsamples,Nneighbors,metric,extensivity,RDCselfconsistently);
			pruneCluster.printPartition();
			pruneCluster.printDistanceVsPruningThreshold();
			break;
			}
		case prgPART:
		case prgPARA:
		case prgPARB:{
			MatrixOfValues mxa(mxofval);
			Partition p;
			switch(analysis){
				case prgPART:{
					p=mxa.cluster();
					break;
					}
				case prgPARA:{
					//MatrixOfValues prunedmx=mxa.pruneEdgesAbove(graph_pruning_thr);
					//p=prunedmx.cluster();
					p=mxa.pruneEdgesAbove(graph_pruning_thr).cluster();
					break;
					}
				case prgPARB:{
					//MatrixOfValues prunedmx=mxa.pruneEdgesBelow(graph_pruning_thr);
					//p=prunedmx.cluster();
					p=mxa.pruneEdgesBelow(graph_pruning_thr).cluster();
					break;
					}
			}
			p.printPartition();
			/*
			*/
			break;
			}
///(For editing partitions)
		case prgPACN:{
		        Partition partition(partitionf1,piformat,cluster1_offset);
			if(!MCLTABF){
				cout<<"Tab file not defined! Did you use option --tab??"<<endl;
				exit(1);
			}
			partition.tabFile(mcltabfile);
			partition.swapLabels();
			partition.printPartition(poformat);
			break;
			}
		case prgPASR:
		case prgPASO:{
		        Partition partition(partitionf1,piformat,cluster1_offset);
		        Partition npart(&partition.clusters,cluster1_offset);
			bool SequentialClusterNames=false;
			switch(analysis){
				case prgPASR:
					if(!QUIET)cout<<"#Using sequential cluster names"<<endl;
					SequentialClusterNames=true;
					break;
			}
			npart.printPartition(poformat,SequentialClusterNames,ClusterPrefix);
			break;
			}
		case prgPAXE:{
		        Partition partition(partitionf1,piformat,cluster1_offset);
			Partition npart=partition.xtractElements(&namelist);
			npart.printPartition(poformat);
			break;
			}
///(For converting between different partition formats)
		case prg2MCL:
		case prg2FRE:
		case prgM2PA:{
		        Partition partition(partitionf1,piformat);
			if(MCLTABF)partition.mclTabFile(mcltabfile);
			partition.printPartition(poformat);
			break;
			}
		case prgA2PA:{
                        MultipleSeqAlign msa;
			Partition part=msa.getPartition(partitionf1);
                        part.printPartition();
			break;
			}
///(For dealing with (fasta) sequence files
		case prgSNOC:{
                        MultipleSeqAlign msa(msaf);
                        if(!QUIET)cout<<"#Dropping cloned sequences ...";
                        if(msafb==NULL){
				if(!QUIET)cout<<endl;
				msa.dropClones().print();
			}else{
				if(!QUIET)cout<<"#against "<<msafb<<endl;
                        	MultipleSeqAlign msab(msafb);
				msa.dropClones(&msab).print();
			}
			break;
			}
///(For analyzing Multiple Sequence Alignments)
		case prgPMSA:
		case prgPMAS:{
                        MultipleSeqAlign msa(msaf);
                        if(!QUIET)cout<<"#avgSeqId "<<msa.averageId()<<endl;
			bool sorted=false;
                        switch(analysis){
				case prgPMAS:
					sorted=true;
					break;
			}
			msa.print(sorted);
			break;
			}
		case prgMSMP:{
                        MultipleSeqAlign msa(msaf);
			Partition partition1(partitionf1,piformat,cluster1_offset);
			msa.printWithClusterLabels(&partition1, msafmt);
			break;
			}
		case prgMSPI:{
                        MultipleSeqAlign msa(msaf);
                        cout<<"#avgSeqId "<<msa.averageId(&positions)<<endl;
                        msa.printPairwiseIds(&positions);
			break;
			}
		case prgMSXP:{
                        MultipleSeqAlign msa(msaf);
                        MultipleSeqAlign nmsa=msa.xtractPositions(&positions);
			nmsa.print();
			break;
			}
		case prgMSDT:
		case prgMSXT:{
                        MultipleSeqAlign msa(msaf);
                        MultipleSeqAlign msab(msafb);
			MultipleSeqAlign nmsa=msa.xtractSequencesHighestId(&msab,Nsamples,&positions);
			nmsa.print();
			break;
			}
		case prgMSDI:
		case prgMSXI:{
                        MultipleSeqAlign msa(msaf);
                        MultipleSeqAlign msab(msafb);
			bool cull=true;
			switch(analysis){
				case prgMSDI:
                        		cull=false;
					break;
			}
			MultipleSeqAlign nmsa=msa.xtractSequencesById(&msab,doublea,doubleb,cull,&positions);;
			nmsa.print();
			break;
			}
		case prgMSDS:
		case prgMSXS:{
                        MultipleSeqAlign msa(msaf);
			bool cull=true;
			switch(analysis){
				case prgMSDS:
                        		cull=false;
					break;
			}
			MultipleSeqAlign nmsa=msa.xtractSequences(&namelist,cull);;
			nmsa.print();
			break;
			}
		case prgMAPI:{
                        MultipleSeqAlign msa(msaf);
                        cout<<"#avgSeqId "<<msa.averageId()<<endl;
                        msa.printAveragePairwiseIds(Id_threshold);
			break;
			}
		case prgMMPI:
		case prgMMAI:{
                        MultipleSeqAlign msa1(msaf);
                        MultipleSeqAlign msa2(msafb);
                        cout<<"#avgSeqId-MSA1 "<<msa1.averageId()<<endl;
                        cout<<"#avgSeqId-MSA2 "<<msa2.averageId()<<endl;
                        switch(analysis){
				case prgMMAI:
					msa1.printAveragePairwiseIds(msa2,Id_threshold,&positions);
					break;
				case prgMMPI:
					msa1.printPairwiseIds(msa2,&positions);
			}
			break;
			}
		case prgMRED:{
                        MultipleSeqAlign msa(msaf);
			if(seed>0)msa.NsamplesRedMSA(seed,Nsamples,Nseq);
			else msa.NsamplesRedMSA(Nsamples,Nseq);
			break;
			}
///(For dealing with -interaction- matrices) 
		case prgMGMX:
		case prgCEMX:
		case prgMMXC:{
			MatrixOfValues mxa(mxofval);
			MatrixOfValues mxb(mxofvalb);
			switch(analysis){
				case prgMGMX:
					mxa.merge(&mxb);
					break;
				case prgMMXC:{
					Partition partition(partitionf1,piformat,cluster1_offset);
					mxa.merge(&mxb,&partition);
					break;
					}
				case prgCEMX:
					mxa.cull(&mxb);
					break;
			}
			/*
			if(analysis==prgMGMX){
				mxa.merge(&mxb);
				exit(0);
			}
			Partition partition(partitionf1,piformat,cluster1_offset);
			mxa.merge(&mxb,&partition);
			*/
			//MX.printMatrix();
			break;
			}
		case prgGPRA:
		case prgGPRB:{
			MatrixOfValues mxa(mxofval);
			//MatrixOfValues prunedmx;
			switch(analysis){
				case prgGPRA:{
					//prunedmx=mxa.pruneEdgesAbove(graph_pruning_thr);
					mxa.pruneEdgesAbove(graph_pruning_thr).printMatrix();
					break;
					}
				case prgGPRB:{
					//prunedmx=mxa.pruneEdgesBelow(graph_pruning_thr);
					mxa.pruneEdgesBelow(graph_pruning_thr).printMatrix();
					break;
					}
			}
			//prunedmx.printMatrix();
			break;
			}
		case prgEDMX:
		case prgEDMP:{
			MatrixOfValues mxa(mxofval);
			switch(analysis){
				case prgEDMP:{
					if(VERBOSE)cout<<"#With partition"<<endl;
					Partition partition(partitionf1,piformat,cluster1_offset);
					mxa.edgeDistribution(&partition);
					break;
					}
				case prgEDMX:
					if(VERBOSE)cout<<"#WithOUT partition"<<endl;
					mxa.edgeDistribution();
					break;
			}
			break;
			}
		default:
			cout<<"ERROR: main: not listed program task : analysis="<<analysis<<endl;
			exit(1);
			//exitWithHelp();
			break;
	}
	//partition.printPartition();
	//partition1.do_vipp(&partition2,threshold);
	if(!QUIET) systemDate();
#ifdef DEBUGDETAILS
	cout<<"#DEBUG MODE\tRUNNING IN DEBUGGING MODE\tQUIET=false (ignoring -q option)"<<endl;
#endif

	return 0;
}
