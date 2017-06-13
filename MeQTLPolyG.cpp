#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <unistd.h> 
#include <armadillo>
#include "Util.h"
#include "PostCal.h"
#include "TopKSNP.h"
#include "MeQTLPolyGModel.h"

using namespace std;


int main( int argc, char *argv[]  ){
	int totalCausalSNP = 2;
	double gamma = 0.01;
	double rho = 0.95;
	bool histFlag = false;
	int oc = 0;	
	string ldFile = "";
	string yFile  = "";
	string outputFileName = "";
	string geneMapFile = "";	
        string weight = "";
        string covariate = "";
        string grm_test ="";
        vector<string> grm_file;
        string X_file="";
        int number;
        int nthread=1;
	while ((oc = getopt(argc, argv, "vhl:t:o:x:p:n:w:g:r:c:G:w:f:")) != -1) {
		switch (oc) {
			case 'v':
				cout << "version 1.0:" << endl;
			case 'h':
				cout << "Options: " << endl;
  				cout << "-h, --help            		show this help message and exit " << endl;
  				cout << "-o OUTFILE, --out=OUTFILE 	specify the output file" << endl;
  				cout << "-l LDFILE, --ld_file=LDFILE  	the ld input file" << endl;
  				cout << "-p yFILE, --y_file=yFILE	phenotype" << endl;
  				cout << "-r RHO, --rho-prob=RHO		set $pho$ probability (default 0.95)" << endl;
				cout << "-g GAMMA, --gamma		set $gamma$ the prior of a SNP being causal (default 0.01)" << endl;
				cout << "-c causal			set the maximum number of causal SNPs" << endl;
                                cout << "-C covariate, --covariate      set the covariate matrix "<<endl;
                                cout << "-G genetic relatedness matrix, --GRM  set the genetic relatedness matrix "<<endl;
				cout << "-f 1				to out the probaility of different number of causal SNP" << endl;
                                cout << "-w Weight file, --weight       set the biological annotation to use" << endl;
                                cout << "-n Number of samples, --number       set the biological annotation to use" << endl;
                                cout << "-x genotype information,       genotype file for explored variants" << endl;
                                cout << "-t threads to use, --nthread       set the threads to use" << endl;
				exit(0);
			case 'l':
				ldFile = string(optarg);
				break;
                        case 'n':
                                number = atoi(optarg);
                                break;
			case 'o':
				outputFileName = string(optarg);
				break;
			case 'p':
				yFile = string(optarg);
				break;
			case 'r':
				rho = atof(optarg);
				break;
                        case 'x':
                                X_file = string(optarg);
                                break;
			case 'c':
				totalCausalSNP = atoi(optarg);
				break;
			case 'g':
				gamma = atof(optarg);
				break;
                        case 'w':
                                weight=string(optarg);
                                break;
                        case 't':
                                nthread=atoi(optarg);
			case 'f':
                                histFlag = true;
                                break;
                        case 'G':
                                grm_test     =string(optarg);
                                grm_file.push_back(grm_test);
                                break;
                        case 'C':
                                covariate=string(optarg);
                                break;
			case ':':
			case '?':
			default:
				cout << "Strange" << endl;
				break;
		}
	}
        cout<<"Getting parameter information is over"<<endl;
        CaviarModel caviar(ldFile, yFile, outputFileName, totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate, grm_file,X_file, number);
	caviar.run();
	caviar.finishUp();		
	return 0;
}
