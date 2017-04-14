#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <limits>
#include <getopt.h>			// arguments

using namespace std;

# define PI 3.14159265358979323846

unsigned static maxProLen = 10000;
unsigned static maxRNAlen = 30000;
unsigned areaBufferSize = 30000;

void make_scale(ifstream &file, map< pair<char, char>, double> &scale){
	char AA, BS;
	for(int i=0; i<20; i++){
		file >> AA;
		file >> scale[make_pair(AA, 'G')];
		file >> scale[make_pair(AA, 'A')];
		file >> scale[make_pair(AA, 'C')];
		file >> scale[make_pair(AA, 'U')];
	}
}

void make_energies(int proLen, int rnaLen, const string &pro, const string &rna, map< pair<char, char>, double> &scale, vector< vector<double> > &energies){	
	char AA, BS;
	for(int i=0; i<proLen; i++){
		AA = pro[i];
		for(int j=0; j<rnaLen; j++){
			BS = rna[j];
			if(scale.find(make_pair(AA,BS)) == scale.end()){ cout << "non-canonical base/AA\n"; std::exit(1); } 
			energies[i][j] = scale[make_pair(AA,BS)];
		}
	}
}

void print_energies(int proLen, int rnaLen, const string &pro, const string &rna, vector< vector<double> > &energies, int x1Min, int x2Min, int y1Min, int y2Min){
	cout << endl << "    ";
	for(int y=0; y<rnaLen; y++){
		printf("      %c", rna[y]);
	}
	cout << endl << endl;
	for(int x=0; x<proLen; x++){
		printf("  %c   ", pro[x]);
		for(int y=0; y<rnaLen; y++){
			if(x>=x1Min && x<=x2Min && y>=y1Min && y<=y2Min){ cout << '|'; }
			else{ cout << ' '; }
			printf("% 2.2f ", energies[x][y]);
		}
		if(x>=x1Min && x<=x2Min && y2Min==rnaLen-1){ cout << '|'; }
		cout << endl;
	}
	cout << endl;
}

void make_buffer(int proLen, int rnaLen, const vector< vector<double> > &energies, vector< vector<double> > &buffer){
	double energy;
	for(int x2=0; x2<proLen; x2++){
	for(int y2=0; y2<rnaLen; y2++){
		if(x2==0 && y2==0){
			energy = energies[x2][y2];
		}
		else if(x2==0){
			energy = energies[x2][y2] + buffer[x2][y2-1];
		}
		else if(y2==0){
			energy = energies[x2][y2] + buffer[x2-1][y2];
		}
		else{
			energy = energies[x2][y2] + buffer[x2-1][y2] + buffer[x2][y2-1] - buffer[x2-1][y2-1];
		}
		buffer[x2][y2] = energy;
	}
	}

}

void next_schedule(ifstream &file, const int proSize, const int rnaSize, bool &b_moreSchedule, int &pStart, int &pStop, int &rStart, int &rStop){
	bool b_num = false;
	char c;
	int numBuf, stage=0;
	string strNumBuf;
	pStart = 0; pStop = proSize; rStart = 0; rStop = rnaSize;
	while(file.good() && file.get(c)){
		if(isdigit(c)){ strNumBuf += c; }
		else if(c==':' || c=='-' || c=='\n'){
			if(strNumBuf.length() > 0){
				b_num = true;
				numBuf = atoi(strNumBuf.c_str()) - 1;	//start at 0 not 1
				strNumBuf.clear();
			}
			else{ b_num = false; }
			if(c==':'){
				if(! b_num){ stage++; }
				if(  b_num && stage==0){ pStart = numBuf; stage++; }
				if(  b_num && stage==2){ rStart = numBuf; stage++; }

			}
			if(c=='-'){
				if(! b_num){ stage = 2; }
				if(  b_num && stage==0){ pStart = numBuf; pStop = numBuf+1; stage = 2; }
				if(  b_num && stage==1){ pStop = numBuf; stage = 2; }
			}
			if(c=='\n'){
				if(  b_num && stage==2){ rStart = numBuf; rStop = numBuf+1; }
				if(  b_num && stage==3){ rStop = numBuf; }
				break;
			}
		}
	}
	file.get(c);
	if(file.good()){
		b_moreSchedule = true;
		file.putback(c);
	}
	else{ b_moreSchedule = false; }
	if(pStart < 0 || pStart >= proSize || pStop < 0 || pStop > proSize || pStart >= pStop || rStart < 0 || rStart >= rnaSize || rStop < 0 || rStop > rnaSize || rStart >= rStop){ cout << "schedule file: range error\n"; exit(1); }
}


int main(int argc, char** argv){
	string usage = 	" Calculates lowest energy regions from RNA/protein scale matrices for given sequences.\n"
			" If no arguments are given, all-to-all calculation is initiated.\n"
			"usage: affinity [options] <db> <scale>\n"
			"  -e print energy matrices\n  -i specify rna/prot pair by index\n"
			"  -n specify rna/prot pair by name (uniprot ID) (overrides -i)\n"
			"  -h display this and exit\n"
			"  -s schedule file\n"
			"Format:\t\t\tprotein ID; protein length; RNA ID; RNA length; protein start index; protein region length; RNA start index; RNA region length; energy of region; boltzmann probability\n"
			"Format schedule:\tInteger slices of protein/RNA sequence indexes with exclusive upper bound divided by a dash. Schedules are separated by newline.\n"
			"Examples for valid schedule formats: 'p:p-r:r' ':-r:r' 'p:p-'  '-r:' ':-:'  '-'\n";
	bool b_enMatrix = false;
	bool b_pairInd = false;
	bool b_pairName = false;
	bool b_schedule = false;
	int b_pairIndPro = -1;
	int b_pairIndRNA = -1;
	string b_pairNamePro, b_pairNameRNA;
	char option;
	ifstream scheduleFile;
	while((option = getopt(argc, argv, "heins:")) != -1){
		switch (option){
		case 'h':
			cout << usage;
			return 0;
		case 'e':
			b_enMatrix = true;
			break;
		case 'i':
			b_pairInd = true;
			b_pairIndPro = atoi(argv[ optind++ ]) - 1;
			b_pairIndRNA = atoi(argv[ optind++ ]) - 1;
			break;
		case 'n':
			b_pairName = true;
			b_pairInd = true;
			b_pairNamePro = argv[ optind++ ];
			b_pairNameRNA = argv[ optind++ ];
			break;
		case 's':
			b_schedule = true;
			scheduleFile.open(optarg);
			if(! scheduleFile.good()){ cout << usage << "schedule file missing!\n"; return 1; }
			break;
		}
	}
	if (argc - optind < 2) { cout << usage << "mandatory argument missing!\n"; return 1; }
	ifstream databFile(argv[optind++]);
	ifstream scaleFile(argv[optind++]);
	if(!databFile || !scaleFile){ cout << usage; return(1); } 

	bool b_moreSchedule = true;
	int i, proLen, rnaLen, x1, x2, y1, y2, x1Min, x2Min, y1Min, y2Min, P_x1Min, P_x2Min, P_y1Min, P_y2Min, minProLen, minRNAlen, areaCount, n, pStart, pStop, rStart, rStop, p, r;
	double energy, minEn, P_minEn, mean, sigma, pValue, sum;
	string name, sequ, ignore, type, pro, rna, proName, rnaName;
	vector< string > proNames;
	vector< string > rnaNames;
	vector< string > proSequs;
	vector< string > rnaSequs;
	vector< double > areaEn ( areaBufferSize, 0 );
	vector< vector<double> > energies ( maxProLen, vector<double> ( maxRNAlen, 0 ) );
	vector< vector<double> > buffer ( maxProLen, vector<double> ( maxRNAlen, 0 ) );
	map< pair<char, char>, double> scale;

	make_scale(scaleFile, scale);
	while(databFile.good()){
		databFile >> name >> sequ >> ignore >> type;
		if(databFile.eof()){ break; }
		if(type.compare("PROT") == 0){
			if(sequ.length() > maxProLen){ printf("protein %s: size limit reached (%d residues)\n", name.c_str(), maxProLen ); return 1; }
			proNames.push_back(name);
			proSequs.push_back(sequ);
		}
		else if(type.compare("RNA") == 0){
			if(sequ.length() > maxRNAlen){ printf("RNA %s: size limit reached (%d bases)\n", name.c_str(), maxRNAlen ); return 1; }
			rnaNames.push_back(name);
			rnaSequs.push_back(sequ);
		}
		else{
			cout << "invalid DB file\n";
			return 1;
		}
	}
	if(b_pairName){
		for(i=0; i<proNames.size(); i++){ if(proNames[i].compare(b_pairNamePro) == 0){ b_pairIndPro = i; }}
		for(i=0; i<rnaNames.size(); i++){ if(rnaNames[i].compare(b_pairNameRNA) == 0){ b_pairIndRNA = i; }}
		if(b_pairIndPro<0 || b_pairIndRNA<0){ cout << "pair names not found\n"; return 1; }
	}
	while(b_moreSchedule){
		next_schedule(scheduleFile, proNames.size(), rnaNames.size(), b_moreSchedule, pStart, pStop, rStart, rStop);
		for(p=pStart; p<pStop; p++){
		for(r=rStart; r<rStop; r++){
			if(b_pairInd){ p = b_pairIndPro; r = b_pairIndRNA; b_moreSchedule = false; }
			pro = proSequs[p];
			rna = rnaSequs[r];
			proName = proNames[p];
			rnaName = rnaNames[r];
			proLen = pro.length();
			rnaLen = rna.length();
			make_energies(proLen, rnaLen, pro, rna, scale, energies);
			make_buffer(proLen, rnaLen, energies, buffer);
	
			minEn = numeric_limits<double>::max();
	#pragma omp parallel for ordered schedule(dynamic) private(x2, y1, y2, energy, P_minEn, P_x1Min, P_x2Min, P_y1Min, P_y2Min) shared(x1, buffer, minEn, x1Min, x2Min, y1Min, y2Min)
		for(x1=0; x1<proLen; x1++){
			P_minEn = numeric_limits<double>::max();
			for(y1=0; y1<rnaLen; y1++){
			for(x2=x1; x2<proLen; x2++){
			for(y2=y1; y2<rnaLen; y2++){
				if(x1==0 && y1==0){
					energy = buffer[x2][y2];
				}
				else if(x1==0){
					energy = buffer[x2][y2] - buffer[x2][y1-1];
				}
				else if(y1==0){
					energy = buffer[x2][y2] - buffer[x1-1][y2];
				}
				else{
					energy = buffer[x2][y2] - buffer[x1-1][y2] - buffer[x2][y1-1] + buffer[x1-1][y1-1];
				}
				if(energy < P_minEn){
					P_minEn = energy;
					P_x1Min = x1; P_y1Min = y1; P_x2Min = x2; P_y2Min = y2;
				}
			}
			}
			}
		#pragma omp ordered
			{
				if(P_minEn < minEn){
					minEn = P_minEn;
					x1Min = P_x1Min; y1Min = P_y1Min; x2Min = P_x2Min; y2Min = P_y2Min;
				}
			}
		}
	
			minProLen = x2Min - x1Min + 1;
			minRNAlen = y2Min - y1Min + 1;
			areaCount = (proLen - minProLen + 1) * (rnaLen - minRNAlen + 1);
			if(areaBufferSize < areaCount){
				areaBufferSize = areaCount;
				areaEn.resize(areaCount, 0);
			}
			n = 0;
			sum = 0;
		for(x1=0; x1 < proLen-minProLen+1; x1++){
		for(y1=0; y1 < rnaLen-minRNAlen+1; y1++){
		x2 = x1 + minProLen - 1;
		y2 = y1 + minRNAlen - 1;
			if(x1==0 && y1==0){
				energy = buffer[x2][y2];
			}
			else if(x1==0){
				energy = buffer[x2][y2] - buffer[x2][y1-1];
			}
			else if(y1==0){
				energy = buffer[x2][y2] - buffer[x1-1][y2];
			}
			else{
				energy = buffer[x2][y2] - buffer[x1-1][y2] - buffer[x2][y1-1] + buffer[x1-1][y1-1];
			}
			areaEn[n++] = energy;
			sum += energy;
		}
		}
			mean = sum / areaCount;
			sum = 0;
			for(i=0; i<areaCount; i++){ sum += pow((areaEn[i] - mean), 2); }
			sigma = sqrt(sum / areaCount);
			pValue = 1 / (sigma * sqrt(2 * PI)) * exp(- pow((minEn - mean), 2) / (2 * sigma * sigma) );
	
			printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%f\t%e\t%f\t%f\t%d\n", proName.c_str(), proLen, rnaName.c_str(), rnaLen, x1Min+1, x2Min-x1Min+1, y1Min+1, y2Min-y1Min+1, minEn, pValue, mean, sigma, areaCount); cout << flush;
			if(b_enMatrix){ print_energies(proLen, rnaLen, pro, rna, energies, x1Min, x2Min, y1Min, y2Min); }
			if(b_pairInd){ return 0; }
		}
		}
	}
}
