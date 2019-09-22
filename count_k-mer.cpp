/*** g++  count_k-mer.cpp -fPIC -shared -o  count_k-mer.so ***/

#include <tr1/unordered_map>
//#include <unordered_map>
#include <vector>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
using namespace std::tr1;
using namespace std;

//
//
////unsigned long power;//parameter, power = 4^(k-1)
//
////receptacle of kmer counts
unordered_map<unsigned long,unsigned long> HashTable;
////
//int seqlength=999999999;
//// seq: save read/genome(A C G T) in each line from file
//char seq[999999999];
//// inverse of seq
//char seq_inverse[999999999];
int ZI = 4;


//
//void printFour(vector<unsigned long> four)
//{
//	cout << "print ";
//	for(int it=0; it<four.size(); it++)
//	{
//		cout << four[it] << "," ;
//	}
//	cout << endl;
//
//}


vector<int> ten2four(unsigned long ten, int k)
{
	vector<int> four (k,0);
	unsigned long tmp = ten;
	for(int currentPos = k-1; currentPos >=0; --currentPos)
	{
		four[currentPos]=tmp%ZI;
		tmp/=ZI;
	}
	//while(tmp>=(ZI-1)) {four[currentPos]=tmp%ZI;tmp/=ZI;currentPos--; }
	//four[currentPos] = tmp;
	return four;
}


int four2ten(vector<int> four, int k)
{
	unsigned long ten = 0;
	for(int currentPos=(k-1); currentPos >= 0; --currentPos)
	{
		int tmp = four[currentPos] * pow(ZI,(k-1 - currentPos));
		ten = ten + tmp;
		//cout << currentPos << " " << ten << endl;
	}
	return ten;

}


// find the reverse compilmentary of a word
vector<int> reverseFour(vector<int> Four)
{
	vector<int> reverseFour(Four.size(), 4);
	for(int revPos = 0; revPos < Four.size(); revPos++)
	{
		reverseFour[revPos] = 3 - Four[Four.size()- 1 - revPos];
	}
	return reverseFour;

}



unsigned long SeqKmerCountSingle(vector<char> seqDNA, int k, unsigned long power)
{
	//..........................
	// int count: The length of char seq
	int count = 0;
	//..........................
	int j=0;
	unsigned long index = 0;
	unsigned long total = 0;//total number of the kmers counted in seq
  //cout << "seqDNA.size()=" << seqDNA.size() << endl;
  for( int i=0; i < seqDNA.size(); i++) //while(seqDNA[i])
	{
    //cout << "seqNDA[i]=" << seqDNA[i] << ",if==A" << seqDNA[i]=='A' << endl;
		//kmer in seq[i,i+1,i+2,...i+k-1] transfered to an index
		//current index = floor(previous index/4^(k-1))*4 + 0 or 1 or 2 or 3
    char seqPos = seqDNA[i];
		if(seqPos == 'A'|| seqPos == 'a') {j++; }
		else if(seqPos == 'C'|| seqPos == 'c') { j++; index++; }
		else if(seqPos == 'G'|| seqPos == 'g') { j++; index+=2; }
		else if(seqPos == 'T'|| seqPos == 't') { j++; index+=3; }
		else { j=0; index=0; }//If seq[i] is ambiguous, reset j and index

    //cout << "index" << index << endl;
		if( j == k )
		{
			HashTable[index]++;
			total++;
			index %= power;// current index = floor(previous index/4^(k-1))
			j--;//the lengh of seq[i+1,i+2,...i+k-1]
		}
		index*=ZI;//current index = floor(previous index/4^(k-1))*4
		count++;
	}


	return total;
}





// compute EuFeature using double strand
void loadToVector(int k, unsigned long total, vector<double>& kmerCount)
{

	//int countWord = 0;
	for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
	{
		//cout << "test2.5" << endl;
		vector<int> currentKmerFour = ten2four(currentKmerTen, k);
		vector<int> currentKmerRevFour = reverseFour(currentKmerFour);

		unsigned long currentKmerRevTen = four2ten(currentKmerRevFour, currentKmerFour.size());
		//cout << currentKmerRevTen << endl;
		//printFour(currentKmerFour);
		//printFour(currentKmerFour);
		//printFour(currentKmerRevFour);
		//cout << "test2.6" << endl;

		if( currentKmerRevTen >= currentKmerTen )
		{
			//kmerTen.push_back(currentKmerTen);
			kmerCount.push_back((HashTable[currentKmerTen] + HashTable[currentKmerRevTen])/double(2 * total));
			//countWord ++;

			//cout << currentKmerTen << "," << (HashTable[currentKmerTen] + HashTable[currentKmerRevTen])/double(2 * total) << endl;
			//featureOutput << currentKmerTen << "," << TransToReal(X_w) << endl;
		}
	}

	//printFour(kmerTen);
	//featureOutput.close();

	//return X_w;
}

const unsigned long KMER_MAX_NUM = 65536;
typedef struct  
{  
    double kmer_freq[KMER_MAX_NUM];  
}kmerList;
typedef kmerList *kmerListPointer; 


// extern interface: countSeqFeatureCpp
extern "C"
{
	//double* countSeqFeatureCpp(char DNA[],  int k);
	//unsigned long kmerNum(int k);
//}
/*
//c++ version
vector<double> countSeqFeatureCpp( vector<char> RseqDNA,  int k) {

	unsigned long power = 1; for( int i = 0; i < k-1; i++) power *= 4;
	HashTable.clear();

	// count kmer
	unsigned long total = SeqKmerCountSingle(RseqDNA, k, power);

	// pair words and output count
	//vector<unsigned long> kmerTen;
	vector<double> kmerCount;
	loadToVector(k, total, kmerCount);
	//Rcout << "\n total:" << total << endl;

	// convert to Rcpp type
	//NumericVector RkmerTen(kmerTen.size());
	//RkmerTen = kmerTen;
	//NumericVector RkmerCount(kmerCount.size());
	//RkmerCount = kmerCount;

	//List ret;
	//ret["kmerTen"] = kmerTen;
	//ret["kmerCount"] = kmerCount;
	return kmerCount;

}
*/
//c version
kmerListPointer countSeqFeatureCpp(char DNA[],  int k) {
  //convert to c++ type
	vector <char> vDNA(DNA,DNA+strlen(DNA));

	unsigned long power = 1; for( int i = 0; i < k-1; i++) power *= 4;
	HashTable.clear();

	// count kmer
	unsigned long total = SeqKmerCountSingle(vDNA, k, power);

	// pair words and output count
	//vector<unsigned long> kmerTen;
	vector<double> kmerCount;
	loadToVector(k, total, kmerCount);
	//Rcout << "\n total:" << total << endl;
	 
	// convert to c type
	//kmerListPointer ret = (kmerListPointer)malloc(sizeof(kmerListPointer));
	kmerList* ret = new kmerList; 
  if (!kmerCount.empty())
	{
			memcpy(ret->kmer_freq, &kmerCount[0], kmerCount.size()*sizeof(double));  
	}

	return ret;
}
//Calculate k-mer number
unsigned long kmerNum(int k)
{
		unsigned long kmer_num;
		unsigned long total = 1; for( int i = 0; i < k; i++) total *= 4;
		if(k%2 != 0)
		{
				kmer_num = total/2;
		}
		else
		{
				unsigned long unique_kmer = 1; for( int i = 0; i < k/2; i++) unique_kmer *= 4;
				kmer_num = unique_kmer + (total-unique_kmer)/2;
		}
		return kmer_num;
}
}

//Test
int main(){
	int k = 6;
	char dna[] = "ATGCCG";
	kmerList* k_mer = countSeqFeatureCpp(dna,k);
	for (int i=0; i < kmerNum(k); i++)
        cout<<k_mer->kmer_freq[i]<<endl;
  return 0;
}

