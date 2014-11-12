/*
Nicholas Botzer
William Botzer
Steven Zamborsky
*/
#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stack>
#include <fstream>
using namespace std;
 
const int NUM_GENES = 1000;
int count;

int longestCommonSequence(unsigned short *matrix, string &gene1, int p1, string &gene2, int p2);
stack <char> walkback(unsigned short *matrix, string &gene1, string &gene2);
void printMatrix(unsigned short *matrix, string& gene1, string& gene2);
void readFile(ifstream& file, string gene[]);
string readGene(ifstream& file);
 
int main (int argc, char* argv[]) {

/******* INIT *******/
	stack <char> walkbackStack;
	ifstream inFile1;

    //Create LCS matrix
	string genes[NUM_GENES] = {""};	// Goes accros the top
	cout << "Start Reading Input\n";
	readFile(inFile1, genes);
	cout << "Done Reading Input\n";
	
	int g1 = 0, g2 = 1;
	string longestSequence = "";
	
/***** End INIT *****/
     

/***** Compute *****/
	clock_t startWatch = clock();	// start the timing
	cout << genes[59].length() << endl;
    //Compute LCS matrix
	for (int loc=0; loc<NUM_GENES-1; loc++) {
		for (int i=loc+1; i<NUM_GENES; i++) {
			cout << "LCS on genes: " << loc << " and " << i << endl ;
			unsigned short* matrix = new unsigned short [ genes[loc].length()*genes[i].length() ];
			memset(matrix, -1, genes[loc].length()*genes[i].length()*sizeof(short));
			
			for (int r=0; r<genes[loc].length(); r++) {
				for(int c=0; c<genes[i].length(); c++) {
					longestCommonSequence(matrix, genes[loc], r, genes[i], c);
				}
			}

			//Walkback through matrix
			walkbackStack = walkback(matrix, genes[loc], genes[i]);

			if (walkbackStack.size() > longestSequence.length()) {
				cout << "New Longest Sequence\n";
				g1 = loc;
				g2 = i;
				longestSequence = "";
				while (walkbackStack.size() > 0) {
						longestSequence += walkbackStack.top();
						walkbackStack.pop();
					}
			}
			free(matrix);
			//delete [] matrix;
		}
	}	
	clock_t stopWatch = clock();	// End the timer
	cout << "\n\n" << double(stopWatch - startWatch) / CLOCKS_PER_SEC << " secs" << endl;
    
	cout << "Longest sequence between: " << g1 << " and " << g2 << endl;
    cout << longestSequence;
	
	
	
	
/***** End Compute *****/
	
 
/***** Results *****/
	

	cout << endl;
	system("PAUSE");
    return 0;
}

int longestCommonSequence(unsigned short *matrix, string &gene1, int p1, string &gene2, int p2){
	int result = 0;
	if (p1 >= 0 && p2 >= 0){
		if ( matrix[p2*gene1.length() +p1] != -1)
			result = matrix[p2*gene1.length() +p1];
		else if (gene1[p1] == gene2[p2])
			result += 1 + longestCommonSequence(matrix, gene1, p1-1, gene2, p2-1);
		else
			result = max(longestCommonSequence(matrix, gene1, p1, gene2, p2-1), longestCommonSequence(matrix, gene1, p1-1, gene2, p2));

		matrix[p2*gene1.length() + p1] = result;
	}

	return result;
}

stack <char> walkback(unsigned short *matrix, string &gene1, string &gene2) {
	stack <char> resultStack;
	bool endWalkback = false;
	int p1 = gene1.length()-1;
	int p2 = gene2.length()-1;
	while (!endWalkback) {
		if (gene1[p1] == gene2[p2]) {
			resultStack.push(gene1[p1]);
		}
		if (p1 > 0 && p2 > 0){
			if ( matrix[p2*gene1.length()+p1] == matrix[p2*gene1.length()+p1-1])
				p1 -= 1;
			else if (matrix[p2*gene1.length()+p1] == matrix[(p2-1)*gene1.length()+p1])
				p2 -= 1;
			else if (matrix[p2*gene1.length()+p1] > matrix[(p2-1)*gene1.length()+p1] && matrix[p2*gene1.length()+p1] > matrix[p2*gene1.length()+p1-1]) {
				p1 -= 1;
				p2 -= 1;
			}
		}
		else if(p1 == 0 && p2 > 0) {
			p2 -= 1;
		}
		else if(p1 > 0 && p2 == 0) {
			p1 -= 1;
		}
		else if (p1 == 0 && p2 == 0) {
			endWalkback = true;
		}
		
	}

	return resultStack;
}

void readFile(ifstream& file, string genes[]) {
	file.open("data1.txt");
	if (file.is_open()) {
		for (int i = 0; i < NUM_GENES; i++) {
			genes[i] = readGene(file);
			if (genes[i].length() > 25000)
				i--;
		}
	}
	else
		cerr << "*****Unable to open " << file << " *****" << endl;

	file.close();
}
string readGene(ifstream& file) {
	string result;
	bool endFlag = false;
	int state = 1;
	char ch = ' ';
	
	while (!endFlag && !file.eof()) {
		switch (state){
		case 1:	//waiting for gene data
			ch = file.get();
			if (ch == '.') 
				state=3;
		break;

		/*case 2:	//check for gene data in next char
			ch = file.get();
			cout << ch << ", ";
			if ((int)ch == 10 || (int)ch == 32){
				state = 1;
			}
			if (ch == 'a' || ch == 'c' || ch == 't' || ch == 'g') {
				result = result + ch;
				state=3;
			}
			else
				state = 1;
		break;*/

		case 3:	//In gene data
			ch = file.get();
			if(ch == '>'){
                endFlag = true;
				count ++;
			}
            else if (ch == 'a' || ch == 'c' || ch == 't' || ch == 'g') {
				result = result + ch;
			}
		break;

		default:
			cout << "***** Error *****\n";
		break;
		}
	}
	return result;
}

void printMatrix(int *matrix, string &gene1, string &gene2){
	int spacing = 3;

	for(int r=0; r<gene2.length(); r++) {
		for (int c=0; c<gene1.length(); c++) {
			cout << setw(spacing) << matrix[r*gene1.length() +c];
		}
		cout << endl;
	}
	cout << endl;
}
