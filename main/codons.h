#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>

using namespace std;

//enum t_nameaa {A='A',C='C',D='D',E='E',F='F',G='G',H='H',I='I',K='K',L='L',M='M',N='N',P='P',Q='Q',R='R',S='S',T='T',V='V',W='W',Y='Y'};
enum t_nameaa {A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y};

struct t_codon {
 char code[3];
 int frequency;
 double weight;
};

struct t_aminoacid {
 t_nameaa name;
 int num_codons;
 t_codon codons[6];
};


void loadtableAA(t_aminoacid &aminoacids);