#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <float.h> 
#include "codons.h"
#include "math.h"
#include "lcs.h"
#include "lrs.h"

using namespace std;

const int NN = 20;              //Size for aminoacid table
//Config parameters
const int TAM = 100;            //Number of solutions to generate
const int FEMALES = 50;         //Number of females
const int PM = 10;              //Mutation percentage 
const int AGE = 10;             //Max dead age

//struct for save info about the aminoacids file 
struct datas
{
    std::string id;             // Chain id
    std::string chain;          // Chain of aminoacids
    int n;                      // Number of CDS 
};

/*struct dataList{
    std::vector<string> amns;
    std::vector<int> genders;
    std::vector<vector<string>> cds_list;
};*/

struct logData{
    int id;                             //Id experimento
    int tPareto;                        //Tamaño en frente de Pareto filtrado
    double tsecC;                       //Tiempos
    double tsecF;                       //
    int newBorns;                       //Numero de nacimientos en el experimento
    int males;                          //Numero de machos
    int females;                        //Numero de hembras
};

struct individual{
    std::vector<vector<string>> cds;    // CDSs
    int gender;                         // 0 -> male, 1 ->female, -1 ->no gender
    int cont=0;                         // Contador soluciones
    std::vector<long double> results;   // 0-> mCAI, 1-> mHD, 2-> MLRCS
    int cub=0;                          // 0->no, 1->cub
    long double crowding;               // Valor para crowding
    int rank;                           // Valor para rank
};


datas dts;                                      //informacion de la proteina obtenida por fichero

std::vector<individual> individualList;         //Structure that store all individuals 

std::vector<individual> paretoFrontier;         //individuals on Pareto frontier

std::vector<individual> paretoFrontierReduced;  //Reduction of Pareto frontier

logData logF[20];                               //Data for log file


//std::vector<datas> d;           // Structure that store the data of file readed

//dataList dList;     // Structure that store all  data generated

void loadtableAA(t_aminoacid *aminoacids) {

    //---------------A
    aminoacids[A].name=A;
    aminoacids[A].num_codons=4;
    
	aminoacids[A].codons[0].code[0]='G';
    aminoacids[A].codons[0].code[1]='C';
    aminoacids[A].codons[0].code[2]='A';
    aminoacids[A].codons[0].frequency=5296;
    aminoacids[A].codons[0].weight=0.390474084;

    aminoacids[A].codons[1].code[0]='G';
    aminoacids[A].codons[1].code[1]='C';
    aminoacids[A].codons[1].code[2]='C';
    aminoacids[A].codons[1].frequency=7223;
    aminoacids[A].codons[1].weight=0.532551795;

    aminoacids[A].codons[2].code[0]='G';
    aminoacids[A].codons[2].code[1]='C';
    aminoacids[A].codons[2].code[2]='G';
    aminoacids[A].codons[2].frequency=1854;
    aminoacids[A].codons[2].weight=0.136695421;

    aminoacids[A].codons[3].code[0]='G';
    aminoacids[A].codons[3].code[1]='C';
    aminoacids[A].codons[3].code[2]='U';
    aminoacids[A].codons[3].frequency=13563;
    aminoacids[A].codons[3].weight=1;

    //---------------C
    aminoacids[C].name=C;
    aminoacids[C].num_codons=2;
 
	aminoacids[C].codons[0].code[0]='U';
    aminoacids[C].codons[0].code[1]='G';
    aminoacids[C].codons[0].code[2]='C';
    aminoacids[C].codons[0].frequency=1234;
    aminoacids[C].codons[0].weight=0.404325033;

    aminoacids[C].codons[1].code[0]='U';
    aminoacids[C].codons[1].code[1]='G';
    aminoacids[C].codons[1].code[2]='U';
    aminoacids[C].codons[1].frequency=3052;
    aminoacids[C].codons[1].weight=1;

    //---------------D
    aminoacids[D].name=D;
    aminoacids[D].num_codons=2;
    
	aminoacids[D].codons[0].code[0]='G';
    aminoacids[D].codons[0].code[1]='A';
    aminoacids[D].codons[0].code[2]='C';
    aminoacids[D].codons[0].frequency=8960;
    aminoacids[D].codons[0].weight=0.703793889;

    aminoacids[D].codons[1].code[0]='G';
    aminoacids[D].codons[1].code[1]='A';
    aminoacids[D].codons[1].code[2]='U';
    aminoacids[D].codons[1].frequency=12731;
    aminoacids[D].codons[1].weight=1;

    //---------------E
    aminoacids[E].name=E;
    aminoacids[E].num_codons=2;
    
	aminoacids[E].codons[0].code[0]='G';
    aminoacids[E].codons[0].code[1]='A';
    aminoacids[E].codons[0].code[2]='A';
    aminoacids[E].codons[0].frequency=19532;
    aminoacids[E].codons[0].weight=1;

    aminoacids[E].codons[1].code[0]='G';
    aminoacids[E].codons[1].code[1]='A';
    aminoacids[E].codons[1].code[2]='G';
    aminoacids[E].codons[1].frequency=6172;
    aminoacids[E].codons[1].weight=0.315994266;

    //---------------F
    aminoacids[F].name=F;
    aminoacids[F].num_codons=2;
    
	aminoacids[F].codons[0].code[0]='U';
    aminoacids[F].codons[0].code[1]='U';
    aminoacids[F].codons[0].code[2]='C';
    aminoacids[F].codons[0].frequency=8251;
    aminoacids[F].codons[0].weight=1;

    aminoacids[F].codons[1].code[0]='U';
    aminoacids[F].codons[1].code[1]='U';
    aminoacids[F].codons[1].code[2]='U';
    aminoacids[F].codons[1].frequency=7773;
    aminoacids[F].codons[1].weight=0.942067628;

    //---------------G
    aminoacids[G].name=G;
    aminoacids[G].num_codons=4;
    
	aminoacids[G].codons[0].code[0]='G';
    aminoacids[G].codons[0].code[1]='G';
    aminoacids[G].codons[0].code[2]='A';
    aminoacids[G].codons[0].frequency=2781;
    aminoacids[G].codons[0].weight=0.177201478;

    aminoacids[G].codons[1].code[0]='G';
    aminoacids[G].codons[1].code[1]='G';
    aminoacids[G].codons[1].code[2]='C';
    aminoacids[G].codons[1].frequency=3600;
    aminoacids[G].codons[1].weight=0.229387027;

    aminoacids[G].codons[2].code[0]='G';
    aminoacids[G].codons[2].code[1]='G';
    aminoacids[G].codons[2].code[2]='G';
    aminoacids[G].codons[2].frequency=1852;
    aminoacids[G].codons[2].weight=0.118006882;

    aminoacids[G].codons[3].code[0]='G';
    aminoacids[G].codons[3].code[1]='G';
    aminoacids[G].codons[3].code[2]='U';
    aminoacids[G].codons[3].frequency=15694;
    aminoacids[G].codons[3].weight=1;

    //---------------H
    aminoacids[H].name=H;
    aminoacids[H].num_codons=2;
    
	aminoacids[H].codons[0].code[0]='C';
    aminoacids[H].codons[0].code[1]='A';
    aminoacids[H].codons[0].code[2]='C';
    aminoacids[H].codons[0].frequency=3288;
    aminoacids[H].codons[0].weight=0.761111111;

    aminoacids[H].codons[1].code[0]='C';
    aminoacids[H].codons[1].code[1]='A';
    aminoacids[H].codons[1].code[2]='U';
    aminoacids[H].codons[1].frequency=4320;
    aminoacids[H].codons[1].weight=1;

    //---------------I
    aminoacids[I].name=I;
    aminoacids[I].num_codons=3;
    
	aminoacids[I].codons[0].code[0]='A';
    aminoacids[I].codons[0].code[1]='U';
    aminoacids[I].codons[0].code[2]='A';
    aminoacids[I].codons[0].frequency=3172;
    aminoacids[I].codons[0].weight=0.26277856;

    aminoacids[I].codons[1].code[0]='A';
    aminoacids[I].codons[1].code[1]='U';
    aminoacids[I].codons[1].code[2]='C';
    aminoacids[I].codons[1].frequency=8251;
    aminoacids[I].codons[1].weight=0.683539061;

    aminoacids[I].codons[2].code[0]='A';
    aminoacids[I].codons[2].code[1]='U';
    aminoacids[I].codons[2].code[2]='U';
    aminoacids[I].codons[2].frequency=12071;
    aminoacids[I].codons[2].weight=1;

    //---------------K
    aminoacids[K].name=K;
    aminoacids[K].num_codons=2;
    
	aminoacids[K].codons[0].code[0]='A';
    aminoacids[K].codons[0].code[1]='A';
    aminoacids[K].codons[0].code[2]='A';
    aminoacids[K].codons[0].frequency=12845;
    aminoacids[K].codons[0].weight=0.846792801;

    aminoacids[K].codons[1].code[0]='A';
    aminoacids[K].codons[1].code[1]='A';
    aminoacids[K].codons[1].code[2]='G';
    aminoacids[K].codons[1].frequency=15169;
    aminoacids[K].codons[1].weight=1;

    //---------------L
    aminoacids[L].name=L;
    aminoacids[L].num_codons=6;
    
	aminoacids[L].codons[0].code[0]='C';
    aminoacids[L].codons[0].code[1]='U';
    aminoacids[L].codons[0].code[2]='A';
    aminoacids[L].codons[0].frequency=4134;
    aminoacids[L].codons[0].weight=0.310150799;

    aminoacids[L].codons[1].code[0]='C';
    aminoacids[L].codons[1].code[1]='U';
    aminoacids[L].codons[1].code[2]='C';
    aminoacids[L].codons[1].frequency=1242;
    aminoacids[L].codons[1].weight=0.093180284;

    aminoacids[L].codons[2].code[0]='C';
    aminoacids[L].codons[2].code[1]='U';
    aminoacids[L].codons[2].code[2]='G';
    aminoacids[L].codons[2].frequency=2852;
    aminoacids[L].codons[2].weight=0.21396954;

    aminoacids[L].codons[3].code[0]='C';
    aminoacids[L].codons[3].code[1]='U';
    aminoacids[L].codons[3].code[2]='U';
    aminoacids[L].codons[3].frequency=3207;
    aminoacids[L].codons[3].weight=0.240603196;

    aminoacids[L].codons[4].code[0]='U';
    aminoacids[L].codons[4].code[1]='U';
    aminoacids[L].codons[4].code[2]='A';
    aminoacids[L].codons[4].frequency=8549;
    aminoacids[L].codons[4].weight=0.64138345;

    aminoacids[L].codons[5].code[0]='U';
    aminoacids[L].codons[5].code[1]='U';
    aminoacids[L].codons[5].code[2]='G';
    aminoacids[L].codons[5].frequency=13329;
    aminoacids[L].codons[5].weight=1;

    //---------------M
    aminoacids[M].name=M;
    aminoacids[M].num_codons=1;
    
	aminoacids[M].codons[0].code[0]='A';
    aminoacids[M].codons[0].code[1]='U';
    aminoacids[M].codons[0].code[2]='G';
    aminoacids[M].codons[0].frequency=7911;
    aminoacids[M].codons[0].weight=1;

    //---------------N
    aminoacids[N].name=N;
    aminoacids[N].num_codons=2;
    
	aminoacids[N].codons[0].code[0]='A';
    aminoacids[N].codons[0].code[1]='A';
    aminoacids[N].codons[0].code[2]='C';
    aminoacids[N].codons[0].frequency=9875;
    aminoacids[N].codons[0].weight=1;

    aminoacids[N].codons[1].code[0]='A';
    aminoacids[N].codons[1].code[1]='A';
    aminoacids[N].codons[1].code[2]='U';
    aminoacids[N].codons[1].frequency=8613;
    aminoacids[N].codons[1].weight=0.872202532;

    //---------------P
    aminoacids[P].name=P;
    aminoacids[P].num_codons=4;
    
	aminoacids[P].codons[0].code[0]='C';
    aminoacids[P].codons[0].code[1]='C';
    aminoacids[P].codons[0].code[2]='A';
    aminoacids[P].codons[0].frequency=8965;
    aminoacids[P].codons[0].weight=1;

    aminoacids[P].codons[1].code[0]='C';
    aminoacids[P].codons[1].code[1]='C';
    aminoacids[P].codons[1].code[2]='C';
    aminoacids[P].codons[1].frequency=1656;
    aminoacids[P].codons[1].weight=0.184718349;

    aminoacids[P].codons[2].code[0]='C';
    aminoacids[P].codons[2].code[1]='C';
    aminoacids[P].codons[2].code[2]='G';
    aminoacids[P].codons[2].frequency=1064;
    aminoacids[P].codons[2].weight=0.11868377;

    aminoacids[P].codons[3].code[0]='C';
    aminoacids[P].codons[3].code[1]='C';
    aminoacids[P].codons[3].code[2]='U';
    aminoacids[P].codons[3].frequency=4575;
    aminoacids[P].codons[3].weight=0.510317903;

    //---------------Q
    aminoacids[Q].name=Q;
    aminoacids[Q].num_codons=2;
    
	aminoacids[Q].codons[0].code[0]='C';
    aminoacids[Q].codons[0].code[1]='A';
    aminoacids[Q].codons[0].code[2]='A';
    aminoacids[Q].codons[0].frequency=10987;
    aminoacids[Q].codons[0].weight=1;

    aminoacids[Q].codons[1].code[0]='C';
    aminoacids[Q].codons[1].code[1]='A';
    aminoacids[Q].codons[1].code[2]='G';
    aminoacids[Q].codons[1].frequency=3312;
    aminoacids[Q].codons[1].weight=0.301447165;

    //----------------R
    aminoacids[R].name=R;
    aminoacids[R].num_codons=6;
    
	aminoacids[R].codons[0].code[0]='A';
    aminoacids[R].codons[0].code[1]='G';
    aminoacids[R].codons[0].code[2]='A';
    aminoacids[R].codons[0].frequency=9784;
    aminoacids[R].codons[0].weight=1;

    aminoacids[R].codons[1].code[0]='A';
    aminoacids[R].codons[1].code[1]='G';
    aminoacids[R].codons[1].code[2]='G';
    aminoacids[R].codons[1].frequency=2175;
    aminoacids[R].codons[1].weight=0.222301717;

    aminoacids[R].codons[2].code[0]='C';
    aminoacids[R].codons[2].code[1]='G';
    aminoacids[R].codons[2].code[2]='A';
    aminoacids[R].codons[2].frequency=489;
    aminoacids[R].codons[2].weight=0.049979558;

    aminoacids[R].codons[3].code[0]='C';
    aminoacids[R].codons[3].code[1]='G';
    aminoacids[R].codons[3].code[2]='C';
    aminoacids[R].codons[3].frequency=658;
    aminoacids[R].codons[3].weight=0.067252657;

    aminoacids[R].codons[4].code[0]='C';
    aminoacids[R].codons[4].code[1]='G';
    aminoacids[R].codons[4].code[2]='G';
    aminoacids[R].codons[4].frequency=342;
    aminoacids[R].codons[4].weight=0.034955029;

    aminoacids[R].codons[5].code[0]='C';
    aminoacids[R].codons[5].code[1]='G';
    aminoacids[R].codons[5].code[2]='U';
    aminoacids[R].codons[5].frequency=3307;
    aminoacids[R].codons[5].weight=0.338000818;

    //---------------S
    aminoacids[S].name=S;
    aminoacids[S].num_codons=6;
    
	aminoacids[S].codons[0].code[0]='A';
    aminoacids[S].codons[0].code[1]='G';
    aminoacids[S].codons[0].code[2]='C';
    aminoacids[S].codons[0].frequency=2623;
    aminoacids[S].codons[0].weight=0.261645885;

    aminoacids[S].codons[1].code[0]='A';
    aminoacids[S].codons[1].code[1]='G';
    aminoacids[S].codons[1].code[2]='U';
    aminoacids[S].codons[1].frequency=3873;
    aminoacids[S].codons[1].weight=0.386334165;

    aminoacids[S].codons[2].code[0]='U';
    aminoacids[S].codons[2].code[1]='C';
    aminoacids[S].codons[2].code[2]='A';
    aminoacids[S].codons[2].frequency=4583;
    aminoacids[S].codons[2].weight=0.457157107;

    aminoacids[S].codons[3].code[0]='U';
    aminoacids[S].codons[3].code[1]='C';
    aminoacids[S].codons[3].code[2]='C';
    aminoacids[S].codons[3].frequency=6403;
    aminoacids[S].codons[3].weight=0.638703242;

    aminoacids[S].codons[4].code[0]='U';
    aminoacids[S].codons[4].code[1]='C';
    aminoacids[S].codons[4].code[2]='G';
    aminoacids[S].codons[4].frequency=2112;
    aminoacids[S].codons[4].weight=0.210673317;

    aminoacids[S].codons[5].code[0]='U';
    aminoacids[S].codons[5].code[1]='C';
    aminoacids[S].codons[5].code[2]='U';
    aminoacids[S].codons[5].frequency=10025;
    aminoacids[S].codons[5].weight=1;

    //---------------T
    aminoacids[T].name=T;
    aminoacids[T].num_codons=4;
    
	aminoacids[T].codons[0].code[0]='A';
    aminoacids[T].codons[0].code[1]='C';
    aminoacids[T].codons[0].code[2]='A';
    aminoacids[T].codons[0].frequency=5037;
    aminoacids[T].codons[0].weight=0.513350999;

    aminoacids[T].codons[1].code[0]='A';
    aminoacids[T].codons[1].code[1]='C';
    aminoacids[T].codons[1].code[2]='C';
    aminoacids[T].codons[1].frequency=6660;
    aminoacids[T].codons[1].weight=0.678760701;

    aminoacids[T].codons[2].code[0]='A';
    aminoacids[T].codons[2].code[1]='C';
    aminoacids[T].codons[2].code[2]='G';
    aminoacids[T].codons[2].frequency=1938;
    aminoacids[T].codons[2].weight=0.197513249;

    aminoacids[T].codons[3].code[0]='A';
    aminoacids[T].codons[3].code[1]='C';
    aminoacids[T].codons[3].code[2]='U';
    aminoacids[T].codons[3].frequency=9812;
    aminoacids[T].codons[3].weight=1;

    //---------------V
    aminoacids[V].name=V;
    aminoacids[V].num_codons=4;
    
	aminoacids[V].codons[0].code[0]='G';
    aminoacids[V].codons[0].code[1]='U';
    aminoacids[V].codons[0].code[2]='A';
    aminoacids[V].codons[0].frequency=3249;
    aminoacids[V].codons[0].weight=0.283953854;

    aminoacids[V].codons[1].code[0]='G';
    aminoacids[V].codons[1].code[1]='U';
    aminoacids[V].codons[1].code[2]='C';
    aminoacids[V].codons[1].frequency=6911;
    aminoacids[V].codons[1].weight=0.604002797;

    aminoacids[V].codons[2].code[0]='G';
    aminoacids[V].codons[2].code[1]='U';
    aminoacids[V].codons[2].code[2]='G';
    aminoacids[V].codons[2].frequency=3700;
    aminoacids[V].codons[2].weight=0.32337004;

    aminoacids[V].codons[3].code[0]='G';
    aminoacids[V].codons[3].code[1]='U';
    aminoacids[V].codons[3].code[2]='U';
    aminoacids[V].codons[3].frequency=11442;
    aminoacids[V].codons[3].weight=1;

    //---------------W
    aminoacids[W].name=W;
    aminoacids[W].num_codons=1;
    
	aminoacids[W].codons[0].code[0]='U';
    aminoacids[W].codons[0].code[1]='G';
    aminoacids[W].codons[0].code[2]='G';
    aminoacids[W].codons[0].frequency=3972;
    aminoacids[W].codons[0].weight=1;

    //---------------Y
    aminoacids[Y].name=Y;
    aminoacids[Y].num_codons=2;
    
	aminoacids[Y].codons[0].code[0]='U';
    aminoacids[Y].codons[0].code[1]='A';
    aminoacids[Y].codons[0].code[2]='C';
    aminoacids[Y].codons[0].frequency=7114;
    aminoacids[Y].codons[0].weight=1;

    aminoacids[Y].codons[1].code[0]='U';
    aminoacids[Y].codons[1].code[1]='A';
    aminoacids[Y].codons[1].code[2]='U';
    aminoacids[Y].codons[1].frequency=5768;
    aminoacids[Y].codons[1].weight=0.810795614;
}

bool cmpSort(LCP i, LCP j) {
    return i.lcp > j.lcp;
}

LCP *createLCP(std::string& s, int *SA) {
    int n = s.size(),
        *rank = new int[n + 1];
    LCP *lcp = new LCP[n + 1];

    for (int i = 0; i < n; i++) {
        rank[SA[i]] = i;
    }

    int i = 0, k = 0, h = 0;

    for (i = 0; i < n; i++) {
        if (h > 0) {
            h--;
        }

        k = SA[rank[i] - 1];

        while(s[i + h] == s[k + h] && (i + h) < n && (k + h) < n) {
            h++;
        }
        
        

        lcp[rank[i]].lcp = h;
        lcp[rank[i]].sa_pos = rank[i];
    }

    delete [] rank;
    return lcp;
}

struct CDSmaxLCS findLCS(std::string str1, std::string str2) {
    std::string str = str1 + '#' + str2, res;
    struct CDSmaxLCS maxLCS;
    int n = str.size(),
        *s = new int[n + 3],
        *SA = new int[n + 3],
        hash_pos = str1.size();

    memset(s, 0, (n + 3) * sizeof(int));
    memset(SA, 0, (n + 3) * sizeof(int));

    for(int i = 0; i < n; i++) {
        s[i] = str.at(i) == '#' ? 0 : (str.at(i) - 'a' + 1);
    }

    suffixArray(s, SA, n, 27);

    LCP *lcp = createLCP(str, SA);
    std::sort(lcp, lcp + n, cmpSort);

    int i;
    for(i = 0; i < n; i++) {
        if ((SA[lcp[i].sa_pos] < hash_pos && SA[lcp[i].sa_pos - 1] > hash_pos) || (SA[lcp[i].sa_pos] > hash_pos && SA[lcp[i].sa_pos - 1] < hash_pos)) {
            break;
        }
    }
    //std::cout << "index start" << lcp[i].sa_pos/8<<'\n';

    if (i < n) {
        for(int j = 0; j < lcp[i].lcp; j++) {
            res += str[SA[lcp[i].sa_pos] + j];
        }

        maxLCS.index=SA[lcp[i].sa_pos];
        maxLCS.maxsubstring=res;
//        std::cout << "str: " <<str<< "    ///   res  "<<res<<'\n';
//        std::cout << "lcp["<<i<<"].lcp " << lcp[i].lcp<<'\n';
//        std::cout << "lcp["<<i<<"].sa_pos " << SA[lcp[i].sa_pos]<<'\n';
    }

    delete [] lcp;
    delete [] SA;
    delete [] s;
    return maxLCS;
}

struct CDSmaxLCS lcs(std::vector<std::string> buf) {

    struct CDSmaxLCS Mlcs;

    // if (buf.size() == 1) {
    //     // search lcs for n-th string and previously found lcs
    //     lcs = findLCS(buf[0], buf[0]);
    //     buf.clear();
    //     buf.push_back(lcs);
    // }
    if (buf.size() > 1) {
        // search lcs for n-th string and previously found lcs
        Mlcs = findLCS(buf[0], buf[1]);
        buf.clear();
        //buf.push_back(Mlcslcs);
    }

    // if (!buf.size() || !buf[0].size()) {
    //     std::cout << "The longest common substring was not found\n";
    // } else {
    //     std::cout << buf[0] << "\n";
    // }

    //return buf[0];
    return Mlcs;
}

/*
* Simple Linear Work Suﬃx Array Construction algorithm

* Juha K¨arkk¨ainen and Peter Sanders
* Max-Planck-Institut f¨ur Informatik
* Stuhlsatzenhausweg 85, 66123 Saarbr¨ucken, Germany
* [juha,sanders]@mpi-sb.mpg.de.
* http://www.cs.helsinki.fi/u/tpkarkka/publications/icalp03.pdf
*/
inline bool leq(int a1, int a2,   int b1, int b2) { // lexic. order for pairs
  return(a1 < b1 || a1 == b1 && a2 <= b2);
}                                                   // and triples

inline bool leq(int a1, int a2, int a3,   int b1, int b2, int b3) {
  return(a1 < b1 || a1 == b1 && leq(a2,a3, b2,b3));
}

// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
static void radixPass(int* a, int* b, int* r, int n, int K){
  // count occurrences
  int* c = new int[K + 1];                          // counter array
  for (int i = 0;  i <= K;  i++) c[i] = 0;         // reset counters
  for (int i = 0;  i < n;  i++) c[r[a[i]]]++;    // count occurences
  for (int i = 0, sum = 0;  i <= K;  i++) { // exclusive prefix sums
     int t = c[i];  c[i] = sum;  sum += t;
  }
  for (int i = 0;  i < n;  i++) b[c[r[a[i]]]++] = a[i];      // sort
  delete [] c;
}

// find the suffix array SA of s[0..n-1] in {1..K}^n
// require s[n]=s[n+1]=s[n+2]=0, n>=2
void suffixArray(int* s, int* SA, int n, int K) {
  int n0=(n+2)/3, n1=(n+1)/3, n2=n/3, n02=n0+n2;
  int* s12  = new int[n02 + 3];  s12[n02]= s12[n02+1]= s12[n02+2]=0;
  int* SA12 = new int[n02 + 3]; SA12[n02]=SA12[n02+1]=SA12[n02+2]=0;
  int* s0   = new int[n0];
  int* SA0  = new int[n0];

  // generate positions of mod 1 and mod  2 suffixes
  // the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
  for (int i=0, j=0;  i < n+(n0-n1);  i++) if (i%3 != 0) s12[j++] = i;

  // lsb radix sort the mod 1 and mod 2 triples
  radixPass(s12 , SA12, s+2, n02, K);
  radixPass(SA12, s12 , s+1, n02, K);
  radixPass(s12 , SA12, s  , n02, K);

  // find lexicographic names of triples
  int name = 0, c0 = -1, c1 = -1, c2 = -1;
  for (int i = 0;  i < n02;  i++) {
    if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2) {
      name++;  c0 = s[SA12[i]];  c1 = s[SA12[i]+1];  c2 = s[SA12[i]+2];
    }
    if (SA12[i] % 3 == 1) { s12[SA12[i]/3]      = name; } // left half
    else                  { s12[SA12[i]/3 + n0] = name; } // right half
  }

  // recurse if names are not yet unique
  if (name < n02) {
    suffixArray(s12, SA12, n02, name);
    // store unique names in s12 using the suffix array
    for (int i = 0;  i < n02;  i++) s12[SA12[i]] = i + 1;
  } else // generate the suffix array of s12 directly
    for (int i = 0;  i < n02;  i++) SA12[s12[i] - 1] = i;

  // stably sort the mod 0 suffixes from SA12 by their first character
  for (int i=0, j=0;  i < n02;  i++) if (SA12[i] < n0) s0[j++] = 3*SA12[i];
  radixPass(s0, SA0, s, n0, K);

  // merge sorted SA0 suffixes and sorted SA12 suffixes
  for (int p=0,  t=n0-n1,  k=0;  k < n;  k++) {
#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
    int i = GetI(); // pos of current offset 12 suffix
    int j = SA0[p]; // pos of current offset 0  suffix
    if (SA12[t] < n0 ?
        leq(s[i],       s12[SA12[t] + n0], s[j],       s12[j/3]) :
        leq(s[i],s[i+1],s12[SA12[t]-n0+1], s[j],s[j+1],s12[j/3+n0]))
    { // suffix from SA12 is smaller
      SA[k] = i;  t++;
      if (t == n02) { // done --- only SA0 suffixes left
        for (k++;  p < n0;  p++, k++) SA[k] = SA0[p];
      }
    } else {
      SA[k] = j;  p++;
      if (p == n0)  { // done --- only SA12 suffixes left
        for (k++;  t < n02;  t++, k++) SA[k] = GetI();
      }
    }
  }
  delete [] s12; delete [] SA12; delete [] SA0; delete [] s0;
}

/* Routine for usual non-domination checking. It will return the
following values:
1 if a dominates b
-1 if b dominates a
0 if both a and b are non-dominated */
int dominate(individual a, individual b){
    int flag1, flag2;

    flag1 = flag2 = 0;

    if (a.results[0] > b.results[0])
        flag1 = 1;
    else if (a.results[0] < b.results[0])
        flag2 = 1;

    if (a.results[1] > b.results[1])
        flag1 = 1;
    else if (a.results[1] < b.results[1])
        flag2 = 1;

    if (a.results[2] < b.results[2])
        flag1 = 1;
    else if (a.results[2] > b.results[2])
        flag2 = 1;

    if (flag1 == 1 && flag2 == 0)
        return 1;
    if (flag1 == 0 && flag2 == 1)
        return -1;
    return 0;
}

/* Calculate lower and upper bounds for the different objective functions */
void calculate_bounds(long double &min_cai,long double &max_cai,long double &min_hd,long double &max_hd,long double &min_mlrcs,long double &max_mlrcs) {
	int i;
    
    min_cai=max_cai=individualList[0].results[0];
    min_hd=max_hd=individualList[0].results[1];
    min_mlrcs=max_mlrcs=individualList[0].results[2];

    for(i=1;i<individualList.size();i++){
        if (individualList[i].results[0] > max_cai) max_cai=individualList[i].results[0];
        else
            if (individualList[i].results[0] < min_cai) min_cai = individualList[i].results[0];
                
        if (individualList[i].results[1] > max_hd) max_hd=individualList[i].results[1];
        else
            if (individualList[i].results[1] < min_hd) min_hd= individualList[i].results[1];
        
        if (individualList[i].results[2] > max_mlrcs) max_mlrcs=individualList[i].results[2];
        else
            if (individualList[i].results[2] < min_mlrcs) min_mlrcs= individualList[i].results[2];     
    }
}

/* Comparison for sorting in ascending order of magnitude by the min CAI */
int compare_cai(const void * a, const void * b) {
	if (((struct individual *)a)->results[0] > ((struct individual *)b)->results[0]) return 1;
	if (((struct individual *)a)->results[0] < ((struct individual *)b)->results[0]) return -1;
	return 0;
}

/* Comparison for sorting in ascending order of magnitude by the min HD */
int compare_hd(const void * a, const void * b) {
	if (((struct individual *)a)->results[1] > ((struct individual *)b)->results[1]) return 1;
	if (((struct individual *)a)->results[1] < ((struct individual *)b)->results[1]) return -1;
	return 0;
}

/* Comparison for sorting in ascending order of magnitude by the MAX LRCS */
int compare_lrcs(const void * a, const void * b) {
	if (((struct individual *)a)->results[2] > ((struct individual *)b)->results[2]) return 1;
	if (((struct individual *)a)->results[2] < ((struct individual *)b)->results[2]) return -1;
	return 0;
}

/* Comparison for sorting in descending order of magnitude by the crowding distance */
int compare_crowding(const void * a, const void * b) {
	if (((struct individual *)a)->crowding > ((struct individual *)b)->crowding) return -1;
	if (((struct individual *)a)->crowding < ((struct individual *)b)->crowding) return 1;
	return 0;
}

/* Fast non-dominated sort and crowding distance assignment and sort */
void rank_and_crowding() {
    int size=individualList.size();
	int p,q,domination,i,max_front,ini;
	int n[size];
	int S_size[size];
	int S[size][size];
	int ranks_size[size+1];
	int rank_pos;
	int rank[size];
	struct individual copy_sol[size];
    long double min_cai,max_cai,min_hd,max_hd,min_mlrcs,max_mlrcs;

	/* Fast non-dominated sort */
	rank_pos=ranks_size[1]=0;
	for (p=0;p < size;p++) {
		S_size[p]=n[p]=0;
		for (q=0;q < size;q++) {
			domination=dominate(individualList[p],individualList[q]);
			if (domination == 1) {
				S[p][S_size[p]]=q;
				S_size[p]++;
			}
			else
				if (domination == -1) n[p]++;
		}

		if (n[p] == 0) {
			individualList[p].rank=1;
			rank[rank_pos]=p;
			rank_pos++;
			ranks_size[1]++;
		}
	}
	p=0;
	i=1;
	while (ranks_size[i] != 0) {
		ranks_size[i+1]=0;
		ini=p;
		for (;p < ini+ranks_size[i];p++) {
			int rp=rank[p];
			for (q=0;q < S_size[rp];q++) {
				n[S[rp][q]]--;
				if (n[S[rp][q]] == 0) {
					individualList[S[rp][q]].rank=i+1;
					rank[rank_pos]=S[rp][q];
					rank_pos++;
					ranks_size[i+1]++;
				}
			}
		}
		i++;
	}

	/* Crowding distance assignment */
	calculate_bounds(min_cai,max_cai,min_hd,max_hd,min_mlrcs,max_mlrcs);
	for (p=0;p < size;p++) {
		copy_sol[p]=individualList[rank[p]];
		copy_sol[p].crowding=0;
	}
	p=0;
	max_front=i;
	for (i=1;i < max_front;i++) {
		qsort(&copy_sol[p],ranks_size[i],sizeof(struct individual),compare_cai);
		copy_sol[p].crowding=copy_sol[p+ranks_size[i]-1].crowding=LDBL_MAX;
		ini=p-1;
		for (p++;p < ini+ranks_size[i];p++)
			copy_sol[p].crowding+=(copy_sol[p+1].results[0] - copy_sol[p-1].results[0])/(double)(max_cai - min_cai);

		p=ini+1;
		qsort(&copy_sol[p],ranks_size[i],sizeof(struct individual),compare_hd);
		copy_sol[p].crowding=copy_sol[p+ranks_size[i]-1].crowding=LDBL_MAX;
		for (p++;p < ini+ranks_size[i];p++)
			if (copy_sol[p].crowding != LDBL_MAX)
				copy_sol[p].crowding+=(copy_sol[p+1].results[1] - copy_sol[p-1].results[1])/(max_hd - min_hd);

        p=ini+1;
		qsort(&copy_sol[p],ranks_size[i],sizeof(struct individual),compare_lrcs);
		copy_sol[p].crowding=copy_sol[p+ranks_size[i]-1].crowding=LDBL_MAX;
		for (p++;p < ini+ranks_size[i];p++)
			if (copy_sol[p].crowding != LDBL_MAX)
				copy_sol[p].crowding+=(copy_sol[p+1].results[2] - copy_sol[p-1].results[2])/(max_mlrcs - min_mlrcs);
            
        //Sorting
		qsort(&copy_sol[ini+1],ranks_size[i],sizeof(struct individual),compare_crowding);
		p=ini+1+ranks_size[i];
	}
    
	for (p=0;p < size;p++){
        individualList[p].cont=copy_sol[p].cont;
        individualList[p].crowding=copy_sol[p].crowding;
        individualList[p].cub=copy_sol[p].cub;
        individualList[p].gender=copy_sol[p].gender;
        individualList[p].rank=copy_sol[p].rank;
        individualList[p].results[0]=copy_sol[p].results[0];
        individualList[p].results[1]=copy_sol[p].results[1];
        individualList[p].results[2]=copy_sol[p].results[2];
        for(int j=0;j<copy_sol[p].cds.size();j++){
            for(int t=0;t<copy_sol[p].cds[j].size();t++){
                individualList[p].cds[j][t] = copy_sol[p].cds[j][t].c_str();
            }    
        }
    } 
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

//Read a file
int readFromFile(string protein){
    std::string line;
    std::string campo[3];
    //std::ifstream file("file.txt");
    std::ifstream file(protein);
    if(file.is_open()){
        printf("Abierto\n");
        //getline (file,line); // lectura cabecera
        while (!file.eof()){
            for(int i=0;i<3;i++){
                getline (file,campo[i],' ');         
            }
            dts.id = campo[0].c_str(); 
            dts.chain = campo[1].c_str();
            dts.n= atoi(campo[2].c_str()); 
        }
        file.close();
    }
     else{
         printf("No existe fichero");
    } 
    printf("Procesado...\n");
    return 0;
}

//Cargar tabla de datos y lectura del fichero
void readAndLoadData(t_aminoacid *am,string protein){
    printf("Cargando tabla de datos...\n");
    loadtableAA(am);
    printf("Leyendo de fichero...\n");
    if(readFromFile(protein) != 0){
        printf("ERROR lectura de fichero\n");
    }
}

//Converts string to enum
t_nameaa convert(const std::string& str){ 
    if(str == "A") return A;
    else if(str == "C") return C;
    else if(str == "D") return D;
    else if(str == "E") return E;
    else if(str == "F") return F;
    else if(str == "G") return G;
    else if(str == "H") return H;
    else if(str == "I") return I;
    else if(str == "K") return K;
    else if(str == "L") return L;
    else if(str == "N") return N;
    else if(str == "M") return M;
    else if(str == "P") return P;
    else if(str == "Q") return Q;
    else if(str == "R") return R;
    else if(str == "S") return S;
    else if(str == "T") return T;
    else if(str == "V") return V;
    else if(str == "W") return W;
    else if(str == "Y") return Y;
}

//-----------------------------------------------------------------

//return num of synonyms
int numSynonyms(t_aminoacid *am,t_nameaa aminoacid){
    return am[aminoacid].num_codons;
}

bool compareCodon(char s1[3], char s2[3]){
    if(s1[0] == s2[0] && s1[1] == s2[1] && s1[2] == s2[2]){
        return true;
    }
    else{
        return false;
    }
}

bool checkCodons(string c1, string c2){
    bool check;
    char s1[3], s2[3];

    s1[0] = c1[0];
    s1[1] = c1[1];
    s1[2] = c1[2];

    s2[0] = c2[0];
    s2[1] = c2[1];
    s2[2] = c2[2];

    check = compareCodon(s1,s2);
        
    return check;  
}

string toLowerCase(std::string s){
    for(int i = 0;i< s.size(); i++){
        s[i] = tolower(s[i]);
    }
    return s;
}

string toUpperCase(std::string s){
    for(int i = 0;i< s.size(); i++){
        s[i] = toupper(s[i]);
    }
    return s;
}

//Converts a cds chain from uppercase to lowercase
std::vector<string> upperToLower(std::vector<string> cds){
    for(int i=0;i<cds.size();i++){
        std::for_each(cds[i].begin(), cds[i].end(), [](char & c){
        c = ::tolower(c);
        });
    }
    return cds;
}

//Returns the CDS vector in a single string 
string vectorToChain(std::vector<string> cds){
    std::string chain;
    chain.clear();
    for(auto s: cds){
        chain.append(s.c_str());
    }
  
    return chain;
}

//Converts a single string into a vector of strings with codons separated
std::vector<string> chainToVector(std::string chain){
    std::vector<string> vector;
    vector.clear();
    char aux[3] ;
    for(int i=0;i<chain.size();i=i+3){
        aux[0] = chain[i];
        aux[1] = chain[i+1];
        aux[2] = chain[i+2];
        vector.push_back(aux);
    }
    return vector;
}

//Converts the chain of aminoacids into a vector with the single aminoacid
std::vector<string> aminoChainToVector(std::string aminoChain){
    std::vector<string> vector;
    for(auto s : aminoChain){
        std::string str(1,s);
        vector.push_back(str); 
    }
    return vector;
}

string calcNewCDS(std::vector<vector<string>> cds_s,int idx,int changePos,std::string codon){
    int pos=0;
    string newCDS;
    for(auto s : cds_s){
        if(pos==idx){
            newCDS= vectorToChain(s);
            newCDS[changePos] = codon[0];
            newCDS[changePos+1] = codon[1];
            newCDS[changePos+2] = codon[2];
        }
        pos++;
    }  
    return newCDS;
}

void calculateNewCDSChain(std::vector<vector<string>> cds_s,std::string newCDS,int pos,std::vector<vector<string>> &newCDS_s){
    std::vector<string> vectorCDS = chainToVector(newCDS);
    newCDS_s.clear();
    int id=0;
    for(auto s : cds_s){
        if(id == pos){
            newCDS_s.push_back(vectorCDS);
        }
        else{
            newCDS_s.push_back(s);
        }   
        id++;
    }  
}

//Searches the weight of a codon
double searchWeight(t_aminoacid *am,t_nameaa amino,char s[3]){
    double w=0.0;
    int ncodons;
    for(int i=0;i<am[amino].num_codons;i++){
        if(am[amino].codons[i].code[0] == s[0] && am[amino].codons[i].code[1] == s[1] && am[amino].codons[i].code[2] == s[2]){
            w=am[amino].codons[i].weight;
        }
    }
    return w;
}

//Generates random Codon 
string randomGeneration(t_aminoacid *am,t_nameaa aminoacid){
    char c[3];
    string codon;
    int num,r1;
    num = am[aminoacid].num_codons;
    r1 = rand() % num; 
    c[0] = am[aminoacid].codons[r1].code[0];
    c[1] = am[aminoacid].codons[r1].code[1];
    c[2] = am[aminoacid].codons[r1].code[2];
    codon.push_back(c[0]);
    codon.push_back(c[1]);
    codon.push_back(c[2]);
    return codon;
}

//Generates random CDS for the given aminoacid
void randomCDSgenerator(t_aminoacid *am,std::string chain, std::vector<string> &cds){
    t_nameaa a; 
    string aux;
    for(int i=0;i<chain.size();i++){
        std::string str(1,dts.chain[i]);
        a = convert(str);
        aux = randomGeneration(am,a);
        cds.push_back(aux);
        aux.clear();    
    }
}

//Generates CDS with a perfect CAI. CAI=1
void perfectCAIgenerator(t_aminoacid *am,std::string chain,std::vector<string> &cds){
    t_nameaa a; 
    char c[3];
    string aux;
    bool exit =false;
    for(int i=0;i<chain.size();i++){
        std::string str(1,dts.chain[i]);
        a = convert(str);
        for(int j=0;j<am[a].num_codons && (!exit);j++){
            if(am[a].codons[j].weight == 1){
                c[0] = am[a].codons[j].code[0];
                c[1] = am[a].codons[j].code[1];
                c[2] = am[a].codons[j].code[2];
                aux.push_back(c[0]);
                aux.push_back(c[1]);
                aux.push_back(c[2]);
                cds.push_back(aux);     
                exit = true;
            }
        }
        aux.clear();
        exit = false;
    }
}

//Generates random codons for a chain of proteins. The number of repetitions is the value readed from file
void randomData(t_aminoacid *am){

    std::vector<string> cds;
    individual aux;
    
    for(int i=0;i<TAM-1;i++){
        aux.gender=-1;
        for(int j=0;j<dts.n;j++){
            randomCDSgenerator(am,dts.chain,cds);
            aux.cds.push_back(cds);  
            cds.clear();   
        }
        individualList.push_back(aux);
        aux.cds.clear(); 
    }
    perfectCAIgenerator(am,dts.chain,cds);
    aux.gender=-1;
    for(int k=0;k<dts.n;k++){
       aux.cds.push_back(cds);
    }
    individualList.push_back(aux);
    cds.clear();
    aux.cds.clear();
}

//----------------------------------------------------------------------------------
//----------------------------Multiobjetive Functions-------------------------------
//----------------------------------------------------------------------------------

//Calculates the total weight of a CDS
long double caiFunction(t_aminoacid *am, std::vector<string> aminoacids,std::vector<string> cds){

    long double result=1.0;
    char s[3];
    t_nameaa prot; 
    for(int i=0;i<aminoacids.size();i++){
        t_nameaa prot = convert(aminoacids[i]);
        s[0]=cds[i][0];
        s[1]=cds[i][1];
        s[2]=cds[i][2];
        result = result * searchWeight(am,prot,s); 
    } 
    return pow(result,(1.0/cds.size()));
}

//Calculates minimum codon adaptation index (CAI)
long double minCai(t_aminoacid *am,std::vector<string> aminoacids,std::vector<vector<string>> cds_s,std::vector<string> &best_cds,int &index){
    long double min = 1.0;
    long double function=0.0;
    std::vector<string> aux;
    int idx =0;
    for(auto s : cds_s){
        function = caiFunction(am,aminoacids,s);
        if(function < min){
            aux.clear();
            min = function;
            aux = s;
            index = idx;
        }
        idx++;
    }
    best_cds = aux;
    return min;
}

//Calculates the total of different nucleotids. If equals = 0.
long double hammingDistance(std::string cds1, std::string cds2){
    long double res=0.0;
    for(int i=0;i<cds1.size();i++){  
        if(cds1[i] != cds2[i]){
            res++;
        }  
    }
    return res;
}

//Calculates the minimun value between pairs
long double minHamming(std::vector<vector<string>> cds_c,std::vector<string> &bestPair,int &index1, int &index2){
    
    long double min = 10000.0; //big number for default
    long double hamming=0.0;
    int count=0;
    std::string chain1,chain2;
    for (int i=0;i<cds_c.size();i++){
        for(int j=i+1;j<cds_c.size();j++){
            chain1 = vectorToChain(cds_c[i]);
            chain2 = vectorToChain(cds_c[j]);
            hamming = hammingDistance(chain1,chain2);
            if(hamming < min){
                min = hamming;
                bestPair.clear();
                bestPair.push_back(vectorToChain(cds_c[i]));
                bestPair.push_back(vectorToChain(cds_c[j]));
                index1 = i;         //indices para conocer la posicion en el vector de los pares
                index2 = j;
            }
            count++;
        }
    }
    return min/(cds_c[0].size()*3);
}

//Returns max LRCS with lcs.cpp
long double maxLRCS(std::vector<string> cds,std::string &subStr,int &posIni, int &posFin,int &index,std::string &out){
    long double max = 0;
    std::string maxSubS;
    struct CDSmaxLCS cdsmaxlcs;
    int idx,cds_size = cds[0].size();
    std::vector<std::string> aux;
    cds = upperToLower(cds);
    for (int i=0;i<cds.size();i++){
        for(int j=i+1;j<cds.size();j++){
            aux.push_back(cds[i]);
            aux.push_back(cds[j]);
            cdsmaxlcs = lcs(aux); 
            if(cdsmaxlcs.maxsubstring.size() > max){
                max = cdsmaxlcs.maxsubstring.size();
                subStr = cdsmaxlcs.maxsubstring;
                idx = cdsmaxlcs.index%cds_size;
                posIni = idx;
                posFin = idx+max-1;
                out = cds[i];
                index = i;
            }
            aux.clear();
        }
    }
    return max;
}

//Returns max LRCS with lrs.cpp
/*long double maxLRCS(std::vector<string> cds,std::string &subStr,int &posIni, int &posFin,int &index,std::string &out){
    long double max = 0;
    std::string maxSubS;
    struct CDSmaxLRS cdsmaxlrs;
    int idx,cds_size = cds[0].size();
    std::string aux;
    //std::vector<std::string> aux;
    cds = upperToLower(cds);
    for (int i=0;i<cds.size();i++){
        aux = cds[i];
        cdsmaxlrs = lrs(cds[i]);
        if(cdsmaxlrs.maxsubstring.size() > max){
            max = cdsmaxlrs.maxsubstring.size();
            subStr = cdsmaxlrs.maxsubstring;
            idx = cdsmaxlrs.index%cds_size;
            posIni = idx;
            posFin = idx+max-1;
            out = cds[i];
            index = i;
        }
        aux.clear();    
    }
    return max;
}*/

//length of repeated or common substrings
long double lrcsFunction(std::vector<vector<string>> cds_s,std::string &subStr,int &posIni, int &posFin,int &index,std::string &out){
    long double max =0.0;
    std::vector<string> chain;
    for(auto s : cds_s){
        chain.push_back(vectorToChain(s));
    }
    max = maxLRCS(chain,subStr,posIni,posFin,index,out);
    return (max/(cds_s[0].size()*3));
}

//----------------------------------------------------------------------------------
//----------------------------------Mutations---------------------------------------
//----------------------------------------------------------------------------------

//Selects new random codon different to the given condon
string newRandCodon(t_aminoacid *am,t_nameaa aminoacid,string codon){
    string newcodon;
    bool exit = false;
    while(!exit){
        newcodon = randomGeneration(am,aminoacid);
        if(!checkCodons(codon,newcodon)){           //Comprueba qe sean distintos, si lo son, exit
            exit=true;
        }
    }  
    return newcodon;
}

//Generates new codon for CAI mutation
string newmCaiCodon(t_aminoacid *am,t_nameaa aminoacid,string codon){
    string newcodon;
    bool exit = false;
    char s1[3], s2[3];
    s1[0] = codon[0];
    s1[1] = codon[1];
    s1[2] = codon[2];
    while(!exit){
        newcodon = randomGeneration(am,aminoacid);
        s2[0] = newcodon[0];
        s2[1] = newcodon[1];
        s2[2] = newcodon[2];
        if(searchWeight(am,aminoacid,s2) > searchWeight(am,aminoacid,s1)){  // si nuevo w es mayor que el antiguo, exit
            exit = true;
        }
    }
    return newcodon;
}

//Random mutation
void randomMutation(t_aminoacid *am,std::vector<string> aminoacids,std::vector<vector<string>> cds_s,int pm,std::vector<vector<string>> &cds_out){ 
    t_nameaa p; 
    string nCodon;
    int r1;
    for(auto s : cds_s){
        for(int i=0;i<s.size();i++){
            p = convert(aminoacids[i]);
            r1 = rand() % 100; 
            if((r1 <= pm) && (numSynonyms(am,p)>1)){
                nCodon = newRandCodon(am,p,s[i]);       
                s[i][0] = nCodon[0];
                s[i][1] = nCodon[1];
                s[i][2] = nCodon[2];
            }
        }
        cds_out.push_back(s);
    }
}

//CAI mutation
void mCaiMutation(t_aminoacid *am,std::vector<string> aminoacids,std::vector<vector<string>> cds_s,int pm,std::vector<vector<string>> &cds_out){
    t_nameaa p; 
    int r1,pos,idx=0;
    char c[3];
    string nCodon;
    std::vector<string> m_cds;
    std::vector<string> aux;
    long double mcai = minCai(am,aminoacids,cds_s,m_cds,pos);
    for(int i=0;i<m_cds.size();i++){
        p = convert(aminoacids[i]);
        r1 = rand() % 100; 
        c[0] = m_cds[i][0];
        c[1] = m_cds[i][1];
        c[2] = m_cds[i][2];
        if((r1 <= pm) &&  searchWeight(am,p,c) != 1){
            nCodon = newmCaiCodon(am,p,m_cds[i]);
            aux.push_back(nCodon.c_str());
        }
        else{
            aux.push_back(m_cds[i]);
        }
    }
    //Copiar al vector salida
    for(auto s : cds_s){
        if(idx == pos){
            cds_out.push_back(aux);
        }
        else{
            cds_out.push_back(s);
        }   
        idx++;
    }
    m_cds.clear();
    aux.clear();
}

//Hamming mutation
void mHDMutation(t_aminoacid *am,std::vector<string> aminoacids,std::vector<vector<string>> cds_s,int pm,std::vector<vector<string>> &cds_out){
    t_nameaa p; 
    int r1,idx1,idx2,idaux,idaux2,synonyms,pos=0;
    bool mutated = false;
    std::vector<string> pair,pairAux,mcAux;
    std::string newCDS,mutatedCodon;
    string c1,c2,bestSyn;
    std::vector<vector<string>> newCDS_s,mSolution;
    long double best_HDpair,best_mHD,new_hd_pair,new_mHD;
    long double c_mH = minHamming(cds_s,pair,idx1,idx2);             //Obtener min hamming y pares
    long double c_hd_pair;  
    for(int i=0;i<pair[0].size();i=i+3){
        p = convert(aminoacids[i/3]);
        r1 = rand() % 100;   
        synonyms = numSynonyms(am,p);     
        if((r1 <= pm) && (synonyms>1)){                         //Comprobar si es posible
            c_hd_pair = hammingDistance(pair[0],pair[1]);       //Calcular distancia hamming entre pares
            c_mH = minHamming(cds_s,pairAux,idaux,idaux2);
            best_HDpair = -1;
            best_mHD = -1;
            for(int j=0;j<synonyms;j++){
                c1[0] = am[p].codons[j].code[0];
                c1[1] = am[p].codons[j].code[1];
                c1[2] = am[p].codons[j].code[2];
                c2[0] = pair[0][i];
                c2[1] = pair[0][i+1];
                c2[2] = pair[0][i+2];
                if(!checkCodons(c1,c2)){
                    newCDS = calcNewCDS(cds_s,idx1,i,c1);
                    new_hd_pair = hammingDistance(newCDS,pair[1]);
                    calculateNewCDSChain(cds_s,newCDS,idx1,newCDS_s);
                    new_mHD = minHamming(newCDS_s,pairAux,idaux,idaux2);
                    if((new_mHD > c_mH) && (new_mHD > best_mHD)){
                        best_mHD = new_mHD;
                        bestSyn[0] = am[p].codons[j].code[0];
                        bestSyn[1] = am[p].codons[j].code[1];
                        bestSyn[2] = am[p].codons[j].code[2];
                    }
                    else{
                        if((best_mHD == -1) && (new_mHD == c_mH) && (new_hd_pair > c_hd_pair) && (new_hd_pair > best_HDpair)){
                            best_HDpair = new_hd_pair;
                            bestSyn[0] = am[p].codons[j].code[0];
                            bestSyn[1] = am[p].codons[j].code[1];
                            bestSyn[2] = am[p].codons[j].code[2];
                        }
                    }
                }
            }
            if((best_mHD != -1) || (best_HDpair != -1)){        //Cambiar solucion
                //printf("bestSyn: %s\n",bestSyn.c_str());
                mutatedCodon = calcNewCDS(cds_s,idx1,i,bestSyn);
                calculateNewCDSChain(cds_s,mutatedCodon,idx1,mSolution);
                mutated=true;
            }
        }
    }
    //Si no muta, se queda igual. Si muta, cambia por la solucion mutada
    if(!mutated){
        cds_out=cds_s;
    }
    else{
        cds_out=mSolution;
    }
}

//LRCS mutation
void lrcsMutation(t_aminoacid *am,std::vector<string> aminoacids,std::vector<vector<string>> cds_s,int pm,std::vector<vector<string>> &cds_out){
    t_nameaa p; 
    int r1,synonyms,idx,idxaux,posI,posF,posIaux,posFaux;
    bool mutated = false;
    std::string subStr,subStraux,cds,cds2,cdsaux,codon,nCodon;
    std::vector<vector<string>> newcds_s;
    long double mlrcs = lrcsFunction(cds_s,subStr,posI,posF,idx,cds);
    long double newmlrcs;
    cds = toUpperCase(cds);
    subStr = toUpperCase(subStr);
    cds2=cds;
    //calcular substring adaptado
    if(posI%3 == 1){
        posI--;
    }
    if(posI%3 == 2){
        posI=posI-2;
    }
    if(posF%3 == 0){
        posF =posF+2;
    }
    if(posF%3==1){
        posF++;
    }
    string subStrAdapt = cds.substr(posI,(posF+1)-posI);
    for(int i=0;i<subStrAdapt.size();i=i+3){ 
        r1 = rand() % 100; 
        p = convert(aminoacids[(posI/3)+(i/3)]);
        if((r1 <= pm) && (numSynonyms(am,p)>1)){
            codon = subStrAdapt.substr(i,3);
            nCodon = newRandCodon(am,p,codon);    
            cds2[posI+i] = nCodon[0];
            cds2[posI+i+1] = nCodon[1];
            cds2[posI+i+2] = nCodon[2];
            calculateNewCDSChain(cds_s,cds2,idx,newcds_s);
            newmlrcs = lrcsFunction(newcds_s,subStraux,posIaux,posFaux,idxaux,cdsaux);
            if(newmlrcs <= mlrcs){
                calculateNewCDSChain(cds_s,cds2,idx,newcds_s);
                mlrcs = newmlrcs;
                mutated = true;    
            }
        }
    }
    if(!mutated){
        cds_out=cds_s;
    }
    else{
        cds_out = newcds_s;
    } 
}

//----------------------------------------------------------------------------------
//-----------------------------------Compute----------------------------------------
//----------------------------------------------------------------------------------

//Distributes gender into population
void genderDistribution(){
    int r,females = FEMALES;
    //Inicializar a valor por defecto
    for(auto d=0;d<individualList.size();d++){
        individualList[d].gender = -1;
    }
    //random females
    for(int i=0;i<individualList.size() && females >0;i++){
        r = rand() % 2;
        if(r==1 && individualList[i].gender != 1){
            individualList[i].gender = r;
            females--;
        }
        if(i==individualList.size()){
            i=0;
        }
    }
    //males
    for(int j=0;j<individualList.size();j++){
        if(individualList[j].gender == -1){
            individualList[j].gender = 0;
        }
    }
}

//Distributes gender into population
//Better distribution for low females population
void genderDistribution2(){
    int r,females = FEMALES;
    //Inicializar a valor por defecto
    for(auto d=0;d<individualList.size();d++){
        individualList[d].gender = -1;
    }
    //females
    while(females > 0){
        r = rand() % 100;
        if(individualList[r].gender == -1){
            individualList[r].gender = 1;
            females--;
        }
    }
    //males
    for(int j=0;j<individualList.size();j++){
        if(individualList[j].gender == -1){
            individualList[j].gender = 0;
        }
    }
}

individual calcObjetives(t_aminoacid *am,std::vector<string> ams,std::vector<vector<string>> cds,individual actual){
    long double mC =0.0,mHD=0.0,MLRCS=0.0;
    std::vector<string> aux;
    std::string subS,output;
    int pos,idx1,idx2,idx3,posIni,posFin;
    individual ind;
    mC = minCai(am,ams,cds,aux,pos);
    mHD = minHamming(cds,aux,idx1,idx2);
    MLRCS = lrcsFunction(cds,subS,posIni,posFin,idx3,output);
    ind.results.push_back(mC);
    ind.results.push_back(mHD);
    ind.results.push_back(MLRCS);
    ind.gender = actual.gender;
    ind.cont = actual.cont;
    ind.cds = cds;
    aux.clear();
    return ind;
}

void calcMultiOjetives(t_aminoacid *am,std::vector<string> ams){
    long double mC =0.0,mHD=0.0,MLRCS=0.0;
    int cont=0;
    std::vector<string> aux;
    std::string subS,output;
    for(int d=0;d<individualList.size();d++){  
        individualList[d].results.clear(); 
        int pos,idx1,idx2,idx3,posIni,posFin;
        mC = minCai(am,ams,individualList[d].cds,aux,pos);
        mHD = minHamming(individualList[d].cds,aux,idx1,idx2);
        MLRCS = lrcsFunction(individualList[d].cds,subS,posIni,posFin,idx3,output);
        individualList[d].results.push_back(mC);
        individualList[d].results.push_back(mHD);
        individualList[d].results.push_back(MLRCS);
        aux.clear();
        subS.clear();
        output.clear();
        cont++;
    }
}

//Copy the rejected solutions to pareto vector, keeps the best solutions
//Filter the added solutions
void saveRejectedSolutions(){
    int equal;
    for(int j=TAM;j<individualList.size();j++){
        if(individualList[j].rank==1){
            paretoFrontier.push_back(individualList[j]);
        } 
        else{
            equal=0;
            for (int k=0;k < paretoFrontier.size();k++){
                if((individualList[j].results[0] == paretoFrontier[k].results[0]) && (individualList[j].results[1] == paretoFrontier[k].results[1]) && (individualList[j].results[2] == paretoFrontier[k].results[2]) ){
                    equal=1;
                    break;
                }
            }
            if(equal==0){
                paretoFrontier.push_back(individualList[j]);
            }
        } 
    }
    for(int i=0;i<TAM;i++){
        individualList.pop_back();
    }  
}

//Computes mutations
void mutations(t_aminoacid *am,int P,int &borns){
    int rF,rM,rM2,dmnt=0;
    int cont=0;                                 //Contador de mutaciones
    bool mutate = false;
    std::vector<vector<string>> cds_Out;
    std::vector<string> aminoacids;
    individual aux;
    std::vector<vector<string>> cub_cds_s;      //Lista cdss auxiliar para los cachorros
    std::vector<string> cub_cds;                //cds auxiliad para los cachorros
    individual cub,teen;                             //cachorros, adolescentes
    aminoacids = aminoChainToVector(dts.chain);

    for(auto d=0;d<individualList.size();d++){
        mutate = false;
        cds_Out.clear();
        if(individualList[d].gender == 1){                      //Si es hembra
            rF = rand() % 3+1;                 //Random para hembra [1-3]
            switch(rF){
                case 1:
                    //printf("mCAI\n");
                    mCaiMutation(am,aminoacids,individualList[d].cds,P,cds_Out);
                    mutate = true;
                    break;
                case 2:
                    //printf("mHD\n");
                    mHDMutation(am,aminoacids,individualList[d].cds,P,cds_Out);
                    mutate = true;
                    break;
                case 3:
                    //printf("MLRCS\n");
                    lrcsMutation(am,aminoacids,individualList[d].cds,P,cds_Out);
                    mutate = true;
                    break;
                default:
                    break;
            }  
        }
        if(individualList[d].gender == 0){                      //Si es macho
            //0,1-> macho adulto, 2-> macho adolescente
            rM = rand() % 2;                    //numero random para machos adolescentes
            switch(rM){
                case 0:
                    rM2 = rand() % 3+1;
                    switch (rM2){
                    case 1:
                        //printf("mCAI\n");
                        mCaiMutation(am,aminoacids,individualList[d].cds,P,cds_Out);
                        mutate = true;
                        break;
                    case 2:
                        //printf("mHD\n");
                        mHDMutation(am,aminoacids,individualList[d].cds,P,cds_Out);
                        mutate = true;
                        break;
                    case 3:
                        //printf("MLRCS\n");
                        lrcsMutation(am,aminoacids,individualList[d].cds,P,cds_Out);
                        mutate = true;
                        break;
                    default:
                        break;
                    }
                    break;
                case 1:
                    randomMutation(am,aminoacids,individualList[d].cds,P,cds_Out);
                    mutate = true;
                    break;
                default:
                    break;
            }
        }
        //Si hay mutacion, comprobar pareto
        if(mutate){
            aux = calcObjetives(am,aminoacids,cds_Out,individualList[d]); 
            dmnt = dominate(individualList[d],aux);
            if(dmnt == 1){
                //printf("Mejor la actual\n");
                individualList[d].cont++;
            }
            if(dmnt == -1){
                //printf("Mejor la nueva\n");
                individualList[d].cont=0;
                for(int j=0;j<cds_Out.size();j++){
                    for(int t=0;t<cds_Out[j].size();t++){
                        individualList[d].cds[j][t] = cds_Out[j][t].c_str();
                    }
                } 
            }
            else{ //Añadir adolescente
                //printf("No dominantes\n");   
                teen.cds.clear();
                if(individualList.size() < 200){   
                    teen.cont=0;
                    teen.cub=0;
                    teen.results.push_back(0.0);
                    teen.results.push_back(0.0);
                    teen.results.push_back(0.0);
                    teen.gender=individualList[d].gender; 
                    for(auto c : cds_Out){
                        teen.cds.push_back(c);
                    }
                    individualList.push_back(teen);
                }
            }  
        }
        //Si llega a 10, muere y genero nueva (cachorro que nace)
        if(individualList[d].cont == AGE){
            paretoFrontier.push_back(individualList[d]); //Añadir solucion perdida
            for(int i=0;i<dts.n;i++){
                randomCDSgenerator(am,dts.chain,cub_cds);
                cub_cds_s.push_back(cub_cds);  
                cub_cds.clear(); 
            }  
            individualList[d].cont=0;                   //Reinicia contador  
            individualList[d].cub=1;                
            for(int j=0;j<cub.cds.size();j++){
                for(int t=0;t<cub.cds[j].size();t++){
                    individualList[d].cds[j][t] = cub.cds[j][t].c_str();
                }
            }     
            borns++;
            cub_cds_s.clear();
        }     
    } 
}

//Compute calculations
void compute(t_aminoacid *am,std::vector<string> ams,int reps,int P,int &borns){
 
    for(int i=0;i<reps;i++){
        mutations(am,P,borns);      //Mutacion
        calcMultiOjetives(am,ams);  //Calcular multiobjetivos
        rank_and_crowding();        //ranking
        saveRejectedSolutions();    //Guardar soluciones deshechadas
        genderDistribution2();      //Distribucion del genero para nueva iteracion
    }
    for(auto p : individualList){   //Copiar al vector de Pareto las soluciones a final
        paretoFrontier.push_back(p);
    }
}

/* Remove the dominated or equal solutions from the file, and calculate the hypervolume */
void filter() {
    int i,j,k,dominated,dmn,equal;
	for (i=0;i < paretoFrontier.size();i++) {
		dominated=0;
		for (j=0;j < paretoFrontier.size();j++){ /* Check if solution i is dominated */
        dmn = dominate(paretoFrontier[i],paretoFrontier[j]);
			if (dmn == -1) {
				dominated=1;
				break;
			}
        }
		if (dominated == 0) {
			equal=0;
			for (j=0;j < paretoFrontierReduced.size();j++){ /* Check if solution i already exists */
				if ((paretoFrontier[i].results[0] == paretoFrontierReduced[j].results[0]) && (paretoFrontier[i].results[1] == paretoFrontierReduced[j].results[1]) && (paretoFrontier[i].results[2] == paretoFrontierReduced[j].results[2])) {
                    equal=1;
					break;
                }                
            }
			if (equal == 0) {
                paretoFrontierReduced.push_back(paretoFrontier[i]);
			}
		}
	}
}

void computeHV(char hv_sol_file_name[20]){
	FILE* hv_sol_file;
	//char hv_sol_file_name[20];
	FILE* better_sol_file;
	char better_sol_file_name[20];
	FILE* hv_value_file;
	char hv_value_file_name[20];
	char hv_command[70];
	double hypervolume;
	FILE* all_hv_file;

    /* Calculate the hypervolume */
    sprintf(hv_value_file_name,"HV_Results.txt");
	sprintf(hv_command,"./hyp_ind hyp_ind_param.txt %s reference_set %s",hv_sol_file_name,hv_value_file_name);
	system(hv_command);
	remove(hv_sol_file_name); /* After using it, remove the file with all the normalized solutions used for computing the hypervolume */
    //Abrir fichero en modo lectura
    hv_value_file=fopen(hv_value_file_name,"r");
	if (hv_value_file == NULL) {
		printf("ERROR: fopen failed to open %s\n",hv_value_file_name);
		exit(-1);
	}
    //leer valor del fichero y guardar en hypervolume
    fscanf(hv_value_file,"%le\n",&hypervolume);
	fclose(hv_value_file);
    //Crear fichero para guardar los resultados
    all_hv_file=fopen("All_HV.txt","a");
	if (all_hv_file == NULL) {
		printf("ERROR: fopen failed to open All_HV.txt\n");
		exit(-1);
	}
    //Guarda en el fichero el id de la proteina y su valor de hipervolumen
    fprintf(all_hv_file,"%s: %lf %%\n",dts.id.c_str(),-hypervolume*100);
	fclose(all_hv_file);
	remove(hv_value_file_name); /* After using it, remove the file with the HV value of this repetition */
}

//Save into file
void saveToFile(string filename){
    int cont=0;
    ofstream fsal;
    fsal.open(filename.c_str());
    if (fsal.is_open()){
        fsal << "ID proteina: " << dts.id << " Cadena: " << dts.chain << " Número de CDSs: " << dts.n << endl;
        fsal << endl;
        for(int d=0;d<individualList.size();d++){
            for(int i=0;i<individualList[d].cds.size();i++){
                fsal << "CDS" << i << ": ";
                for(int j=0;j<individualList[d].cds[i].size();j++){
                    fsal << individualList[d].cds[i][j].c_str();
                }
                fsal << endl;
            }
            fsal << "Individual "<< cont << " || mCAI: " << individualList[d].results[0] << " mHamming: " << individualList[d].results[1] 
            << " LRCS: " << individualList[d].results[2] << " Gender: " << individualList[d].gender << " Cub: " << individualList[d].cub 
            << " Cont: " << individualList[d].cont << " Rank: "<< individualList[d].rank <<endl;
            fsal << endl;
            cont++;
            }
        fsal.close();
    }
}

//Save into file Pareto individuals
void saveFilePareto(string filename){
    int cont=0;
     ofstream fsal;
     fsal.open(filename.c_str());
    if (fsal.is_open()){
        fsal << "ID proteina: " << dts.id << " Cadena: " << dts.chain << " Número de CDSs: " << dts.n << endl;
        fsal << endl;
        for(int i=0;i<paretoFrontierReduced.size();i++){
            for(int j=0;j<paretoFrontierReduced[i].cds.size();j++){
                fsal << "CDS" << j << ": ";
                for(int k=0;k<paretoFrontierReduced[i].cds[j].size();k++){
                    fsal << paretoFrontierReduced[i].cds[j][k].c_str();
                }
                fsal << endl;
            }
            fsal << "Individual "<< cont << " || mCAI: " << paretoFrontierReduced[i].results[0] << " mHamming: " << paretoFrontierReduced[i].results[1] 
            << " LRCS: " << paretoFrontierReduced[i].results[2] << " Gender: " << paretoFrontierReduced[i].gender << " Cub: " << paretoFrontierReduced[i].cub 
            << " Cont: " << paretoFrontierReduced[i].cont << " Rank: "<< paretoFrontierReduced[i].rank <<endl;
            fsal << endl;
            cont++;
        }
        fsal.close();
    }
}

void saveFileParetoFrontier(string filename){
    int cont=0;
    ofstream fsal;
    fsal.open(filename.c_str());
    if (fsal.is_open()){
        fsal << "ID proteina: " << dts.id << " Cadena: " << dts.chain << " Número de CDSs: " << dts.n << endl;
        fsal << endl;
        for(int i=0;i<paretoFrontier.size();i++){
            for(int j=0;j<paretoFrontier[i].cds.size();j++){
                fsal << "CDS" << j << ": ";
                for(int k=0;k<paretoFrontier[i].cds[j].size();k++){
                    fsal << paretoFrontier[i].cds[j][k].c_str();
                }
                fsal << endl;
            }
            fsal << "Individual "<< cont << " || mCAI: " << paretoFrontier[i].results[0] << " mHamming: " << paretoFrontier[i].results[1] 
            << " LRCS: " << paretoFrontier[i].results[2] << " Gender: " << paretoFrontier[i].gender << " Cub: " << paretoFrontier[i].cub 
            << " Cont: " << paretoFrontier[i].cont << " Rank: "<< paretoFrontier[i].rank <<endl;
            fsal << endl;
            cont++;
        }
        fsal.close();
    }
}

//File for HV  
//Hamming distance divided by 0.4
void saveToHVfile(string filename){
    ofstream fsal;
    fsal.open(filename.c_str());
    if (fsal.is_open()){
        for(int d=0;d<paretoFrontierReduced.size();d++){
            fsal << paretoFrontierReduced[d].results[0] << " " << paretoFrontierReduced[d].results[1]/0.4 << " " << paretoFrontierReduced[d].results[2] << endl;
        }
        fsal.close();
    }
}

//Save info into log.txt
void saveLog(int exp,int reps, double tTotal){
    double parteDecimalC,parteEnteraC,parteDecimalF,parteEnteraF,parteDecimal2C,parteEntera2C,
    parteDecimal2F,parteEntera2F,parteDecimalT,parteEnteraT,parteDecimalT2,parteEnteraT2;
    int p1C,p2C,p1F,p2F,pt,pt2,cont;
    string filename = "log.txt";
    ofstream fsal;
    fsal.open(filename.c_str());
    if (fsal.is_open()){ 
        fsal << "ID Proteina: " << dts.id << endl << "Cadena: "  <<dts.chain << endl << "Número de CDSs: " << dts.n << endl;
        fsal << endl;
        fsal << "Numero de experimentos: " << exp << endl;
        fsal << "Numero de repeticiones: " << reps << " | Porcentaje de mutacion: " << PM << "%" << " | Numero Hembras inicial: " << FEMALES << " | Edad maxima: " << AGE << endl;
        fsal << endl;
        for (int i =0;i<exp;i++){
            //Formatear de segundos a horas
            parteDecimalC = modf((logF[i].tsecC/60),&parteEnteraC);     //Minutos
            parteDecimal2C = modf((logF[i].tsecC/60)/60,&parteEntera2C);//Horas
            p1C = (parteDecimalC*0.6)*100;                              //Segundos
            p2C = (parteDecimal2C*0.6)*100;                             //Minutos
            parteDecimalF = modf((logF[i].tsecF/60),&parteEnteraF);     
            parteDecimal2F = modf((logF[i].tsecF/60)/60,&parteEntera2F);
            p1F = (parteDecimalF*0.6)*100;                              
            p2F = (parteDecimal2F*0.6)*100;  
            fsal << "Iteracion: " << logF[i].id << " --> N. nacimientos: " << logF[i].newBorns << " | N. Hembras: " << logF[i].females 
            << " | N. Machos: " << logF[i].males << " | Individuos en Pareto: " << logF[i].tPareto <<  endl;
            fsal << "Tiempo en horas Calculos: " << parteEntera2C << "h " << p2C << "m " << p1C << "s" << endl; 
            fsal << "Tiempo en horas Filtrado: " << parteEntera2F << "h " << p2F << "m " << p1F << "s" << endl;
            fsal << endl;
        }
        fsal << endl;
        parteDecimalT = modf((tTotal/60),&parteEnteraT);
        parteDecimalT2= modf(((tTotal/60)/60),&parteEnteraT2);
        pt = (parteDecimalT*0.6)*100;
        pt2 = (parteDecimalT2*0.6)*100;
        //fsal << "Tiempo total del experimento en minutos: " << parteEnteraT << "m " << pt << "s" << endl;
        fsal << "Tiempo total del experimento en horas: " << parteEnteraT2 << "h " << pt2 << "m " << pt << "s" << endl;
        fsal.close();
    }
}

//Cleans all data structures for each new iteration
void flushStructures(){
    individualList.clear();
    paretoFrontier.clear();
    paretoFrontierReduced.clear();
}

//Show information on terminal
void showInformation(){
    int c=0;
    for(auto d : individualList){
        printf("Num: %d, Cont: %d, Gender: %d\n",c,d.cont,d.gender);
        for(auto s : d.cds){
            printf("CDS: ");
            for(auto c : s){
                printf("%s",c.c_str());
            }
            printf("\n");
        } 
        printf("mCAI: %Lf, mHD: % Lf, MLRCS: %Lf\n",d.results[0],d.results[1],d.results[2]);
        printf("\n");
        c++;
    }
}

//Print the protein data on terminal
void printData(){    
    printf("ID: %s\n",dts.id.c_str());
    printf("CHAIN: %s\n",dts.chain.c_str());
    printf("N CDS: %d\n",dts.n);   
}

// argv[1] -> numero de repeticiones (epocas)
// argv[2] -> numero de experimentos
// argv[3] -> proteina
int main(int argc, char const *argv[]) 
{   
    t_aminoacid *am;
    am=(t_aminoacid *)malloc(sizeof(t_aminoacid)*NN);
    srand(time(NULL));
    string protein;
    std::vector<string> aminoacids;
    double total=0.0;
    time_t start, end,t_start,t_end;
    char hv_filename[20];
    int nBorns,m,f;
    int pm = PM;

    if(argc < 1 || argc > 4){
	 	printf("Error al iniciar. Indique valores de entrada correctos\n");
        printf("N repeticiones, N experimentos, Proteina\n");
	}
	else{
        switch(atoi(argv[3])){                  //Eleccion de la proteina a procesar indicada por parametro
            case 1:
                protein = "B7KHU9.txt";
                break;
            case 2:
                protein = "Q88X33.txt";
                break;
            case 3:
                protein = "A6L9J9.txt";
                break;
            case 4:
                protein = "Q89BP2.txt";
                break;
            case 5:
                protein = "Q91X51.txt";
                break;
            case 6:
                protein = "B4TWR7.txt";
                break;
            case 7:
                protein = "B3LS90.txt";
                break;
            case 8:
                protein = "A4Y1B6.txt";
                break;
            case 9:
                protein = "Q5VZP5.txt";
                break;
            default:
                break;   
        }
        readAndLoadData(am,protein);                                    //Cargar datos de fichero
        aminoacids = aminoChainToVector(dts.chain);             
        printData();                                                    //Imprimir en pantalla los datos cargados
        printf("Número de experimentos: %d\n",atoi(argv[2]));
        time(&t_start);

        for(int i=0;i<atoi(argv[2]);i++){
            nBorns=0;
            logF[i].id = i+1;
            randomData(am);                                             //Genera CDS aleatorios
            printf("Inicio de la prueba %d\n",i+1);
            genderDistribution2();                                       //Distribuir el género entre los distintos CDS
            calcMultiOjetives(am,aminoacids);                           //Calcular los resultados
            saveToFile("salida " + std::to_string(i+1) + ".txt");       //Guardar datos base en el fichero
            time(&start);
            compute(am,aminoacids,atoi(argv[1]),pm,nBorns);             //Realizar calculos
            time(&end);
            //saveFileParetoFrontier("salidaPF " + std::to_string(i+1) + ".txt");
            logF[i].newBorns = nBorns;
            logF[i].tsecC = (end-start);
            time(&start);
            filter();                                                   //Realizar filtrado final
            time(&end);
            logF[i].tsecF = (end-start);
            logF[i].tPareto = paretoFrontierReduced.size();
            f=0;
            m=paretoFrontierReduced.size();
            for(auto c : paretoFrontierReduced){
                if(c.gender ==1){
                    f++;
                }
            }
            logF[i].females=f;
            logF[i].males =m-f;
            saveToFile("salidaMutada " + std::to_string(i+1) + ".txt"); //Guardar datos mutaciones
            saveFilePareto("salidaPareto " + std::to_string(i+1) + ".txt");//Guarda los individuos en Pareto en fichero
            saveToHVfile("HV_file.txt");                                //Guardar resultados para calcular HV
            strcpy(hv_filename, "HV_file.txt");
            computeHV(hv_filename);                                     //Realizar calculos de hipervolumen
            flushStructures();                                          //Limpiar estructuras de datos para cada iteracion
            printf("FIN de la prueba %d\n",i+1);
        }
        time(&t_end);
        total = (double)(t_end-t_start);
        saveLog(atoi(argv[2]),atoi(argv[1]),total);                     //Guardar informacion de ejecucion en el fichero log 
    }   
    free(am);
    return 0;
}