#ifndef __ROUTE_H__
#define __ROUTE_H__

#include "lib_io.h"
#include "iostream" 
#include "vector"
#include<string>
#include<math.h> 
#include<algorithm>
#include<time.h>
using namespace std;

#define random(a,b) (rand()%((b)-(a)+1)+(a))

void predict_server(char * info[MAX_INFO_NUM], char * data[MAX_DATA_NUM], int data_num, char * filename);

class CleanDate{
public:
  int day;
  int nums;
  int need=0;
};

class INPUT{
public:
  class Server{
  public:
    string name;
    int kernel_num=0;
    int memory=0;
    int harddisk=0;
  };
  class Flaver{
  public:
    string name;
    int flaver_label;
    int memory;
    int kernel;
    double mem_cpu_rate;
  };
  class PreDate{
  public:
    int month;
    int day;
  };
  Server physicalserver[3];
  Flaver flaver[160];
  int flaver_num;
  bool allottype;
  int server_num;
  PreDate prestart;
  PreDate preend;
  
};

class TrainData{
public:
  class Date{
  public:
    int month;
    int day;
  };
  Date nowdate;
  string name;
  int flaver_label=999;
  int flaver_num[160]={0};
  int seven_num[160]={0};
  int delta_num[160]={0};
};

class Flaver{
public:
  int flaver_type;
  double flaver_w[100]={0};
  double flaver_ans[100]={0};
  
};
class PUT_Flaver{
public:
   int kernel=0;
   int memory=0;
   int numbers=0;
   int type;
   string name;
};

///////////////////////////////////////////////////////
 

struct gene{
    public:
        vector<int> seq;
        float score;
        int bagsize[3];
        int leftdiff;
        int lastbag_sum;
        

    gene(vector<int>& seqls,float scorels,int bag0,int bag1,int bag2,int lastbag_sumls,int leftdiffls)
    {
        seq=seqls;
        score=scorels;
        bagsize[0]=bag0;
        bagsize[1]=bag1;
        bagsize[2]=bag2;
        lastbag_sum=lastbag_sumls;
        leftdiff=leftdiffls;

    }
};

struct linknode
{
	int bag_num;
	linknode *next;
	linknode *prior;
    linknode()
    {
    	bag_num=-1;
    	next= nullptr;
    	prior= nullptr;
    }
};




class binpacking{
	public:
        int iyear;


		int totalCPU,totalMEM;
		int serverCPU[3],serverMEM[3];
		int N;  //number kinds of flavors
		int Nall;//number: all flavors,length of gene.seq
        int* cpu;
	    int* mem;
	    int* nums;
        int cpu_mem_ratio[18]; //value=(mem[i]/cpu[i]/2)   0:high performance   1:general  2:high memory
	    int min_cpu;
	    int min_mem;
        int n_serval_type;//how many serval types can be chosen
        int general_fla_cs;
        int belong_server[3];


		vector<int> iid;   //look up table for FF
		
        vector<int> ans_seq; //answer: the gene seq
        float ans_score;
        int ans_bagnum[3];//answer: the bags
        int ans_ff[3][1000][18]={{{0}}};//answer: solution of the binpacking in "high performance"

    float GeneticA_adaptive(int years,int ind_num);

    gene new_first_fit(vector<int>&seq);
    void new_ans_first_fit(); //FF for the answer print
    vector<int> generate_new_seq(vector<int>&seq);

    
   
    binpacking(int*cpuls,int*memls,int*numsls,int *sCPU,int *sMEM,int Nfla)
    {
        N=Nfla;
        cpu=cpuls;
        mem=memls;
        nums=numsls;
        for (int i=0;i<Nfla;i++)
        {
            cpu_mem_ratio[i]=mem[i]/cpu[i];
            cpu_mem_ratio[i]/=2;      //0:high performance   1:general  2:high memory
        }
        totalCPU=0;
        totalMEM=0;



        min_cpu=1000;
        min_mem=1000;
        Nall=0;
        for (int i=0;i<N;i++)
        {
            if (min_cpu>cpu[i]) min_cpu=cpu[i];
            if (min_mem>mem[i]) min_mem=mem[i];
            Nall+=nums[i];
        }
        
        int nums_acc[18];
        for (int i=0;i<N;i++) { nums_acc[i]=nums[i];} 
    	for(int i=1;i<N;i++)
    	{
    		nums_acc[i]+=nums_acc[i-1];
    	}
    	vector<int> iidls(Nall,0);
    	iid=iidls;
    	for(int i=1;i<N;i++)
    		for(int j=nums_acc[i-1];j<nums_acc[i];j++)
    			iid[j]=i;
        
        vector<int> ans_seqls(Nall,0);
        ans_seq=ans_seqls;
        ans_score=0;
        ans_bagnum[0]=-1;//answer: the bags
        ans_bagnum[1]=-1;//answer: the bags
        ans_bagnum[2]=-1;

        for(int k=0;k<N;k++) {
            totalCPU+=cpu[k]*nums[k];
            totalMEM+=mem[k]*nums[k];
        }

        //method for 1-2 server
        belong_server[0]=0;
        belong_server[1]=1;
        belong_server[2]=2;

        float server_cpu_ratio[3];
        server_cpu_ratio[0]=(sCPU[0]==0)? 100:(float(sMEM[0])/sCPU[0]);
        server_cpu_ratio[1]=(sCPU[1]==0)? 100:(float(sMEM[1])/sCPU[1]);
        server_cpu_ratio[2]=(sCPU[2]==0)? 100:(float(sMEM[2])/sCPU[2]);
        int ls;
        float lsf;
        for (int i=0;i<2;i++)
            for (int j=i+1;j<3;j++)
                if (server_cpu_ratio[i]>server_cpu_ratio[j])
                {
                    ls=sCPU[i];sCPU[i]=sCPU[j];sCPU[j]=ls;
                    ls=sMEM[i];sMEM[i]=sMEM[j];sMEM[j]=ls;
                    lsf=server_cpu_ratio[i];server_cpu_ratio[i]=server_cpu_ratio[j];server_cpu_ratio[j]=lsf;
                    ls=belong_server[i];belong_server[i]=belong_server[j];belong_server[j]=ls;

                }

        serverCPU[0]=sCPU[0];
        serverMEM[0]=sMEM[0];
        serverCPU[1]=sCPU[1];
        serverMEM[1]=sMEM[1];
        serverCPU[2]=sCPU[2];
        serverMEM[2]=sMEM[2];
        n_serval_type=3;
        general_fla_cs=-1;
        if (server_cpu_ratio[2]>99)
        {
            if (server_cpu_ratio[1]>99)  //one kind of server: serverCPU[0],serverMEM[0]
            {
                n_serval_type=1;
            }
            else    //two kind of server:serverCPU[0],serverMEM[0];serverCPU[1],serverMEM[1]
            {
                n_serval_type=2;
                if (server_cpu_ratio[0]>1.5)      //0 is general
                {
                    general_fla_cs=0;
                    if (server_cpu_ratio[1]<2.5)  //1 is better general
                        { general_fla_cs=1;}
                }
                else
                    general_fla_cs=1;

            }

        }

    }
};



 
 

#endif