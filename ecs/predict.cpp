#include "predict.h"
#include <stdio.h>
#include <string.h> 
#include<algorithm>
#include<fstream>
#include <math.h> 
#define MAX_ROWS 50000
INPUT S_input;
Flaver cal_w[160],cal_seven_w[160];
int total_days;
int out_sum_flavers;
int out_sum_flavers_temp;
int out_bags;
int serverCPUtemp[3]={0},serverMEMtemp[3]={0};
int out_bags_fusai[3];
string out_bags_fusai_name[3];
int max_label;
TrainData T_data[20000];
TrainData O_data[20000];
int all_flaver[160]={0};
PUT_Flaver p_flaver[160],out_flaver[160];
int s_bag[3][1000][1000];
int s_bag_temp[3][1000][1000];
int final_bag[3][1000][1000];
int train_label_sum[160]={0};
int train_label_avg[160]={0};
int type_avg[100]={0};
int cpu_all_use=0,mem_all_use=0;
int cpu_all_use_temp=0,mem_all_use_temp=0;
CleanDate cleandata[100][200];
int delta_num_avg[160]={0};
int all_bag_delete=0;
int ifdelete=0;
int type_avg_end[100];
int sum_day[1000];
int all_can_CPU=0,all_can_MEM=0;
////////////////////

void s_bag_zero()
{
  for(int i=0;i<3;i++)
    for(int j=0;j<1000;j++)
      for(int k=0;k<1000;k++)
	s_bag[i][j][k]=0;
}
bool cmp_population(const gene &x,const gene &y)
{
    if (x.score==y.score)
    {
        if (x.lastbag_sum==y.lastbag_sum)
        {
            return x.leftdiff<y.leftdiff;
        }
        return x.lastbag_sum>y.lastbag_sum;
    }
    return x.score>y.score;

}

float binpacking::GeneticA_adaptive(int years,int ind_num)
{
    vector<int> seq(Nall,0);
    vector<gene> population;
    clock_t start_time,end_time;
    //start_time = clock();

    for (int i=0;i<Nall;i++) { seq[i]=i;}

    for (int i=0;i<(3*ind_num);i++)
    {
        std::random_shuffle(seq.begin(),seq.end());
        population.push_back(new_first_fit(seq));
    }

    sort(population.begin(),population.end(),cmp_population);
    cout<<"init-finished"<<endl;


    //adaptive stage1： previous 1/3 years, one father only generate two child
    int father_num=ind_num;
    int years_stage1=years/3;
    for (iyear=0;iyear<years_stage1;iyear++)
    {
        for (int i=0;i<father_num;i++)
        {

            seq=generate_new_seq(population[i].seq);
            population[i+father_num]=new_first_fit(seq);

            seq=generate_new_seq(population[i].seq);
            population[i+father_num*2]=new_first_fit(seq);

        }

        sort(population.begin(),population.end(),cmp_population);
        if (iyear % 200==0)
        {
            cout << "***" << iyear << "***" << population[0].score*50 << "  " << population[0].bagsize[0] << "  "
                 << population[0].bagsize[1] << "  " << population[0].bagsize[2]
                 << "  "<<population[0].lastbag_sum<<"  "<<population[0].leftdiff<<endl;

            //end_time = clock();
            //cout << "Running Time : " << (double) (end_time - start_time) / CLOCKS_PER_SEC << endl;
        }


    }

    //adaptive stage2： latter 2/3 years, one father generate more than two child
    father_num=ind_num/5;
    vector<gene>::iterator end_stage2=population.begin();
    for (int i=0;i<father_num*10;i++) {end_stage2++;}

    for (iyear=years_stage1;iyear<years;iyear++)
    {


        for (int i=0;i<father_num;i++)
        {

            seq=generate_new_seq(population[i].seq);
            population[i+father_num]=new_first_fit(seq);

            seq=generate_new_seq(population[i].seq);
            population[i+father_num*2]=new_first_fit(seq);

            seq=generate_new_seq(population[i].seq);
            population[i+father_num*3]=new_first_fit(seq);

            seq=generate_new_seq(population[i].seq);
            population[i+father_num*4]=new_first_fit(seq);

            seq=generate_new_seq(population[i].seq);
            population[i+father_num*5]=new_first_fit(seq);

            seq=generate_new_seq(population[i].seq);
            population[i+father_num*6]=new_first_fit(seq);

            seq=generate_new_seq(population[i].seq);
            population[i+father_num*7]=new_first_fit(seq);

            seq=generate_new_seq(population[i].seq);
            population[i+father_num*8]=new_first_fit(seq);

            seq=generate_new_seq(population[i].seq);
            population[i+father_num*9]=new_first_fit(seq);
        }

        sort(population.begin(),end_stage2,cmp_population);
        if (iyear % 200==0)
            cout << "***" << iyear << "***" << population[0].score*50 << "  " << population[0].bagsize[0] << "  "
                 << population[0].bagsize[1] << "  " << population[0].bagsize[2]
                 << "  "<<population[0].lastbag_sum<<"  "<<population[0].leftdiff<<endl;



    }


    //cout<<"-----"<<years<<"years finished"<<endl;
    ans_seq=population[0].seq;
    ans_score=population[0].score*50;
    ans_bagnum[belong_server[0]]=population[0].bagsize[0];
    ans_bagnum[belong_server[1]]=population[0].bagsize[1];
    ans_bagnum[belong_server[2]]=population[0].bagsize[2];
    new_ans_first_fit();


    return ans_score;
}

gene binpacking::new_first_fit(vector<int>&seq)
{
    int i,j,k,cs;
    vector<int>left_cpu0(800,serverCPU[0]);
    vector<int>left_mem0(800,serverMEM[0]);
    vector<int>left_cpu1(800,serverCPU[1]);
    vector<int>left_mem1(800,serverMEM[1]);
    vector<int>left_cpu2(800,serverCPU[2]);
    vector<int>left_mem2(800,serverMEM[2]);

    int bag_size0=0,bag_size1=0,bag_size2=0;//number of bags(three kinds sever)
    linknode bag_node0[800],bag_node1[800],bag_node2[800];
    linknode *bag_list0=nullptr,*bag_list1=nullptr,*bag_list2=nullptr;
    linknode *tail0,*tail1,*tail2;
    linknode *ptr0,*ptr1,*ptr2;
    linknode *temp;

    float tmp_left0,tmp_left1,tmp_left2;

    bool flag;
    bool flag_serval[3];//weather this kind of  serval can be put into this fla
    unsigned long seqsize=seq.size();
    for(int s=0;s<seqsize;s++)
    {
        i=iid[seq[s]];//seq_element mapping to kind of flavor
        flag_serval[0]=false;
        ptr0=bag_list0;
        while(ptr0!=nullptr)
        {
            j=ptr0->bag_num;
            if(left_cpu0[j]>=cpu[i]&&left_mem0[j]>=mem[i])
            {
                flag_serval[0]=true;
                break;
            }
            ptr0=ptr0->next;
        }

        flag_serval[1]=false;
        ptr1=bag_list1;
        while(ptr1!=nullptr)
        {
            j=ptr1->bag_num;
            if(left_cpu1[j]>=cpu[i]&&left_mem1[j]>=mem[i])
            {
                flag_serval[1]=true;
                break;
            }
            ptr1=ptr1->next;
        }

        flag_serval[2]=false;
        ptr2=bag_list2;
        while(ptr2!=nullptr)
        {
            j=ptr2->bag_num;
            if(left_cpu2[j]>=cpu[i]&&left_mem2[j]>=mem[i])
            {
                flag_serval[2]=true;
                break;
            }
            ptr2=ptr2->next;
        }
        flag=flag_serval[0]||flag_serval[1]||flag_serval[2];

        //choose cs to put

        if(flag)
        {
            if(n_serval_type==3) {
                if (flag_serval[cpu_mem_ratio[i]]) cs = cpu_mem_ratio[i];
                else
                {
                    if(cpu_mem_ratio[i]!=1)//fix here to decide if to use method 1
                    {
                        if (flag_serval[1]) cs = 1;
                        else if (flag_serval[2]) cs = 2;
                        else cs = 0;
                    }
                    else//method1 choose the best place to put 1:2 type flavor
                    {
                        if(flag_serval[0]&&flag_serval[2])
                        {
                            j=ptr0->bag_num;
                            tmp_left0=(left_mem0[j]-mem[i]+serverMEM[0])/float(left_cpu0[j]-cpu[i]+serverCPU[0]);
                            j=ptr2->bag_num;
                            tmp_left2=(left_mem2[j]-mem[i]+serverMEM[2])/float(left_cpu2[j]-cpu[i]+serverCPU[2]);
                            tmp_left0=tmp_left0*serverCPU[0]/serverMEM[0];
                            tmp_left2=tmp_left2*serverCPU[2]/serverMEM[2];
                            if (tmp_left0>1) tmp_left0=1/tmp_left0;
                            if (tmp_left2>1) tmp_left2=1/tmp_left2;
                            if(tmp_left0>tmp_left2) cs=0;
                            else cs=2;

                        }
                        else
                        {
                            if(flag_serval[0]) cs=0;
                            else cs=2;
                        }
                    }
                }

            }
            else if(n_serval_type==1)
                cs=0;
            else
            {
                if(cpu_mem_ratio[i]==0)
                {
                    if(flag_serval[0])
                        cs=0;
                    else
                        cs=1;
                }
                else if (cpu_mem_ratio[i]==2)
                {
                    if(flag_serval[1])
                        cs=1;
                    else
                        cs=0;
                }
                else
                {
                    if(flag_serval[general_fla_cs])
                        cs=general_fla_cs;
                    else
                        cs=1-general_fla_cs;
                }
            }



            if(cs==0)
            {
                j=ptr0->bag_num;
                left_cpu0[j]-=cpu[i];
                left_mem0[j]-=mem[i];
                if(left_cpu0[j]<min_cpu||left_mem0[j]<min_mem)
                {
                    if(bag_list0==ptr0)
                        bag_list0=ptr0->next;
                    else if(tail0==ptr0)
                    {
                        ptr0->prior->next=nullptr;
                        tail0=ptr0->prior;
                    }
                    else
                    {
                        ptr0->prior->next=ptr0->next;
                        ptr0->next->prior=ptr0->prior;
                    }
                }
            }
            else if(cs==1)
            {
                j=ptr1->bag_num;
                left_cpu1[j]-=cpu[i];
                left_mem1[j]-=mem[i];
                if(left_cpu1[j]<min_cpu||left_mem1[j]<min_mem)
                {
                    if(bag_list1==ptr1)
                        bag_list1=ptr1->next;
                    else if(tail1==ptr1)
                    {
                        ptr1->prior->next=nullptr;
                        tail1=ptr1->prior;
                    }
                    else
                    {
                        ptr1->prior->next=ptr1->next;
                        ptr1->next->prior=ptr1->prior;
                    }
                }
            }
            else
            {
                j=ptr2->bag_num;
                left_cpu2[j]-=cpu[i];
                left_mem2[j]-=mem[i];
                if(left_cpu2[j]<min_cpu||left_mem2[j]<min_mem)
                {
                    if(bag_list2==ptr2)
                        bag_list2=ptr2->next;
                    else if(tail2==ptr2)
                    {
                        ptr2->prior->next=nullptr;
                        tail2=ptr2->prior;
                    }
                    else
                    {
                        ptr2->prior->next=ptr2->next;
                        ptr2->next->prior=ptr2->prior;
                    }
                }
            }
        }
        else
        {
            if(n_serval_type==3)
                cs=cpu_mem_ratio[i];
            else if(n_serval_type==1)
                cs=0;
            else
            {
                if(cpu_mem_ratio[i]==2)
                    cs=1;
                else if(cpu_mem_ratio[i]==0)
                    cs=0;
                else
                {
                    cs=general_fla_cs;
                }

            }


            if(cs==0)
            {
                if(bag_list0==nullptr)
                {
                    bag_list0=&bag_node0[bag_size0];
                    tail0=bag_list0;
                    tail0->bag_num=bag_size0;
                }
                else
                {
                    tail0->next=&bag_node0[bag_size0];
                    temp=tail0;
                    tail0=tail0->next;
                    tail0->bag_num=bag_size0;
                    tail0->prior=temp;
                }

                left_cpu0[bag_size0]-=cpu[i];
                left_mem0[bag_size0]-=mem[i];
                bag_size0+=1;
            }
            else if(cs==1)
            {
                if(bag_list1==nullptr)
                {
                    bag_list1=&bag_node1[bag_size1];
                    tail1=bag_list1;
                    tail1->bag_num=bag_size1;
                }
                else
                {
                    tail1->next=&bag_node1[bag_size1];
                    temp=tail1;
                    tail1=tail1->next;
                    tail1->bag_num=bag_size1;
                    tail1->prior=temp;
                }

                left_cpu1[bag_size1]-=cpu[i];
                left_mem1[bag_size1]-=mem[i];
                bag_size1+=1;
            }
            else
            {
                if(bag_list2==nullptr)
                {
                    bag_list2=&bag_node2[bag_size2];
                    tail2=bag_list2;
                    tail2->bag_num=bag_size2;
                }
                else
                {
                    tail2->next=&bag_node2[bag_size2];
                    temp=tail2;
                    tail2=tail2->next;
                    tail2->bag_num=bag_size2;
                    tail2->prior=temp;
                }

                left_cpu2[bag_size2]-=cpu[i];
                left_mem2[bag_size2]-=mem[i];
                bag_size2+=1;
            }
        }
    }


    int total_cpu_ser=bag_size0*serverCPU[0]+bag_size1*serverCPU[1]+bag_size2*serverCPU[2];
    int total_mem_ser=bag_size0*serverMEM[0]+bag_size1*serverMEM[1]+bag_size2*serverMEM[2];
    float score=(totalCPU/float(total_cpu_ser)+totalMEM/float(total_mem_ser));

    //left_diff
    int leftdiff_sum=0;
    int ls1=0;
    for(i=0;i<bag_size0;i++)
    {
        ls1=left_cpu0[i]-left_mem0[i];
        ls1*=ls1;
        leftdiff_sum+=ls1;
    }
    for(i=0;i<bag_size1;i++)
    {
        ls1=left_cpu1[i]-left_mem1[i];
        ls1*=ls1;
        leftdiff_sum+=ls1;
    }
    for(i=0;i<bag_size2;i++)
    {
        ls1=left_cpu2[i]-left_mem2[i];
        ls1*=ls1;
        leftdiff_sum+=ls1;
    }


    //lastbag_sum
    int lastbag_sum=0;

    j=-1;
    ptr0=bag_list0;
    while(ptr0!=nullptr)  { j=ptr0->bag_num;ptr0=ptr0->next; }
    if (j!=-1)
    {
        lastbag_sum+=left_cpu0[j]*left_mem0[j];
        leftdiff_sum-=(left_cpu0[j]-left_mem0[j])^2;  // left_diff should be subtract the last_bag
    }
    j=-1;
    ptr1=bag_list1;
    while(ptr1!=nullptr)  { j=ptr1->bag_num;ptr1=ptr1->next; }
    if (j!=-1)
    {
        lastbag_sum+=left_cpu1[j]*left_mem1[j];
        leftdiff_sum-=(left_cpu1[j]-left_mem1[j])^2;  // left_diff should be subtract the last_bag
    }
    j=-1;
    ptr2=bag_list2;
    while(ptr2!=nullptr)  { j=ptr2->bag_num;ptr2=ptr2->next; }
    if (j!=-1)
    {
        lastbag_sum+=left_cpu2[j]*left_mem2[j];
        leftdiff_sum-=(left_cpu2[j]-left_mem2[j])^2;  // left_diff should be subtract the last_bag
    }
    gene newFFgene(seq,score,bag_size0,bag_size1,bag_size2,lastbag_sum,leftdiff_sum);
    return newFFgene;
}



vector<int> binpacking::generate_new_seq(vector<int>& seq)
{
   
    
    vector<int> newseq(seq);
    int i;
    int index1,index2;
    
    //exchange several point
    
    int exchange_bit;
    //srand((unsigned)time(NULL));
    for (i=0;i<2;i++)
    {
        index1=random(0,Nall-1);
        index2=random(0,Nall-1);
        exchange_bit=newseq[index1];
        newseq[index1]=newseq[index2];
        newseq[index2]=exchange_bit;
    }
    

   
  //print for check
/* 
    cout<<"ex_len="<<exchange_length<<"  "<<"index1="<<index1<<"  index2="<<index2<<endl;
    cout<<"orgin_seq="<<endl;cout<<"  ";
    for (int j=0;j<Nall;j++) {cout<<seq[j]<<' ';}
    cout<<endl;
    cout<<"new_seq="<<endl;cout<<"  ";
    for (int j=0;j<Nall;j++) {cout<<newseq[j]<<' ';}
    cout<<endl; 
    cout<<endl;
*/

  //newseq_chek
  /*
    vector<int> checkseq(Nall,0);
    for (i=0;i<Nall;i++){checkseq[newseq[i]]++;}
    for (i=0;i<Nall;i++)
    {
        if (checkseq[i]!=1) 
        {
            cout<<"This Newseq is  illegal! "<<endl;
            cout<<"len="<<reverse_length<<"  "<<"index1="<<index1<<"  index2="<<index2<<endl;
            cout<<"orgin_seq="<<endl;cout<<"  ";
            for (int j=0;j<Nall;j++) {cout<<seq[j]<<' ';}
            cout<<endl;
            cout<<"new_seq="<<endl;cout<<"  ";
            for (int j=0;j<Nall;j++) {cout<<newseq[j]<<' ';}
            cout<<endl; 
        }
    }
    */
    return newseq;
        
}


void binpacking::new_ans_first_fit()
{
    int i,j,k,cs;
    vector<int> seq;
    seq=ans_seq;
    vector<int>left_cpu0(800,serverCPU[0]);
    vector<int>left_mem0(800,serverMEM[0]);
    vector<int>left_cpu1(800,serverCPU[1]);
    vector<int>left_mem1(800,serverMEM[1]);
    vector<int>left_cpu2(800,serverCPU[2]);
    vector<int>left_mem2(800,serverMEM[2]);

    int bag_size0=0,bag_size1=0,bag_size2=0;//number of bags(three kinds sever)
    linknode bag_node0[800],bag_node1[800],bag_node2[800];
    linknode *bag_list0=nullptr,*bag_list1=nullptr,*bag_list2=nullptr;
    linknode *tail0,*tail1,*tail2;
    linknode *ptr0,*ptr1,*ptr2;
    linknode *temp;
    for(k=0;k<3;k++)
        for (i=0;i<ans_bagnum[k];i++)
            for (j=0;j<N;j++)
                ans_ff[k][i][j]=0;

    float tmp_left0,tmp_left2;

    bool flag;
    bool flag_serval[3];//weather this kind of  serval can be put into this fla
    unsigned long seqsize=seq.size();
    for(int s=0;s<seqsize;s++)
    {
        i=iid[seq[s]];//seq_element mapping to kind of flavor
        flag_serval[0]=false;
        ptr0=bag_list0;
        while(ptr0!=nullptr)
        {
            j=ptr0->bag_num;
            if(left_cpu0[j]>=cpu[i]&&left_mem0[j]>=mem[i])
            {
                flag_serval[0]=true;
                break;
            }
            ptr0=ptr0->next;
        }

        flag_serval[1]=false;
        ptr1=bag_list1;
        while(ptr1!=nullptr)
        {
            j=ptr1->bag_num;
            if(left_cpu1[j]>=cpu[i]&&left_mem1[j]>=mem[i])
            {
                flag_serval[1]=true;
                break;
            }
            ptr1=ptr1->next;
        }

        flag_serval[2]=false;
        ptr2=bag_list2;
        while(ptr2!=nullptr)
        {
            j=ptr2->bag_num;
            if(left_cpu2[j]>=cpu[i]&&left_mem2[j]>=mem[i])
            {
                flag_serval[2]=true;
                break;
            }
            ptr2=ptr2->next;
        }
        flag=flag_serval[0]||flag_serval[1]||flag_serval[2];

        //choose cs to put

        if(flag)
        {
            if(n_serval_type==3) {
                if (flag_serval[cpu_mem_ratio[i]]) cs = cpu_mem_ratio[i];
                else
                {
                    if(cpu_mem_ratio[i]!=1)//fix here to decide if to use method 1
                    {
                        if (flag_serval[1]) cs = 1;
                        else if (flag_serval[2]) cs = 2;
                        else cs = 0;
                    }
                    else//method1 choose the best place to put 1:2 type flavor
                    {
                        if(flag_serval[0]&&flag_serval[2])
                        {
                            j=ptr0->bag_num;
                            tmp_left0=(left_mem0[j]-mem[i]+serverMEM[0])/float(left_cpu0[j]-cpu[i]+serverCPU[0]);
                            j=ptr2->bag_num;
                            tmp_left2=(left_mem2[j]-mem[i]+serverMEM[2])/float(left_cpu2[j]-cpu[i]+serverCPU[2]);
                            tmp_left0=tmp_left0*serverCPU[0]/serverMEM[0];
                            tmp_left2=tmp_left2*serverCPU[2]/serverMEM[2];
                            if (tmp_left0>1) tmp_left0=1/tmp_left0;
                            if (tmp_left2>1) tmp_left2=1/tmp_left2;
                            if(tmp_left0>tmp_left2) cs=0;
                            else cs=2;

                        }
                        else
                        {
                            if(flag_serval[0]) cs=0;
                            else cs=2;
                        }
                    }
                }

            }
            else if(n_serval_type==1)
                cs=0;
            else
            {
                if(cpu_mem_ratio[i]==0)
                {
                    if(flag_serval[0])
                        cs=0;
                    else
                        cs=1;
                }
                else if (cpu_mem_ratio[i]==2)
                {
                    if(flag_serval[1])
                        cs=1;
                    else
                        cs=0;
                }
                else
                {
                    if(flag_serval[general_fla_cs])
                        cs=general_fla_cs;
                    else
                        cs=1-general_fla_cs;
                }
            }



            if(cs==0)
            {
                j=ptr0->bag_num;
                left_cpu0[j]-=cpu[i];
                left_mem0[j]-=mem[i];
                ans_ff[belong_server[cs]][j][i]+=1;
                if(left_cpu0[j]<min_cpu||left_mem0[j]<min_mem)
                {
                    if(bag_list0==ptr0)
                        bag_list0=ptr0->next;
                    else if(tail0==ptr0)
                    {
                        ptr0->prior->next=nullptr;
                        tail0=ptr0->prior;
                    }
                    else
                    {
                        ptr0->prior->next=ptr0->next;
                        ptr0->next->prior=ptr0->prior;
                    }
                }
            }
            else if(cs==1)
            {
                j=ptr1->bag_num;
                left_cpu1[j]-=cpu[i];
                left_mem1[j]-=mem[i];
                ans_ff[belong_server[cs]][j][i]+=1;
                if(left_cpu1[j]<min_cpu||left_mem1[j]<min_mem)
                {
                    if(bag_list1==ptr1)
                        bag_list1=ptr1->next;
                    else if(tail1==ptr1)
                    {
                        ptr1->prior->next=nullptr;
                        tail1=ptr1->prior;
                    }
                    else
                    {
                        ptr1->prior->next=ptr1->next;
                        ptr1->next->prior=ptr1->prior;
                    }
                }
            }
            else
            {
                j=ptr2->bag_num;
                left_cpu2[j]-=cpu[i];
                left_mem2[j]-=mem[i];
                ans_ff[belong_server[cs]][j][i]+=1;
                if(left_cpu2[j]<min_cpu||left_mem2[j]<min_mem)
                {
                    if(bag_list2==ptr2)
                        bag_list2=ptr2->next;
                    else if(tail2==ptr2)
                    {
                        ptr2->prior->next=nullptr;
                        tail2=ptr2->prior;
                    }
                    else
                    {
                        ptr2->prior->next=ptr2->next;
                        ptr2->next->prior=ptr2->prior;
                    }
                }
            }
        }
        else
        {
            if(n_serval_type==3)
                cs=cpu_mem_ratio[i];
            else if(n_serval_type==1)
                cs=0;
            else
            {
                if(cpu_mem_ratio[i]==2)
                    cs=1;
                else if(cpu_mem_ratio[i]==0)
                    cs=0;
                else
                {
                    cs=general_fla_cs;
                }

            }


            if(cs==0)
            {
                if(bag_list0==nullptr)
                {
                    bag_list0=&bag_node0[bag_size0];
                    tail0=bag_list0;
                    tail0->bag_num=bag_size0;
                }
                else
                {
                    tail0->next=&bag_node0[bag_size0];
                    temp=tail0;
                    tail0=tail0->next;
                    tail0->bag_num=bag_size0;
                    tail0->prior=temp;
                }

                left_cpu0[bag_size0]-=cpu[i];
                left_mem0[bag_size0]-=mem[i];
                ans_ff[belong_server[cs]][bag_size0][i]+=1;
                bag_size0+=1;
            }
            else if(cs==1)
            {
                if(bag_list1==nullptr)
                {
                    bag_list1=&bag_node1[bag_size1];
                    tail1=bag_list1;
                    tail1->bag_num=bag_size1;
                }
                else
                {
                    tail1->next=&bag_node1[bag_size1];
                    temp=tail1;
                    tail1=tail1->next;
                    tail1->bag_num=bag_size1;
                    tail1->prior=temp;
                }

                left_cpu1[bag_size1]-=cpu[i];
                left_mem1[bag_size1]-=mem[i];
                ans_ff[belong_server[cs]][bag_size1][i]+=1;
                bag_size1+=1;
            }
            else
            {
                if(bag_list2==nullptr)
                {
                    bag_list2=&bag_node2[bag_size2];
                    tail2=bag_list2;
                    tail2->bag_num=bag_size2;
                }
                else
                {
                    tail2->next=&bag_node2[bag_size2];
                    temp=tail2;
                    tail2=tail2->next;
                    tail2->bag_num=bag_size2;
                    tail2->prior=temp;
                }

                left_cpu2[bag_size2]-=cpu[i];
                left_mem2[bag_size2]-=mem[i];
                ans_ff[belong_server[cs]][bag_size2][i]+=1;
                bag_size2+=1;
            }
        }
    }
}




//////////////////////



void read_info(char * info[MAX_INFO_NUM])
{
  const char * split = " "; 
  const char * split_1 = "-"; 
  S_input.server_num=atoi(info[0]);
  char  * server;
  for(int i=0;i<S_input.server_num;i++)
  {
    server = strtok(info[i+1],split);
    S_input.physicalserver[i].name=server;
    cout<<"name: "<<server<<endl;
    server = strtok(NULL,split);  
    S_input.physicalserver[i].kernel_num = stoi(server);
    cout<<" kernel: "<<stoi(server)<<endl;
    server = strtok(NULL,split);  
    S_input.physicalserver[i].memory = stoi(server);
    cout<<" mem: "<<atoi(server)<<endl;
    server = strtok(NULL,split); 
    S_input.physicalserver[i].harddisk=stoi(server);
  }
  S_input.flaver_num = atoi(info[2+S_input.server_num]);
  for(int i=0;i<S_input.flaver_num;i++)
  {
    char * flaver;
    flaver = strtok(info[3+i+S_input.server_num],split);   
    S_input.flaver[i].name =flaver;
    string str=flaver;
    if(str.length()==7){
      str=str.substr(str.length()-1,1);
      S_input.flaver[i].flaver_label=stoi(str);
    }
    else if(str.length()==9){
      str=str.substr(str.length()-3,3);
      S_input.flaver[i].flaver_label=stoi(str);
    }
    else if(str.length()==10){
      str=str.substr(str.length()-4,4);
      S_input.flaver[i].flaver_label=stoi(str);
    }
    else{
      str=str.substr(str.length()-2,2);
      S_input.flaver[i].flaver_label=stoi(str); 
    }

    flaver= strtok(NULL,split); 
    S_input.flaver[i].kernel=atoi(flaver);
    flaver=strtok(NULL,split); 
    S_input.flaver[i].memory=atoi(flaver)/1024;
    S_input.flaver[i].mem_cpu_rate = S_input.flaver[i].memory/S_input.flaver[i].kernel;
    cout<< S_input.flaver[i].name<<" cpu:"<<S_input.flaver[i].kernel<<" mem:"<<S_input.flaver[i].memory <<" rate: "<<S_input.flaver[i].mem_cpu_rate<<endl;
  }
/*cout<<"[1]"<<int((info[S_input.flaver_num+4][1]))<<endl;
  cout<<"[3]"<<int(info[S_input.flaver_num+4][3])<<endl;
  cout<<"[4]"<<int(info[S_input.flaver_num+4][4])<<endl;
  cout<<strlen(cmp_1)<<":"<<strlen(cmp_2)<<endl;*/
 
  char * starttime;
  starttime=strtok(info[S_input.flaver_num+4+S_input.server_num],split);
  cout<<starttime<<endl;
  char * startdate;
  startdate=strtok(starttime,split_1);
  startdate=strtok(NULL,split_1);
  S_input.prestart.month=atoi(startdate);
  startdate=strtok(NULL,split_1);
  S_input.prestart.day=atoi(startdate);
//  cout<<S_input.prestart.month<<":"<<S_input.prestart.day<<endl;
  
  char * endtime;
  endtime=strtok(info[S_input.flaver_num+5+S_input.server_num],split); 
  char * enddate;
  enddate=strtok(endtime,split_1);
  enddate=strtok(NULL,split_1);
  S_input.preend.month=atoi(enddate);
  enddate=strtok(NULL,split_1);
  S_input.preend.day=atoi(enddate);
//  cout<<S_input.preend.month<<":"<<S_input.preend.day<<endl;
}

void read_data(char * data[MAX_DATA_NUM],int linenum)
{
  const char * split = "\t"; 
  const char * split_1 = "-";
  const char * split_2 = " "; 
  int j=0,temp=0;
  cout<<"line::"<<linenum<<endl;
  for(int i=0;i<linenum;i++){
    char * flaver;
   // cout<<data[i];
    flaver = strtok(data[i],split);
    flaver = strtok(NULL,split);
    string str=flaver;
    flaver = strtok(NULL,split);
    char * f_date;
    char * f_time;
    f_time=strtok(flaver,split_2);
    f_date=strtok(f_time,split_1);
    f_date=strtok(NULL,split_1);
    O_data[i].nowdate.month=atoi(f_date);
    f_date=strtok(NULL,split_1);
    O_data[i].nowdate.day=atoi(f_date);
    if((O_data[i].nowdate.day==O_data[i-1].nowdate.day)){
      j=j;
    }
    else{
      j++;
    }
    
    T_data[j].name = str;
    if(str.length()==7){
      str=str.substr(str.length()-1,1);
      T_data[j].flaver_label=stoi(str);
    }
    else{
      str=str.substr(str.length()-2,2);
      T_data[j].flaver_label=stoi(str); 
    }
    
    int label=T_data[j].flaver_label;
    if(label>temp)
    {
      temp=label;
    }
    T_data[j].flaver_num[label]+=1;
    
    T_data[j].flaver_num[0]=0;
    T_data[j].nowdate.month=O_data[i].nowdate.month;    
    T_data[j].nowdate.day=O_data[i].nowdate.day;
 //   cout<<i<<"  "<<j<<"--"<<T_data[j].nowdate.month<<":"<<T_data[j].nowdate.day<<" label:"<<T_data[j].flaver_label<<" num:"<<T_data[j].flaver_num[label]<<endl;   
   }
   total_days=j;
   max_label=temp;
   cout<<"total_days : "<<total_days<<endl;
   cout<<"max label : "<<max_label<<endl;
   cout<<"end-day: "<<T_data[total_days].nowdate.month<<"-"<<T_data[total_days].nowdate.day<<endl;
  
}

int train_test_time()
{
  int month_days[13]={0,31,28,31,30,31,30,31,31,30,31,30,31};
  if(S_input.prestart.month==T_data[total_days].nowdate.month)
  {
    return  S_input.prestart.day-T_data[total_days].nowdate.day;
  }
  else{
    return S_input.prestart.day+month_days[T_data[total_days].nowdate.month]-T_data[total_days].nowdate.day;
  }
}

int pre_times()
{
  int month_days[13]={0,31,28,31,30,31,30,31,31,30,31,30,31};
  if(S_input.preend.month==S_input.prestart.month)
  {
    return S_input.preend.day-S_input.prestart.day;
  }
  else{
    return S_input.preend.day+month_days[S_input.prestart.month]-S_input.prestart.day;
  }
  
}
void view_nocal_train()
{
   for(int i=1;;i++){
     if(T_data[i].flaver_label==999)
       return;
     cout<<i<<"::";
     for(int j=0;j<max_label+1;j++)
     {
       cout<<T_data[i].flaver_num[j]<<" ";
     }
     cout<<endl;
   }
}
void view_train()
{
   for(int i=1;;i++){
     if(T_data[i].flaver_label==999)
       return;
     cout<<i<<"::";
     for(int j=0;j<max_label+1;j++)
     {
       train_label_sum[j]+=T_data[i].flaver_num[j];
       cout<<T_data[i].flaver_num[j]<<" ";
     }
     cout<<endl;
   }
}
void get_delta_train()
{
  int delta_times = pre_times();
  int start=total_days-delta_times+1-7;
  for(int i=1;i<=total_days-delta_times+2-start;i++)
  {
    cout<<i<<":delta:";
    for(int j=0;j<=max_label;j++)
    {
      for(int k=0;k<delta_times;k++)
      {
	T_data[i].delta_num[j]+=T_data[i+k-1+start].flaver_num[j];
      }
      cout<<T_data[i].delta_num[j]<<" ";
      delta_num_avg[j]+=T_data[i].delta_num[j];
    }
    cout<<endl;    
  }
  for(int i=0;i<max_label;i++)
  {
    delta_num_avg[i]=delta_num_avg[i]/(total_days-delta_times+2-start);
    cout<<" "<<i<<" avg: "<<delta_num_avg[i]<<" "<<endl;
  }
  cout<<endl;
}






void out_csv_file()
{
   ofstream write;
   write.open("/home/hou/Windy/DEVcloud/data/train1512.csv");
   write<<"date";
    for(int j=1;j<max_label+1;j++)
    {
      write<<","<<"flavor"<<j;
    }
    write<<'\n';
    for(int i=1;;i++){
     if(T_data[i].flaver_label==999)
       return;
     write<<T_data[i].nowdate.month<<'-'<<T_data[i].nowdate.day;
     for(int j=1;j<max_label+1;j++)
     {
       write<<","<<T_data[i].flaver_num[j];
     }
     write<<'\n';
   }
   write.close();
}

void avg_cal()
{
  int delta_time = pre_times(); 
  for(int i=1;i<max_label+1;i++){
     train_label_avg[i]=train_label_sum[i]*delta_time/total_days;
     cout<<train_label_avg[i]<<" = ";
  }
  cout<<endl;
}


void FreeData(double **dat, double *d, int count)
{
    int i;
    free(d);
    for (i = 0; i < count; i ++)
        free(dat[i]);
    free(dat);
} 
int LinearEquations(double *data, int count, double *Answer)
{
    int j, m, n;
    double tmp, **dat, *d = data;
    dat = (double**)malloc(count * sizeof(double*));
    for (m = 0; m < count; m ++, d += (count + 1))
    {
        dat[m] = (double*)malloc((count + 1) * sizeof(double));
        memcpy(dat[m], d, (count + 1) * sizeof(double));
    }
    d = (double*)malloc((count + 1) * sizeof(double));
    for (m = 0; m < count - 1; m ++)
    { 
        for (n = m + 1; n < count && dat[m][m] == 0.0; n ++)
        {
            if ( dat[n][m] != 0.0)
            {
                memcpy(d, dat[m], (count + 1) * sizeof(double));
                memcpy(dat[m], dat[n], (count + 1) * sizeof(double));
                memcpy(dat[n], d, (count + 1) * sizeof(double));
            }
        }
        if (dat[m][m] == 0.0)
        {
            FreeData(dat, d, count);
            return -1;
        } 
        for (n = m + 1; n < count; n ++)
        {
            tmp = dat[n][m] / dat[m][m];
            for (j = m; j <= count; j ++)
                dat[n][j] -= tmp * dat[m][j];
        }
    }
    for (j = 0; j < count; j ++)
        d[j] = 0.0; 
    Answer[count - 1] = dat[count - 1][count] / dat[count - 1][count - 1]; 
    for (m = count - 2; m >= 0; m --)
    {
        for (j = count - 1; j > m; j --)
            d[m] += Answer[j] * dat[m][j];
        Answer[m] = (dat[m][count] - d[m]) / dat[m][m];
    }
    FreeData(dat, d, count);
    return 0;
}

 
int MultipleRegression(double *data, int rows, int cols, double *Answer, double *SquarePoor)
{
    int m, n, i, count = cols - 1;
    double *dat, *p, a, b;
    if (data == 0 || Answer == 0 || rows < 2 || cols < 2)
        return -1;
    dat = (double*)malloc(cols * (cols + 1) * sizeof(double));
    dat[0] = (double)rows;
    for (n = 0; n < count; n ++)                     // n = 0 to cols - 2
    {
        a = b = 0.0;
        for (p = data + n, m = 0; m < rows; m ++, p += cols)
        {
            a += *p;
            b += (*p * *p);
        }
        dat[n + 1] = a;                            
        dat[(n + 1) * (cols + 1)] = a;              
        dat[(n + 1) * (cols + 1) + n + 1] = b;      
        for (i = n + 1; i < count; i ++)             
        {
            for (a = 0.0, p = data, m = 0; m < rows; m ++, p += cols)
                a += (p[n] * p[i]);
            dat[(n + 1) * (cols + 1) + i + 1] = a;   
            dat[(i + 1) * (cols + 1) + n + 1] = a;  
        }
    }
    for (b = 0.0, m = 0, p = data + n; m < rows; m ++, p += cols)
        b += *p;
    dat[cols] = b;                                  
    for (n = 0; n < count; n ++)
    {
        for (a = 0.0, p = data, m = 0; m < rows; m ++, p += cols)
            a += (p[n] * p[count]);
        dat[(n + 1) * (cols + 1) + cols] = a;       
    }
    n = LinearEquations(dat, cols, Answer);          
    
    if (n == 0 && SquarePoor)
    {
        b = b / rows;                               
        SquarePoor[0] = SquarePoor[1] = 0.0;
        p = data;
        for (m = 0; m < rows; m ++, p ++)
        {
            for (i = 1, a = Answer[0]; i < cols; i ++, p ++)
                a += (*p * Answer[i]);              
            SquarePoor[0] += ((a - b) * (a - b));    
            SquarePoor[1] += ((*p - a) * (*p - a)); 
        }
        SquarePoor[2] = SquarePoor[0] / count;      
	if (rows - cols > 0.0)
	  SquarePoor[3] = SquarePoor[1] / (rows - cols); 
	else
	  SquarePoor[3] = 0.0;
    }
    free(dat);
    return n;
}


/*
void sum_day_flaver()
{
  int temp_sum=0,temp_avg_a=0;
  for(int i=0;i<=max_label;i++)
  {
    for(int j=total_days-14;j<=total_days;j++)
    {
      sum_day[j] += T_data[j].flaver_num[i];
      temp_sum+=T_data[j].flaver_num[i];
    } 
  }
    for(int j=0;j<=total_days;j++)
    {
      cout<<"sum_day   "<<sum_day[j]<<endl;
    } 
  temp_avg_a=temp_sum/14;
  cout<<"temp_avg "<<temp_avg_a<<endl;
  for(int i=0;i<=max_label;i++)
  {
    cout<<i<<" "<<type_avg[i]<<endl;
  }
  for(int j=total_days-14;j<=total_days;j++)
  {
    if(sum_day[j]>4*temp_avg_a){
      for(int i=0;i<=max_label;i++)
      {
	T_data[j].flaver_num[i]=type_avg_end[i];
      }
    }
  }
}
*/

void sum_day_flaver()
{
    int temp_sum,temp_avg_a,today_sum;
    int threshold_deal_noise=3;
    //cyxadd
    
    for(int day=14;day<total_days;day++)
    {
      //calculate the avg of  14days before the day
      temp_sum=0;
      today_sum=0;
      for(int i=0;i<=max_label;i++)
      {
	  for(int j=day-14;j<day;j++)
	  {
	      temp_sum+=T_data[j].flaver_num[i];
	  }
	  today_sum+=T_data[day].flaver_num[i];
      }
      temp_sum/=14;
      cout<<"day"<<day<<' '<<" avg is "<<temp_sum<<' '<<today_sum<<endl;
      
      //fix the noise day to the avg of 14days before
      if(today_sum>threshold_deal_noise*temp_sum)
      {
	  for(int i=0;i<=max_label;i++)
	  {
	    temp_avg_a=0;
	    for(int j=day-14;j<day;j++)
		temp_avg_a+=T_data[j].flaver_num[i];
	    temp_avg_a/=14;
	    cout<<"fix the day "<<day<<" origin"<<i<<' '<<T_data[day].flaver_num[i]<<"to"<<temp_avg_a<<endl;
	    T_data[day].flaver_num[i]=temp_avg_a;
	  }
      }
      
    }
}

void clean_data(){

  for(int i=0;i<=max_label;i++)
  {
    for(int j=0;j<=total_days;j++)
    {
      cleandata[i][j].nums = T_data[j].flaver_num[i];
      cleandata[i][j].day = j;
      cleandata[i][j].need = 0; 
    } 
  }
  for(int i=0;i<=max_label;i++)
  {
    for(int j=0;j<=total_days;j++)
    {
      if(cleandata[i][j].nums>(type_avg[i]*3))
      {
	cout<<"type "<<i<<" day "<<j<<" be changed"<<endl;
	cleandata[i][j].need=1;
	T_data[j].flaver_num[i]=type_avg[i];
      }
    }
  }
}


void multi_cal(int type,int n_w,int val)
{
  int total_days_new=total_days-val;
  double weight[n_w+1],SquarePoor[4];
  n_w=n_w+1;
  cal_w[type].flaver_type=type;
  double mult_data[MAX_ROWS][n_w];
//create rows*n_w
  int row=0;
  while(row<=total_days_new-n_w)
  {
  
    for(int i=0;i<n_w;i++){
      mult_data[row][i]=T_data[i+row+1].flaver_num[type];
     cout<<mult_data[row][i]<<" ";
    }
    row++;
    cout<<endl;
  }
 // cout<<"row::"<<row<<endl;
  MultipleRegression((double*)mult_data,row,n_w,weight,SquarePoor);
  for(int i=0;i<n_w;i++){
    cal_w[type].flaver_w[i]=weight[i];
    cout<<type<<" weight "<<cal_w[type].flaver_w[i]<<" ";
  }
  cout<<endl;
}



void predict_cal(int type,int n_w,int val)
{
  int total_days_new = total_days-val;
  int cal_x[100]={0};
  int delta_times=pre_times();
  double temp_sum=0;
 // cout<<"type: "<<type<<endl;
  for(int i=0;i<n_w;i++)
  {
    cal_x[i]=T_data[total_days-n_w+i+1].flaver_num[type];
  }
  for(int i=0;i<delta_times;i++)
  {
    cout<<temp_sum<<":::";
    for(int j=0;j<n_w;j++)
    {
      temp_sum+=cal_w[type].flaver_w[j+1]*cal_x[j];
     // cout<<cal_x[j]<<"*";
      cout<<cal_x[j]<<"*"<<cal_w[type].flaver_w[j+1]<<" = "<<temp_sum<<"---";
    }
   cout<<" [0] "<<cal_seven_w[type].flaver_w[0]<<"  ";
    temp_sum+=cal_w[type].flaver_w[0];
    
    for(int k=0;k<n_w-1;k++)
    {
      double right_move= cal_x[k+1];
      cal_x[k] = right_move;
    }
    //cout<<temp_sum<<endl;
    if(temp_sum>=0)
    {
      if(temp_sum>3*type_avg[type])
      {
	temp_sum=type_avg[type];
	cal_x[n_w-1]=temp_sum;
        cal_w[type].flaver_ans[i]=int(temp_sum);
	temp_sum=0;
      }
      else{
        cal_x[n_w-1]=temp_sum;
        cal_w[type].flaver_ans[i]=int(temp_sum);
	temp_sum=0;
      }
      
    }
    else
    {
      cal_x[n_w-1]=0;
      cal_w[type].flaver_ans[i]=0;  
      temp_sum=0;
    }
    
  }
}



void view_weight(int n_w)
{
 for(int i=1;i<=max_label;i++)
 {
   cout<<"type: "<<cal_w[i].flaver_type<<endl;
   for(int j=0;j<n_w;j++)
   {
     cout<<" "<<cal_w[i].flaver_w[j];
   }
   cout<<endl;
 }
}

void view_ans(int n_w)
{
  
  int delta_times=pre_times();
  for(int j=1;j<=max_label;j++)
  {
    int type=j;
    for(int i=0;i<delta_times;i++)
    {
      cout<<cal_w[type].flaver_ans[i]<<" - ";
    }
    cout<<endl;
  }
}


void get_flaver_back()
{ 
  int type_avg_all[100]={0};
  int train_test_delta = train_test_time();
  cout<<"train_test: "<<train_test_delta<<endl;
  int temp_all=0;
  int delta_times=pre_times();
  cout<<"delta time:::"<<delta_times<<endl;
  double uprate=0;
  uprate = (double)(delta_times+train_test_delta)/(double)(delta_times);
  cout<<"uprate: "<<uprate<<endl;
  for(int j=1;j<=max_label;j++)
  {    
    type_avg_all[j]=type_avg[j]*delta_times;
    all_flaver[j]=0;
    int type=j;
    for(int i=0;i<delta_times;i++)
    {
      all_flaver[j]+=cal_w[type].flaver_ans[i];
    }
   cout<<"("<<all_flaver[j]<<"<->"<<delta_num_avg[j]<<"<->"<<type_avg_all[j]<<")";
    
  }
  cout<<endl;
  for(int j=1;j<=max_label;j++)
  {
    //use delta_avg to filter
    if(all_flaver[j]>delta_num_avg[j]*3)
    {
      all_flaver[j]=delta_num_avg[j]*2;
    }
    else if(all_flaver[j]<delta_num_avg[j]/2)
    {
      if(all_flaver[j]==0)
      {
	all_flaver[j]=delta_num_avg[j]/2;
      }
      else
      {
	all_flaver[j]=delta_num_avg[j]/1.2;
      }
      
    }
    else{
      temp_all=all_flaver[j];
      all_flaver[j]=(temp_all+delta_num_avg[j]+type_avg_all[j])/3;
      all_flaver[j]=all_flaver[j]*(uprate);
    }
  }
//



}

void get_delta_flaver(){
  int delta_time=pre_times();
  for(int j=1;j<=max_label;j++)
  {
    if(delta_num_avg[j]>=delta_time*3)
    {
      all_flaver[j] = delta_num_avg[j]*1.5;
    }
    else if(delta_num_avg[j]<delta_time*3&&delta_num_avg[j]>=delta_time)
    {
      all_flaver[j] = delta_num_avg[j]*1.2;
    }
    else
    {
      all_flaver[j] = delta_num_avg[j];
    }
    
  }
}



void view_allans()
{
  for(int j=1;j<=max_label;j++)
  {
    cout<<"("<<j<<":"<<all_flaver[j]<<")";
  }
  cout<<endl;
  
}
int get_stop()
{
  int ans=0;
  for(int i=0;i<=max_label;i++){
    ans+=p_flaver[i].numbers;
  }
  return ans;
}
/*
void view_bag_xiaoxi(int bags)
{
  int all_bag[100]={0};
  for(int i=0;i<bags;i++){
    cout<<"bags: "<<i+1<<": ";
    int temp=0;
    int other=0;
    for(int j=0;j<=max_label;j++)
    {
      all_bag[j]+=s_bag[i][j];
      temp+=s_bag[i][j]*p_flaver[j].memory;
      other+=s_bag[i][j]*p_flaver[j].kernel;
      cout<<s_bag[i][j]<<" ";
    }
    cout<<" mem: "<<temp <<" cpu: "<<other;
    cout<<endl;
  }

  for(int j=0;j<=max_label;j++)
  {
    cout<<all_bag[j]<<" ";
  }
  cout<<endl;
}
*/

void view_bag_xiaoxi_fusai()
{
  int all_bag[100]={0};
  
  for (int i = 0; i <out_bags_fusai[0]; i++) {
    cout<<"High-"<<i+1<<" ";
    int tempcpu=0,tempmem=0;
    for (int j = 0; j < max_label; j++) {
      all_bag[j]+=s_bag[0][i][j];
      tempcpu+=s_bag[0][i][j]*p_flaver[j].kernel;
      tempmem+=s_bag[0][i][j]*p_flaver[j].memory;
      cout<<s_bag[0][i][j]<<" ";  
    }
    cout<<" mem: "<<tempmem <<" cpu: "<<tempcpu;
    cout<<endl;
  }
  cout<<endl;
  for (int i = 0; i < out_bags_fusai[1]; i++) {
    cout<<"Generl-"<<i+1<<" ";
    int tempcpu=0,tempmem=0;
    for (int j = 0; j < max_label; j++) {
      all_bag[j]+=s_bag[1][i][j];
      tempcpu+=s_bag[1][i][j]*p_flaver[j].kernel;
      tempmem+=s_bag[1][i][j]*p_flaver[j].memory;
      cout<<s_bag[1][i][j]<<" ";  
    }
    cout<<" mem: "<<tempmem <<" cpu: "<<tempcpu;
    cout<<endl;
  }
  cout<<endl;
  for (int i = 0; i < out_bags_fusai[2]; i++) {
    cout<<"Large-"<<i+1<<" ";
    int tempcpu=0,tempmem=0;
    for (int j = 0; j < max_label; j++) {
      all_bag[j]+=s_bag[2][i][j];
      tempcpu+=s_bag[2][i][j]*p_flaver[j].kernel;
      tempmem+=s_bag[2][i][j]*p_flaver[j].memory;
      cout<<s_bag[2][i][j]<<" ";  
    }
    cout<<" mem: "<<tempmem <<" cpu: "<<tempcpu;
    cout<<endl;
  }
  cout<<endl;
  for(int j=0;j<=max_label;j++)
  {
    cout<<all_bag[j]<<" ";
  }
  cout<<endl;
}


void avg_data_end(){

  for(int i=0;i<=max_label;i++)
  {
    type_avg_end[i]=0;
    for(int j=total_days-14;j<=total_days;j++)
    {
      
      type_avg_end[i]+=T_data[j].flaver_num[i];
    }
    type_avg_end[i]= type_avg_end[i]/15;
    cout<<" type_avg_end: "<<i<<" avg: "<<type_avg_end[i]<<endl;
  } 
}



void avg_data(){

  for(int i=0;i<=max_label;i++)
  {
    for(int j=0;j<=total_days;j++)
    {
       
      type_avg[i]+=T_data[j].flaver_num[i];
    }
    type_avg[i]= type_avg[i]/total_days+1;
    cout<<" avg_flaver: "<<i<<" avg: "<<type_avg[i]<<endl;
  } 
}


int cal_sum()
{
  int sum=0;
  for(int i=0;i<=S_input.flaver_num;i++)
  {
    sum+=p_flaver[i].numbers;
  }
  return sum;
}

int put_predict_xiaoxi_fusai(){
  int serverCPU[3]={0},serverMEM[3]={0};
  
  string serverName;
  int numls[18]={0},cpuls[18]={0},memls[18]={0},Nfla=0;
  Nfla = S_input.flaver_num;
  int sum_test=0;
  for(int i=0;i<S_input.flaver_num;i++)
  {
    p_flaver[i].type=S_input.flaver[i].flaver_label;
    p_flaver[i].name=S_input.flaver[i].name;
    p_flaver[i].kernel=S_input.flaver[i].kernel;
    p_flaver[i].memory=S_input.flaver[i].memory;
    p_flaver[i].numbers = all_flaver[S_input.flaver[i].flaver_label];
  }
  for(int i=0;i<S_input.flaver_num;i++){
    numls[i] = all_flaver[S_input.flaver[i].flaver_label];
    sum_test+=numls[i];
    cpuls[i] = S_input.flaver[i].kernel;
    memls[i] = S_input.flaver[i].memory;
  }
  for(int i=0;i<3;i++)
  {
    if(S_input.physicalserver[i].name=="High-Performance")
    {
      out_bags_fusai_name[0]="High-Performance";
      cout<<"High-Performance"<<endl;
      serverCPU[0]=S_input.physicalserver[i].kernel_num;
      serverMEM[0]=S_input.physicalserver[i].memory;
      cout<<"server: "<<serverCPU[0]<<" "<<serverMEM[0]<<endl;
    }
    else if(S_input.physicalserver[i].name=="General")
    {
      out_bags_fusai_name[1]="General";
      cout<<"General"<<endl;
      serverCPU[1]=S_input.physicalserver[i].kernel_num;
      serverMEM[1]=S_input.physicalserver[i].memory;
      cout<<"server: "<<serverCPU[1]<<" "<<serverMEM[1]<<endl;
    }
    else if(S_input.physicalserver[i].name=="Large-Memory")
    {
      out_bags_fusai_name[2]="Large-Memory";
      cout<<"Large-Memory"<<endl;
      serverCPU[2]=S_input.physicalserver[i].kernel_num;
      serverMEM[2]=S_input.physicalserver[i].memory;
      cout<<"server: "<<serverCPU[2]<<" "<<serverMEM[2]<<endl;
    }
  }
  for(int i=0;i<3;i++){
    serverCPUtemp[i]=serverCPU[i];
    serverMEMtemp[i]=serverMEM[i];
  }
  binpacking bp1(cpuls,memls,numls,serverCPU,serverMEM,Nfla);
  cout<<"bp1.totalCPU: "<<bp1.totalCPU<<" bp1.totalMEM: "<<bp1.totalMEM<<endl;
  if(bp1.Nall<1000)
	cout<<"USING GeneticAlgorithm -------- Score= "<<bp1.GeneticA_adaptive(6000,150)<<" --------"<<endl;
  else
	cout<<"USING GeneticAlgorithm -------- Score= "<<bp1.GeneticA_adaptive(2000,150)<<" --------"<<endl;  

  cout << endl << "High performance" << endl;
  for (int i = 0; i < bp1.ans_bagnum[0]; i++) {
    cout << i + 1 << ' ';
    for (int j = 0; j < 18; j++) {
      s_bag_temp[0][i][j]=0;
      if (bp1.ans_ff[0][i][j] > 0) {
	  s_bag_temp[0][i][j] = bp1.ans_ff[0][i][j];
	  cout << "flavor" << j + 1 << ' ' << bp1.ans_ff[0][i][j] << ' ';
      }
    }
    cout << endl;
  }
  cout << endl << "General" << endl;
  for (int i = 0; i < bp1.ans_bagnum[1]; i++) {
    cout << i + 1 << ' ';
    for (int j = 0; j < 18; j++) {
      s_bag_temp[1][i][j]=0;
      if (bp1.ans_ff[1][i][j] > 0) {
	s_bag_temp[1][i][j] = bp1.ans_ff[1][i][j];
	cout << "flavor" << j + 1 << ' ' << bp1.ans_ff[1][i][j] << ' ';
      }
    }
    cout << endl;
  }
  cout << endl << "Large memory" << endl;
  for (int i = 0; i < bp1.ans_bagnum[2]; i++) {
    cout << i + 1 << ' ';
    for (int j = 0; j < 18; j++) {
      s_bag_temp[2][i][j]=0;
      if (bp1.ans_ff[2][i][j] > 0) {
	s_bag_temp[2][i][j] = bp1.ans_ff[2][i][j];
	cout << "flavor" << j + 1 << ' ' << bp1.ans_ff[2][i][j] << ' ';
      }
    }
    cout << endl;
  }
  cpu_all_use_temp = bp1.totalCPU;
  mem_all_use_temp = bp1.totalMEM;  
  out_sum_flavers_temp = cal_sum();
  for(int i=0;i<3;i++){
    out_bags_fusai[i]=bp1.ans_bagnum[i];
    cout<<"outbag "<<i<<" num: "<<out_bags_fusai[i]<<" "<<serverCPUtemp[i]<<" "<<serverMEMtemp[i]<<endl;
  }

  all_can_CPU = serverCPUtemp[0]*out_bags_fusai[0]+serverCPUtemp[1]*out_bags_fusai[1]+serverCPUtemp[2]*out_bags_fusai[2];
  all_can_MEM = serverMEMtemp[0]*out_bags_fusai[0]+serverMEMtemp[1]*out_bags_fusai[1]+serverMEMtemp[2]*out_bags_fusai[2];
  cout<<"all_can "<<all_can_CPU<<" "<<all_can_MEM<<endl;
  return bp1.ans_bagnum[0]+bp1.ans_bagnum[1]+bp1.ans_bagnum[2];
  
}
/*
int put_predict_xiaoxi(){
  int serverCPU=0,serverMEM=0;
  string caretarget;
  int sum_test=0;
  int numls[18]={0},cpuls[18]={0},memls[18]={0},Nfla=0;
  for(int i=0;i<S_input.flaver_num;i++)
  {
    p_flaver[i].type=S_input.flaver[i].flaver_label;
    p_flaver[i].name=S_input.flaver[i].name;
    p_flaver[i].kernel=S_input.flaver[i].kernel;
    p_flaver[i].memory=S_input.flaver[i].memory;
    p_flaver[i].numbers = all_flaver[S_input.flaver[i].flaver_label];
  }
  for(int i=0;i<S_input.flaver_num;i++){
    numls[i] = all_flaver[S_input.flaver[i].flaver_label];
    sum_test+=numls[i];
    cpuls[i] = S_input.flaver[i].kernel;
    memls[i] = S_input.flaver[i].memory;
  }
  Nfla = S_input.flaver_num;
  serverCPU = S_input.physicalserver[0].kernel_num;
  serverMEM = S_input.physicalserver[0].memory;
  if(S_input.allottype==0)
  { caretarget="CPU"; }
  else
  { caretarget="MEM"; }
  binpacking bp1(cpuls,memls,numls,serverCPU,serverMEM,Nfla,caretarget);
  cout<<bp1.totalCPU<<' '<<bp1.totalMEM<<endl;
  cout<<"USING GeneticAlgorithm -------- Need "<<bp1.GeneticA(1000,100)<<" servers--------"<<endl;
  
  for (int i=0;i<bp1.ans_bagnum;i++)
  {
//      cout<<i+1<<' ';
      for (int j=0;j<S_input.flaver_num;j++)
      {
	s_bag_temp[i][j]=0;
	if (bp1.ans_ff[i][j]>0)
	{
  //	  cout<<S_input.flaver[j].name<<' '<<bp1.ans_ff[i][j]<<' ';
	    s_bag_temp[i][j] = bp1.ans_ff[i][j];
	} 
  }
//  cout<<endl;
  }
  
  cpu_all_use_temp = bp1.totalCPU;
  mem_all_use_temp = bp1.totalMEM;
  out_sum_flavers_temp = cal_sum();
 // cout<<"sum_flavers: "<<sum_test<<endl;
//  view_bag_xiaoxi(bp1.ans_bagnum);
  return bp1.ans_bagnum;
} 
*/
void select_change_xiaoxi(){
  int temp_out_bags;
  int iters=1;
  int nums_fla_inlastbag=0;
  out_bags=99999;
  int min_fla_inlastbag=999999;
  while(iters)
  {
    nums_fla_inlastbag=0;
    srand(iters*111+2018);  
    iters--;
    temp_out_bags = put_predict_xiaoxi_fusai();
    
    if(temp_out_bags<=out_bags){
      cout<<"iters: "<<iters<<" temp_bags: "<<temp_out_bags<<endl<<endl;
      out_bags=temp_out_bags;
      cpu_all_use = cpu_all_use_temp;
      mem_all_use = mem_all_use_temp;
      out_sum_flavers = out_sum_flavers_temp; 
      
      for(int k=0;k<3;k++)
      {
	for (int i = 0; i <out_bags_fusai[k]; i++) {
	  for (int j = 0; j < 18; j++) {
	    s_bag[k][i][j]=0;
	    if (s_bag_temp[k][i][j] > 0) {
		s_bag[k][i][j] = s_bag_temp[k][i][j];
	    }
	  }
	}
	
      }
       
      
      //cyxadd
 /*
      for (int j=0;j<S_input.flaver_num;j++)
      {
	  nums_fla_inlastbag+=s_bag_temp[temp_out_bags-1][j];
      }
      cout<<"in last bag has "<<nums_fla_inlastbag<<endl;
      if(nums_fla_inlastbag<min_fla_inlastbag)
      {
	cout<<"fix min_fla_inlastbag from "<<min_fla_inlastbag<<" to "<<nums_fla_inlastbag<<endl;
	min_fla_inlastbag=nums_fla_inlastbag;
	
	out_bags=temp_out_bags;
	cpu_all_use = cpu_all_use_temp;
	mem_all_use = mem_all_use_temp;
	out_sum_flavers = out_sum_flavers_temp;
	for (int i=0;i<out_bags;i++)
	{
          
	    for (int j=0;j<S_input.flaver_num;j++)
	    {
	      s_bag[i][j]=0;
	      if (s_bag_temp[i][j]>0)
	      {
		  s_bag[i][j] = s_bag_temp[i][j];
	      } 
	    }
	}
 
      }*/
    }
  }
  
  
//  view_bag_xiaoxi(out_bags);
  cout<<endl;
}


void view_rate()
{
  double cpu_rate,mem_rate;
  cpu_rate = (double)cpu_all_use/(double)(all_can_CPU);
  mem_rate = (double)mem_all_use/(double)(all_can_MEM);
  cout<<"cpu_use: "<<cpu_all_use<<" mem_use: "<<mem_all_use<<" cpu_rate:: "<<cpu_rate<<" mem_rate::"<<mem_rate<<endl;
}

void update_flaver()
{
  double cpu_rate=0,mem_rate=0,temp_cpu,temp_mem;
  int iters=1;
  s_bag_zero();
  while(iters)
  {
    iters--;
    view_allans();
    select_change_xiaoxi(); 
    temp_cpu=cpu_rate;
    temp_mem=mem_rate;
    cpu_rate = (double)cpu_all_use/(double)(all_can_CPU);
    mem_rate = (double)mem_all_use/(double)(all_can_MEM);
    cout<<" cpu_rate:: "<<cpu_rate<<" mem_rate::"<<mem_rate<<endl;
    view_bag_xiaoxi_fusai();
  }
  float baserate=1;
//delete last one
  if(out_bags_fusai[0])//High-Performance
  {
    int last_bag=0,last_mem=0,last_cpu=0;
    last_bag = out_bags_fusai[0]-1;
    for(int j=0;j<S_input.flaver_num;j++)
    {
      last_mem+=s_bag[0][last_bag][j]*p_flaver[j].memory;
      last_cpu+=s_bag[0][last_bag][j]*p_flaver[j].kernel;
    }
    float rate=0;
    rate = (float)(last_cpu+last_mem)/(float)(serverCPUtemp[0]+serverMEMtemp[0]);
    if(rate<baserate)
    {      
      for(int i=0;i<=max_label;i++)
      {
	if(s_bag[0][last_bag][i]!=0)
	{
	  all_bag_delete+=s_bag[0][last_bag][i];
	}
      }
      out_bags_fusai[0]-=1;
      ifdelete=1;
    }
    mem_all_use-=last_mem;
    cpu_all_use-=last_cpu;
  }
  if(out_bags_fusai[1])//General
  {
    int last_bag=0,last_mem=0,last_cpu=0;
    last_bag = out_bags_fusai[1]-1;
    for(int j=0;j<S_input.flaver_num;j++)
    {
      last_mem+=s_bag[1][last_bag][j]*p_flaver[j].memory;
      last_cpu+=s_bag[1][last_bag][j]*p_flaver[j].kernel;
    }
    float rate=0;
    rate = (float)(last_cpu+last_mem)/(float)(serverCPUtemp[1]+serverMEMtemp[1]);
    if(rate<baserate)
    {      
      for(int i=0;i<=max_label;i++)
      {
	if(s_bag[1][last_bag][i]!=0)
	{
	  all_bag_delete+=s_bag[1][last_bag][i];
	}
      }
      out_bags_fusai[1]-=1;
      ifdelete=1;
    }
    mem_all_use-=last_mem;
    cpu_all_use-=last_cpu;
  }
  if(out_bags_fusai[2])//Large-Memory
  {
    int last_bag=0,last_mem=0,last_cpu=0;
    last_bag = out_bags_fusai[2]-1;
    for(int j=0;j<S_input.flaver_num;j++)
    {
      last_mem+=s_bag[2][last_bag][j]*p_flaver[j].memory;
      last_cpu+=s_bag[2][last_bag][j]*p_flaver[j].kernel;
    }
    float rate=0;
    rate = (float)(last_cpu+last_mem)/(float)(serverCPUtemp[2]+serverMEMtemp[2]);
    if(rate<baserate)
    {      
      for(int i=0;i<=max_label;i++)
      {
	if(s_bag[2][last_bag][i]!=0)
	{
	  all_bag_delete+=s_bag[2][last_bag][i];
	}
      }
      out_bags_fusai[2]-=1;
      ifdelete=1;
    }
    mem_all_use-=last_mem;
    cpu_all_use-=last_cpu;    
  }
//update all CPU MEM  
  all_can_CPU = serverCPUtemp[0]*out_bags_fusai[0]+serverCPUtemp[1]*out_bags_fusai[1]+serverCPUtemp[2]*out_bags_fusai[2];
  all_can_MEM = serverMEMtemp[0]*out_bags_fusai[0]+serverMEMtemp[1]*out_bags_fusai[1]+serverMEMtemp[2]*out_bags_fusai[2];

/*
  int last_bag=0,last_mem=0,last_cpu=0;
  last_bag = out_bags-1;
  for(int j=0;j<S_input.flaver_num;j++)
  {
    last_mem+=s_bag[last_bag][j]*p_flaver[j].memory;
    last_cpu+=s_bag[last_bag][j]*p_flaver[j].kernel;
  }
  if(1) 
  {
    
    float rate=0;
    rate = float(last_cpu+last_mem)/float(S_input.physicalserver[0].kernel_num+S_input.physicalserver[0].memory);
    if(rate<0.6)
    {
      for(int i=0;i<=max_label;i++)
      {
	if(s_bag[out_bags-1][i]!=0)
	{
	  all_bag_delete+=s_bag[out_bags-1][i];
	}
      }
      out_bags-=1;
      ifdelete=1;
    }
  }
 
  mem_all_use-=last_mem;
  cpu_all_use-=last_cpu;
*/
}

/*
void out_file(char * filename)
{
  
  ofstream write;
  write.open(filename);
  write<<out_sum_flavers-all_bag_delete<<'\n';
  for(int i=0;i<S_input.flaver_num;i++)
  {
    if(ifdelete)
    {
      cout<<s_bag[out_bags][i]<<endl;
      write<<p_flaver[i].name<<" "<<p_flaver[i].numbers-s_bag[out_bags][i]<<'\n';
     }
    else
    {
      write<<p_flaver[i].name<<" "<<p_flaver[i].numbers<<'\n';
    }
  }
  write<<'\n';
  write<<S_input.physicalserver[0].name<<' '<<out_bags<<'\n';
  for(int i=1;i<=out_bags;i++)
  {
    write<<S_input.physicalserver[0].name<<'-'<<i;
    for(int j=0;j<max_label;j++)
    {
      if(s_bag[i-1][j]!=0)
      {
	write<<" ";
      }
      if(s_bag[i-1][j]!=0)
      {
	write<<p_flaver[j].name<<" "<<s_bag[i-1][j];
      }
    }
    write<<'\n';
  }
  
  write.close();
 
}
*/


void out_file_fusai(char * filename)
{
  ofstream write;
  write.open(filename);
  write<<out_sum_flavers-all_bag_delete<<'\n';
  for(int i=0;i<S_input.flaver_num;i++)
  {
    if(ifdelete)
    {
      cout<<s_bag[0][out_bags_fusai[0]][i]<<" "<<s_bag[1][out_bags_fusai[1]][i]<<" "<<s_bag[2][out_bags_fusai[2]][i]<<endl;
      write<<p_flaver[i].name<<" "<<p_flaver[i].numbers-s_bag[0][out_bags_fusai[0]][i]-s_bag[1][out_bags_fusai[1]][i]-s_bag[2][out_bags_fusai[2]][i]<<'\n';
    }
    else
    {
      cout<<s_bag[0][out_bags_fusai[0]][i]<<" "<<s_bag[1][out_bags_fusai[1]][i]<<" "<<s_bag[2][out_bags_fusai[2]][i]<<endl;
      write<<p_flaver[i].name<<" "<<p_flaver[i].numbers<<'\n';
    }
  }
  write<<'\n';
  if(out_bags_fusai[1])
  {
    write<<out_bags_fusai_name[1]<<' '<<out_bags_fusai[1]<<'\n';
    for(int i=1;i<=out_bags_fusai[1];i++)
    {
      write<<out_bags_fusai_name[1]<<'-'<<i;
      for(int j=0;j<max_label;j++)
      {
	if(s_bag[1][i-1][j]!=0)
	{
	  write<<" ";
	}
	if(s_bag[1][i-1][j]!=0)
	{
	  write<<p_flaver[j].name<<" "<<s_bag[1][i-1][j];
	}
      }
      write<<'\n';
    }
    write<<'\n';
  }
  if(out_bags_fusai[2])
  {
    write<<out_bags_fusai_name[2]<<' '<<out_bags_fusai[2]<<'\n';
    for(int i=1;i<=out_bags_fusai[2];i++)
    {
      write<<out_bags_fusai_name[2]<<'-'<<i;
      for(int j=0;j<max_label;j++)
      {
	if(s_bag[2][i-1][j]!=0)
	{
	  write<<" ";
	}
	if(s_bag[2][i-1][j]!=0)
	{
	  write<<p_flaver[j].name<<" "<<s_bag[2][i-1][j];
	}
      }
      write<<'\n';
    }
    write<<'\n'; 
  }
  if(out_bags_fusai[0])
  {
    write<<out_bags_fusai_name[0]<<' '<<out_bags_fusai[0]<<'\n';
    for(int i=1;i<=out_bags_fusai[0];i++)
    {
      write<<out_bags_fusai_name[0]<<'-'<<i;
      for(int j=0;j<max_label;j++)
      {
	if(s_bag[0][i-1][j]!=0)
	{
	  write<<" ";
	}
	if(s_bag[0][i-1][j]!=0)
	{
	  write<<p_flaver[j].name<<" "<<s_bag[0][i-1][j];
	}
      }
      write<<'\n';
    }
  }
  write.close();
}

void predict_server(char * info[MAX_INFO_NUM], char * data[MAX_DATA_NUM], int data_num, char * filename)
{
  read_info(info);
  read_data(data,data_num);
  view_train();
  avg_data();
  avg_data_end();
  sum_day_flaver();
  clean_data();
  get_delta_train();
  for(int i=1;i<=max_label;i++){
    multi_cal(i,10,0);
  }
  for(int i=1;i<=max_label;i++){
    predict_cal(i,10,0); 
  }
  view_ans(10);
  
  get_flaver_back();  
  
  update_flaver();
  out_file_fusai(filename); 
}

