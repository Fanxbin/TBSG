#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <time.h>
#include <stdlib.h>
#include <queue>
#include "TBSG.h"
using namespace std;
#define dist_t float
void load_fvecs(string filename,float*&data,int&num,int&dim){
    ifstream fin(filename,ios::binary);
    if(!fin.is_open()){
        cout<<"open "<<filename<<" error\n";
        exit(-1);
    }
    fin.read((char*)&dim,4);
    fin.seekg(0,ios::end);
    streampos sp=fin.tellg();
    num=(long long)(sp)/(dim+1)/4;
    data=new float[num*dim];
    fin.seekg(0,ios::beg);
    for(int i=0;i<num;i++){
        fin.seekg(4,ios::cur);
        fin.read((char*)(data+i*dim),dim*sizeof(float));
    }
    fin.close();
}
void load_ivecs(string filename,vector<vector<int> >&data){
    ifstream fin(filename,ios::binary);
    if(!fin.is_open()){
        cout<<"open "<<filename<<" error\n";
        exit(-1);
    }
    int num,dim;
    fin.read((char*)&dim,4);
    fin.seekg(0,ios::end);
    streampos sp=fin.tellg();
    num=(long long)(sp)/(dim+1)/4;
    data.resize(num,vector<int>(dim,-1));
    int nums[dim];
    fin.seekg(0,ios::beg);
    for(int i=0;i<num;i++){
        fin.seekg(4,ios::cur);
        fin.read((char*)nums,dim*4);
        for(int j=0;j<dim;j++) data[i][j]=nums[j];
    }
    fin.close();
}
int main(int argc,char**argv){
    if(argc!=7){
        cout<<"Error with parameters"<<endl;
        exit(-1);
    }
    int vecsize,vecdim,qsize,k;
    float* data;
    load_fvecs(argv[1],data,vecsize,vecdim);
    std::cout<<"Data_path = "<<argv[1]<<std::endl;
    float* query;
    load_fvecs(argv[2],query,qsize,vecdim);
    std::cout<<"Query_path = "<<argv[2]<<std::endl;
    vector<float*>massQ(qsize,NULL);
    for(int i=0;i<qsize;i++) massQ[i]=query+i*vecdim;
    vector<vector<int> >massQA;
    load_ivecs(argv[3],massQA);
    std::cout<<"GT_path = "<<argv[3]<<std::endl;
    k=massQA[0].size();
    vector<std::priority_queue<std::pair<dist_t,int> > >result;
    TBSG * tbsg=new TBSG(data,vecdim,argv[4]);
    int step=atoi(argv[5]);
    vector<pair<float,float> >record;
    int fss;
    int qk=1;
    for(int i=0;i<=100;i++){
        //fss=k+i*step;
        if(i<5) fss=qk+i*step;
        else if(i<10) fss+=2*step;
        else fss+=2*step;
        cout<<"fss_="<<fss<<endl;
        auto s = std::chrono::high_resolution_clock::now();
        tbsg->SearchKnn(massQ,k,result,fss);
        auto e = std::chrono::high_resolution_clock::now();
 
        std::chrono::duration<double> diff = e - s;
        float prec=precision_calc(massQA,result,qk);
        float sec_cost=diff.count();
        printf("TBSG precision=%.6f,costing %.6f seconds, QPS=%.6f\n",prec,sec_cost,1.0f*qsize/sec_cost);
        record.push_back(make_pair(prec,1.0f*qsize/sec_cost));
        //record.push_back(make_pair(prec,vecsize/(cal_cnt/qsize)));
        record_output(argv[6],record);
    }
    return 0;
}
