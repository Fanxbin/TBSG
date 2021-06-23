#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <time.h>
#include <stdlib.h>
#include "TBSG.h"
using namespace std;
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
int main(int argc,char** argv){
    if(argc!=7){
         std::cout<<"Error with parameters"<<std::endl;
         exit(-1);
    }
    int vecsize,vecdim;
    float*data;
    load_fvecs(argv[1],data,vecsize,vecdim);
    vector<float*>mass(vecsize,NULL);
    for(int i=0;i<vecsize;i++) mass[i]=data+i*vecdim;
    int M,S;
    string nnfile;
    float MC=0.5;
    M=atoi(argv[2]);
    S=atoi(argv[3]);
    MC=atof(argv[4]);
    nnfile=argv[5];
    std::cout<<"Data_path = "<<argv[1]<<" vecsize = "<<vecsize<<" vecdim = "<<vecdim<<endl;
    std::cout<<"M = "<<M<<std::endl;
    std::cout<<"S = "<<S<<std::endl;
    std::cout<<"MC = "<<MC<<std::endl;
    std::cout<<"NNfile = "<<nnfile<<std::endl;
    std::cout<<"Save_path = "<<argv[6]<<std::endl;
    TBSG* tbsg=new TBSG(mass,data,M,vecdim,S,MC,nnfile);
    tbsg->TBSG_save(argv[6]);
    return 0;
}
