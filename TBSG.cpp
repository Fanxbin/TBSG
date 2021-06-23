#include <iostream>
#include <fstream>
#include <queue>
#include "cover_tree.h"
#include <set>
#include <time.h>
#include <stdio.h>
#include <thread>
#include <stack>
#include <map>
#include <stdlib.h>
#include <boost/dynamic_bitset.hpp>
#include <bitset>
#include <chrono>
#include <sstream>
#include <random>
#include "space_l2.h"
#include "TBSG.h"
#include <algorithm>
#define Node CoverTree::Node
#define dist_t float
#define LL long long
#define PI 3.1415926
#define Tiny 1e-10
using namespace std;
using namespace hnswlib;
atomic<long long>cal_cnt;
DISTFUNC<float> fstdistfunc_;
void* dist_func_param_;
hnswlib::L2Space* l2space;
float l2dist(float*a,float*b,int dim){
    return fstdistfunc_(a,b,dist_func_param_);
    cal_cnt.fetch_add(1);
    float ret=0;
    for(int i=0;i<dim;i++){
        float t=a[i]-b[i];
        ret+=t*t;
    }
    return ret;
}


void GenRandom(std::mt19937& rng, unsigned* addr, unsigned size, unsigned N) {
  for (unsigned i = 0; i < size; ++i) {
    addr[i] = rng() % (N - size);
  }
  std::sort(addr, addr + size);
  for (unsigned i = 1; i < size; ++i) {
    if (addr[i] <= addr[i - 1]) {
      addr[i] = addr[i - 1] + 1;
    }
  }
  unsigned off = rng() % N;
  for (unsigned i = 0; i < size; ++i) {
    addr[i] = (addr[i] + off) % N;
  }
}

void TBSG::covertree_dfs(Node* par){
    if(par==NULL) return;
    for(int i=0;i<par->dataId.size();i++) nodelist[par->dataId[i]]=par;
    for(int i=0;i<par->children.size();i++) covertree_dfs(par->children[i]);
}
int TBSG::Did(Node* self){
    return self->dataId[0];
}
vector<int> TBSG::all_ids(Node* self){
    return self->dataId;
}
dist_t TBSG::Dist(int fir,int sec){
    return l2dist(data_+fir*dim,data_+sec*dim,dim);
}
dist_t TBSG::Dist(int id,float* p){
    return l2dist(data_+id*dim,p,dim);
}
float TBSG::cover_prob(float a,float b){
    return 1-acos(sin(2*a+b)/sin(a+b+Tiny)/2)/PI;
}
float cover_prob_2(float a,float b, float l,float r){
    return 1-acos(sin(2*a+b)/sin(a+b+Tiny)/2*l/r)/PI;
}
bool TBSG::occludeByNeighbor(dist_t neigh_dist,dist_t curr_dist,dist_t cand_dist,float min_cover){
    if(neigh_dist==0||cand_dist==0) return false;
    if(curr_dist==0) return true;
    if(neigh_dist<=cand_dist&&curr_dist<=cand_dist){
        float alpha=acos((neigh_dist+cand_dist-curr_dist)/2/sqrt(neigh_dist*cand_dist));
        float theta=acos((curr_dist+cand_dist-neigh_dist)/2/sqrt(curr_dist*cand_dist));
        if(cover_prob(alpha,theta)>=min_cover) return true;
    }
    return false;
}
bool TBSG::occludeByNeighbor(dist_t neigh_dist,dist_t curr_dist,dist_t cand_dist,float min_cover,float r){
    if(neigh_dist==0||cand_dist==0) return false;
    if(curr_dist==0) return true;
    if(neigh_dist<=cand_dist&&curr_dist<=cand_dist){
        float alpha=acos((neigh_dist+cand_dist-curr_dist)/2/sqrt(neigh_dist*cand_dist));
        float theta=acos((curr_dist+cand_dist-neigh_dist)/2/sqrt(curr_dist*cand_dist));
        if(cover_prob_2(alpha,theta,sqrt(cand_dist),r)>=min_cover) return true;
    }
    return false;
}
void TBSG::NeighborSelection(int self,vector<int>candidate,dist_t* Neigh_dist,float min_c){
    vector<pair<float,int> >candidates;
    for(int i=0;i<candidate.size();i++) candidates.push_back({Dist(self,candidate[i]),candidate[i]});
    sort(candidates.begin(),candidates.end());
    vector<pair<float,int>>res;
    for(int i=0;i<candidates.size();i++){
        bool exclude=false;
        if(res.size()==Max_neigh) break;
        for(int j=0;j<res.size();j++){
             if(occludeByNeighbor(res[j].first,Dist(res[j].second,candidates[i].second),candidates[i].first,min_c/*,radius[candidates[i].second]*/)){
                 exclude=true;
                 break;
             }
        }
        if(!exclude) res.push_back(candidates[i]);
    }
    nodelist[self]->mut.lock();
    int* neighbor=neighbors[self];
    neighbor[0]=res.size();
    for(int i=0;i<res.size();i++){
        neighbor[i+1]=res[i].second;
        *(Neigh_dist+self*Max_neigh+i)=res[i].first;
    }
    nodelist[self]->mut.unlock();
}
void TBSG::add_neigh(int self,int new_neigh,dist_t *Neigh_dist,float min_c){
    if(self==new_neigh) return;
    nodelist[self]->mut.lock();
    int* neighbor=neighbors[self];
    int nsize=neighbor[0];
    neighbor++;
    for(int i=0;i<nsize;i++) if(neighbor[i]==new_neigh){
        nodelist[self]->mut.unlock();
        return;
    }
    dist_t new_dist=Dist(self,new_neigh);
    dist_t* ndist=Neigh_dist+self*Max_neigh;
    int temp=nsize;
    bool add=true;
    for(int i=0;i<nsize;i++){
        dist_t neigh_dist=ndist[i];
        if(neigh_dist<=new_dist){
            if(occludeByNeighbor(neigh_dist,Dist(neighbor[i],new_neigh),new_dist,min_c)){
                add=false;
                break;
            }
        }
        else{
            temp=i;
            break;
        }

    }
    if(!add){
        nodelist[self]->mut.unlock();
        return;
    }
    vector<std::pair<int,dist_t> >result;
    for(int i=0;i<temp;i++) result.push_back(std::make_pair(neighbor[i],ndist[i]));
    result.push_back(std::make_pair(new_neigh,new_dist));
    for(int i=temp;i<nsize;i++){
        if(occludeByNeighbor(new_dist,Dist(neighbor[i],new_neigh),ndist[i],min_c)) continue;
        result.push_back(std::make_pair(neighbor[i],ndist[i]));
    }
    while(result.size()>Max_neigh) result.pop_back();
    neighbor=neighbors[self];
    neighbor[0]=result.size();
    for(int i=0;i<result.size();i++){
        neighbor[i+1]=result[i].first;
        ndist[i]=result[i].second;
    }
    nodelist[self]->mut.unlock();
}
void TBSG::select_candidates(Node* self,Node* par,std::priority_queue<std::pair<dist_t,int> >&candidates,int cand_size,boost::dynamic_bitset<> &flags){
    while(!candidates.empty()) candidates.pop();
    if(par==NULL) return;
    dist_t dist_bound=std::numeric_limits<dist_t>::max();
    priority_queue<std::pair<dist_t,int>,vector<std::pair<dist_t,int> >,greater<std::pair<dist_t,int> > >neigh_set;
    float* self_point=self->data;
    neigh_set.emplace(Dist(Did(par),self_point),Did(par));
    flags[Did(par)]=true;
    while(!neigh_set.empty()){
        int curr_neigh=neigh_set.top().second;
        dist_t curr_dist=neigh_set.top().first;
        neigh_set.pop();
        if(curr_dist>dist_bound) break;
        vector<Node*>children=nodelist[curr_neigh]->children;
        int csize=children.size();
        int ch_id[csize+1];
        ch_id[csize]=0;
        for(int i=0;i<csize;i++) ch_id[i]=Did(children[i]);
#ifdef USE_SSE
        _mm_prefetch(data_+ch_id[0]*dim, _MM_HINT_T0);;
#endif
        for(int i=0;i<csize;i++){
            if(children[i]==self) continue;
#ifdef USE_SSE
            _mm_prefetch(data_+ch_id[i+1]*dim, _MM_HINT_T0);
#endif
            dist_t ch_dist=l2dist(data_+ch_id[i]*dim,self_point,dim);
            if(ch_dist>dist_bound) continue;
            candidates.emplace(ch_dist,ch_id[i]);
            if(candidates.size()>cand_size){
                dist_bound=candidates.top().first;
                candidates.pop();
            }
        }
        int* neigh_neigh=neighbors[curr_neigh];
        int next_size=neigh_neigh[0];
        int next_neigh[next_size+1];
        for(int i=0;i<next_size;i++) next_neigh[i]=neigh_neigh[i+1];
        next_neigh[next_size]=0;
#ifdef USE_SSE
        _mm_prefetch(data_+next_neigh[0]*dim, _MM_HINT_T0);
#endif
        for(int i=0;i<next_size;i++){
            int nn_id=next_neigh[i];
            if(nodelist[nn_id]->children.size()==0) continue;
            if(flags[nn_id]==true) continue;
#ifdef USE_SSE
        _mm_prefetch(data_+next_neigh[i+1]*dim, _MM_HINT_T0);
#endif
            flags[nn_id]=true;
            neigh_set.emplace(l2dist(data_+nn_id*dim,self_point,dim),nn_id);
        }
    }
}


void Load(const char *filename,int**&data,int&num,int&dim) {
  std::ifstream in(filename, std::ios::binary);
  unsigned k;
  in.read((char*)&k,4);
  dim=k;
  in.seekg(0,std::ios::end);
  std::ios::pos_type ss = in.tellg();
  size_t fsize = (size_t)ss;
  num = fsize / ((size_t)k + 1) / 4;
  in.seekg(0,std::ios::beg);
  data=new int*[num];
  for(int i=0;i<num;i++) data[i]=new int[dim];
  for(size_t i = 0; i < num; i++){
    in.seekg(4,std::ios::cur);
    in.read((char*)data[i], k * sizeof(unsigned));
  }
  in.close();
}
void TBSG::insert(int* nn,int num,int curr,float mc,dist_t* Neigh_dist){
    for(int i=0;i<num;i++){
        add_neigh(curr,nn[i],Neigh_dist,mc);
        add_neigh(nn[i],curr,Neigh_dist,mc);
    }
}
void TBSG::TBSG_from_covertree(){
    int cnt=0;
    vector<Node*>prev,next;
    for(int i=0;i<dataNum;i++) neighbors[i][0]=0;
    dist_t* Neigh_dist=new dist_t[dataNum*Max_neigh];
    prev.push_back(covertree->root);
    long long prev_cnt=0;
    cal_cnt=0;
    /*
    while(true){
        if(prev.empty()) break;
        #pragma omp parallel
        {
            boost::dynamic_bitset<>flags{dataNum,0};
            #pragma omp for schedule(dynamic,100)
            for(int i=0;i<prev.size();i++){
                flags.reset();
                Node* temp=prev[i];
                priority_queue<std::pair<dist_t,int> >candidates;
                select_candidates(temp,temp->parent,candidates,S1,flags);
                vector<int>result;
                while(!candidates.empty()){
                    result.push_back(candidates.top().second);
                    candidates.pop();
                }
                for(int j=result.size()-1;j>=0;j--) add_neigh(Did(temp),result[j],Neigh_dist,min_cover_);
                //for(int j=result.size()-1;j>=0;j--) add_neigh(result[j],Did(temp),Neigh_dist,0.5);
                /*while(candidates.size()>Max_neigh) candidates.pop();
                int* neighbor=neighbors[Did(temp)];
                neighbor[0]=candidates.size();
                while(candidates.size()){
                    neighbor[candidates.size()]=candidates.top().second;
                    Neigh_dist[Did(temp)*Max_neigh+candidates.size()-1]=candidates.top().first;
                    candidates.pop();
                }
            }
        }
        next.clear();
        for(int i=0;i<prev.size();i++){
            for(int j=0;j<prev[i]->children.size();j++) next.push_back(prev[i]->children[j]);
        }
        prev=next;
    }*/
    
    
    int num,nndim,**nngraph;
    Load(nnfile.c_str(),nngraph,num,nndim);
    radius= new float[dataNum];
    for(int i=0;i<dataNum;i++) radius[i]=sqrt(Dist(i,nngraph[i][0]));
    vector<vector<int> >bknng(dataNum,vector<int>(0));
    for(int i=0;i<dataNum;i++){
        for(int j=0;j<S;j++){
            bknng[i].push_back(nngraph[i][j]);
            bknng[nngraph[i][j]].push_back(i);
        }
        for(int j=0;j<nodelist[i]->children.size();j++) bknng[i].push_back(Did(nodelist[i]->children[j]));
    }

    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,100)
        for(int i=0;i<dataNum;i++){  
            NeighborSelection(i,bknng[i],Neigh_dist,min_cover_);
        }
    }
    
    
    delete[] Neigh_dist;
    for(int i=0;i<num;i++) delete[]nngraph[i];
    delete[]nngraph;
}
void TBSG::CoverTree_test(){
    stack<Node*>nstack;
    nstack.push(covertree->root);
    int ncnt=covertree->N;
    cout<<"ncnt: "<<ncnt<<endl;
    while(!nstack.empty()){
        Node* temp=nstack.top();
        nstack.pop();
        for(int i=0;i<temp->children.size();i++) nstack.push(temp->children[i]);
        ncnt--;
    }
    if(ncnt>0) cout<<"error ncnt\n";
}
void TBSG::CoverTree_save(string filename){
    ofstream fout(filename.c_str(),std::ios::binary);
    stack<Node*>nstack;
    nstack.push(covertree->root);
    int ncnt=covertree->N;
    fout.write((char*)&ncnt,sizeof(int));
    while(!nstack.empty()){
        ncnt--;
        Node* temp=nstack.top();
        nstack.pop();
        for(int i=0;i<temp->children.size();i++) nstack.push(temp->children[i]);
        fout.write((char*)&(temp->ID),sizeof(unsigned));
        fout.write((char*)&(temp->level),sizeof(int));
        int size0=temp->children.size();
        fout.write((char*)&(size0),sizeof(int));
        for(int i=0;i<temp->children.size();i++) fout.write((char*)&(temp->children[i]->ID),sizeof(unsigned));
        size0=temp->dataId.size();
        fout.write((char*)&(size0),sizeof(int));
        for(int i=0;i<temp->dataId.size();i++) fout.write((char*)&(temp->dataId[i]),sizeof(int));
        if(temp->parent!=NULL) fout.write((char*)&(temp->parent->ID),sizeof(unsigned));
        else {
            int parentid=-1;
            fout.write((char*)&parentid,sizeof(int));
        }
        fout.write((char*)&(temp->maxdistUB),sizeof(double));
    }
    if(ncnt!=0) cout<<"covertree node num error\n";
    fout.close();
}
CoverTree* TBSG::CoverTree_load(string filename,vector<float*> plist){
    ifstream fin(filename.c_str(),std::ios::binary);
    CoverTree* covertree=new CoverTree(-1);
    int ncnt=0;
    fin.read((char*)&ncnt,sizeof(int));
    covertree->N=ncnt;
    covertree->D=dim;
    vector<Node*>nlist(ncnt,NULL);
    vector<vector<int> >children_ids(ncnt,vector<int>(0));
    vector<int>par_id(ncnt,0);
    for(int i=0;i<ncnt;i++){
        Node* temp=new Node;
        if(i==0) covertree->root=temp;
        fin.read((char*)&(temp->ID),sizeof(unsigned));
        nlist[temp->ID]=temp;
        fin.read((char*)&(temp->level),sizeof(int));
        int childnum;
        fin.read((char*)&(childnum),sizeof(int));
        int children[childnum];
        fin.read((char*)(children),sizeof(int)*childnum);
        for(int j=0;j<childnum;j++) children_ids[temp->ID].push_back(children[j]);
        int datanum;
        fin.read((char*)&datanum,sizeof(int));
        int dataid[datanum];
        fin.read((char*)dataid,sizeof(int)*datanum);
        for(int j=0;j<datanum;j++) temp->dataId.push_back(dataid[j]);
        temp->data=plist[dataid[0]];
        int parentid;
        fin.read((char*)&(parentid),sizeof(unsigned));
        par_id[temp->ID]=parentid;
        fin.read((char*)&(temp->maxdistUB),sizeof(double));
    }
    fin.close();
    for(int i=0;i<nlist.size();i++){
        Node* temp=nlist[i];
        if(temp==NULL) continue;
        for(int j=0;j<children_ids[temp->ID].size();j++){
            int id=children_ids[temp->ID][j];
            temp->children.push_back(nlist[id]);
        }
        if(par_id[temp->ID]>=0) temp->parent=nlist[par_id[temp->ID]];
        else temp->parent=NULL;
    }
    return covertree;
}
void TBSG::TBSG_save(string filename){
    ofstream fout(filename.c_str(),std::ios::binary);
    fout.write((char*)&dataNum,sizeof(int));
    fout.write((char*)&ep_,sizeof(int));
    fout.write((char*)&Max_neigh,sizeof(int));
    for(int i=0;i<dataNum;i++){
        int nsize=neighbors[i][0];
        fout.write((char*)neighbors[i],sizeof(int)*(nsize+1));
    }
    fout.close();
}
void TBSG::TBSG_load(string filename){
    if(filename[filename.size()-1]=='2'){
        ifstream fin(filename.c_str(),std::ios::binary);
        int eps_size;
        fin.read((char*)&eps_size,4);
        eps.resize(eps_size,0);
        for(int i=0;i<eps_size;i++) fin.read((char*)&eps[i],4);
        ep_=eps[0];
        fin.read((char*)&dataNum,4); 
        fin.read((char*)&Max_neigh,4);
        cout<<dataNum<<' '<<Max_neigh<<endl;
        neighbors=new int*[dataNum];
        for(int i=0;i<dataNum;i++) neighbors[i]=new int[Max_neigh+1];
        for(int i=0;i<dataNum;i++){
            fin.read((char*)neighbors[i],4);
            fin.read((char*)(neighbors[i]+1),4*neighbors[i][0]);
        }
        return;
    }
    ifstream fin(filename.c_str(),std::ios::binary);
    fin.read((char*)&dataNum,sizeof(int));
    fin.read((char*)&ep_,sizeof(int));
    fin.read((char*)&Max_neigh,sizeof(int));
    neighbors=new int*[dataNum];
    for(int i=0;i<dataNum;i++) neighbors[i]=new int[Max_neigh+1];
    for(int i=0;i<dataNum;i++){
        fin.read((char*)neighbors[i],sizeof(int));
        fin.read((char*)(neighbors[i]+1),sizeof(int)*neighbors[i][0]);
    }
    fin.close();
}
TBSG::TBSG(vector<float*> pList,float* data,int m,int Dim,int s,float min_c,string nnfile1){
    dataNum=pList.size();
    dim=Dim;
    data_=data;
    l2space =new hnswlib::L2Space(dim);
    fstdistfunc_ =l2space->get_dist_func();
    dist_func_param_ = l2space->get_dist_func_param();
    cout<<"CoverTree Construction starts\n";
    auto start = std::chrono::high_resolution_clock::now();
    covertree=new CoverTree(pList,dim,-1);
    cout<<"root->level "<<covertree->root->level<<endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    printf("CoverTree Construction finishes, costing %f seconds\n",(float)diff.count());
    nodelist.resize(dataNum,NULL);
    covertree_dfs(covertree->root);
    ep_=Did(covertree->root);
    Max_neigh=m;
    S=s;
    min_cover_=min_c;
    nnfile=nnfile1;
    neighbors=new int*[dataNum];
    for(int i=0;i<dataNum;i++) neighbors[i]=new int[Max_neigh+1];
    for(int i=0;i<dataNum;i++) neighbors[i][0]=0;
    start = std::chrono::high_resolution_clock::now();
    printf("TBSG Construction starts\n");
    cal_cnt=0;
    TBSG_from_covertree();
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    printf("TBSG Construction finishes, costing %f seconds\n",(float)diff.count());
    int avg=0;
    for(int i=0;i<dataNum;i++) avg+=neighbors[i][0];
    cout<<"avg neigh = "<<avg/dataNum<<endl;
}
TBSG::TBSG(float* data,int Dim,string tbsg_file){
    data_=data;
    dim=Dim;
    l2space =new hnswlib::L2Space(dim);
    fstdistfunc_ =l2space->get_dist_func();
    dist_func_param_ = l2space->get_dist_func_param();
    TBSG_load(tbsg_file);
    int avg=0;
    for(int i=0;i<dataNum;i++) avg+=neighbors[i][0];
    cout<<"avg neigh = "<<avg/dataNum<<endl;
}
void TBSG::free_mem(){
    for(int i=0;i<dataNum;i++) delete []neighbors[i];
    delete[]neighbors;
}
void TBSG::SearchInGraph(bool internal,float* p,priority_queue<std::pair<dist_t,int> >&top_candidate,boost::dynamic_bitset<> &flags,int fss_){
    int curr_id=ep_;
    dist_t curr_dist=Dist(curr_id,p);
    dist_t dist_bound=std::numeric_limits<dist_t>::max();
    priority_queue<std::pair<dist_t,int>,vector<std::pair<dist_t,int> >,greater<std::pair<dist_t,int> > >candidate_set;
    candidate_set.emplace(curr_dist,curr_id);
    flags[curr_id]=true;
    while(!top_candidate.empty()) top_candidate.pop();
    top_candidate.emplace(curr_dist,curr_id);
    if(eps.size()){
        int L = max(int(eps.size()),fss_);
        std::vector<unsigned> init_ids(L);
       // initializer_->Search(query, nullptr, L, parameter, init_ids.data());
       std::mt19937 rng(rand());
       GenRandom(rng, init_ids.data(), L , (unsigned)dataNum);
       for(int i=0;i<eps.size();i++) init_ids[i]=eps[i];
        for(int i=1;i<init_ids.size();i++){
            dist_t dist=Dist(init_ids[i],p);
            candidate_set.emplace(dist,init_ids[i]);
            top_candidate.emplace(dist,init_ids[i]);
            flags[init_ids[i]]=true;
        }
    }
    while(!candidate_set.empty()){
        curr_dist=candidate_set.top().first;
        curr_id=candidate_set.top().second;
        candidate_set.pop();
        if(curr_dist>dist_bound) break;
        if(internal) nodelist[curr_id]->mut.lock();
        int* neighbor=neighbors[curr_id];
        int nsize=neighbor[0];
        int next[nsize+1];
        for(int i=0;i<nsize;i++) next[i]=neighbor[i+1];
        if(internal) nodelist[curr_id]->mut.unlock();
        next[nsize]=0;
#ifdef USE_SSE
        //_mm_prefetch(data_+next[0]*dim, _MM_HINT_T0);
#endif
        for(int i=0;i<nsize;i++){
            int next_id=next[i];
            if(flags[next_id]) continue;
            flags[next_id]=true;
#ifdef USE_SSE
        //_mm_prefetch(data_+next[i+1]*dim, _MM_HINT_T0);
#endif
            dist_t next_dist=Dist(next_id,p);
            if(next_dist>dist_bound) continue;
            candidate_set.emplace(next_dist,next_id);
            top_candidate.emplace(next_dist,next_id);
            while(top_candidate.size()>fss_){
                dist_bound=top_candidate.top().first;
                top_candidate.pop();
            }
        }
    }
}
void TBSG::SearchKnn(vector<float* >qlist,int k,vector<priority_queue<std::pair<dist_t,int> > >&result,int fss_){
    std::priority_queue<std::pair<dist_t,int> > que;
    result.clear();
    cal_cnt=0;
    boost::dynamic_bitset<> flags{dataNum, 0};
    for(int i=0;i<qlist.size();i++){
        flags.reset();
        SearchInGraph(false,qlist[i],que,flags,fss_);
        while(que.size()>k) que.pop();
        result.push_back(que);
    }
    cout<<cal_cnt/qlist.size()<<' '<<dataNum/(cal_cnt/qlist.size())<<endl;
}
float precision_calc(vector<vector<int> >gt,vector<std::priority_queue<std::pair<dist_t,int> > >result,int k){
    if(gt.size()!=result.size()){
        cout<<"gt.size()!=result.size()\n";
        return 0;
    }
    int total=0;
    for(int i=0;i<gt.size();i++) total+=k;
    int correct=0;
    set<int>ret_set;
    for(int i=0;i<gt.size();i++){
        ret_set.clear();
        while(result[i].size()>k) result[i].pop();
        while(!result[i].empty()){
            ret_set.insert(result[i].top().second);
            result[i].pop();
        }
        for(int j=0;j<k;j++) if(ret_set.find(gt[i][j])!=ret_set.end()) correct++;
    }
    return 1.0*correct/total;
}
void test(vector<vector<int> >gt,vector<std::priority_queue<std::pair<dist_t,int> > >result,int k){
    set<int>ret_set,gt_set;
    int cnt0=0,cnt1=0,cnt2=0,cnt3=0;
    int cnt4=0;
    for(int i=0;i<gt.size();i++){
        ret_set.clear();
        gt_set.clear();
        for(int j=0;j<k;j++) gt_set.insert(gt[i][j]);
        while(result[i].size()>k) result[i].pop();
        while(!result[i].empty()){
            ret_set.insert(result[i].top().second);
            result[i].pop();
        }
        cnt4+=ret_set.size();
        for(set<int>::iterator ite=gt_set.begin();ite!=gt_set.end();ite++){
            if(ret_set.find(*ite)==ret_set.end()) cnt0++;
            else cnt1++;
        }
        for(set<int>::iterator ite=ret_set.begin();ite!=ret_set.end();ite++){
            if(gt_set.find(*ite)==gt_set.end()) cnt2++;
            else cnt3++;
        }
    }
    cout<<cnt0<<' '<<cnt1<<' '<<cnt2<<' '<<cnt3<<' '<<cnt4<<endl;
}
void search_result_save(string filename,vector<priority_queue<std::pair<dist_t,int> > >result){
    ofstream fout(filename);
    char buff[100];
    for(int i=0;i<result.size();i++){
        fout<<result[i].size()<<' ';
        while(!result[i].empty()){
            std::pair<dist_t,int> temp=result[i].top();
            result[i].pop();
            sprintf(buff,"%d %.2lf, ",temp.second,temp.first);
            fout<<buff;
        }
        fout<<endl;
    }
}
void gt_result_save(string filename,vector<vector<int> >gt,vector<float*>massQ,vector<float*>mass,int dim){
    ofstream fout(filename);
    char buff[100];
    for(int i=0;i<gt.size();i++){
        fout<<i<<endl;
        for(int j=0;j<gt[i].size();j++){
           sprintf(buff,"%d %.2f, ",gt[i][j],l2dist(mass[gt[i][j]],massQ[i],dim));
           fout<<buff;
        }
        fout<<endl;
    }
}
void record_output(string filename,vector<pair<float,float> >record){
    ofstream fout(filename);
    fout<<record.size()<<endl;
    for(int i=0;i<record.size();i++){
        fout<<record[i].first<<' '<<record[i].second<<endl;
    }
    fout.close();
}
