#ifndef _TBSG_h_
#define _TBSG_h_
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
#include "space_l2.h"
#define Node CoverTree::Node
#define dist_t float
using namespace hnswlib;
using namespace std;
float precision_calc(vector<vector<int> >gt,vector<std::priority_queue<std::pair<dist_t,int> > >result,int k);
void record_output(string filename,vector<pair<float,float> >record);
class TBSG{
private:
    CoverTree* covertree;
    int** neighbors;
    vector<Node*>nodelist;
    int Max_neigh;
    int S;
    int dataNum;
    int dim;
    float *radius;
    float* data_;
    float min_cover_;
    int ep_;
    vector<int>eps;
    void covertree_dfs(Node* par);
    int Did(Node* self);
    vector<int> all_ids(Node* self);
    dist_t Dist(int fir,int sec);
    dist_t Dist(int id,float* p);
    float cover_prob(float a,float b);
    bool occludeByNeighbor(dist_t neigh_dist,dist_t curr_dist,dist_t cand_dist,float min_cover,float r);
    void NeighborSelection(int self,vector<int>candidate,dist_t* Neigh_dist,float min_c);
    string nnfile;
    void insert(int* nn,int num,int curr,float mc,dist_t* Neigh_dist);
    bool occludeByNeighbor(dist_t neigh_dist,dist_t curr_dist,dist_t cand_dist,float min_cover);
    void add_neigh(int self,int new_neigh,dist_t *Neigh_dist,float min_c);
    void select_candidates(Node* self,Node* par,std::priority_queue<std::pair<dist_t,int> >&candidates,int cand_size,boost::dynamic_bitset<> &flags);
    void Insert(int self,dist_t *Neigh_dist,boost::dynamic_bitset<> &flags);
    void TBSG_from_covertree();
    void CoverTree_test();
    void CoverTree_save(string filename);
    CoverTree* CoverTree_load(string filename,vector<float*> plist);
    void free_mem();
    void SearchInGraph(bool internal,float* p,priority_queue<std::pair<dist_t,int> >&top_candidate,boost::dynamic_bitset<> &flags,int fss_);
public:
    TBSG(vector<float*> pList,float* data,int m,int Dim,int s,float min_c,string nnfile1);
    TBSG(float* data,int Dim,string tbsg_file);
    void TBSG_save(string filename);
    void TBSG_load(string filename);
    void SearchKnn(vector<float* >qlist,int k,vector<priority_queue<std::pair<dist_t,int> > >&result,int fss_);
};
#endif
