/*
 * Copyright (c) 2017 Manzil Zaheer All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "cover_tree.h"
#include "utils.h"
#include <stack>
#include <numeric>
double* CoverTree::compute_pow_table()
{
    double* powdict = new double[2048];
    for (int i = 0; i<2048; ++i)
        powdict[i] = pow(CoverTree::base, i - 1024);
    return powdict;
}

double* CoverTree::powdict = compute_pow_table();

bool CoverTree::insert(CoverTree::Node* current, float* p,int p_id)
{
    bool result = false;
#ifdef DEBUG
    if (current->dist(p) > current->covdist())
        throw std::runtime_error("Internal insert got wrong input!");
    if (truncateLevel > 0 && current->level < maxScale - truncateLevel)
    {
        std::cout << maxScale;
        std::cout << " skipped" << std::endl;
        return false;
    }
#endif
    if (truncate_level > 0 && current->level < max_scale-truncate_level)
        return false;

    // acquire read lock
    current->mut.lock_shared();

    // Sort the children
    unsigned num_children = current->children.size();
    std::vector<int> idx(num_children);
    std::iota(std::begin(idx), std::end(idx), 0);
    std::vector<float> dists(num_children);
    for (unsigned i = 0; i < num_children; ++i)
        dists[i] = current->children[i]->dist(p, D);
    auto comp_x = [&dists](int a, int b) { return dists[a] < dists[b]; };
    std::sort(std::begin(idx), std::end(idx), comp_x);

    bool flag = true;
    for (const auto& child_idx : idx)
    {
        Node* child = current->children[child_idx];
        float dist_child = dists[child_idx];
        if (dist_child <= 0.0)
        {
            // release read lock then enter child
            current->mut.unlock_shared();
	        child->mut.lock();
	        child->dataId.push_back(p_id);
            child->mut.unlock();
            flag = false;
            result=true;
            //std::cout << "Duplicate entry!!!" << std::endl;
            break;
        }
        else if (dist_child <= child->covdist())
        {
            // release read lock then enter child
            if (child->maxdistUB < dist_child)
                child->maxdistUB = dist_child;
            current->mut.unlock_shared();
            result = insert(child, p , p_id);
            flag = false;
            break;
        }
    }

    if (flag)
    {
        // release read lock then acquire write lock
        current->mut.unlock_shared();
        current->mut.lock();
        // check if insert is still valid, i.e. no other point was inserted else restart
        if (num_children==current->children.size())
        {
            int new_id = N++;
            current->setChildp(p, p_id, new_id);//changed by Fan;
            result = true;
            current->mut.unlock();
            int local_min = min_scale.load();
            while( local_min > current->level - 1){
                min_scale.compare_exchange_weak(local_min, current->level - 1, std::memory_order_relaxed, std::memory_order_relaxed);
                local_min = min_scale.load();
            }
        }
        else
        {

            current->mut.unlock();
            result = insert(current, p, p_id);
        }
        //if (min_scale > current->level - 1)
        //{
            //min_scale = current->level - 1;
            ////std::cout << minScale << " " << maxScale << std::endl;
        //}
    }
    return result;
}

bool CoverTree::insert(float* p,int p_id){
    bool result = false;
    id_valid = false;
    global_mut.lock_shared();
    if (root->dist(p, D) > root->covdist())
    {
        global_mut.unlock_shared();
        std::cout<<"Entered case 1: " << root->dist(p, D) << " " << root->covdist() << " " << root->level <<std::endl;
        std::cout<<"Requesting global lock!" <<std::endl;
        global_mut.lock();
        while (root->dist(p, D) > base * root->covdist()/(base-1))
        {
            CoverTree::Node* current = root;
            CoverTree::Node* parent = NULL;
            while (current->children.size()>0)
            {
                parent = current;
                current = current->children.back();
            }
            if (parent != NULL)
            {
                parent->children.pop_back();
                current->level = root->level + 1;
                //current->maxdistUB = powdict[current->level + 1025];
                current->children.push_back(root);
                root = current;
            }
            else
            {
                root->level += 1;
                //root->maxdistUB = powdict[root->level + 1025];
            }
        }
        CoverTree::Node* temp = new CoverTree::Node;
        temp->level = root->level + 1;
        temp->parent = NULL;
        temp->dataId.push_back(p_id);//added by Fan
        temp->ID=N++;
        temp->data=p;
        //temp->maxdistUB = powdict[temp->level+1025];
        temp->children.push_back(root);
        root->parent = temp;
        root = temp;
        max_scale = root->level;
        result = true;
        //std::cout << "Upward: " << minScale << " " << maxScale << std::endl;
        global_mut.unlock();
        global_mut.lock_shared();
    }
    else
    {
        //root->tempDist = root->dist(p);
        result = insert(root, p, p_id);
    }
    global_mut.unlock_shared();
    return result;
}


/****************************** Internal Constructors of Cover Trees *************************************/

// constructor: NULL tree
CoverTree::CoverTree(const int truncate /*=-1*/ )
    : root(NULL)
    , min_scale(1000)
    , max_scale(0)
    , truncate_level(truncate)
    , id_valid(false)
    , N(0)
    , D(0)
{
}

// constructor: needs at least 1 point to make a valid cover-tree
CoverTree::CoverTree(float* p,int dimensions, int truncate /*=-1*/)
    : min_scale(1000)
    , max_scale(0)
    , truncate_level(truncate)
    , id_valid(false)
    , N(1)
    , D(dimensions)
{
    root = new CoverTree::Node;
    root->data=p;
    root->level = 0;
    root->maxdistUB = 0;
}


void CoverTree::free_memory(){
    std::stack<CoverTree::Node*>nstack;
    nstack.push(root);
    while(!nstack.empty()){
        CoverTree::Node* temp=nstack.top();
        nstack.pop();
        for(int i=0;i<temp->children.size();i++) nstack.push(temp->children[i]);
        delete temp;
    }
}

// constructor: cover tree using points in the list between begin and end
CoverTree::CoverTree(std::vector<float*>pList, int dimensions, int truncateArg /*= 0*/)
{
    int nd=pList.size();
    int dim=dimensions;
    float *mx=new float[dim];
    for(int i=0;i<dim;i++) mx[i]=0;
    for(int i=0;i<nd;i++){
        for(int j=0;j<dim;j++) mx[j]+=pList[i][j];
    }
    for(int i=0;i<dim;i++) mx[i]/=nd;

    float *dists=new float[nd];

    for(int i=0;i<nd;i++) dists[i]=sqrt(l2dist(pList[i], mx, dim));
    //3. argort the distance to find approximate mediod
    std::vector<int> idx(nd);
    std::iota(std::begin(idx), std::end(idx), 0);
    auto comp_x = [&dists](int a, int b) { return dists[a] > dists[b]; };
    std::sort(std::begin(idx), std::end(idx), comp_x);
    std::cout<<"Max distance: " << dists[idx[0]] << std::endl;

    //4. Compute distance of every point from the mediod
    delete [] mx;
    mx = pList[idx[0]];
    for(int i=0;i<nd;i++) dists[i]=sqrt(l2dist(pList[i],mx,dim));
    float max_dist=0;
    for(int i=0;i<nd;i++) max_dist=std::max(max_dist,dists[i]);

    delete []dists;
    int scale_val = std::ceil(std::log(max_dist)/std::log(base));
    std::cout<<"Scale chosen: " << scale_val << std::endl;
    float* temp = pList[idx[0]];
    min_scale = scale_val; //-1000;
    max_scale = scale_val; //-1000;
    truncate_level = truncateArg;
    N = 1;
    D = dim;

    root = new CoverTree::Node;
    root->level = scale_val; //-1000;
    root->maxdistUB = powdict[scale_val+1024];
    root->dataId.push_back(idx[0]);
    root->ID=0;
    root->data=temp;
    root->parent=NULL;

    int run_till = 50000 < nd ? 50000 : nd;
    for (int i = 1; i < run_till; ++i){
        utils::progressbar(i, run_till);
        if(!insert(pList[idx[i]],idx[i]))
            std::cout << "Insert failed!!!" << std::endl;
    }
    utils::progressbar(run_till, run_till);
    std::cout<<std::endl;

    std::cout << dimensions<< ", " << nd << std::endl;

    utils::parallel_for_progressbar(50000,nd,[&](int i)->void{
    //for (int i = 50000; i < end; ++i){
        //utils::progressbar(i, end-50000);
        if(!insert(pList[idx[i]],idx[i]))
            std::cout << "Insert failed!!!" << std::endl;
    });
}

// constructor: cover tree using points in the list between begin and end
// CoverTree::CoverTree(Eigen::Map<Eigen::MatrixXd>& pMatrix, int begin, int end, int truncateArg /*= 0*/)
// {
    // pointType temp = pMatrix.col(begin);

    // min_scale = 7;
    // max_scale = 7;
    // truncate_level = truncateArg;
    // N = 1;
    // D = temp.rows();

    // root = new CoverTree::Node;
    // root->_p = temp;
    // root->level = 7;
    // root->maxdistUB = powdict[7+1024];

    // //utils::parallel_for(begin+1,end,[&](int i)->void{
    // for (int i = begin + 1; i < end; ++i){
        // //std::cout<<i<<std::endl;
        // if(i%1000==0)
            // std::cout << i << std::endl;
        // insert(pMatrix.col(i));
    // }//);
// }

// destructor: deallocating all memories by a post order traversal

/******************************************* Auxiliary Functions ***************************************************/

// get root level == max_level
int CoverTree::get_level()
{
    return root->level;
}


/******************************************* Functions to remove ***************************************************/


