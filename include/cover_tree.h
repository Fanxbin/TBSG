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

# ifndef _COVER_TREE_H
# define _COVER_TREE_H

//#define DEBUG
#include <atomic>
#include <fstream>
#include <iostream>
#include <stack>
#include <map>
#include <vector>
#include <shared_mutex>
#include <mutex>
#include <math.h>
#include <thread>

#define SHARED_MUTEX_TYPE shared_mutex
/*
#ifdef __clang__
#define SHARED_MUTEX_TYPE shared_mutex
#else
#define SHARED_MUTEX_TYPE shared_timed_mutex
#endif
*/

//typedef pointType::Scalar dtype;
float l2dist(float* a,float* b, int dim);
class CoverTree
{
/************************* Internal Functions ***********************************************/
protected:
    /*** Base to use for the calculations ***/
    static constexpr double base = 1.3;
    static double* compute_pow_table();
    static double* powdict;
    friend class TBSG;

public:
    /*** structure for each node ***/
    struct Node
    {
        std::vector<Node*> children;        // list of children
        int level;                          // current level of the node
        double maxdistUB;                   // upper bound of distance to any of descendants
        unsigned ID;                   // unique ID of current node
        std::vector<int>dataId;
        Node* parent;                       // parent of current node
        float* data;
        mutable std::SHARED_MUTEX_TYPE mut;// lock for current node

        /*** Node modifiers ***/
        double covdist()                    // covering distance of subtree at current node
        {
            return powdict[level + 1024];
        }
        double sepdist()                    // separating distance between nodes at current level
        {
            return powdict[level + 1023];
        }

        float dist(Node* n,int dim) const
        {
            return sqrt(l2dist(data,n->data,dim));
        }
        float dist(float* p,int dim) const
        {
            return sqrt(l2dist(data,p,dim));
        }
        Node* setChildp(float* p,int p_id, int new_id=-1)   // insert a new child of current node with point pIns
        {
            Node* temp = new Node;
            temp->level = level - 1;
            temp->maxdistUB = 0; // powdict[level + 1024];
            temp->ID = new_id;
            temp->dataId.push_back(p_id);//added by Fan
            temp->parent = this;
            temp->data=p;
            children.push_back(temp);
            return temp;
        }
        /*** erase child ***/
        void erase(size_t pos)
        {
            children[pos] = children.back();
            children.pop_back();
        }

        void erase(std::vector<Node*>::iterator pos)
        {
            *pos = children.back();
            children.pop_back();
        }

        /*** Iterator access ***/
        inline std::vector<Node*>::iterator begin()
        {
            return children.begin();
        }
        inline std::vector<Node*>::iterator end()
        {
            return children.end();
        }
        inline std::vector<Node*>::const_iterator begin() const
        {
            return children.begin();
        }
        inline std::vector<Node*>::const_iterator end() const
        {
            return children.end();
        }
    };
    // mutable std::map<int,std::atomic<unsigned>> dist_count;
    std::map<int,unsigned> level_count;

protected:
    Node* root;                         // Root of the tree
    std::atomic<int> min_scale;         // Minimum scale
    std::atomic<int> max_scale;         // Minimum scale
    //int min_scale;                    // Minimum scale
    //int max_scale;                    // Minimum scale
    int truncate_level;                 // Relative level below which the tree is truncated
    bool id_valid;

    std::atomic<unsigned> N;            // Number of points in the cover tree
    //unsigned N;                       // Number of points in the cover tree
    unsigned D;                         // Dimension of the points

    std::SHARED_MUTEX_TYPE global_mut;  // lock for changing the root

    bool insert(Node* current, float* p,int p_id);
    //bool insert(Node* current, Node* p,int p_id);


public:
    /*** Internal Contructors ***/
    /*** Constructor: needs at least 1 point to make a valid cover-tree ***/
    // NULL tree
    explicit CoverTree(int truncate = -1);
    // cover tree with one point as root
    CoverTree(float* p,int dimensions, int truncate = -1);
    // cover tree using points in the list between begin and end
    CoverTree(std::vector<float*> pList, int dimensions, int truncate = -1);
    /*** Destructor ***/
    /*** Destructor: deallocating all memories by a post order traversal ***/
    ~CoverTree();
    void free_memory();
/************************* Public API ***********************************************/
public:

    /*** Insert point p into the cover tree ***/
    bool insert(float* p,int p_id);

    int get_level();

};

#endif //_COVER_TREE_H

