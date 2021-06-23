TBSG: Tree-based Search Graph for Approximate Nearest Neighbor Search.
------
TBSG is a graph-based algorithm for ANNS based on Cover Tree, which is also an approximation of Monotonic Search Network (MSNET). TBSG is very efficient with high precision.

Benchmark datasets
------
Datasets    |   No. of base   |  dimension   |  No. of query  |    download link  
Sift        |   1,000,000      |   128        |    10,000     |    (http://corpus-texmex.irisa.fr/)  
Gist        |   1,000,000      |   300        |    1,000      |    (http://corpus-texmex.irisa.fr/)  
Glove       |   1,183,514      |   100        |    10,000     |    (http://downloads.zjulearning.org.cn/data/glove-100.tar.gz)  
Crawl       |   1,989,995      |   300        |    10,000     |    (http://commoncrawl.org/)  

### How to use TBSG
### 1) compile

* Prerequisite : openmp, cmake, eigen3
```bash
$ cd /path/to/project  
$ cmake . && make  
```
### 2) build an approximate kNNG

We use [efanna\_graph](https://github.com/ZJULearning/efanna\_graph) to build the kNNG.

### 3) create a TBSG index
  

```bash  
$ cd /path/to/project/  
$ ./TBSG_index data_path M S MP nnfile save_path  
```

**data_path** is the path of base data.  
**M** is the maximum of size of neighbors.  
**S** is the candidate set size to build TBSG.  
**MP** is the minimum of min_prob.  
**nnfile** is the file of k nearest neighbor graph.  
**save_path** is the path to save the index.  

### 4) search with TBSG index

```bash
$ cd /path/to/project/
$ ./TBSG_search data_path query_path groundtruth_path save_path step
```
**data_path** is the path of base data.  
**query_path** is the path of query data.  
**groundtruth** is the path of groundtruth data.  
**save_path** is the path to save the index.  
**step** is the step size to expand the search pool.  

### Parameters used for four datasets

### parameters for building kNNG

|  Dataset  |  K  |  L  | iter |  S |  R  |
|:---------:|:---:|:---:|:----:|:--:|:---:|
| Sift      | 200 | 200 |  12  | 10 | 100 |
| Gist      | 400 | 400 |  12  | 15 | 100 |
| Glove     | 400 | 420 |  12  | 20 | 300 |
| Crawl     | 400 | 420 |  12  | 20 | 100 |

### parameters for building index

|Datasets   |  M  |  S  |   MP   |
|:---------:|:---:|:---:|:------:|
|Sift       |  50 | 100 |  0.53  |
|Gist       |  70 | 200 |  0.515 |
|Glove      |  80 | 300 |  0.53  |
|Crawl      |  50 | 200 |  0.53  |


