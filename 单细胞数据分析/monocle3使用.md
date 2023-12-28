

## monocle包

- 在许多单细胞研究中，单个细胞以不同步的方式执行基因表达程序。实际上，每个细胞都是所研究的转录程序的快照。**Monocle包 提供了用于分析单细胞表达的工具**。
- Monocle 引入了在伪时间中对单个细胞进行排序的策略，通过利用单个细胞的异步进程，将它们沿着与生物过程（例如细胞分化）相对应的轨迹放置。
- **Monocle 通过利用先进的机器学习技术（反向图嵌入）从单细胞基因组数据中学习明确的主图来对细胞进行排序**，从而稳健而准确地解决复杂的生物过程。
- Monocle 还执行聚类（即使用 t-SNE 和密度峰值聚类）。然后，Monocle 进行差异基因表达测试，从而识别在不同状态、生物过程以及替代细胞命运之间差异表达的基因。

## monocle3介绍

monocle3是一个用于单细胞RNA测序数据伪时序建模和分析的R包。它是monocle软件包的最新版本，2019年再nature上发表研究 **哺乳动物器官发生的单细胞转录景观** （之前常用的是monocle2）

![](https://files.mdnice.com/user/23696/7ce53eaf-37fd-4d1f-b43f-2f95270bdef2.png)

### monocle3 用途

- `拓扑学建模`：monocle3基于单细胞数据的拓扑结构，将细胞状态建模成一个有向无环图（DAG），其中节点代表细胞状态，边代表状态之间的转换关系。
- `伪时序建模`：monocle3可以根据单细胞RNA测序数据构建伪时序，**即确定细胞在发育、分化等动态过程中的顺序**。
- `动态系统模型`：monocle3使用动态系统的原理来描述和模拟细胞状态之间的转变。用于研究细胞在发育过程中的动态变化，并识别关键的分支点和状态转换。
- `分支分析`：monocle3可以识别和分析细胞发育树中的分支点，帮助理解细胞发育过程中的重要节点和分支。
- `差异表达分析`：monocle3可以进行单细胞差异表达分析，**用于识别在不同状态或时间点下哪些基因的表达发生了显著变化**。
- `可视化工具`：monocle3 提供了丰富的可视化工具，用于可视化单细胞数据的伪时序、分支结构、差异表达等信息。

### monocle3更新了什么？
- 支持UMAP算法初始化轨迹推断。
- 支持具有多个根的轨迹。
- 学习具有环路或收敛点的轨迹的方法。
- 使用“近似图抽象”的思想自动划分单元以学习不相交或并行轨迹的算法。
- 对具有轨迹依赖性表达的基因进行新的统计测试。**这取代了旧 differentialGeneTest()函数 和BEAM()**.
- 保存和加载 Monocle 对象和转换模型。
- fit_models 的混合负二项分布。
- 用于可视化轨迹和基因表达的 3D 界面。



![monocle2](https://files.mdnice.com/user/23696/049973c4-9eaf-4db8-ad93-227c245eddb0.png)


![monocle3](https://files.mdnice.com/user/23696/6d39013f-b019-4a25-972d-8b3229f14a04.png)


### monocle3的原理

- `拓扑数据结构`： Monocle使用拓扑学的概念构建了一个表示细胞状态之间关系的拓扑数据结构。这被称为可达性图（reachable graph）。在这个图中，细胞是节点，边表示细胞之间的可达性，即它们在基因表达模式上相似。
- `状态定义`： 在可达性图中，每个细胞都被分配到不同的状态。这些状态代表了细胞在发育轨迹上的位置或其他生物学过程中的特定状态。Monocle使用拓扑学的概念来定义这些状态，并识别在细胞群中表达相似基因的亚群。
- `基于流形的降维`： Monocle使用非线性降维技术t-SNE和UMAP，将高维的基因表达数据映射到低维的空间。这有助于可视化细胞状态之间的关系，尤其是在发育轨迹上。
- `动态时间规整`（DDRTree）： Monocle中的DDRTree方法是一种基于流形学习的技术，用于在发育过程中识别动态的时间点。这有助于在时间序列数据中揭示细胞的动态变化。
- `基因表达的演变模型`： Monocle建立了基因表达的演变模型，该模型描述了基因在细胞状态转变过程中的表达变化。


## monocle3的使用

monocle3的使用和分析流程，**前面部分处理可以使用seurat包实现。**

![](https://files.mdnice.com/user/23696/b00b8f66-a9c3-49cb-8f9c-2c6a28ed2f0d.png)



#### 包的安装和加载
```r
# monocle3的安装
# devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)
```

#### 数据下载
**这里给一个示例数据（标准化之后的单细胞表达矩阵）**，有需要的可以下载一下。
百度云链接：https://pan.baidu.com/s/17woJ6RVdlReNjYH7Zw8h0w 提取码：ke2p

#### 先进行基础的Seurat流程
由于原数据进行过标准化，因此这里就不再标准化了。**关于Seurat的基础操作可以看之前的推文**:
> [单细胞分析基础（一）| 详解数据类型和预处理](https://mp.weixin.qq.com/s?__biz=Mzg2NjYzNjQ4Ng==&mid=2247486603&idx=1&sn=0dfa6227358de081ccb45717987cf723&chksm=ce468b22f9310234fe5d4ecd478e69bb8f3cff0a64da151c080c1d34ae69d94b2de47c78c5d8&token=1105195915&lang=zh_CN#rd)              
> [单细胞分析基础（二）| 数据降维和细胞类型注释](https://mp.weixin.qq.com/s?__biz=Mzg2NjYzNjQ4Ng==&mid=2247486644&idx=1&sn=7be5e4d55eeea6e65ff1f2e8e3f9f457&chksm=ce468b1df931020b82e7aa5c55569b75c10833a4e0e1765a35b546e8903ee8ff75249d51387e&token=1105195915&lang=zh_CN#rd)

```r
seu_kidney <- CreateSeuratObject(sparse_mat)
seu_kidney <- ScaleData(seu_kidney)
seu_kidney <- FindVariableFeatures(seu_kidney, selection.method = "vst", nfeatures = 1000)
seu_kidney <- RunPCA(seu_kidney, features = VariableFeatures(object = seu_kidney))
ElbowPlot(seu_kidney)
seu_kidney <- FindNeighbors(seu_kidney, dims = 1:10)
seu_kidney <- FindClusters(seu_kidney, resolution = 0.05)
seu_kidney <- RunUMAP(seu_kidney, dims = 1:10)
DimPlot(seu_kidney, reduction = "umap")
```
##### 这里选择1:10作为非线性降维的输入
![](https://files.mdnice.com/user/23696/6dfbd50c-8522-47df-8c6b-3102c2e1cc83.png)

##### 细胞被分为四类
![](https://files.mdnice.com/user/23696/63832cc3-514e-4f85-8f77-0929fc20af0f.png)

#### 构建cell_data_set类
cell_data_set 类是 monocle3 包中用于存储单细胞数据的主要数据结构之一。它包含了单细胞RNA测序数据以及与之相关的元数据信息。

```r
data <- GetAssayData(seu_kidney,assay = "RNA",slot = "counts")
cell_metadata <- seu_kidney@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)
```
![](https://files.mdnice.com/user/23696/710dc253-0027-4249-a8be-8752c6a36d2d.png)

- `rowRanges`：包含了关于基因的信息，例如基因的名称、ID等。
- `colData`：包含了每个单细胞样本的元数据信息，例如细胞类型、时间点等。
- `assays`：包含了单细胞数据的表达矩阵，可以包括原始的、标准化的、降维后的数据等。
- `reduce_dim_aux`：包含了降维后的数据，如PCA、t-SNE等。
- `principal_graph`：包含了单细胞数据的图结构，通常用于构建细胞状态的拓扑关系。
- `......`

#### 降维及可视化
这里用了seurat聚类的标签，也可以不使用，或者根据自己注释的真实细胞类型来标注。
```r
cds <- preprocess_cds(cds, num_dim = 5)
cds <- reduce_dimension(cds)
plot_cells(cds,color_cells_by = "seurat_clusters")
```

![](https://files.mdnice.com/user/23696/5a548c9e-2093-4734-918d-88353e8c210a.png)

##### 可以使用plot_cells()函数可视化单个基因的变化
```r
plot_cells(cds, genes=c("PTN", "GPM6B", "MT3", "CLU"))
plot_cells(cds, genes=c("KLF6", "DUSP1", "SGK1", "RGS2"))
```

![](https://files.mdnice.com/user/23696/ea8dfff4-d1cf-4045-a4d4-457b54f63632.png)

![](https://files.mdnice.com/user/23696/02ad5afa-0fef-4475-a79d-0a700909ec2c.png)

#### 将细胞聚集分区
Monocle能够通过其聚类过程了解何时应将细胞放置在同一轨迹中，而不是放置在单独的轨迹中。运行cluster_cells()，每个细胞不仅分配给一个簇，还分配给一个分区。
```r
cds <- cluster_cells(cds,cluster_method = "louvain")
plot_cells(cds, color_cells_by = "partition")
```
- 用 "louvain" 方法对单细胞数据进行聚类。"Louvain" 是一种基于图的聚类算法，通常用于发现社区结构，即在数据中寻找具有相似特征的细胞群。

![](https://files.mdnice.com/user/23696/40aec29c-43cd-457d-8f4c-6c9ef3d28ef1.png)

#### 轨迹图
Monocle 并不假设数据集中的所有细胞都源自共同的转录“祖先”。在许多实验中，实际上可能存在多个不同的轨迹。
```r
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition",label_groups_by_cluster=FALSE, label_leaves=FALSE,label_branch_points=FALSE)
```
**在两个不同分区分别绘制了轨迹图路径。**

![](https://files.mdnice.com/user/23696/d8fad18e-b208-4863-95be-8b945755b006.png)

#### 以伪时间对细胞进行排序
**为了将细胞按顺序排列，我们需要告诉Monocle生物过程的“起点”在哪里**。我们通过选择图表中标记为轨迹“根”的区域来实现这一点。通常**可以通过在UMAP空间中查找早期时间点的细胞占据的点来实现**：

```r
colData(cds)["partition"] <- cds@clusters$UMAP$partitions
plot_cells(cds, color_cells_by = "partition",label_cell_groups=FALSE,label_leaves=TRUE,label_branch_points=TRUE,graph_label_size=1.5)
```
**下图就是所有可能成为起点的位置**，接下来可以选择手动决定起点或者通过函数进行计算获得。

![](https://files.mdnice.com/user/23696/e6c08084-8104-4bf6-b0dd-e87d2bf6d98b.png)

#### 手动选择起点
现在已经了解了早期细胞落在哪里，我们可以使用order_cells()函数计算每个细胞在伪时间内落在哪里。**order_cells()需要指定根节点的轨迹图**。如果不提供它们作为参数，**它将启动一个图形用户界面来选择一个或多个根节点**。

```r
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE, graph_label_size=1.5)
```

![](https://files.mdnice.com/user/23696/80314100-f69b-43d2-9cd2-6a171b46462a.png)

#### 计算选择起点
下面的函数首先根据单元最接近的轨迹图节点对单元进行分组。然后，它计算每个节点上的单元有多少部分来自最早的时间点。
```r

get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == 2)
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds, color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,
           label_branch_points=FALSE,graph_label_size=1.5)
```

![](https://files.mdnice.com/user/23696/6e7cc3c3-5f67-4343-af11-3a3e80eee39d.png)

#### 寻找随伪时间变化的基因
**识别随着细胞沿着轨迹前进而变化的基因是此类分析的核心目标**。了解基因运行和关闭的顺序可以为新的发育模型提供信息。
```r
## 伪时序相关基因
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=8)
pr_deg_ids <- ciliated_cds_pr_test_res %>% filter(q_value < 0.05) %>% arrange(-morans_I)
plot_cells(cds, genes=c("CLDN11","WSB1","MACF1","MMP16"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
```

![](https://files.mdnice.com/user/23696/a2d70123-d613-4bfc-b2d9-d68b1570deac.png)

##### Monocle提供了另一种绘图功能，可以更清晰地了解基因沿单一路径的动态

```r
plot_genes_in_pseudotime(cds[c("CLDN11","WSB1","MMP16","MACF1"),],
                         color_cells_by="seurat_clusters",
                         min_expr=0.5)
```

![](https://files.mdnice.com/user/23696/9dd7bfab-ea8c-461b-8b4a-58891d800c95.png)









