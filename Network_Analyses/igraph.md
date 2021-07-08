
Use igraph on the adjacency matrix created by WGCNA


For MSI
```bash
tmux new -s wgcna
srun -N 1 -n 1 -c 1 --mem=150gb -t 12:00:00 -p interactive --pty bash
module load R/3.6.0
cd /scratch.global/edittmar/WGCNA
R

```

For sapelo
```bash
tmux new -s igraph
cd /scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses_Inbred
#srun --pty -p batch --mem=100G --nodes=1 --ntasks-per-node=1  --constraint=EPYC --time=12:00:00 --job-name=qlogin /bin/bash -l
srun --pty  -p inter_p  --mem=150G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l # crashed at 100
#srun --pty  -p highmem_q  --mem=250G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l # 
module load R/4.0.0-foss-2019b
R
```

```R
library(WGCNA)
library(igraph)

load("Consensus_Network/ConsensusTOM_signed.RData")
# consensusTOM
length(consensusTOM[,1]) # 34408 x 34408 matrix

load("multiExpr.RData")
checkSets(multiExpr)

length(colnames(multiExpr[[1]]$data)) # gene identities

# assign to consensus TOM:
rownames(consensusTOM) <- colnames(multiExpr[[1]]$data)
colnames(consensusTOM) <- colnames(multiExpr[[1]]$data)

# diagonal values?
consensusTOM[1,1] #1
consensusTOM[2,2] # 1
consensusTOM[2,1:4]

graph <- graph_from_adjacency_matrix(
consensusTOM,
mode = "undirected",
weighted = TRUE,
diag = FALSE,
add.colnames = NULL,
add.rownames = NA
)
format(object.size(graph), units = "auto") #22.1 Gb
save(graph, file="igraph_object.RData")


graph_simp <- simplify(graph, remove.multiple = F, remove.loops = T) # removes loops in the graph

##############

### attributes:
E(graph)       # The edges of the "graph" object

V(graph)       # The vertices of the "graph" object

E(graph)$type  # Edge attribute "type"

V(graph)$media # Vertex attribute "media"

####

### plotting:
plot(net, edge.arrow.size=.4,vertex.label=NA)
net <- simplify(net, remove.multiple = F, remove.loops = T) # removes loops in the graph

```

