if (!require("igraph")) {install.packages("igraph")}
if (!require("NetIndices")) {install.packages("NetIndices")}
library(reshape)
library(ggplot2)
setwd("/Users/trvinh/work/thesis/microsporidia/koAnnotation")

pathIDs <- c("00010", "00020")#, "00030", "00230", "00240", "90001")  # glycolysis, gluconeogenesis, the Krebs cycle, pentose phosphate pathway, purine and pyrimidine metabolism, and amino acid metabolism
# id <- "00010"
networkProperty <- data.frame()
for(id in pathIDs){
  file <- paste0("/Users/trvinh/work/R_projects/keggcxn/data/keggcxn/",id,".cxn")
  inputDf <- as.data.frame(read.table(file, sep='\t',header=F,check.names=FALSE,comment.char=""))
  colnames(inputDf) <- c("from","to","interaction")
  
  # get all nodes
  sourceDf <- data.frame(inputDf$from,stringsAsFactors = FALSE)
  colnames(sourceDf) <- "id"
  targetDf <- data.frame(inputDf$to,stringsAsFactors = FALSE)
  colnames(targetDf) <- "id"
  nodeDf <- unique(rbind(sourceDf,targetDf))
  
  # filter nodes based on annotated proteins
  annoDf <- as.data.frame(read.table("/Users/trvinh/work/R_projects/keggcxn/data/koMappedProteins.list", sep='\t',header=TRUE,check.names=FALSE,comment.char="",stringsAsFactors=FALSE))
  colnames(annoDf)[2] <- "id"
  
  # map annotated nodes to reference net
  joinedDf <- merge(nodeDf,annoDf, by="id", all.x=TRUE)
  joinedDf$group <- joinedDf$source
  
  # remove un-annotated nodes if requested
  joinedDf <- joinedDf[!is.na(joinedDf$geneID),]
  
  ### change group type for reference nodes
  # if(length(unique(complete.cases(joinedDf))) > 1){
  if(nrow(joinedDf[is.na(joinedDf$geneID),]) > 0){
    joinedDf[is.na(joinedDf$geneID),]$group <- "reference"
  }

  for(source in levels(as.factor(joinedDf$source))){
    # node data
    nodeDf <- data.frame("id" = as.character(unique(joinedDf[joinedDf$source == source,"id"])))
    
    ### edge data
    edgeDf <- inputDf[,(1:2)]
    edgeDf <- edgeDf[edgeDf$from %in% nodeDf$id & edgeDf$to %in% nodeDf$id,]

    # create igraph object
    graph <- graph_from_data_frame(edgeDf, directed = FALSE, vertices = nodeDf)
    
    # use GenInd function from NetIndices package to output network properties
    graph.adj<-get.adjacency(graph,sparse=FALSE)
    graph.properties<-GenInd(graph.adj)
    
    # get degree for all nodes
    all.deg.graph <- as.data.frame(degree(graph,v=V(graph),mode="all"))
    
    # network properties
    netProp <- data.frame(
      "Nodes" = graph.properties$N,  # number of nodes
      "Edges" = graph.properties$Ltot/2, # number of links
      "Avg_Degree" = graph.properties$LD, # same as mean(all.deg.graph[,1]); link density (average # of links per node)
      "Max_Degree" = max(all.deg.graph[,1]),
      "Avg_Path_Len" = average.path.length(graph, unconnected=TRUE),
      "Max_Path_Len_Diameter" = diameter(graph)
    )
    
    netProp$PathID <- id
    netProp$Source <- source
    networkProperty <- rbind(networkProperty,netProp)
  }
}

### output network properties
networkProperty <- networkProperty[,c("PathID","Source","Nodes","Edges","Avg_Degree","Max_Degree","Avg_Path_Len","Max_Path_Len_Diameter")]
networkProperty$PathID[networkProperty$PathID == "00010"] <- "Glycolysis / Gluconeogenesis"
networkProperty$PathID[networkProperty$PathID == "00020"] <- "TCA cycle"
networkProperty$PathID[networkProperty$PathID == "00030"] <- "Pentose phosphate pathway"
networkProperty$PathID[networkProperty$PathID == "00230"] <- "Purine metabolism"
networkProperty$PathID[networkProperty$PathID == "00240"] <- "Pyrimidine metabolism"
networkProperty$PathID[networkProperty$PathID == "90001"] <- "Amino acid metabolism"
colnames(networkProperty) <- c("Pathway","Source","Nodes","Edges","Avg_Degree","Max_Degree","Avg_Path_Len","Max_Path_Len_Diameter")

# head(networkProperty)
write.table(networkProperty,"networkProperty.txt",sep="\t",row.names = FALSE,quote = FALSE)

### properties statistics
meltedNetworkProp <- melt(networkProperty, id.vars = c("Pathway","Source"))
colnames(meltedNetworkProp) <- c("Pathway","Source","Property","Value")

# networkStat <- data.frame("source" = as.character(), "prop" = as.character(), "mean" = as.numeric(), "sd" = as.numeric(), stringsAsFactors = FALSE)
# i = 1
# head(networkProperty)
# 
# for(source in levels(as.factor(networkProperty$Source))){
#   for(prop in levels(as.factor(meltedNetworkProp$Property))){
#     sd <- sd(networkProperty[networkProperty$Source == source,][,prop])
#     mean <- mean(networkProperty[networkProperty$Source == source,][,prop])
#     networkStat[i,] <- c(source,prop,mean,sd)
#     i <- i+1
#   }
# }

meltedNetworkProp$Property <- as.character(meltedNetworkProp$Property)
meltedNetworkProp$Property[meltedNetworkProp$Property == "Avg_Degree"] <- "Avg. degree"
meltedNetworkProp$Property[meltedNetworkProp$Property == "Avg_Path_Len"] <- "Avg. path length"
meltedNetworkProp$Property[meltedNetworkProp$Property == "Max_Path_Len_Diameter"] <- "Diameter"

meltedNetworkProp$Source[meltedNetworkProp$Source == "01_LCA"] <- "LCA_Microsporidia"
meltedNetworkProp$Source[meltedNetworkProp$Source == "02_Encephalitozoon cuniculi"] <- "E.cuniculi"
meltedNetworkProp$Source[meltedNetworkProp$Source == "03_Encephalitozoon hellem"] <- "E.hellem"
meltedNetworkProp$Source[meltedNetworkProp$Source == "04_Encephalitozoon intestinalis"] <- "E.intestinalis"
meltedNetworkProp$Source[meltedNetworkProp$Source == "05_Nosema ceranae"] <- "N.ceranae"
meltedNetworkProp$Source[meltedNetworkProp$Source == "06_Saccharomyces cerevisiae"] <- "S.cerevisiae"

# data$geneID <- factor(data$geneID, levels = unique(data$geneID))
meltedNetworkProp$Source <- factor(meltedNetworkProp$Source, levels = c("LCA_Microsporidia","E.cuniculi","E.hellem","E.intestinalis","N.ceranae","S.cerevisiae"))
meltedNetworkProp$Pathway <- factor(meltedNetworkProp$Pathway, levels = c("Glycolysis / Gluconeogenesis","TCA cycle","Pentose phosphate pathway","Purine metabolism","Pyrimidine metabolism","Amino acid metabolism"))

### plot properties stat
textSize <- 15

### for avg. degree, avg. path length and max path length (density)
meltedNetworkPropSub <- meltedNetworkProp[!(meltedNetworkProp$Property %in% c("Nodes","Edges","Max_Degree")),]

library(dplyr)
meltedNetworkPropSubSummary <- meltedNetworkPropSub %>%
  group_by(Property,Source) %>%
  summarize(mean = mean(Value))

# g <- ggplot(meltedNetworkPropSub, aes(x=Property, y=Value, fill=Property)) +
#   facet_wrap( ~ Source) +
#   geom_violin() + 
#   geom_point(data = meltedNetworkPropSubSummary, aes(y = mean), color = "black", size = 2) +
#   theme(axis.title = element_blank(),axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         strip.text.x = element_text(size = textSize),
#         legend.position = "bottom", legend.text = element_text(size=textSize),legend.title = element_text(size=textSize)) +
#   scale_fill_brewer(palette = "Set2")
# g

# levels(as.factor(meltedNetworkPropSub$Property))
# meltedNetworkPropSub[meltedNetworkPropSub$Property=="Avg. path length",]

g <- ggplot(meltedNetworkPropSub, aes(x=Source, y=Value, fill=Source)) +
  facet_wrap( ~ Property) +
  geom_violin() + 
  geom_point(data = meltedNetworkPropSubSummary, aes(y = mean), color = "black", size = 2) +
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=textSize),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = textSize),
        legend.position = "top", legend.text = element_text(size=textSize),legend.title = element_blank()) +
  scale_fill_brewer(palette = "Set2")
g

# ggsave("network_stat.pdf", width = 30, height = 15, units = "cm")

### for nodes and edges
meltedNetworkPropSub2 <- meltedNetworkProp[meltedNetworkProp$Property %in% c("Nodes","Edges"),]
p <- ggplot(meltedNetworkPropSub2, aes(x=Pathway, y=Value,fill=Source)) +
  facet_wrap( ~ Property) +
  geom_point(aes(col=Source)) +
  coord_flip() +
  labs(y="Count") +
  theme(axis.title.x = element_text(size=textSize), axis.text.x = element_text(size=textSize),
        axis.title.y = element_text(size=textSize), axis.text.y = element_text(size=textSize),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = textSize),
        legend.position = "top", legend.text = element_text(size=textSize),legend.title = element_text(size=textSize)) +
  scale_fill_brewer(palette = "Set2")
p
# ggsave("network_node_edge.pdf", width = 30, height = 10, units = "cm")
