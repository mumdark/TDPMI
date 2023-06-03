#' For time-dependent progression-related gene set enrichment analysis of myocardial infarction
#'
#' @param genes          a vector of a list of genes
#' @param pvalueFilter   pvalue filter for enrichment
#' @param qvalueFilter   qvalue filter for enrichment
#' @param highcol        color for high statistical significance
#' @param lowcol         color for low statistical significance
#' @param useFilter      Whether to filter the results based on the pvalue and qvalue
#' @param typefigure     The type of the output figure, including "barplot", "bubble" and "all"
#' @param organism       the organism of input genes (using "hsa" or "mmu")
#'
#' @return              a dataframe of TDPMI enrichment result
#' @export
#'
#' @import ggplot2
#' @import msigdbr
#' @import dplyr
#' @import clusterProfiler
#' @import enrichplot
#' @import stringr
#'
#' @examples
#' result=TDPMIenrich(genes=c("Aqp1","Cxcl3","Gm26870","Hba-a1","Hba-a2","Hbb-bs","Hbb-bt","Mb","mt-Atp6","mt-Atp8","mt-Co1","mt-Co2"),organism = "mmu",pvalueFilter=0.05,qvalueFilter=0.05,highcol="grey",lowcol="red3",useFilter="TRUE",typefigure="bubble")
TDPMIenrich=function(genes,
                     pvalueFilter=0.05,
                     qvalueFilter=0.05,
                     highcol="grey",
                     lowcol="red3",
                     useFilter="TRUE",
                     typefigure="bubble", #or "barplot" or "all"
                     organism = "hsa"  # or mmu
){
 # library("msigdbr")
 # library("dplyr")
 # library("clusterProfiler")
 # library("enrichplot")
 # library("ggplot2")
 # library("stringr")
  width <- options()$width*0.91
  cat('[', paste0(rep('#', 1/10*width), collapse=''),
      paste0(rep('-', width - 1/10*width), collapse=''),
      ']',
      round(1/10*100),'%')
 # load("data/cluster.rda")
  if (organism == "hsa") {
    m_df_sub = msigdbr(species = "Rattus norvegicus")
    m_df_sub = m_df_sub[,c(4,7,5,8)] # select only the gene ids necessary
    m_df_sub = distinct(m_df_sub) # remove non unqique rows
    m_df_sub=as.data.frame(na.omit(m_df_sub))
    genes=(m_df_sub[,1])[which(as.vector(m_df_sub[,2]) %in% genes)]
  }else{
    if (organism != "mmu") {
      stop("Incorrect parameter used, the organism should be one of hsa or mmu")
    }
  }
  cat('[', paste0(rep('#', 3/10*width), collapse=''),
      paste0(rep('-', width - 3/10*width), collapse=''),
      ']',
      round(3/10*100),'%')
  cat('[', paste0(rep('#', 4/10*width), collapse=''),
      paste0(rep('-', width - 7/10*width), collapse=''),
      ']',"Id converted successfully")
  egmt <- enricher(genes, TERM2GENE=cluster,pvalueCutoff = 1,qvalueCutoff = 1)
  afrt=as.data.frame(egmt)
  if (length(rownames(afrt))<1){
    next
  }
  if (useFilter %in% c("T","TRUE")) {
    afrt=afrt[(afrt$pvalue<pvalueFilter & afrt$qvalue<qvalueFilter),]
  }else{
    if (useFilter %in% c("F","FALSE")) {
      afrt=afrt
    }else{
      stop("Incorrect command (useFilter) used")
    }
  }
  cat('[', paste0(rep('#', 1/2*width), collapse=''),
      paste0(rep('-', width - 1/2*width), collapse=''),
      ']',
      round(1/2*100),'%')
  if (length(rownames(afrt))<1){
    next
  }
  afrt$Ontology=as.character(as.data.frame(strsplit(afrt$ID,"[.]"))[1,])
  afrt$cluster=as.character(as.data.frame(strsplit(afrt$ID,"[.]"))[2,])
  #write.table(afrt,file=paste0(i,"_","enrich.xls"),sep="\t",quote=F,row.names = F)

  rt=afrt[,c(10,11,11,9,3,5,7)]
  names(rt)=c("Ontology","ID","Term","Count","Ratio","pvalue","qvalue")

  split_b<-str_split(rt$Ratio,"/")
  b<-sapply(split_b,"[",1)
  c<-sapply(split_b,"[",2)
  rt$Ratio=as.numeric(rt$Count)/as.numeric(c[1])

  labels=rt[order(rt$Ratio),"Term"]
  cat('[', paste0(rep('#', 6/10*width), collapse=''),
      paste0(rep('-', width - 6/10*width), collapse=''),
      ']',
      round(6/10*100),'%')
  p = ggplot(rt,aes(Ratio, Term)) +
    geom_point(aes(size=Count, color=qvalue))
  p1 = p +
    scale_colour_gradient(high=highcol, low = lowcol) +
    labs(color="qvalue",size="Count",x="Gene ratio",y="Cluster")+
    theme(axis.text.x=element_text(color="black", size=10),axis.text.y=element_text(color="black", size=10)) +
    scale_size_continuous(range=c(4,9))+
    theme_bw()+
    facet_grid( Ontology~. ,scales="free")+
    theme(strip.text.x = element_text(size = 10,colour = "black"))+
    theme(strip.background.x = element_rect(fill = c("white"), colour = "grey"),strip.background.y = element_rect(fill = c("white"), colour = "grey"))
  #ggsave("bubble.pdf", width=5.5, height=9)

  rt=afrt[,c(10,11,11,9,3,5,7)]
  names(rt)=c("Ontology","ID","Term","Count","Ratio","pvalue","qvalue")
  cat('[', paste0(rep('#', 4/5*width), collapse=''),
      paste0(rep('-', width - 4/5*width), collapse=''),
      ']',
      round(4/5*100),'%')
  labels=rt[order(rt$qvalue,decreasing =T),"Term"]

  p=ggplot(data=rt)+geom_bar(aes(x=Term, y=Count, fill=qvalue), stat='identity')+
    coord_flip() + scale_fill_gradient(high=highcol, low = lowcol) +
    xlab("Cluster") + ylab("Gene count") +
    theme(axis.text.x=element_text(color="black", size=10),axis.text.y=element_text(color="black", size=10)) +
    scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0,0))+
    theme_bw()+
    facet_grid( Ontology~. ,scales="free")+
    theme(strip.text.x = element_text(size = 10,colour = "black"))+
    theme(strip.background.x = element_rect(fill = c("white"), colour = "grey"),strip.background.y = element_rect(fill = c("white"), colour = "grey"))
  #ggsave("barplot.pdf", width=5.5, height=9)
  cat('[', paste0(rep('#', 9/10*width), collapse=''),
      paste0(rep('-', width - 9/10*width), collapse=''),
      ']',
      round(9/10*100),'%')
  if (typefigure=="barplot") {
    print(p)
  }else{
    if (typefigure=="bubble") {
      print(p1)
    }else{
      print(p1+p)
    }

  }
  cat('[', paste0(rep('#', 1*width), collapse=''),
      paste0(rep('-', width - 1*width), collapse=''),
      ']',
      round(1*100),'%')
  return(afrt)
}
