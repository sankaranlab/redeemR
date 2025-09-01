#' Function to plot the mito depth summary
#'
#' This function allows you to plot both position-wise and cell-wise mito depth summary
#' @param ob The redeemR object
#' @param name The plot name shown on top (optional, for backward compatibility)
#' @param w the Width of the plot, default=10 (optional, for backward compatibility)
#' @param h the height of the plot default=3 (optional, for backward compatibility)
#' @return A list containing both the combined patchwork plot and individual plots
#' @examples
#' depth_plots <- plot_depth(redeemR_obj)
#' # Use the combined plot:
#' depth_plots$combined
#' # Or use individual plots:
#' depth_plots$position_depth | depth_plots$cell_depth
#' @export
#' @import ggplot2
#' @importFrom patchwork plot_annotation
plot_depth<-function(ob,name="",w=10,h=3){
require(patchwork)
depth=ob@DepthSummary
d1<-depth[[1]]
d2<-depth[[2]]
names(d1)<-c("pos","meanCov")
names(d2)<-c("cell","meanCov")

# Create position-wise depth plot
p1<-ggplot(d1, aes(pos, meanCov)) + 
    geom_point() + 
    theme_bw() + 
    labs(title = ifelse(name != "", name, "Position-wise Depth"),
         x = "Position", 
         y = "Mean Coverage")

# Create cell-wise depth plot  
p2<-ggplot(d2, aes("cell", meanCov)) + 
    geom_violin() + 
    geom_boxplot() + 
    theme_bw() + 
    labs(title = "Cell-wise Depth",
         x = "", 
         y = "Mean Coverage")

# Combine plots using patchwork
combined_plot <- (p1 | p2) +
    patchwork::plot_annotation(
        title = ifelse(name != "", name, "Depth Summary"),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(family = "sans"))
    ) +patchwork::plot_layout(widths = c(8, 1))

# Return both combined and individual plots
return(list(
    combined = combined_plot,
    position_depth = p1,
    cell_depth = p2
))
}

#' Legacy Function to plot the mito depth summary
#'
#' This function allows you to plot both position-wise and cell-wise mito depth summary
#' @param depth The .depth file by function DepthSummary
#' @param name The plot name shown on top
#' @param w the Width of the plot, default=10
#' @param h the height of the plot default=3
#' @return directly out put the plot
#' @examples
#' plot_depth(DN1CD34_1.depth$Total,"Total")
#' @export
#' @import ggplot2
plot_depth_legacy<-function(depth=DN1CD34_1.depth,name="",w=10,h=3){
d1<-depth[[1]]
d2<-depth[[2]]
names(d1)<-c("pos","meanCov")
names(d2)<-c("cell","meanCov")
options(repr.plot.width=w, repr.plot.height=h)
p1<-ggplot(d1)+aes(pos,meanCov)+geom_point()+theme_bw()
p2<-ggplot(d1)+aes("cell",meanCov)+geom_violin()+geom_boxplot()+theme_bw()
grid.arrange(p1,p2,layout_matrix=rbind(c(1,1,1,1,1,1,2)),top=name)
}

#' Function to plot bulk level mutation signatures
#'
#' This function allows you to plot the mito mutation signatures
#' @param cell_variants a vector of variants formated as c('93_A_G''103_G_A''146_T_C'
#' @return p from ggplot2
#' @examples
#' MutationProfile.bulk(DN1CD34_1.Variants.feature.lst[[name]]$Variants
#' @export
#' @import ggplot2
MutationProfile.bulk<-function(cell_variants){  ## cell_variants look like c('93_A_G''103_G_A''146_T_C''150_C_T''152_T_C''182_C_T')
# Annotate with called variants
data(ref_all_long)
called_variants <- strsplit(cell_variants,"_") %>% sapply(.,function(x){paste(x[1],x[2],">",x[3],sep="")})
ref_all_long$called <- ref_all_long$variant %in% called_variants
# Compute changes in expected/observed
total <- dim(ref_all_long)[1]
total_called <- sum(ref_all_long$called)
prop_df <- ref_all_long %>% group_by(three_plot, group_change, strand) %>%
  dplyr::summarize(observed_prop_called = sum(called)/total_called, expected_prop = n()/total, n = n()) %>%
  mutate(fc_called = observed_prop_called/expected_prop)
prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)
# Visualize
p1 <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.title.x=element_blank(),
        axis.text.x =element_blank())+
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +theme_bw()+
  theme(legend.position = "bottom",axis.text.x  = element_blank(),axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide", y = "Substitution Rate")+
   facet_grid(.~group_change,scales = "free",space="free")
return(p1)
}

#' Function to plot variant metrics
#'
#' This function allows you to plot the mito mutation metrics
#' For each category(stringency),
#' p1: Variant allele frequency(VAF);
#' p2: Heteroplasmy histogram
#' p3: CellN(Number of caells that carry the variants) VS maxcts( The number of variant counts in the highest cell)
#' p4: Histogram to show the distribution of the number of variant per cell
#' @param ob The redeemR object
#' @param p4xlim the p4 xlim(number of variant per cell), default is 50
#' @param QualifyCellCut median coverage for qualified cells, default is 10
#' @return A list containing both the combined patchwork plot and individual plots
#' @export
#' @import ggplot2 ggExtra
#' @importFrom gridExtra grid.arrange
plot_variant<-function(ob,p4xlim=50,QualifyCellCut=10){
options(repr.plot.width=20, repr.plot.height=6)
require(ggExtra)
feature<-ob@V.fitered
GTSummary<-ob@GTsummary.filtered 
depth<-ob@DepthSummary
qualifiedCell<-subset(depth[[2]],meanCov>=QualifyCellCut)[,1,drop=T]  ## Filter Qualified cell based on total depth
#p1
p1<-ggplot(feature)+aes(log10(totalVAF),log10(CV),color=HomoTag)+geom_point(size=0.2)+theme_bw()+scale_color_brewer(palette = "Set1")+theme(axis.text=element_text(size=20,color="black"))
p1<-ggMarginal(p1, type = "histogram",)
#p2
p2<-GTSummary %>% subset(.,!Variants %in% ob@HomoVariants) %>% .$hetero %>% data.frame(HeteroPlasmy=.) %>% ggplot()+aes(log10(HeteroPlasmy))+geom_histogram(binwidth = 0.05)+theme_bw()+theme(axis.text=element_text(size=20,color="black"))
#p3
QualifiedV<-subset(feature,HomoTag=="Hetero")
p3title=paste(nrow(QualifiedV),"Qualified Hetero Variants\nMedian Cell # per V:",median(QualifiedV$CellN),"\nVariants # maxcts>=3:",length(which(QualifiedV$maxcts>=3)))
p3<-ggplot(feature)+aes(log2(CellN),log2(maxcts))+geom_jitter(color="grey80")+geom_point(data=subset(feature,maxcts>=2 & CellN>=2 & HomoTag=="Hetero"),color="black",size=1)+theme_classic()+ggtitle(p3title)
p3<-ggMarginal(p3, type = "histogram",)
#p4
CellVar.Sum<-subset(GTSummary,Variants %in% QualifiedV$Variants & Cell %in% qualifiedCell) %>% group_by(Cell) %>% dplyr::summarise(VN=n(), maxcts=max(Freq),mediancts=median(Freq))
p4title<-paste("Qualified Cell number:",length(qualifiedCell),"\nMedian V number is",median(CellVar.Sum$VN),"\n",CountVperCell(CellVar.Sum$VN,c,CellN=nrow(CellVar.Sum)))
p4<-ggplot(CellVar.Sum)+aes(VN)+geom_histogram(binwidth = 1,color="black",fill="white")+xlim(0,p4xlim)+ggtitle(p4title)+theme(axis.text=element_text(size=20))+geom_vline(xintercept = median(CellVar.Sum$VN),linetype=2)
# Return both combined and individual plots
gridExtra::grid.arrange(p1,p2,p3,p4,ncol=4)
return(list(
    p1 = p1, p2 = p2, p3 = p3, p4 = p4
))
}


#' Legacy Function to plot variant metrics
#'
#' This function works with CW_mgatk.read and Vfilter_v3 
#' This allows you to plot the mito mutation metrics
#' This legacy function is useful to look at all threadhold simultaneous
#' For each category(stringency),
#' p1: Variant allele frequency(VAF);
#' p2: Heteroplasmy histogram
#' p3: CellN(Number of caells that carry the variants) VS maxcts( The number of variant counts in the highest cell)
#' p4: Histogram to show the distribution of the number of variant per cell
#' @param GTSummary GTSummary from CW_mgatk.read
#' @param feature.list feature.list from Vfilter_v3
#' @param cat The catogories, it can be cat = c("Total", "VerySensitive", "Sensitive", "Specific") or a subset
#' @param p4xlim the p4 xlim(number of variant per cell), default is 50
#' @param QualifyCellCut median coverage for qualified cells, default is 10
#' @return no returns, directly plot
#' @examples
#' plot_variant(DN1CD34_1.VariantsGTSummary,DN1CD34_1.Variants.feature.lst,depth=DN1CD34_1.depth,cat=c("Total","VerySensitive","Sensitive","Specific"),p4xlim = 30)
#' @export
#' @import ggplot2 ggExtra
#' @importFrom gridExtra grid.arrange
plot_variant_legacy<-function(GTSummary, feature.list, depth, cat = c("Total", "VerySensitive", "Sensitive", "Specific"), p4xlim = 50, QualifyCellCut = 10) 
{
    options(repr.plot.width = 20, repr.plot.height = 6)
    require(gridExtra)
    require(ggExtra)
    qualifiedCell <- subset(depth[["Total"]][[2]], meanCov >= QualifyCellCut)[, 1, drop = T]
    for (c in cat) {
        p1 <- ggplot(feature.list[[c]]) + aes(log10(VAF), log10(CV), 
            color = HomoTag) + geom_point(size = 0.2) + theme_bw() + 
            scale_color_brewer(palette = "Set1") + theme(axis.text = element_text(size = 20, 
            color = "black"))
        p1 <- ggMarginal(p1, type = "histogram", )
        p2 <- subset(GTSummary[[c]], depth > 20) %>% .$hetero %>% 
            data.frame(HeteroPlasmy = .) %>% ggplot() + aes(HeteroPlasmy) + 
            geom_histogram(binwidth = 0.05) + theme_bw() + theme(axis.text = element_text(size = 20, 
            color = "black"))
        QualifiedV <- subset(feature.list[[c]], maxcts >= 2 & CellN >= 2 & HomoTag == "Hetero")
        p3title = paste(c, ":", nrow(QualifiedV), "Qualified Hetero Variants\nMedian Cell # per V:", 
            median(QualifiedV$CellN), "\nVariants # maxcts>=3:", 
            length(which(QualifiedV$maxcts >= 3)))
        p3 <- ggplot(feature.list[[c]]) + aes(log2(CellN), log2(maxcts)) + 
            geom_jitter(color = "grey80") + geom_point(data = subset(feature.list[[c]], 
            maxcts >= 2 & CellN >= 2 & HomoTag == "Hetero"), 
            color = "black", size = 1) + theme_classic() + ggtitle(p3title)
        p3 <- ggMarginal(p3, type = "histogram", )
        CellVar.Sum <- subset(GTSummary[[c]], Variants %in% QualifiedV$Variants & 
            Cell %in% qualifiedCell) %>% group_by(Cell) %>% dplyr::summarise(VN = n(), 
            maxcts = max(Freq), mediancts = median(Freq))
        p4title<-paste("Qualified Cell number:",length(qualifiedCell),"\nMedian V number is",median(CellVar.Sum$VN),"\n",CountVperCell(CellVar.Sum$VN,c,CellN=nrow(CellVar.Sum)))
        p4<-ggplot(CellVar.Sum)+aes(VN)+geom_histogram(binwidth = 1,color="black",fill="white")+xlim(0,p4xlim)+ggtitle(p4title)+theme(axis.text=element_text(size=20))+geom_vline(xintercept = median(CellVar.Sum$VN),linetype=2)
        return(list(p1,p2,p3,p4))
    }
}

#' Internal function in plot_variant
#'
#' @param x CellVar.Sum$VN
#' @param name c
#' @param CellN nrow(CellVar.Sum)
#' @examples
#' CountVperCell(CellVar.Sum$VN,c,CellN=nrow(CellVar.Sum))
CountVperCell<-function(x,name,CellN){
s<-c(CellN-length(x),length(which(x==1)),length(which(x>=2 &x<=5)),length(which(x>=6 &x<=10)),length(which(x>10)))
names(s)<-c("0","1","2-5","6-10",">10")
paste(paste(s,"Cells have",names(s),"Variants"),collapse="\n")
}


#' Function to plot consensus mtDNA mutation benchmark 
#'
#' This function allows you to plot the mito mutation consensus levels
#' It will print out Quantiles of UMI family size; Quantile of consensus score; Percentage of R1/R2 overlaped mutation detections
#' It will also plot random N mutations as examples to show consensus metrics
#' @param ob The redeemR object
#' @param N number of example variants to show, default is 25
#' @export
#' @import dplyr ggplot2
#' @importFrom gridExtra grid.arrange
Show_Consensus<-function(ob,N=25){
require(dplyr)
require(ggplot2)
if(ob@para["Threhold"]=="T"){
    RawGenotypes<-read.table(paste(ob@attr$path,"/RawGenotypes.Total.StrandBalance",sep=""))
}else if(ob@para["Threhold"]=="LS"){
    RawGenotypes<-read.table(paste(ob@attr$path,"/RawGenotypes.VerySensitive.StrandBalance",sep=""))
}else if(ob@para["Threhold"]=="S"){
    RawGenotypes<-read.table(paste(ob@attr$path,"/RawGenotypes.Sensitive.StrandBalance",sep=""))
}else if(ob@para["Threhold"]=="VS"){
    RawGenotypes<-read.table(paste(ob@attr$path,"/RawGenotypes.Specific.StrandBalance",sep=""))
}    
colnames(RawGenotypes)<-c("ID","Cell","Pos","Call","Alt","Ref","FamSize","GT_cts","CSS","DB_cts","SG_cts","Strand0","Strand1","SensitiveQualifyCTS")    
Example_V<-subset(ob@V.fitered,CellN>25 & HomoTag=="Hetero")$Variants %>% sample(.,N)
RawGenotypes_example<-subset(RawGenotypes,Call %in% Example_V)
###p1
library(reshape2)
RawGenotypes_example$Consensus<-with(RawGenotypes_example,GT_cts/FamSize)
p1<-ggplot(RawGenotypes_example)+aes(Call,log2(FamSize))+geom_boxplot(width=0.5)+coord_flip()+theme_bw()+
    labs(y="Consensus family size")+theme(axis.title.y=element_blank())
###p2
TagC<-function(x){
    tag<-c()
    tag[which(x==1)]<-"100%"
    tag[which(x<1 &x>0.8)]<-"80%"
    tag[which(x<0.8)]<-"75%"
    tag<-factor(tag,levels = c("100%","80%","75%"))
    return(tag)
}
RawGenotypes_example$tag<-TagC(RawGenotypes_example$Consensus)
p2<-RawGenotypes_example[,c("Call","tag")] %>% table %>% as.data.frame %>% ggplot()+aes(Call,Freq,fill=tag)+geom_bar(stat="identity",position="fill")+
    coord_flip()+theme_bw()+scale_fill_brewer(palette="Set2")+theme(axis.title.y=element_blank() )+labs(y="Consensus score")+theme(axis.title.y=element_blank())
##p3
p3<-RawGenotypes_example %>% group_by(Call) %>% dplyr::summarise(DoubleCall=sum(DB_cts),SingleCall=sum(SG_cts)) %>% melt %>% 
    ggplot()+aes(Call,value,fill=variable)+geom_bar(stat="identity",position="fill")+coord_flip()+theme_bw()+scale_fill_brewer(palette="Set1")+theme(axis.title.y=element_blank())+
    theme(axis.title.y=element_blank() )+labs(y="Overlap sequencing(Double call)")+theme(axis.title.y=element_blank())
#p4
p4<-RawGenotypes_example %>% group_by(Call) %>% dplyr::summarise(Strandratio=sum(Strand0)/sum(Strand1)) %>% 
    ggplot()+aes(Call,Strandratio)+geom_point()+ylim(0,2)+coord_flip()+theme_bw()+theme(axis.title.y=element_blank() )+
    theme(axis.title.y=element_blank() )+labs(y="Strand ratio(plus vs minus)")+theme(axis.title.y=element_blank())
##Draw
options(repr.plot.width=18, repr.plot.height=6,repr.plot.res=150,warn=-1)
print("Quantile of UMI family size")
print(quantile(RawGenotypes$FamSize))
print("Quantile of consensus score")
print(quantile(RawGenotypes$CSS))
print("Percentage of R1/R2 overlaped mutation detections")
print(nrow(subset(RawGenotypes,SG_cts==0))/nrow(RawGenotypes))
gridExtra::grid.arrange(p1,p2,p3,p4,nrow=1)
}


#' run_redeem_qc
#'
#' This function generate qc plot to assess filtering strategy
#' @param redeem The redeemR object
#' @param homosets the homoplasmic mutation set
#' @param hotcall make sure hotcall mutations are not included
#' @export
run_redeem_qc <- function(redeem, homosets, hotcall= c("310_T_C","9979_G_A","3109_T_C")){
    print("Make sure add_raw_fragment has been run  after clean_redeem")
    require(ggplot2)   
    require(gridExtra)
    require(igraph)    
    ## check the position distribution</span>
    p_pos<-subset(redeem@raw.fragment.uniqV, !variant %in% c(homosets,hotcall))  %>% ggplot()+aes(rel_position)+geom_histogram(bins=100,fill="white",color="black")+theme_cw1()+ggtitle("1+molecule")
    pos_1mol <- subset(redeem@raw.fragment.uniqV, !variant %in% c(homosets, hotcall)) %>% subset(., Freq==1)  %>% 
    ggplot()+aes(rel_position)+geom_histogram(bins=100,fill="white",color="black")+theme_cw1()+ggtitle("1 mol")
    grid.arrange(p_pos,pos_1mol)

    ## Add transition and transversion information
    redeem@V.fitered <- redeem@V.fitered %>% mutate(changes=add_changes(Variants)) %>% mutate(types=add_types(changes))
    redeem@raw.fragment.uniqV<-redeem@raw.fragment.uniqV %>% mutate(changes=add_changes(variant)) %>% mutate(types=add_types(changes)) 

    p_cell_maxcts <- ggplot(redeem@V.fitered) + aes(log2(CellN), log2(maxcts), color=types) +geom_point()+scale_color_brewer(palette="Set1")+theme_cw1()
    p_cell_meancts <- ggplot(redeem@V.fitered) + aes(log2(CellN), log2(PositiveMean_cts), color=types) +geom_point()+scale_color_brewer(palette="Set1")+theme_cw1()
    grid.arrange(p_cell_maxcts, p_cell_meancts, nrow=1)

    ## report transversion rate
    raw.fragment.uniqV_allhetero <- subset(redeem@raw.fragment.uniqV, !variant %in% homosets)
    raw.fragment.uniqV_1mol<-subset(redeem@raw.fragment.uniqV, Freq==1 & !variant %in% homosets)
    raw.fragment.uniqV_2mol<-subset(redeem@raw.fragment.uniqV, Freq==2 & !variant %in% homosets)
    raw.fragment.uniqV_3molplus<-subset(redeem@raw.fragment.uniqV, Freq>=3 & !variant %in% homosets)
    transversion_rate<-c(unweight=sum(redeem@V.fitered$types=="transversion")/nrow(redeem@V.fitered),
                         weight_freq_all=sum(raw.fragment.uniqV_allhetero$types=="transversion")/nrow(raw.fragment.uniqV_allhetero),
                         weight_freq_1mol=sum(raw.fragment.uniqV_1mol$types=="transversion")/nrow(raw.fragment.uniqV_1mol),
                         weight_freq_2mol=sum(raw.fragment.uniqV_2mol$types=="transversion")/nrow(raw.fragment.uniqV_2mol),
                         weight_freq_3molplus=sum(raw.fragment.uniqV_3molplus$types=="transversion")/nrow(raw.fragment.uniqV_3molplus))
    print(transversion_rate)

    GTsummary.filtered<- redeem@GTsummary.filtered %>% subset(.,!Variants %in% c(homosets,hotcall))

    ##check number of total mutations; mutation per cell, number of cells connected </span>

    Cts.Mtx <- dMcast(GTsummary.filtered,Cell~Variants,value.var = "Freq")
    Cts.Mtx.bi <- Cts.Mtx
    Cts.Mtx.bi[Cts.Mtx.bi>=1]<-1

    number_of_total_mut <- ncol(Cts.Mtx.bi)
    number_of_total_cells <- nrow(Cts.Mtx.bi)
    mut_per_cell <- rowSums(Cts.Mtx.bi)
    redeem_n0.adj.mtx<-CountOverlap_Adj(Cts.Mtx.bi,n=0)
    diag(redeem_n0.adj.mtx)<-0
    number_cell_connected<-rowSums(redeem_n0.adj.mtx) 
    mc <- graph_from_adjacency_matrix(redeem_n0.adj.mtx) %>% igraph::components(mode = "weak") 

    mc <- graph_from_adjacency_matrix(redeem_n0.adj.mtx) %>% igraph::components(mode = "weak") 
    fraction_in_component<-max(mc$csize)/nrow(redeem_n0.adj.mtx)

    print(paste0("number of cells: ", as.character(number_of_total_cells)))
    print(paste0("number of total mutations: ", as.character(number_of_total_mut)))
    print(paste0("number of mutations per cell: ", as.character(median(mut_per_cell))))
    print(paste0("number of cells connected: ", as.character(median(number_cell_connected))))
    print(paste0("fraction of cells in component: ", as.character(median(fraction_in_component))))

    report_metric<-list(number_of_total_mut=number_of_total_mut, 
                        number_of_total_cells=number_of_total_cells,
                        mut_per_cell=mut_per_cell, 
                        fraction_in_component=fraction_in_component,
                        number_cell_connected=number_cell_connected )

    plots <- list(p_pos=p_pos,pos_1mol=pos_1mol, p_cell_maxcts=p_cell_maxcts,p_cell_meancts=p_cell_meancts)

    return(list(plots=plots, transversion_rate=transversion_rate, report_metric=report_metric))
}




#' plot_tree_heatmap_teng_cw
#'
#' Draw a composite phylogenetic tree + cell-by-variant heatmap for mtDNA data.
#'
#' This implementation is derived from the mitodrift visualization by Teng
#' ("plot_redeem_heatmap"), adapted and extended to support multiple stacked
#' annotation bars and consolidated legends.
#'
#' @param gtree A phylo or tbl_graph object representing the clonal tree.
#' @param df_var A data.frame with at least columns: cell, variant, a (alt), d (depth),
#'   or a precomputed vaf column in [0,1].
#' @param draw_depthbranch_width Numeric scaling for branch width drawing.
#' @param root_edge Logical; whether to include a root edge.
#' @param dot_size Numeric size for node/tip dots.
#' @param ylim Optional y limits.
#' @param clade_annot Optional data.frame of clade annotations (columns: cell, clade[, score]).
#' @param tip_annot Optional data.frame of tip annotations (columns: cell, annot).
#' @param title Optional plot title.
#' @param label_site Logical; label mutational sites on the tree.
#' @param cell_annot Either a single data.frame (columns: cell, annot) or a named list
#'   of such data.frames; each element produces a stacked annotation bar.
#' @param tip_lab,node_lab Logical; show tip/node labels.
#' @param layered Logical; layered rendering for annotation bars.
#' @param annot_bar_height,clade_bar_height Numeric heights for stacked bars.
#' @param het_max Numeric VAF max for color scale.
#' @param conf_min,conf_max Numeric limits for node confidence fill.
#' @param conf_label Logical; render numeric confidence labels.
#' @param branch_length Logical; use branch lengths if available.
#' @param node_conf Logical; display internal node confidence.
#' @param annot_scale Optional scaling for annotation bars.
#' @param annot_legend,label_group Logical flags controlling annotation display.
#' @param annot_legend_title Character legend title used when cell_annot is a single table
#'   (named list elements use their own names as titles).
#' @param text_size,label_size Numeric text sizes.
#' @param mut Optional column name in node data to color nodes by.
#' @param post_max,mark_low_cov,facet_by_group,flip Additional display flags.
#'
#' @return A patchwork/ggplot object combining tree, heatmap, and optional bars.
#'
#' @export
#' @import ggplot2 ggraph ggtree patchwork
plot_tree_heatmap_teng_cw <- function (gtree, df_var, draw_depth, branch_width = 0.25, root_edge = TRUE, 
    dot_size = 1, ylim = NULL, clade_annot = NULL, tip_annot = NULL, 
    title = NULL, label_site = FALSE, cell_annot = NULL, tip_lab = FALSE, 
    node_lab = FALSE, layered = FALSE, annot_bar_height = 0.1, 
    clade_bar_height = 1, het_max = 0.1, conf_min = 0.5, conf_max = 0.5, 
    conf_label = FALSE, branch_length = TRUE, node_conf = FALSE, 
    annot_scale = NULL, annot_legend = FALSE, label_group = FALSE, 
    annot_legend_title = "", text_size = 3, label_size = 1, mut = NULL, 
    post_max = FALSE, mark_low_cov = FALSE, facet_by_group = FALSE, 
    flip = FALSE) 
{
    if (inherits(gtree, "tbl_graph")) {
        phylo = to_phylo_reorder(gtree)
    }
    else {
        phylo = gtree
        node_conf = FALSE
    }
    if (!branch_length) {
        phylo$edge.length = NULL
    }
    p_tree = phylo %>% ggtree(ladderize = TRUE, linewidth = branch_width, 
        right = flip) + theme_bw() + theme(plot.margin = margin(0, 
        0, 0, 0, unit = "mm"), axis.title.x = element_blank(), 
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_text(), 
        axis.line.y = element_line(), axis.ticks.y = element_line(), 
        axis.ticks.length.x = unit(0, "pt"), panel.background = element_rect(fill = "transparent", 
            colour = NA), plot.background = element_rect(fill = "transparent", 
            color = NA), panel.grid = element_blank()) + coord_flip() + 
        scale_x_reverse(expand = expansion(mult = 0.05)) + scale_y_continuous(expand = expansion(add = 1)) + 
        ggtitle(title)

    # initialize composition containers BEFORE adding optional components
    heights <- c(1)
    plot_components <- list(p_tree)
    if (!is.null(mut)) {
        dat = gtree %>% activate(nodes) %>% as.data.frame() %>% 
            select(all_of(c("name", p_v = mut)))
        p_tree = p_tree %<+% dat + ggraph::geom_node_point(aes(color = p_v), 
            size = dot_size, stroke = 0) + scale_color_gradient(limits = c(0, 
            het_max), oob = scales::oob_squish)
    }
    if (node_conf) {
        dat = gtree %>% activate(nodes) %>% mutate(isRoot = node_is_root()) %>% 
            as.data.frame() %>% select(any_of(c("name", "isRoot", 
            "conf")))
        if ("conf" %in% colnames(dat)) {
            p_tree = p_tree %<+% dat + geom_nodepoint(aes(fill = conf, 
                subset = !isTip & !isRoot, x = branch), size = dot_size, 
                pch = 22, stroke = 0) + scale_fill_gradient(low = "white", 
                high = "firebrick", limits = c(conf_min, conf_max), 
                oob = scales::oob_squish)
            if (conf_label) {
                p_tree = p_tree + geom_text2(aes(label = round(conf, 
                  2), x = branch, subset = !isTip & !isRoot), 
                  size = text_size, hjust = 0, vjust = 0.25)
            }
        }
    }
    if (!is.null(tip_annot)) {
        dat = tip_annot %>% select(any_of(c(name = "cell", "annot")))
        p_tree = p_tree %<+% dat + geom_tippoint(aes(color = annot, 
            subset = !is.na(annot)), size = dot_size, pch = 19, 
            show.legend = FALSE)
    }
    if (tip_lab) {
        p_tree = p_tree + geom_tiplab(size = label_size, vjust = 0.5, 
            hjust = 1.2, angle = 90) + scale_x_reverse(expand = expansion(mult = c(0.15, 
            0.05)))
    }
    if (node_lab) {
        p_tree = p_tree + geom_text2(aes(label = label, subset = !isTip, 
            x = branch), size = label_size, vjust = 0.5, hjust = 0.5)
    }
    if (!"vaf" %in% colnames(df_var)) {
        df_var = df_var %>% mutate(vaf = a/d)
    }
    cell_order = p_tree$data %>% filter(isTip) %>% arrange(y) %>% 
        pull(label)
    mut_order = order_muts(cell_order, df_var)
    df_var = df_var %>% filter(cell %in% cell_order) %>% mutate(variant = factor(variant, 
        rev(mut_order))) %>% mutate(cell = factor(as.integer(factor(cell, 
        cell_order)), 1:length(cell_order)))
    p_heatmap = df_var %>% ggplot(aes(x = cell, y = variant, 
        fill = vaf)) + geom_raster() + theme_bw() + theme(plot.margin = margin(0, 
        1, 0, 0, unit = "mm"), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = text_size), panel.background = element_rect(fill = "transparent", 
            colour = NA), plot.background = element_rect(fill = "transparent"), 
        panel.grid = element_blank()) + scale_x_discrete(expand = expansion(add = 1), 
        drop = F) + scale_y_discrete(expand = expansion(mult = 0.01)) + 
        scale_fill_gradient(low = "white", high = "red", limits = c(0, 
            het_max), oob = scales::oob_squish) + guides(fill = guide_colorbar(title = "VAF")) + 
        xlab(paste0("cells (n=", length(cell_order), ")")) + 
        ylab(paste0("variants (n=", length(unique(df_var$variant)), 
            ")"))
    if (facet_by_group) {
        p_heatmap = p_heatmap + facet_grid(group ~ ., scales = "free_y", 
            space = "free_y") + theme(panel.spacing.y = unit(0, 
            "mm"))
    }
    if (mark_low_cov) {
        p_heatmap = p_heatmap + geom_point(data = df_var %>% 
            filter(d < 5), pch = 4, size = dot_size)
    }
    if (!is.null(clade_annot)) {
        if ("score" %in% colnames(clade_annot)) {
            show.legend = TRUE
        }
        else {
            clade_annot = clade_annot %>% mutate(score = 1)
            show.legend = FALSE
        }
        p_clade = clade_annot %>% mutate(cell = factor(as.integer(factor(cell, 
            cell_order)), 1:length(cell_order))) %>% ggplot(aes(x = cell, 
            y = factor(clade), fill = score)) + geom_raster(show.legend = show.legend) + 
            theme_bw() + theme(axis.text.x = element_blank(), 
            axis.text.y = element_text(size = text_size), axis.title = element_blank(), 
            axis.ticks = element_blank(), panel.grid = element_blank(), 
            plot.margin = margin(1, 0, 0, 0, unit = "mm"), ) + 
            scale_x_discrete(expand = expansion(add = 1), drop = F) + 
            scale_y_discrete(expand = expansion(add = 1))
        plot_components <- c(plot_components, list(p_clade))
        heights <- c(heights, clade_bar_height)
    }
    if (!is.null(cell_annot)) {
        .make_bar <- function(df, title_label) {
            if (!is.data.frame(df) || !all(c("cell","annot") %in% names(df))) {
                stop("cell_annot items must be data.frames with columns 'cell' and 'annot'")
            }
            df %>%
              dplyr::filter(cell %in% phylo$tip.label) %>%
              dplyr::mutate(cell = factor(cell, cell_order)) %>%
              annot_bar(
                legend       = TRUE,
                label_group  = label_group,
                label_size   = text_size/2,
                annot_scale  = annot_scale,
                legend_title = title_label,
                layered      = layered
              )
        }

        if (is.list(cell_annot)) {
            bar_names <- names(cell_annot)
            if (is.null(bar_names)) bar_names <- rep(annot_legend_title, length(cell_annot))
            bar_list <- mapply(
              function(df, nm) .make_bar(df, nm),
              cell_annot, bar_names,
              SIMPLIFY = FALSE
            )
            plot_components <- c(plot_components, bar_list)
            heights <- c(heights, rep(annot_bar_height, length(bar_list)))
        } else {
            p_bar <- .make_bar(cell_annot, annot_legend_title)
            plot_components <- c(plot_components, list(p_bar))
            heights <- c(heights, annot_bar_height)
        }
    }
    plot_components <- c(plot_components, list(p_heatmap))
    heights <- c(heights, 2)
    wrap_plots(plot_components) + plot_layout(heights = heights, guides = "collect") &
        theme(legend.position = "right")
}



