getAbundanceStats <- function(x, k = "meta20", group_by = "condition", shape_by = NULL){
  
  CATALYST:::.check_sce(x, TRUE)
  k <- CATALYST:::.check_k(x, k)
  CATALYST:::.check_cd_factor(x, group_by)
  
  ns <- table(cluster_id = cluster_ids(x, k), sample_id = sample_ids(x))
  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)
  
  #get the positions of each sample in x
  m <- match(df$sample_id, x$sample_id)
  
  for (i in c(shape_by, group_by)) df[[i]] <- x[[i]][m]
  
  df
}


plotAbundancesWithStats_ggarrange <- function(df, by="cluster_id", group_by = "condition", 
                                              outfile="abundances.pdf", shape_by = NULL,
                                              group_by_palette = CATALYST:::.cluster_cols,
                                              group_by_stats_ref_group = NULL, 
                                              group_by_stats_comparisons = NULL, 
                                              group_by_stats_method = "t.test",
                                              nrows=3,width_single=40,height_single=60){
  
  # df <- sce_all_day1_abundances
  # group_by <- "Tx_Group"
  # group_by_palette <-Tx_color_palette
  # comparisonList <- list(c("Placebo","PolyICLC"),c("Placebo","Resiquimod"))
  
  cluster_list <- levels(df %>% pull(.,var=by))
  i=1
  panels <- list()
  
  if(is.null(shape_by)){
    columnToSelect <- c(group_by,"Freq")
  }else{
    columnToSelect <- c(group_by,shape_by,"Freq")
  }
  
  for(i in 1:length(cluster_list)){
    
    subtable <- df %>% dplyr::filter(cluster_id == cluster_list[i]) %>%
      dplyr::select(all_of(columnToSelect))
    
    subtable <- subtable %>% relocate(Freq, .after = last_col())
    
    y1 = ggboxplot(subtable, 
                   x = group_by, 
                   y = "Freq",
                   color = group_by, 
                   palette = group_by_palette,
                   size=1.5, outlier.shape = NA)
    
    if(is.null(shape_by)){      
      y1 = y1 + geom_point(aes_(fill=as.name(group_by)), 
                           size = 2, shape=21, 
                           position = position_jitterdodge())
    }else{
      y1 = y1 + geom_point(aes_(fill=as.name(group_by),shape=as.name(shape_by)), 
                           size = 2, col="black", position = position_jitterdodge())+
        scale_shape_manual(values=21:22)
    }
    
    y1 <- y1 + theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=18),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      plot.title=element_text(size=18))+
      
      ggtitle(cluster_list[i])
    
    if(!is.null(group_by_stats_ref_group)){
      
      y1<- y1 + scale_y_continuous(expand=expansion(mult=c(0.1,0.1))) + 
        stat_compare_means(ref.group = group_by_stats_ref_group, 
                           size=5, hide.ns = T, method = group_by_stats_method, label = "p.format",vjust=-1)
      
    }else{
      
      if(!is.null(group_by_stats_comparisons)){
        y1<- y1 + scale_y_continuous(expand=expansion(mult=c(0.1,0.1))) + 
          stat_compare_means(comparisons = group_by_stats_comparisons, 
                             size=5, hide.ns = T, method = group_by_stats_method,label = "p.format",vjust=-1)
      }else{
        ## default comparison (can be pairwise or group)
        y1<- y1 + scale_y_continuous(expand=expansion(mult=c(0.1,0.1))) + stat_compare_means(size=5, hide.ns = T, vjust=-1)  
      }
    }
    
    panels[[i]]=y1
  }
  
  ncols <- floor(length(cluster_list)/nrows)+1
  y <- ggarrange(plotlist=panels,ncol=ncols,nrow=nrows)
  
  ggsave(y, filename=outfile,
         width=width_single*ncols,
         height=height_single*nrows,
         units="mm")
}
