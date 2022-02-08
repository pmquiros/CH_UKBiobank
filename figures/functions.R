
"%ni%" <- Negate("%in%")

not_all_na <- function(x) any(!is.na(x))

dup_str <- function(x, na.rm=T) {
  x2 <- str_split(x, ",")
  x3 <- paste(unique(unlist(x2)), collapse = ',')
  return(x3)
}

pheno_seeker <- function(name="all", list_icd10=NULL, list_icd9 = NULL, list_self_cancer=NULL, list_self_illness=NULL, list_operation=NULL, list_opcs4=NULL, keep = FALSE) {
  
  # list_icd10 = Data-Coding 19
  # list_icd9 = Data-Coding 87
  # list_self_cancer = Data-Coding 3
  # list_self_illness = Data-Coding 6
  # list_operation = Data-Coding 5
  # list_opcs4 = Data-Coding 240
  
  res <- data.frame(matrix(ncol = 1, nrow = 0))
  colnames(res) <- "eid"
  
  
  
  if (!is.null(list_icd10)) {
    
    res_dis <- pheno_icd10_all %>%
      mutate(icd10 = case_when(disease_code %in% list_icd10 ~ diag_date)) %>%
      select(eid, icd10) %>%
      group_by(eid) %>% arrange(icd10) %>% slice(1) %>% ungroup() %>%
      mutate(icd10 = lubridate::ymd(icd10)) %>%
      drop_na
    
    res_death <- pheno_death %>%
      mutate(icd10_death = case_when(cause_icd10 %in% list_icd10 ~ dod)) %>%
      select(eid, icd10_death) %>%
      group_by(eid) %>% arrange(icd10_death) %>% slice(1) %>% ungroup() %>%
      mutate(icd10_death = lubridate::ymd(icd10_death)) %>%
      drop_na
    
    res_icd10 <- full_join(res_dis, res_death, by=("eid"))
    res <- res %>% full_join(res_icd10, by=("eid"))
  }
  
  if (!is.null(list_icd9)) {
    
    res_dis_icd9 <- pheno_icd9_all %>%
      mutate(icd9 = case_when(disease_code %in% list_icd9 ~ diag_date)) %>%
      select(eid, icd9) %>%
      group_by(eid) %>% arrange(icd9) %>% slice(1) %>% ungroup() %>%
      mutate(icd9 = lubridate::ymd(icd9)) %>%
      drop_na
    
    res <- res %>% full_join(res_dis_icd9, by=("eid"))
  }
  
  if (!is.null(list_self_cancer)) {
    
    res_self_cancer <- pheno_self_cancer %>%
      mutate(self_cancer = case_when(disease_code %in% list_self_cancer ~ diag_date)) %>%
      select(eid, self_cancer) %>%
      group_by(eid) %>% arrange(self_cancer) %>% slice(1) %>% ungroup() %>%
      mutate(self_cancer = lubridate::ymd(self_cancer)) %>%
      drop_na
    
    res <- res %>% full_join(res_self_cancer, by=("eid"))
    
  }
  
  if (!is.null(list_self_illness)) {
    
    res_self_illness <- pheno_self_cancer %>%
      mutate(self_illness = case_when(disease_code %in% list_self_illness ~ diag_date)) %>%
      select(eid, self_illness) %>%
      group_by(eid) %>% arrange(self_illness) %>% slice(1) %>% ungroup() %>%
      mutate(self_illness = lubridate::ymd(self_illness)) %>%
      drop_na
    
    res <- res %>% full_join(res_self_illness, by=("eid"))
    
  }
  
  if (!is.null(list_operation)) {
    
    res_operation <- pheno_operations %>%
      mutate(operation = case_when(disease_code %in% list_operation ~ diag_date)) %>%
      select(eid, operation) %>%
      group_by(eid) %>% arrange(operation) %>% slice(1) %>% ungroup() %>%
      mutate(operation = lubridate::ymd(operation)) %>%
      drop_na
    
    res <- res %>% full_join(res_operation, by=("eid"))
    
  }
  
  if (!is.null(list_opcs4)) {
    
    res_opcs4 <- pheno_opcs4_all %>%
      mutate(opcs4 = case_when(disease_code %in% list_opcs4 ~ diag_date)) %>%
      select(eid, opcs4) %>%
      group_by(eid) %>% arrange(opcs4) %>% slice(1) %>% ungroup() %>%
      mutate(opcs4 = lubridate::ymd(opcs4)) %>%
      drop_na
    
    res <- res %>% full_join(res_opcs4, by=("eid"))
    
  }
  
  if(keep == TRUE) {
    
    res <- res %>%
      mutate(all = select(.,-matches("eid")) %>% 
               mutate_all(as.numeric) %>% reduce(.,na.rm=T, pmin)) %>%
      mutate(all = lubridate::as_date(as.numeric(all))) %>%
      rename(!!name:= all)
    
  } else {
    
    res <- res %>%
      mutate(all = select(.,-matches("eid")) %>% 
               mutate_all(as.numeric) %>% reduce(.,na.rm=T, pmin)) %>%
      mutate(all = lubridate::as_date(as.numeric(all))) %>%
      select(eid, !!name:= all)
    
  }
  
  
  return(res)
  
}


paste3 <- function(..., sep=",") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}



# test probability that a mutation is true. Binomial test (x=ALT, n=DP, p=0.5). Filter pval <= 0.001
vaf_bin <- function(ref, alt, chr, sex, p=.5){
  #ref = Ref allele counts
  #alt = Alt allele counts
  #sex c("Male", "Female")
  #chr: chromosome chr1-22, chrX, chrY
  
  sexchr = c("chrX", "chrY")
  
  if(sex == "Female") {
    pval <- binom.test(as.numeric(alt), as.numeric(ref)+as.numeric(alt), p=p, alternative = "less")$p.value
    return(pval)  
  } else if (sex == "Male") {
    if (chr %in% sexchr ) {
      pval <- binom.test(as.numeric(alt), as.numeric(ref)+as.numeric(alt), p=.95, alternative = "less")$p.value
      return(pval) 
    } else {
      pval <- binom.test(as.numeric(alt), as.numeric(ref)+as.numeric(alt), p=p, alternative = "less")$p.value
      return(pval) 
    }
  }
  
}



# this tidier uses a Wald method for calculating the CIs
tidy_ci_exp <- function(x, exponentiate =  TRUE, conf.level = 0.95, ...) {
  if (exponentiate == TRUE) {
    ci <- as_tibble(exp(confint.default(x, conf.level = conf.level)), rownames = "term")
  } else {
    ci <- as_tibble(confint.default(x, conf.level = conf.level), rownames = "term")
  }
  names(ci) <- names(ci) <- c("term", "conf.low", "conf.high")
  full_join(
    broom::tidy(x, exponentiate = exponentiate, conf.int = FALSE),
    ci, by="term"
  )
  
}


######## forestplot - new update - log scale and grouped data

forestplot_from_tidy <- function (x, type = c('Hazard ratio', 'Odds ratio', 'Risk ratio'), 
                                  pval = NULL, cutoff = 0.05, color = TRUE,
                                  title = NULL, shadows=TRUE,
                                  min=NULL, max=NULL, table = TRUE,
                                  order=NULL, log = FALSE,
                                  group = FALSE) {
  
  if (is.null(pval)) pval = "nominal"
  if (type == "Hazard ratio") {
    type_name = "HR [95% CI]"
  } else if (type == "Risk ratio") {
    type_name = "RR [95% CI]"
  } else if (type == "Odds ratio") {
    type_name = "OR [95% CI]"
  } else { stop("type required ('Hazard ratio', 'Odds ratio' or 'Risk ratio'") }
  
  if (table == T & group == T) {stop("Table and group are incompatible options")}
  
  x <- x %>% 
    mutate(color = ifelse(p_value <= cutoff, "blue", "black")) %>%
    mutate(p_value = formatC(p_value, format = "e", digits = 2)) %>% 
    mutate(estimate_ci = paste0(format(estimate, digits=2, nsmall=2, trim=T), " [", format(conf_low, digits=2, nsmall=2, trim=T), "-", format(conf_high, digits=2, nsmall=2, trim=T), "]"))
  
  if (pval == "fdr") {
    x <- x %>% 
      mutate(color = ifelse(fdr <= cutoff, "blue", "black")) %>%
      mutate(fdr = formatC(fdr, format = "e", digits = 2))
  }
  
  
  yint = 1
  
  
  
  if (is.null(max)) {
    if(max(x$conf_high) < 10) {max_y = ceiling(max(x$conf_high))}
    if(max(x$conf_high) >= 10) {max_y = ceiling(max(x$conf_high)/10)*10}
  } else { 
    max_y = max
  }
  
  
  
  if (is.null(min)) {
    min_y = 0
  } else { 
    min_y = min
  }
  
  
  # Reorder names
  
  if (is.null(order)) {
    x <- x %>% 
      mutate(name = fct_reorder(name, estimate))
  } else if (is.vector(order)) {
    x <- x %>% 
      mutate(name = factor(name, levels = rev(order)))
  } else {print("Error in order")}
  
  n_samples=nrow(x)
  
  if (group == T){
    p1 <- ggplot(x, aes(y= estimate, x = name, shape=group))
  } else {
    p1 <- ggplot(x, aes(y= estimate, x = name))
  }
  
  
  
  if (shadows == TRUE) {
    
    if (n_samples > 1) {sel <- seq(2, n_samples, 2)} else {sel = 0}
    
    p1 <- p1 +
      geom_rect(data = x[sel,],
                xmin=sel-0.5, xmax=sel+.5,
                ymin = -Inf, ymax = Inf,
                fill = 'grey', alpha = 0.5,
                inherit.aes = FALSE)
  }
  
  if (color == T & group == F) {
    
    p1 <- p1 + 
      geom_point(size=2, color=x$color, fill=x$color) +
      geom_errorbar(aes(ymin= conf_low, ymax=conf_high), width=.1, color=x$color)
    
  } else if (color == T & group == T) {
    
    p1 <- p1 + 
      geom_point(size=2, color=x$color, fill=x$color, position = position_dodge(.8)) +
      scale_shape_manual(values=c(19, 15))+
      geom_errorbar(aes(ymin= conf_low, ymax=conf_high), width=.1, color=x$color,
                    position = position_dodge(.8)) 
    
    
  } else if (color == F & group == T) {
    
    p1 <- p1 + 
      geom_point(size=2, position = position_dodge(.8)) +
      scale_shape_manual(values=c(19, 15))+
      geom_errorbar(aes(ymin= conf_low, ymax=conf_high), width=.1,
                    position = position_dodge(.8))
    
    
  } else {
    
    p1 <- p1 + 
      geom_point(size=2) +
      geom_errorbar(aes(ymin= conf_low, ymax=conf_high), width=.1)
  }
  
  #Add to set number of breaks = 5
  seq_lim <- seq(min_y, max_y)
  if(length(seq_lim) < 5){
    if(5-length(seq_lim) == 1) {max_y = max_y + seq_lim[2]-seq_lim[1]}
    if(5-length(seq_lim) == 2) {
      if (min_y == 0) {
        max_y = max_y + 2*seq_lim[2]-seq_lim[1]
      } else {
        max_y = max_y + (seq_lim[2]-seq_lim[1])
        min_y = min_y - (seq_lim[2]-seq_lim[1])
      }
    }
  }
  
  p0 <- p1 + 
    scale_y_continuous(expand = c(0, 0), n.breaks = 5, limits = c(min_y,max_y)) +
    coord_flip() 
  
  
  if (table == TRUE) {
    
    b = ggplot_build(p0)$layout$panel_params[[1]]$x$breaks
    b=b[!is.na(b)]
    seg = (max(b)-min(b))/(length(b)-1)
    max_y_table=max(b)+5*seg
    pos_est = max(b) + 2*seg
    pos_p = max(b) + 4*seg
    lb = b
    
    if(min_y < yint) { 
      p1 <- p1 + 
        geom_hline(yintercept = yint, linetype=2)}
    
    p1 <- p1 +
      scale_x_discrete(expand = expansion(add = 0.5)) +
      coord_flip(clip = "off") + 
      theme_c2()+
      theme(axis.title.y = element_blank(),
            axis.line.x = element_blank(),
            axis.title.x = element_text(hjust = 0.2),
            plot.title = element_text(size = 18),
            plot.margin = unit(c(2,3,1,1), "lines")) +
      geom_segment(aes(x=0.5,xend=0.5,y=min_y,yend=max_y))+
      labs(y = type_name, title = title) +
      
      geom_text(aes(label = estimate_ci, x = name, y = pos_est), size = 5) +
      annotate("text", x = n_samples+.5, y = pos_est, label = type_name, size = 5, vjust=-.5)
    
    if(log == T){
      p1 <- p1 + scale_y_continuous(expand = c(0, 0), limits = c(min_y,max_y_table), breaks = b, labels = lb, trans = 'log2')
    } else {
      p1 <- p1 + scale_y_continuous(expand = c(0, 0), limits = c(min_y,max_y_table), breaks = b, labels = lb)
    }
    
    if (pval == "nominal") {
      p1 <- p1 +
        geom_text(aes(label = p_value, x = name, y = pos_p), size = 5) +
        annotate("text", x = n_samples+.5, y = pos_p, label = "p-val", size = 5, vjust=-.5)
      
    } else if (pval == "fdr") {
      p1 <- p1 +
        geom_text(aes(label = fdr, x = name, y = pos_p), size = 5) +
        annotate("text", x = n_samples+.5, y = pos_p, label = "FDR", size = 5, vjust=-.5)
      
    }
    
  } else {
    
    p1 <- p1 + 
      geom_hline(yintercept = yint, linetype=2) +
      scale_x_discrete(expand = expansion(add = 0.5)) +
      coord_flip() + 
      theme_c2()+
      theme(axis.title.y = element_blank(),
            plot.title = element_text(size = 18),
            plot.margin = unit(c(2,3,1,1), "lines")) +
      labs(y = type_name, title = title) 
    
    if (log == TRUE) {
      p1 <- p1 + scale_y_continuous(expand = c(0, 0), n.breaks = 5, limits = c(min_y,max_y), trans = 'log2')#, labels = no_zero)
    } else { 
      p1 <- p1 + scale_y_continuous(expand = c(0, 0), n.breaks = 5, limits = c(min_y,max_y))
    }
    
    if (group == T) {
      p1 <- p1 + theme(legend.title = element_blank(),
                       legend.position = 'bottom', legend.direction = 'horizontal',
                       legend.text = element_text(size=12)) +
        guides(shape = guide_legend(reverse = T))
    }
    
  }
  
  
  return(p1)
  
}



data_summary <- function(x) {
  m <- mean(x); 
  ymin <- m-sd(x); 
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

data_summary_CI <- function(x) {
  m <- mean(x, na.rm=TRUE)
  sem <-sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))
  ymin<-m-1.96*sem
  ymax<-m+1.96*sem
  return(c(y=m,ymin=ymin,ymax=ymax))
}

theme_c2 <- function () {
  theme_classic() %+replace% 
    theme(
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18)
    )
}

clean_title <- function(s, caps = FALSE){
  n=str_count(s)
  if (caps == TRUE) {res=toupper(s)}
  else {
    res=paste0(toupper(str_sub(s,1,1)), str_sub(s,2))
    res=str_replace_all(res, "\\.|_", " ")
  }
  return(res)
}



reverselog_trans <- function(base = exp(1)) {     
  trans <- function(x) -log(x, base)     
  inv <- function(x) base^(-x)     
  scales::trans_new(paste0("reverselog-", format(base)), 
                    trans, 
                    inv,
                    scales::log_breaks(base = base),
                    domain = c(1e-100, Inf)) 
}


waterfall_plot <- function (df, top=NULL, log=NULL) {
  
  samples <- df %>%
    select(eid, symbol, mutation_type) %>%
    mutate(eid=as.character(eid))
  
  genes <- df %>%
    group_by(symbol) %>%
    count() %>%
    ungroup() %>%
    mutate(perc = round(100*n / sum(n),1)) %>%
    arrange(-n)
  
  if(!is.null(top)) {
    genes <- head(genes, top)
    samples <- samples %>% filter(symbol %in% genes$symbol)
  }
  
  sort_genes <- genes %>% pull(symbol)
  
  #select only one mutation type per eid/gene
  samples <- samples %>% 
    mutate(symbol = fct_relevel(symbol, sort_genes)) %>%
    arrange(symbol, mutation_type) %>% 
    group_by(eid, symbol) %>% 
    slice(1) %>% 
    ungroup()
  
  sort_samples <- samples %>% 
    arrange(symbol, mutation_type) %>% select(eid) %>% distinct %>% pull(eid) 
  
  plot_var <- ggplot(samples, aes(x=eid, y=symbol, fill=factor(mutation_type))) + 
    # geom_rect(aes(xmin = 0.5, xmax = length(unique(sort_samples))+.5, ymin = 0.5, ymax = length(sort_genes)+.5), 
    #           fill = "grey92")+
    # geom_vline(xintercept=seq(0.5, 1+length(unique(sort_samples))-0.5, 1), 
    #            lwd=.5, colour="white")+
    geom_tile()+
    scale_y_discrete(limits=rev(sort_genes))+
    scale_x_discrete(limits=sort_samples) +
    scale_fill_manual(values=mut_type_colors)+
    theme(axis.title= element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(size = 18, color='black', face="italic"),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.title = element_blank(), legend.text = element_text(size = 10))+
    geom_hline(yintercept=seq(0.5, 1+length(sort_genes)-0.5, 1), 
               lwd=.5, colour="white") 
  
  
  tmp_all <- ggplot_gtable(ggplot_build(plot_var))
  leg_all <- which(sapply(tmp_all$grobs, function(x) x$name) == "guide-box")
  legend_all <- tmp_all$grobs[[leg_all]]
  
  plot_var2 <- plot_var + theme(legend.position="none")+
    theme(plot.margin = unit(c(1,1,1,1), "pt"))
  
  
  lim_y_genes = max(genes$n)+0.2*max(genes$n)
  lim_y_genes_log= 10^(ceiling(log10(max(genes$n))))
  
  
  
  plot_genes <- ggplot(genes, aes(x=reorder(symbol,n), y=n))+
    geom_bar(stat="identity",fill="steelblue")+
    theme_c2()+ coord_flip(clip = "off") + 
    theme(axis.text.x = element_text(), axis.title.x = element_text(), 
          axis.title.y = element_blank(), axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid = element_blank(), panel.border = element_blank()) +
    labs( x = "Gene", y = "No of mutations")
  
  if (is.null(log)){
    
    plot_genes <- plot_genes +
      scale_y_reverse(limits = c(lim_y_genes, 0), expand = c(0, 0)) +
      geom_text(aes(label=paste0(perc, "%")), hjust=1.1, color="black", size=4)
  } else if (log == TRUE) {
    
    plot_genes <- plot_genes +
      scale_y_continuous(trans=reverselog_trans(10), 
                         limits = c(lim_y_genes_log, 1), expand = c(0, 0),
                         labels=scales::trans_format('log10',scales::math_format(10^.x))
                         ) +
      theme(plot.margin = unit(c(1, 15, 1, 12), "pt"))+
      
      geom_text(aes(label=paste0(perc, "%")), hjust=-0.2, color="white", size=4) 
    
  }
  
  
  plot_vaf = df %>% 
    filter(symbol %in% genes$symbol) %>% 
    mutate(symbol = factor(symbol, levels=genes$symbol)) %>% 
    ggplot(., aes(x=fct_rev(symbol), y=vaf))+
    geom_violin(draw_quantiles = 0.5, fill = '#999999')+
    stat_summary(fun.data=data_summary, size=.6) +
    theme_c2()+ 
    coord_flip(clip = "off")+ 
    theme(axis.text.x = element_text(), axis.title.x = element_text(), 
          axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid = element_blank(), panel.border = element_blank(),
          #plot.margin = unit(c(11, 5.5, 5.5, 5.5), "pt"),
          legend.position = 'none') +
    scale_y_continuous(expand = c(0, 0), breaks=c(0,.25, .5, .75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    labs(y = "VAF") 
  
  
  ##
  
  
  plot_grid(plot_genes, plot_var, plot_vaf, nrow = 1, rel_widths = c(1,2,1),
            align = "h", axis = "lrtb")
  
  
  
}

inverse_normal <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}


tidy_crr <- function(mod) {
  estimate = summary(mod)$conf.int[, 1]
  std.error = summary(mod)$coef[, 3]
  statistic = summary(mod)$coef[, 4]
  #calculate manually the p-values
  p.value = 2 * pnorm(-abs(mod$coef)/sqrt(diag(mod$var)))
  conf.low = summary(mod)$conf.int[, 3]
  conf.high = summary(mod)$conf.int[, 4]
  term = row.names(summary(mod)$coef)
  res <- tibble(term, estimate, std.error, statistic, p.value, conf.low, conf.high)
  return(res)
}


