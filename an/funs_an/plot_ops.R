## Egor Ananyev
## 2019-11-15

## plot_ops.R

## Plotting for eb1/2-an

outpdf = function(n='plot', w=3.5, h=2.5, s=F, j=T){  # name, width, height, suppress warn
  # pdf(paste0('plots/', n, '.pdf'), width=w, height=h)
  pdf(paste0(n, '.pdf'), width=w, height=h)
  if(s){
    suppressWarnings(plot(p))
  } else {
    plot(p)
  }
  x = dev.off()
  if(j){
    # png(paste0('plots/', n, '.png'), width=w, height=h, units='in', res=300)
    png(paste0(n, '.png'), width=w, height=h, units='in', res=600)
    if(s){
      suppressWarnings(plot(p))
    } else {
      plot(p)
    }
    x = dev.off()
  }
}

plot_themefy = function(p_in){
  # Colorblind-friendly palette:
  cb_palette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
                 "#CC79A7")
  p_out = p_in + theme_bw() + scale_fill_manual(values=cb_palette) +
    scale_color_manual(values=cb_palette) + 
    theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(),
       axis.text=element_text(size=8), axis.title=element_text(size=9),
       legend.text=element_text(size=8), legend.title=element_text(size=9),
       legend.key = element_blank(), legend.margin=margin(t=-.04, unit='in'),
       legend.background = element_rect(fill='transparent'),
       plot.title=element_text(face='bold'))
  return(p_out)
}

# Plotting box plot data across CTOAs
plot_box_ctoa = function(sss, y_lab, plot_title, out, x_str='ctoa_lowbound', eb='eb', l=T){
  sss$x = sss[,x_str]  # this allows plotting non-binned CTOAs; x_str needs to change to a diff var
  p = ggplot(data=sss, aes(x=factor(x), y=RT, fill=Cond)) + 
      geom_boxplot() + xlab('Cue-to-Target Asynchrony (ms)') + ylab(y_lab) +
      theme(panel.grid.minor=element_blank(), plot.title = element_text(hjust = 0.5)) +
      ggtitle(plot_title) + guides(alpha=F)
  if(eb == 'eb'){
    if('CuePredFull' %in% colnames(sss)){  # means eb1
      p = p + facet_grid(.~CuePredFull)
      eb_string = 'eb1_ctoa'
    } else {
      eb_string = 'eb2_ctoa'
    }
  } else {
    eb_string = paste0(eb, '_ctoa')
  }
  p = plot_themefy(p)
  if(!l){  # if legend is unnecessary
    p = p + theme(legend.position='none')
  }
  plot(p)
  if(out){
    pdf(paste(eb_string, ' ', plot_title, '.pdf', sep=''), width=8.5, height=2.5)
    plot(p)
    dev.off()
  }
}

# Bar plots for single-subject visualization:
plot_bar_ctoa = function(sss, y_lab, plot_title, out=F, x_str='targ_soa', eb = 'eb'){
  sss$x = sss[,x_str]  # this allows plotting non-binned CTOAs; x_str needs to change to a diff var
  p = ggplot(data=sss, aes(x=factor(x), y=RT, fill=Cond)) + 
      geom_bar(stat="identity", position=position_dodge(), colour='black') +
      xlab('Cue-to-Target Asynchrony (ms)') + ylab(y_lab) +
      theme(panel.grid.minor=element_blank(), plot.title=element_text(hjust=0.5)) +
      ggtitle(plot_title) + guides(alpha=F)
  if(eb == 'eb'){
    if('CuePredFull' %in% colnames(sss)){  # means eb1
      p = p + facet_grid(.~CuePredFull)
      eb_string = 'eb1_ctoa'
    } else {
      eb_string = 'eb2_ctoa'
    }
  } else {
    eb_string = paste0(eb, '_ctoa')
  }
  p = plot_themefy(p)
  plot(p)
  if(out){
    pdf(paste(eb_string, ' ', plot_title, '.pdf', sep=''), width=8.5, height=2.5)
    plot(p)
    dev.off()
  }
}

plot_acc = function(acc_ss, cue_type, plot_title, out=F,
                    x_lab='Cue-to-Target Asynchrony (ms)', y_lab='Accuracy'){
  p = ggplot(data=acc_ss[acc_ss$CueValid==cue_type, ], aes(x=factor(targ_soa), y=acc, fill=Cond)) +
      geom_bar(stat="identity", position=position_dodge(), colour='black') +
      xlab(x_lab) + ylab(y_lab) +
      theme(panel.grid.minor=element_blank(), plot.title=element_text(hjust=0.5)) +
      ggtitle(plot_title)
  p = plot_themefy(p)
  p = p + theme(legend.position='none')
  plot(p)
  if(out){
    pdf(paste(eb_string, ' ', plot_title, '.pdf', sep=''), width=3.5, height=2.5)
    plot(p)
    dev.off()
  }
}
