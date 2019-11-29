## Egor Ananyev
## 2019-11-15

## plot_ops.R

## Plotting for eb1/2-an

outpdf = function(n='plot.pdf', w=3.5, h=2.5, s=F, j=T){  # name, width, height, suppress warn
  pdf(paste0('plots/', n, '.pdf'), width=w, height=h)
  if(s){
    suppressWarnings(plot(p))
  } else {
    plot(p)
  }
  x = dev.off()
  if(j){
    png(paste0('plots/', n, '.png'), width=w, height=h, units='in', res=300)
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
    scale_color_manual(values=cb_palette)
  return(p_out)
}

# Plotting box plot data across CTOAs
plot_box_ctoa = function(sss, y_lab, plot_title, out){
  p = ggplot(data=sss, aes(x=factor(ctoa_lowbound), y=RT, fill=Cond)) + 
      geom_boxplot() + xlab('Cue-to-Target Asynchrony (ms)') + ylab(y_lab) +
      theme(panel.grid.minor=element_blank(), plot.title = element_text(hjust = 0.5)) +
      ggtitle(plot_title) + guides(alpha=F)
  if('CuePredFull' %in% colnames(sss)){  # means eb1
    p = p + facet_grid(.~CuePredFull)
    eb_string = 'eb1_ctoa'
  } else {
    eb_string = 'eb2_ctoa'
  }
  p = plot_themefy(p)
  plot(p)
  if(out){
    pdf(paste(eb_string, ' ', plot_title, '.pdf', sep=''), width=8.5, height=2.5)
    plot(p)
    dev.off()
  }
}

# Plotting smooth line plots
