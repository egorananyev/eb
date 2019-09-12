
# Visualize the data:
library(ggplot2)
plot_qc = function(samples, trials){
    ylims = c(-10, 10)
    p = ggplot(samples[samples$trial==cur_trial,], aes(x=time))
    p = p + geom_line(aes(y=xr, colour='Gaze'))
    # p = p + geom_line(aes(y=yr-432, colour='Gaze'), linetype='longdash')
    p = p + geom_line(aes(y=((psr-mean(psr))/100), colour='Pupil'))  # colour='Pupil'))
    p = p + scale_y_continuous(sec.axis = sec_axis(~.*5, name='Pupil Diameter (a.u.)'),limits=ylims)
    # p = p + xlim(xlims)
    p = p + scale_colour_manual(values=c('blue', 'magenta'))
    # number of blanks in this trial:
    numof_blanks = sum(blanks$trial==cur_trial)
    if(numof_blanks >= 1){
        cur_blanks = blanks[blanks$trial==cur_trial,]
        for(cur_blank in 1:numof_blanks){
            p = p + geom_rect(data=cur_blanks[cur_blank,], inherit.aes=F, 
                              aes(xmin=blank_time_beg, xmax=blank_time_end, 
                                  ymin=ylims[1], ymax=ylims[2]),
                              color='transparent', fill='black', alpha=.3)
        }
    }
    # Marking the cue:
    p = p + geom_rect(data=trials[as.character(cur_trial),], inherit.aes=F, 
                      aes(xmin=cue_time-trial_time_beg, xmax=cue_time+.2-trial_time_beg, 
                          ymin=ylims[1], ymax=ylims[2]),
                      color='transparent', fill='orange', alpha=.3)
    # Marking the shutter on/off times [2019-06-28] and space presses:
    if(cond == 'cond-a'){
        p = p + geom_rect(data=trials[as.character(cur_trial),], inherit.aes=F,
                          aes(xmin=shutter_time_beg-trial_time_beg,
                              xmax=shutter_time_end-trial_time_beg, ymin=ylims[1], ymax=ylims[2]),
                          color='transparent', fill='green', alpha=.3)
        p = p + geom_vline(data=trials[as.character(cur_trial),], 
                           aes(xintercept=cue_resp_time-trial_time_beg), linetype='dotted')
    }
    # Indicating the start of key monitoring & trial end times (optional/for debug):
    p = p + geom_vline(data=trials[as.character(cur_trial),], 
                       aes(xintercept=targ_time-trial_time_beg), linetype='longdash')
    # Marking the response time to the cue in the experimental conditions:
    if(cond != 'cond-m'){
        p = p + geom_vline(data=trials[as.character(cur_trial),], 
                           aes(xintercept=resp_time-trial_time_beg))
        p = p + geom_hline(aes(yintercept=-8), color='grey')
        p = p + geom_hline(aes(yintercept=8), color='grey')
    }
    p = p + labs(y='Gaze X-Position (deg)', x='Time (s)', colour='Parameter')
    p = p + ggtitle(paste('Trial #', as.character(cur_trial)))
    # Dealing with exclusions:
    if(cur_trial %in% all_excl$trial){
        excl_cells = all_excl$excl_reason[all_excl$trial==cur_trial]
        if(length(excl_cells) == 1){
            excl_reasons = excl_cells
        } else {
            excl_reasons = excl_cells[1]
            for(cur_excl_cell in 1:(length(excl_cells)-1)){
                excl_reasons = paste(excl_reasons, excl_cells[cur_excl_cell+1], sep='\n')
            }
        }
        excl_text = paste('excluded:', excl_reasons, sep='\n')
        p = p + annotate('text', x=0, y=9, label=excl_text, colour='red', hjust=0, vjust=1)
    } else{
        p = p + theme_bw()
    }
    p = p + theme(axis.title.y.right = element_text(color='magenta'),
              axis.title.y.left = element_text(color='blue'),
              legend.position = 'none')
    suppressWarnings(print(p))
}
