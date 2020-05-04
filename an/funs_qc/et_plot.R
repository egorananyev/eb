
# Visualize the data:
library(ggplot2)
plot_qc = function(trial_samples, trial_blanks, this_trial, trial_sacc){
    ylims = c(-10, 10)
    p = ggplot(trial_samples, aes(x=time))
    p = p + geom_line(aes(y=xr, colour='Gaze'))
    # p = p + geom_line(aes(y=yr-432, colour='Gaze'), linetype='longdash')
    p = p + geom_line(aes(y=((psr-mean(psr))/100), colour='Pupil'))  # colour='Pupil'))
    p = p + scale_y_continuous(sec.axis = sec_axis(~.*5, name='Pupil Diameter (a.u.)'),limits=ylims)
    # p = p + xlim(xlims)
    p = p + scale_colour_manual(values=c('blue', 'chocolate'))
    
    # Assigning the number of blanks to be processed:
    if(nrow(trial_blanks) == 1){  # if there's a single row in <blanks>
        if(is.na(trial_blanks[1,1])){
            numof_blanks = 0  # if there is a blank table
        } else {
            numof_blanks = 1
        }
    } else {
        numof_blanks = nrow(trial_blanks)  # number of blanks in this trial:
    }
    
    # Marking blanks:
    if(!is.na(numof_blanks) & numof_blanks >= 1){
        for(cur_blank in 1:numof_blanks){
            this_blank_df = trial_blanks[cur_blank, ]
            this_alpha = .35
            if(!this_blank_df$blank_post_cue){this_alpha = .65}
            if(this_blank_df$blank_post_targ){this_alpha = .25}
            p = p + geom_rect(data=this_blank_df, inherit.aes=F, 
                              aes(xmin=blank_time_beg, xmax=blank_time_end, 
                                  ymin=ylims[1], ymax=ylims[2]),
                              color='transparent', fill='black', alpha=this_alpha)
        }
    }
    # Marking saccades:
    numof_sacc = nrow(trial_sacc)
    if(!is.na(numof_sacc) & numof_sacc>=1){
        for(cur_sacc in 1:numof_sacc){
            this_sacc_df = trial_sacc[cur_sacc,]
            if(this_sacc_df$bef_targ){this_alpha=.35} else {this_alpha=.15}
            p = p + geom_rect(data=this_sacc_df, inherit.aes=F, 
                              aes(xmin=sacc_beg_time, xmax=sacc_end_time, 
                                  ymin=ylims[1], ymax=ylims[2]),
                              color='transparent', fill='blue', alpha=this_alpha)
        }
    }
    # Marking the cue:
    p = p + geom_rect(data=this_trial, inherit.aes=F, 
                      aes(xmin=cue_time-trial_time_beg, xmax=cue_time+.1-trial_time_beg, 
                          ymin=ylims[1], ymax=ylims[2]),
                      color='transparent', fill='orange', alpha=.3)
    # Marking the shutter on/off times:
    if(cond == 'cond-a'){
        p = p + geom_rect(data=this_trial, inherit.aes=F,
                          aes(xmin=shutter_time_beg-trial_time_beg,
                              xmax=shutter_time_end-trial_time_beg, ymin=ylims[1], ymax=ylims[2]),
                          color='transparent', fill='green', alpha=.3)
    }
    # [2019-06-28] and space presses AB condition:
    # [2020-04-07] unfortunately, the space bar is recorded in beh file, but not in et ds, for NB
    if(cond %in% c('cond-a')){
        p = p + geom_vline(data=this_trial, aes(xintercept=cue_resp_time-trial_time_beg),
                           linetype='dotted')
    }
    # Indicating the start of key monitoring & trial end times (optional/for debug):
    p = p + geom_vline(data=this_trial, aes(xintercept=targ_time-trial_time_beg), linetype='longdash')
    # Marking the response time to the cue in the experimental conditions:
    if(cond != 'cond-m'){
        p = p + geom_vline(data=this_trial, aes(xintercept=resp_time-trial_time_beg))
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
    p = p + theme(axis.title.y.right = element_text(color='chocolate'),
              axis.title.y.left = element_text(color='blue'),
              legend.position = 'none')
    suppressWarnings(print(p))
}
