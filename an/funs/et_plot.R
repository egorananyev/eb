
# Visualize the data:
library(ggplot2)
plot_qc = function(samples, trials, xlims){
    ylims = c(0, 1000)
    p = ggplot(samples[samples$trial==cur_trial,], aes(x=time))
    p = p + geom_line(aes(y=yr, colour='Gaze'))
    p = p + geom_line(aes(y=psr/5, colour='Pupil'))
    p = p + scale_y_continuous(sec.axis = sec_axis(~.*5, name='Pupil Diameter (a.u.)'),
                               limits=ylims)
    p = p + xlim(xlims)
    p = p + scale_colour_manual(values=c('blue', 'red'))
    p = p + theme(legend.position = c(0.15, 0.2))
    # number of blanks in this trial:
    numof_blanks = sum(blanks$trial==cur_trial)
    if(numof_blanks >= 1){
        cur_blanks = blanks[blanks$trial==cur_trial,]
        for(cur_blank in 1:numof_blanks){
            p = p + geom_rect(data=cur_blanks[cur_blank,], inherit.aes=F, 
                              aes(xmin=blank_time_beg, xmax=blank_time_end, 
                                  ymin=ylims[1], ymax=ylims[2]),
                              color='transparent', fill='purple', alpha=.3)
        }
    }
    # Marking the cue:
    p = p + geom_rect(data=trials[cur_trial,], inherit.aes=F, 
                      aes(xmin=cue_time-trial_time_beg, xmax=cue_time+.2-trial_time_beg, 
                          ymin=ylims[1], ymax=ylims[2]),
                      color='transparent', fill='orange', alpha=.3)
    # Marking the supposed blink time window:
    p = p + geom_rect(data=trials[cur_trial,], inherit.aes=F, 
                      aes(xmin=blink_window_time-trial_time_beg-.26, 
                          xmax=blink_window_time+.3-trial_time_beg, 
                          ymin=ylims[1], ymax=ylims[2]),
                      color='transparent', fill='green', alpha=.3)
    # Indicating the start of key monitoring & trial end times (optional/for debug):
    # p = p + geom_vline(data=trials[cur_trial,], aes(xintercept=resp_time-trial_time_beg))
    # p = p + geom_vline(data=trials[cur_trial,], aes(xintercept=trial_time_end-trial_time_beg))
    p = p + labs(y='Gaze Y-Position', x='Time', colour='Parameter')
    p = p + ggtitle(paste('Trial #', as.character(cur_trial))) + theme_bw()
    print(p)
}
