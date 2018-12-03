####!/usr/bin/arch -i386 /usr/bin/python # -*- coding: utf-8 -*-
from __future__ import print_function

"""
The influence of eye blinks on attentional cueing.
Creator: Egor Ananyev
Original date: 2018-11-15
"""

# todo: test if I need to use iohub for keyboard, or can get away with default psychopy

from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import visual, core, data, event, gui, sound
from psychopy.constants import *  # things like STARTED, FINISHED
from psychopy.iohub import (EventConstants, EyeTrackerConstants, getCurrentDateTimeString,
                            ioHubExperimentRuntime)
import numpy as np
import pandas as pd
from datetime import datetime
import shutil
import copy

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))

# ====================================================================================
## Initial variables.
# Experiment variables:
exp_name = 'eb1'
# General timing variables:
fixation_dur = .5  # the time frame before the cue, if any
iti_dur_min = .3
iti_dur_max = .7
# Display dimensions:
fix_cross_sz = .3  # length of the two cross lines, in dva
fix_cross_thick = 3  # pixels
ds = 58  # distance to screen in cm
dr = (1152, 864)  # display resolution in px
dd = (29.5, 16.6)  # display dimensions in cm
# TODO update the display dimensions once known
# Cue:
cue_off_y = .8
cue_dur = .25
cue_valid = .7  # validity of the cue
beep_dur = .1
# Shutter:
shutter_dur_min = .25
shutter_dur_max = .45
# Target:
targ_off_x = 3.5
targ_diam = .5

# ====================================================================================
# Converter functions:
def cm2px(cm, dr_=dr, dd_=dd):
    px = int(cm * (dr_[0] / dd_[0]))
    return px

def px2cm(px, dr_=dr, dd_=dd):
    cm = px / (dr_[0] / dd_[0])
    return cm

def cm2dg(cm, ds_=ds):
    dg = np.degrees(np.arctan(cm / ds_))
    return dg

def dg2cm(dg, ds_=ds):
    cm = ds_ * np.tan(np.radians(dg))
    return cm

def px2dg(px, cm2dg_=cm2dg, px2cm_=px2cm):
    dg = cm2dg_(px2cm_(px))
    return dg

def dg2px(dg, cm2px_=cm2px, dg2cm_=dg2cm):
    px = int(cm2px_(dg2cm_(dg)))
    return px

# ====================================================================================
# Converting win dimensions to pixels
fix_cross_sz = dg2px(fix_cross_sz)
cue_off_y = dg2px(cue_off_y)
targ_diam = dg2px(targ_diam)

# ====================================================================================
## Getting user info about the experiment session:
exp_info = {u'experiment': exp_name, u'exp_run': '0', u'participant': u'',
            u'cond': u'C', u'session': u''}
exp_name = exp_info['experiment']
dlg = gui.DlgFromDict(dictionary=exp_info, title=exp_name)  # dialogue box
if not dlg.OK:
    core.quit()  # user pressed cancel
exp_info['time'] = datetime.now().strftime('%Y-%m-%d_%H%M')
end_exp_now = False  # flag for 'escape' or other condition => quit the exp

# ====================================================================================
## Setting up the display window:
win = visual.Window(size=dr, fullscr=False, screen=1)  #, allowGUI=False,
                    # allowStencil=False, color='grey', blendMode='avg', useFBO=True,
                    # units='pix')
frame_rate = win.getActualFrameRate()
print('frame_rate=' + str(frame_rate))

# ====================================================================================
# Output files:
file_name = '%s_ER%s_P%s_%s_sess%s_%s' % (exp_name, exp_info['exp_run'], exp_info['participant'],
                                          exp_info['cond'], exp_info['session'], exp_info['time'])
file_path = '..' + os.sep + 'data' + os.sep + file_name
print('file path is ' + file_path)

# Condition-related file paths (input):
if exp_info['exp_run'] == '1':  # if this is an experimental run (not a trial run):
    conditionsFilePath = 'cond-files' + os.sep + 'cond-' + exp_name + '.csv'
else:  # otherwise, setting the condition file to the trial run:
    conditionsFilePath = 'cond-files' + os.sep + 'cond_' + exp_name + '_tr' + '.csv'
print(conditionsFilePath)
os.chdir(_thisDir)

# ====================================================================================
## Initialize the stimuli and instructions
trial_clock = core.Clock()
# TODO Do I need this ISI period definition?
isi = core.StaticPeriod(win=win, screenHz=frame_rate, name='isi')
fix_cross = visual.ShapeStim(win, vertices=((0, -fix_cross_sz), (0, fix_cross_sz), (0, 0),
                                            (-fix_cross_sz, 0), (fix_cross_sz, 0)),
                             lineWidth=2, closeShape=False, lineColor='white')

# Auditory tone (for blink condition):
print('Using %s (with %s) for sounds' % (sound.audioLib, sound.audioDriver))
beep = sound.Sound('C', octave=4, sampleRate=44100, secs=beep_dur, stereo=True)
beep.setVolume(.1)

# Cue:
arrow_vert = [(.5, 0), (0, .3), (0, .1), (-.5, .1), (-.5, -.1), (0, -.1), (0, -.3)]
cue_arrow = visual.ShapeStim(win, vertices=arrow_vert, fillColor='black', size=30, lineColor='black',
                             pos=(0, cue_off_y))

# Target:
targ = visual.Circle(win, radius=targ_diam/2, edges=32, pos=(50, 0), fillColor='white')

# Target response:
key_arrow = event.BuilderKeyResponse()  # create an object of type keyresponse
# todo: we need multiple instruction screens in the TR (not ER)
# instr_text = visual.TextStim(win, text='press any key\n     to start', font='Cambria',
#                              height=int(dg2px(.5)), wrapWidth=int(dg2px(4)),
#                              color='white', alignHoriz='center')

# ====================================================================================
# Initiating the trial loop

nTrialsDone = 0

while nTrialsDone < 1:
    nTrialsDone += 1
    t = 0
    trial_clock.reset()  # clock
    frame_n = -1
    targ_resp_given = False

    # update component parameters for each repeat
    # keep track of which components have finished
    trialComponents = [fix_cross, key_arrow, isi, cue_arrow]
    for thisComponent in trialComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED

    # -------Start Routine "trial"-------
    continue_routine = True
    while continue_routine:
        t_prev = t
        t = trial_clock.getTime()  # get current time
        frame_n = frame_n + 1  # number of completed frames (0 is the first frame)
        print('frame_n=' + str(frame_n) + ' took ' + str(t - t_prev))
        beep.play()
        fix_cross.draw()
        cue_arrow.draw()
        targ.draw()

        # print('time=' + str(trial_clock.getTime()))
        # *key_arrow* updates for target responses:
        if key_arrow.status == NOT_STARTED:  # record key strokes only after trial end
            # keep track of start time/frame for later
            key_arrow.tStart = t  # underestimates by a little under one frame
            key_arrow.frameNStart = frame_n  # exact frame index
            key_arrow.status = STARTED
            # keyboard checking is just starting
            key_arrow.clock.reset()  # now t=0
            event.clearEvents(eventType='keyboard')
            #kb_device.clearEvents()

        # registering response at the end of the trial
        if key_arrow.status == STARTED:
            theseKeys = event.getKeys(keyList=['left', 'right'])
            if len(theseKeys) > 0:
                if 'left' in theseKeys:
                    print('response: left')
                    behRespTrial = -1  # the first number is negative
                    targ_resp_given = True
                elif 'right' in theseKeys:
                    print('response: right')
                    behRespTrial = 1  # the second number is positive
                    targ_resp_given = True
                # if targ_resp_given:  # this is overwritten every time any key is pressed
                #     rt = t - targTon
                #     if behRespTrial == thisTargXoff:
                #         corrResp = 1  # correct dir resp
                #     else:
                #         corrResp = 0  # this is the result from both 0 and wrong dir resp
                #     if thisContr <= -2:
                #         corrResp = 0

        # pause text and data exporting
        if targ_resp_given:
            # thisStair.addResponse(corrResp)
            print('finished trial')
            continue_routine = False

        # # *isi* period
        # if isi.status == NOT_STARTED:
        #     # keep track of start time/frame for later
        #     isi.tStart = t  # underestimates by a little under one frame
        #     isi.frameNStart = frame_n  # exact frame index
        #     # isi.start(isi_duration)
        # # one frame should pass before updating params and completing
        # elif isi.status == STARTED:
        #     isi.complete()  # finish the static period
        #     continue_routine = False

        # check if all components have finished
        # a component has requested a forced-end of Routine:
        if not continue_routine:
            break
        # will revert to True if at least one component still running
        continue_routine = False
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "status") and \
                    thisComponent.status != FINISHED:
                continue_routine = True
                break  # at least one component has not yet finished

        # check for quit (the Esc key)
        if end_exp_now or event.getKeys(keyList=["escape"]):
            core.quit()

        print('time=' + str(trial_clock.getTime()))
        # refresh the screen
        # don't flip if this routine is over or we'll get a blank screen
        if continue_routine:
            flip_time = win.flip()
            print('flip_time=' + str(flip_time))
        print('time=' + str(trial_clock.getTime()))

    # -------Ending Routine "trial"-------
    for thisComponent in trialComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)

    # thisExp.nextEntry()

print("finished the experiment")

win.close()
core.quit()

if __name__ == "__main__":
    import os
    from psychopy.iohub import module_directory


    def main(configurationDirectory):
        #-- Sol's comments:
        """
        Creates an instance of the ExperimentRuntime class, gets the eye tracker
        the user wants to use for the demo, and launches the experiment logic.
        """
        # The following code merges a iohub_config file called iohub_config.yaml.part,
        # that has all the iohub_config settings, other than those for the eye tracker, with
        # the eye tracker configs in the yaml files in the eyetracker_configs dir. {...}
        #
        # The merged result is saved as iohub_config.yaml so it can be picked up
        # by the Experiment _runtime as normal.
        #--
        eye_tracker_config_file = 'eyetracker_configs/eyelink_config.yaml'
        base_config_file = os.path.normcase(os.path.join(configurationDirectory,
                                                         'eyetracker_configs/iohub_config.yaml.part'))

        eyetrack_config_file = os.path.normcase(os.path.join(configurationDirectory, 'eyetracker_configs',
                                                             eye_tracker_config_file))

        combined_config_file_name = os.path.normcase(os.path.join(configurationDirectory,
                                                                  'eyetracker_configs/iohub_config.yaml'))

        ExperimentRuntime.mergeConfigurationFiles(base_config_file, eyetrack_config_file,
                                                  combined_config_file_name)

        runtime = ExperimentRuntime(configurationDirectory, "experiment_config.yaml")
        runtime.start(('SR Research EyeLink',))


    #-- Get the current directory, using a method that does not rely on __FILE__
    # or the accuracy of the value of __FILE__.
    configurationDirectory = module_directory(main)

    #-- Run the main function, which starts the experiment runtime
    main(configurationDirectory)