####!/usr/bin/arch -i386 /usr/bin/python # -*- coding: utf-8 -*-
from __future__ import print_function

"""
The influence of eye blinks on attentional cueing.
Creator: Egor Ananyev
Original date: 2018-11-15
"""

# todo: test if I need to use iohub for keyboard, or can get away with default psychopy

from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import visual, core, data, event, gui, sound, monitors
# from psychopy.constants import *  # things like STARTED, FINISHED
from psychopy.iohub import (EventConstants, EyeTrackerConstants, getCurrentDateTimeString,
                            ioHubExperimentRuntime)
import numpy as np
from psychopy.data import TrialHandler, importConditions
# import pandas as pd
from datetime import datetime
# import shutil
# import copy
import os


# Ensure that relative paths start from the same directory as this script
# _thisDir = os.path.dirname(os.path.abspath(__file__))

class ExperimentRuntime(ioHubExperimentRuntime):
    # -
    """
    Create an experiment using psychopy and the ioHub framework by extending the ioHubExperimentRuntime class.
    At minimum all that is needed in the __init__ for the new class, here called ExperimentRuntime, is the a
    call to the ioHubExperimentRuntime __init__ itself.
    """

    def run(self, *args):

        # ====================================================================================
        ## Initial variables.
        # experiment variables:
        exp_name = 'eb1'
        debug = True  # are we using the eye tracker in the session? must be False for expt
        # general timing variables:
        fix_1_dur = .4  # the time frame before the cue, if any
        fix_2_min = 240 # these are in ms, because we need a random _integer_ in this range
        fix_2_max = 500
        # the time window for the blink - quite conservative - should include the whole blink,
        # but is independent of the blink start/end
        blink_time_window = .3
        post_blink_dur = .05
        # display dimensions:
        fix_cross_sz = .3  # length of the two cross lines, in dva
        fix_cross_thick = 3  # pixels
        if debug:
            dr = (576, 432)
            ds = 61
            dd = (17.2, 12.9)  # display dimensions in cm
        else:
            dr = (1152, 864)  # display resolution in px
            ds = 58  # distance to screen in cm
            dd = (34.4, 25.8)  # display dimensions in cm
        # cue:
        cue_off_y = .8
        cue_dur = .2
        cue_valid = .7  # validity of the cue
        beep_dur = .05
        # target:
        targ_off_x = 3.5
        targ_diam = .5

        # ====================================================================================
        ## Conversion functions:

        def cm2px(cm, dr_=dr, dd_=dd):
            px = int(cm * (dr_[0] / dd_[0]))
            return px

        def dg2cm(dg, ds_=ds):
            cm = ds_ * np.tan(np.radians(dg))
            return cm

        def dg2px(dg, cm2px_=cm2px, dg2cm_=dg2cm):
            px = int(cm2px_(dg2cm_(dg)))
            return px

        # Converting stimulus dimensions to pixels


        # ====================================================================================
        ## getting user info about the experiment session:
        exp_info = {u'expt': exp_name, u'subj': u'', u'cond': u'c', u'sess': u''}
        # conditions: 't'=training, 'c'=control, 'a'=artificial blink, 'v'=voluntary blink
        exp_name = exp_info['expt']
        dlg = gui.DlgFromDict(dictionary=exp_info, title=exp_name)  # dialogue box
        if not dlg.OK:
            core.quit()  # user pressed cancel
        exp_info['time'] = datetime.now().strftime('%Y-%m-%d_%H%M')
        end_exp_now = False  # flag for 'escape' or other condition => quit the exp

        # ====================================================================================
        ## Input and output

        # condition file:
        exp_conditions = importConditions('cond-files/cond_' + exp_name + '_' +
                                          exp_info['cond'] + '.xlsx')
        trials = TrialHandler(exp_conditions, 1)

        # - Inform the ioDataStore that the experiment is using a TrialHandler. The ioDataStore will
        # create a table which can be used to record the actual trial variable values (DV or IV)
        # in the order run / collected:
        self.hub.createTrialHandlerRecordTable(trials)

        # output file:
        out_file_name = '%s_subj-%s_cond-%s_sess-%s_%s' % (exp_name, exp_info['subj'], exp_info['cond'],
                                                           exp_info['sess'], exp_info['time'])
        out_file_path = '..' + os.sep + 'data' + os.sep + out_file_name
        print('output file path is ' + out_file_path)

        # ====================================================================================
        ## Handy shortcuts:
        tracker = self.hub.devices.tracker
        display = self.hub.devices.display
        kb = self.hub.devices.keyboard
        mouse = self.hub.devices.mouse

        # ====================================================================================
        ## EyeLink & screen setup

        # - Start by running the eye tracker default setup procedure.
        selected_eyetracker_name = args[0]
        if not debug:
            tracker.runSetupProcedure()

        if debug:
            mon = monitors.Monitor('Dell', width=dd[0], distance=ds)
            window = visual.Window(size=dr, monitor=mon, fullscr=True, screen=1, units='deg')
        else:
            # - Create a psychopy window, full screen resolution, full screen mode...
            res = display.getPixelResolution()
            window = visual.Window(res, monitor='Screen1', units=display.getCoordinateType(),
                                   fullscr=True, allowGUI=False, waitBlanking=False,
                                   screen=display.getIndex())

        display_coord_type = display.getCoordinateType()
        print('unit type: ', display_coord_type)
        frame_rate = window.getActualFrameRate()
        # print('frame rate: ' + str(frame_rate))

        # ====================================================================================
        ## Initialize the stimuli and instructions

        instruction_text = "press space key to start the experiment"
        # instructions_text_stim = visual.TextStim(window, text=instruction_text, pos=[0, 0], height=24,
        #                                          color=[-1, -1, -1], colorSpace='rgb', alignHoriz='center',
        #                                          alignVert='center', wrapWidth=window.size[0] * .9)
        instructions_text_stim = visual.TextStim(window, text=instruction_text, height=.4)
        # fix_cross = visual.shapestim(window, vertices=((0, -fix_cross_sz), (0, fix_cross_sz), (0, 0),
        #                                             (-fix_cross_sz, 0), (fix_cross_sz, 0)),
        #                              linewidth=2, closeshape=False, linecolor='white')
        # todo: potentially adjust the 'height' for better cross size
        fix_cross = visual.TextStim(window, text='+', bold='True', pos=[0, 0], rgb=1, height=.3)

        # auditory tone (for blink condition):
        print('using %s (with %s) for sounds' % (sound.audioLib, sound.audioDriver))
        beep = sound.Sound('A', octave=5, sampleRate=44100, secs=beep_dur, stereo=True)
        # todo: test if this is too soft on the machine
        beep.setVolume(.1)  # this isn't needed in Aaron's paradigm -- the sound is softer

        # cue:
        arrow_vert = [(.5, 0), (0, .3), (0, .1), (-.5, .1), (-.5, -.1), (0, -.1), (0, -.3)]
        cue_arrow = visual.ShapeStim(window, vertices=arrow_vert, fillColor='black', size=.3,
                                     lineColor='black', pos=(0, cue_off_y))

        # target:
        targ = visual.Circle(window, radius=targ_diam / 2, edges=32, pos=(50, 0), fillColor='white')

        # target response:
        key_arrow = event.BuilderKeyResponse()

        # ====================================================================================
        ## Initiating the iohub routines

        # wait until a key event occurs after the instructions are displayed
        self.hub.clearEvents('all')

        # - Send some information to the ioHub DataStore as experiment messages
        # including the eye tracker being used for this session.
        self.hub.sendMessageEvent(text="IO_HUB EXPERIMENT_INFO START")
        self.hub.sendMessageEvent(text="ioHub Experiment started {0}".format(getCurrentDateTimeString()))
        self.hub.sendMessageEvent(
            text="Experiment ID: {0}, Session ID: {1}".format(self.hub.experimentID, self.hub.experimentSessionID))
        self.hub.sendMessageEvent(
            text="Stimulus Screen ID: {0}, Size (pixels): {1}, CoordType: {2}".format(display.getIndex(),
                                                                                      display.getPixelResolution(),
                                                                                      display.getCoordinateType()))
        self.hub.sendMessageEvent(
            text="Calculated Pixels Per Degree: {0} x, {1} y".format(*display.getPixelsPerDegree()))
        self.hub.sendMessageEvent(text="Eye Tracker being Used: {0}".format(selected_eyetracker_name))
        self.hub.sendMessageEvent(text="IO_HUB EXPERIMENT_INFO END")

        self.hub.clearEvents('all')
        print('initiated the iohub routines')

        # ====================================================================================
        ## Handy routines:

        # Frame-skipping check:
        def frame_skip_check(elapsed_t, elapsed_frames):
            # The number of elapsed frames should match the time:
            print('time=%.3f  frames=%d  rate=%.4f' % (elapsed_t, elapsed_frames,
                                                       (elapsed_t / elapsed_frames)))

        # This is done at every frame update, regardless of trial phase, so predefining:
        def frame_routine():
            window.flip()
            fix_cross.draw()
            # Checking for quit (the Esc key)
            if event.getKeys(keyList=['escape']):
                exit_routine()
            # Measuring time elapsed since the start of the trial:
            trial_elapsed_t_ = trial_clock.getTime() - trial_t_start
            return trial_elapsed_t_

        # Also no variation across frames, but only available upon call, which is made only in key
        # registering phase.
        def exit_routine():
            window.close()
            core.quit()

        # ====================================================================================
        # initiating the trial loop

        n_trials_done = 0
        trial_clock = core.Clock()
        if debug:
            it_clock = core.Clock()  # inter-trial clock, for measuring unintended delays

        for trial in trials:

            # ================================================================================
            ## First trial initiates instructions and sends the expt initiation message to the
            # eye tracker:
            if n_trials_done == 0:
                # - Update the instruction screen text...
                instructions_text_stim.setText(instruction_text)
                instructions_text_stim.draw()
                flip_time = window.flip()
                self.hub.sendMessageEvent(text="EXPERIMENT_START", sec_time=flip_time)
                # wait until a space key event occurs after the instructions are displayed
                kb.waitForPresses(keys=' ')

            # ================================================================================
            ## Trial components:

            # Trial components pertaining to time, frames, and trial number:
            n_trials_done += 1
            print('======TRIAL#' + str(n_trials_done) + '======')
            trial_clock.reset()  # clock

            if debug:
                if n_trials_done > 1:
                    print(n_trials_done)
                    iti_dur = it_clock.getTime() - iti_end_trial
                    print('inter-trial duration: %.3f' % iti_dur)

            # Randomize the duration of the post-cue fixation:
            fix_2_dur = np.random.randint(fix_2_min, fix_2_max)/1000  # randomizing & converting to sec
            if debug:
                print('fix2dur = %.3f' % fix_2_dur)

            # Randomize target location:
            this_targ_loc = np.random.randint(0, 1) * 2 - 1
            print('this target location is ' + str(this_targ_loc))

            # Trial components pertaining to eye blinks:
            blinked = False  # blink detection
            t_blink_start = 0  # blink start & end time
            t_blink_end = 0


            # ================================================================================
            ## Starting the eye-tracking recording:
            flip_time = window.flip()
            trial['session_id'] = self.hub.getSessionID()
            trial['trial_id'] = n_trials_done
            trial['TRIAL_START'] = flip_time
            self.hub.sendMessageEvent(text="TRIAL_START", sec_time=flip_time)
            self.hub.clearEvents('all')
            if not debug:
                tracker.setRecordingState(True)

            # ================================================================================
            ## Starting the frame cycle & waiting for the response:
            trial_t_start = trial_clock.getTime()
            if debug:
                trial_elapsed_frames = 0  # counting frames for frame skip test

            # ================================================================================
            ## Fixation cross:
            fix_1_frames = int(fix_1_dur*frame_rate)
            for fix_1_frame in range(fix_1_frames):
                frame_routine()
                if debug:
                    trial_elapsed_frames += 1

            # ================================================================================
            ## The brief period with beep & cue:
            beep_cue_frames = int(beep_dur*frame_rate)
            for beep_cue_frame in range(beep_cue_frames):
                frame_routine()
                beep.play()
                cue_arrow.draw()
                if debug:
                    trial_elapsed_frames += 1

            # ================================================================================
            ## The rest of the period without the beep, but with the cue:
            cue_only_frames = int((cue_dur-beep_dur)*frame_rate)
            for cue_only_frame in range(cue_only_frames):
                frame_routine()
                cue_arrow.draw()
                if debug:
                    trial_elapsed_frames += 1

            # ================================================================================
            ## Fixation 2 + blink period, i.e., the fixation period after the cue:
            blink_time_period_frames = int((fix_2_dur+blink_time_window)*frame_rate)
            for blink_time_period_frame in range(blink_time_period_frames):
                frame_routine()
                if debug:
                    trial_elapsed_frames += 1

            # ================================================================================
            ## Behavioural response: measuring the reaction time:

            # Trial components pertaining to behavioural response:
            targ_resp_given = False
            rt_start = trial_clock.getTime()

            # Displaying the target and measuring the reaction time.
            while not targ_resp_given:

                # Measuring the time it takes for the behavioural response:
                trial_elapsed_t = frame_routine()

                # Drawing the target:
                targ.draw()

                if debug:
                    trial_elapsed_frames += 1

                # Monitoring for key presses:
                arrow_keys = event.getKeys(keyList=['left', 'right'])
                if len(arrow_keys) > 0:
                    if 'left' in arrow_keys:
                        print('response: left')
                        beh_resp = -1
                        targ_resp_given = True
                    elif 'right' in arrow_keys:
                        print('response: right')
                        beh_resp = 1
                        targ_resp_given = True
                    if targ_resp_given:  # this is overwritten every time any key is pressed
                        rt = trial_clock.getTime() - rt_start
                        if beh_resp == this_targ_loc:
                            corr_resp = 1  # correct location response
                        else:
                            corr_resp = 0  # incorrect location response
                        print('RT=%.2f correct?=%s' % (rt, corr_resp))
                        if debug:  # in debug mode, check if the frame rate looks okay
                            frame_skip_check(trial_elapsed_t, trial_elapsed_frames)
                            iti_end_trial = it_clock.getTime()

if __name__ == "__main__":
    import os
    from psychopy.iohub import module_directory

    def main(base_dir):
        # -- Sol's comments:
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
        # --

        ##
        # I revamped the code so that all iohub config files are stored in the 'iohub_configs'
        # subdirectory. All the .yaml files should be stored there.
        #
        # Also, since a single eye-tracker is used, I removed all references to other
        # eye trackers beside the EyeLink and removed the gui prompt asking for the type of eye
        # tracker.

        eye_config_dir = os.path.normcase(os.path.join(base_dir, 'iohub_configs'))
        base_config_file = os.path.normcase(os.path.join(eye_config_dir, 'iohub_config.yaml.part'))
        eyetracker_config_file = os.path.normcase(os.path.join(eye_config_dir, 'eyelink_config.yaml'))
        combined_config_file_name = os.path.normcase(os.path.join(eye_config_dir, 'iohub_config.yaml'))
        ExperimentRuntime.mergeConfigurationFiles(base_config_file, eyetracker_config_file,
                                                  combined_config_file_name)
        runtime = ExperimentRuntime(eye_config_dir, 'experiment_config.yaml')
        runtime.start(('SR Research EyeLink',))


    # -- Get the current directory, using a method that does not rely on __FILE__
    # or the accuracy of the value of __FILE__.
    base_dir = module_directory(main)

    # -- Run the main function, which starts the experiment runtime
    main(base_dir)
