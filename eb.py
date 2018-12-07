####!/usr/bin/arch -i386 /usr/bin/python # -*- coding: utf-8 -*-
from __future__ import print_function

"""
The influence of eye blinks on attentional cueing.
Creator: Egor Ananyev
Original date: 2018-11-15
"""

# todo: test if I need to use iohub for keyboard, or can get away with default psychopy

from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import visual, core, event, gui, sound, monitors
from psychopy.iohub import (EventConstants, EyeTrackerConstants, getCurrentDateTimeString,
                            ioHubExperimentRuntime)
import numpy as np
from psychopy.data import TrialHandler, importConditions
import pandas as pd
from datetime import datetime


class ExperimentRuntime(ioHubExperimentRuntime):
    # -
    """
    Create an experiment using psychopy and the ioHub framework by extending the ioHubExperimentRuntime class.
    At minimum all that is needed in the __init__ for the new class, here called ExperimentRuntime, is the a
    call to the ioHubExperimentRuntime __init__ itself.
    """

    def run(self, *args):

        ## Initial variables.
        # experiment variables:
        exp_name = 'eb1'
        trial_n = 15  # trials per row; 15 gives 300 trials
        debug = True  # are we using the eye tracker in the session? must be False for expt general timing variables:
        fix_1_dur = .4  # the time frame before the cue, if any
        blink_latency_min = 240  # these are in ms, because we need a random _integer_ in this range
        blink_latency_max = 500
        # the time window for the blink - quite conservative - should include the whole blink, but is independent of
        # the blink start/end
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
        beep_dur = .05
        # target:
        targ_off_x = 3
        targ_diam = .5

        ## getting user info about the experiment session:
        exp_info = {u'expt': exp_name, u'subj': u'', u'cond': u'd', u'sess': u''}
        # conditions: 't'=training, 'c'=control, 'a'=artificial blink, 'v'=voluntary blink, 'd'=debug
        exp_name = exp_info['expt']
        dlg = gui.DlgFromDict(dictionary=exp_info, title=exp_name)  # dialogue box
        if not dlg.OK:
            core.quit()  # user pressed cancel
        exp_info['time'] = datetime.now().strftime('%Y-%m-%d_%H%M')
        if exp_info['cond'] == 'd':
            trial_n = 1
        end_exp_now = False  # flag for 'escape' or other condition => quit the exp

        ## Input and output

        # condition file:
        exp_conditions = importConditions('cond-files/cond_' + exp_name + '_' + exp_info['cond'] + '.xlsx')
        trials = TrialHandler(exp_conditions, trial_n, extraInfo=exp_info)

        # output file:
        out_file_name = '%s_subj-%s_cond-%s_sess-%s_%s' % (exp_name, exp_info['subj'], exp_info['cond'],
                                                           exp_info['sess'], exp_info['time'])
        out_file_path = '..' + os.sep + 'data' + os.sep + out_file_name

        # output matrix:
        output_mat = {}

        ## Handy shortcuts:
        tracker = self.hub.devices.tracker
        display = self.hub.devices.display
        kb = self.hub.devices.keyboard

        ## EyeLink & screen setup

        # - Start by running the eye tracker default setup procedure.
        selected_eyetracker_name = args[0]
        if not debug:
            tracker.runSetupProcedure()

        if debug:
            mon = monitors.Monitor('Dell', width=dd[0], distance=ds)
            window = visual.Window(size=dr, monitor=mon, fullscr=False, screen=1, units='deg')
        else:
            # - Create a psychopy window, full screen resolution, full screen mode...
            res = display.getPixelResolution()
            window = visual.Window(res, monitor='Screen1', units=display.getCoordinateType(), fullscr=True,
                                   allowGUI=False, waitBlanking=False, screen=display.getIndex())

        # display_coord_type = display.getCoordinateType()
        # print('unit type: ', display_coord_type)
        if debug:
            frame_rate = 60
        else:
            frame_rate = window.getActualFrameRate()
            print('frame rate: ' + str(frame_rate))
            if frame_rate < 100:
                print('WARNING! The measured frame rate is lower than expected')

        ## Initialize the stimuli and instructions
        instruction_text = "Please press the ''space'' key\nto start the experiment"
        # instructions_text_stim = visual.TextStim(window, text=instruction_text, pos=[0, 0], height=24,
        #                                          color=[-1, -1, -1], colorSpace='rgb', alignHoriz='center',
        #                                          alignVert='center', wrapWidth=window.size[0] * .9)
        instructions_text_stim = visual.TextStim(window, text=instruction_text, height=.4)
        # todo: potentially adjust the 'height' for better cross size
        fix_cross = visual.TextStim(window, text='+', bold='True', pos=[0, 0], rgb=1, height=.3)

        # auditory tone (for blink condition):
        print('using %s (with %s) for sounds' % (sound.audioLib, sound.audioDriver))
        beep = sound.Sound('A', octave=5, sampleRate=44100, secs=beep_dur, stereo=True)
        # todo: test if this is too soft on the machine
        beep.setVolume(.1)  # this isn't needed in Aaron's paradigm -- the sound is softer

        # cue:
        arrow_vert = [(.5, 0), (0, .3), (0, .1), (-.5, .1), (-.5, -.1), (0, -.1), (0, -.3)]
        cue_arrow = visual.ShapeStim(window, vertices=arrow_vert, fillColor='black', size=.3, lineColor='black',
                                     pos=(0, cue_off_y))

        # target:
        targ = visual.Circle(window, radius=targ_diam / 2, edges=32, pos=(10, 0), fillColor='white')

        # target response:
        key_arrow = event.BuilderKeyResponse()

        ## Initiating the iohub routines

        # wait until a key event occurs after the instructions are displayed
        self.hub.clearEvents('all')

        # - Send some information to the ioHub DataStore as experiment messages
        # including the eye tracker being used for this session.
        self.hub.sendMessageEvent(text="IO_HUB EXPERIMENT_INFO START")
        self.hub.sendMessageEvent(text="ioHub Experiment started {0}".format(getCurrentDateTimeString()))
        self.hub.sendMessageEvent(text="Experiment ID: {0}, Session ID: {1}".format(self.hub.experimentID,
                                                                                    self.hub.experimentSessionID))
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

        ## Handy routines:

        # Frame-skipping check:
        def frame_skip_check(elapsed_t, elapsed_frames):
            # The number of elapsed frames should match the time:
            print('time=%.3f  frames=%d  rate=%.4f' % (elapsed_t, elapsed_frames, (elapsed_t / elapsed_frames)))

        # This is done at every frame update, regardless of trial phase, so predefining:
        def frame_routine():
            flip_time_ = window.flip()
            fix_cross.draw()
            # Checking for quit (the Esc key)
            if event.getKeys(keyList=['escape']):
                exit_routine()
            return flip_time_

        # Also no variation across frames, but only available upon call, which is made only in key registering phase.
        def exit_routine():
            window.close()
            core.quit()

        ## Initiating the trial loop

        n_trials_done = 0

        for trial in trials:

            ## First trial initiates instructions and sends the expt initiation message to the eye tracker:
            if n_trials_done == 0:
                # - Update the instruction screen text...
                instructions_text_stim.setText(instruction_text)
                instructions_text_stim.draw()
                flip_time = window.flip()
                self.hub.sendMessageEvent(text="EXPERIMENT_START", sec_time=flip_time)
                # wait until a space key event occurs after the instructions are displayed
                kb.waitForPresses(keys=' ')

            ## Trial components:

            # Trial components pertaining to time, frames, and trial number:
            n_trials_done += 1
            print('======TRIAL#' + str(n_trials_done) + '======')

            if debug:
                if n_trials_done > 1:
                    # noinspection PyUnboundLocalVariable
                    iti_dur = flip_time - iti_end_trial
                    print('inter-trial duration: %.3f' % iti_dur)

            # Randomize the duration of the post-cue fixation & converting to sec:
            blink_latency = np.random.randint(blink_latency_min,
                                              blink_latency_max + 1) / 1000  # max value has to be one up
            if debug:
                print('blink_latency = %.3f' % blink_latency)

            # Target location:
            # this_targ_loc = np.random.randint(2) * 2 - 1
            this_targ_loc = trial['targ_right'] * 2 - 1
            if this_targ_loc > 0:
                print('target location: Right')
            else:
                print('target location: Left')
            targ.pos = (targ_off_x * this_targ_loc, 0)

            # Cue validity:
            cue_dir = (trial['cue_valid'] * 2 - 1) * this_targ_loc
            # Logic: First, the cue validity is converted from binary [0, 1] to [-1, 1].
            # It is then multiplied by the target location, which is either left [-1] or right [1].
            # E.g., if the cue is valid for a target that appears on the right, cue direction is 1*1, rightward.
            # If the cue is invalid for such a target, cue direction is -1*1, leftward.
            cue_arrow.ori = 90 * (cue_dir - 1)
            # Logic: If cue_dir == 1, 90 * 0 = 0, rightward orientation. If cue_dir == -1, 90 * (-2) = -180, leftward.
            if trial['cue_valid']:
                print('valid cue')
            else:
                print('invalid cue')

            ## Starting the eye-tracking recording.

            # Trial components pertaining to eye blinks:
            blinked = False  # blink detection
            t_blink_start = 0  # blink start & end time
            t_blink_end = 0

            # Recording trial characteristics in the trial output:
            flip_time = window.flip()

            # Starting the recording:
            self.hub.sendMessageEvent(text="TRIAL_START", sec_time=flip_time)
            self.hub.clearEvents('all')
            if not debug:
                tracker.setRecordingState(True)

            ## Starting the frame cycle & waiting for the response.
            trial_t_start = flip_time
            if debug:
                trial_elapsed_frames = 0  # counting frames for frame skip test

            ## Fixation cross:
            fix_1_frames = int(fix_1_dur * frame_rate)
            for fix_1_frame in range(fix_1_frames):
                frame_routine()
                if debug:
                    # noinspection PyUnboundLocalVariable
                    trial_elapsed_frames += 1

            ## The rest of the period without the beep, but with the cue:
            cue_frames = int(cue_dur * frame_rate)
            for cue_frame in range(cue_frames):
                frame_routine()
                cue_arrow.draw()
                if debug:
                    # noinspection PyUnboundLocalVariable
                    trial_elapsed_frames += 1

            ## The brief period with the beep:
            beep_frames = int(beep_dur * frame_rate)
            for beep_frame in range(beep_frames):
                frame_routine()
                beep.play()
                if debug:
                    # noinspection PyUnboundLocalVariable
                    trial_elapsed_frames += 1

            ## Fixation 2 + blink period, i.e., the fixation period after the beep:
            blink_time_period_frames = int((blink_latency + blink_time_window) * frame_rate)
            for blink_time_period_frame in range(blink_time_period_frames):
                flip_time = frame_routine()
                if debug:
                    # noinspection PyUnboundLocalVariable
                    trial_elapsed_frames += 1

            ## Behavioural response: measuring the reaction time:

            # Trial components pertaining to behavioural response:
            targ_resp_given = False
            rt_start = flip_time

            # Displaying the target and measuring the reaction time.
            while not targ_resp_given:

                # Measuring the time it takes for the behavioural response:
                flip_time = frame_routine()

                # Measuring time elapsed since the start of the trial:
                trial_elapsed_t = flip_time - trial_t_start

                # Drawing the target:
                targ.draw()

                if debug:
                    # noinspection PyUnboundLocalVariable
                    trial_elapsed_frames += 1

                ## Monitoring for key presses:
                arrow_keys = event.getKeys(keyList=['left', 'right'])
                if len(arrow_keys) > 0:
                    if 'left' in arrow_keys:
                        print('subject response: Left')
                        beh_resp = -1
                        targ_resp_given = True
                    else:
                        print('subject response: Right')
                        beh_resp = 1
                        targ_resp_given = True
                    if targ_resp_given:  # this is overwritten every time any key is pressed
                        rt = flip_time - rt_start
                        if beh_resp == this_targ_loc:
                            corr_resp = 1  # correct location response
                        else:
                            corr_resp = 0  # incorrect location response
                        print('RT=%.2f correct?=%s' % (rt, corr_resp))
                        if debug:  # in debug mode, check if the frame rate looks okay
                            # noinspection PyUnboundLocalVariable
                            frame_skip_check(trial_elapsed_t, trial_elapsed_frames)
                            iti_end_trial = flip_time

            if not debug:
                tracker.setRecordingState(False)

            ## Recording the data
            output_mat[n_trials_done - 1] = {'exp_name': exp_name,
                                             'subj': exp_info['subj'],
                                             'cond': exp_info['cond'],
                                             'sess': exp_info['sess'],
                                             'sess_id': self.hub.getSessionID(),
                                             'trial_id': n_trials_done,
                                             'targ_right': trial['targ_right'],
                                             'cue_valid': trial['cue_valid'],
                                             'blink_latency': blink_latency,
                                             'trial_start': trial_t_start,
                                             'trial_end': flip_time,
                                             'corr_resp': corr_resp,
                                             'rt': rt}

            # Passing messages to the eye tracker before trial termination:
            self.hub.sendMessageEvent(text="TRIAL_END %d" % n_trials_done, sec_time=flip_time)
            self.hub.clearEvents('all')

        ## Data output:
        data_columns = ['exp_name', 'subj', 'cond', 'sess', 'sess_id', 'trial_id', 'targ_right', 'cue_valid',
                        'blink_latency', 'trial_start', 'trial_end', 'corr_resp', 'rt']
        pd.DataFrame.from_dict(output_mat, orient='index').to_csv(out_file_path + '.csv', index=False,
                                                                  columns=data_columns)
        # .to_csv(out_file_path + '.csv', index=False, columns=data_columns)
        print('output file path is ' + out_file_path)

        ## Termination procedures post-trial
        # So the experiment is done, all trials have been run. Clear the screen and show an 'experiment  done' message
        # using the instructionScreen state. Wait for the trigger to exit that state (i.e. the space key was pressed).

        # Disconnect the eye tracking device:
        tracker.setConnectionState(False)

        # Update the instruction screen text:
        window.flip()
        instructions_text_stim.setText('      Completed!\nPress any key to exit')
        instructions_text_stim.draw()
        flip_time = window.flip()
        self.hub.sendMessageEvent(text="EXPERIMENT COMPLETED", sec_time=flip_time)

        # wait until any key is pressed
        kb.waitForPresses()

        self.hub.sendMessageEvent(text='EXPERIMENT_COMPLETE')


if __name__ == "__main__":
    import os
    from psychopy.iohub import module_directory


    def main(base_dir):
        # -- Sol's comments:
        """
        Creates an instance of the ExperimentRuntime class, gets the eye tracker the user wants to use for the demo, and
        launches the experiment logic.
        """
        # The following code merges a iohub_config file called iohub_config.yaml.part, that has all the iohub_config
        # settings, other than those for the eye tracker, with the eye tracker configs in the yaml files in the
        # eyetracker_configs dir. {...}
        #
        # The merged result is saved as iohub_config.yaml so it can be picked up by the Experiment _runtime as normal.
        # --

        ##
        # I revamped the code so that all iohub config files are stored in the 'iohub_configs' subdirectory. All the
        # .yaml files should be stored there.
        #
        # Also, since a single eye-tracker is used, I removed all references to other eye trackers beside the EyeLink
        # and removed the gui prompt asking for the type of eye tracker.

        eye_config_dir = os.path.normcase(os.path.join(base_dir, 'iohub_configs'))
        base_config_file = os.path.normcase(os.path.join(eye_config_dir, 'iohub_config.yaml.part'))
        eyetracker_config_file = os.path.normcase(os.path.join(eye_config_dir, 'eyelink_config.yaml'))
        combined_config_file_name = os.path.normcase(os.path.join(eye_config_dir, 'iohub_config.yaml'))
        ExperimentRuntime.mergeConfigurationFiles(base_config_file, eyetracker_config_file, combined_config_file_name)
        runtime = ExperimentRuntime(eye_config_dir, 'experiment_config.yaml')
        runtime.start(('SR Research EyeLink',))


    # -- Get the current directory, using a method that does not rely on __FILE__
    # or the accuracy of the value of __FILE__.
    base_dir = module_directory(main)

    # -- Run the main function, which starts the experiment runtime
    main(base_dir)
