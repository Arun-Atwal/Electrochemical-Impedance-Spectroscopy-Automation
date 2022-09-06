"""
This user interface automates variable temperature impedance measurements from a selected range of Linkam heating stages
and impedance analysers. It was adapted as part of a FUSE internship (Faraday Institute) at the University of
Cambridge, Department of Physics by Arun Atwal (asa73@cam.ac.uk) in the summer of 2022.
This work was adapted and expanded on to add functionality with various stages and analysers,
from the user interface code of Adam Alderton, the previous FUSE intern working
on the project. The original author header is below.
"""

# This is user interface to perform variable temperature impedance measurements
# using a heating stage and an impedance analyser. This was created as part of 
# the FUSE studentship (Faraday institute) during the summer of 2020. It is 
# written by Adam Alderton (aa816@exeter.ac.uk) with supervision and assistance
# from Josh Tuffnell (jmt83@cam.ac.uk). Please consult the README before using 
# this software.

import sys, os, time, datetime, json, shutil, ctypes
from ctypes import byref

import numpy as np
import threading
import visa
import tkinter as tk
from tkinter import ttk
from tkinter import font as tkfont

# For better user interface, a third party theme is used. Documentation can be found at https://github.com/rdbende/Sun-Valley-ttk-theme
import sv_ttk

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.animation as animation
from matplotlib import style
from matplotlib.figure import Figure
style.use("ggplot")

from hardware import virtual
from hardware import linkam
from hardware import biologic
from hardware import solartron

os.add_dll_directory("{}".format(os.getcwd()))


class Ramp():
    # This class is to store data relevant to a temperature ramp that may be required of the stage.
    # Sequential ramps should be stored in a list of Ramp() instances.
    
    def __init__(self, num, start_temp, end_temp, rate, interval, min_holdtime, voltage, fmin, fmax, ppd, num_sweeps, sweep_delay):
        self.num = num                          # Integer: Stores the index of the ramp for future reference.
        self.interval = interval                # Float: The temperature by which to change in an interval
        self.start_temp = start_temp            # Float: Starting temperature, (°C)
        self.end_temp = end_temp                # Float: Temperature to reach in num_intervals intervals, (°C)
        self.rate = rate                        # Float: The rate at which to change temperature, (°C/min)
        self.min_holdtime = min_holdtime        # Float: The minimum hold time to hold at each interval BEFORE measuring impedance (seconds). Temperature is held automatically while impedances are measured.
        self.voltage = voltage                  # Float: The A.C voltage amplitude (mV)
        self.fmin = fmin                        # Float: Minimum frequency for frequency sweep (Hz)
        self.fmax = fmax                        # Float: Maximum frequency for frequency sweep (Hz)
        self.ppd = ppd                          # Integer: Points per decade at which to measure impedance
        self.num_sweeps_at_T = num_sweeps       # Integer: Number of impedance measurements to take at each temperature interval
        self.sweep_delay = sweep_delay          # Integer: Time (seconds) between each impedance sweep

        # Calculate "up", describes direction of ramp
        if self.start_temp <= self.end_temp:
            self.up = True
        else:
            self.up = False
        
        # If ramp is downwards, set cooling rate maximum of 15 C/min
        # if self.up == False and self.rate > 15:
            # self.rate = 15
        
        # Calculate num_intervals, Number of intervals in which to change temperature. At each interval, hold by holdtime.
        self.num_intervals = abs(int((self.end_temp - self.start_temp) / self.interval))
        if self.up == False:
            self.interval = abs(self.interval) * -1
        
        # Calculate temperatures (°C)
        self.temps = []
        for i in range(self.num_intervals + 1):
            self.temps.append(round(self.start_temp + (i * self.interval), 2))

        # Calculate numpoints and frequency range
        self.numpoints = int(self.ppd * (np.log10(self.fmax) - np.log10(self.fmin)))
        frange_full = np.logspace(*np.log10( (self.fmin, self.fmax) ), num = self.numpoints)
        self.frange = list(frange_full[::-1])


class vtipy2(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        sv_ttk.set_theme("dark")
        self.check_hardware_config()
        self._frame = None
        self.switch_frame(StartPage)
        global overall_ramp_counter
        overall_ramp_counter = 0  # running count of how many ramp profiles have been created

    def switch_frame(self, frame_class):
        new_frame = frame_class(self)
        if self._frame is not None:
            self._frame.destroy()
        self._frame = new_frame
        self._frame.pack()
    
    def check_hardware_config(self):
        # Search for analysers by attempting to initialise the allowed options
        global available_analysers, stage_name, stage, analyser
        available_analysers = [("Virtual",1)]   # Allow "virtual" analyser by default
        analyser = None
        #analyser = biologic.SP_200()
        print("Detecting Analyser...")
        try:
            print("Searching for Solartron...")
            analyser = solartron.solartron1260()
            available_analysers.append(("Solartron 1260A",3))
            print("Solartron detected")
        except:
            print("Solartron not found.")
        try:
            print("Searching for SP_200...")
            analyser = biologic.SP_200()
            available_analysers.append(("SP_200",2))
            print("SP_200 detected")
        except:
            print("SP_200 not found.")
        print("Available analysers: {}".format(available_analysers))

        try:
            stage = linkam._GeneralStage()
            stage_name= ctypes.create_string_buffer(16)
            stage._process_msg(linkam.Msg.GetStageName, byref(stage_name), 16)
            print(stage_name.value)
        except:
            pass


class StartPage(ttk.Frame):
    def __init__(self, master):
        ttk.Frame.__init__(self, master)
        ttk.Label(self, text = "Arun Atwal - 2022").pack(anchor = tk.NW)
        
        ttk.Label(self, text = "\n EIS Application\n", font = "Helvetica 16 bold").pack()

        ttk.Label(self, text = "Ensure only one stage is connected. If not, close the program and re-run.\n\n"
                              "Select from Available Stages: ", justify = "center",
                 font = tkfont.Font(size = 12)).pack()

        global stage_choice, analyser_choice, returned
        returned = False  # So as to not re-write file headers upon editing ramp data
        stage_choice= tk.IntVar()
        stage_choice.set(1)
        ttk.Radiobutton(self, text = "Virtual", variable= stage_choice, value = 1).pack()
        try:
            ttk.Radiobutton(self, text=stage_name.value, variable=stage_choice, value=2).pack()
        except:
            pass
        ttk.Label(self, text = "\n Select from Available Analysers:", font = tkfont.Font(size=12)).pack()

        analyser_choice = tk.IntVar()
        analyser_choice.set(1)
        for i, val in available_analysers:
            ttk.Radiobutton(self, text = i, variable=analyser_choice, value = val).pack()

        
        ttk.Label(self, text = "\n To begin, please enter a name for the experiment:\n", font = tkfont.Font(size = 12)).pack()
        
        ent_exp_name = ttk.Entry(self, width = 30, font = tkfont.Font(size = 12), justify = "center")
        ent_exp_name.insert(0, "Name")
        ent_exp_name.pack()
        
        start_button = ttk.Button(self, text = "Start", command = lambda: self.start_button_func(master, ent_exp_name, start_button))
        start_button.pack(pady = 15)
        
        ttk.Label(self, text = "Please consult the README before using this software.").pack()
    
    def start_button_func(self, master, entry, start_button):
        
        # Get and store experiment name
        global experiment_name, stage, analyser
        experiment_name = entry.get()
        

        # Set Stage and Analyser choice
        if stage_choice.get() == 1:
            stage = virtual.stage()
        elif stage_choice.get() == 2:
            if stage_name.value == "HFS350EV-PB4    ".encode():
                stage = linkam.HFS350()
            elif stage_name.value == "TS1000-PB4      ".encode():
                stage = linkam.TS1000()

        if analyser_choice.get() == 1:
            analyser = virtual.analyser()
        elif analyser_choice.get() == 2:
            # analyser = biologic.SP_200()
            print(analyser.name)
        elif analyser_choice.get() == 3:
            analyser = solartron.solartron1260()

        print(stage, analyser)

        # Create directory to store experiment files and data
        if not os.path.isdir("experiments"):
            os.mkdir("experiments")
        duplicate_num = 1 # This is to handle duplicate names of experiment
        while True:
            if not os.path.isdir("experiments\\" + experiment_name):
                os.mkdir("experiments\\" + experiment_name)
                break
            else:
                if experiment_name.split("_")[-1] == str(duplicate_num - 1):
                    experiment_name = experiment_name.split("_")[0] + "_{}".format(duplicate_num)
                else:
                    experiment_name += "_{}".format(duplicate_num)
                duplicate_num += 1
                
        
        # Store initial experiment details in a json file for easy reading
        json_handler = lambda obj: (obj.isoformat() if isinstance(obj, (datetime.datetime)) else None)
        with open("experiments\\" + experiment_name + "\\details.json", "w") as jf:
            jsondict = {"experiment_name" : experiment_name,
                        "datetime" : datetime.datetime.now(),
                        "complete" : False,
                        "stage" : stage.name,
                        "analyser" : analyser.name
                        }
            jf.write(json.dumps(jsondict, indent = 4, default = json_handler))
        
        # Switch frame to RampInputPage
        master.switch_frame(RampInputPage)


class RampInputPage(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        try:
            test_name= ctypes.create_string_buffer(16)
            stage._process_msg(linkam.Msg.GetStageName, byref(test_name), 16)
            print(test_name.value)
        except:
            pass
        self.frm_form = tk.Frame()
        self.frm_form.pack()
        
        # Store ramps in list of Ramp instances (Global for now, consider changing)
        global ramps, returned, overall_ramp_counter

        if returned:
            overall_ramp_counter += len(ramps)

        # print(overall_ramp_counter)

        self.ramp_count = overall_ramp_counter
        ramps = []

        
        self.ramp_attributes = [
            "Starting Temperature (°C):",
            "End Temperature (°C):",
            "Rate (°C/min):",
            "Temperature Interval (°C):",
            "Minimum Holdtime (s):",
            "Voltage Amplitude (mV):",
            "Minimum Frequency (Hz):",
            "Maximum Frequency (Hz):",
            "Points per Decade:",
            "Sweeps at T:",
            "Sweep Delay (s):"
        ]
        
        default_entries = [
            "30.0",     # start_temp
            "60.0",     # end_temp
            "5.0",      # rate
            "10",       # interval
            "300",       # min_holdtime
            "50",       # voltage amplitude
            "1e-1",     # fmin
            "3e6",      # fmax
            "10",       # ppd
            "1",        # num_sweeps_at_T
            "300"        # sweep_delay
        ]
        
        # Organise labels and entries
        self.entries = []
        for i, attr in enumerate(self.ramp_attributes):
            label = ttk.Label(master = self.frm_form, text = attr, font = tkfont.Font(size = 12), anchor = tk.W, justify = tk.LEFT)
            entry = ttk.Entry(master = self.frm_form, width = 10, font = tkfont.Font(size = 12), justify = tk.RIGHT)
            self.entries.append(entry)
            entry.insert(0, default_entries[i])
            label.grid(row = i, column = 0, sticky = tk.W, pady = 3)
            entry.grid(row = i, column = 2, sticky = tk.E, pady = 3)
        
        # Organise ranges
        tk.Label(master = self.frm_form, text = "{:.1f} \u2192 {:.1f}".format(stage.min_temp, stage.max_temp), font = tkfont.Font(size = 12), width = 12).grid(row = 0, column = 1, pady = 3) # starting temperature range
        tk.Label(master = self.frm_form, text = "{:.1f} \u2192 {:.1f}".format(stage.min_temp, stage.max_temp), font = tkfont.Font(size = 12), width = 12).grid(row = 1, column = 1, pady = 3) # end temperature range
        tk.Label(master = self.frm_form, text = "{:.1f} \u2192 {:.1f}".format(stage.min_rate, stage.max_rate), font = tkfont.Font(size = 12), width = 12).grid(row = 2, column = 1, pady = 3) # heater rate (C/min)
        tk.Label(master = self.frm_form, text = "\u2265 1", font = tkfont.Font(size = 12), width = 12).grid(row = 3, column = 1, pady = 3) # temperature interval
        tk.Label(master = self.frm_form, text = "10 \u2192 3e7", font = tkfont.Font(size = 12), width = 12).grid(row = 4, column = 1, pady = 3) # minimum holdtime
        tk.Label(master = self.frm_form, text = "{} \u2192 {}".format(analyser.min_V, analyser.max_V), font = tkfont.Font(size = 12), width = 12).grid(row = 5, column = 1, pady = 3) # minimum holdtime
        tk.Label(master = self.frm_form, text = "\u2265 {:.2e}".format(analyser.min_frequency), font = tkfont.Font(size = 12), width = 12).grid(row = 6, column = 1, pady = 3) # minimum frequency
        tk.Label(master = self.frm_form, text = "\u2264 {:.2e}".format(analyser.max_frequency), font = tkfont.Font(size = 12), width = 12).grid(row = 7, column = 1, pady = 3) # maximum frequency
        tk.Label(master = self.frm_form, text = "\u2265 1", font = tkfont.Font(size = 12), width = 12).grid(row = 8, column = 1, pady = 3) # ppd
        tk.Label(master = self.frm_form, text = "\u2265 0", font = tkfont.Font(size = 12), width = 12).grid(row = 9, column = 1, pady = 3) # num sweeps at T
        tk.Label(master = self.frm_form, text = "\u2265 5", font = tkfont.Font(size = 12), width = 12).grid(row = 10, column = 1, pady = 3) # sweep delay

        self.frm_buttons = tk.Frame()
        self.frm_buttons.pack(fill = tk.X, padx = 5, pady = 5)

        # Finish button
        btn_finish = ttk.Button(master = self.frm_buttons, text = "Finish", command = lambda: self.save_ramps(master))
        btn_finish.pack(side = tk.RIGHT, padx = 10)
        
        # Add ramp button
        btn_ramp = ttk.Button(master = self.frm_buttons, text = "Create Ramp", command = self.input_ramp)
        btn_ramp.pack(side=tk.RIGHT, padx = 10)
        
        # Clear all button
        btn_clear = ttk.Button(master = self.frm_buttons, text = "Clear All", command = self.clear_all)
        btn_clear.pack(side = tk.LEFT, padx = 10)
        
        # Clear previous button
        btn_clear_prev = ttk.Button(master = self.frm_buttons, text = "Clear Prev", command = self.clear_prev)
        btn_clear_prev.pack(side = tk.LEFT, padx = 10)
        
        # Ramp counter
        self.lbl_counter = ttk.Label(master = self.frm_buttons, text = "Total Ramps = " + str(self.ramp_count), font = tkfont.Font(size = 12))
        self.lbl_counter.pack(side = tk.RIGHT, padx = 10)

    def increase_ramp_counter(self):
        self.ramp_count += 1
        self.lbl_counter["text"] = "Total Ramps = " + str(self.ramp_count)

    def clear_all(self):
        # for entry in self.entries:
            # entry.delete(0, tk.END)
        
        ramps.clear()
        
        self.ramp_count = 0
        self.lbl_counter["text"] = "Total Ramps = " + str(self.ramp_count)

    def clear_prev(self):
        ramps.pop()
        self.ramp_count -= 1
        self.lbl_counter["text"] = "Total Ramps = " + str(self.ramp_count)

    def sanity_check(self, ramp):
        global stage, analyser
        
        if not stage.min_temp <= ramp.start_temp <= stage.max_temp:
            return False, "Invalid Start Temp"
        if not stage.min_temp <= ramp.end_temp <= stage.max_temp:
            return False, "Invalid End Temp"
        if not stage.min_rate <= ramp.rate <= stage.max_rate:
            return False, "Invalid Rate"
        if stage.name == "TS1000" and ramp.rate > 50:
            tk.messagebox.showwarning(title="Rate Warning", message="""It is not recommended to regularly exceed a rate
            of 50°C/min on the TS1000 Stage. Proceed with care.""")
        if abs((ramp.end_temp - ramp.start_temp) / ramp.interval) % 1 != 0: # If interval does not divide in \Delta{T}
            return False, "Temperature interval does not divide range"
        if abs(ramp.interval) < 1:
            return False, "Interval"
        if not 10 <= ramp.min_holdtime <= 3e7:
            return False, "Minimum Holdtime"
        if not analyser.min_V <= ramp.voltage <= analyser.max_V:
            return False, "Voltage"
        if ramp.fmin >= ramp.fmax:
            return False, "fmin >= fmax"
        if ramp.fmin < analyser.min_frequency:
            return False, "Minimum Frequency"
        if ramp.fmax > analyser.max_frequency:
            return False, "Maximum Frequency"
        if ramp.ppd < 1:
            return False, "PPD must be greater than 1"
        if not isinstance(ramp.ppd, int):
            return False, "PPD not int"
        if ramp.num_sweeps_at_T < 0:
            return False, "No. of Sweeps"
        if not isinstance(ramp.num_sweeps_at_T, int):
            return False, "No. of Sweeps"
        if ramp.sweep_delay < 5:
            return False, "Sweep Delay must be greater than 5"
        
        return True, None

    def input_ramp(self):
        entries_types_valid = False
        try:
            ramp = Ramp(
                num = self.ramp_count + 1, # Such that first ramp is ramp 1 etc
                start_temp =    float(self.entries[0].get()),
                end_temp =      float(self.entries[1].get()),
                rate =          float(self.entries[2].get()),
                interval =      int(self.entries[3].get()),
                min_holdtime =  int(self.entries[4].get()),
                voltage =       float(self.entries[5].get()),
                fmin =          float(self.entries[6].get()),
                fmax =          float(self.entries[7].get()),
                ppd =           int(self.entries[8].get()),
                num_sweeps =    int(self.entries[9].get()),
                sweep_delay =   int(self.entries[10].get())
            )
            entries_types_valid = True
        except:
            tk.messagebox.showwarning(title = "Input Error", message = "Input(s) invalid.\nPlease check inputs!")
        if entries_types_valid:
            entries_types_valid, fault = self.sanity_check(ramp)
            
            if entries_types_valid: # if input values are valid, proceed

                # Swap start_temp and end_temp entries such that user easily inputs reverse of ramp just entered
                if len(ramps) == 0:
                    n_ramp = ramp
                    new_start_temp = n_ramp.end_temp
                    new_end_temp = n_ramp.start_temp
                else:
                    n_ramp = ramps[-1]
                    new_start_temp = n_ramp.start_temp
                    new_end_temp = n_ramp.end_temp

                for entry in self.entries:
                    entry.delete(0, tk.END)
                
                self.entries[0].insert(0, str(new_start_temp))
                self.entries[1].insert(0, str(new_end_temp))
                self.entries[2].insert(0, str(n_ramp.rate))
                self.entries[3].insert(0, str(abs(n_ramp.interval)))
                self.entries[4].insert(0, str(n_ramp.min_holdtime))
                self.entries[5].insert(0, str(n_ramp.voltage))
                self.entries[6].insert(0, "{:.2e}".format(n_ramp.fmin))
                self.entries[7].insert(0, "{:.2e}".format(n_ramp.fmax))
                self.entries[8].insert(0, str(n_ramp.ppd))
                self.entries[9].insert(0, str(n_ramp.num_sweeps_at_T))
                self.entries[10].insert(0, str(n_ramp.sweep_delay))

                # Store ramp
                ramps.append(ramp)
                self.increase_ramp_counter()
                
            else: # if inputs are not valid
                tk.messagebox.showwarning(title = "Input Error", message = "Input(s) invalid.\nFAULT: {}\n".format(fault))

    def save_ramps(self, master):
        global returned, overall_ramp_counter
        if self.ramp_count != overall_ramp_counter:
            # Open existing details, and allocate space to store ramp data
            with open("experiments\\{}\\details.json".format(experiment_name), "r") as jf:
                jsondict = json.load(jf)

            if not returned:    # initialise
                jsondict["ramps"] = []

            water_warning_shown = False
            with open("experiments\\{}\\details.json".format(experiment_name), "w") as jf:
                for ramp in ramps:
                    if not water_warning_shown:
                        if ramp.end_temp >= 200 or ramp.start_temp >= 200:
                            water_warning_shown = True
                            tk.messagebox.showwarning(title = "Cooling Water Warning", message = "Temperature will go above 200°C.\nCheck cooling water before proceeding!")
                    ramp_dict = ramp.__dict__.copy()
                    del(ramp_dict["frange"]) # don't store frequency range, as fmin, fmax and ppd are already stored
                    jsondict["ramps"].append(ramp_dict)
                jf.write(json.dumps(jsondict, indent = 4))

            # Destroy frames and open details file
            self.frm_buttons.destroy()
            self.frm_form.destroy()
            master.switch_frame(ReviewRampsPage)
        else:
            tk.messagebox.showwarning(title="Input Error", message="Must input at least one ramp.")


class ReviewRampsPage(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        
        self.process_ramps(master)
        
        tk.Label(self, text = "\nPlease see a schematic of the experiment above.\n\n",
                 font = tkfont.Font(size = 12)).pack()

        tk.Label(self, text = "Estimated Time = {} (H:M:S)\n".format(estimated_time),
                 font = tkfont.Font(size = 18)).pack()
        
        frm_buttons = tk.Frame(self)
        frm_buttons.pack()
        
        ttk.Button(frm_buttons, text = "Begin Experiment",command = lambda: self.begin_experiment_button_func(master)).pack(side = tk.RIGHT, padx = 5, pady = 15)
    
        ttk.Button(frm_buttons, text = "Start Over", command = lambda: self.start_over(master)).pack(side = tk.RIGHT, padx = 5, pady = 15)
    
        self.two_hr_warning_shown = False
        
    def start_over(self, master):
        self.plot_widget.destroy()
        master.switch_frame(RampInputPage)
        
    def begin_experiment_button_func(self, master):
        global estimated_time
        if self.two_hr_warning_shown == True or estimated_time.total_seconds() <= 120 * 60: # two hours in seconds
            self.plot_widget.destroy()
            estimated_time = str(estimated_time)
            master.switch_frame(RunningPage)
        else:
            self.two_hr_warning_shown = True
            tk.messagebox.showwarning(title = "2 Hour Warning", message = "Estimated time is above 2 hours.\nEnsure stage is not held at a high temperature for more than 2 hours at a time over the temperature profile!")
    
    def process_ramps(self, master):
        # Plots an estimated schematic of the temperature profile over time, as well as calculates the estimated time for the whole experiment
        # Note that estimated time is also calculated in the ramp input page, but it is messy to 
        # The time taken to take an impedance scan is roughly t(f) = 1/f + 2 for the solartron, as was measured manually and fit to by a curve.
        # For the biologic, the time was fit to a polynomial
        # These times are then summed over for the entire frequency sweep
            
        def dt_sweep(ramp):
            t = 0
            frange = ramp.frange
            if analyser.name == "Virtual Analyser" or analyser.name == "Solartron 1260" or analyser.name == "Biologic SP-200":
                t += round(1.8047 * ramp.numpoints)
                for f in frange:
                    t += round(1.0025 / f)
            # elif analyser.name == "Biologic SP-200":
            #     t += round(1.76e1 * ramp.numpoints)
            #     for f in frange:
            #         t += (-3.12e-8 * f**3) + (5.37e-5 * f**2) + (-3e-2 * f) + (2.36e1 * f**-1) + (-3.11e1 * f**-2) + (1.84e1 * f**-3)
            else:
                raise Exception("Analyser type unrecognised when calculating estimated time.")
            return abs((ramp.min_holdtime + (ramp.num_sweeps_at_T - 1) * ramp.sweep_delay + ramp.num_sweeps_at_T * t) / 60) # Convert to minutes
        
        global estimated_time
        estimated_time = 0 
        
        fig = plt.figure()
        ax = plt.axes()
        
        # Take into account first heating
        start_temp = stage.get_temperature()
        dT = abs(ramps[0].temps[0] - start_temp)
        dt = abs(dT / ramps[0].rate)
        ax.plot([0, dt], [start_temp, ramps[0].temps[0]], color = "red" if ramps[0].temps[0] >= start_temp else "blue", linewidth = 2)
        estimated_time += dt
        
        # Plot first sweep and hold time for first ramp
        first_hold_time = dt_sweep(ramps[0])
        ax.plot([estimated_time, estimated_time + first_hold_time], [ramps[0].temps[0], ramps[0].temps[0]], ":", color = "black", linewidth = 2)
        estimated_time += first_hold_time
        
        for ramp in ramps:
            temps = ramp.temps
            rate = ramp.rate
            
            num_intervals = ramp.num_intervals
            
            # Time taken at each interval to hold and sweep
            ramp_dt_sweep = dt_sweep(ramp)
            
            for i in range(num_intervals):
                dT = abs(temps[i + 1] - temps[i])
                dt = abs(dT / rate) # in minutes
                
                # Plot temperature change line
                ax.plot([estimated_time, estimated_time + dt], [temps[i], temps[i + 1]], color = "red" if temps[i + 1] >= temps[i] else "blue", linewidth = 2)
                estimated_time += dt
                # Plot hold and sweep line
                ax.plot([estimated_time, estimated_time + ramp_dt_sweep], [temps[i + 1], temps[i + 1]], ":", color = "black", linewidth = 2)
                
                estimated_time += ramp_dt_sweep 
        
        ax.set_xlabel("Time [mins]")
        ax.set_ylabel("Temperature [°C]")
        
        canvas = FigureCanvasTkAgg(fig, master = master)
        self.plot_widget = canvas.get_tk_widget()
        self.plot_widget.pack()
        fig.canvas.draw()
        
        # Convert estimated time to hours:minutes:seconds, rounding up to the nearest minute
        estimated_time = datetime.timedelta(minutes = estimated_time)
        estimated_time += datetime.timedelta(seconds = np.floor(60 - (estimated_time.seconds % 60)))
        estimated_time -= datetime.timedelta(microseconds = estimated_time.microseconds)


class RunningPage(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        
        # Info frame
        self.frm_info_master = tk.Frame(self)
        self.frm_info_master.grid(row = 0, column = 0)
        self.setup_info(self.frm_info_master)
        self.setup_status_box(self.frm_info_master)
        
        # Live plot frame
        self.frm_plot_master = tk.Frame(self)
        self.frm_plot_master.grid(row = 0, column = 1, rowspan = 2)
        self.setup_plot(self.frm_plot_master)
        
        # Stop / Finish Button
        self.frm_stop = tk.Frame(self, pady = 20)
        self.frm_stop.grid(row = 1, column = 0)
        self.paused = False
        self.btn_pause = tk.Button(self.frm_stop, text = "PAUSE", font = "Helvetica 22 bold", fg = "green", bg = "white",
                                   width = 12, height = 1, command = lambda: self.pause_button_func(master))
        self.btn_pause.pack()
        self.btn_return = tk.Button(self.frm_stop, text="EDIT RAMPS", font="Helvetica 22 bold", fg="green", bg = "white",
                                    width=12, height=1, command=lambda: self.return_func(master))
        self.btn_return.pack()
        self.btn_stop = tk.Button(self.frm_stop, text = "STOP", font = "Helvetica 22 bold", bg = "red", fg = "white",
                                  width = 12, height = 1, command = lambda: self.stop_button_func(master))
        self.btn_stop.pack()
        
        # Begin main experiment thread
        self.experiment_thread = threading.Thread(target = self.experiment_mainloop, args = (master,))
        self.experiment_thread.daemon = True
        self.experiment_thread.start()
        
        # Save the current time to measure the time elapsed
        global time_elapsed
        if not returned:
            time_elapsed = 0
        self.start_time = time.perf_counter() - time_elapsed

        if not returned:
            with open("experiments\\" + experiment_name + "\\details.json", "r") as jf:
                self.start_time_datetime = datetime.datetime.strptime(json.load(jf)["datetime"], '%Y-%m-%dT%H:%M:%S.%f')

    def return_func(self, master):
        global returned, ramping, time_elapsed, holding, currently_returned
        time_elapsed = time.perf_counter()
        returned = True
        currently_returned = True
        self.stop_thread = True
        if ramping or holding:     # on ramp or hold i.e. only if a sweep has not been started
            if not self.paused:
                if ramping:
                    stage.toggle_hold()
                self.paused = True
            self.frm_info_master.destroy()
            self.frm_plot_master.destroy()
            self.frm_stop.destroy()
            master.switch_frame(RampInputPage)


    def pause_button_func(self, master):
        global ramping, holding
        if ramping or holding:     # on ramp or hold i.e. only if a sweep has not been started
            if not self.paused:
                # self.stop_thread = True
                if ramping:
                    stage.toggle_hold()
                self.update_status("Paused")
                self.btn_pause.configure(text = "RESUME")
                self.paused = True
            else:
                self.btn_pause.configure(text="PAUSE")
                # self.stop_thread = False
                self.paused = False
                if ramping:
                    stage.toggle_hold()

        
    def stop_button_func(self, master):
        # stage.stop_heating()
        
        with open("experiments\\" + experiment_name + "\\details.json", "r") as jf:
            json_dict = json.load(jf)
            json_dict["time_elapsed"] = str(self.info_data["Time Elapsed"].cget("text"))
        
        with open("experiments\\" + experiment_name + "\\details.json", "w") as jf:
            jf.write(json.dumps(json_dict, indent = 4))
        
        try:
            stage.move_log(experiment_name)
        except:
            pass
        
        self.stop_thread = True
        time.sleep(0.5)
            
        self.save_plot()
        
        app.destroy()
        """POTENTIAL WORRY:
        Stopping the stage heating at the start of this function did not work (without warning)
        when in the middle of a ramp.
        Doing it here after the other commands seems to prevent this issue.
        This is potentially dangerous and care should be taken."""
        stage.stop_heating()
        sys.exit()
    
    def finish_button_func(self, master):
        app.destroy()
        sys.exit()

    def setup_status_box(self, master):
        # Simple way to create a simple gap info frame and status frame
        lbl_running = tk.Label(master)
        lbl_running.pack()
        
        frm_status_box = tk.Frame(master)
        frm_status_box.pack()
        
        lbl_status_title = ttk.Label(master = frm_status_box, text = "CURRENT STATUS", font = "Helvetica 13 bold")
        lbl_status_title.pack()
        
        self.lbl_status = ttk.Label(master = frm_status_box, font = tkfont.Font(size = 12), text = "Experiment is Starting.", width = 45, relief = tk.SUNKEN, borderwidth = 5)
        self.lbl_status.pack(pady = 20)

    def setup_info(self, master):
        frm_info = tk.Frame(master)
        frm_info.pack()

        # Store data for easy updating
        self.info_labels = {}
        self.info_data = {}

        def setup_row(frame, label, row):
            self.info_labels[label] = tk.Label(master = frame, text = label + " =", font = tkfont.Font(size = 12), anchor = tk.W, justify = tk.LEFT, width = 25)
            self.info_data[label] = tk.Label(master = frame, text = "-", font = tkfont.Font(size = 12), anchor = tk.E, justify = tk.RIGHT, width = 20)
            
            self.info_labels[label].grid(row = row, column = 0, sticky = tk.W, pady = 5)
            self.info_data[label].grid(row = row, column = 1, sticky = tk.E, pady = 5)

        # Frame 1: Experiment Name, Estimated Total Time and Time Elapsed
        frm_1 = ttk.Frame(master = frm_info, relief = tk.SUNKEN, borderwidth = 5)
        frm_1.pack(pady = 1)
        setup_row(frm_1, "Experiment Name", 0)
        setup_row(frm_1, "Estimated Total Time", 1)
        setup_row(frm_1, "Time Elapsed", 2)
        
        # Space between frame 1 and frame 2
        frm_space_1 = tk.Label(master = frm_info, pady = 1)
        frm_space_1.pack()

        # Frame 2: Temperature
        frm_2 = ttk.Frame(master = frm_info, relief = tk.SUNKEN, borderwidth = 5)
        frm_2.pack(pady = 1,)
        setup_row(frm_2, "Temperature (°C)", 0)
        
        # Space between frame 2 and frame 3
        frm_space_2 = tk.Label(master = frm_info, pady = 1)
        frm_space_2.pack()
        
        # Frame 3: Ramp, Temperature Interval, Impedance Sweep Number
        frm_3 = ttk.Frame(master = frm_info, relief = tk.SUNKEN, borderwidth = 5)
        frm_3.pack(pady = 1,)
        setup_row(frm_3, "Ramp", 0)
        
        # Set up custom row displaying ramp details
        self.info_labels["Ramp Details"] = tk.Label(master = frm_3, font = tkfont.Font(size = 12), width = 35, relief = tk.GROOVE, borderwidth = 2)
        self.info_labels["Ramp Details"].grid(row = 1, column = 0, columnspan = 2, pady = 5, padx = 5)
        
        # Continue Frame 3
        setup_row(frm_3, "Temperature Point", 2)
        setup_row(frm_3, "Impedance Sweep", 3)
        setup_row(frm_3, "Impedance Scan", 4)
        
        # Insert experiment name and estimated time
        self.info_data["Experiment Name"]["text"] = experiment_name
        self.info_data["Estimated Total Time"]["text"] = estimated_time

    def setup_plot(self, master):
        fig_ani = Figure(figsize = (10, 6))
        ax_ani = fig_ani.add_subplot(111)
        
        with open("experiments\\" + experiment_name + "\\details.json") as jf:
            ramps = json.load(jf)["ramps"]
        
        terminating_temps = []
        for ramp in ramps:
            terminating_temps.append(ramp["start_temp"])
            terminating_temps.append(ramp["end_temp"])
        
        min_temp = min(terminating_temps) - 1
        max_temp = max(terminating_temps) + 1
        
        def animate(i):
            if os.path.exists("experiments\\" + experiment_name + "\\temperature_data.txt"):
                lines = open("experiments\\" + experiment_name + "\\temperature_data.txt").readlines()[1:]
                data = np.asarray([line.strip().split(',') for line in lines[0:]]).astype(np.float)

                temp_data = data[:,1]
                time_data = [seconds / 60 for seconds in data[:,0]] # change to minutes
                
                # If more than 30 minutes have pased, only plot last 30 minutes
                time_data = time_data[-1800:]
                temp_data = temp_data[-1800:]

                ax_ani.clear()
                ax_ani.plot(time_data, temp_data)
                
                ax_ani.set_xlabel('Time [mins]')
                ax_ani.set_ylabel('Temperature [$^o$C]')
                
                ax_ani.set_xlim(right = max(5, time_data[-1] + 1)) # axis range will never be below 5 minutes
                
                y_low = min_temp
                y_high = max_temp
                
                highest_temp = max(temp_data)
                lowest_temp = min(temp_data)
                
                if highest_temp > y_high:
                    y_high = highest_temp
                if lowest_temp < y_low:
                    y_low = lowest_temp

                y_low *= 0.9
                y_high *= 1.1   # add margins of 10% for clarity on boundaries.

                ax_ani.set_ylim(y_low, y_high)
                
            else:
                time.sleep(4)

        global returned
        if not returned:
            with open("experiments\\" + experiment_name + "\\temperature_data.txt", "a") as f:
                f.write("seconds, temperature\n")
 
        canvas = FigureCanvasTkAgg(fig_ani, master = master)
        canvas.draw()
        canvas.get_tk_widget().pack()
        canvas._tkcanvas.pack()
        
        app.ani = animation.FuncAnimation(fig_ani, animate, frames = None, repeat = False)

    def save_plot(self):
        fig, ax = plt.subplots(1,1,figsize=(10,6))

        lines = open("experiments\\" + experiment_name + "\\temperature_data.txt").readlines()[1:]
        data = np.asarray([line.strip().split(',') for line in lines[0:]]).astype(np.float)

        temp_data = data[:,1]
        time_data = [seconds / 60 for seconds in data[:,0]]

        ax.plot(time_data, temp_data)
        ax.set_xlabel('Time [mins]')
        ax.set_ylabel('Temperature [$^o$C]')

        fig.canvas.set_window_title(experiment_name)

        plt.savefig("experiments\\" + experiment_name + "\\temp_profile.png", dpi = 200)

    def update_measurement_info(self, temp, ramp, stage, temp_interval_num, sweep_num):
        global overall_ramp_counter
        self.info_data["Ramp"]["text"] = "{} / {}".format(ramp.num, overall_ramp_counter + len(ramps))
        self.info_labels["Ramp Details"]["text"] = "[ {:.2f} \u2192 {:.2f} ] °C, \u0394T = {:.2f} K".format(ramp.start_temp, ramp.end_temp, ramp.interval)
        self.info_data["Temperature Point"]["text"] = "{} °C    {} / {}".format(round(temp, 2), temp_interval_num, ramp.num_intervals + 1)
        self.info_data["Impedance Sweep"]["text"] = "{} / {}".format(sweep_num, ramp.num_sweeps_at_T)
        self.info_data["Impedance Scan"]["text"] = "{} Hz    {} / {}".format("-", "-", len(ramp.frange))

    def _update_temperature_and_time(self, stage):
        while not self.stop_thread:
            temperature = round(stage.get_temperature(), 2)
            seconds = np.floor(time.perf_counter() - self.start_time)
            self.info_data["Time Elapsed"]["text"] = str(datetime.timedelta(seconds = seconds))
            self.info_data["Temperature (°C)"]["text"] = str(temperature)

            with open("experiments\\" + experiment_name + "\\temperature_data.txt", "a") as f:
                f.write("{:.2f},{:.2f}\n".format(seconds, temperature))
                
            time.sleep(1)

    def update_status(self, status):
        self.lbl_status["text"] = status

    def experiment_mainloop(self, master):
        global currently_returned
        currently_returned = False
        # Spawn subprocesses to update temperature and time every second
        self.stop_thread = False    # this corresponds to updating measurement info
        self.complete = False
        self.temp_update_process = threading.Thread(target = self._update_temperature_and_time, args = (stage,))
        self.temp_update_process.daemon = True
        self.temp_update_process.start()

        complete_count = len(ramps) # to check if all ramps were completed

        # for thread in threading.enumerate():
            # print(thread.name)

        for ramp in ramps:
            try:
                # Create directory to store ramp data
                os.mkdir("experiments//" + experiment_name + "//ramp_{}_{}-{}".format(ramp.num, ramp.start_temp, ramp.end_temp))
                # Measure impedances and varying temperatures along the ramp
                self.ramp_and_measure(ramp, stage, analyser)
            except:
                pass

        time.sleep(0.3)   # wait for ramps to clear in event of editing ramps

        # print(ramps)
        # print(complete_count)
        if len(ramps) == complete_count:
            self.complete = True

        # print(self.complete)

        try:    # if paused to edit ramps, this may fail
            if self.complete:
                self.update_status("Experiment Complete")

                # Change experiment to being complete in details.json
                with open("experiments\\{}\\details.json".format(experiment_name), "r") as jf:
                    jsondict = json.load(jf)
                    jsondict["complete"] = True
                    jsondict["time_elapsed"] = self.info_data["Time Elapsed"].cget("text")

                with open("experiments\\{}\\details.json".format(experiment_name), "w") as jf:
                    jf.write(json.dumps(jsondict, indent = 4))

                # Stop hardware and threads, and change STOP button to FINISH
                stage.stop_heating()
                try:
                    stage.move_log(experiment_name)
                except:
                    pass

                self.stop_thread = True
                time.sleep(0.5)

                self.save_plot()

                self.btn_stop.configure(text = "FINISH", font = "Helvetica 22 bold", bg = "green",
                                        width = 12, height = 1, command = lambda: self.finish_button_func(master))
        except:
            pass

    def ramp_and_measure(self, ramp, stage, analyser):
        # print(self.paused)
        global ramping, holding, currently_returned
        ramping = False
        holding = False
        for temp_interval_num, temp in enumerate(ramp.temps):
            self.update_measurement_info(temp, ramp, stage, temp_interval_num + 1, 0)
            self.btn_pause.configure(fg="green")
            self.btn_return.configure(fg="green")
            ramping = True

            # Heat up to desired temperature
            #first_attempt = True
            while True:
                #if not first_attempt:
                #    stage.stop_heating()
                #first_attempt = False
                if not self.paused:
                    stage._start_heating(temp, ramp.rate, ramp.min_holdtime)
                    time.sleep(3) # if important as if within a few degrees, stage may be returning that it's holding but we actually need to keep heating

                    self.update_status("Heating/Cooling to {} °C".format(temp))

                    # Holdtime is 0.0 while stage is heating/cooling
                    i = 0
                    while stage.get_holdtime_remaining() == 0.0 and not self.paused:
                        # INEFFICIENT but implemented due to safety concerns:
                        # Pausing at the wrong time could cause indefinite heating
                        # compares heating direction every 10 seconds
                        i += 1
                        if i == 1:
                            t1 = stage.get_temperature()
                        elif i == 10:
                            t2 = stage.get_temperature()
                            i = 0
                        if i == 0:
                            sign = t2 - t1
                            if sign < 0 and ramp.up and not ramp.start_temp < t2: # and not in process of decreasing to start temp
                                stage._start_heating(temp, ramp.rate, ramp.min_holdtime)
                            elif sign > 0 and not ramp.up and not ramp.start_temp > t2: # and not in process of increasing to start temp
                                stage._start_heating(temp, ramp.rate, ramp.min_holdtime)
                            else:
                                # print(temp, ramp.rate)
                                pass
                        if currently_returned:
                            break
                        time.sleep(1)

                    # If temperature has not been reached, try again.
                    if temp - 1 <= stage.get_temperature() <= temp + 1:
                        break
                elif currently_returned:   # finish and close the thread without measuring further
                    break
                else:
                    time.sleep(1)
            ramping = False

            holding = True
            # Hold at starting temperature
            while True:
                if currently_returned:
                    break
                holdtime_remaining = stage.get_holdtime_remaining()
                if not self.paused:
                    self.update_status("Holding at T = {:.2f} °C, Holdtime Remaining = {:.2f}(s)".format(temp, round(holdtime_remaining)))
                time.sleep(1)
                
                if holdtime_remaining <= 0.0 and not self.paused:
                    break
            holding = False

            self.btn_pause.configure(fg="red")
            self.btn_return.configure(fg="red")

            # stage.stop_heating()
            # holdtime_remaining = ramp.min_holdtime
            # while holdtime_remaining > 0.0:
            #     self.update_status("Holding at T = {:.2f} °C, Holdtime Remaining = {:.2f}(s)".format(temp, round(holdtime_remaining)))
            #     time.sleep(1)
            #     holdtime_remaining -= 1
            
            if ramp.num_sweeps_at_T != 0 and not currently_returned:
                # Toggle indefinite temperature hold ON to begin impedance measurements
                stage.toggle_hold()
                
                # If temperature has not stabilised to within 0.01 C, keep holding
                if not isinstance(stage, virtual.stage):
                    stage.update_tolerance(temp)
                    while not (temp - stage.tolerance <= stage.get_temperature() <= temp + stage.tolerance):
                        self.update_status("Temperature not T = {:.2f} °C, holding...".format(temp))
                        time.sleep(1)
                
                # Make desired measurements
                for sweep_num in range(ramp.num_sweeps_at_T):
                    self.update_measurement_info(temp, ramp, stage, temp_interval_num + 1, sweep_num + 1)
                    
                    # If not the first sweep, wait sweep_delay seconds
                    if sweep_num != 0:
                        for seconds_elapsed in range(ramp.sweep_delay):
                            self.update_status("Waiting to perform next impedance sweep, t = {}s".format(ramp.sweep_delay - seconds_elapsed))
                            time.sleep(1)
                    
                    # Get temperature as measured by the stage to include in the metadata of the sweep
                    Tcell = stage.get_temperature()
                    
                    # Begin impedance measurement subprocess
                    self.update_status("Beginning sweep {} / {}.".format(sweep_num + 1, ramp.num_sweeps_at_T))
                    Z_measurement_process = analyser.measure_impedance(self.info_data["Impedance Scan"], ramp, sweep_num, Tcell, experiment_name, temp_interval_num)
                    time.sleep(1) # Wait for process to start
                    
                    # Obtain temperature every second until impedance sweeps are complete
                    self.update_status("Performing sweep {} / {}.".format(sweep_num + 1, ramp.num_sweeps_at_T))
                    while Z_measurement_process.is_alive():
                        time.sleep(1)
                
                # Now that impedance measurements are complete, stop holding at temperature and proceed
                stage.toggle_hold()
                
            else:
                self.update_status("No sweeps to be carried during this ramp.")
        
            self.update_measurement_info(temp, ramp, stage, temp_interval_num + 1, sweep_num + 1)

def on_quit():
    try:
        stage.stop_heating()
    except:
        pass
    app.destroy()
    sys.exit()

if __name__ == "__main__":
    global app
    app = vtipy2()
    app.title("EIS asa73")
    app.protocol("WM_DELETE_WINDOW", on_quit)
    app.mainloop()
    sys.exit()