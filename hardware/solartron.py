"""
This is an adaptation of work carried out by Jack Hodkinson as part of his PhD thesis:
https://github.com/jackhodkinson/vtipy
"""

import os, time
import visa
import numpy as np, sys, datetime
import threading


class solartron1260():
    def __init__(self,
                 Vdc = 0.0,             # Float: Optional, D.C bias.
                 integration_time = 1,  # Integer: Optional, Integration time (s).
                 gpib_address = 8):     # Integer: Optional, GPIB address as determined by switches on the Solartron1260.

        self.name = "Solartron 1260"
        self.min_frequency = 1e-5 # from manual
        self.max_frequency = 3.2e7
        self.min_V = 1e-6
        self.max_V = 1e3
        
        # Open solartron connection using pyvisa (using NI-VISA backend)
        ADDRESS = u"GPIB0::{}::INSTR".format(gpib_address)
        rm = visa.ResourceManager()
        self.analyser = rm.open_resource(ADDRESS)
        
        # Set read and write termination characters, as determined by switches on solartron (off, off)
        self.analyser.read_termination = "\r\n"
        self.analyser.write_termination = "\r\n"
        
        self.send("*RST")   # Resets analyser
        self.send("*CLS")   # Clears registers and error queues
        self.send("TT 2")   # Resets analyser
        
        time.sleep(2)       # Allow enough time for reset to take place

        command_list = ["OS 0",                         # Command seperator is comma
                        "OT 0",                         # Command terminator is "\r\n"
                        "OP 2,1",                       # Send all data in ASCII format
                        "CZ 1",                         # Display impedance coordinates on solartron
                        "UW 1",                         # Display phase normally
                        "IP 1,1",                       # Input V1 is single
                        "OU 1,0",                       # Input V1 is floating
                        "IP 2,1",                       # Input V2 is single
                        "OU 2,0",                       # Input V2 is floating
                        "VB " + str(Vdc),               # D.C bias
                        "DC 1,0",                       # Input V1 coupling is False
                        "DC 3,0",                       # Input current coupling is False
                        "RA 1,0",                       # Input V1 range is 0
                        "IS " + str(integration_time)]  # Integration time

        for command in command_list:
            self.send(command)

    def send(self, message):
        # Send a message to the analyser and delay slightly
        self.analyser.write(message)
        time.sleep(0.1)

    def measure_frequency(self, frequency):

        # Sends an instruction to the Solarton to measures impedance at the specified frequency.
        # Returns impedance data formatted as a comma separated string.
        # First three values are frequency, magnitude and argument.

        fm = "{:.1E}".format(frequency)
        self.send('FR ' + fm)       # Set frequency
        self.send("SI")             # Take single measurement
        
        if frequency < 3:
            time.sleep(1.0 / (frequency + 2))
 
        self.send('SO 1,3')         # Set Display source to 'Z1=V1/I' to convert last measurement to Z, theta   
        self.send('DO')             # Output last result
        
        time.sleep(0.05)
        
        # Read bytes from solartron as ASCII string, and strip trailing zeroes
        result = self.analyser.read().strip()
        return result

    def measure_impedance_process(self, label, ramp, sweep_num, Tcell, experiment_name, temp_interval_num):
        # Send command specifying voltage amplitude for this ramp
        self.send("VA " + str(float(ramp.voltage)))
        
        # If directory for this temperature interval does not exist, create it
        dir_name = "experiments//{}//{}//{}".format(experiment_name, "ramp_{}_{}".format(ramp.num, "up" if ramp.up else "down"), "interval_{}".format(temp_interval_num + 1))
        if not os.path.isdir(dir_name):
            os.mkdir(dir_name)
        
        # Create data file for sweep and write sweep file header
        file_path = "{}//{}".format(dir_name, "sweep_{}.txt".format(sweep_num + 1))
        now = datetime.datetime.now()
        with open(file_path, 'w') as f:
            f.write('sweep_num, date, time, Tcell, setT, ramp_direction\n')
            f.write( '{}, {}, {}, {}, {}, {}\n'.format(sweep_num, now.strftime("%Y-%m-%d"), now.strftime("%H:%M:%S"), Tcell, ramp.end_temp, "up" if ramp.up else "down"))
        
        for i, f in enumerate(ramp.frange):
            label["text"] = "{:.2e} Hz    {} / {}".format(f, i, ramp.numpoints)
            try:
                result = self.measure_frequency(f)
            except:
                result = self.measure_frequency(1.)
                result = self.measure_frequency(f)
            
            # Write result of measurement to file, discarding trailing zeroes
            result = ",".join(result.split(",")[:4])
            with open(file_path, "a") as f:
                f.write(result + "\n")

    def measure_impedance(self, label, ramp, sweep_num, Tcell, experiment_name, temp_interval_num):
        # Begin python thread to measure impedances
        process = threading.Thread(
            target = self.measure_impedance_process,
            args=(label, ramp, sweep_num, Tcell, experiment_name, temp_interval_num),
            daemon = True
            )
        process.start()
        return process
