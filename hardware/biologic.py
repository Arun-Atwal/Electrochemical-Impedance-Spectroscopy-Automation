import ctypes
import sys
import time
import threading
import os
import datetime
from kbio.kbio_api import KBIO_api
from kbio.kbio_tech import ECC_parm, make_ecc_parm, make_ecc_parms
import numpy as np

class eclib(object):

    # Constants
    UNITS_NB = 16

    # EC-Lab DLL (private)
    path = os.path.join(os.path.dirname(__file__), os.pardir)
    ecpath = os.path.join(path, "EC-Lab_Files")
    __dll = ctypes.CDLL(os.path.join(ecpath, "EClib64.dll"))

    # Device Info Structure
    class DeviceInfoType(ctypes.Structure):
        _fields_ = [("DeviceCode", ctypes.c_int32),
                    ("RAMSize", ctypes.c_int32),
                    ("CPU", ctypes.c_int32),
                    ("NumberOfChannels", ctypes.c_int32),
                    ("NumberOfSlots", ctypes.c_int32),
                    ("FirmwareVersion", ctypes.c_int32),
                    ("FirmwareDate_yyyy", ctypes.c_int32),
                    ("FirmwareDate_mm", ctypes.c_int32),
                    ("FirmwareDate_dd", ctypes.c_int32),
                    ("HTdisplayOn", ctypes.c_int32),
                    ("NbOfConnectedPC", ctypes.c_int32)]

    # Current Values Type
    class CurrentValuesType(ctypes.Structure):
        _fields_ = [("State", ctypes.c_int32),
                    ("MemFilled", ctypes.c_int32),
                    ("TimeBase", ctypes.c_float),
                    ("Ewe", ctypes.c_float),
                    ("EweRangeMin", ctypes.c_float),
                    ("EweRangeMax", ctypes.c_float),
                    ("Ece", ctypes.c_float),
                    ("EceRangeMin", ctypes.c_float),
                    ("EceRangeMax", ctypes.c_float),
                    ("Eoverflow", ctypes.c_int32),
                    ("I", ctypes.c_float),
                    ("IRange", ctypes.c_int32),
                    ("Ioverflow", ctypes.c_int32),
                    ("ElapsedTime", ctypes.c_float),
                    ("Freq", ctypes.c_float),
                    ("Rcomp", ctypes.c_float),
                    ("Saturation", ctypes.c_int32),
                    ("OptErr", ctypes.c_int32),
                    ("OptPos", ctypes.c_int32)]

    # Data Information Type
    class DataInfosType(ctypes.Structure):
        _fields_ = [("IRQskipped", ctypes.c_int32),
                    ("NbRows", ctypes.c_int32),
                    ("NbCols", ctypes.c_int32),
                    ("TechniqueIndex", ctypes.c_int32),
                    ("TechniqueID", ctypes.c_int32),
                    ("ProcessIndex", ctypes.c_int32),
                    ("loop", ctypes.c_int32),
                    ("StartTime", ctypes.c_double),
                    ("MuxPad", ctypes.c_int32)]

    # Data buffer Type
    DataBufferType = ctypes.c_uint32 * 1000

    # ECC parameter structure
    class EccParamType(ctypes.Structure):
        _fields_ = [("ParamStr", 64 * ctypes.c_byte),
                    ("ParamType", ctypes.c_int32),
                    ("ParamVal", ctypes.c_uint32),
                    ("ParamIndex", ctypes.c_int32)]

    # ECC parameters structure
    class EccParamsType(ctypes.Structure):
        _fields_ = [("len", ctypes.c_int32),
                    ("pParams", ctypes.c_void_p)]

    # Array of units
    UnitsType = ctypes.c_byte * UNITS_NB

    # Array of results
    ResultsType = ctypes.c_int32 * UNITS_NB

    # Error Enumeration
    class ErrorCodeEnum(object):
        ERR_NOERROR = 0

    # Technique Parameter Type Enumeration
    class ParamTypeEnum(object):
        PARAM_INT = 0
        PARAM_BOOLEAN = 1
        PARAM_SINGLE = 2

    BL_GetUSBdeviceinfos = __dll["BL_GetUSBdeviceinfos"]
    BL_GetUSBdeviceinfos.restype= ctypes.c_bool

    # ErrorCode BL_ConvertNumericIntoSingle(int num, ref float psgl)
    BL_ConvertNumericIntoSingle = __dll["BL_ConvertNumericIntoSingle"]
    BL_ConvertNumericIntoSingle.restype = ctypes.c_int

    # ErrorCode BL_Connect(string server, byte timeout, ref int connection_id, ref DeviceInfo pInfos)
    BL_Connect = __dll["BL_Connect"]
    BL_Connect.restype = ctypes.c_int
    
    # ErrorCode BL_Disconnect(int ID)
    BL_Disconnect = __dll["BL_Disconnect"]
    BL_Disconnect.restype = ctypes.c_int
    
    # ErrorCode BL_TestConnection(int ID)
    BL_TestConnection = __dll["BL_TestConnection"]
    BL_TestConnection.restype = ctypes.c_int

    # ErrorCode BL_LoadFirmware(int ID, byte[] pChannels, int[] pResults, byte Length, bool ShowGauge, bool ForceReload, string BinFile, string XlxFile)
    BL_LoadFirmware = __dll["BL_LoadFirmware"]
    BL_LoadFirmware.restype = ctypes.c_int

    # bool BL_IsChannelPlugged(int ID, byte ch)
    BL_IsChannelPlugged = __dll["BL_IsChannelPlugged"]
    BL_IsChannelPlugged.restype = ctypes.c_bool

    # ErrorCode BL_GetChannelsPlugged(int ID, byte[] pChPlugged, byte Size)
    BL_GetChannelsPlugged = __dll["BL_GetChannelsPlugged"]
    BL_GetChannelsPlugged.restype = ctypes.c_int

    # ErrorCode BL_GetMessage(int ID, byte ch, [MarshalAs(UnmanagedType.LPArray)] byte[] msg, ref int size)
    BL_GetMessage = __dll["BL_GetMessage"]
    BL_GetMessage.restype = ctypes.c_int

    # ErrorCode BL_LoadTechnique(int ID, byte channel, string pFName, EccParams pparams, bool FirstTechnique, bool LastTechnique, bool DisplayParams)
    BL_LoadTechnique = __dll["BL_LoadTechnique"]
    BL_LoadTechnique.restype = ctypes.c_int

    # ErrorCode BL_DefineSglParameter(string lbl, float value, int index, IntPtr pParam)
    BL_DefineSglParameter = __dll["BL_DefineSglParameter"]
    BL_DefineSglParameter.restype = ctypes.c_int

    # ErrorCode BL_DefineSglParameter(string lbl, int value, int index, IntPtr pParam)
    BL_DefineIntParameter = __dll["BL_DefineIntParameter"]
    BL_DefineIntParameter.restype = ctypes.c_int

    # ErrorCode BL_DefineSglParameter(string lbl, bool value, int index, IntPtr pParam)
    BL_DefineBoolParameter = __dll["BL_DefineBoolParameter"]
    BL_DefineBoolParameter.restype = ctypes.c_int

    # ErrorCode BL_TestCommSpeed(int ID, byte channel, ref int spd_rcvt, ref int spd_kernel)
    BL_TestCommSpeed = __dll["BL_TestCommSpeed"]
    BL_TestCommSpeed.restype = ctypes.c_int

    # ErrorCode BL_StartChannel(int ID, byte channel)
    BL_StartChannel = __dll["BL_StartChannel"]
    BL_StartChannel.restype = ctypes.c_int
    
    # ErrorCode BL_StopChannel(int ID, byte channel)
    BL_StopChannel = __dll["BL_StopChannel"]
    BL_StopChannel.restype = ctypes.c_int

    # ErrorCode BL_GetData(int ID, byte channel, [MarshalAs(UnmanagedType.LPArray, SizeConst=1000)] int[] buf, ref DataInfos pInfos, ref CurrentValues pValues)
    BL_GetData = __dll["BL_GetData"]
    BL_GetData.restype = ctypes.c_int

class SP_200(eclib):
    def __init__(self):
        self.name = "Biologic SP-200"
        
        # Hardware specifications
        self.min_frequency = 1e-5
        self.max_frequency = 3e6
        self.min_V = 1e-6 #micro volt, although data will likely be terrible at this amplitude
        self.max_V = 1e4 # 10V (10,000 mV), absolute maximum but not entirely recomended
        
        # Configuration
        self.cfg_conn_ip = b"USB0"
        self.cfg_conn_timeout = 10
        self.cfg_channel = 0
        self.cfg_debug_enabled = False

        # global variables (in SP_300 scope)
        self.glob_firmware_loaded = False
        self.glob_stop = False
        self.glob_conn_id = ctypes.c_int(-1)

        # Connect to instrument (call BL_CONNECT)
        device_info = self.DeviceInfoType()
        error = self.BL_Connect(ctypes.c_char_p(self.cfg_conn_ip), ctypes.c_byte(self.cfg_conn_timeout), ctypes.byref(self.glob_conn_id), ctypes.byref(device_info))
        if error != self.ErrorCodeEnum.ERR_NOERROR:
            raise Exception("Error connecting to instrument. Errcode = {}".format(error))
        
        # Get connected channels
        units = self.UnitsType()
        error = eclib.BL_GetChannelsPlugged(self.glob_conn_id, ctypes.byref(units), ctypes.c_ubyte(self.UNITS_NB))
        if error != self.ErrorCodeEnum.ERR_NOERROR:
            raise Exception("Error retrieving connected channel(s). Errcode = {}".format(error))
        
        # Load firmware
        results = self.ResultsType()
        error = self.BL_LoadFirmware(self.glob_conn_id, units, ctypes.byref(results), ctypes.c_ubyte(self.UNITS_NB), False, True, None, None)
        if error != self.ErrorCodeEnum.ERR_NOERROR:
            raise Exception ("Error loading firmware. Errcode = {}".format(error))
        self.glob_firmware_loaded = True
    
    def disconnect(self):
        error = self.BL_StopChannel(ctypes.c_int32(0), 0)
        if error != self.ErrorCodeEnum.ERR_NOERROR:
            raise Exception ("Error closing channel. Errcode = {}".format(error))
        error = self.BL_Disconnect(ctypes.c_int32(0))
        if error != self.ErrorCodeEnum.ERR_NOERROR:
            raise Exception ("Error disconnecting. Errcode = {}".format(error))

    def define_parameter(self, p_type, p_index, label, value, index):
        if p_type == "float":
            if not isinstance(value, float):
                return False
            error = self.BL_DefineSglParameter(ctypes.c_char_p(label), ctypes.c_float(value), ctypes.c_int32(index), ctypes.byref(self.EccParamArray[p_index]))
            if error != self.ErrorCodeEnum.ERR_NOERROR:
                return False
            return True
        elif p_type == "int":
            if not isinstance(value, int):
                return False
            error = self.BL_DefineIntParameter(ctypes.c_char_p(label), ctypes.c_int32(value), ctypes.c_int32(index), ctypes.byref(self.EccParamArray[p_index]))
            if error != self.ErrorCodeEnum.ERR_NOERROR:
                return False
            return True
        elif p_type == "bool":
            if not isinstance(value, bool):
                return False
            error = self.BL_DefineBoolParameter(ctypes.c_char_p(label), ctypes.c_bool(value), ctypes.c_int32(index), ctypes.byref(self.EccParamArray[p_index]))
            if error != self.ErrorCodeEnum.ERR_NOERROR:
                return False
            return True
        else:
            return False

    def measure_impedance_process(self, label, ramp, sweep_num, Tcell, experiment_name, temp_interval_num):
        # If directory for this temperature interval does not exist, create it
        dir_name = "experiments//{}//{}//{}".format(experiment_name, "ramp_{}_{}-{}".format(ramp.num, ramp.start_temp, ramp.end_temp), "interval_{}°C".format(round(Tcell)))
        if not os.path.isdir(dir_name):
            os.mkdir(dir_name)
        
        # Create data file for sweep and write sweep file header
        file_path = "{}//{}".format(dir_name, "sweep_{}_{}°C.txt".format(sweep_num + 1, Tcell))
        now = datetime.datetime.now()
        # All header data as previously included below. Causes ZView to skip the first lines of data due to format.
        # with open(file_path, 'w') as f:
        #     f.write('sweep_num, date, time, Tcell, setT, ramp_direction\n')
        #     f.write( '{}, {}, {}, {}, {}, {}\n'.format(sweep_num, now.strftime("%Y-%m-%d"), now.strftime("%H:%M:%S"), Tcell, ramp.end_temp, "up" if ramp.up else "down"))
        with open(file_path, 'w') as f:
            f.write('Experiment name: {}\n'.format(experiment_name))
            f.write('{}\n'.format(now.strftime("%Y-%m-%d")))
            f.write('{}\n'.format(now.strftime("%H:%M:%S")))
            f.write('Temperature interval number {} in ramp from {}°C to {}°C\n'.format(temp_interval_num+1, ramp.start_temp, ramp.end_temp))
            f.write('Sweep number: {}\n'.format(sweep_num+1))
            f.write('Current Temperature: {}\n'.format(Tcell))
            f.write('freq/Hz\tRe(Z)/Ohm\tIm(Z)/Ohm\n')
        
        # Set up parameters
        params = [
            ["float", 0, b"Final_frequency", ramp.fmin, 0],
            ["float", 1, b"Initial_frequency", ramp.fmax, 0],
            ["bool", 2, b"sweep", False, 0], #logarithmic sweep
            ["float", 3, b"Amplitude_Voltage", ramp.voltage, 0], 
            ["int", 4, b"Frequency_number", ramp.numpoints, 0],
            ["int", 5, b"Average_N_times", 1, 0],
            ["bool", 6, b"Correction", False, 0],
            ["float", 7, b"Wait_for_steady", 0.1, 0]
        ]

        path = os.path.join(os.path.dirname(__file__), os.pardir)
        ecpath = os.path.join(path, "EC-Lab_Files")
        dll_path = os.path.join(ecpath, "EClib64.dll")
        api = KBIO_api(dll_path)
        
        """IMPORTANT: Below is a result of issues with ctypes and the 64 bit environment required for the Linkam stages.
        Here, the module layering over the dll provided by biologic, kbio_api, and related scripts are used to handle the ctypes and other intricacies involved with different bit regimes
        However, the entirety of the rest of the code, aside from the loading of ecc techniques, worked fine doing this manually.
        Perhaps with a more careful handling of ctypes here, the issue could have been avoided. However, this set of scripts certainly works.
        Instead of entirely rewriting the code, I have opted to use the api only for this function, hence the vast incosistency.
        If trying to improve the code, rewriting it to entirely use the api script would be a good place to start. 
        However, for the time being, it does not seem too difficult to do everything as is done here, aside from loading any other techniques with the api as below."""

        api_params = {
            'ramp.fmin' : ECC_parm("Final_frequency", float),
            'ramp.fmax' : ECC_parm("Initial_frequency",float),
            'sweep' : ECC_parm("sweep", bool),
            'ramp.voltage' : ECC_parm("Amplitude_Voltage", float),
            'ramp.numpoints' : ECC_parm("Frequency_number", int),
            'average' : ECC_parm("Average_N_times", int),
            'correction' : ECC_parm("Correction", bool),
            'wait' : ECC_parm("Wait_for_steady", float),
        }

        p_finalf = make_ecc_parm(api, api_params['ramp.fmin'], ramp.fmin)
        p_initf = make_ecc_parm(api, api_params['ramp.fmax'], ramp.fmax)
        p_sweep = make_ecc_parm(api, api_params['sweep'], False)
        p_voltage = make_ecc_parm(api, api_params['ramp.voltage'], ramp.voltage)
        p_numpoints = make_ecc_parm(api, api_params['ramp.numpoints'], ramp.numpoints)
        p_average = make_ecc_parm(api, api_params['average'], 1)
        p_correction = make_ecc_parm(api, api_params['correction'], False)
        p_wait = make_ecc_parm(api, api_params['wait'], 0.1)
        ecc_parms = make_ecc_parms(api, p_finalf, p_initf, p_sweep, p_voltage, p_numpoints, p_average, p_correction, p_wait)
        
        params_nb = len(params)
        EccParamArrayType = self.EccParamType * params_nb
        self.EccParamArray = EccParamArrayType()
        self.EccParams = self.EccParamsType()
        self.EccParams.len = ctypes.c_int32(params_nb)
        self.EccParams.pParams = ctypes.cast(self.EccParamArray, ctypes.c_void_p)
        
        # Pass parameters
        for param in params:
            success = self.define_parameter(param[0], param[1], param[2], param[3], param[4])
            if success == False:
                raise Exception("Error defining parameter: {}".format(param[2]))
            time.sleep(0.2)

        # Load technique
        # peis = ctypes.create_string_buffer(b"peis4.ecc")
        # error = self.BL_LoadTechnique(self.glob_conn_id, ctypes.c_ubyte(self.cfg_channel), ctypes.byref(peis), self.EccParams, True, True, True)

        api.LoadTechnique(self.glob_conn_id, self.cfg_channel, "peis4.ecc", ecc_parms, first=True, last=True, display=False)
        # if error != self.ErrorCodeEnum.ERR_NOERROR:
        #     raise Exception("BL_LoadTechnique error. Errcode = {}".format(error))
        
        # Start channel
        error = self.BL_StartChannel(self.glob_conn_id, ctypes.c_ubyte(self.cfg_channel))
        if error != self.ErrorCodeEnum.ERR_NOERROR:
            raise Exception("BL_StartChannel error. Errcode = {}".format(error))

        # Main loop
        for i in range(ramp.numpoints):
            label["text"] = "{:.2e} Hz    {} / {}".format(ramp.frange[i], i, ramp.numpoints)
            result, PI = self.get_result()
            with open(file_path, "a") as f:
                if PI == 1: # if data is from process index 1 (process of interest), store the data
                    f.write("\t".join(str(x) for x in result) + "\n")

    
    def get_result(self):
        # Process and return data
        data = []
        while not data: # if data is empty, try again
            # Retrieve data
            buffer = self.DataBufferType()
            infos = self.DataInfosType()
            values = self.CurrentValuesType()
            error = self.BL_GetData(self.glob_conn_id, ctypes.c_ubyte(self.cfg_channel), ctypes.byref(buffer), ctypes.byref(infos), ctypes.byref(values))
            if error != self.ErrorCodeEnum.ERR_NOERROR:
                raise Exception("BL_GetData error. Errcode = {}".format(error))
            
            if infos.ProcessIndex == 1:
                for i in range(infos.NbRows * infos.NbCols):
                    receive_data = ctypes.c_float(0.0)
                    error = self.BL_ConvertNumericIntoSingle(buffer[i], ctypes.byref(receive_data))
                    if error != self.ErrorCodeEnum.ERR_NOERROR:
                        raise Exception("BL_ConvertNumericIntoSingle error. Errcode = {}".format(error))
                    data.append(receive_data.value)

        # below: freq, Re = V/I * cos(arg(z)), Im = V/I * sin(arg(z))
        result = [data[0],abs(data[1]) / abs(data[2]) * np.cos(data[3]),abs(data[1]) / abs(data[2]) * np.sin(data[3])]
        # Conversion to degrees included - depends on what the expected data is
        return result, int(infos.ProcessIndex)
        
    def measure_impedance(self, label, ramp, sweep_num, Tcell, experiment_name, temp_interval_num):
        process = threading.Thread(
            target = self.measure_impedance_process,
            args = (label, ramp, sweep_num, Tcell, experiment_name, temp_interval_num),
            daemon = True
        )
        process.start()
        return process