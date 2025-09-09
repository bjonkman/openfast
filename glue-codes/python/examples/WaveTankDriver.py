
from ctypes import (
    CDLL,
    POINTER,
    create_string_buffer,
    byref,
    c_byte,
    c_int,
    c_double,
    c_float, 
    c_char,
    c_char_p, 
    c_bool
)

import numpy as np
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
from scipy.spatial.transform import Rotation

from pyOpenFAST.interface_abc import OpenFASTInterfaceType
from pyOpenFAST.tdmslib import parse_tdms

project_root = '../../../'
library_path = project_root + '/build-Single-Debug/glue-codes/labview/libwavetanktestinglib.dylib'

class WaveTankLib(OpenFASTInterfaceType):

    #--------------------------------------
    # Error levels
    #--------------------------------------
    error_levels: Dict[int, str] = {
        0: "None",
        1: "Info",
        2: "Warning",
        3: "Severe Error",
        4: "Fatal Error"
    }

    #--------------------------------------
    # Constants
    #--------------------------------------
    # NOTE: The length of the error message in Fortran is determined by the
    #       ErrMsgLen variable in the NWTC_Base.f90 file. If ErrMsgLen is modified,
    #       the corresponding size here must also be updated to match.
    ERROR_MESSAGE_LENGTH: int = 8197
    DEFAULT_STRING_LENGTH: int = 1025

    def __init__(self, library_path: str, input_file_names: dict):
        """

        Args:
            library_path (str): Path to the compile wavetank interface shared library
            input_file_names (dict): Map of file names for each included module:
                - WT_InputFile
                - MD_InputFile
                - SS_InputFile
                - AD_InputFile
                - IfW_InputFile
        """
        super().__init__(library_path)

        # Create C-compatible string buffers for input file names
        self.input_file_names = {
            k: create_string_buffer(str(Path(v).absolute()).encode('utf-8'), self.IntfStrLen)
            for k,v in input_file_names.items()
        }

        self._initialize_routines()

        self.ended = False   # For error handling at end
        self.print_error_level = 1

        # Create buffers for class data
        # These will generally be overwritten by the Fortran code
        # self.ss_output_channel_names = []
        # self.ss_output_channel_units = []
        # self.ss_output_values = None

        self.md_output_channel_names = []
        self.md_output_channel_units = []
        self.md_output_values = None

        self.adi_output_channel_names = []
        self.adi_output_channel_units = []
        self.adi_output_values = None

        # Error handling setup
        self.abort_error_level = 4
        self.error_status_c = c_int(0)
        self.error_message_c = create_string_buffer(self.ERROR_MESSAGE_LENGTH)

    def _initialize_routines(self):
        self.WaveTank_Init.argtypes = [
            POINTER(c_char),        #  intent(in   ) :: WT_InputFile_c(IntfStrLen)
            POINTER(c_char),        #  intent(in   ) :: MD_InputFile_c(IntfStrLen)
            POINTER(c_char),        #  intent(in   ) :: SS_InputFile_c(IntfStrLen)
            POINTER(c_char),        #  intent(in   ) :: AD_InputFile_c(IntfStrLen)
            POINTER(c_char),        #  intent(in   ) :: IfW_InputFile_c(IntfStrLen)
            POINTER(c_int),         #  intent(  out) :: ErrStat_C
            POINTER(c_char),        #  intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.WaveTank_Init.restype = c_int

        self.WaveTank_CalcOutput.argtypes = [
            POINTER(c_double),      # real(c_double) :: time
            POINTER(c_float),       # real(c_float),          intent(in   ) :: positions_x
            POINTER(c_float),       # real(c_float),          intent(in   ) :: positions_y
            POINTER(c_float),       # real(c_float),          intent(in   ) :: positions_z
            POINTER(c_float),       # real(c_float),          intent(in   ) :: floater_rotation_matrix(9)
            POINTER(c_float),       # real(c_float),          intent(in   ) :: blade_rotation_matrix(9)
            POINTER(c_float),       # real(c_float),          intent(  out) :: MD_Forces_C(1,6)
            POINTER(c_float),       # real(c_float),          intent(  out) :: ADI_MeshFrc_C(NumMeshPts,6)
            POINTER(c_float),       # real(c_float),          intent(  out) :: ADI_HHVel_C(3)
            POINTER(c_float),       # real(c_float),          intent(  out) :: md_outputs(MD_NumChannels_C)
            POINTER(c_float),       # real(c_float),          intent(  out) :: adi_outputs(ADI_NumChannels_C)
            POINTER(c_int),         # integer(c_int),         intent(  out) :: ErrStat_C
            POINTER(c_char),        # character(kind=c_char), intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.WaveTank_CalcOutput.restype = c_int

        self.WaveTank_End.argtypes = [
            POINTER(c_int),         # integer(c_int),         intent(  out) :: ErrStat_C
            POINTER(c_char),        # character(kind=c_char), intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.WaveTank_End.restype = c_int

        self.WaveTank_SetWaveFieldPointer.argtypes = [
            POINTER(c_int),         # integer(c_int),         intent(  out) :: ErrStat_C
            POINTER(c_char),        # character(kind=c_char), intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.WaveTank_SetWaveFieldPointer.restype = c_int

        self.WaveTank_Sizes.argtypes = [
            POINTER(c_int),
            POINTER(c_int),
            POINTER(c_int),
        ]
        self.WaveTank_Sizes.restype = c_int


    def check_error(self) -> None:
        """Checks for and handles any errors from the Fortran library.

        Raises:
            RuntimeError: If a fatal error occurs in the Fortran code
        """
        # If the error status is 0, return
        if self.error_status_c.value == 0:
            return

        # Get the error level and error message
        error_level = self.error_levels.get(
            self.error_status_c.value,
            f"Unknown Error Level: {self.error_status_c.value}"
        )
#FIXME: losing the second line of the message here!!!!
        error_msg = self.error_message_c.value.decode('utf-8') #.strip()
        message = f"WaveTank library {error_level}: {error_msg}"
        # If the error level is fatal, call WaveTank_End() and raise an error
        if self.error_status_c.value >= self.abort_error_level:
            print(f"Fatal error occured.")
            print(message)
            raise RuntimeError(message)
            try:
                self.end()
            except Exception as e:
                message += f"\nAdditional error during cleanup: {e}"
#FIXME: why isn't this triggered???

            raise RuntimeError(message)
        else:
            print(message)


    def init(self):
        _error_message = create_string_buffer(self.ERROR_MSG_C_LEN)

        # # Convert the initial positions array into c_float array
        # init_positions_c = (c_float * 6)(0.0, )
        # for i, p in enumerate(platform_init_pos):
        #     init_positions_c[i] = c_float(p)

        self.WaveTank_Init(
            self.input_file_names["WaveTank"],
            self.input_file_names["MoorDyn"],
            self.input_file_names["SeaState"],
            self.input_file_names["AeroDyn"],
            self.input_file_names["InflowWind"],
            byref(self.error_status_c),             # OUT <- error status code
            self.error_message_c                    # OUT <- error message buffer
        )
        self.check_error()

        # self.output_channel_names = [n.decode('UTF-8') for n in _channel_names.value.split()] 
        # self.output_channel_units = [n.decode('UTF-8') for n in _channel_units.value.split()] 
        # self.output_values = np.zeros( self.num_outs_c.value, dtype=c_float, order='C' )

    def calc_output(
        self,
        time: float,
        positions_x: float,
        positions_y: float,
        positions_z: float,
        floater_rotation_matrix: np.ndarray,
        blade_rotation_matrix: np.ndarray,
        md_loads: np.ndarray,
        ad_loads: np.ndarray,
        hub_height_velocities: np.ndarray,
    ):

        self.WaveTank_CalcOutput(
            byref(c_double(time)),
            byref(c_float(positions_x)),
            byref(c_float(positions_y)),
            byref(c_float(positions_z)),
            floater_rotation_matrix.ctypes.data_as(POINTER(c_float)),
            blade_rotation_matrix.ctypes.data_as(POINTER(c_float)),
            md_loads.ctypes.data_as(POINTER(c_float)),
            ad_loads.ctypes.data_as(POINTER(c_float)),
            hub_height_velocities.ctypes.data_as(POINTER(c_float)),
            self.md_output_values.ctypes.data_as(POINTER(c_float)),
            self.adi_output_values.ctypes.data_as(POINTER(c_float)),
            byref(self.error_status_c),             # OUT <- error status code
            self.error_message_c                    # OUT <- error message buffer
        )
        self.check_error()

    def end(self) -> None:
        _error_message = create_string_buffer(self.ERROR_MSG_C_LEN)

        self.WaveTank_End(
            byref(self.error_status_c),             # OUT <- error status code
            self.error_message_c                    # OUT <- error message buffer
        )
        self.check_error()

    def allocate_outputs(self):
        ss_numouts = c_int(0)
        md_numouts = c_int(0)
        adi_numouts = c_int(0)
        self.WaveTank_Sizes(
            byref(ss_numouts),
            byref(md_numouts),
            byref(adi_numouts),
        )

        # self.ss_output_values = np.zeros(ss_numouts.value, dtype=np.float32, order='C')
        self.md_output_values = np.zeros(md_numouts.value, dtype=np.float32, order='C')
        self.adi_output_values = np.zeros(adi_numouts.value, dtype=np.float32, order='C')
        # self.ss_output_channel_names = [b""] * ss_numouts.value
        # self.ss_output_channel_units = [b""] * ss_numouts.value
        # self.md_output_channel_names = [b""] * md_numouts.value
        # self.md_output_channel_units = [b""] * md_numouts.value
        # self.adi_output_channel_names = [b""] * adi_numouts.value
        # self.adi_output_channel_units = [b""] * adi_numouts.value


if __name__=="__main__":

    floater_motions = parse_tdms("Full_10N_Wrench_1.tdms")
    n_timesteps = len(floater_motions["x"])

    wavetanklib = WaveTankLib(
        library_path,
        {
            "WaveTank":   "wavetankconfig.in",
            "MoorDyn":    "../../../reg_tests/r-test/glue-codes/openfast/MHK_RM1_Floating/MHK_RM1_Floating_MoorDyn.dat",
            "SeaState":   "../../../reg_tests/r-test/glue-codes/openfast/MHK_RM1_Floating/SeaState.dat",
            "AeroDyn":    "../../../reg_tests/r-test/glue-codes/openfast/MHK_RM1_Floating/MHK_RM1_Floating_AeroDyn.dat",
            "InflowWind": "../../../reg_tests/r-test/glue-codes/openfast/MHK_RM1_Floating/MHK_RM1_Floating_InflowWind.dat",
        },
    )
    wavetanklib.init()

    rotation_matrix = np.eye(3, 3, dtype=np.float32)
    md_loads = np.zeros((1,6), dtype=np.float32, order='C')
    ad_loads = np.zeros((2,6), dtype=np.float32, order='C')
    hub_height_velocities = np.zeros((3,1), dtype=np.float32, order='C')

    wavetanklib.allocate_outputs()

    blade_dcm = np.zeros((2*9), dtype=np.float32, order='C')

    for i in range(n_timesteps):
        R = Rotation.from_euler(
            "xyz",
            (
                floater_motions["phi"][i],      # roll
                floater_motions["theta"][i],    # pitch
                floater_motions["psi"][i],      # yaw
            )
        )
        floater_dcm = R.as_matrix().flatten()

        # Create the rotation matrix for the blades using the loop index as a rotation angle
        # The second blade is rotated 180 degrees from the first
        blade_dcm[0:9] = Rotation.from_euler("xyz", (np.deg2rad(i), 0.0, 0.0)).as_matrix().flatten()
        blade_dcm[9:18] = Rotation.from_euler("xyz", (np.deg2rad(i + 180), 0.0, 0.0)).as_matrix().flatten()

        wavetanklib.calc_output(
            time=i,
            positions_x=floater_motions["x"][i],
            positions_y=floater_motions["y"][i],
            positions_z=floater_motions["z"][i],
            floater_rotation_matrix=floater_dcm,
            blade_rotation_matrix=blade_dcm,
            md_loads=md_loads,
            ad_loads=ad_loads,
            hub_height_velocities=hub_height_velocities
        )
        # print(hub_height_velocities)

        # print(wavetanklib.md_output_values)
    # print(wavetanklib.md_output_values)

    wavetanklib.end()
