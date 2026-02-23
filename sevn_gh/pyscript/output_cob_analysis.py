#!/usr/bin/env python
"""
SCRIPT: output_cob_analysis.py
PURPOSE: Take all the output of a sevn simulation filter compact object binaries
and save the following files:
	- BHBH.csv: containing all the bound BHBH binaries
	- BHBHm.csv: containing all the bound BHBH binaries that will merge within an Hubble time (14000 My)
	- NSNS.csv: containing all the bound NSNS binaries
	- NSNSm.csv: containing all the bound NSNS binaries that will merge within an Hubble time (14000 My)
	- BHNS.csv: containing all the bound BHNS binaries
	- BHNSm.csv: containing all the bound BHNSm binaries that will merge within an Hubble time (14000 My)

There are two ways to use this script,
****1- Analyse a single sevn output folder*********
	USAGE;  ./output_cob_analysis.py [-h] [-n NPROC] [-o OUT] input
	e.g. ./output_cob_analysis.py sevn_output -n 2
	positional arguments:
	  input                 path to the sevn folder with the input to analyse

	optional arguments:
	  -h, --help            show this help message and exit
	  -n NPROC, --nproc NPROC
							Number of parallel processes to use [2]
	  -o OUT, --output OUT  Output folder where to store the files, if None use the same input folder [None]
*******************************************************************

****2- Analyse multiple sevn output folders*********
	usage: ./output_cob_analysis.py -m [-h] [-n NPROC]   [--froot FROOT] [--subfolder SUBFOLDER] input
	e.g. ./output_cob_analysis.py -m   -n 2  --froot sevn_output --subfolder 0 main_folder
		In this case, the script will analyse all the data folder starting with  name starting with
		sevn_output inside the main_folder, the output data are inside the subfolder with name 0, e.g.
		sevn_output_Z02/0, the output are automatically saved inside the root sevn output folders
		(e.g. sevn_output_Z02)

	positional arguments:
	  input                 root folder containing all the sevn output folders to analyse

	optional arguments:
	  -h, --help            show this help message and exit
	  -n NPROC, --nproc NPROC
							Number of parallel processes to use [2]
	  -m, --multifolders    If true analyse a list of output folders inside the path specified
	  --froot FROOT         common root for multiple sevn output analysis, used only if option -m is enabled [sevn_output]
	  --subfolder SUBFOLDER
							Subfolder containg the SEVN output, used only if option -m is enabled [None]on input and with folder root specified with
							the froot parameter [False]

REQUIRED PACKAGE:
    -numpy
    -pandas

V 1.0: 04/01/22 Giuliano Iorio: First version
V 1.1: 05/01/22 Giuliano Iorio: Added the option to analyse multiple sevn output folders
V 1.2: 13/04/22 Giuliano Iorio: Changed the function get_history to have  more information in the summary labels
				+ added the possibility to use proper tGW estimate if the pyblack module is available
V 1.3: 14/04/22 Giuliano Iorio: Changed the function get_history to have  more information in the summary labels for kicks
V 1.4: 12/07/22 Giuliano Iorio: Update the estimate of tgw based on the update of pyblack
V 1.5: 28/07/22 Giuliano Iorio: Added the namespace-class channelify and added channel information in the output
V 1.6: 07/02/24 Giuliano Iorio: Many changes
								- changed the name of the method get_history in get_loginfo since now it returns more than just the history
							    - the method get_longinfo now returns a DataFrame containing the list of the events but also details about the SN explosions
								- more robust regex pattern identification in get_loginfo
								- added a try..exception check when running the analysis
"""

import numpy as np
import pandas as pd
import re
import glob
from itertools import repeat
from multiprocessing import Pool
import argparse
import os
import sys

try:
	from pyblack.gw.gw_time import estimate_tgw
	tgw_loaded=True
except:
	print("Warning module pyblack not available, the delay time will be estimated using the GWtime column from the SEVN ouput. pyblack can be installed with pip install pyblack.")
	tgw_loaded=False

class Channel:

	_available_channels    = (-1,0,1,2,3,33,4,44,5)
	_available_subchannels = ("n","s","c","d")

	@staticmethod
	def summarise_channels(df: pd.DataFrame) -> dict:
		"""
        Create a dictonary containing a summary of the statiscs of different channels subchannels and  other
        cases (i.e stripped SN, ultra stripped SN, kollisions)
        :param df: Pandas dataframe, it has to containt the column Eventplus, if it does not contain the column
        channel and subchannel they will be created using the function channellify.
        :return: A dictionary with format:
            - channel : (Ntot,fraction, dictsubchannel, dictspecial)
            where
                - Ntot is the total number of object in the specific channel
                - fraction is the fraction of these obeject with respect the total number  of object in the dataframe
                - dictsubchannel: a dictionary with format
                    - subchannel : frcation of the objcet in the subchannel with respect to the total number of object i the channel
                - dictspecial: a dictionary with format
                    - St: fraction of stripped (no Hydrogen envelope) SN explosion among all the SN in the channel
                    - Su: fraction of ultrastripped (no Hydrogen/Helium envelope) SN explosion among all the SN in the channel
                    - K: fraction of systems that trigger at leat a collision at periastron  among all the binary in the channel
        """
		Ntot=len(df)
		outdict={}
		for channel in Channel._available_channels:
			dft=Channel.filter_channel(df,channel=channel)
			Nch=len(dft)
			fch= 0 if Ntot==0 else Nch/Ntot
			subch={}
			for subchannel in Channel._available_subchannels:
				dftt = Channel.filter_channel(dft, subchannel=subchannel)
				Nsch=len(dftt)
				fsch= 0 if Nch==0 else Nsch/Nch
				subch[subchannel]=fsch
			extra_dic={}
			if "NSt" in dft.columns:  extra_dic["St"] = 0 if Nch==0 else np.sum(dft["NSt"])/(2*Nch) #2 times because there are two explosions per binary
			if "NSu" in dft.columns:  extra_dic["Su"] = 0 if Nch==0 else np.sum(dft["NSu"])/(2*Nch) #2 times because there are two explosions per binary
			if "collisions" in df.columns: extra_dic["K"] = 0 if Nch==0  else np.sum(dft.collisions>0)/Nch
			outdict[channel] = (Nch, fch, subch,extra_dic)

		return outdict


	@staticmethod
	def filter_channel(df: pd.DataFrame, channel: int=None, subchannel: str=None) -> pd.DataFrame:
		"""
        Filter a dataframe based on channels and subchannels
        :param df:pandas dataframe containing the column EventsPlus, it should contain also
        the columns channel and subchannel, but if they are not present they will be automatically added
        :param channel: formation channel to filter
            0- No interaction at all
            1- Only SMT before the first SN explosion
            2- Only SMT
            3- At least 1 CE before the first SN explosion, at the first SN explosion one star is stripped by the envelope
                - 33 Exactly 1 CE happens before the first SN between a star with a developed core and a MS
            4- At least 1 CE before the first SN explosion, at the first SN explosion both stars are stripped by the envelope
                - 44 Exactly 1 CE happens before the first SN, both stars in the CE have a core
            5- No interactions before the first SN explotion
            -1- All the other classes (Btw, the class before should contains all the possible cases)
        :param subchannel: formation subchannel to filter
            n- No interaction after the first SN explosion
            s- Only stable mass transfers after the first SN explosion
            c- Exactly 1 CE after the first SN explosion
            d- More than 1 CE after the first SN explosion
        :return: a filtered dataframe, if the columns channel and subchannel were not present the dataframe
        will contains also these new columns plus additional one (see the funciton channellify). The only exception
        is when both channel and subchannel is None, in this case the original daframe will be returned.
        """

		# If both filters None just return df
		if channel is None and subchannel is None:
			return df
		elif channel  is not  None and channel not in Channel._available_channels:
			raise ValueError(f"Channel {channel} not in the list of available channels {Channel._available_channels}")
		elif subchannel  is not None and subchannel not in Channel._available_subchannels:
			raise ValueError(f"Subchannel {subchannel} not in the list of available channels {Channel._available_subchannels}")


		# Check columns
		# If channel and subchannel not in the dataframe craete the columns
		if "channel"  not in df.columns or "subchannel"  not in df.columns:
			df=Channel.channellify(df,inplace=False)

		if channel is not None:
			# If channel 3 consider also the 33 case
			if channel==3: idxch=(df.channel==3) | (df.channel==33)
			# If channel 4 consider also the 44 case
			elif channel==4: idxch=(df.channel==4) | (df.channel==44)
			else: idxch=df.channel==channel
		else:
			idxch=np.array([True,]*len(df)) #No filters

		if subchannel is not None:
			idxsc=df.subchannel==subchannel
		else:
			return df[idxch]

		return df[idxch&idxsc]

	@staticmethod
	def channellify(df: pd.DataFrame,inplace: bool = True) -> pd.DataFrame:
		"""
        Add columns to the dataframe storing the channel, the number of collisions and number of CE events
        :param df:  pandas dataframe containing the column EventsPlus
        :param inplace:  if False make a copy of the dataframe
        :return: the input (or a copy, depening on parameter in place) dataframe with the following new columns:
            channel: formation channel (-1,0,1,3,4,5)
                        0- No interaction at all
                        1- Only SMT before the first SN explosion
                        2- Only SMT
                        3- At least 1 CE before the first SN explosion, at the first SN explosion one star is stripped by the envelope
                            - 33 Exactly 1 CE happens before the first SN between a star with a developed core and a MS
                        4- At least 1 CE before the first SN explosion, at the first SN explosion both stars are stripped by the envelope
                            - 44 Exactly 1 CE happens before the first SN, both stars in the CE have a core
                        5- No interactions before the first SN explotion
                        -1- All the other classes (Btw, the class before should contains all the possible cases)
            subchannel: sub channel based on the history after the second SN explosion
                        n- No interaction after the first SN explosion
                        s- Only stable mass transfers after the first SN explosion
                        c- Exactly 1 CE after the first SN explosion
                        d- More than 1 CE after the first SN explosion
            collisions: number of triggered collisions at periastron
            CE: total number of triggered CE events
            CEb: number of triggered CE events before the explosion of the first SN
            CEa: number of triggered CE events after the explosion of the first SN
            SMTa: if 0 not stable mass transfers after the first SN, otherwise 1
            NSt: number of stripped SN explosions (without Hydrogen envelope)
            NSu: nuber of ultra-stripped SN explosions (bare CO-core stars)
        """

		if "EventsPlus" not in df.columns:
			raise ValueError("column EventsPlus not in dataframe")

		if inplace: df=df.copy(deep=True)


		df["NSt"] = df["EventsPlus"].str.count("St") # number of stripped SN explosions (without Hydrogen envelope)
		df["NSu"] = df["EventsPlus"].str.count("Su") # nuber of ultra-stripped SN explosions (bare CO-core stars)

		dd = Channel.prepare_history(df) # Series containing a simplified version of EventsPlus

		df["channel"] = dd.apply(Channel.chf)
		df["subchannel"] = dd.apply(Channel.subchf)
		df["collisions"] = dd.str.count("K")
		df["CE"] = dd.str.count("C")
		df["CEb"] = dd.apply(Channel.ceb)
		df["CEa"] = dd.apply(Channel.cea)
		df["SMTa"] = dd.apply(Channel.smta)
		df["EventsSimple"] = dd

		return df

	@staticmethod
	def chf(x: str) -> int:
		"""
        Find channel. The channel are based mainly on the history before the first SN explosion
        :param x: short history with events separated by :
        :return: formation channel:
            0- No interaction at all
            1- Only SMT before the first SN explosion
            2- Only SMT
            3- At least 1 CE before the first SN explosion, at the first SN explosion one star is stripped by the envelope
                - 33 Exactly 1 CE happens before the first SN between a star with a developed core and a MS
            4- At least 1 CE before the first SN explosion, at the first SN explosion both stars are stripped by the envelope
                - 44 Exactly 1 CE happens before the first SN, both stars in the CE have a core
            5- No interactions before the first SN explotion
            -1- All the other classes (Btw, the class before should contains all the possible cases)
        """
		# Notice, the if order is important, do not change it!

		try:
			# Channel0
			if ("C" not in x and "K" not in x and "M" not in x):
				return 0

			# ChannelII
			if ("C" not in x and "K" not in x and "M" in x):
				return 2

			subs0, subs1 = x.split("S")[:2]

			# ChannelV
			if ("C" not in subs0 and "K" not in subs0 and "M" not in subs0):
				return 5

			# ChannelI
			if ("C" not in subs0 and "K" not in subs0 and "M" in subs0):
				return 1

			# ChannelIII+IV
			if ("C" in subs0 or "K" in subs0):
				if subs0.count("C") == 1 and "M" not in subs0 and "ee" in subs0:
					return 44 # To follow the definition of channel 4 in Broekgaarden+21
				elif "ee" in subs0:
					return 4
				elif subs0.count("C") == 1:
					return 33 # To follow the definition of channel 3 in Broekgaarden+21
				elif subs0.count("C") > 1:
					return 3
		except Exception as err:
			print(f"Error raised in method chf  with input {x}. Rerasing the erro",flush=True)
			raise err

		return -1

	@staticmethod
	def subchf(x: str) -> str:
		"""
        Find sub channels. Sub channels are bases on the history after the first SN explosion
        :param x: short history with events separated by :
        :return: formation sub-channel:
            n- No interaction after the first SN explosion
            s- Only stable mass transfers after the first SN explosion
            c- Exactly 1 CE after the first SN explosion
            d- More than 1 CE after the first SN explosion
        """

		try:
			SMTa=Channel.smta(x) # At least 1 mass transfer after the first SN?
			CEa=Channel.cea(x) # How many CE after the first SN

			if CEa==1:
				return "c"
			elif  CEa>1:
				return "d"
			elif SMTa==0:
				return "n"
			elif SMTa==1:
				return "s"
			else:
				raise ValueError(f"Impossible to find subchannels for history {x}")
		except Exception as err:
			print(f"Error raised in method chf  with input {x}. Rerasing the erro",flush=True)
			raise err

		return ""

	@staticmethod
	def smta(x: str) -> int:
		"""
        Check if after the first SN there have been only stable mass transfer
        :param x: short history with events separated by :
        :return: 1 if there were only smt episodes after the first SN
        """

		try:
			if ("M" not in x):
				return 0

			subs = x.split("S")[1]

			if ("C"  in subs):
				return 0
			elif ("M" in subs):
				return 1

			return 0
		except Exception as err:
			print(f"Error raised in method chf  with input {x}. Rerasing the erro",flush=True)
			raise err

		return -1



	@staticmethod
	def ceb(x: str) -> int:
		"""
        Number of CE events before the explosion of the first SN
        :param x: short history with events separated by :
        :return: Number of CE events before the explosion of the first SN
        """

		try:
			if ("C" not in x):
				return 0

			subs = x.split("S")[0]

			return subs.count("C")
		except Exception as err:
			print(f"Error raised in method chf  with input {x}. Rerasing the erro",flush=True)
			raise err

		return -1

	@staticmethod
	def cea(x: str) -> int:
		"""
        Number of CE events after the explosion of the first SN
        :param x: short history with events separated by :
        :return: Number of CE events before the explosion of the first SN
        """

		try:
			if ("C" not in x):
				return 0

			subs = x.split("S")[1]

			return subs.count("C")
		except Exception as err:
			print(f"Error raised in method chf  with input {x}. Rerasing the erro",flush=True)
			raise err

		return -1

	@staticmethod
	def prepare_history(df: pd.DataFrame) -> pd.Series:
		"""
        Simplify history in EventsPlus so that it can be used for the other methods
        of this class.
        """
		if "EventsPlus" not in df.columns:
			raise ValueError("column EventsPlus not in dataframe")

		dd = df["EventsPlus"].str.replace("St|Su", "S", regex=True)
		dd = dd.str.replace("RC:", "", regex=True)
		dd = dd.str.replace("RB:RE", "M", regex=True)
		dd = dd.str.replace("RB:S:RE", "M:S", regex=True)
		dd = dd.str.replace("RB:C[s|d]?[e]?:RE", "C", regex=True)
		dd = dd.str.replace("RB:S", "M:S", regex=True)
		dd = dd.str.replace("C[s|d]?[e]?:RE", "C", regex=True)
		dd = dd.str.replace("K:C[s|d]?[e]?", "K:C", regex=True)
		dd = dd.str.replace("RB:", "M:", regex=True)
		dd = dd.str.replace("RE:", "", regex=True)

		return dd
class SEVN_output:

	def __init__(self,outfolder):

		self.outfolder = outfolder

	@staticmethod
	def get_loginfo(logfile):
		"""
		This method extract information from the log files
		Extract information about the Binary event for each system
			:param logfile: Path to the logfile to analyse
			:return: a pandas dataframe with columns: name, ID, Events, EventsPlus, EventsAll, SNtime, SNMej, SNvcom, SNcalpha, SNtype, SNvkick
			* The SN columns are replicated two times addind the suffix _0 or _1 (for the SN from the object with ID 0 or 1).
				The columns reports:
				- SNtime: explosion time of the SN in Myr
				- SNMej: total mass ejected in the SN (including neutrinos) in Myr
				- SNvcom: velocity of the binary centre of mass after the kick in (km/s) assuming that the previous Vcom=0 (it includes 
							also the effect of the Blauw kick)
				- SNcalpha: cosine of the angle between the normal to the old orbital plane (before the SN) and 
							the normal to the new orbital plane (after the SN)
				- SNtype: SN type (0-Unknown, 1- ECSN, 2- CCSN, 3- PPISN, 4-PISN, 5-Ia)
				- SNvkick: effective magnitude of the SN kick velocity in km/s
			* the Events column contains string with all the keyword of the events separated by a :,
				RB= RLO BEGIN
				RE= RLO END
				RC= CIRC
				C= Common Envelope
				K= Collision
				M= Merger
				W= Swallowed
				S= Supernova Explosion (notice that this will appear only if the system is still a binary)
				X= A Common envelope ending with a merger
			the EventsPlus column add some additional information to the tag
				- S, Supernova,
					- two prefixes are addes before S,
						-hh: before the explosion the binary was composed by two Hydrogen stars
						-he: before the explosion the binary was composed by an Hydrogen stars (the star that explodes) and a pureHe star
						-hr: before the explosion the binary was composed by an Hydrogen stars (the star that explodes) and a remnant
						-eh: before the explosion the binary was composed by a pureHe stars (the star that explodes) and an Hydrogen star
						-ee: before the explosion the binary was composed by two pureHe stars
						-er: before the explosion the binary was composed by a pureHe stars (the star that explodes)and a remnant
					- two suffixes can be added after S
						t= Stripped SN, the preSN star was almost or totally stripped by the Hydrogen envelope
							(Mtot-MHE)<0.02Mtot
						u= Ultra stripped SN, the preSN star was almost or totally stripped by the Hydrogen and and
							Helium envelope (Mtot-MHCO)<0.02Mtot
					e.g. ehSt -> before the explosion the sytems was composed by a pureHe star and an Hydrogen star,
								the pureHe star explosed as a stripped SN.
				- C and X, Common envelope:
					a first suffix is added after C or X:
						s= single core CE, a CE where just one of the two star has a clear core-envelope separation
						d= double core CE, a CE where both stars have a clear core-envelope separation
					an optional suffix can be present after C,X and s,d
						e= the star that started the CE is a pureHE star
					e.g. Cde
			the EventsAll columns contains the same labels of the EventPlus column and some additional labels
			for the RB events
				- RB, Roche Lobe Begin:
					- two number are added after RB:
						i0: the first number is 0 or 1 and indicates what star (primary=0, secondary=1) triggered
							the RLO
						i1: the second number (between 0 and 7) indicates the SEVN Phase of the star triggered the
							the RLO
					- An optional label is added
						e: if the star that triggered the RLO id a pureHe star


		"""
		dic_map={"RLO_BEGIN":"RB",
				 "RLO_BEGIN0":"RB0",
				 "RLO_BEGIN1":"RB1",
				 "RLO_BEGINhe0":"RBhe0",
				 "RLO_BEGINhe1":"RBhe1",
				 "RLO_BEGIN00":"RB00",
				 "RLO_BEGIN01":"RB01",
				 "RLO_BEGIN02":"RB02",
				 "RLO_BEGIN03":"RB03",
				 "RLO_BEGIN04":"RB04",
				 "RLO_BEGIN05":"RB05",
				 "RLO_BEGIN06":"RB06",
				 "RLO_BEGIN07":"RB07",
				 "RLO_BEGINhe00":"RB00e",
				 "RLO_BEGINhe01":"RB01e",
				 "RLO_BEGINhe02":"RB02e",
				 "RLO_BEGINhe03":"RB03e",
				 "RLO_BEGINhe04":"RB04e",
				 "RLO_BEGINhe05":"RB05e",
				 "RLO_BEGINhe06":"RB06e",
				 "RLO_BEGINhe07":"RB07e",
				 "RLO_BEGIN10":"RB10",
				 "RLO_BEGIN11":"RB11",
				 "RLO_BEGIN12":"RB12",
				 "RLO_BEGIN13":"RB13",
				 "RLO_BEGIN14":"RB14",
				 "RLO_BEGIN15":"RB15",
				 "RLO_BEGIN16":"RB16",
				 "RLO_BEGIN17":"RB17",
				 "RLO_BEGINhe10":"RB10e",
				 "RLO_BEGINhe11":"RB11e",
				 "RLO_BEGINhe12":"RB12e",
				 "RLO_BEGINhe13":"RB13e",
				 "RLO_BEGINhe14":"RB14e",
				 "RLO_BEGINhe15":"RB15e",
				 "RLO_BEGINhe16":"RB16e",
				 "RLO_BEGINhe17":"RB17e",
				 "RLO_END":"RE","CE":"C","CEs":"Cs","CEd":"Cd",
				 "CIRC": "RC",
				 "CEhe":"Ce","CEhes":"Cse","CEhed":"Cde",
				 "XE":"X","XEs":"Xs","XEd":"Xd",
				 "XEhe":"Xe","XEhes":"Xse","XEhed":"Xde",
				 "COLLISION":"K",
				 "MERGER":"M","SWALLOWED":"W",
				 "BSN":"S",
				 "BSNS":"St",
				 "BSNU":"Su",
				 "hhBSN":"hhS",
				 "heBSN":"heS",
				 "hrBSN":"hrS",
				 "ehBSNS":"ehSt",
				 "eeBSNS":"eeSt",
				 "erBSNS":"erSt",
				 "ehBSNU":"ehSu",
				 "eeBSNU":"eeSu",
				 "erBSNU":"erSu",}

		list_map = lambda x: ":".join([dic_map[_x] for _x in x])

		#Match patterns
		matchexp=r"[+|-]?[0-9]+\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan)"
		matchname=r"(?:[0-9|A-Za-z]*\_)?[0-9]*"
		matchid=r"[0-9]+"
		matchtype=r"[+|-]?\d+"


		with open(logfile,"r") as fo:
			i=fo.read()
			ma = re.findall(f"B;({matchname});({matchid});(.*);.*;",i)
		na = np.array(ma)
		nasimple = na.copy()


		#More information for SN explosion
		with open(logfile,"r") as fo:
			i=fo.read()
			sa = re.findall(f"B;({matchname});({matchid});BSN;({matchexp});({matchid}):({matchexp}):({matchexp}):({matchexp}):({matchtype}):"
							f"{matchtype}:\d:({matchexp}):({matchexp}):({matchexp}):({matchtype}):(?:{matchtype}):(?:{matchexp}):(?:{matchexp}):(?:{matchexp}):(:?{matchexp}):({matchexp}):({matchexp})",i)
		if (len(sa)>0):
			sa=np.array(sa)
			sa_float=np.array(sa[:,4:],dtype=float)
			sn_arr=np.empty(shape=(len(sa),3),dtype=object)
			sn_arr[:,0]=sa[:,0]
			sn_arr[:,1]=sa[:,1]
			sn_arr[:,2]="BSN"

			idx1h=(sa_float[:,0]>sa_float[:,1]) & (sa_float[:,3]<7)
			idx1he=(sa_float[:,0]==sa_float[:,1])
			idx2h=(sa_float[:,4]>sa_float[:,5]) & (sa_float[:,7]<7)
			idx2he=(sa_float[:,4]==sa_float[:,5])
			idx2r=sa_float[:,7]==7
			sn_arr[(sa_float[:,0]-sa_float[:,1])<0.02*sa_float[:,0],2]="BSNS"
			sn_arr[(sa_float[:,0]-sa_float[:,2])<0.02*sa_float[:,0],2]="BSNU"

			sn_arr[idx1h&idx2h,2]="hh"+sn_arr[idx1h&idx2h,2]
			sn_arr[idx1he&idx2h,2]="eh"+sn_arr[idx1he&idx2h,2]
			sn_arr[idx1h&idx2he,2]="he"+sn_arr[idx1h&idx2he,2]
			sn_arr[idx1he&idx2he,2]="ee"+sn_arr[idx1he&idx2he,2]
			sn_arr[idx1h&idx2r,2]="hr"+sn_arr[idx1h&idx2r,2]
			sn_arr[idx1he&idx2r,2]="er"+sn_arr[idx1he&idx2r,2]

		na[na[:,2]=="BSN",2]=sn_arr[:,2]


		# Info about SN (Sntime, SNvcom, SNcalpha)
		if (len(sa)>0):
			SNsid = np.array(sa[:,3],dtype=int)
			idx_SNid0 = SNsid==0
			idx_SNid1 = SNsid==1
			SNname_sid0 	 =  np.array(sa[idx_SNid0,0],dtype=object)
			SNname_sid1 	 =  np.array(sa[idx_SNid1,0],dtype=object)
			SNID_sid0 	 	 =  np.array(sa[idx_SNid0,1],dtype=object)
			SNID_sid1 		 =  np.array(sa[idx_SNid1,1],dtype=object)
			# SN explosion time
			SNtime_sid0 	 =  np.array(sa[idx_SNid0,2],dtype=float)
			SNtime_sid1 	 =  np.array(sa[idx_SNid1,2],dtype=float)
			# New centre of mass velocity after the kick (considering the initial Vcom=0) in km/s
			SNvcom_sid0 	 =  np.array(sa[idx_SNid0,-1],dtype=float)
			SNvcom_sid1 	 =  np.array(sa[idx_SNid1,-1],dtype=float)
			# cosine of the angle between the normal to the old orbital plane (before the event) and the normal to the new orbital plane (after the event)
			SNcalpha_sid0    =  np.array(sa[idx_SNid0,-2],dtype=float)
			SNcalpha_sid1	 =  np.array(sa[idx_SNid1,-2],dtype=float)
			dfSN_sid0 = pd.DataFrame({"name":SNname_sid0, "ID":SNID_sid0, "SNtime_0":SNtime_sid0, "SNvcom_0":SNvcom_sid0, "SNcalpha_0":SNcalpha_sid0})
			dfSN_sid1 = pd.DataFrame({"name":SNname_sid1, "ID":SNID_sid1, "SNtime_1":SNtime_sid1, "SNvcom_1":SNvcom_sid1, "SNcalpha_1":SNcalpha_sid1})
			dfSN = dfSN_sid0.merge(right=dfSN_sid1, on=("name","ID"),how="outer")
			# Why how=outer? Because in some rare case when two stars have exactly the same mass, a single BSN event is triggered
			# even if a compact binary object is created. So here we just merge anything and take into account that one of the two set of values can be set to nan
			del dfSN_sid0
			del dfSN_sid1
		else:
			dfSN = None

		## Complete the SN information with the event SN
		with open(logfile,"r") as fo:
			i=fo.read()
			ssesna = re.findall(f"S;({matchname});({matchid});SN;(?:{matchexp});({matchexp}):(?:{matchexp}):(?:{matchexp}):({matchexp}):(?:{matchtype}):({matchtype}):(?:{matchexp}):({matchexp})",i)

		if (len(ssesna)>0):
			ssesna=np.array(ssesna)
			# Additional info about SN (SNtype, SNvkick)
			SSESNsid = np.array(ssesna[:,1],dtype=int)
			Mej = np.array(ssesna[:,2],dtype=float) - np.array(ssesna[:,3],dtype=float)
			idx_SSESNid0 = SSESNsid==0
			idx_SSESNid1 = SSESNsid==1
			SSESNname_sid0 	     =  np.array(ssesna[idx_SSESNid0,0],dtype=object)
			SSESNname_sid1 	 	 =  np.array(ssesna[idx_SSESNid1,0],dtype=object)
			SSESNMej_sid0 	 	 =  np.array(Mej[idx_SSESNid0],dtype=float)
			SSESNMej_sid1 	 	 =  np.array(Mej[idx_SSESNid1],dtype=float)
			SSESNtype_sid0 	 	 =  np.array(ssesna[idx_SSESNid0,4],dtype=int)
			SSESNtype_sid1 	 	 =  np.array(ssesna[idx_SSESNid1,4],dtype=int)
			SSESNvkick_sid0 	 =  np.array(ssesna[idx_SSESNid0,5],dtype=float)
			SSESNvkick_sid1 	 =  np.array(ssesna[idx_SSESNid1,5],dtype=float)
			dfSSESN_sid0 = pd.DataFrame({"name":SSESNname_sid0, "SNMej_0":SSESNMej_sid0, "SNtype_0":SSESNtype_sid0, "SNvkick_0":SSESNvkick_sid0})
			dfSSESN_sid1 = pd.DataFrame({"name":SSESNname_sid1, "SNMej_1":SSESNMej_sid1, "SNtype_1":SSESNtype_sid1, "SNvkick_1":SSESNvkick_sid1})
			dfSSESN = dfSSESN_sid0.merge(right=dfSSESN_sid1, on=("name"),how="inner")
			# Here instead we use how=inner because to create a binary compact object we have triggered for sure two SN events
			del dfSSESN_sid0
			del dfSSESN_sid1
		else:
			dfSSESN = None

		if (dfSN is  None or dfSSESN is None):
			pass
		else:
			# Combine all the SN info
			dfSN = dfSSESN.merge(right=dfSN, on=("name"),how="inner")[["name","ID","SNtime_0","SNMej_0","SNvcom_0","SNcalpha_0","SNtype_0","SNvkick_0",
																 "SNtime_1","SNMej_1","SNvcom_1","SNcalpha_1","SNtype_1","SNvkick_1"]]
			# How=inner so that, we remove the single BSN event not combined with two single SN events


		#More on CE
		with open(logfile,"r") as fo:
			i=fo.read()
			ca = re.findall(f"B;({matchname});({matchid});CE;.*;\d:({matchexp}):({matchexp}):({matchexp}):{matchtype}:"
							f"{matchtype}:\d:({matchexp}):({matchexp}):({matchexp}):.*:(\d)",i)
		if (len(ca)>0):
			ca= np.array(ca,dtype=object)
			ca_float=np.array(ca[:,2:],dtype=float)
			cn_arr=np.empty(shape=(len(ca),3),dtype=object)
			cn_arr[:,0]=ca[:,0]
			cn_arr[:,1]=ca[:,1]
			idx_destroyed=ca_float[:,-1]>0
			cn_arr[:,2]="CE"
			cn_arr[idx_destroyed,2]="XE"
			idx_single=  ( (ca_float[:,3]==ca_float[:,4]) &  (ca_float[:,5]==0)) | (ca_float[:,4]==0)
			idx_he = (ca_float[:,0]==ca_float[:,1])
			cn_arr[idx_he,2]=cn_arr[idx_he,2]+"he"
			cn_arr[idx_single,2]=cn_arr[idx_single,2]+"s"
			cn_arr[~idx_single,2]=cn_arr[~idx_single,2]+"d"
			na[na[:,2]=="CE",2]=cn_arr[:,2]
		namedium=na.copy()

		#More on RLO
		with open(logfile,"r") as fo:
			i=fo.read()
			ra = re.findall(f"B;({matchname});({matchid});RLO_BEGIN;.*;(\d):({matchexp}):({matchexp}):(?:{matchexp}):({matchtype}):.*",i)
		if (len(ra)>0):
			ra= np.array(ra,dtype=object)
			ra_float=np.array(ra[:,2:],dtype=float)
			ra_arr=np.empty(shape=(len(ra),3),dtype=object)
			ra_arr[:,0]=ra[:,0]
			ra_arr[:,1]=ra[:,1]
			ra_arr[:,2]="RLO_BEGIN"
			idx_pureHe = ra_float[:,1]==ra_float[:,2]
			ra_arr[idx_pureHe,2]=ra_arr[idx_pureHe,2]+"he"
			ra_arr[:,2]=ra_arr[:,2]+ra[:,2]
			ra_arr[:,2]=ra_arr[:,2]+ra[:,-1]
			na[na[:,2]=="RLO_BEGIN",2]=ra_arr[:,2]

		df = pd.DataFrame(na,columns=["name","ID","EventsAll"])
		df["Events"] = nasimple[:,2]
		df["EventsPlus"] = namedium[:,2]
		df = df.groupby(["ID","name"]).agg(list_map).reset_index()

		# Add SN info
		if (dfSN is not None):
			df = df.merge(right=dfSN, on=("name","ID"), how="inner")

		return df.astype({"ID":"int64"})

	@staticmethod
	def get_COBs(outfile):
		"""
		Filter all the bound Compact object binaries (BHBH BHNS NSNS) from an outputfile.
		Notice this function assume that the simulation is stopped when both stars are remnant, therefore it gives more than
		one result for systems if the evolution is not stopped at the remnant formation.
		:param outfile: Path to the output file to analyse
		:return: a pandas dataframe containing the all COB binaries
		"""
		df = pd.read_csv(outfile)
		idx_BBH = (df.RemnantType_0 == 6) & (df.RemnantType_1 == 6)
		idx_BNS = ((df.RemnantType_0 == 4) | (df.RemnantType_0 == 5)) & ((df.RemnantType_1 == 4) | (df.RemnantType_1 == 5))
		idx_BNB = ((df.RemnantType_0 == 6) & ((df.RemnantType_1 == 4) | (df.RemnantType_1 == 5))) | (
				(df.RemnantType_1 == 6) & ((df.RemnantType_0 == 4) | (df.RemnantType_0 == 5)))
		idx = (idx_BBH | idx_BNS | idx_BNB) & (df.Semimajor > 0)
		dff = df[idx]

		return dff.astype({"name":object})

	@staticmethod
	def get_IC(outfile):
		"""

		:param outfile: Path to the evolved file to analyse
		:return: a pandas dataframe with the column:
			- ID: system id
			- name: sytem name
			- Mzams_0: initial zams mass of the first star
			- Mzams_1: initial zams mass of the second star
			- a: initial semimajor axis
			- e: initial eccentricity
			- Z: metallicity
		"""
		df = pd.read_csv(outfile, sep="\s+")
		df = df.rename(columns={"#ID": "ID",
								"Mass_0": "Mzams_0",
								"Mass_1": "Mzams_1",
								"a": "Semimajor_ini",
								"e": "Eccentricity_ini",
								"Z_0": "Z"})

		return df

	@staticmethod
	def make_COB_analysis_single(file_id, outpath):
		"""
		Return a pandas dataframe containg all the bound compact object binaries with information from the output, evolved and logfile
		:param file_id: number of thread that generate the output to analyse
		:param outpath: Path to the folder containing the sevn output
		:return:
		"""

		try:
			logfile = f"{outpath}/logfile_{file_id}.dat"
			outfile = f"{outpath}/output_{file_id}.csv"
			evolvedfile = f"{outpath}/evolved_{file_id}.dat"
			# History
			dfh = SEVN_output.get_loginfo(logfile)
			# Evolution
			dfe = SEVN_output.get_COBs(outfile)

			# Filter columns
			cols_in_file=dfe.columns
			wanted_columns=["ID", "name", "BWorldtime", "Mass_0", "Radius_0", "Zams_0", "Phase_0", "RemnantType_0", "Xspin_0","Mass_1", "Radius_1",
							"Phase_1", "Zams_1", "RemnantType_1", "Xspin_1","Semimajor", "Eccentricity", "GWtime"]
			effective_columns = [col for col in wanted_columns if col in cols_in_file]
			dfe=dfe[effective_columns]
			dfe["name"] = dfe["name"].astype('str')
			# IC
			dfi = SEVN_output.get_IC(evolvedfile)[["ID", "name", "Mzams_0", "Mzams_1", "Semimajor_ini", "Eccentricity_ini", "Z"]]
			dfi["name"] = dfi["name"].astype('str')

			dff = dfe.merge(dfh, on=["ID", "name"], how="left")
			dff = dff.merge(dfi, on=["ID", "name"], how="left")
		except Exception as err:
			print(f"Error when analysing sevn output ID {file_id} in folder {outpath}. See error details below",flush=True)
			raise err 
		
		return dff

	def make_COB_analysis(self,nproc=4):
		"""

		:param nproc: number of parallel process to use
		:return:
		"""
		files = glob.glob(self.outfolder + "/logfile_*")
		idlist = [int(x.split("logfile_")[1].split(".")[0]) for x in files]

		if nproc==1:
			dfl=[]
			for fileid in idlist: 
				dfl.append(SEVN_output.make_COB_analysis_single(fileid,self.outfolder))
		else:
			with Pool(nproc) as p:
				dfl = p.starmap(SEVN_output.make_COB_analysis_single, zip(idlist, repeat(self.outfolder)))

		return pd.concat(dfl)

	@staticmethod
	def split_COBs(df):

		idxbbh = (df.RemnantType_0==6) & (df.RemnantType_1==6) & (df.Semimajor.notnull())
		idxbns = ( (df.RemnantType_0==4) | (df.RemnantType_0==5) ) &  ((df.RemnantType_1==4) | (df.RemnantType_1==5) ) & (df.Semimajor.notnull())
		idxbhns = ( ( ((df.RemnantType_0==4) | (df.RemnantType_0==5)) & (df.RemnantType_1==6)  ) \
					| ( ((df.RemnantType_1==4) | (df.RemnantType_1==5)) & (df.RemnantType_0==6)  )  ) & (df.Semimajor.notnull())

		dfbbh  = df[idxbbh]
		dfbns  = df[idxbns]
		dfbhns = df[idxbhns]

		return dfbbh, dfbns, dfbhns

	@staticmethod
	def merging(df,tshold=14000):

		if tgw_loaded:
			tmerge = estimate_tgw(df.Semimajor,df.Eccentricity,df.Mass_0,df.Mass_1,method="combined",nproc=1,a_Rsun=True)
		else:
			tmerge = df.GWtime

		idxmerging = df.BWorldtime + tmerge < tshold

		return df[idxmerging]



def analyse_sevn_output(input_folder, nproc=2, output_older=None):
	"""
	Analyse the file of a sevn output filtering the compact objectes binaries mixing information from output, evolved an logfiles.
	6 csv files will be saved:
	- BHBH.csv: containing all the bound BHBH binaries
	- BHBHm.csv: containing all the bound BHBH binaries that will merge within an Hubble time (14000 My)
	- NSNS.csv: containing all the bound NSNS binaries
	- NSNSm.csv: containing all the bound NSNS binaries that will merge within an Hubble time (14000 My)
	- BHNS.csv: containing all the bound BHNS binaries
	- BHNSm.csv: containing all the bound BHNSm binaries that will merge within an Hubble time (14000 My)
	:param input_folder: Name of sevn output folder containt the output, evolved and logfiles
	:param nproc: number of parallel processes to use
	:param output_older: name of the folder where to store the output file, if None the same input folder will be used
	:return: 0
	"""
	if output_older is None:
		output_older = input_folder

	if not os.path.isdir(output_older):
		os.makedirs(output_older)

	so = SEVN_output(input_folder)
	df = so.make_COB_analysis(nproc)
	df = Channel.channellify(df) #Add channel information
	dfbbh, dfbns, dfbhns = so.split_COBs(df)

	dfbbh.to_csv(output_older + "/BHBH.csv", index=False)
	dfbhbhm = so.merging(dfbbh)
	dfbhbhm.to_csv(output_older+ "/BHBHm.csv", index=False)

	dfbns.to_csv(output_older + "/NSNS.csv", index=False)
	dfbnsm = so.merging(dfbns)
	dfbnsm.to_csv(output_older + "/NSNSm.csv", index=False)

	dfbhns.to_csv(output_older + "/BHNS.csv", index=False)
	dfbhnsm = so.merging(dfbhns)
	dfbhnsm.to_csv(output_older + "/BHNSm.csv", index=False)

	return 0

def analyse_multi_sevn_output(main_folder,  nproc=2, sevn_root="sevn_output", subfolder_name=None):
	"""
	Analyse the file of multiple sevn outputs folder having the same root name.
	For each simulation folder the output files are analysed filtering the compact objectes binaries mixing information from output, evolved an logfiles.
	6 csv files will be saved:
	- BHBH.csv: containing all the bound BHBH binaries
	- BHBHm.csv: containing all the bound BHBH binaries that will merge within an Hubble time (14000 My)
	- NSNS.csv: containing all the bound NSNS binaries
	- NSNSm.csv: containing all the bound NSNS binaries that will merge within an Hubble time (14000 My)
	- BHNS.csv: containing all the bound BHNS binaries
	- BHNSm.csv: containing all the bound BHNSm binaries that will merge within an Hubble time (14000 My)
	:param main_folder: path to the folder containing the multiple sevn output folders
	:param nproc: number of parallel processes to use
	:param sevn_root: common root name of the sevn output folders, e.g. if we have sevn_out1 sevn_out2,
	the root name should be sevn_out
	:param subfolder_name: Sometime the outputs are not contained directly in the main sevn output folder but thera a number
	of subfolders to consider, this is the name of the subfolder we want to analyse.
	:return: 0
	"""

	folders = glob.glob(main_folder+"/"+sevn_root+"*")

	for folder in folders:
		print("Analysing folder ",folder, flush=True)
		fname=folder
		if subfolder_name is not None:  fname = folder+"/"+subfolder_name
		try:
			analyse_sevn_output(fname, nproc=nproc, output_older=folder)
			print("Done with folder ",folder, flush=True)
		except Exception as e:
			print("Error analyse folder ",folder," with message: ",file=sys.stderr)
			print(e,file=sys.stderr)


	return 0

if __name__=="__main__":

	info = "	Analyse the file of a sevn output filtering the compact objectes binaries mixing information from output, evolved an logfiles.\n"
	info += "6 csv files will be saved:\n"
	info += "- BHBH.csv: containing all the bound BHBH binaries\n"
	info += "- BHBHm.csv: containing all the bound BHBH binaries that will merge within an Hubble time (14000 My)\n"
	info += "- NSNS.csv: containing all the bound NSNS binaries\n"
	info += "- NSNSm.csv: containing all the bound NSNS binaries that will merge within an Hubble time (14000 My)\n"
	info += "- BHNS.csv: containing all the bound BHNS binaries\n"
	info += "- BHNSm.csv: containing all the bound BHNS binaries that will merge within an Hubble time (14000 My)"


	parser = argparse.ArgumentParser(description=info)
	parser.add_argument('input', type=str, help='path to the sevn folder with the input to analyse, if option -m is enabled, this'
												'is the root folder containing all the sevn output folders to analyse')
	parser.add_argument('-n','--nproc', dest='nproc', default=2, type=int, help="Number of parallel processes to use [2]")
	parser.add_argument('-o', '--output', dest='out', default=None, type=str, help="Output folder where to store the files, if None the input folder is used [None]")
	parser.add_argument('-m', '--multifolders', dest='multifolders', action='store_true', default=False, help="If true analyse a list of output folders inside the path specified")
	parser.add_argument('--froot', dest='froot', default="sevn_output", type=str, help="common root for multiple sevn output analysis, used only if option -m is enabled [sevn_output]")
	parser.add_argument('--subfolder', dest='subfolder', default=None, type=str, help="Subfolder containg the SEVN output, used only if option -m is enabled [None]"
																					  "on input and with folder root specified with the froot parameter [False]")
	args = parser.parse_args()

	if (args.multifolders):
		analyse_multi_sevn_output(args.input,args.nproc,args.froot,args.subfolder)
	else:
		analyse_sevn_output(args.input,args.nproc,args.out)
