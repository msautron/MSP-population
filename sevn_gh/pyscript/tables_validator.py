#!/usr/bin/env python
"""
SCRIPT: tables_validator
PURPOSE: Make a series of check on a given tables set. For each check it returns:
		- OK if there are no problems
		- WARNING, if a given feature is not  standard for SEVN, but SEVN can handle it.
		- ERROR,  a given feature is not standard for SEVN, this make the table not usable in SEVN

USAGE; ./tables_valitor <TABLES_PATH>  --plot[OPTIONAL] --purehe[Optional]
REQUIRED PACKAGE:
    -numpy
    -matplotlib

V 1.0: 06/04/21 Giuliano Iorio
"""

import numpy as np
import sys
import glob
import matplotlib.pyplot as plt
import matplotlib
import os

CHECK_FALSE    = 0
CHECK_WARNING  = 1
CHECK_TRUE     = 2

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

class FontColors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKCYAN = '\033[96m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'

"""
Decorator to handle the messages in thest.
"""
def check_message(message=""):
	def decorator(f):
		def wrap(*args, **kwargs):
			print(FontColors.HEADER+"*"+message+FontColors.ENDC+"...",end="",flush=True)
			ret = f(*args, **kwargs)
			if ret.check==CHECK_TRUE:
				print(FontColors.OKGREEN+FontColors.BOLD+"OK"+FontColors.ENDC,flush=True)
			elif ret.check==CHECK_WARNING:
				print(FontColors.WARNING+FontColors.BOLD+"WARNING"+FontColors.ENDC,flush=True)
				print("WRN message: %s"%ret.msg)
			elif ret.check==CHECK_FALSE:
				print(FontColors.FAIL+FontColors.BOLD+"ERROR"+FontColors.ENDC,flush=True)
				print("ERR message: %s"%ret.msg)
			return ret
		return wrap
	return decorator

class SEVNTablesValidator:

	def __init__(self,folder_path):

		self._fundamental_tables = ("mass.dat", "mhe.dat", "mco.dat", "radius.dat", "lumi.dat", "time.dat", "phase.dat")
		self._semifundamental_tables = ("inertia.dat","rco.dat","rhe.dat")
		self._sup_tables = ("csup.dat","hesup.dat","hsup.dat","nsup.dat","osup.dat")
		self._envconv_tables = ("depthconv.dat","qconv.dat","tconv.dat")
		self._all_tables = self._fundamental_tables + self._semifundamental_tables + self._sup_tables + self._envconv_tables

		self.folders = np.sort(glob.glob(folder_path + "/*"))
		self.summary_dict={ folder : {} for folder in self.folders }

		self.len_phases=14

		##RETURN TYPE
		self.CHECK_FALSE    = CHECK_FALSE
		self.CHECK_WARNING  = CHECK_WARNING
		self.CHECK_TRUE     = CHECK_TRUE

	class CheckReturn:
		__slots__ = ['check', 'msg']
		def __init__(self,check,msg=""):

			if (check!=CHECK_FALSE and check!=CHECK_TRUE and check!=CHECK_WARNING):
				raise ValueError("Unkown check value")

			self.check = check
			self.msg = msg

	def check_all(self):

		for folder in self.folders:
			Z=folder.split("/")[-1]
			Z=float(Z[0]+"."+Z[1:])
			self.check_all_folder(folder)
			self.summary_dict[folder]["Z"]=Z

		print("***** FINAL SUMMARY ******")
		for element in self.summary_dict.values():
			print(FontColors.HEADER+"Z=%s"%(element["Z"]) + FontColors.ENDC +  "-> %s (NWarnings=%i, NErrors=%i)"%(element["result"], element["Nwarning"],element["Nerror"]))
		print("***************************")

	def check_all_folder(self,folder):

		all_ret=[]
		print("------------------------------------------------", flush=True)
		print("Checking folder \"%s\"" % folder, flush=True)
		print("------------------------------------------------", flush=True)
		all_ret.append(self.check_metallicity_folder_name(folder))
		all_ret.append(self.check_fundamental_tables(folder))
		all_ret.append(self.check_tables_name(folder))
		all_ret.append(self.check_table_dimensions(folder))
		all_ret.append(self.check_mass(folder))
		all_ret.append(self.check_radius(folder))
		all_ret.append(self.check_times(folder))
		all_ret.append(self.check_he_core_time(folder))
		all_ret.append(self.check_he_core_mass(folder))
		all_ret.append(self.check_rhe_time(folder))
		all_ret.append(self.check_rhe(folder))
		all_ret.append(self.check_co_core_time(folder))
		all_ret.append(self.check_co_core_mass(folder))
		all_ret.append(self.check_rco_time(folder))
		all_ret.append(self.check_rco(folder))
		print("***** Check Summary *****")
		final_ret = np.array([x.check for x in all_ret])
		NWarning = np.sum(final_ret==CHECK_WARNING)
		NError = np.sum(final_ret==CHECK_FALSE)
		if (NError>0):
			result = FontColors.FAIL+"Not usable in SEVN"+FontColors.ENDC
		elif (NWarning>0):
			result = FontColors.WARNING+"Usable in SEVN, but be careful of Warnings"+FontColors.ENDC
		else:
			result = FontColors.OKGREEN+"Usable in SEVN"+FontColors.ENDC
		print(FontColors.WARNING+"NWarning"+FontColors.ENDC,NWarning)
		print(FontColors.FAIL+"NError"+FontColors.ENDC,NError)
		print(FontColors.HEADER+"Results"+FontColors.ENDC,result)
		self.summary_dict[folder]["Nwarning"]=NWarning
		self.summary_dict[folder]["Nerror"]=NError
		self.summary_dict[folder]["result"]=result
		print("***** INFO *****")
		self.info_zams(folder)


	@check_message("Checking folder name")
	def check_metallicity_folder_name(self,folder_path):
		name = folder_path.split("/")[-1]

		err_msg = ""
		if not name.isdigit() or name[0] != "0":
			msg_metallicity = "Metallicity folder name must contains only digits and start with 0"
			err_msg = "Invalid folder name %s. %s" % (name, msg_metallicity)
			return self.CheckReturn(self.CHECK_FALSE,err_msg)

		return self.CheckReturn(self.CHECK_TRUE)

	@check_message("Checking fundamental tables")
	def check_fundamental_tables(self,folder_path):

		missing_tables = []
		for table_name in self._fundamental_tables:
			if not os.path.isfile(folder_path + "/" + table_name):
				missing_tables.append(table_name)

		if len(missing_tables) > 0:
			err_msg = "Fundamental tables are missing: %s" % (", ".join(missing_tables))
			return self.CheckReturn(self.CHECK_FALSE,err_msg)

		return self.CheckReturn(self.CHECK_TRUE)

	@check_message("Checking unexpected  tables/files")
	def check_tables_name(self, folder_path):

		cumulated_warnings=[]
		#Read tables
		tables = self.read_tables_onlyname(folder_path)
		for table in tables:
			if  not table   in self._all_tables:
				cumulated_warnings.append(" - File: %s"%table)

		if len(cumulated_warnings)>0:
			wrn_msg = "Unknown file(s) found:\n%s" %("\n".join(cumulated_warnings))
			return self.CheckReturn(self.CHECK_WARNING,wrn_msg)

		return self.CheckReturn(self.CHECK_TRUE)

	@check_message("Checking tables dimension")
	def check_table_dimensions(self,folder_path):

		tables = self.read_tables(folder_path)
		tables_name=[]
		nrows_list=[]
		ncols_list=[]
		for table in tables:
			table_name=table.split("/")[-1]
			tables_name.append(table_name)
			ncols_list_tmp=[]
			with open(table,"r") as fo:
				ncols = [len(line.split()) for line in fo.readlines()]
				nrows = len(ncols)
				nrows_list.append(nrows)
				ncols_list_tmp.append(ncols)
				self.summary_dict[folder_path][table_name]={"nrows": nrows, "ncols":ncols}
			ncols_list.append(ncols_list_tmp)

		tables_name=np.array(tables_name)
		nrows_list=np.array(nrows_list)
		ncols_list=np.array(ncols_list,dtype=object)

		#Check rows
		if len(np.unique(nrows_list))>1:
			idx=np.argsort(nrows_list)
			n_row_str=[" - %s nrows=%i"%("{:<15}".format(tname),nrows) for tname,nrows in zip(tables_name[idx],nrows_list[idx])]
			n_row_str="\n".join(n_row_str)
			err_msg="Number of row mismatch in tables: \n"+n_row_str
			return self.CheckReturn(self.CHECK_FALSE,err_msg)

		#Rows are ok check columns
		cumulative_err_msg=[]
		for row in range(nrows_list[0]): #Here we are sure that all the files have the same number of rows
			ncols_list_at_row = []

			for tname in tables_name:
				if tname!="phase.dat":
					ncols_list_at_row.append(self.summary_dict[folder_path][tname]["ncols"][row])
				else:
					ncols_phase_list_at_row=self.summary_dict[folder_path][tname]["ncols"][row]

			ncols_list_at_row = np.array(ncols_list_at_row)
			ncols_phase_list_at_row = np.array(ncols_phase_list_at_row)


			#Check
			if ncols_phase_list_at_row!=self.len_phases:
				err_msg = "\n   -Row: %i -> Number of columns in table phases is %i instead of  14"%(row,ncols_phase_list_at_row)
				cumulative_err_msg.append(err_msg)


			if (len(np.unique(ncols_list_at_row))>1):
				err_msg = "\n  -Row: %i -> Number of columns mismatch in tables \n" % row
				tables_name_nophase=tables_name[tables_name!="phase.dat"]
				idx=np.argsort(ncols_list_at_row)
				for tname, ncols in zip(tables_name_nophase[idx],ncols_list_at_row[idx]):
					err_msg += "    - %s ncols=%i \n"%("{:<15}".format(tname),ncols)
				cumulative_err_msg.append(err_msg)


		if (len(cumulative_err_msg)>0):

			return self.CheckReturn(self.CHECK_FALSE,"\n".join(cumulative_err_msg))



		return self.CheckReturn(self.CHECK_TRUE)

	@check_message("Checking He core mass time")
	def check_he_core_time(self,folder_path):

		fmass=open(folder_path+"/mass.dat","r")
		fhe=open(folder_path+"/mhe.dat","r")
		ftime=open(folder_path+"/time.dat","r")
		fphase=open(folder_path+"/phase.dat","r")

		cumulative_err_msg=[]

		row=1
		for mass, hemasses, times, phases in zip(self.read_first(fmass), self.read_lines(fhe), self.read_lines(ftime), self.read_phases(fphase)):

			col_list=np.arange(len(times))

			times_with_he = times[hemasses != 0]
			times_without_he = times[hemasses == 0]
			t_he_start = phases[2]

			he_wrongly_nonzero = times_with_he <= t_he_start
			t_he_wrongly_nonzero = times_with_he[he_wrongly_nonzero]
			he_wrongly_zero = times_without_he > t_he_start
			t_he_wrongly_zero = times_without_he[he_wrongly_zero]


			if (len(t_he_wrongly_nonzero)>0 or len(t_he_wrongly_zero)>0):

				msg="\nProblem for \033[96m Zams  %s \033[0m (row %i, tphase_2=%f): \n"%(mass,row,t_he_start)

				for  i,t_error in enumerate(t_he_wrongly_nonzero):
					msg+=" - MHE>0 before tphase_2 at t=%f phase=%i\n"%(t_error,self.find_phase(t_error,phases))

				for  i,t_error in enumerate(t_he_wrongly_zero):
					msg+=" - MHE=0 after tphase_2 at t=%f phase=%i\n"%(t_error,self.find_phase(t_error,phases))

				cumulative_err_msg.append(msg[:-1]) #-1 to remove the last/n


			row+=1

		fmass.close()
		fhe.close()
		ftime.close()
		fphase.close()


		if len(cumulative_err_msg)>0:
			return self.CheckReturn(self.CHECK_FALSE,"\n".join(cumulative_err_msg))

		return  self.CheckReturn(self.CHECK_TRUE)

	@check_message("Checking He core mass")
	def check_he_core_mass(self,folder_path):

		fmass=open(folder_path+"/mass.dat","r")
		fhe=open(folder_path+"/mhe.dat","r")
		ftime=open(folder_path+"/time.dat","r")
		fphase=open(folder_path+"/phase.dat","r")

		cumulative_err_msg=[]

		row=1
		for masses, hemasses, times, phases in zip(self.read_lines(fmass), self.read_lines(fhe), self.read_lines(ftime), self.read_phases(fphase)):

			problematic_idx = hemasses>=masses
			col_list=np.arange(len(masses))

			if (np.sum(problematic_idx)>0):
				msg="\nProblem for \033[96m Zams  %s \033[0m (row %i): \n"%(masses[0],row)

				for time,mass,mhe,col in zip(times[problematic_idx],masses[problematic_idx],hemasses[problematic_idx],col_list[problematic_idx]):
					msg+=" - MHE(%f)>=Mass(%f) at t=%f (col=%i/%i) phase=%i\n"%(mhe,mass,time,col,len(masses)-1,self.find_phase(time,phases))


				cumulative_err_msg.append(msg[:-1]) #-1 to remove the last/n


			row+=1

		fmass.close()
		fhe.close()
		ftime.close()
		fphase.close()


		if len(cumulative_err_msg)>0:
			return self.CheckReturn(self.CHECK_WARNING,"\n".join(cumulative_err_msg))

		return  self.CheckReturn(self.CHECK_TRUE)

	@check_message("Checking CO core mass time")
	def check_co_core_time(self,folder_path):

		fmass=open(folder_path+"/mass.dat","r")
		fco=open(folder_path+"/mco.dat","r")
		ftime=open(folder_path+"/time.dat","r")
		fphase=open(folder_path+"/phase.dat","r")

		cumulative_err_msg=[]

		row=1
		for mass, comasses, times, phases in zip(self.read_first(fmass), self.read_lines(fco), self.read_lines(ftime), self.read_phases(fphase)):


			times_with_co = times[comasses != 0]
			times_without_co = times[comasses == 0]
			t_co_start = phases[5]

			co_wrongly_nonzero = times_with_co <= t_co_start
			t_co_wrongly_nonzero = times_with_co[co_wrongly_nonzero]
			co_wrongly_zero = times_without_co > t_co_start
			t_co_wrongly_zero = times_without_co[co_wrongly_zero]


			if (len(t_co_wrongly_nonzero)>0 or len(t_co_wrongly_zero)>0):

				msg="\nProblem for \033[96m Zams  %s \033[0m (row %i, tphase_5=%f): \n"%(mass,row,t_co_start)

				for  i,t_error in enumerate(t_co_wrongly_nonzero):
					msg+=" - MCO>0 before tphase_5 at t=%f phase=%i\n"%(t_error,self.find_phase(t_error,phases))

				for  i,t_error in enumerate(t_co_wrongly_zero):
					msg+=" - MCO=0 after tphase_5 at t=%f phase=%i\n"%(t_error,self.find_phase(t_error,phases))

				cumulative_err_msg.append(msg[:-1]) #-1 to remove the last/n


			row+=1

		fmass.close()
		fco.close()
		ftime.close()
		fphase.close()


		if len(cumulative_err_msg)>0:
			return self.CheckReturn(self.CHECK_FALSE,"\n".join(cumulative_err_msg))

		return  self.CheckReturn(self.CHECK_TRUE)

	@check_message("Checking CO core mass")
	def check_co_core_mass(self,folder_path):

		fmass=open(folder_path+"/mass.dat","r")
		fco=open(folder_path+"/mco.dat","r")
		fhe=open(folder_path+"/mhe.dat","r")
		ftime=open(folder_path+"/time.dat","r")
		fphase=open(folder_path+"/phase.dat","r")

		cumulative_err_msg=[]

		row=1
		for masses, comasses,hemasses, times, phases in zip(self.read_lines(fmass),self.read_lines(fco), self.read_lines(fhe), self.read_lines(ftime), self.read_phases(fphase)):

			problematic_idx_decreasing = (np.diff(comasses,prepend=comasses[0])/comasses)<-0.01
			problematic_idx_tmass = (comasses>=masses) & (problematic_idx_decreasing==False)
			problematic_idx_hemass = (comasses>=hemasses) & (problematic_idx_tmass==False) & (problematic_idx_decreasing==False) & (comasses!=0)


			col_list=np.arange(len(masses))

			if (np.sum(problematic_idx_tmass)>0 or np.sum(problematic_idx_hemass)>0 or np.sum(problematic_idx_decreasing)) :
				msg="\nProblem for \033[96m Zams  %s \033[0m (row %i): \n"%(masses[0],row)

				for time,mco,col in zip(times[problematic_idx_decreasing],comasses[problematic_idx_decreasing],col_list[problematic_idx_decreasing]):
					msg+=" - Decreasing MCO (more than 1pcent): MCO[i+1](%.8f)<MCO[i](%.8f) at t=%f (col=%i) phase=%i\n"%(comasses[col],comasses[col-1],time,col,self.find_phase(time,phases))


				for time,mass,mco,col in zip(times[problematic_idx_tmass],masses[problematic_idx_tmass],comasses[problematic_idx_tmass],col_list[problematic_idx_tmass]):
					msg+=" - MCO(%f)>=Mass(%f) at t=%f (col=%i/%i) phase=%i\n"%(mco,mass,time,col,len(masses)-1,self.find_phase(time,phases))

				for time,mhe,mco,col in zip(times[problematic_idx_hemass],hemasses[problematic_idx_hemass],comasses[problematic_idx_hemass],col_list[problematic_idx_hemass]):
					msg+=" - MCO(%f)>=MHE(%f) at t=%f (col=%i/%i) phase=%i\n"%(mco,mhe,time,col,len(masses)-1,self.find_phase(time,phases))


				cumulative_err_msg.append(msg[:-1]) #-1 to remove the last/n


			row+=1

		fmass.close()
		fco.close()
		fhe.close()
		ftime.close()
		fphase.close()


		if len(cumulative_err_msg)>0:
			return self.CheckReturn(self.CHECK_WARNING,"\n".join(cumulative_err_msg))

		return  self.CheckReturn(self.CHECK_TRUE)

	@check_message("Checking radius")
	def check_radius(self,folder_path):

		fmass=open(folder_path+"/mass.dat","r")
		fradius=open(folder_path+"/radius.dat","r")
		ftime=open(folder_path+"/time.dat","r")


		cumulative_err_msg=[]

		row=1
		for zams, radii, times in zip(self.read_first(fmass),self.read_lines(fradius),self.read_lines(ftime)):

			problematic_idx = radii<=0.0

			if np.sum(problematic_idx)>0 :
				msg="\nProblem for \033[96m Zams  %s \033[0m (row %i): \n"%(zams,row)


				col_list=np.arange(len(radii))[problematic_idx]

				for rad,time,col in zip(radii[problematic_idx],times[problematic_idx],col_list):
					msg+=" - R(%f)<=0.0 at t=%f (col=%i/%u) \n"%(rad,time,col,len(radii)-1)

				cumulative_err_msg.append(msg[:-1]) #-1 to remove the last/n


			row+=1

		if len(cumulative_err_msg)>0:
			return self.CheckReturn(self.CHECK_ERROR,"\n".join(cumulative_err_msg))

		fmass.close()
		fradius.close()
		ftime.close()

		return  self.CheckReturn(self.CHECK_TRUE)

	@check_message("Checking HE core radius times")
	def check_rhe_time(self,folder_path):

		if not os.path.isfile(folder_path+'/rhe.dat'):

			return  self.CheckReturn(self.CHECK_TRUE)

		else:

			fmass=open(folder_path+"/mass.dat","r")
			fhe=open(folder_path+"/rhe.dat","r")
			ftime=open(folder_path+"/time.dat","r")
			fphase=open(folder_path+"/phase.dat","r")

			cumulative_err_msg=[]

			row=1
			for mass, hemasses, times, phases in zip(self.read_first(fmass), self.read_lines(fhe), self.read_lines(ftime), self.read_phases(fphase)):

				times_with_he = times[hemasses != 0]
				times_without_he = times[hemasses == 0]
				t_he_start = phases[2]

				he_wrongly_nonzero = times_with_he <= t_he_start
				t_he_wrongly_nonzero = times_with_he[he_wrongly_nonzero]
				he_wrongly_zero = times_without_he > t_he_start
				t_he_wrongly_zero = times_without_he[he_wrongly_zero]


				if (len(t_he_wrongly_nonzero)>0 or len(t_he_wrongly_zero)>0):

					msg="\nProblem for \033[96m Zams  %s \033[0m (row %i, tphase_2=%f): \n"%(mass,row,t_he_start)

					for  i,t_error in enumerate(t_he_wrongly_nonzero):
						msg+=" - RHE>0 before tphase_2 at t=%f phase=%i\n"%(t_error,self.find_phase(t_error,phases))

					for  i,t_error in enumerate(t_he_wrongly_zero):
						msg+=" - RHE=0 after tphase_2 at t=%f phase=%i\n"%(t_error,self.find_phase(t_error,phases))

					cumulative_err_msg.append(msg[:-1]) #-1 to remove the last/n


				row+=1

			fmass.close()
			fhe.close()
			ftime.close()
			fphase.close()

			if len(cumulative_err_msg)>0:
				return self.CheckReturn(self.CHECK_FALSE,"\n".join(cumulative_err_msg))


		return  self.CheckReturn(self.CHECK_TRUE)

	@check_message("Checking He core radius")
	def check_rhe(self, folder_path):

		if not os.path.isfile(folder_path+'/rhe.dat'):
			return  self.CheckReturn(self.CHECK_TRUE)

		fmass=open(folder_path+"/mass.dat","r")
		fmhe=open(folder_path+"/mhe.dat","r")
		fhe=open(folder_path+"/rhe.dat","r")
		frad=open(folder_path+"/radius.dat","r")
		ftime=open(folder_path+"/time.dat","r")
		fphase=open(folder_path+"/phase.dat","r")

		cumulative_err_msg=[]

		row=1
		for masses, hemasses, rads, rhes, times, phases in zip(self.read_lines(fmass), self.read_lines(fmhe), self.read_lines(frad), self.read_lines(fhe) , self.read_lines(ftime), self.read_phases(fphase)):

			problematic_idx_typea = (rhes>=rads) & (hemasses!=masses)
			problematic_idx_typeb = (rhes!=rads) & (hemasses==masses)
			col_list=np.arange(len(masses))

			if (np.sum(problematic_idx_typea)>0 or np.sum(problematic_idx_typeb)>0):
				msg="\nProblem for \033[96m Zams  %s \033[0m (row %i): \n"%(masses[0],row)

				for time,rad,rhe,mass,mhe,col in zip(times[problematic_idx_typea],rads[problematic_idx_typea],rhes[problematic_idx_typea],masses[problematic_idx_typea],hemasses[problematic_idx_typea],col_list[problematic_idx_typea]):
					msg+=" - RHE(%f)>=Radius(%f) when MHE(%f)!=Mass(%f) at t=%f (col=%i/%i) phase=%i\n"%(rhe,rad,mhe,mass,time,col,len(masses)-1,self.find_phase(time,phases))

				for time,rad,rhe,mass,mhe,col in zip(times[problematic_idx_typeb],rads[problematic_idx_typeb],rhes[problematic_idx_typeb],masses[problematic_idx_typeb],hemasses[problematic_idx_typeb],col_list[problematic_idx_typeb]):
					msg+=" - RHE(%f)!=Radius(%f) when MHE(%f)==Mass(%f) at t=%f (col=%i/%i) phase=%i\n"%(rhe,rad,mhe,mass,time,col,len(masses)-1,self.find_phase(time,phases))




				cumulative_err_msg.append(msg[:-1]) #-1 to remove the last/n


			row+=1

		frad.close()
		fmass.close()
		fmhe.close()
		fhe.close()
		ftime.close()
		fphase.close()


		if len(cumulative_err_msg)>0:
			return self.CheckReturn(self.CHECK_WARNING,"\n".join(cumulative_err_msg))

		return  self.CheckReturn(self.CHECK_TRUE)

	@check_message("Checking CO core radius times")
	def check_rco_time(self,folder_path):


		if not os.path.isfile(folder_path+'/rco.dat'):

			return  self.CheckReturn(self.CHECK_TRUE)

		fmass=open(folder_path+"/mass.dat","r")
		fco=open(folder_path+"/rco.dat","r")
		ftime=open(folder_path+"/time.dat","r")
		fphase=open(folder_path+"/phase.dat","r")

		cumulative_err_msg=[]

		row=1
		for mass, comasses, times, phases in zip(self.read_first(fmass), self.read_lines(fco), self.read_lines(ftime), self.read_phases(fphase)):


			times_with_co = times[comasses != 0]
			times_without_co = times[comasses == 0]
			t_co_start = phases[5]

			co_wrongly_nonzero = times_with_co <= t_co_start
			t_co_wrongly_nonzero = times_with_co[co_wrongly_nonzero]
			co_wrongly_zero = times_without_co > t_co_start
			t_co_wrongly_zero = times_without_co[co_wrongly_zero]


			if (len(t_co_wrongly_nonzero)>0 or len(t_co_wrongly_zero)>0):

				msg="\nProblem for \033[96m Zams  %s \033[0m (row %i, tphase_5=%f): \n"%(mass,row,t_co_start)

				for  i,t_error in enumerate(t_co_wrongly_nonzero):
					msg+=" - RCO>0 before tphase_5 at t=%f phase=%i\n"%(t_error,self.find_phase(t_error,phases))

				for  i,t_error in enumerate(t_co_wrongly_zero):
					msg+=" - RCO=0 after tphase_5 at t=%f phase=%i\n"%(t_error,self.find_phase(t_error,phases))

				cumulative_err_msg.append(msg[:-1]) #-1 to remove the last/n


			row+=1

		fmass.close()
		fco.close()
		ftime.close()
		fphase.close()


		if len(cumulative_err_msg)>0:
			return self.CheckReturn(self.CHECK_FALSE,"\n".join(cumulative_err_msg))


		return self.CheckReturn(self.CHECK_TRUE)

	@check_message("Checking CO core radius")
	def check_rco(self, folder_path):

		if not os.path.isfile(folder_path+'/rco.dat'):
			return  self.CheckReturn(self.CHECK_TRUE)

		if not os.path.isfile(folder_path+'/rhe.dat'):
			fhe=None
		else:
			fhe=open(folder_path+"/rhe.dat","r")

		fmass=open(folder_path+"/mass.dat","r")
		fmhe=open(folder_path+"/mhe.dat","r")
		fmco=open(folder_path+"/mco.dat","r")
		frco=open(folder_path+"/rco.dat","r")
		frad=open(folder_path+"/radius.dat","r")
		ftime=open(folder_path+"/time.dat","r")
		fphase=open(folder_path+"/phase.dat","r")

		cumulative_err_msg=[]
		row=1

		if fhe is not None:
			for masses, comasses, hemasses, rads, rcos, rhes, times, phases in zip(self.read_lines(fmass), self.read_lines(fmco), self.read_lines(fmhe), self.read_lines(frad), self.read_lines(frco), self.read_lines(fhe) , self.read_lines(ftime), self.read_phases(fphase)):

				problematic_idx_typea = (rcos>=rads) & (comasses!=masses)
				problematic_idx_typeb = (rcos!=rads) & (comasses==masses)
				problematic_idx_typec = (rcos>=rhes) & (comasses!=hemasses) & (problematic_idx_typea==False)
				problematic_idx_typed = (rcos!=rhes) & (comasses==hemasses) & (problematic_idx_typeb==False)

				col_list=np.arange(len(masses))

				if (np.sum(problematic_idx_typea)>0 or np.sum(problematic_idx_typeb)>0 or np.sum(problematic_idx_typec) or np.sum(problematic_idx_typed)):
					msg="\nProblem for \033[96m Zams  %s \033[0m (row %i): \n"%(masses[0],row)

					for time,rad,rco,mass,mco,col in zip(times[problematic_idx_typea],rads[problematic_idx_typea],rcos[problematic_idx_typea],masses[problematic_idx_typea],comasses[problematic_idx_typea],col_list[problematic_idx_typea]):
						msg+=" - RCO(%f)>=Radius(%f) when MCO(%f)!=Mass(%f) at t=%f (col=%i/%i) phase=%i\n"%(rco,rad,mco,mass,time,col,len(masses)-1,self.find_phase(time,phases))

					for time,rad,rco,mass,mco,col in zip(times[problematic_idx_typeb],rads[problematic_idx_typeb],rcos[problematic_idx_typeb],masses[problematic_idx_typeb],comasses[problematic_idx_typeb],col_list[problematic_idx_typeb]):
						msg+=" - RCO%f)!=Radius(%f) when MCO(%f)==Mass(%f) at t=%f (col=%i/%i) phase=%i\n"%(rco,rad,mco,mass,time,col,len(masses)-1,self.find_phase(time,phases))

					for time,rhe,rco,mhe,mco,col in zip(times[problematic_idx_typec],rhes[problematic_idx_typec],rcos[problematic_idx_typec],hemasses[problematic_idx_typec],comasses[problematic_idx_typec],col_list[problematic_idx_typec]):
						msg+=" - RCO(%f)>=RHE(%f) when MCO(%f)!=MHE(%f) at t=%f (col=%i/%i) phase=%i\n"%(rco,rhe,mco,mhe,time,col,len(masses)-1,self.find_phase(time,phases))

					for time,rhe,rco,mhe,mco,col in zip(times[problematic_idx_typed],rhes[problematic_idx_typed],rcos[problematic_idx_typed],hemasses[problematic_idx_typed],comasses[problematic_idx_typed],col_list[problematic_idx_typed]):
						msg+=" - RCO(%f)!=Radius(%f) when MCO(%f)==Mass(%f) at t=%f (col=%i/%i) phase=%i\n"%(rco,rhe,mco,mhe,time,col,len(masses)-1,self.find_phase(time,phases))


					cumulative_err_msg.append(msg[:-1]) #-1 to remove the last/n

				row+=1

		else:
			for masses, comasses, hemasses, rads, rcos, rhes, times, phases in zip(self.read_lines(fmass), self.read_lines(fmco), self.read_lines(fmhe), self.read_lines(frad), self.read_lines(frco), self.read_lines(fhe) , self.read_lines(ftime), self.read_phases(fphase)):

				problematic_idx_typea = (rcos>=rads) & (comasses!=masses)
				problematic_idx_typeb = (rcos!=rads) & (comasses==masses)

				col_list=np.arange(len(masses))

				if (np.sum(problematic_idx_typea)>0 or np.sum(problematic_idx_typeb)>0):
					msg="\nProblem for \033[96m Zams  %s \033[0m (row %i): \n"%(masses[0],row)

					for time,rad,rco,mass,mco,col in zip(times[problematic_idx_typea],rads[problematic_idx_typea],rcos[problematic_idx_typea],masses[problematic_idx_typea],comasses[problematic_idx_typea],col_list[problematic_idx_typea]):
						msg+=" - RCO(%f)>=Radius(%f) when MCO(%f)!=Mass(%f) at t=%f (col=%i/%i) phase=%i\n"%(rco,rad,mco,mass,time,col,len(masses)-1,self.find_phase(time,phases))

					for time,rad,rco,mass,mco,col in zip(times[problematic_idx_typeb],rads[problematic_idx_typeb],rcos[problematic_idx_typeb],masses[problematic_idx_typeb],comasses[problematic_idx_typeb],col_list[problematic_idx_typeb]):
						msg+=" - RCO%f)!=Radius(%f) when MCO(%f)==Mass(%f) at t=%f (col=%i/%i) phase=%i\n"%(rco,rad,mco,mass,time,col,len(masses)-1,self.find_phase(time,phases))

					cumulative_err_msg.append(msg[:-1]) #-1 to remove the last/n

			row+=1


		frad.close()
		fmass.close()
		fmhe.close()
		frco.close()
		if fhe!=None: fhe.close()
		ftime.close()
		fphase.close()


		if len(cumulative_err_msg)>0:
			return self.CheckReturn(self.CHECK_WARNING,"\n".join(cumulative_err_msg))

		return  self.CheckReturn(self.CHECK_TRUE)

	@check_message("Checking masses")
	def  check_mass(self,folder_path):

		fmass=open(folder_path+"/mass.dat","r")

		cumulative_msg=[]
		for row, masses in enumerate(self.read_lines(fmass)):

			cumulative_msg_inner=[]
			for col in range(1,len(masses),1):
				if masses[col]>masses[col-1]:
					cumulative_msg_inner.append(" - M[%i](%f)>M[%i](%f)"%(col,masses[col],col-1,masses[col-1]))

			if len(cumulative_msg_inner)>0:
				msg="\nProblem for \033[96m Zams  %s \033[0m (row %i), increasing mass in given track:"%(masses[0],row+1)
				cumulative_msg = [msg, ] + cumulative_msg_inner

			if row==0:
				old_zams=masses[0]
			else:
				if masses[0]<=old_zams:
					msg="\nProblem for \033[96m Zams  %s \033[0m (row %i), Zams are not sorted:"%(masses[0],row+1)
					cumulative_msg.append(msg)
					cumulative_msg.append(" - ZAMS[%i](%f)>ZAMS[%i](%f)"%(row-1,old_zams,row,masses[0]))

				old_zams=masses[0]

		fmass.close()


		if len(cumulative_msg)>0:
			return self.CheckReturn(self.CHECK_FALSE,"\n".join(cumulative_msg))



		return  self.CheckReturn(self.CHECK_TRUE)

	@check_message("Checking times")
	def check_times(self,folder_path):

		ftimes=open(folder_path+"/time.dat","r")
		fmass=open(folder_path+"/mass.dat","r")
		fphase=open(folder_path+"/phase.dat","r")

		cumulative_msg=[]
		row=1
		for  times,zams,phases in zip(self.read_lines(ftimes),self.read_first(fmass),self.read_phases(fphase)):

			col_list = np.arange(1,len(times),1)
			idx_time_error = times[1:]<=times[:-1]
			last_time_error=times[-1]<=phases[6]

			if (np.sum(idx_time_error)>0 or last_time_error):
				cumulative_msg.append("Problem for \033[96m Zams  %s \033[0m (row %i):"%(zams,row))

				if (last_time_error):
					cumulative_msg.append(" - T_last(%f)<=Tphase[6](%f)"%(times[-1],phases[6]))


				for col in col_list[idx_time_error]:
					cumulative_msg.append(" - t[%i](%f)>=t[%i](%f)"%(col,times[col],col-1,times[col-1]))




			row+=1

		fphase.close()
		ftimes.close()
		fmass.close()

		if len(cumulative_msg)>0:

			return self.CheckReturn(self.CHECK_FALSE,"\n"+"\n".join(cumulative_msg))

		return self.CheckReturn(self.CHECK_TRUE)


	def get_all_zams(self,folder_path):

		return np.atleast_1d(np.loadtxt(folder_path+"/mass.dat",usecols=0))

	def info_zams(self,folder_path):

		lzams = self.get_all_zams(folder_path)
		Dzams=lzams[1:]-lzams[:-1]
		print("Min Zams=%.3f Max Zams=%.3f  NZams=%i"%(np.min(lzams),np.max(lzams),len(lzams)))
		if len(lzams)>1:
			print("Shortest zams interval=%.3f  Longest zams interval=%.3f"%(np.min(Dzams),np.max(Dzams)))

	def markdown_table(self,ouput_file=None):

		msg=[]
		msg.append("|   Z  | Mzams_min | Mzams_max | Ntrack |")
		msg.append("|:----:|:---------:|:---------:|:------:|")
		for folder in self.folders:
			Z=folder.split("/")[-1]
			Z=float(Z[0]+"."+Z[1:])
			lzams = self.get_all_zams(folder)
			msg.append("| %f |     %.3f     |     %.3f      |   %i    |"%(Z,np.min(lzams),np.max(lzams),len(lzams)))

		final_msg="\n".join(msg)

		if ouput_file is not None:
			fo=open(ouput_file,"w")
			print(final_msg,file=fo)
			fo.close()

		return final_msg

	def plot_all_final_masses(self,output_folder=None):


		Nz=len(self.folders)
		fig, axl = self.plot_organizer(Nz)

		i=0
		for folder,ax in zip(self.folders,axl):
			Z = folder.split("/")[-1]
			ax.set_title("Z=%s"%(Z[0]+"."+Z[1:]))
			self.plot_final_masses(folder,ax)
			if i==0:
				ax.legend()

			i+=1


		for ax in axl[i:]:
			plt.sca(ax)
			plt.axis('off')

		fig.tight_layout()

		if output_folder is None:
			fig.savefig("All_mfinal.png")
		else:
			fig.savefig(output_folder+"/All_mfinal.png")

	def plot_all_zams(self,output_folder=None):

		fig, ax = plt.subplots(1,1, figsize=(8,8))
		plt.sca(ax)

		Zmin=1e30
		Zmax=1e-30
		for i,folder in enumerate(self.folders):
			Z=folder.split("/")[-1]
			Z=float(Z[0]+"."+Z[1:])
			if (Z<Zmin): Zmin=Z
			if (Z>Zmax): Zmax=Z
			lzams = self.get_all_zams(folder)
			plt.scatter(np.ones_like(lzams)*Z,lzams,c="k")

		plt.xlabel("Z")
		plt.ylabel("Mzams [Msun]")
		plt.xlim(Zmin-0.1*Zmin,Zmax+0.1*Zmax)
		plt.xscale("log")
		plt.yscale("log")

		fig.tight_layout()

		if output_folder is None:
			fig.savefig("All_Zams.png")
		else:
			fig.savefig(output_folder+"/All_Zams.png")

	def plot_all_tracks(self,output_folder=None):

		Nz=len(self.folders)
		fig, axl = self.plot_organizer(Nz)

		i=0
		for folder,ax in zip(self.folders,axl):
			Z = folder.split("/")[-1]
			ax.set_title("Z=%s"%(Z[0]+"."+Z[1:]))
			self.plot_tracks(folder,ax)

			i+=1


		for ax in axl[i:]:
			plt.sca(ax)
			plt.axis('off')

		fig.tight_layout()

		if output_folder is None:
			fig.savefig("All_tracks.png")
		else:
			fig.savefig(output_folder+"/All_tracks.png")

	def plot_all_tracks_phases(self,output_folder=None):

		if not os.path.exists(output_folder):
			os.makedirs(output_folder)

		i=0
		for folder in self.folders:
			self.plot_track_phases(folder,output_folder)

			i+=1

	def plot_all_tlife(self,ouput_folder=None):

		Nz=len(self.folders)
		fig, axl = self.plot_organizer(Nz)

		i=0
		for folder,ax in zip(self.folders,axl):
			Z = folder.split("/")[-1]
			ax.set_title("Z=%s"%(Z[0]+"."+Z[1:]))
			self.plot_tlife(folder,ax)
			if i==0:
				ax.legend()

			i+=1


		for ax in axl[i:]:
			plt.sca(ax)
			plt.axis('off')

		fig.tight_layout()

		if output_folder is None:
			fig.savefig("All_tlife.png")
		else:
			fig.savefig(output_folder+"/All_tlife.png")


	def plot_final_masses(self,folder_path,ax):


		lzams=self.get_all_zams(folder_path)
		fmass=open(folder_path+"/mass.dat","r")
		fco=open(folder_path+"/mco.dat","r")
		fhe=open(folder_path+"/mhe.dat","r")


		farray = np.array([ [zams, massf, mhef, mcof] for zams, massf, mhef, mcof in zip(lzams,self.read_last(fmass),self.read_last(fhe),self.read_last(fco))])
		plt.sca(ax)

		if len(farray)==1:
			plt.scatter(farray[:,0],farray[:,1],label="Mtot")
			plt.scatter(farray[:,0],farray[:,2],label="MHE",ls="dashed")
			plt.scatter(farray[:,0],farray[:,3],label="MCO",ls="dotted")
			xmin=farray[0,0]*0.9
			xmax=farray[0,0]*1.1
			ymin=farray[0,3]*0.9
			ymax=farray[0,3]*1.1
		else:
			plt.plot(farray[:,0],farray[:,1],label="Mtot")
			plt.plot(farray[:,0],farray[:,2],label="MHE",ls="dashed")
			plt.plot(farray[:,0],farray[:,3],label="MCO",ls="dotted")
			xmin=np.min(farray[:,0])
			xmax=np.max(farray[:,0])
			ymin=np.min(farray[:,3])
			ymax=np.max(farray[:,0])

		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)
		plt.xlabel("Mzams [Msun]")
		plt.ylabel("Mfinal [Msun]")
		plt.xscale("log")
		plt.yscale("log")



		fmass.close()
		fco.close()
		fhe.close()

	def plot_track_phases(self, folder_path, output_folder=None):

		Z = folder_path.split("/")[-1]

		lzams=self.get_all_zams(folder_path)
		frad=open(folder_path+"/radius.dat","r")
		flumi=open(folder_path+"/lumi.dat","r")
		ftime=open(folder_path+"/time.dat")
		fphase=open(folder_path+"/phase.dat")

		Nplot=16
		dzams=np.ceil(len(lzams)/Nplot)
		idx_plot=np.arange(0,len(lzams)+dzams,dzams)


		SigmaStefBoltz=7.144796315707217e-17 #LSun^3 RSun^-2 K^-4
		#cmap = matplotlib.cm.get_cmap('plasma_r')
		#norm = matplotlib.colors.LogNorm(vmin=minm, vmax=maxm)


		fig, axl = self.plot_organizer(Nplot)

		color_phase={0: "brown", 1: "gold" , 2: "orange", 3: "red", 4: "darkgreen", 5: "blue", 6: "pink"}

		i=0
		row=0
		for mass, radii, lums, times, dphase in zip(lzams, self.read_lines(frad), self.read_lines(flumi),self.read_lines(ftime),self.read_phases(fphase)):

			if row  in idx_plot:
				plt.sca(axl[i])
				temp =  (lums/(4*np.pi*radii*radii*SigmaStefBoltz))**0.25
				phases=self.find_phase(times,dphase)
				idx_phase=phases>0

				for phase in np.sort(np.unique(phases[idx_phase])):
					cphase=phases==phase
					if phase==5:
						plt.plot(temp[cphase],lums[cphase],color=color_phase[phase],label=phase,zorder=10000)
					else:
						plt.plot(temp[cphase],lums[cphase],color=color_phase[phase],label=phase)

				if i==0:
					plt.legend()

				plt.title("M$_{ZAMS}$=%f"%mass)
				plt.xlabel("T [K]")
				plt.ylabel("L [Lsun]")

				plt.gca().invert_xaxis()
				plt.xscale("log")
				plt.yscale("log")

				i+=1

				#plt.plot(temp,lums,color=cmap(norm(mass)),lw=0.5)
				#plt.scatter(temp[idx_phase],lums[idx_phase],c=phases[idx_phase])
			row+=1

		for ax in axl[i:]:
			plt.sca(ax)
			plt.axis('off')


		fig.tight_layout()

		frad.close()
		flumi.close()
		ftime.close()
		fphase.close()

		if output_folder is None:
			fig.savefig(Z+"_All_tracks_phases.png")
		else:
			fig.savefig(output_folder+"/"+Z+"_All_tracks_phases.png")

	def plot_tracks(self,folder_path,ax):

		lzams=self.get_all_zams(folder_path)
		frad=open(folder_path+"/radius.dat","r")
		flumi=open(folder_path+"/lumi.dat","r")
		if (len(lzams)==1):
			minm=0.9*lzams[0]
			maxm=1.1*lzams[0]
		else:
			minm=np.min(lzams)
			maxm=np.max(lzams)

		SigmaStefBoltz=7.144796315707217e-17 #LSun^3 RSun^-2 K^-4
		cmap = matplotlib.cm.get_cmap('plasma_r')
		norm = matplotlib.colors.LogNorm(vmin=minm, vmax=maxm)
		plt.sca(ax)

		for mass, radii, lums in zip(lzams, self.read_lines(frad), self.read_lines(flumi)):
			temp =  (lums/(4*np.pi*radii*radii*SigmaStefBoltz))**0.25
			plt.plot(temp,lums,color=cmap(norm(mass)),lw=0.5)

		plt.scatter(np.ones_like(lzams)*np.mean(temp),np.ones_like(lzams)*np.mean(lums),c=lzams,cmap="plasma_r",s=0,vmin=minm,vmax=maxm,norm=matplotlib.colors.LogNorm())
		cb=plt.colorbar(pad=0)
		cb.set_label('Mzams [Msun]')



		plt.xlabel("T [K]")
		plt.ylabel("L [Lsun]")

		plt.gca().invert_xaxis()
		plt.xscale("log")
		plt.yscale("log")

		frad.close()
		flumi.close()

	def plot_tlife(self, folder_path,ax):

		lzams=self.get_all_zams(folder_path)
		ftimes=open(folder_path+"/time.dat","r")


		farray = np.array([ [zams, times] for zams, times in zip(lzams,self.read_last(ftimes))])
		plt.sca(ax)

		if len(farray)==1:
			plt.scatter(farray[:,0],farray[:,1])
			xmin=farray[0,0]*0.9
			xmax=farray[0,0]*1.1
			ymin=farray[0,1]*0.9
			ymax=farray[0,1]*1.1
		else:
			plt.plot(farray[:,0],farray[:,1])
			xmin=np.min(farray[:,0])
			xmax=np.max(farray[:,0])
			ymin=np.min(farray[:,1])
			ymax=np.max(farray[:,1])

		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)
		plt.xlabel("Mzams [Msun]")
		plt.ylabel("Tlife [Myr]")
		plt.xscale("log")
		plt.yscale("log")

		ftimes.close()


	def plot_organizer(self,nplot):

		nx = int(np.sqrt(nplot))
		ny = nplot//nx + nplot%nx

		fig, ax = plt.subplots(nx,ny,figsize=(4*ny,4*nx))
		ax = np.atleast_1d(ax)


		return fig, ax.flatten()

	def read_tables(self,folder_path):
		return np.sort(glob.glob(folder_path + "/*"))

	def read_tables_onlyname(self,folder_path):
		return [table.split("/")[-1] for table in self.read_tables(folder_path)]

	def read_lines(self,file):

		for line in file.readlines():

			yield self.make_vector(line)

	def read_first(self,file):

		for line in file.readlines():
			vec = self.make_vector(line)
			yield vec[0]

	def read_last(self,file):

		for line in file.readlines():
			vec = self.make_vector(line)
			yield vec[-1]

	def read_phases(self,file):

		for line in file.readlines():

			yield self.make_dphase_vec(line)

	def make_vector(self,vector_str):

		return np.fromstring(vector_str,sep="    ")

	def make_dphase_vec(self,dphase_str):

		dphase=self.make_vector(dphase_str)

		return dict((int(id), float(tt)) for id,tt in zip(dphase[1::2],dphase[:-1:2]))

	def find_phase(self,t,dphase):

		return np.searchsorted(list(dphase.values()),t,side="right")-1

class SEVNTablesValidatorPureHE (SEVNTablesValidator) :

	def __init__(self,folder_path):

		super().__init__(folder_path)
		self.len_phases=6

	@check_message("Checking He core mass time")
	def check_he_core_time(self,folder_path):

		fmass=open(folder_path+"/mass.dat","r")
		fhe=open(folder_path+"/mhe.dat","r")
		ftime=open(folder_path+"/time.dat","r")
		fphase=open(folder_path+"/phase.dat","r")

		cumulative_err_msg=[]

		row=1
		for mass, hemasses, times, phases in zip(self.read_first(fmass), self.read_lines(fhe), self.read_lines(ftime), self.read_phases(fphase)):

			col_list=np.arange(len(times))

			times_without_he = times[hemasses == 0]
			t_he_start = 0

			he_wrongly_zero = times_without_he > t_he_start
			t_he_wrongly_zero = times_without_he[he_wrongly_zero]


			if (len(t_he_wrongly_zero)>0):

				msg="\nProblem for \033[96m Zams  %s \033[0m (row %i, tphase_2=%f): \n"%(mass,row,t_he_start)

				for  i,t_error in enumerate(t_he_wrongly_zero):
					msg+=" - MHE=0 after tphase_2 at t=%f phase=%i\n"%(t_error,self.find_phase(t_error,phases))

				cumulative_err_msg.append(msg[:-1]) #-1 to remove the last/n


			row+=1

		fmass.close()
		fhe.close()
		ftime.close()
		fphase.close()


		if len(cumulative_err_msg)>0:
			return self.CheckReturn(self.CHECK_FALSE,"\n".join(cumulative_err_msg))

		return  self.CheckReturn(self.CHECK_TRUE)


	@check_message("Checking HE core radius times")
	def check_rhe_time(self,folder_path):

		if not os.path.isfile(folder_path+'/rhe.dat'):

			return  self.CheckReturn(self.CHECK_TRUE)

		else:

			fmass=open(folder_path+"/mass.dat","r")
			fhe=open(folder_path+"/rhe.dat","r")
			ftime=open(folder_path+"/time.dat","r")
			fphase=open(folder_path+"/phase.dat","r")

			cumulative_err_msg=[]

			row=1
			for mass, hemasses, times, phases in zip(self.read_first(fmass), self.read_lines(fhe), self.read_lines(ftime), self.read_phases(fphase)):

				times_without_he = times[hemasses == 0]
				t_he_start = 0


				he_wrongly_zero = times_without_he > t_he_start
				t_he_wrongly_zero = times_without_he[he_wrongly_zero]


				if ( len(t_he_wrongly_zero)>0):

					msg="\nProblem for \033[96m Zams  %s \033[0m (row %i, tphase_2=%f): \n"%(mass,row,t_he_start)

					for  i,t_error in enumerate(t_he_wrongly_zero):
						msg+=" - RHE=0 after tphase_2 at t=%f phase=%i\n"%(t_error,self.find_phase(t_error,phases))

					cumulative_err_msg.append(msg[:-1]) #-1 to remove the last/n


				row+=1

			fmass.close()
			fhe.close()
			ftime.close()
			fphase.close()

			if len(cumulative_err_msg)>0:
				return self.CheckReturn(self.CHECK_FALSE,"\n".join(cumulative_err_msg))


		return  self.CheckReturn(self.CHECK_TRUE)


	@check_message("Checking He core mass")
	def check_he_core_mass(self,folder_path):

		fmass=open(folder_path+"/mass.dat","r")
		fhe=open(folder_path+"/mhe.dat","r")
		ftime=open(folder_path+"/time.dat","r")
		fphase=open(folder_path+"/phase.dat","r")

		cumulative_err_msg=[]

		row=1
		for masses, hemasses, times, phases in zip(self.read_lines(fmass), self.read_lines(fhe), self.read_lines(ftime), self.read_phases(fphase)):

			problematic_idx = hemasses>masses
			col_list=np.arange(len(masses))

			if (np.sum(problematic_idx)>0):
				msg="\nProblem for \033[96m Zams  %s \033[0m (row %i): \n"%(masses[0],row)

				for time,mass,mhe,col in zip(times[problematic_idx],masses[problematic_idx],hemasses[problematic_idx],col_list[problematic_idx]):
					msg+=" - MHE(%f)>Mass(%f) at t=%f (col=%i/%i) phase=%i\n"%(mhe,mass,time,col,len(masses)-1,self.find_phase(time,phases))


				cumulative_err_msg.append(msg[:-1]) #-1 to remove the last/n


			row+=1

		fmass.close()
		fhe.close()
		ftime.close()
		fphase.close()


		if len(cumulative_err_msg)>0:
			return self.CheckReturn(self.CHECK_WARNING,"\n".join(cumulative_err_msg))

		return  self.CheckReturn(self.CHECK_TRUE)

	def find_phase(self,t,dphase):

			return (np.searchsorted(list(dphase.values()),t,side="right")-1)+4

###### Global Variables
global_fundamental_tables=("mass.dat","mhe.dat","mco.dat","radius.dat","lumi.dat","time.dat","phase.dat")

@check_message("Checking folder name")
def check_metallicity_folder_name(folder_path):

	name = folder_path.split("/")[-1]

	err_msg=""
	if not name.isdigit() or name[0]!="0":
		msg_metallicity = "Metallicity folder name must contains only digits and start with 0"
		err_msg="Invalid folder name %s. %s"%(name,msg_metallicity)
		return False, err_msg

	return True,err_msg

@check_message("Checking fundamental tables")
def check_fundamental_tables(folder_path):

	err_msg = ""

	missing_tables=[]
	for table_name in global_fundamental_tables:
		if not os.path.isfile(folder_path+"/"+table_name):
			missing_tables.append(table_name)

	if len(missing_tables)>0:
		err_msg = "Fundamental tables are missing: %s"%(", ".join(missing_tables))
		return False, err_msg

	return True,err_msg

def make_check(check_function,message):
	print("*"+message+"...", end="", flush=True)

if __name__=="__main__":

	if __name__ == "__main__":

		if (len(sys.argv) == 1):
			raise ValueError("path to the ouput folder missing")

		input_folder = sys.argv[1]
		if input_folder[-1]=="/":
			input_folder=input_folder[:-1]

		make_plot=False
		if "--plot" in sys.argv:
			make_plot=True

		pureHE=False
		if "--purehe" in sys.argv:
			pureHE=True


		if ("/" in input_folder):
			output_folder = "tables_analysis_"+input_folder.split("/")[-1]
		else:
			output_folder = "tables_analysis_"+input_folder

		if not os.path.exists(output_folder):
			os.makedirs(output_folder)

		if pureHE:
			svt = SEVNTablesValidatorPureHE(input_folder)
		else:
			svt = SEVNTablesValidator(input_folder)
		svt.check_all()
		if make_plot:
			svt.plot_all_final_masses(output_folder)
			svt.plot_all_zams(output_folder)
			svt.plot_all_tracks(output_folder)
			svt.plot_all_tlife(output_folder)
			svt.plot_all_tracks_phases(output_folder+"/plot_phases")
		mkdown_table=svt.markdown_table(output_folder+"/markdown_info_table.txt")
		print()
		print(mkdown_table)
