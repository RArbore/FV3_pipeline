from argparse import ArgumentParser
from netCDF4 import Dataset
import numpy as np
import json

argparser = ArgumentParser(description="Destruct FV3 ADM MODTRAN output into a single netCDF file.")
argparser.add_argument("input", help="Path to directory containing MODTRAN's output.")
argparser.add_argument("output", help="Output file path.")
argparser.add_argument("cases", help="Number of cases.")
argparser.add_argument("szas", help="Number of szas.")
args = argparser.parse_args()

CASES = int(args.cases)
SZAS = int(args.szas)

szas_list = []
prof_pressure_list = []
prof_temperature_list = []
prof_h2o_list = []
prof_co2_list = []
prof_o3_list = []
prof_n2o_list = []
prof_co_list = []
prof_ch4_list = []
toa_rad_list = []
ref_sol_list = []
sol_at_obs_list = []
depth_list = []
obszen_list = []
azinp_list = []
vert_path_list = []
los_paths_list = []
geometries_vert_path_list = []
geometries_los_paths_list = []
day_list = []
temp_list = []
toa_radiances_list = []
upward_diffuse_list = []
downward_diffuse_list = []
direct_solar_list = []

flux_unit = ""
toa_radiance_unit = "W CM-2 SR-1 / CM-1"

for case in range(CASES):
	for sza in range(SZAS):
		print("Handling case {} and SZA {}.".format(case+1, sza+1))
		file_name = args.input + "/Libera_ADM_FV3_Input_Case_{:02d}_SZA_{:02d}".format(case+1, sza+1)

		json_data = json.loads(open(file_name+".json", "r").read())
		if case == 0:
			szas_list.append(json_data["MODTRAN"][0]["MODTRANINPUT"]["GEOMETRY"]["PARM2"])
		profiles = json_data["MODTRAN"][0]["MODTRANINPUT"]["ATMOSPHERE"]["PROFILES"]
		prof_pressure = np.array(profiles[1]["PROFILE"])
		prof_temperature = np.array(profiles[2]["PROFILE"])
		prof_h2o = np.array(profiles[3]["PROFILE"])
		prof_co2 = np.array(profiles[4]["PROFILE"])
		prof_o3 = np.array(profiles[5]["PROFILE"])
		prof_n2o = np.array(profiles[6]["PROFILE"])
		prof_co = np.array(profiles[7]["PROFILE"])
		prof_ch4 = np.array(profiles[8]["PROFILE"])
		radiance = json_data["MODTRAN"][0]["MODTRANOUTPUT"]["SPECTRA"]["RADIANCE"]
		toa_rad = np.array(radiance["TOA_RAD"])
		ref_sol = np.array(radiance["REF_SOL"])
		sol_at_obs = np.array(radiance["SOL_AT_OBS"])
		depth = np.array(radiance["DEPTH"])
		if case == 0 and sza == 0:
			try:
				xflag = json_data["MODTRAN"][0]["MODTRANINPUT"]["SPECTRAL"]["XFLAG"]
				if xflag == "N":
					toa_radiance_unit = "uW CM-2 SR-1 / NM"
				elif xflag == "M":
					toa_radiance_unit = "W CM-2 SR-1 / uM"
				elif xflag == "W":
					toa_radiance_unit = "W CM-2 SR-1 / CM-1"
				if xflag == "N" or xflag == "M" or xflag == "W":
					pass
					#print("Found TOA radiances unit " + toa_radiance_unit)
			except:
				pass
		mlos = json_data["MODTRAN"][0]["MODTRANINPUT"]["GEOMETRY"]["MLOS"]
		obszen = np.array(list(map(lambda x: x["OBSZEN"], mlos)))
		azinp = np.array(list(map(lambda x: x["AZ_INP"], mlos)))

		tp6_data = open(file_name+".tp6", "r").read()	

		vert_path = tp6_data.split("TOTAL COLUMN ABSORBER AMOUNTS FOR A VERTICAL PATH FROM GROUND TO SPACE:")[1].split("SINGLE SCATTER SOLAR PATH GEOMETRY TABLE FOR MULTIPLE SCATTERING VERTICAL GROUND-TO-SPACE PATH")[0]
		vert_path = list(filter(lambda x: "." in x and not "(" in x and not ")" in x, vert_path.split("\n")[1:]))
		vert_path = np.concatenate(list(map(lambda x: np.fromstring(x, sep=" "), vert_path)))

		los_paths = tp6_data.split("TOTAL COLUMN ABSORBER AMOUNTS FOR LINE-OF-SIGHT PATH")[1:]
		los_paths = list(map(lambda x: x.split("Line-of-Sight No.")[0].split("\n"), los_paths))
		los_paths = list(map(lambda x: list(filter(lambda y: "." in y and not "(" in y and not ")" in y, x)), los_paths))
		los_paths = np.stack(list(map(lambda x: np.concatenate(list(map(lambda y: np.fromstring(y, sep=" "), x))), los_paths)))

		geometries_vert_path = tp6_data.split("SINGLE SCATTER SOLAR PATH GEOMETRY TABLE FOR MULTIPLE SCATTERING VERTICAL GROUND-TO-SPACE PATH")[1].split("SOLAR ILLUMINATION PATH COLUMN")[0]
		geometries_vert_path = geometries_vert_path.split("[DEG]")[-1]
		geometries_vert_path = list(filter(lambda x: "." in x and not any([y.islower() for y in x]), geometries_vert_path.split("\n")))
		geometries_vert_path = np.stack(list(map(lambda x: np.fromstring(x, sep=" "), geometries_vert_path)))

		geometries_los_paths = tp6_data.split("SINGLE SCATTER SOLAR PATH GEOMETRY TABLE FOR LINE-OF-SIGHT PATH")[1:]
		geometries_los_paths = list(map(lambda x: x.split("SOLAR ILLUMINATION PATH COLUMN")[0].split("\n")[1:], geometries_los_paths))
		geometries_los_paths = list(map(lambda x: list(filter(lambda y: "." in y and not any([z.islower() for z in y]), x)), geometries_los_paths))
		geometries_los_paths = np.stack(list(map(lambda x: np.stack(list(map(lambda y: np.fromstring(y, sep=" ")[:7], x))), geometries_los_paths)))

		day = list(filter(lambda x: "DAY OF THE YEAR" in x, tp6_data.split("\n")))[0]
		day = int(day.split(" ")[-1])

		temp = list(filter(lambda x: "AREA-AVERAGED GROUND TEMPERATURE" in x, tp6_data.split("\n")))[0]
		temp = float(temp.split(" ")[-1])

		psc_data = open(file_name+".psc", "r").read()
		num_toa_radiances = psc_data.count("$")
		old_psc_data = psc_data.split("\n")[:-1]
		psc_data = []
		for line_num in range(0, len(old_psc_data), int(len(old_psc_data) / num_toa_radiances)):
			psc_data += old_psc_data[line_num:(line_num + int(len(old_psc_data) / num_toa_radiances - 1))]
		psc_data = "".join(psc_data)
		toa_radiances = np.transpose(np.reshape(np.fromstring(psc_data, sep=" "), (num_toa_radiances, -1, 2)), axes=(0, 2, 1))
		
		flx_data = open(file_name+".flx", "r").read()	
		fluxes = flx_data.split("WAVLEN")[-1].split("=")[0]
		fluxes = fluxes.split("\n")[:-1]
		offset = len(fluxes[3].split(" ")[1]) + 2
		fluxes = list(map(lambda x: x[offset:], fluxes[3:]))
		fluxes = np.stack(list(map(lambda x: np.fromstring(x, sep=" "), fluxes)))
		upward_diffuse = np.reshape(np.stack((fluxes[:, 0], fluxes[:, 3]), axis=1), (-1))
		downward_diffuse = np.reshape(np.stack((fluxes[:, 1], fluxes[:, 4]), axis=1), (-1))
		direct_solar = np.reshape(np.stack((fluxes[:, 2], fluxes[:, 5]), axis=1), (-1))
		upward_diffuse = np.reshape(upward_diffuse, (toa_radiances.shape[2], -1))[::-1, :]
		downward_diffuse = np.reshape(downward_diffuse, (toa_radiances.shape[2], -1))[::-1, :]
		direct_solar = np.reshape(direct_solar, (toa_radiances.shape[2], -1))[::-1, :]
		if case == 0 and sza == 0:
			unit_flx_data = list(filter(lambda x: "SPECTRAL VERTICAL FLUX TABLE" in x,flx_data.split("\n")))[0]
			flux_unit = unit_flx_data.split("(")[1].split(")")[0]		
			#print("Found flux unit " + flux_unit)

		prof_pressure_list.append(prof_pressure) # millibar
		prof_temperature_list.append(prof_temperature) # kelvin
		prof_h2o_list.append(prof_h2o) # ppmv
		prof_co2_list.append(prof_co2) # ppmv
		prof_o3_list.append(prof_o3) # ppmv
		prof_n2o_list.append(prof_n2o) # ppmv
		prof_co_list.append(prof_co) # ppmv
		prof_ch4_list.append(prof_ch4) # ppmv
		toa_rad_list.append(toa_rad)
		ref_sol_list.append(ref_sol) # W cm-2
		sol_at_obs_list.append(sol_at_obs) # W cm-2
		depth_list.append(depth)
		obszen_list.append(obszen)
		azinp_list.append(azinp)
		vert_path_list.append(vert_path)
		los_paths_list.append(los_paths)
		geometries_vert_path_list.append(geometries_vert_path)
		geometries_los_paths_list.append(geometries_los_paths)
		day_list.append(day)
		temp_list.append(temp)
		toa_radiances_list.append(toa_radiances)
		upward_diffuse_list.append(upward_diffuse)
		downward_diffuse_list.append(downward_diffuse)
		direct_solar_list.append(direct_solar)

szas = np.stack(szas_list)
prof_pressure = np.stack(prof_pressure_list)
prof_temperature = np.stack(prof_temperature_list)
prof_h2o = np.stack(prof_h2o_list)
prof_co2 = np.stack(prof_co2_list)
prof_o3 = np.stack(prof_o3_list)
prof_n2o = np.stack(prof_n2o_list)
prof_co = np.stack(prof_co_list)
prof_ch4 = np.stack(prof_ch4_list)
toa_rad = np.stack(toa_rad_list)
ref_sol = np.stack(ref_sol_list)
sol_at_obs = np.stack(sol_at_obs_list)
depth = np.stack(depth_list)
obszen = np.stack(obszen_list)
azinp = np.stack(azinp_list)
vert_path = np.stack(vert_path_list)
los_paths = np.stack(los_paths_list)
geometries_vert_path = np.stack(geometries_vert_path_list)
geometries_los_paths = np.stack(geometries_los_paths_list)
day = np.stack(day_list)
temp = np.stack(temp_list)
toa_radiances = np.stack(toa_radiances_list)
upward_diffuse = np.stack(upward_diffuse_list)
downward_diffuse = np.stack(downward_diffuse_list)
direct_solar = np.stack(direct_solar_list)

prof_pressure = np.reshape(prof_pressure, (CASES, SZAS) + prof_pressure.shape[1:])
prof_temperature = np.reshape(prof_temperature, (CASES, SZAS) + prof_temperature.shape[1:])
prof_h2o = np.reshape(prof_h2o, (CASES, SZAS) + prof_h2o.shape[1:])
prof_co2 = np.reshape(prof_co2, (CASES, SZAS) + prof_co2.shape[1:])
prof_o3 = np.reshape(prof_o3, (CASES, SZAS) + prof_o3.shape[1:])
prof_n2o = np.reshape(prof_n2o, (CASES, SZAS) + prof_n2o.shape[1:])
prof_co = np.reshape(prof_co, (CASES, SZAS) + prof_co.shape[1:])
prof_ch4 = np.reshape(prof_ch4, (CASES, SZAS) + prof_ch4.shape[1:])
toa_rad = np.reshape(toa_rad, (CASES, SZAS) + toa_rad.shape[1:])
ref_sol = np.reshape(ref_sol, (CASES, SZAS) + ref_sol.shape[1:])
sol_at_obs = np.reshape(sol_at_obs, (CASES, SZAS) + sol_at_obs.shape[1:])
depth = np.reshape(depth, (CASES, SZAS) + depth.shape[1:])
obszen = np.reshape(obszen, (CASES, SZAS) + obszen.shape[1:])
azinp = np.reshape(azinp, (CASES, SZAS) + azinp.shape[1:])
vert_path = np.reshape(vert_path, (CASES, SZAS) + vert_path.shape[1:])
los_paths = np.reshape(los_paths, (CASES, SZAS) + los_paths.shape[1:])
geometries_vert_path = np.reshape(geometries_vert_path, (CASES, SZAS) + geometries_vert_path.shape[1:])
geometries_los_paths = np.reshape(geometries_los_paths, (CASES, SZAS) + geometries_los_paths.shape[1:])
day = np.reshape(day, (CASES, SZAS))
temp = np.reshape(temp, (CASES, SZAS))
toa_radiances = np.reshape(toa_radiances, (CASES, SZAS) + toa_radiances.shape[1:])
upward_diffuse = np.reshape(upward_diffuse, (CASES, SZAS) + upward_diffuse.shape[1:])
downward_diffuse = np.reshape(downward_diffuse, (CASES, SZAS) + downward_diffuse.shape[1:])
direct_solar = np.reshape(direct_solar, (CASES, SZAS) + direct_solar.shape[1:])

def spectrally_integrate(array, lam, dim):
	b_shape = np.ones((len(array.shape),), dtype=np.int)
	b_shape[dim] = lam.shape[0]
	lam = np.reshape(lam, b_shape)
	lam_2 = lam * lam
	integral = array / lam_2
	return np.trapz(integral, lam, axis=dim)

toa_radiances_wavelens_test = toa_radiances[:, :, :, 0, :]
toa_radiances_wavelens_test = np.reshape(toa_radiances_wavelens_test, (-1, toa_radiances.shape[4],))
for i in range(toa_radiances_wavelens_test.shape[0] - 1):
	assert (toa_radiances_wavelens_test[i, :] == toa_radiances_wavelens_test[i + 1, :]).all()
vis_cutoff = np.reshape(np.where(toa_radiances_wavelens_test[0, :] <= 700.0), (-1,))[-1]
nir_cutoff = np.reshape(np.where(toa_radiances_wavelens_test[0, :] >= 700.0), (-1,))[0]
vis_lam = toa_radiances_wavelens_test[0, :(vis_cutoff+1)]
nir_lam = toa_radiances_wavelens_test[0, nir_cutoff:]
toa_radiances_vis = spectrally_integrate(toa_radiances[:, :, :, 1, :(vis_cutoff+1)], vis_lam, 3)
toa_radiances_nir = spectrally_integrate(toa_radiances[:, :, :, 1, nir_cutoff:], nir_lam, 3)
upward_diffuse_vis = spectrally_integrate(upward_diffuse[:, :, :(vis_cutoff+1), :], vis_lam, 2)
upward_diffuse_nir = spectrally_integrate(upward_diffuse[:, :, nir_cutoff:, :], nir_lam, 2)
downward_diffuse_vis = spectrally_integrate(downward_diffuse[:, :, :(vis_cutoff+1), :], vis_lam, 2)
downward_diffuse_nir = spectrally_integrate(downward_diffuse[:, :, nir_cutoff:, :], nir_lam, 2)
direct_solar_vis = spectrally_integrate(direct_solar[:, :, :(vis_cutoff+1), :], vis_lam, 2)
direct_solar_nir = spectrally_integrate(direct_solar[:, :, nir_cutoff:, :], nir_lam, 2)

netcdf = Dataset(args.output, "w", format="NETCDF4")

netcdf.createDimension("num_cases", CASES)
netcdf.createDimension("num_szas", SZAS)
netcdf.createDimension("profiles_len", prof_pressure.shape[2])
netcdf.createDimension("num_wavenumbers", toa_rad.shape[2])
netcdf.createDimension("num_input_geometries", obszen.shape[2])
netcdf.createDimension("num_path_geometries", geometries_vert_path.shape[2])
netcdf.createDimension("num_toa_radiances", num_toa_radiances)
netcdf.createDimension("num_wavelengths", toa_radiances.shape[4])
netcdf.createDimension("fluxes_per_index", upward_diffuse.shape[3])

szas_var = netcdf.createVariable("szas", "f8", ("num_szas",))
prof_pressure_var = netcdf.createVariable("prof_pressure", "f8", ("num_cases", "num_szas", "profiles_len",))
prof_temperature_var = netcdf.createVariable("prof_temperature", "f8", ("num_cases", "num_szas", "profiles_len",))
prof_h2o_var = netcdf.createVariable("prof_h2o", "f8", ("num_cases", "num_szas", "profiles_len",))
prof_co2_var = netcdf.createVariable("prof_co2", "f8", ("num_cases", "num_szas", "profiles_len",))
prof_o3_var = netcdf.createVariable("prof_o3", "f8", ("num_cases", "num_szas", "profiles_len",))
prof_n2o_var = netcdf.createVariable("prof_n2o", "f8", ("num_cases", "num_szas", "profiles_len",))
prof_co_var = netcdf.createVariable("prof_co", "f8", ("num_cases", "num_szas", "profiles_len",))
prof_ch4_var = netcdf.createVariable("prof_ch4", "f8", ("num_cases", "num_szas", "profiles_len",))
toa_rad_var = netcdf.createVariable("toa_rad", "f8", ("num_cases", "num_szas", "num_wavenumbers",))
ref_sol_var = netcdf.createVariable("ref_sol", "f8", ("num_cases", "num_szas", "num_wavenumbers",))
sol_at_obs_var = netcdf.createVariable("sol_at_obs", "f8", ("num_cases", "num_szas", "num_wavenumbers",))
depth_var = netcdf.createVariable("depth", "f8", ("num_cases", "num_szas", "num_wavenumbers",))
obszen_var = netcdf.createVariable("obszen", "f8", ("num_cases", "num_szas", "num_input_geometries",))
azinp_var = netcdf.createVariable("azinp", "f8", ("num_cases", "num_szas", "num_input_geometries",))
day_var = netcdf.createVariable("day", "u8", ("num_cases", "num_szas",))
temp_var = netcdf.createVariable("temp", "f8", ("num_cases", "num_szas",))
toa_radiances_wavelengths_var = netcdf.createVariable("toa_radiances_wavelengths", "f8", ("num_cases", "num_szas", "num_toa_radiances", "num_wavelengths",))
toa_radiances_var = netcdf.createVariable("toa_radiances", "f8", ("num_cases", "num_szas", "num_toa_radiances", "num_wavelengths",))
toa_radiances_var_vis = netcdf.createVariable("toa_radiances_vis", "f8", ("num_cases", "num_szas", "num_toa_radiances",))
toa_radiances_var_nir = netcdf.createVariable("toa_radiances_nir", "f8", ("num_cases", "num_szas", "num_toa_radiances",))
upward_diffuse_var = netcdf.createVariable("upward_diffuse", "f8", ("num_cases", "num_szas", "num_wavelengths", "fluxes_per_index",))
upward_diffuse_var_vis = netcdf.createVariable("upward_diffuse_vis", "f8", ("num_cases", "num_szas", "fluxes_per_index",))
upward_diffuse_var_nir = netcdf.createVariable("upward_diffuse_nir", "f8", ("num_cases", "num_szas", "fluxes_per_index",))
downward_diffuse_var = netcdf.createVariable("downward_diffuse", "f8", ("num_cases", "num_szas", "num_wavelengths", "fluxes_per_index",))
downward_diffuse_var_vis = netcdf.createVariable("downward_diffuse_vis", "f8", ("num_cases", "num_szas", "fluxes_per_index",))
downward_diffuse_var_nir = netcdf.createVariable("downward_diffuse_nir", "f8", ("num_cases", "num_szas", "fluxes_per_index",))
direct_solar_var = netcdf.createVariable("direct_solar", "f8", ("num_cases", "num_szas", "num_wavelengths", "fluxes_per_index",))
direct_solar_var_vis = netcdf.createVariable("direct_solar_vis", "f8", ("num_cases", "num_szas", "fluxes_per_index",))
direct_solar_var_nir = netcdf.createVariable("direct_solar_nir", "f8", ("num_cases", "num_szas", "fluxes_per_index",))

szas_var.units = "DEGREES"
prof_pressure_var.units = "MBAR"
prof_temperature_var.units = "KELVIN"
prof_h2o_var.units = "PPMV"
prof_co2_var.units = "PPMV"
prof_o3_var.units = "PPMV"
prof_n2o_var.units = "PPMV"
prof_co_var.units = "PPMV"
prof_ch4_var.units = "PPMV"
toa_rad_var.units = "W CM-2 SR-1 / CM-1"
ref_sol_var.units = "W CM-2"
sol_at_obs_var.units = "W CM-2"
depth_var.units = "UNITLESS"
obszen_var.units = "DEGREES"
azinp_var.units = "DEGREES"
day_var.units = "DAYS"
temp_var.units = "KELVIN"
toa_radiances_wavelengths_var.units = "NM"
toa_radiances_var.units = toa_radiance_unit
upward_diffuse_var.units = flux_unit
downward_diffuse_var.units = flux_unit
direct_solar_var.units = flux_unit

szas_var[:] = szas
prof_pressure_var[:] = prof_pressure
prof_temperature_var[:] = prof_temperature
prof_h2o_var[:] = prof_h2o
prof_co2_var[:] = prof_co2
prof_o3_var[:] = prof_o3
prof_n2o_var[:] = prof_n2o
prof_co_var[:] = prof_co
prof_ch4_var[:] = prof_ch4
toa_rad_var[:] = toa_rad
ref_sol_var[:] = ref_sol
sol_at_obs_var[:] = sol_at_obs
depth_var[:] = depth
obszen_var[:] = obszen
azinp_var[:] = azinp
day_var[:] = day
temp_var[:] = temp
toa_radiances_wavelengths_var[:] = toa_radiances[:, :, :, 0, :]
toa_radiances_var[:] = toa_radiances[:, :, :, 1, :]
upward_diffuse_var[:] = upward_diffuse
downward_diffuse_var[:] = downward_diffuse
direct_solar_var[:] = direct_solar
toa_radiances_var_vis[:] = toa_radiances_vis
upward_diffuse_var_vis[:] = upward_diffuse_vis
downward_diffuse_var_vis[:] = downward_diffuse_vis
direct_solar_var_vis[:] = direct_solar_vis
toa_radiances_var_nir[:] = toa_radiances_nir
upward_diffuse_var_nir[:] = upward_diffuse_nir
downward_diffuse_var_nir[:] = downward_diffuse_nir
direct_solar_var_nir[:] = direct_solar_nir

vert_path_group = netcdf.createGroup("vert_path")
los_paths_group = netcdf.createGroup("los_paths")

los_paths_group.createDimension("num_los_paths", los_paths.shape[2])

vert_path_names = ["HNO3", "O3 UV", "CNTMSLF1", "CNTMSLF2", "CNTMFRN", "N2 CONT", "MOL SCAT 550NM", "TOTAL AER 550NM", "AER 1 550NM", "AER 2 550NM", "AER 3 550NM", "AER 4 550NM", "CIRRUS 550NM", "WAT DROP", "ICE PART", "MEAN AER RH", "MOL SCAT 459.45NM", "TOT AER 459.45NM", "AER 1 459.45NM", "AER 2 459.45NM", "AER 3 459.45NM", "AER 4 459.45NM", "CIRRUS 459.45NM", "H2O", "O3", "CO2", "CO", "CH4", "N2O", "O2", "NH3", "NO", "NO2", "SO2", "F11", "F12", "CCl3F", "CF4", "F22", "F113", "F114", "R115", "ClONO2", "HNO4", "CHCl2F", "CCl4", "N2O5", "H2-H2", "H2-HE", "H2-CH4", "CH4-CH4"]
vert_path_units = ["ATM CM", "ATM CM", "MOL CM-2", "MOL CM-2", "MOL CM-2", "ATM^2 KM", "550 NM EXTINCTION", "550 NM EXTINCTION", "550 NM EXTINCTION", "550 NM EXTINCTION", "550 NM EXTINCTION", "550 NM EXTINCTION", "550 NM EXTINCTION", "KM GM/M3", "KM GM/M3", "PERCENT", "459.45 NM EXTINCTION", "459.45 NM EXTINCTION", "459.45 NM EXTINCTION", "459.45 NM EXTINCTION", "459.45 NM EXTINCTION", "459.45 NM EXTINCTION", "459.45 NM EXTINCTION", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM^2 CM", "ATM^2 CM", "ATM^2 CM", "ATM^2 CM"]
for n, u, i in zip(vert_path_names, vert_path_units, range(len(vert_path_names))):
        vert_path_var = vert_path_group.createVariable(n, "f8", ("num_cases", "num_szas"))
        vert_path_var[:] = vert_path[:, :, i]
        vert_path_var.units = u

los_paths_names = ["HNO3", "O3 UV", "CNTMSLF1", "CNTMSLF2", "CNTMFRN", "N2 CONT", "MOL SCAT 550NM", "TOTAL AER 550NM", "AER 1 550NM", "AER 2 550NM", "AER 3 550NM", "AER 4 550NM", "CIRRUS 550NM", "WAT DROP", "ICE PART", "MEAN AER RH", "MOL SCAT 459.45NM", "TOT AER 459.45NM", "AER 1 459.45NM", "AER 2 459.45NM", "AER 3 459.45NM", "AER 4 459.45NM", "CIRRUS 459.45NM", "H2O", "O3", "CO2", "CO", "CH4", "N2O", "O2", "NH3", "NO", "NO2", "SO2", "F11", "F12", "CCl3F", "CF4", "F22", "F113", "F114", "R115", "ClONO2", "HNO4", "CHCl2F", "CCl4", "N2O5", "H2-H2", "H2-HE", "H2-CH4", "CH4-CH4"]
los_paths_units = ["ATM CM", "ATM CM", "MOL CM-2", "MOL CM-2", "MOL CM-2", "ATM^2 KM", "550 NM EXTINCTION", "550 NM EXTINCTION", "550 NM EXTINCTION", "550 NM EXTINCTION", "550 NM EXTINCTION", "550 NM EXTINCTION", "550 NM EXTINCTION", "KM GM/M3", "KM GM/M3", "PERCENT", "459.45 NM EXTINCTION", "459.45 NM EXTINCTION", "459.45 NM EXTINCTION", "459.45 NM EXTINCTION", "459.45 NM EXTINCTION", "459.45 NM EXTINCTION", "459.45 NM EXTINCTION", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM CM", "ATM^2 CM", "ATM^2 CM", "ATM^2 CM", "ATM^2 CM"]
for n, u, i in zip(los_paths_names, los_paths_units, range(len(los_paths_names))):
        los_paths_var = los_paths_group.createVariable(n, "f8", ("num_cases", "num_szas", "num_los_paths"))
        los_paths_var[:] = los_paths[:, :, :, i]
        los_paths_var.units = u

geometries_names = ["SCATTER_ALTITUDE", "SUBTENDED_ANGLE", "SOLAR_ZENITH", "PATH_ZENITH", "RELATIVE_AZIMUTH", "SCATTER_ANGLE"]
geometries_units = ["KM", "DEG", "DEG", "DEG", "DEG", "DEG"]
for n, u, i in zip(geometries_names, geometries_units, range(len(geometries_names))):
	geometries_vert_path_var = vert_path_group.createVariable(n, "f8", ("num_cases", "num_szas", "num_path_geometries"))
	geometries_vert_path_var[:] = geometries_vert_path[:, :, :, i]
	geometries_vert_path_var.units = u

	geometries_los_paths_var = los_paths_group.createVariable(n, "f8", ("num_cases", "num_szas", "num_los_paths", "num_path_geometries"))
	geometries_los_paths_var[:] = geometries_los_paths[:, :, :, :, i]
	geometries_los_paths_var.units = u

netcdf.close()
