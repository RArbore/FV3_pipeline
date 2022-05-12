from argparse import ArgumentParser
import scipy.interpolate
import multiprocessing
import numpy as np
import json
import os

argparser = ArgumentParser(description="Add structure to JSON files for MODTRAN input (specific to FV3 OSSE inputs).")
argparser.add_argument("input", help="Input file path.")
argparser.add_argument("input_sap", help="Input SAP file path.")
argparser.add_argument("output", help="Output file path.")
args = argparser.parse_args()

case_num = 0

MODTRAN_MINIMUM_ZENITH = 115.8

def write_modtran6_json_file(fname,
                             output_fname,
                             observer_azimuth_angles,
                             observer_zenith_angles,
                             sza,
                             case_name,
                             case_description,
                             nstr,
                             alts,
                             press,
                             temps,
                             h2o_vals,
                             co2_vals,
                             o3_vals,
                             n2o_vals,
                             co_vals,
                             ch4_vals,
                             cloud_alts,
                             cliq_prof,
                             cice_prof,
                             sap_file_path):
    global case_num
    out_json = [
                  { "MODTRANINPUT" :
                     { "NAME" : case_name,
                        "DESCRIPTION" : case_description,
                        "CASE" : case_num,
                        "RTOPTIONS" :
                        {
                            "MODTRN" : "RT_MODTRAN",
                            "LYMOLC" : False,
                            "T_BEST" : False,
                            "IEMSCT" : "RT_SOLAR_AND_THERMAL",
                            "IMULT" : "RT_DISORT", # or RT_DISORT_AT_OBS
                            "DISALB" : True,
                            "NSTR" : nstr,
                            "SOLCON" : 0.0
                        },
                        "ATMOSPHERE" : {
                            "MODEL" : "ATM_USER_ALT_PROFILE",
                            "MDEF" : 1,
                            "CO2MX" : 330.0,
                            "HMODEL" : "New Atm Profile",
                            "NPROF" : 9, # possibly change this
                            "NLAYERS" : alts.shape[0],
                            "PROFILES" : [
                                {
                                    "TYPE" : "PROF_ALTITUDE",
                                    "UNITS" : "UNT_KILOMETERS",
                                    "PROFILE" : list(map(lambda x: round(x, 1), alts.tolist()))
                                },
                                {
                                    "TYPE" : "PROF_PRESSURE",
                                    "UNITS" : "UNT_PMILLIBAR",
                                    "PROFILE" : list(map(lambda x: round(x, 5), press.tolist()))
                                },
                                {
                                    "TYPE" : "PROF_TEMPERATURE",
                                    "UNITS" : "UNT_TKELVIN",
                                    "PROFILE" : list(map(lambda x: round(x, 5), temps.tolist()))
                                },
                                {
                                    "TYPE" : "PROF_H2O",
                                    "UNITS" : "UNT_DPPMV",
                                    "PROFILE" : list(map(lambda x: round(x, 5), h2o_vals.tolist()))
                                },
                                {
                                    "TYPE" : "PROF_CO2",
                                    "UNITS" : "UNT_DPPMV",
                                    "PROFILE" : list(map(lambda x: round(x, 5), co2_vals.tolist()))
                                },
                                {
                                    "TYPE" : "PROF_O3",
                                    "UNITS" : "UNT_DPPMV",
                                    "PROFILE" : list(map(lambda x: round(x, 5), o3_vals.tolist()))
                                },
                                {
                                    "TYPE" : "PROF_N2O",
                                    "UNITS" : "UNT_DPPMV",
                                    "PROFILE" : list(map(lambda x: round(x, 5), n2o_vals.tolist()))
                                },
                                {
                                    "TYPE" : "PROF_CO",
                                    "UNITS" : "UNT_DPPMV",
                                    "PROFILE" : list(map(lambda x: round(x, 5), co_vals.tolist()))
                                },
                                {
                                    "TYPE" : "PROF_CH4",
                                    "UNITS" : "UNT_DPPMV",
                                    "PROFILE" : list(map(lambda x: round(x, 5), ch4_vals.tolist()))
                                }
                            ],
                        },
                        "AEROSOLS" : {
                            "IHAZE" : "AER_RURAL",
                            "ICLD" : "CLOUD_CUMULUS",
                            "VIS" : 5.0, # change this
                            "ARUSS": "SAP",
                            "SAPFILE": sap_file_path,
                            "CLDALT" : {
                                "NCRALT" : cloud_alts.shape[0],
                                "ZPCLD" : list(map(lambda x: round(x, 9), cloud_alts.tolist())),
                                "CLD" : list(map(lambda x: round(x, 9), cliq_prof.tolist())),
                                "CLDICE" : list(map(lambda x: round(x, 9), cice_prof.tolist())),
                                "RR" : [0.0] * cloud_alts.shape[0]
                            }
                        },
                        "GEOMETRY" : {
                            "IDAY" : 16,
                            "IPARM" : 12,
                            "PARM1" : 180.0,
                            "PARM2" : round(sza, 1),
                            "GMTIME" : 15.5,
                            "TRUEAZ" : 270.0,
                            "NLOS" : observer_azimuth_angles.shape[0],
                            "MLOS" : list(map(lambda angles: {
                                                         "H1ALT" : 705.0,
                                                         "OBSZEN" : round(angles[0], 1),
                                                         "AZ_INP" : round(angles[1], 1),
                                                         "LENN" : 1}, zip(observer_zenith_angles.tolist(), observer_azimuth_angles.tolist())))
                        },
                        "SURFACE" : {
                            "SURFTYPE" : "REFL_LAMBER_MODEL",
                            "SURREF" : 0.4,
                            "NSURF" : 2,
                            "SURFA" : { "CSALB" : "LAMB_OCEAN_WATER" },
                            "SURFNLOS" : observer_azimuth_angles.shape[0],
                            "SURFLOS" : [{ "CSALB" : "LAMB_OCEAN_WATER" }] * (observer_azimuth_angles.shape[0])
                        },
                        "SPECTRAL" : {
                            "V1" : 300.0,
                            "V2" : 5000.0,
                            "DV" : 5.0,
                            "FWHM" : 40.0,
                            "YFLAG" : "R",
                            "XFLAG" : "N",
                            "FLAGS" : "NTA   T",
                            "LBMNAM" : "T",
                            "BMNAME" : "15_2013"
                        },
                        "FILEOPTIONS" : {
                            "MSGPRNT" : "MSG_DEBUG",
                            "JSONOPT" : "WRT_ALL",
                            "NOPRNT" : 1,
                            "JSONPRNT" : output_fname + ".json"
                        }
                     }
                  }
               ]
    case_num += 1
    return out_json

observer_azimuth_angles = 5 * np.array(range(73), dtype=float)
observer_zenith_angles = 120 * np.ones(observer_azimuth_angles.shape, dtype=float)

alts = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0])
press = np.array([1013.25165, 898.7637, 795.0062, 701.212, 616.59973, 540.48126, 472.17044, 411.04807, 356.50916, 308.00012, 264.98944, 226.99094, 193.99063, 165.7897, 141.6997, 121.11009, 103.5195, 88.49714, 75.6517, 64.67406, 55.293205, 47.289112, 40.47482, 34.668102, 29.717024, 25.491879, 17.42908, 11.969957, 8.257732, 5.7459197, 4.041395, 2.8714008, 2.0596993, 1.4910045, 1.0885973, 0.79778993, 0.42524916, 0.21957971, 0.109290056, 0.052209016, 0.023881052, 0.010523967, 0.0044568186, 0.0018359008, 7.596596E-4, 3.201086E-4, 1.4477066E-4, 7.1041926E-5, 4.0096143E-5, 2.538202E-5])
temps = np.array([288.15, 281.651, 275.154, 268.659, 262.166, 255.676, 249.187, 242.7, 236.215, 229.733, 223.252, 216.774, 216.65, 216.65, 216.65, 216.65, 216.65, 216.65, 216.65, 216.65, 216.65, 217.581, 218.574, 219.567, 220.56, 221.552, 224.032, 226.509, 229.587, 236.513, 243.434, 250.35, 257.26, 264.164, 270.65, 270.65, 260.771, 247.021, 233.292, 219.585, 208.399, 198.639, 188.893, 186.87, 188.42, 195.08, 208.84, 240.0, 300.0, 360.0])
h2o_vals = np.array([7745.0, 6071.0, 4631.0, 3182.0, 2158.0, 1397.0, 925.4, 572.0, 366.7, 158.3, 69.96, 36.13, 19.06, 10.85, 5.927, 5.0, 3.95, 3.85, 3.825, 3.85, 3.9, 3.975, 4.065, 4.2, 4.3, 4.425, 4.575, 4.725, 4.825, 4.9, 4.95, 5.025, 5.15, 5.225, 5.25, 5.225, 5.1, 4.75, 4.2, 3.5, 2.825, 2.05, 1.33, 0.85, 0.54, 0.4, 0.34, 0.28, 0.24, 0.2])
co2_vals = np.array([330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 328.0, 320.0, 310.0, 270.0, 195.0, 110.0, 60.0, 40.0, 35.0])
o3_vals = np.array([0.0266, 0.02931, 0.03237, 0.03318, 0.03387, 0.03768, 0.04112, 0.05009, 0.05966, 0.09168, 0.1313, 0.2149, 0.3095, 0.3846, 0.503, 0.6505, 0.8701, 1.187, 1.587, 2.03, 2.579, 3.028, 3.647, 4.168, 4.627, 5.118, 5.803, 6.553, 7.373, 7.837, 7.8, 7.3, 6.2, 5.25, 4.1, 3.1, 1.8, 1.1, 0.7, 0.3, 0.25, 0.3, 0.5, 0.7, 0.7, 0.4, 0.2, 0.05, 0.005, 5.0E-4])
n2o_vals = np.array([0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.3195, 0.3179, 0.314, 0.3095, 0.3048, 0.2999, 0.2944, 0.2877, 0.2783, 0.2671, 0.2527, 0.2365, 0.2194, 0.2051, 0.1967, 0.1875, 0.1756, 0.1588, 0.1416, 0.1165, 0.09275, 0.06693, 0.04513, 0.02751, 0.01591, 0.009378, 0.004752, 0.003, 0.002065, 0.001507, 0.001149, 8.89E-4, 7.056E-4, 5.716E-4, 4.708E-4, 3.932E-4, 3.323E-4, 2.837E-4, 2.443E-4, 2.12E-4, 1.851E-4])
co_vals = np.array([0.15, 0.145, 0.1399, 0.1349, 0.1312, 0.1303, 0.1288, 0.1247, 0.1185, 0.1094, 0.09962, 0.08964, 0.07814, 0.06374, 0.05025, 0.03941, 0.03069, 0.02489, 0.01966, 0.01549, 0.01331, 0.01232, 0.01232, 0.01307, 0.014, 0.01498, 0.01598, 0.0171, 0.0185, 0.02009, 0.0222, 0.02497, 0.02824, 0.03241, 0.03717, 0.04597, 0.06639, 0.1073, 0.1862, 0.3059, 0.6375, 1.497, 3.239, 5.843, 10.13, 16.92, 24.67, 33.56, 41.48, 50.0])
ch4_vals = np.array([1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.699, 1.697, 1.693, 1.685, 1.675, 1.662, 1.645, 1.626, 1.605, 1.582, 1.553, 1.521, 1.48, 1.424, 1.355, 1.272, 1.191, 1.118, 1.055, 0.987, 0.9136, 0.83, 0.746, 0.6618, 0.5638, 0.4614, 0.3631, 0.2773, 0.21, 0.165, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.14, 0.13, 0.12, 0.11, 0.095, 0.06, 0.03])

#cloud_alts = np.array([0.0, 0.09, 0.11, 2.49, 2.51])
#cliqwp_vals = np.array([0.0, 0.0, 0.68, 0.68, 0.0])
#cicewp_vals = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

input_file_path = args.input
with open(input_file_path, "r") as input_file:
    json_data = json.loads(input_file.read())

lat = np.array(json_data["lat"])
lon = np.array(json_data["lon"])
q_plev = np.array(json_data["q_plev"])
t_plev = np.array(json_data["t_plev"])
qi_plev = np.array(json_data["qi_plev"])
ql_plev = np.array(json_data["ql_plev"])
ts = np.array(json_data["ts"])
ps = np.array(json_data["ps"])
plev = np.array(json_data["plev"])
cldc = np.array(json_data["cldc"])
SZA = np.array(json_data["SZA"])
VZA = np.array(json_data["VZA"])
RAA = np.array(json_data["RAA"])

json_data = { "MODTRAN" : [] }

for i in range(lat.shape[0]):
    for j in range(SZA.shape[0]):
        #print("Starting case " + str(i) + " w/ SZA " + str(j))
        out_file_path = "FV3_modtran6_json_mine/Libera_ADM_FV3_Input_Case_{:02d}_SZA_{:02d}.json".format(i+1, j+1)

        case_name = "Libera_ADM_FV3_Input_Case_{:02d}_SZA_{:02d}".format(i+1, j+1)
        case_description = "Libera ADM Calculation from FV3 output for Profile {} and SZA = {}. Created by Daniel Feldman for Libera Science Team in cooperation with Jake Gristey, Sebastian Schmidt, Hong Chen, Maria Hakuba, Bruce Kindel, and Peter Pilewskie".format(i+1, round(SZA[j], 2))

        press_fv3 = plev[::-1]
        alts_fv3 = -7 * np.log(press_fv3 / 1000)
        temps_fv3 = t_plev[:, i][::-1]
        h2o_vals_fv3 = q_plev[:, i][::-1] * 28 / 18 * 1E6
        co2_vals_fv3 = 430 * np.ones(alts_fv3.shape)

        n2o_vals_fv3 = scipy.interpolate.interp1d(press, n2o_vals, kind="linear")(press_fv3)
        o3_vals_fv3 = scipy.interpolate.interp1d(press, o3_vals, kind="linear")(press_fv3)
        co_vals_fv3 = scipy.interpolate.interp1d(press, co_vals, kind="linear")(press_fv3)
        ch4_vals_fv3 = scipy.interpolate.interp1d(press, ch4_vals, kind="linear")(press_fv3)

        fv3_observer_zenith_angles = np.zeros((VZA.shape[0] * RAA.shape[0]))
        fv3_observer_azimuth_angles = np.zeros((VZA.shape[0] * RAA.shape[0]))

        for k in range(RAA.shape[0]):
            fv3_observer_zenith_angles[k * VZA.shape[0] : (k + 1) * VZA.shape[0]] = 90 + VZA
            fv3_observer_azimuth_angles[k * VZA.shape[0] : (k + 1) * VZA.shape[0]] = RAA[k]

        cloud_alts = np.array([0.0, 0.09, 0.11, 2.49, 2.51])
        cliqwp_vals = np.array([0.0, 0.0, 0.001, 0.001, 0.0])
        cicewp_vals = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

        thresh_ql_plev = ql_plev[:, i]
        thresh_ql_plev[thresh_ql_plev < 1E-8] = 0.0
        thresh_qi_plev = qi_plev[:, i]
        thresh_qi_plev[thresh_qi_plev < 1E-8] = 0.0

        cloud_alts_fv3 = []
        cliqwp_vals_fv3 = []
        cicewp_vals_fv3 = []

        alt_last = 0.0
        liq_last = 0.0
        ice_last = 0.0
        for (alt, liq, ice) in zip(alts_fv3, thresh_ql_plev, thresh_qi_plev):
            if (liq > 0.0 and liq_last == 0.0) or (ice > 0.0 and ice_last == 0.0):
                cloud_alts_fv3 += [alt - 0.01, alt]
                cliqwp_vals_fv3 += [0.0 if liq_last == 0.0 else liq, liq]
                cicewp_vals_fv3 += [0.0 if ice_last == 0.0 else ice, ice]
                alt_last = alt
                liq_last = liq
                ice_last = ice
            elif (liq == 0.0 and liq_last > 0.0) or (ice == 0.0 and ice_last > 0.0):
                cloud_alts_fv3 += [alt_last + 0.01, alt]
                cliqwp_vals_fv3 += [0.0 if liq == 0.0 else liq_last, liq]
                cicewp_vals_fv3 += [0.0 if ice == 0.0 else ice_last, ice]
                alt_last = alt
                liq_last = liq
                ice_last = ice
            elif not liq == liq_last or not ice == ice_last:
                cloud_alts_fv3 += [alt]
                cliqwp_vals_fv3 += [liq]
                cicewp_vals_fv3 += [ice]
                alt_last = alt
                liq_last = liq
                ice_last = ice

        cloud_alts_fv3 = np.array(cloud_alts_fv3) * 1.0
        cliqwp_vals_fv3 = np.array(cliqwp_vals_fv3) * 1.0
        cicewp_vals_fv3 = np.array(cicewp_vals_fv3) * 1.0

        case_json = write_modtran6_json_file(out_file_path,
                                             case_name,
                                             fv3_observer_azimuth_angles,
                                             np.maximum(fv3_observer_zenith_angles, MODTRAN_MINIMUM_ZENITH),
                                             SZA[j],
                                             case_name,
                                             case_description,
                                             16,
                                             alts_fv3,
                                             press_fv3,
                                             temps_fv3,
                                             h2o_vals_fv3,
                                             co2_vals_fv3,
                                             o3_vals_fv3,
                                             n2o_vals_fv3,
                                             co_vals_fv3,
                                             ch4_vals_fv3,
                                             cloud_alts_fv3,
                                             cliqwp_vals_fv3,
                                             cicewp_vals_fv3,
                                             os.path.abspath(args.input_sap))

        #temp_data = { "MODTRAN" : [] }
        #temp_data["MODTRAN"] = case_json
        #json.dump(temp_data, open(out_file_path, "w"), indent=4)
        json_data["MODTRAN"] += case_json

json.dump(json_data, open(args.output, "w"), indent=4)
