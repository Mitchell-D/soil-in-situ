""" """
from pathlib import Path
import numpy as np
import json
from datetime import datetime,timedelta
import pickle as pkl
from multiprocessing import Pool

from ismn.interface import ISMN_Interface as ISMN

def _mp_preproc_station_data(args):
    return _preprocess_station_data(**args)

def _preprocess_station_data(station_dict:dict, var_mapping:dict,
        station_pkl_dir:Path):
    """
    for each station, make a dict containing all information for each sensor,
    including parsed value and flag data, and store it in a pkl file
    """
    stn = station_dict
    ssr_dict = {}
    for ssr in stn["sensors"]:
        cur_var = var_mapping[ssr["variable"]]
        if cur_var not in ssr_dict.keys():
            ssr_dict[cur_var] = []
        ## combine sensor meta dict and parsed sensor data
        ssr_dict[cur_var].append({
            **ssr,
            **{k:v for k,v in zip(
                ["header", "values", "flags", "datetimes", "etimes"],
                parse_ismn_stm(
                    stm_path=ismn_stations_path.joinpath(ssr["file"]),
                    return_epoch_times=True,
                    debug=False,
                    ))
                },
            })
    ## Determine the universal minimum and maximum time for all sensors
    print(stn["network"], stn["station"])
    min_times, max_times = [],[]
    for vk in ssr_dict.keys():
        for ssr in ssr_dict[vk]:
            min_times.append(min(ssr["datetimes"]))
            max_times.append(max(ssr["datetimes"]))
    tmin,tmax = min(min_times),max(max_times)

    ## Create a dict mapping each hourly time to a unique index
    t,ix_fulltime = tmin,0
    all_times = {}
    while t <= tmax:
        all_times[t.strftime("%Y%m%d%H")] = ix_fulltime
        ix_fulltime += 1
        t += timedelta(hours=1)

    ## make a dict of consistent-time arrays of data and flags per sensor
    duplicates = {}
    data_dict = {}
    for vk in ssr_dict.keys():
        data_dict[vk] = []
        for ix_ssr,ssr in enumerate(ssr_dict[vk]):
            #hdict,vals,datetimes,etimes = ssr["data"]
            vals = ssr["values"]
            datetimes = ssr["datetimes"]
            ismn_flags,src_flags = ssr["flags"]
            tmp_array = np.full(ix_fulltime+1, np.nan)
            tmp_iflags = np.full(ix_fulltime+1, "-")
            if src_flags is None:
                tmp_sflags = None
            else:
                tmp_sflags = np.full(ix_fulltime+1, "-")
            for i,(v,t) in enumerate(zip(vals,datetimes)):
                tk = t.strftime("%Y%m%d%H")
                if np.isfinite(tmp_array[all_times[tk]]):
                    if vk not in duplicates.keys():
                        duplicates[vk] = {}
                    duplicates[vk][ix_ssr] = all_times[tk]
                tmp_array[all_times[tk]] = v
                tmp_iflags[all_times[tk]] = ismn_flags[i]
                if not src_flags is None:
                    tmp_sflags[all_times[tk]] = src_flags[i]
            data_dict[vk].append((tmp_array, (ismn_flags,src_flags)))

    ## collect all depths and sensors per depth into a single dict
    adata = []
    alabels = []
    amasks = []
    adepths = []
    aduplicates = []
    sensors = {}
    times = list(sorted(all_times.keys()))
    for vk in ssr_dict.keys():
        ssr_combos = []
        valid_depths = set(tuple(ssr["depth"]) for ssr in ssr_dict[vk])
        depths_dict = {d:[
            (six,ssr) for six,ssr in enumerate(ssr_dict[vk])
            if tuple(ssr["depth"])==d
            ] for d in valid_depths}
        ## sort by depths
        for dix,(depth,depth_members) in enumerate(sorted(
                depths_dict.items(), key=lambda s:s[0])):
            ## sort by instrument name per depth
            depth_members = list(sorted(
                depth_members, key=lambda d:d[1]["instrument"]
                ))
            adepths.append(depth)
            for nix,(six,ssr) in enumerate(depth_members):
                skey = f"{vk}_{dix:02}_{nix:02}"
                alabels.append(skey)
                adata.append(data_dict[vk][six][0])
                amasks.append(data_dict[vk][six][1])
                sensors[skey] = stn["sensors"][six]
                if vk in duplicates.keys():
                    if six in duplicates[vk].keys():
                        aduplicates.append(duplicates[vk][six])

    network_name = stn["network"].replace(".","")
    station_name = stn["station"].replace(".","")
    station_pkl_dict = {
            ## string network and station ID
            "network":network_name,
            "station":station_name,
            ## dict of file, depth, etc info on each sensor
            "sensors":sensors,
            ## soil, climate, etc meta-info on the station
            "station_meta":stn["station_meta"],
            ## (lat, lon, elevation) of the station
            "location":stn["location"],
            ## all valid depths ordered by the sensor labels' second field
            "depths":adepths,
            ## YYYYMMDDHH times associated with the first axis of data
            "times":times,
            ## all ordered according to the sensor labels labels
            "labels":alabels,
            "data":np.stack(adata, axis=-1),
            "masks":amasks,
            "duplicate_times":aduplicates,
            }
    print(network_name, station_name, alabels)
    pkl_path = station_pkl_dir.joinpath(
            f"station_{network_name}_{station_name}.pkl")
    pkl.dump(station_pkl_dict, pkl_path.open("wb"))
    return pkl_path

def parse_ismn_stm(stm_path:Path, return_epoch_times=False, debug=False):
    """
    Extract instrument file time series from files following the ISMN format

    :@param stm_path: Path to a sensor stm file from ismn.bafg.de
    :@param return_epoch_times: If True, returns an array of float epoch times
        rather than a list of datetime objects corresponding to each timestep
    :@param debug: If True, return information on sample count and time
        interval consistency for the file.

    :@return: 3-tuple (header_dict, data_values, times) Where header_dict
        contains the CSE, network, and station id, station location, and
        sensor
    """
    stm_path = Path(stm_path)
    assert stm_path.exists()

    ## header labels
    hlabels = ["cseid", "network", "station", "lat", "lon", "elev", "d0", "df"]
    ## header fields that should be converted to float
    ffields = ["lat", "lon", "elev", "d0", "df"]
    ## parse the header as a dict
    lines = stm_path.open("r").readlines()
    hline = [v for v in lines[0].replace("\n","").split(" ") if v]
    sensor_info = " ".join(hline[len(hlabels):])
    hline = hline[:len(hlabels)]
    hdict = {l:v if l not in ffields else float(v)
            for l,v in zip(hlabels,hline)}
    hdict["sensor_info"] = sensor_info

    ## read the lines, skipping the header
    date,time,values,*flags = zip(*map(
        lambda l:l.replace("\n","").split(" "), lines[1:]
        ))

    ## separate the flag columns if there are more than 1
    if len(flags) == 2:
        flags_ismn,flags_provider = flags
    elif len(flags) == 1:
        flags_ismn = flags
        flags_provider = None
    else:
        raise ValueError(f"Why are there more than 3 flag columns? >:(", flags)

    ## extract datetimes and convert to epoch times
    datetimes = [datetime.strptime(" ".join(tt), "%Y/%m/%d %H:%M")
            for tt in zip(date,time) ]
    etimes = np.asarray([float(t.strftime("%s")) for t in datetimes])

    ## if debugging print sample size and time interval consistency
    if debug:
        dt = np.diff(etimes)
        print(f"{etimes.size = } {np.amin(dt)}, {np.amax(dt)}, {np.average(dt)}, {np.count_nonzero(dt < 3599)}, {np.count_nonzero(dt > 3601)}")

    ## convert values to a float array
    values = np.asarray([float(v) for v in values])

    return [hdict,values,(flags_ismn,flags_provider),datetimes] + \
            [[],[etimes]][return_epoch_times]

if __name__=="__main__":
    #proj_root_dir = Path("/Users/mtdodson/desktop/soilm-in-situ")
    proj_root_dir = Path("/rhome/mdodson/soil-in-situ")
    proj_data_dir = Path("/rstor/mdodson/in-situ/ismn")
    ismn_stations_path = proj_data_dir.joinpath("station-data")
    ismn_available_json = proj_root_dir.joinpath("data/ismn-datamenu.json")
    station_pkl_dir = proj_data_dir.joinpath("station-pkls")

    ## subset of sensor networks to utilize
    extract_networks = ["TxSON", "SOILSCAPE", "SNOTEL", "SCAN", "RISMA",
            "FLUXNET-AMERIFLUX", "CW3E", "ARM"]
    ## mapping to variable shorthand
    var_mapping = {
            "snow_depth":"snod", "surface_temperature":"tsfc",
            "precipitation":"prcp", "soil_temperature":"tsoil",
            "air_temperature":"tair", "soil_moisture":"soilm",
            "snow_water_equivalent":"swe",
            }
    keep_meta = ["clay_fraction", "sand_fraction", "silt_fraction",
            "saturation", "climate_KG", "climate_insitu", "elevation",
            "instrument", "organic_carbon"]

    ## Load the metadata via the interface tool
    '''
    interface = ISMN(
            ismn_stations_path,
            parallel=True,
            meta_path=proj_root_dir.joinpath("data/ismn_meta"),
            temp_root=proj_root_dir.joinpath("ismn_tmp")
            )

    ismn_dm = {}
    for network in interface:
        ismn_dm[network.name] = {}
        for station in network:
            ismn_dm[network.name][station.name] = {
                    "location":(station.lat, station.lon, station.elev),
                    "station_meta":{
                        k:v for k,v in station.metadata.to_dict().items()
                        if k in keep_meta
                        },
                    "time_range":[None if t is None else int(t.timestamp())
                        for t in station.get_min_max_obs_timestamp()],
                    "sensors":[
                        {"name":s.name, "variable":s.variable,
                            "depth":(s.depth.start, s.depth.end),
                            "instrument":s.instrument,
                            "file":s.filehandler.posix_path.as_posix(),
                            }
                        for s in station.sensors.values()
                        ]
                    }
    json.dump(ismn_dm, ismn_available_json.open("w"), indent=2)
    '''

    ## fix each station's data to a consistent time range/interval and store
    ## the data in a pkl file alongside its metadata
    #'''
    ## Collect a list of station dicts alsongside network and station info
    ismn_dm = json.load(ismn_available_json.open("r"))
    stations = []
    for ntw in ismn_dm.keys():
        if ntw not in extract_networks:
            continue
        for stn in ismn_dm[ntw].keys():
            stations.append({
                "network":ntw,
                "station":stn,
                "sensors":ismn_dm[ntw][stn]["sensors"],
                "location":ismn_dm[ntw][stn]["location"],
                "station_meta":ismn_dm[ntw][stn]["station_meta"]
                })

    args = [{
        "station_dict":s,
        "var_mapping":var_mapping,
        "station_pkl_dir":station_pkl_dir,
        } for s in stations
        ]
    nworkers = 15
    with Pool(nworkers) as pool:
        for ppath in pool.imap_unordered(_mp_preproc_station_data, args):
            print(f"generated {ppath.name}")
