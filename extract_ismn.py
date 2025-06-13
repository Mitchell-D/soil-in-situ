""" """
from pathlib import Path
import numpy as np
import json
from datetime import datetime,timedelta
import pickle as pkl

from ismn.interface import ISMN_Interface as ISMN

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

    return [hdict,values,datetimes] + [[],[etimes]][return_epoch_times]

if __name__=="__main__":
    #proj_root_dir = Path("/Users/mtdodson/desktop/soilm-in-situ")
    proj_root_dir = Path("/rhome/mdodson/soil-in-situ")
    proj_data_dir = Path("/rstor/mdodson/in-situ/ismn")
    ismn_stations_path = proj_data_dir.joinpath("station-data")
    ismn_available_json = proj_root_dir.joinpath("data/ismn-datamenu.json")

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

    ismn_dm = json.load(ismn_available_json.open("r"))
    variables = set()
    stations = []
    for ntw in ismn_dm.keys():
        for stn in ismn_dm[ntw].keys():
            stations.append({
                "network":ntw,
                "station":stn,
                "sensors":ismn_dm[ntw][stn]["sensors"],
                "location":ismn_dm[ntw][stn]["location"],
                "station_meta":ismn_dm[ntw][stn]["station_meta"]
                })

    for stn in stations:
        data = {}
        for ssr in stn["sensors"]:
            cur_var = var_mapping[ssr["variable"]]
            if cur_var not in data.keys():
                data[cur_var] = []
            data[cur_var].append({
                **ssr,
                "data":parse_ismn_stm(
                    stm_path=ismn_stations_path.joinpath(ssr["file"]),
                    return_epoch_times=True,
                    debug=False,
                    )
                })
        print()
        print(stn["network"], stn["station"])
        unq_times = set()
        for vk in data.keys():
            for ssr in data[vk]:
                _,vals,datetimes,etimes = ssr["data"]
                unq_times.update(set(datetimes))
            ## strangely common to have skipped and doubled hourly timesteps
            ## separated by 3205 hours exactly
            '''
            dt = np.diff(etimes)
            m_dt_under = dt<3599
            m_dt_over = dt>3601
            print(k, vals.size, np.count_nonzero(m_dt_under),
                    np.count_nonzero(m_dt_over))
            print(np.where(m_dt_under)[0])
            print(np.where(m_dt_over)[0])
            if np.count_nonzero(m_dt_under) == np.count_nonzero(m_dt_over)+1:
                print(np.where(m_dt_over)[0]-np.where(m_dt_under)[0][1:])
            if np.count_nonzero(m_dt_under) == np.count_nonzero(m_dt_over):
                print(np.where(m_dt_over)[0]-np.where(m_dt_under)[0])
            print()
            '''

        unq_times = sorted(list(unq_times))
        print(len(unq_times), unq_times[0], unq_times[-1])
        for t0,tf in zip(unq_times[:-1],unq_times[1:]):
            if abs((tf-t0).seconds)-3600 > 60:
                print(t0)
