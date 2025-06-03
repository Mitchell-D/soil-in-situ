
import requests
from pathlib import Path
import numpy as np
import json
from datetime import datetime
import pickle as pkl

base_url = "https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1"

def get_data_menu():
    r = requests.get(base_url + "/reference-data")
    try:
        assert r.status_code == 200
    except:
        print(f"{r.status_code = }")
        print(r)
        raise ValueError(f"Unable to complete request")
    rj = r.json()
    return {
            "locales":{d["code"]:d["name"] for d in rj["dcos"]},
            "feats":{d["code"]:d for d in rj["elements"]},
            "periods":{d["code"]:d for d in rj["forecastPeriods"]},
            #"instruments":rj["instruments"],
            "networks":{d["code"]:d for d in rj["networks"]},
            "states":{d["code"]:d for d in rj["states"]},
            "units":{d["code"]:d for d in rj["units"]},
            }

def get_stations_menu():
    r = requests.get(base_url + "/stations")
    try:
        assert r.status_code == 200
    except:
        print(f"{r.status_code = }")
        print(r)
        raise ValueError(f"Unable to complete request")
    rj = r.json()
    return {d["stationTriplet"]:d for d in rj}

def get_station_data(
        station_triplet:str, features:list,  begin_date:datetime,
        end_date:datetime, duration="DAILY", return_flags=True,
        **other_params,
        ):
    """
    api reference: https://wcc.sc.egov.usda.gov/awdbRestApi/
    """
    query = base_url + "/data?"
    params = {
            "stationTriplets":station_triplet,
            "elements":",".join([f"{f}:*:*" for f in features]),
            "beginDate":begin_date.strftime("%Y-%m-%d %H:%M"),
            "endDate":end_date.strftime("%Y-%m-%d %H:%M"),
            "duration":duration.upper(),
            "returnFlags":["false","true"][return_flags],
            **other_params,
            }
    query += "&".join([f"{k}={v}" for k,v in params.items()])
    print(query)
    r = requests.get(query)
    try:
        assert r.status_code == 200
    except:
        print(f"{r.status_code = }")
        print(r)
        raise ValueError(f"Unable to complete request")
    if len(r.json()) == 0:
        raise ValueError(f"No data found: {station_triplet}")
    if len(r.json()) > 1:
        raise ValueError(f"Multiple stations found: {station_triplet}")
    station = r.json()[0]
    sdata = {}
    for d in station["data"]:
        se = d["stationElement"]
        sdkey = ":".join(map(str,[
            se["elementCode"], se.get("heightDepth",""), se.get("ordinal", "")
            ]))
        etimes,dvals,m_qc = [],[],[]
        try:
            for v in d["values"]:
                t = datetime.strptime(v["date"], "%Y-%m-%d %H:%M")
                etimes.append(int(t.strftime("%s")))
                dvals.append(float(v["value"]))
                m_qc.append(v["qcFlag"])
        except:
            print(v)
            raise ValueError(f"Error extracting data: {station_triplet}")
        sdata[sdkey] = {
                "finfo":d["stationElement"],
                "data":np.asarray(dvals),
                "etimes":etimes,
                "flags":np.asarray(m_qc),
                }
    return sdata

if __name__=="__main__":
    proj_root_dir = Path("/Users/mtdodson/desktop/soilm-in-situ")
    station_json_path = proj_root_dir.joinpath("scan-stations.json")
    datamenu_json_path = proj_root_dir.joinpath("scan-datamenu.json")
    pkl_dir = proj_root_dir.joinpath("station-pkls")

    ## determine which features to extract if available
    extract_feats = [
            "TAVG", "PRES", "DPTP", "PREC", "PRCP", "LWINV", "SWINV", "WTEQ",
            "WDIR", "WSPDV", "SMS", "SMV", "SMX", "SMN",

            #"TMAX", "TMIN", "HFTV", "EVAP", "TGSV",
            #"NTRDV", "NTRDX", "NTRDN", "PTEMP",
            #"LWOTV", "SWOTV", "REST", "RESC", "SRDOO", "SNDN", "SNWD", "SNWDV",
            #"SNWDX", "SNWDN", "SNRR", "STV", "STX", "STN", "SRAD", "SRADV",
            #"SRADX", "SRADN", "SRMV", "SRMX", "SRMN", "WSPDX", "WSPDN",
            ]
    begin_date = datetime(2018, 1, 1)
    end_date = datetime(2024, 1, 1)
    skip_existing = True

    ## get the REST API parameters for stations and data types
    if not station_json_path.exists():
        smenu = get_stations_menu()
        json.dump(smenu, station_json_path.open("w"))
    else:
        smenu = json.loads(station_json_path.open("r").read())
    if not datamenu_json_path.exists():
        dmenu = get_data_menu()
        json.dump(dmenu, datamenu_json_path.open("w"))
    else:
        dmenu = json.load(datamenu_json_path.open("r"))

    ## print full names of selected features that are available
    for k in extract_feats:
        print(f"{k:12} : {dmenu['feats'][k]['name']}")

    ## filter for
    scan = []
    for s in smenu.keys():
        if smenu[s]["networkCode"] == "SCAN":
            scan.append(s)
    smenu = {k:v for k,v in smenu.items() if smenu[k]["networkCode"]=="SCAN"}
    for s in smenu.keys():
        pkl_path = pkl_dir.joinpath(f"scandata_{s.replace(':','-')}.pkl")
        if pkl_path.exists() and skip_existing:
            print(f"Skipping existing file: {pkl_path.as_posix()}")
            continue
        try:
            sdata = get_station_data(
                    station_triplet=s,
                    features=extract_feats,
                    begin_date=begin_date,
                    end_date=end_date,
                    duration="HOURLY",
                    return_flags=True,
                    )
        except ValueError as e:
            print(e)
            continue
        pkl.dump(sdata, pkl_path.open("wb"))
        print(f"Generated {pkl_path.as_posix()}")
