"""
Extract from USCRN fixed-column text files into pkl files per station across
years, identify contiguous data sequences having all valid needed data, and
extract all such data from all stations into a combined pkl file
"""
import requests
from pathlib import Path
import numpy as np
from pprint import pprint
import json
from datetime import datetime
import pickle as pkl

base_url = "https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1"

fields = [
        ("wbanno", slice(0,5)), ## station number
        ("utc-datetime", slice(6,19)), ## "YYYYmmdd HHMM"
        ("lon", slice(41,48)),
        ("lat", slice(49,56)),
        ("temp-final", slice(57,64)), ## C
        ("temp-mean", slice(65,72)), ## C
        ("precip", slice(89,96)), ## mm/hr
        ("dswrf", slice(97,103)), ## w/m^2
        ("qf-dswrf", slice(104,105)),
        ("sfctemp", slice(126,133)), ## C
        ("qt-sfctemp", slice(124,125)), ## type
        ("qf-sfctemp", slice(134,135)), ## flag
        ("rh", slice(156,161)), ## %
        ("qf-rh", slice(162,163)), ## flag
        ("vsm-5", slice(164,171)), ## m^3/m^3 @ 5cm
        ("vsm-10", slice(172,179)), ## m^3/m^3 @ 10cm
        ("vsm-20", slice(180,187)), ## m^3/m^3 @ 20cm
        ("vsm-50", slice(188,195)), ## m^3/m^3 @ 50cm
        ("vsm-100", slice(196,203)), ## m^3/m^3 @ 100cm
        ("tsoil-5", slice(204,211)), ## C @ 5cm
        ("tsoil-10", slice(212,219)), ## C @ 10cm
        ("tsoil-20", slice(220,227)), ## C @ 20cm
        ("tsoil-50", slice(228,235)), ## C @ 50cm
        ("tsoil-100", slice(236,243)), ## C @ 100cm
        ]

str_fields = ["wbanno", "qt-sfctemp", "qf-dswrf"]

nan_values = {
        "sfctemp":"-9999.0",
        "dswrf":"-99999",
        "precip":"-9999.0",
        "temp-mean":"-9999.0",
        "temp-final":"-9999.0",
        "vsm-5":"-9999.0",
        "vsm-10":"-9999.0",
        "vsm-20":"-9999.0",
        "vsm-50":"-9999.0",
        "vsm-100":"-9999.0",
        "tsoil-5":"-9999.0",
        "tsoil-10":"-9999.0",
        "tsoil-20":"-9999.0",
        "tsoil-50":"-9999.0",
        "tsoil-100":"-9999.0",
        }

def extract_uscrn_file(text_file:Path, skip_existing=True):
    """ """
    print(f"Extracting {text_file.name}")
    values = {f[0]:[] for f in fields}
    for l in text_file.open("r").readlines():
        tmpl = l.replace("\n","")
        for f,s in fields:
            values[f].append(tmpl[s])
    for k,v in values.items():
        if k=="utc-datetime":
            values[k] = [
                    int(datetime.strptime(t,"%Y%m%d %H%M").strftime("%s"))
                    for t in v
                    ]
        elif k in str_fields:
            continue
        else:
            try:
                values[k] = np.asarray(list(map(
                    float, [[p,np.nan][p==nan_values.get(k,np.nan)] for p in v]
                    )))
            except Exception as e:
                print(k, values[k])
                raise e

    return values

if __name__=="__main__":
    proj_root_dir = Path("/Users/mtdodson/desktop/soilm-in-situ")
    uscrn_file_dir = proj_root_dir.joinpath("uscrn")
    pkl_dir = proj_root_dir.joinpath("uscrn-pkls-new")
    combined_pkl_dir = proj_root_dir.joinpath("uscrn-pkls-combined")
    init_idx_json = proj_root_dir.joinpath("uscrn-valid-init-idxs_24hr.json")
    skip_existing = False
    contiguous_window_min_size = 24 ## must have 24 contiguous hours
    min_spacing_hours = 7

    ## extract pkls from text files
    '''
    for f in uscrn_file_dir.iterdir():
        year_state,*location = f.stem.split("_")
        _,year,state = year_state.split("-")
        locale = "-".join(location).lower().replace(".","")
        pkl_path = pkl_dir.joinpath(
                f"uscrn_{state.lower()}_{locale}_{year}.pkl")
        if pkl_path.exists() and skip_existing:
            print(f"Exists: {pkl_path.as_posix()}")
            continue
        try:
            pkl.dump(extract_uscrn_file(f), pkl_path.open("wb"))
        except Exception as e:
            print(e)
    '''

    ## combine pickles per locale
    '''
    pkl_paths = list(pkl_dir.iterdir())
    pdict = {}
    for p in pkl_paths:
        _,state,locale,year = p.stem.split("_")
        if state not in pdict.keys():
            pdict[state] = {}
        if locale not in pdict[state].keys():
            pdict[state][locale] = []
        pdict[state][locale].append(p)
        pdict[state][locale] = list(sorted(pdict[state][locale]))
    pprint(pdict)

    fkeys = [f[0] for f in fields]
    for s in pdict.keys():
        for l in pdict[s].keys():
            all_data = {}
            years = [p.stem.split("_")[-1] for p in pdict[s][l]]
            new_pkl_path = combined_pkl_dir.joinpath(
                    f"uscrn_{s}_{l}_{years[0]}-{years[-1]}.pkl")
            if skip_existing and new_pkl_path.exists():
                print(f"Skipping: {new_pkl_path.name}")
                continue
            for p in pdict[s][l]:
                tmpd = pkl.load(p.open("rb"))
                for k in fkeys:
                    if k in all_data.keys():
                        if k in ("utc-datetime", *str_fields):
                            all_data[k] = all_data[k] + tmpd[k]
                        else:
                            all_data[k] = np.concatenate(
                                    [all_data[k], tmpd[k]], axis=0)
                    else:
                        all_data[k] = tmpd[k]
            pkl.dump(all_data, new_pkl_path.open("wb"))
            print(f"Generated: {new_pkl_path.as_posix()}")
    '''

    ## identify contiguous strings of entirely valid data
    cwdw = contiguous_window_min_size
    '''
    fkeys = [f[0] for f in fields]
    valid_init_ixs = {}
    for tmpp in combined_pkl_dir.iterdir():
        pd = pkl.load(tmpp.open("rb"))
        masks = []
        for k in fkeys:
            if k=="utc-datetime":
                tdiffs = np.diff(pd[k])
                m_contig = (tdiffs-3600)**2 < 100
            elif k in str_fields:
                continue
            masks.append(np.isfinite(pd[k]))
            m_all = np.all(np.stack(masks, axis=1), axis=1)
        print(f"{np.count_nonzero(m_all):<6}",tmpp.as_posix())
        valid_init_ixs[tmpp.name] = []
        for ix in range(0, m_contig.size-cwdw):
            if np.all(m_contig[ix:ix+cwdw]) & np.all(m_all[ix:ix+cwdw+1]):
                valid_init_ixs[tmpp.name].append(ix)
    json.dump(valid_init_ixs, init_idx_json.open("w"))
    '''

    ## Extract valid sequences to time series (N,W,F) for each station
    valid_init_ixs = json.load(init_idx_json.open("r"))
    cpkl_names = list(valid_init_ixs.keys())
    strdata = []
    fdata = []
    times = []
    sflag = []
    for pi,pk in enumerate(cpkl_names):
        pd = pkl.load(combined_pkl_dir.joinpath(pk).open("rb"))
        fkeys = [fk for fk in pd.keys()
                 if fk not in ("utc-datetime",*str_fields)]

        ## identify starting indeces that are sufficiently spaced
        vixs = sorted(valid_init_ixs[pk]) ## valid initial indeces
        svixs = [] ## spaced-out valid initial indeces
        cur_ix = -9999
        for ix in vixs:
            if ix - cur_ix >= min_spacing_hours:
                svixs.append(slice(ix, ix+cwdw))
                cur_ix = ix

        tmp_fdata = []
        tmp_times = []
        tmp_strdata = []
        tmp_sflag = []
        for s in svixs:
            tmp_fdata.append(np.stack([pd[fk][s] for fk in fkeys], axis=-1))
            tmp_times.append([
                    datetime.fromtimestamp(t) for t in pd["utc-datetime"][s]
                    ])
            tmp_strdata = list(zip(*[pd[sk][s] for sk in str_fields]))
            tmp_sflag.append(np.full(s.stop-s.start, pi))

        fdata.append(np.stack(tmp_fdata, axis=0)) ## (N, W, F)
        times.append(np.stack(tmp_times, axis=0)) ## (N, W)
        strdata.append(tmp_strdata) ## (N, W, Fs)
        sflag.append(tmp_sflag) ## (N,)
    fdata = np.concatenate(fdata, axis=0)
    times = np.concatenate(times, axis=0)
    sflag = np.concatenate(sflag, axis=0)

    labels = (fkeys, str_fields, cpkl_names)
    data = (fdata, strdata, sflag, times)
    pkl_path = proj_root_dir.joinpath(
            f"uscrn_samples_{cwdw}h_{min_spacing_hours}p.pkl")
    pkl.dump((labels,data), pkl_path.open("wb"))
