"""
Acquire SCAN data from USDA's REST API and store each station's data as a
separate pkl file alongside quality flags and feature information
"""
from pathlib import Path
import numpy as np
import json
from datetime import datetime,timedelta
import pickle as pkl

if __name__=="__main__":
    #proj_root_dir = Path("/Users/mtdodson/desktop/soilm-in-situ")
    proj_root_dir = Path("/rhome/mdodson/soil-in-situ")
    station_json_path = proj_root_dir.joinpath("data/scan-stations.json")
    datamenu_json_path = proj_root_dir.joinpath("data/scan-datamenu.json")
    pkl_dir = proj_root_dir.joinpath("data/scan/scan-pkls")
    init_idx_json = proj_root_dir.joinpath(
            "data/scan-valid-init-idxs_48hr.json")
    contiguous_window_min_size = 48
    min_spacing_hours = 17

    ## require temp, humidity, precipitation, wind, soilm at 2,4,8,20,40
    require_feats = [
            "DPTP::1", "PRCP::1",
            #"PREC::1",
            "SMS:-2:1", "SMS:-4:1", "SMS:-8:1", "SMS:-20:1", "SMS:-40:1",
            "TAVG::1", "WSPDV::1"
            ]

    valid_init_ixs = {}
    cwdw = contiguous_window_min_size
    for spp in pkl_dir.iterdir():
        pd = pkl.load(spp.open("rb"))
        ## only continue if all required features are present
        if not all(k in pd.keys() for k in require_feats):
            continue

        ## Make sure time axes align and collect quality masks for all featuers
        masks = []
        times = None
        valid_times = None
        feat_times = {}
        feat_ixs = {}
        for k in require_feats:
            ## Get times associated with valid points for this feature
            tmp_ixs,tmp_times = zip(*[
                    (i,datetime.fromtimestamp(t))
                    for i,(t,f) in enumerate(zip(
                        pd[k]["etimes"], pd[k]["flags"]))
                    if f=="V"
                    ])
            ## Round the times and create a serializable string for each
            tmp_times = [
                    t.strftime("%Y%m%d%H") if t.minute < 45
                    else (t+timedelta(hours=1)).strftime("%Y%m%d%H")
                    for t in tmp_times
                    ]
            ## Keep valid timesteps for this feature for mask creation and
            ## determine intersection of all valid timesteps.
            feat_times[k] = tmp_times
            feat_ixs[k] = tmp_ixs
            if valid_times is None:
                valid_times = set(tmp_times)
            else:
                valid_times = valid_times.intersection(set(tmp_times))
        '''
        valid_times = sorted([
            datetime.strptime(t, "%Y%m%d%H")
            for t in valid_times
            ])
        '''
        time_masks = []
        for k in require_feats:
            fixs = np.array([
                i for i,t in zip(feat_ixs[k],feat_times[k])
                if t in valid_times
                ])
            m_valid = np.full(len(pd[k]["etimes"]), False)
            m_valid[fixs] = True
            time_masks.append(np.array(pd[k]["etimes"])[m_valid])
            print(k, np.count_nonzero(m_valid))
        time_masks = np.stack(time_masks, axis=-1)
        print(time_masks)

    '''
        if times is None:
            times = np.array(pd[k]["etimes"])
        else:
            assert np.all(np.abs(times-np.array(pd[k]["etimes"])) < 60)
        masks.append([c=="V" for c in pd[k]["flags"]])

        ## Collect indeces that are valid in all features and contiguous with
        ## the following timestep
        m_valid = np.all(np.stack(masks, axis=-1), axis=-1)
        m_contig = (np.diff(times)-3600)**2 < 100
        valid_init_ixs = []
        for ix in range(0, m_contig.size-cwdw):
            if np.all(m_contig[ix:ix+cwdw]) & np.all(m_valid[ix:ix+cwdw+1]):
                valid_init_ixs[spp.name].append(ix)
        print(f"{valid_init_ixs[spp.name]} valid times in {spp.stem}")
    json.dump(valid_init_ixs, init_idx_json.open("w"))
    '''
