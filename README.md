# soil-in-situ

## dependencies

- [ISMN parser package][1]

## extract\_ismn.py

Given a directory tree of ISMN networks from the [web portal][2],
extracts the metadata for all network/station combos, and dumps them
into a common JSON file.

Then, for each station, extracts a dictionary containing all data
and flag values on a consistent hourly grid for each of the sensors
supported by that station.

All the extracted data is dumped into a unique pickle file for each
network/station combo, which includes lists of duplicate times
and other ancillary information.

Extracted metadata includes the location, soil and sensor
characteristics, and valid observation time range.

Also, data quality flags include both those from the individual
networks (as described in the network metadata), as well as the
post-processed quality flags from ISMN.

**before executing**

- Acquire and unzip data from the ISMN database
- Modify `ismn_stations_path` under `if __name__=="__main__:"` to
  point to the downloaded ISMN data.
- Modify `ismn_available_json` and `station_pkl_dir` to point to new
  paths where station metadata and extracted data will be deposited,
  respectively.
- Modify `extract_networks` to specify which networks to consider,
  `var_mapping` to shorten data variable names, and `keep_meta` to
  specify a subset of metadata fields that are relevant to extract.

[1]:https://pypi.org/project/ismn/
[2]:https://ismn.earth/en/dataviewer
