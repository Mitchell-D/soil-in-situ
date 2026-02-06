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

#### before executing

- Acquire and unzip data from the ISMN database
- Modify `ismn_stations_path` under `if __name__=="__main__":` to
  point to the downloaded ISMN data.
- Modify `ismn_available_json` and `station_pkl_dir` to point to new
  paths where station metadata and extracted data will be deposited,
  respectively.
- Modify `extract_networks` to specify which networks to consider,
  `var_mapping` to shorten data variable names, and `keep_meta` to
  specify a subset of metadata fields that are relevant to extract.
- Set `nworkers` to the maximum number of parallel processes to
  dispatch for extracting

#### output

The dictionary stored in the pkl file associated with each station
has the following fields:

- **network**: Unique string name for this station's network

- **station**: String name for this station, unique within network.

- **location**: 3-tuple (lat, lon, elevation (meters))

- **depths**: List of 2-tuples (min, max) depths in meters measured
  by sensors.  The second element of the sensor 3-tuples in 'labels'
  provides the index in this list of that sensor's depth profile.

- **times**: List of size 10 strings like "%Y%m%d%H" (YYYYmmddHH)
  indicating the hourly time of each data point in the array.

- **labels**: List of 3-tuples uniquely identifying each sensor
  implemented by this station like: (feat:str, level:int, count:int)
  where *feat* indicates what is measured, *level* is an index
  mapping to the measurement location in 'depths', and *count*
  is an integer incremented when the station has multiple sensors
  measuring the same feature at the same depth.

- **data**: Numpy array of data values shaped like (T, L) for T
  timesteps (corresponding to T string entries in 'times'), and L
  sensor features (corresponding to L entries in 'labels').

- **masks**: Length L list of 2-tuples (ismn\_flags, network\_flags)
  for each of the sensors in 'labels'. Each of the sets of flags is
  a 1d array of characters. Data features definitely have an
  ismn\_flags array, but not all networks independently provide
  flags. For info on ISMN flag values see [this text file][3], and
  for network-specific flag values see [this file][4].

- **duplicate_times**: Sometimes a timestep is repeated, though this
  is rare. This list contains the time strings of the duplicates.

- **texture**:

- **climate**:

- **saturation**:

- **sensor_meta**: Dictionary mapping each 3-tuple from the 'labels'
  field to a dictionary of metadata for that sensor. Most of this
  metadata is redundant to other extracted fields, and is mainly
  here to help refer to the original data format.

- **station_meta**: Dictionary mapping original metadata fields to
  their data, which includes climate class, full sensor names, soil
  properties, etc.

[1]:https://pypi.org/project/ismn/
[2]:https://ismn.earth/en/dataviewer
[3]:info/ISMN_qualityflags_description.txt
[4]:info/ISMN_network_flags_descriptions.txt
