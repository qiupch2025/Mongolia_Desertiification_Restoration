import os
import numpy as np
from netCDF4 import Dataset

# Data directory
datadir = '/groups/lzu_public/home/qiupch2023/data_wrfout_no_mongolia/'
filelist = sorted([f for f in os.listdir(datadir) if f.startswith('wrfout_d01')])

k = len(filelist)

# Initialize empty arrays
dust = edust = pm10 = tsk = wind_speed = wind_stress = None
time = 0
rho_air = 1.2  # Air density (kg/m³)

print(f'Processing {k} files...')
for i, filename in enumerate(filelist, start=1):
    ncFilePath = os.path.join(datadir, filename)
    print(f'Reading file {i}/{k}: {ncFilePath}')

    with Dataset(ncFilePath, 'r') as nc:
        # Read variables
        dust1 = np.array(nc.variables['DUSTLOAD_1'][:])
        dust2 = np.array(nc.variables['DUSTLOAD_2'][:])
        dust3 = np.array(nc.variables['DUSTLOAD_3'][:])
        dust4 = np.array(nc.variables['DUSTLOAD_4'][:])
        dust5 = np.array(nc.variables['DUSTLOAD_5'][:])
        edust1 = np.array(nc.variables['EDUST1'][:])
        edust2 = np.array(nc.variables['EDUST2'][:])
        edust3 = np.array(nc.variables['EDUST3'][:])
        edust4 = np.array(nc.variables['EDUST4'][:])
        edust5 = np.array(nc.variables['EDUST5'][:])
        pm10_data = np.array(nc.variables['PM10'][:, 0, :, :])
        tsk_data = np.array(nc.variables['TSK'][:])
        u10 = np.array(nc.variables['U10'][:])
        v10 = np.array(nc.variables['V10'][:])
        ust = np.array(nc.variables['UST'][:])  # Friction velocity

        # Calculate variables
        dust_sum = dust1 + dust2 + dust3 + dust4 + dust5
        edust_sum = edust1 + edust2 + edust3 + edust4 + edust5
        wind_speed_data = np.sqrt(u10**2 + v10**2)
        wind_stress_data = rho_air * ust * wind_speed_data

        time1 = dust5.shape[0]
        time += time1

        # Concatenate data
        if dust is None:
            dust = dust_sum
            edust = edust_sum
            pm10 = pm10_data
            tsk = tsk_data
            wind_speed = wind_speed_data
            wind_stress = wind_stress_data
        else:
            dust = np.concatenate((dust, dust_sum), axis=0)
            edust = np.concatenate((edust, edust_sum), axis=0)
            pm10 = np.concatenate((pm10, pm10_data), axis=0)
            tsk = np.concatenate((tsk, tsk_data), axis=0)
            wind_speed = np.concatenate((wind_speed, wind_speed_data), axis=0)
            wind_stress = np.concatenate((wind_stress, wind_stress_data), axis=0)

    print(f'File {i} processed, cumulative timesteps: {time}')

# Read latitude and longitude data
with Dataset(ncFilePath, 'r') as nc:
    lon = np.array(nc.variables['XLONG'][0, :, :])
    lat = np.array(nc.variables['XLAT'][0, :, :])

lat_l, lon_l = lon.shape
time_l = dust.shape[0]

output_file = 'dust_shao04_ideal_new.nc'
print(f'Creating NetCDF file: {output_file}')

# Create NetCDF file
with Dataset(output_file, 'w', format='NETCDF4') as nc_out:
    nc_out.createDimension('lon', lon_l)
    nc_out.createDimension('lat', lat_l)
    nc_out.createDimension('time', time_l)

    nc_out.createVariable('lon', 'f4', ('lat', 'lon'))[:, :] = lon
    nc_out.createVariable('lat', 'f4', ('lat', 'lon'))[:, :] = lat
    nc_out.createVariable('dust', 'f4', ('time', 'lat', 'lon'))[:, :, :] = dust
    nc_out.createVariable('edust', 'f4', ('time', 'lat', 'lon'))[:, :, :] = edust
    nc_out.createVariable('pm10', 'f4', ('time', 'lat', 'lon'))[:, :, :] = pm10
    nc_out.createVariable('tsk', 'f4', ('time', 'lat', 'lon'))[:, :, :] = tsk
    nc_out.createVariable('wind_speed', 'f4', ('time', 'lat', 'lon'))[:, :, :] = wind_speed
    nc_out.createVariable('wind_stress', 'f4', ('time', 'lat', 'lon'))[:, :, :] = wind_stress

print('✅ NetCDF file created successfully!')
