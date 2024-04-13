import numpy as np
import netCDF4
import xarray as xr
import warnings
import matplotlib.pyplot as plt
import glob
import os
import pickle

# Función para cargar datos
def load_data():
    try:
        with open('data_dict.pickle', 'rb') as handle:
            # Load the entire dictionary stored in the pickle file
            dic_container = pickle.load(handle)

        # Remove the dictionaries you no longer need
        dic_container.pop('dic1', None)  # Remove 'dic1' if it exists
        dic_container.pop('dic2', None)  # Remove 'dic2' if it exists

        return dic_container
    except (FileNotFoundError, EOFError):
        # If there's no file, return an empty dictionary
        return {}


# Intenta cargar el diccionario al inicio del script
try:
    data_dict_f = load_data()
except (FileNotFoundError, EOFError):
    data_dict_f = {}  # Si no existe, comienza con un diccionario vacío


#%%
#lvt_ds = netCDF4.Dataset(r"C:\Users\gusta\Downloads\OR_ABI-L2-LVTPC-M6_G16_s20220020101173_e20220020103546_c20220020106347.nc")
#cod_ds = netCDF4.Dataset(r"D:\GOES 16\Database\2018\001\01\OR_ABI-L2-CODC-M3_G16_s20180010107197_e20180010109569_c20180010111261.nc")
#tpw_ds = netCDF4.Dataset(r"C:\Users\gusta\Downloads\OR_ABI-L2-TPWC-M6_G16_s20210030201173_e20210030203546_c20210030205408.nc")
#lst_ds = netCDF4.Dataset(r"C:\Users\gusta\Downloads\OR_ABI-L2-LSTC-M6_G16_s20220020101173_e20220020103546_c20220020105278.nc")
#%%
def file_reading(dataset_path,size,dataset_name)->dict:
    center_lat_lon = (25.6866, -100.3161)  # Monterrey's geographical coordinates
    crop_size_lat_lon = (4, 4)  # Degrees of latitude and longitude for cropping
    nc_files = glob.glob(os.path.join(dataset_path, '*.nc'))
    variables = ['COD', 'LVT', 'TPW', 'LST']
    mty_lat_lon = (25.6866, -100.3161)  # Monterrey's geographical coordinates
    orientation = 'north'  # Orientation for the final image
    # Itera sobre cada archivo en la carpeta
    data_dict = {var: {} for var in variables}

    for file_path in nc_files:
        #ds = xr.open_dataset(file_path,decode_times=False)
        ds=netCDF4.Dataset(file_path,decodetimes=False)
        time_coverage_start = ds.getncattr('time_coverage_start')


        # Procesa y almacena los datos para cada variable
        for var in variables:
            if var in ds.variables:
                if var == 'COD':
                    processed_data = center_crop(ds, mty_lat_lon, size, var)
                    data_dict[var][time_coverage_start] = processed_data
                    save_data(data_dict,dataset_name)

                    data_dict[var][time_coverage_start] = processed_data
                if var == 'TPW':
                    processed_data = center_crop(ds, mty_lat_lon, size, var)
                    data_dict[var][time_coverage_start] = processed_data
                    save_data(data_dict,dataset_name)

                    data_dict[var][time_coverage_start] = processed_data
                if var == 'LST':
                    processed_data = center_crop(ds, mty_lat_lon, size, var)
                    data_dict[var][time_coverage_start] = processed_data
                    save_data(data_dict,dataset_name)

                    data_dict[var][time_coverage_start] = processed_data
                if var == 'LVT':
                        pressures=[0,25,50,75,100]
                        for pressure_level in pressures:
                            processed_data = center_crop(ds, mty_lat_lon, size, var, pressure_level)
                            data_dict[var][time_coverage_start] = processed_data
                            data_dict[var][time_coverage_start] = processed_data
                            save_data(data_dict,dataset_name)
                print(f'Processed {var} for {time_coverage_start}')
        
        ds.close()
    print('Data processed and saved.', data_dict)
    return data_dict

def save_data(partial_data, dataset_name):
    try:
        with open('data_dict.pickle', 'rb') as handle:
            existing_data = pickle.load(handle)
    except (FileNotFoundError, EOFError):
        existing_data = {}
    
    # Asegúrate de que el diccionario para este conjunto de datos exista.
    if dataset_name not in existing_data:
        existing_data[dataset_name] = {}
    
    # Para cada variable en los nuevos datos.
    for var_key, timestamps in partial_data.items():
        # Si la variable aún no existe en el diccionario, añádela completa.
        if var_key not in existing_data[dataset_name]:
            existing_data[dataset_name][var_key] = timestamps
        else:
            # Si la variable ya existe, actualiza solo las nuevas marcas de tiempo.
            for timestamp, value in timestamps.items():
                # Solo añade la marca de tiempo si no existe; no actualiza las existentes.
                if timestamp not in existing_data[dataset_name][var_key]:
                    existing_data[dataset_name][var_key][timestamp] = value
    
    # Guarda los datos, asegurándote de que 'dic1' y 'dic2' no estén incluidos.
    data_to_save = {k: v for k, v in existing_data.items() if k not in ['dic1', 'dic2']}
    
    with open('data_dict.pickle', 'wb') as handle:
        pickle.dump(data_to_save, handle, protocol=pickle.HIGHEST_PROTOCOL)


    with open('data_dict.pickle', 'wb') as handle:
        pickle.dump(data_to_save, handle, protocol=pickle.HIGHEST_PROTOCOL)


def lat_lon_square(center, size):
    half_size = size / 2
    lat_min = center[0] - half_size
    lat_max = center[0] + half_size
    lon_min = center[1] - half_size
    lon_max = center[1] + half_size

    return (lat_min, lat_max), (lon_min, lon_max)
#%%
# Function to get the x and y coordinates from a latitude and longitude square
def get_xy_from_latlon(ds, lats, lons)->tuple:
    lat1, lat2 = lats
    lon1, lon2 = lons

    lat = ds.lat.data
    lon = ds.lon.data

    x = ds.x.data
    y = ds.y.data

    x,y = np.meshgrid(x,y)

    x = x[(lat >= lat1) & (lat <= lat2) & (lon >= lon1) & (lon <= lon2)]
    y = y[(lat >= lat1) & (lat <= lat2) & (lon >= lon1) & (lon <= lon2)]

    return ((min(x), max(x)), (min(y), max(y)))
#%%
def calc_latlon(ds)->xr.Dataset:
    x = ds.x
    y = ds.y
    goes_imager_projection = ds.goes_imager_projection

    x,y = np.meshgrid(x,y)

    r_eq = goes_imager_projection.attrs["semi_major_axis"]
    r_pol = goes_imager_projection.attrs["semi_minor_axis"]
    l_0 = goes_imager_projection.attrs["longitude_of_projection_origin"] * (np.pi/180)
    h_sat = goes_imager_projection.attrs["perspective_point_height"]
    H = r_eq + h_sat

    a = np.sin(x)**2 + (np.cos(x)**2 * (np.cos(y)**2 + (r_eq**2 / r_pol**2) * np.sin(y)**2))
    b = -2 * H * np.cos(x) * np.cos(y)
    c = H**2 - r_eq**2

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        r_s = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)

    s_x = r_s * np.cos(x) * np.cos(y)
    s_y = -r_s * np.sin(x)
    s_z = r_s * np.cos(x) * np.sin(y)

    lat = np.arctan((r_eq**2 / r_pol**2) * (s_z / np.sqrt((H-s_x)**2 +s_y**2))) * (180/np.pi)
    lon = (l_0 - np.arctan(s_y / (H-s_x))) * (180/np.pi)

    ds = ds.assign_coords({
        "lat":(["y","x"],lat),
        "lon":(["y","x"],lon)
    })
    ds.lat.attrs["units"] = "degrees_north"
    ds.lon.attrs["units"] = "degrees_east"

    return ds
#%%
def center_crop(ds, center, size, var, pressure_level = 0)->np.ndarray:
    xds = xr.open_dataset(xr.backends.NetCDF4DataStore(ds),decode_times=False)

    xds.coords["x"] = xds.x
    xds.coords["y"] = xds.y
    xds.coords["goes_imager_projection"] = xds.goes_imager_projection

    xds = calc_latlon(xds)

    lat_min = xds.coords['lat'].min().values
    lat_max = xds.coords['lat'].max().values
    lon_min = xds.coords['lon'].min().values
    lon_max = xds.coords['lon'].max().values

    lats, lons = lat_lon_square(center, size)
    ((x1,x2), (y1, y2)) = get_xy_from_latlon(xds, lats, lons)
    subset = xds.sel(x=slice(x1, x2), y=slice(y2, y1))

    if var == 'COD':
        data_array = subset.COD
    elif var == 'TPW':
        data_array = subset.TPW
    elif var == 'LST':
        data_array = subset.LST
    elif var == 'LVT':
        data_array = subset.LVT.isel(pressure=pressure_level)
    else:
        print('Not a valid variable.')

    return data_array.values
#%%
def plot_center_crops(ds,var):
    mty_lat_lon = (25.6866, -100.3161)  # Monterrey's geographical coordinates
    size = 6 # Size of the cropping square in degrees
    if var == 'LVT':
        matrix = center_crop(ds, mty_lat_lon, size, 'LVT', pressure_level = 0)
        plt.imshow(matrix, cmap='viridis')
        plt.colorbar()
        plt.show()
        
        matrix = center_crop(ds, mty_lat_lon, size, 'LVT', pressure_level = 25)
        plt.imshow(matrix, cmap='viridis')
        plt.colorbar()
        plt.show()
        
        matrix = center_crop(ds, mty_lat_lon, size, 'LVT', pressure_level = 50)
        plt.imshow(matrix, cmap='viridis')
        plt.colorbar()
        plt.show()
        
        matrix = center_crop(ds, mty_lat_lon, size, 'LVT', pressure_level = 75)
        plt.imshow(matrix, cmap='viridis')
        plt.colorbar()
        plt.show()
        
        matrix = center_crop(ds, mty_lat_lon, size, 'LVT', pressure_level = 100)
        plt.imshow(matrix, cmap='viridis')
        plt.colorbar()
        plt.show()
    if var == 'COD':
        matrix = center_crop(ds, mty_lat_lon, size, 'COD')
        plt.imshow(matrix, cmap='viridis')
        plt.colorbar()
        plt.show()
    if var == 'TPW':    
        matrix = center_crop(ds, mty_lat_lon, size, 'TPW')
        plt.imshow(matrix, cmap='viridis')
        plt.colorbar()
        plt.show()
    if var == 'LST':    
        matrix = center_crop(ds, mty_lat_lon, size, 'LST')
        plt.imshow(matrix, cmap='viridis')
        plt.colorbar()
        plt.show()

if __name__ == '__main__':
    existing_data = load_data()
    print("Diccionario previamente guardado", existing_data)
    dataset_path = r"F:/dataSetsperProduct5"
    size1, size2 = 4, 6
    # Process and save data for size 4
    file_reading(dataset_path, size1, 'dataset_size_4')
    # Process and save data for size 6
    file_reading(dataset_path, size2, 'dataset_size_6')
    print("Diccionario guardado", existing_data)

