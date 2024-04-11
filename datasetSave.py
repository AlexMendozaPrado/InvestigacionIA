from pandas import date_range, Timedelta
from goes2go.data import goes_timerange
import os
import requests

start = "2019-12-05 18:00"
end = "2022-12-31 23:00"

save_directory = 'F:\\dataSetsperProduct4'
if not os.path.exists(save_directory):
    os.makedirs(save_directory)

dates = date_range(start=start, end=end, freq='H')
products = ["ABI-L2-CODC"]

for product in products:
    for current in dates:
        start_str = current.strftime("%Y-%m-%d %H:%M")
        end_str = (current + Timedelta(hours=1)).strftime("%Y-%m-%d %H:%M")
        
        try:
            datos = goes_timerange(
                start=start_str,
                end=end_str,
                satellite="goes16",
                product=product,
                return_as='filelist',
                download=False
            )
            
            if not datos.empty:
                file_to_download = datos.iloc[0]['file']
                file_to_download = file_to_download.split('noaa-goes16/', 1)[-1]  # Nombre del archivo
                file_name = file_to_download.split('/')[-1]# Nombre del archivo
                save_path = os.path.join(save_directory, file_name)
                
                # Check if the file already exists to avoid re-downloading it
                if not os.path.isfile(save_path):
                    download_url = f"https://noaa-goes16.s3.amazonaws.com/{file_to_download}"
                    
                    
                    # Descargar y guardar el archivo
                    response = requests.get(download_url)
                    with open(save_path, 'wb') as file:
                        file.write(response.content)
                    print(f"Archivo descargado y guardado en: {save_path}")
                    print("archivo de fecha: ", start_str + " - " + end_str)
                else:
                    print(f"El archivo ya existe: {save_path}")
                    
        except FileNotFoundError:
            print(f"No se encontró archivo para: {start_str} - {end_str} con el producto {product}.")
            continue
