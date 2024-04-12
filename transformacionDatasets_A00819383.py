import pandas as pd
import glob
import os

path = "C:\\Users\\Alejandro\\Desktop\\transformacionDatasets"

# Asegúrate de que tu archivo 'merged_monterrey3.csv' esté en el mismo directorio que el script,
# o proporciona la ruta completa al archivo.
df_master = pd.read_csv('C:\\Users\\Alejandro\\Desktop\\InvestigacionIA\\merged_monterrey3.csv')
df_master['Fecha'] = df_master['Fecha'].astype(str)
df_master['Estacion'] = df_master['Estacion'].astype(str)

file_paths = glob.glob(os.path.join(path, '*.csv'))
 
for file_path in file_paths:
    # Leer solo las primeras líneas para obtener la estación y el parámetro
    with open(file_path, 'r', encoding='ISO-8859-1') as file:
        first_line = file.readline()
        parts = first_line.split(',')
        station_name = parts[1].strip()
        parameter = parts[3].strip()
    
    # Ahora lee el archivo CSV omitiendo las primeras líneas con los metadatos
    df = pd.read_csv(file_path, skiprows=3, encoding='ISO-8859-1')
    
    # Actualizar los valores en df_master con los de df
    for index, row in df.iterrows():
        if index < len(df_master):
            fecha = row['Fecha']
            valor = row['Valor']
            df_master.at[index, 'Fecha'] = fecha
            df_master.at[index, parameter] = valor
            df_master.at[index, 'Estacion'] = station_name
            print('Actualizando fila', index, 'con', parameter, '=', valor, 'de', station_name, 'en', fecha)
        else:
            # Si df tiene más filas que df_master, detener el bucle
            break
       
# Después de procesar todos los archivos CSV, puedes querer guardar tu df_master:
df_master.to_csv('df_master.csv', index=False)
print('Proceso terminado')
print('El archivo df_master.csv ha sido creado')
