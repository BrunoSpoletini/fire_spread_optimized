import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Cargar los datos desde los archivos CSV
df_atom = pd.read_csv('resultadosComp3Atom.csv')
df_local = pd.read_csv('resultadosComp3Local.csv')

# Unir (merge) los dataframes por las columnas comunes
df_merged = pd.merge(
    df_atom, df_local,
    on=['Compilador', 'Optimizador', 'Landscape'],
    suffixes=('_Atom', '_Local')
)

# Calcular el promedio de los dos sistemas para cada fila
df_merged['Promedio'] = df_merged[['Celdas quemadas por microsegundo_Atom',
                                   'Celdas quemadas por microsegundo_Local']].mean(axis=1)

# Filtrar solo el landscape '2015_50'
df_filtered = df_merged[df_merged['Landscape'] == '2015_50']

# Crear una tabla pivote: índice = Landscape, columnas = Optimizador, valores = Promedio
pivot = df_filtered.pivot(index='Landscape', columns='Optimizador', values='Promedio')

# Configurar la gráfica de barras agrupadas
fig, ax = plt.subplots(figsize=(8, 6))
landscapes = pivot.index.tolist()  # Será una lista con un solo elemento: ['2015_50']
optimizers = pivot.columns.tolist()

x = np.arange(len(landscapes))   # Posiciones para cada Landscape (solo 1)
width = 0.2                      # Ancho de cada barra

# Para centrar las barras en el grupo, calculamos un offset según la cantidad de optimizadores
n_opts = len(optimizers)
offset = (n_opts - 1) / 2.0

# Dibujar una barra para cada optimizador
for i, opt in enumerate(optimizers):
    values = pivot[opt].values
    ax.bar(x + (i - offset) * width, values, width, label=str(opt))

# Configurar etiquetas y título
ax.set_xticks(x)
ax.set_xticklabels(landscapes, rotation=0)
ax.set_ylabel('Celdas quemadas por microsegundo (Promedio)')
ax.set_title('Comparación por optimizador en Landscape 2015_50')
ax.legend(title='Optimizador')

plt.tight_layout()
plt.grid(True)

# Guardar la gráfica en un archivo PNG
plt.savefig("grafica3.png")


