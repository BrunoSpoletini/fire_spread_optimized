import pandas as pd
import matplotlib.pyplot as plt
import os

# Crear carpeta si no existe
output_dir = "./comp3"
#output_dir = "./comp3Atom"
os.makedirs(output_dir, exist_ok=True)

# Lista de archivos con sus etiquetas
archivos = {
    "Atom": "resultadosComp3Atom.csv",
    "Local1": "resultadosComp3Local1.csv",
    "Local2": "resultadosComp3Local2.csv"
    #"Reduced": "resultadosComp3AtomReduced.csv",
    #"Original": "resultadosComp3AtomReduced.csv",
}

# Diccionario para almacenar los datos
data_frames = {}

# Leer los archivos y agregarlos al diccionario
for etiqueta, archivo in archivos.items():
    df = pd.read_csv(archivo)
    df["Fuente"] = etiqueta  # Agregar columna para identificar la fuente
    data_frames[etiqueta] = df

# Unir todos los dataframes en uno solo
df_total = pd.concat(data_frames.values())

# Eliminar el prefijo "./data/" en Landscape
df_total["Landscape"] = df_total["Landscape"].str.replace(r"^\.\/data\/", "", regex=True)

# Reemplazar valores vacíos en la columna Optimizador con "base"
df_total["Optimizador"] = df_total["Optimizador"].fillna("base")

# Obtener los nombres únicos de los landscapes
landscapes = df_total["Landscape"].unique()

# Crear un gráfico para cada Landscape
for landscape in landscapes:
    df_landscape = df_total[df_total["Landscape"] == landscape]
    
    # Crear tabla pivote: filas = Optimizador, columnas = Fuente, valores = promedio de Celdas quemadas
    pivot = df_landscape.pivot_table(index="Optimizador", columns="Fuente", values="Celdas quemadas por microsegundo", aggfunc="mean")
    
    # Ordenar las optimizaciones para que "base" esté primero
    pivot = pivot.reindex(sorted(pivot.index, key=lambda x: (x != "base", x)))
    
    # Ordenar las fuentes de datos según la más rápida (menor promedio)
    pivot = pivot[pivot.mean().sort_values().index]
    
    # Graficar
    ax = pivot.plot(kind="bar", figsize=(10, 6))
    ax.set_xlabel("Optimizador")
    ax.set_ylabel("Celdas incendiadas por microsegundo")
    ax.set_title(f"{landscape}\nCeldas incendiadas por microsegundo vs Parámetros")
    ax.legend(title="Fuente de datos", loc="upper left")
    
    # Agregar etiquetas de valores en las barras
    for container in ax.containers:
        ax.bar_label(container, fmt="%.2f")
    
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    # Guardar la imagen
    filename = os.path.join(output_dir, f"{landscape}.png")
    plt.savefig(filename)
    plt.close()
