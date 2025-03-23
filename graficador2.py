import pandas as pd
import matplotlib.pyplot as plt

def plot_comparacion_celdas(csv_path):
    # Leer el archivo CSV
    df = pd.read_csv(csv_path)
    
    # Calcular el promedio de "Celdas quemadas por microsegundo" para cada optimizador
    promedio_opt = df.groupby("Optimizador")["Celdas quemadas por microsegundo"].mean()
    # Ordenar de mayor a menor (la mayor cantidad implica mejor rendimiento)
    optimizer_order = promedio_opt.sort_values(ascending=True).index.tolist()
    
    # Crear la tabla dinámica: filas = Optimizador, columnas = Landscape,
    # valores = Celdas quemadas por microsegundo
    pivot = df.pivot(index="Optimizador", columns="Landscape", values="Celdas quemadas por microsegundo")
    
    # Eliminar el prefijo "./data/" de los nombres de las columnas
    pivot.columns = pivot.columns.str.replace(r"^\.\/data\/", "", regex=True)
    
    # Reordenar los optimizadores según el promedio (más rápido primero)
    pivot = pivot.reindex(optimizer_order)
    
    # Ahora se define el orden deseado sin el prefijo
    orden_landscapes = ["2000_8", "1999_27j_S", "2015_50"]
    pivot = pivot[orden_landscapes]
    
    # Crear el gráfico de barras agrupadas
    ax = pivot.plot(kind='bar', figsize=(10, 6))
    ax.set_xlabel("Opción de optimización")
    ax.set_ylabel("Celdas quemadas por microsegundo")
    ax.set_title("Comparación de celdas quemadas por microsegundo usando g++ -Ofast")
    ax.legend(title="Landscape", loc="lower left")
    
    # Modificar las etiquetas del eje x para que no muestren "-Ofast"
    nuevos_labels = [opt.replace("-Ofast", "").strip() for opt in pivot.index]
    
    nuevos_labels = ["base" if label == "" else label for label in nuevos_labels]
    ax.set_xticklabels(nuevos_labels)
    
    # Añadir etiquetas con el valor encima de cada barra
    for container in ax.containers:
        ax.bar_label(container, fmt="%.2f")
    
    plt.tight_layout()
    plt.savefig("grafica2v2.png")

# Ejemplo de uso:
plot_comparacion_celdas("resultadosComp2.csv")
