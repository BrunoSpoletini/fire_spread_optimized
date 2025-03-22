import pandas as pd
import matplotlib.pyplot as plt

def plot_comparacion_celdas(csv_path):
    # Leer el archivo CSV
    df = pd.read_csv(csv_path)
    
    # Crear una tabla dinámica (pivot):
    #   - Índice: Optimizador
    #   - Columnas: Landscape
    #   - Valores: Celdas quemadas por microsegundo
    pivot = df.pivot(index="Optimizador", columns="Landscape", values="Celdas quemadas por microsegundo")
    
    # Reordenar los optimizadores si se desea un orden específico
    orden_optimizadores = ["-Ofast", "-Ofast -march=native", "-Ofast -mtune=native", "-Ofast -flto", "-Ofast -funroll-loops"]
    pivot = pivot.reindex(orden_optimizadores)
    
    # Crear el gráfico de barras comparativo
    ax = pivot.plot(kind='bar', figsize=(10, 6))
    ax.set_xlabel("Optimizador")
    ax.set_ylabel("Celdas quemadas por microsegundo")
    ax.set_title("Comparación de celdas quemadas por microsegundo")
    ax.legend(title="Landscape")
    
    # Añadir etiquetas de valor encima de cada barra
    for container in ax.containers:
        ax.bar_label(container, fmt='%.3f')
    
    plt.tight_layout()
    plt.savefig("grafica2.png")

# Ejemplo de uso:
plot_comparacion_celdas("resultadosComp2.csv")
