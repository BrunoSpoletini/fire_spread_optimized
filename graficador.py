import pandas as pd
import matplotlib.pyplot as plt


# Leer los datos desde un archivo CSV llamado "datos.csv"
df = pd.read_csv("resultados.csv")

# Obtener los diferentes landscapes presentes en los datos
landscapes = df['Landscape'].unique()

for landscape in landscapes:
    # Filtrar datos por landscape
    df_l = df[df['Landscape'] == landscape]

    # Construir una tabla dinámica:
    #   - Filas: Optimizador
    #   - Columnas: Compilador
    #   - Valores: Celdas quemadas por segundo
    pivot = df_l.pivot(index='Optimizador', columns='Compilador', values='Celdas quemadas por segundo')

    # Crear el gráfico de barras comparativo
    ax = pivot.plot(kind='bar', figsize=(8, 6), rot=0)
    ax.set_title(f'Comparación de celdas quemadas por segundo\nLandscape: {landscape}')
    ax.set_xlabel('Optimizador')
    ax.set_ylabel('Celdas quemadas por segundo')
    ax.legend(title='Compilador')

    print("Test")

    # Ajustar la visualización y mostrar el gráfico
    plt.tight_layout()
    plt.show()
    plt.savefig("nombre.png")
    input("Presiona Enter para finalizar...")