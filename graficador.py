import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv("resultados.csv")

landscapes = df['Landscape'].unique()

for landscape in landscapes:
    partes = landscape.strip().split('/')
    name = partes[-1]

    df_l = df[df['Landscape'] == landscape]

    pivot = df_l.pivot(index='Optimizador', columns='Compilador', values='Celdas quemadas por segundo')

    ax = pivot.plot(kind='bar', figsize=(8, 6), rot=0)
    ax.set_title(f'Comparación de celdas quemadas por segundo\nLandscape: {landscape}')
    ax.set_xlabel('Optimizador')
    ax.set_ylabel('Celdas quemadas por segundo')
    ax.legend(title='Compilador')

    plt.tight_layout()
    plt.savefig(name + ".png")