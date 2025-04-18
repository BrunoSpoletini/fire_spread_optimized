import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv("resultadosComp1.csv")

landscapes = df['Landscape'].unique()

for landscape in landscapes:
    partes = landscape.strip().split('/')
    name = partes[-1]

    df_l = df[df['Landscape'] == landscape]

    pivot = df_l.pivot(index='Optimizador', columns='Compilador', values='Celdas quemadas por segundo')

    ax = pivot.plot(kind='bar', figsize=(8, 6), rot=0)
    ax.set_title(f'{name}\nCeldas incendiadas por microsegundo vs Compilador y Optimizador')
    ax.set_xlabel('Optimizador')
    ax.set_ylabel('Celdas incendiadas por microsegundo')
    ax.legend(title='Compilador')

    plt.tight_layout()
    plt.savefig("./comp1/" + name + ".png")