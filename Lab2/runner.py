import subprocess
import os
import csv

compiladores = ["g++"]
optimizadores = ["-O1", "-O2", "-Ofast"]
landscapes = [("../data/2000_8", 16), ("../data/1999_27j_S", 4), ("../data/2015_50", 4)]

def maxCeldasPorSeg(res):
    sumaCeldas = []
    for linea in res.stderr.splitlines():
        if "celdas incendiadas por microsegundo" in linea:
            partes = linea.strip().split()
            sumaCeldas.append(float(partes[-1]))
    return max(sumaCeldas)

def createCsvFile(filename, src, id=""):
    with open(filename + "_" + id + ".csv", mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(["Compilador", "Optimizador", "Landscape", "Celdas quemadas por microsegundo"])
        for landscape in landscapes:
            for comp in compiladores:
                for opt in optimizadores:
                    os.system("cd .. && make clean && make SRC=\"" + src + "\" CXX=\"" + comp + "\" CXXOPT=\"" + opt + "\"")
                    comandoPerf = ["perf", "stat", "-r", str(landscape[1]), "../graphics/fire_animation_data" , landscape[0], "0"]
                    resPerf = subprocess.run(comandoPerf, stderr=subprocess.PIPE, text=True)
                    writer.writerow([comp, opt, landscape[0], maxCeldasPorSeg(resPerf)])

def main():
    # El ultimo parametro es para generar resultados con las distintas optimizaciones progresivas
    # (Cada optimizacion contiene a las anteriores)
    createCsvFile("resultadosVect_1", "srcVect", 1)


if __name__ == "__main__":
    main()