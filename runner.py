import subprocess
import os
import csv

#compiladores = ["g++", "clang++ -std=c++17"]
compiladores = ["g++", "clang++ -std=c++17", "icpx", "nvc++"]
optimizadores = ["-O0", "-O1", "-O2", "-O3", "-Ofast"]
landscapes = [("./data/2000_8", 16), ("./data/1999_27j_S", 16), ("./data/2015_50", 4)]

def maxCeldasPorSeg(res):
    sumaCeldas = []
    for linea in res.stderr.splitlines():
        if "celdas incendiadas por microsegundo" in linea:
            partes = linea.strip().split()
            sumaCeldas.append(float(partes[-1]))
    return max(sumaCeldas)

def main():

    resultados = dict()

    with open('resultados.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(["Compilador", "Optimizador", "Landscape", "Celdas quemadas por microsegundo"])
        for landscape in landscapes:
            for comp in compiladores:

                for opt in optimizadores:
                    os.system("make clean && make CXX=" + comp + " CXXOPT=" + opt)
                    comandoPerf = ["/usr/lib/linux-tools/5.15.0-134-generic/perf", "stat", "-r", str(landscape[1]), "./graphics/fire_animation_data" , landscape[0], "0"]
                    resPerf = subprocess.run(comandoPerf, stderr=subprocess.PIPE, text=True)

                    writer.writerow([comp, opt, landscape[0], maxCeldasPorSeg(resPerf)])

if __name__ == "__main__":
    main()