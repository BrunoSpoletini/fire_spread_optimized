import subprocess
import os
import csv

#optimizadores = ["-O1", "-O2", "-O3","-Ofast"]
#landscapes = [("../data/2000_8", 16), ("../data/1999_27j_S", 4)]#, ("../data/2015_50", 1)]
landscapes = [("../data/2000_8", 32), ("../data/1999_27j_S", 16), ("../data/2015_50", 4)]
#landscapes = [("../data/1999_27j_S", 4)]

def maxCeldasPorSeg(res):
    sumaCeldas = []
    for linea in res.stderr.splitlines():
        if "celdas incendiadas por microsegundo" in linea:
            partes = linea.strip().split()
            sumaCeldas.append(float(partes[-1]))
    return max(sumaCeldas)

def createCsvFile(filename, src, compiladores=None, optimizadores=None):
    with open(filename + ".csv", mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(["Landscape", "Celdas quemadas por microsegundo"])
        for landscape in landscapes:
            if src == "srcBase":
                os.system("cd .. && make -f MakefileBase clean && make -f MakefileBase")
            elif src == "srcVect":
                os.system("cd .. && make -f MakefileVect clean && make -f MakefileVect")
            elif src == "srcParalel":
                os.system("cd .. && make -f MakefileParalel clean && make -f MakefileParalel")

            # Para vos loquito
            # comandoPerf = ["perf", "stat", "-r", str(landscape[1]), "../graphics/fire_animation_data" , landscape[0], "0"]
            comandoPerf = ["/usr/lib/linux-tools/5.15.0-139-generic/perf", "stat", "-r", str(landscape[1]), "../graphics/fire_animation_data" , landscape[0], "0"]
            resPerf = subprocess.run(comandoPerf, stderr=subprocess.PIPE, text=True)
            writer.writerow([landscape[0], maxCeldasPorSeg(resPerf)])

def main():
    createCsvFile("Base", "srcBase")
    createCsvFile("Vectorizado", "srcVect")
    createCsvFile("Paralelo", "srcParalel")

if __name__ == "__main__":
    main()