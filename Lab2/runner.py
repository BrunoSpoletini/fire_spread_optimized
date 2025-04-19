import subprocess
import os
import csv

#optimizadores = ["-O1", "-O2", "-O3","-Ofast"]
landscapes = [("../data/2000_8", 16), ("../data/1999_27j_S", 4)]#, ("../data/2015_50", 1)]
#landscapes = [("../data/2000_8", 32), ("../data/1999_27j_S", 16), ("../data/2015_50", 4)]

def maxCeldasPorSeg(res):
    sumaCeldas = []
    for linea in res.stderr.splitlines():
        if "celdas incendiadas por microsegundo" in linea:
            partes = linea.strip().split()
            sumaCeldas.append(float(partes[-1]))
    return max(sumaCeldas)

def createCsvFile(filename, src, id="", compiladores=None, optimizadores=None):
    with open(filename + "_" + id + ".csv", mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(["Compilador", "Optimizador", "Landscape", "Celdas quemadas por microsegundo"])
        for landscape in landscapes:
            for comp in compiladores:
                for opt in optimizadores:
                    os.system("cd .. && make clean && make SOURCE=\"" + src + "\" CXX=\"" + comp + "\" CXXOPT=\"" + opt + "\"")
                    comandoPerf = ["perf", "stat", "-r", str(landscape[1]), "../graphics/fire_animation_data" , landscape[0], "0"]
                    resPerf = subprocess.run(comandoPerf, stderr=subprocess.PIPE, text=True)
                    writer.writerow([comp, opt, landscape[0], maxCeldasPorSeg(resPerf)])

def main():
    compiladoresOptiAutomatica = ["g++"]
    optimizadoresAuto = [   "-O0 -march=native -flto -ftree-vectorize",
                            "-O1 -march=native -flto -ftree-vectorize", #deberia agregar -flto
                            "-O2 -march=native -flto -ftree-vectorize", 
                            "-O3 -march=native -flto -ftree-vectorize",
                            "-Ofast -march=native -flto -ftree-vectorize"]

    compiladoresVectManual = ["g++"]
    optimizadoresVectManual = ["-O0 -march=native -flto -fopenmp-simd", 
                               "-O1 -march=native -flto -fopenmp-simd", 
                               "-O2 -march=native -flto -fopenmp-simd", 
                               "-O3 -march=native -flto -fopenmp-simd",
                                "-Ofast -march=native -flto -fopenmp-simd"]


    createCsvFile("resultadosBase", "srcBase", "1", compiladoresOptiAutomatica, optimizadoresAuto)

    #createCsvFile("resultadosSrc", "srcVect", "6", compiladoresVectManual, optimizadoresVectManual)


if __name__ == "__main__":
    main()