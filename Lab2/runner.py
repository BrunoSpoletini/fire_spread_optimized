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

def createCsvFile(filename, src, compiladores=None, optimizadores=None):
    with open(filename + "_" + ".csv", mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(["Compilador", "Optimizador", "Landscape", "Celdas quemadas por microsegundo"])
        for landscape in landscapes:
            for comp in compiladores:
                for opt in optimizadores:
                    if src == "srcBase":
                        os.system("cd .. && make -f MakefileBase clean && make -f MakefileBase CXX=\"" + comp + "\" CXXOPT=\"" + opt + "\"")
                    else:
                        os.system("cd .. && make -f MakefileVect clean && make -f MakefileVect CXX=\"" + comp + "\" CXXOPT=\"" + opt + "\"")

                    comandoPerf = ["perf", "stat", "-r", str(landscape[1]), "../graphics/fire_animation_data" , landscape[0], "0"]
                    resPerf = subprocess.run(comandoPerf, stderr=subprocess.PIPE, text=True)
                    writer.writerow([comp, opt, landscape[0], maxCeldasPorSeg(resPerf)])


def get_vectorization_info():
    optims = ["-O0", "-O1", "-O2", "-O3", "-Ofast"]
    for opt in optims:
        os.system("cd .. && make -f MakefileVect clean && make -f MakefileVect CXX=g++ CXXOPT=\"" + opt +" -fopt-info-vec=vector_srcVect_" + opt + ".log\"")

def get_vectorization_info_not_vectorized():
    optims = ["-O0", "-O1", "-O2", "-O3", "-Ofast"]
    for opt in optims:
        os.system("cd .. && make -f MakefileVect clean && make -f MakefileVect CXX=g++ CXXOPT=\"" + opt +" -fopt-info-missed=not_vectorized_srcVect_" + opt + ".log\"")


def main():
    compiladores = ["g++"]
    optimizadores = ["-O0", "-O1", "-O2", "-O3","-Ofast"]

    optimizadoresBase = list (map(lambda x: x + " -flto", optimizadores))
    optimizadoresVectAuto = list (map(lambda x: x + " -flto -ftree-vectorize", optimizadores))
    optimizadoresVectManual = list (map(lambda x: x + " -flto", optimizadores))



    # createCsvFile("resultadosBase", "srcBase", compiladores, optimizadoresBase)
    # createCsvFile("resultadosVectAuto", "srcBase", compiladores, optimizadoresVectAuto)
    # createCsvFile("resultadosVectManual_SinFlagDeSimd", "srcVect", compiladores, optimizadoresVectManual)
    createCsvFile("resultadosVectManual_SOA", "srcVect", compiladores, optimizadoresVectManual)

    # get_vectorization_info_not_vectorized()
    # get_vectorization_info()

if __name__ == "__main__":
    main()