import subprocess
import os
import csv

#compiladores = ["g++"] #["g++", "clang++"]
compiladores = ["g++", "clang++"]
#compiladores = ["clang++"]
#optimizadores = ["-O0"]
optimizadores = ["-O0", "-O1", "-O2", "-O3", "-Ofast"]
#landscapes = [("./data/2000_8", 1)]
landscapes = [("./data/2000_8", 16), ("./data/1999_27j_S", 16), ("./data/2015_50", 4)]

# Filtra el tiempo de usuario promedio de todas las ejecuciones
def filtrarPerf(res):
    tiempo_promedio = None
    for linea in res.stderr.splitlines():
        if "seconds time elapsed" in linea:
            partes = linea.strip().split()
            tiempo_promedio = partes[0]
            break
    return tiempo_promedio

# Numero maximo de celdas quemadas por segundo de varias ejecuciones
def maxCeldasPorSeg(res):
    lines = res.stderr.splitlines()
    sumaCeldas = []
    for line in lines:
        partes = line.strip().split()
        num = partes[-1]
        sumaCeldas.append(float(num))
    return max(sumaCeldas)

def main():

    resultados = dict()

    # comando = ["/usr/lib/linux-tools/5.15.0-134-generic/perf", "stat", "-r", "1", "./graphics/fire_animation_data" , "./data/2000_8", "0"]
    # res = subprocess.run(comando, stderr=subprocess.PIPE, text=True)
    # tiempo_promedio_user = filtrarPerf(res)
    # print("Tiempo medio de ejecución:", tiempo_promedio_user, "segundos")
    # # Cuando usariamos esto? Si la idea es usar la metrica.


    # comando = ["/usr/lib/linux-tools/5.15.0-134-generic/perf", "stat", "-r", "5", "--quiet", "./graphics/fire_animation_data" , "./data/2000_8", "0"]
    # res = subprocess.run(comando, stderr=subprocess.PIPE, text=True)
    # print("Tiempo medio de celdas por segundo: ", mediaCeldasPorSeg(res))

    with open('resultados.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(["Compilador", "Optimizador", "Landscape", "Tiempo promedio de ejecuciones", "Celdas quemadas por segundo"])
        for landscape in landscapes:
            for comp in compiladores:

                for opt in optimizadores:
                    os.system("make clean && make CXX=" + comp + " CXXOPT=" + opt)
                    comandoPerf = ["/usr/lib/linux-tools/5.15.0-134-generic/perf", "stat", "-r", str(landscape[1]), "./graphics/fire_animation_data" , landscape[0], "0"]
                    resPerf = subprocess.run(comandoPerf, stderr=subprocess.PIPE, text=True)

                    comandoCells = ["/usr/lib/linux-tools/5.15.0-134-generic/perf", "stat", "-r", str(landscape[1]), "--quiet", "./graphics/fire_animation_data" , landscape[0], "0"]
                    resCells = subprocess.run(comandoCells, stderr=subprocess.PIPE, text=True)

                    # exportar datos a csv:
                    writer.writerow([comp, opt, landscape[0], filtrarPerf(resPerf), maxCeldasPorSeg(resCells)])

if __name__ == "__main__":
    main()