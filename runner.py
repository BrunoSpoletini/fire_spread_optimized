import subprocess
import os

compiladores = ["g++", "clang"]
optimizadores = ["-O0"]
#optimizadores = ["-O0", "-O1", "-O2", "-O3", "-Ofast"]
landscapes = [("./data/2000_8", 16)]
#landscapes = [("./data/2000_8", 16), ("./data/1999_27j_S", 16), ("./data/2015_50", 4)]

def main():

    resultados = dict()

    for landscape in landscapes:
        for comp in compiladores:
            resultados[comp] = []
            for opt in optimizadores:
                os.system("make clean && make CXX=" + comp + " CXXOPT=" + opt)

                comando = ["perf", "stat", "-r", str(landscape[1]), "./grafics/fire_animation_data" , landscape[0], "0"]
                resultado = subprocess.run(comando, stderr=subprocess.PIPE, text=True)

                print(resultado)
    #             for linea in resultado.stderr.split("\n"):
    #                 if "time elapsed" in linea:
    #                     resultados[comp].insert(re.findall(r"[\d.]+", linea)[0])

    #             print(resultados)

    # print(resultados)


if __name__ == "__main__":
    main()