import os
import matplotlib.pyplot as plt

# # inputFile = '/data/user/borzari/cmssw/pixeltrack-standalone/output.txt'
eventsStr = ''
# inputFile = 'output_' + eventsStr + '.txt'
inputFile = 'output.txt'

with open(inputFile, "r") as f:
    linhas = f.readlines()

numLinhas = 22

if len(linhas) != numLinhas:
    raise ValueError(f"O arquivo deve conter exatamente {numLinhas} linhas")

def parse_linha(linha):
    return [float(valor.strip()) for valor in linha.strip().split(",") if valor.strip()]

y1 = parse_linha(linhas[0])
x1 = parse_linha(linhas[1])
y2 = parse_linha(linhas[2])
x2 = parse_linha(linhas[3])
y3 = parse_linha(linhas[4])
x3 = parse_linha(linhas[5])
y4 = parse_linha(linhas[6])
x4 = parse_linha(linhas[7])
y5 = parse_linha(linhas[8])
x5 = parse_linha(linhas[9])
y6 = parse_linha(linhas[10])
x6 = parse_linha(linhas[11])
y7 = parse_linha(linhas[12])
x7 = parse_linha(linhas[13])
y8 = parse_linha(linhas[14])
x8 = parse_linha(linhas[15])
y9 = parse_linha(linhas[16])
x9 = parse_linha(linhas[17])
y10 = parse_linha(linhas[18])
x10 = parse_linha(linhas[19])
y11 = parse_linha(linhas[20])
x11 = parse_linha(linhas[21])

plt.figure(figsize=(6, 6))
plt.scatter(x1, y1, label="Pt", color='red')
plt.xlabel("Pt")
plt.xscale("log")
plt.ylabel("Eff")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('testEffPt' + eventsStr + '.pdf')


plt.figure(figsize=(6, 6))
plt.scatter(x2, y2, label="Eta", color='green')
plt.xlabel("Eta")
plt.ylabel("Eff")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('testEffEta' + eventsStr + '.pdf')


plt.figure(figsize=(6, 6))
plt.scatter(x3, y3, label="Phi", color='blue')
plt.xlabel("Phi")
plt.ylabel("Eff")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('testEffPhi' + eventsStr + '.pdf')

plt.figure(figsize=(6, 6))
plt.scatter(x4, y4, label="Pt", color='red')
plt.xlabel("Pt")
plt.xscale("log")
plt.ylabel("Fake")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('testFakePt' + eventsStr + '.pdf')


plt.figure(figsize=(6, 6))
plt.scatter(x5, y5, label="Eta", color='green')
plt.xlabel("Eta")
plt.ylabel("Fake")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('testFakeEta' + eventsStr + '.pdf')


plt.figure(figsize=(6, 6))
plt.scatter(x6, y6, label="Phi", color='blue')
plt.xlabel("Phi")
plt.ylabel("Fake")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('testFakePhi' + eventsStr + '.pdf')

plt.figure(figsize=(6, 6))
plt.scatter(x7, y7, label="Pt", color='blue')
plt.xlabel("Eta")
plt.ylabel("Res")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('testResPt' + eventsStr + '.pdf')

plt.figure(figsize=(6, 6))
plt.scatter(x8, y8, label="Eta", color='blue')
plt.xlabel("Eta")
plt.ylabel("Res")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('testResEta' + eventsStr + '.pdf')

plt.figure(figsize=(6, 6))
plt.scatter(x9, y9, label="Phi", color='blue')
plt.xlabel("Eta")
plt.ylabel("Res")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('testResPhi' + eventsStr + '.pdf')

plt.figure(figsize=(6, 6))
plt.scatter(x10, y10, label="D0", color='blue')
plt.xlabel("Eta")
plt.ylabel("Res")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('testResD0' + eventsStr + '.pdf')

plt.figure(figsize=(6, 6))
plt.scatter(x11, y11, label="Dz", color='blue')
plt.xlabel("Eta")
plt.ylabel("Res")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('testResDz' + eventsStr + '.pdf')

