import re

input_bibfile = "/home/marc/LaTex/cGAN-islands/egbibsample.bib"
input_tex = "/home/marc/LaTex/cGAN-islands/main.tex"

refs = {}

with open(input_tex, "r") as f:
    _lines = f.readlines()
    lines = []
    for l in _lines:
        lines += l.split("}")
    expr = re.compile(r"\\cite{(.*)$")
    for line in lines:
        citations = expr.findall(line)
        for _cites in citations:
            cites = _cites.split(",")
            for cite in cites:
                refs[cite] = ""

with open(input_bibfile, "r") as f:
    lines = f.readlines()
    currentRead = None
    for line in lines:
        if currentRead is None:
            for author in refs:
                if author in line:
                    currentRead = author
                    break
        if currentRead:
            refs[currentRead] += line

        if line.strip() == "}":
            currentRead = None

for ref in refs:
    print(refs[ref])



