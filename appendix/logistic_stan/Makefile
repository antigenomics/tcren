# Build the distribution-aware Bayesian logistic appendix.
# build.py needs PyMC+arviz+sklearn: pass PY=<python with those>, e.g. the manuscript results/eda/.venv.
PY ?= python

all: logistic_stan.pdf

.PHONY: data figs clean
data:
	$(PY) build.py

figs: data
	gnuplot roc.gp
	gnuplot pr.gp
	gnuplot forest.gp

logistic_stan.pdf: logistic_stan.tex figs
	latexmk -pdf -interaction=nonstopmode logistic_stan.tex

clean:
	latexmk -C 2>/dev/null || true
	rm -f figures/*.dat figures/*.pdf base.txt metrics.tex forest_labels.tex
