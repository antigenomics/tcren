# Build the Shuffled2026 BN appendix: train/eval (build.py) -> gnuplot figures -> latexmk.
all: shuffle_bn.pdf

.PHONY: data figs clean
data:
	python build.py

figs: data
	gnuplot roc.gp
	gnuplot pr.gp
	gnuplot marginals.gp

shuffle_bn.pdf: shuffle_bn.tex figs
	latexmk -pdf -interaction=nonstopmode shuffle_bn.tex

clean:
	latexmk -C 2>/dev/null || true
	rm -f figures/*.dat figures/*.pdf figures/*.dot metrics.tex marginal_labels.tex base.txt
