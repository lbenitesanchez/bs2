.PHONY: install check quick smoke full bmd aggregate figures benchmark freeze clean

install:
	Rscript scripts/00_install.R

check:
	Rscript scripts/00_check_environment.R

quick:
	CORES=$${CORES:-2} BACKEND=$${BACKEND:-multisession} Rscript scripts/00_run_all_quick.R

smoke:
	Rscript tests/smoke_test.R

full:
	./scripts/00_run_all_full_ubuntu.sh

bmd:
	Rscript scripts/08_bmd_analysis.R --B=499 --B-sign=99999 --out-dir=results/tables

aggregate:
	Rscript scripts/04_aggregate_results.R --raw-est=results/raw_estimation --raw-tests=results/raw_tests --raw-boot=results/raw_bootstrap --out-dir=results/tables

figures:
	Rscript scripts/07_make_figures_publication.R --fig-dir=figures --formats=eps,pdf
	Rscript scripts/09_bmd_bs2_lognormal_contours.R --fig-dir=figures --formats=pdf,png

benchmark:
	Rscript scripts/10_algorithmic_benchmark.R --R=200 --cores=$${CORES:-2} --backend=$${BACKEND:-multisession}

freeze:
	Rscript scripts/00_freeze_renv.R

clean:
	rm -rf results/raw_* results/tables_debug results/tables_quick results/tables_power_debug results/tables_gof_debug logs/*.log figures/*
