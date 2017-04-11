help:
	@echo 'Make for some simple commands            '
	@echo '                                         '
	@echo ' Usage:                                  '
	@echo '     make lint    flake8 the codebase    '
	@echo '     make test    run unit tests         '
	@echo '     make asv     run performance tests  '

lint:
	flake8 ./blotter

test:
	pytest blotter/tests -v --cov=blotter --cov-report term-missing

asv:
	asv run -v
