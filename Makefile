help:
	@echo 'Make for some simple commands        '
	@echo '                                     '
	@echo ' Usage:                              '
	@echo '     make lint    flake8 the codebase'
	@echo '     make test    run unit tests     '

lint:
	flake8 ./blotter

test:
	python -m unittest -v blotter.tests.test_holdings
	python -m unittest -v blotter.tests.test_event
	python -m unittest -v blotter.tests.test_marketdata
	python -m unittest -v blotter.tests.test_blotter
asv:
	asv run -v
