# Developer entry points: make test / make lint / make all (see docs/CONTRIBUTING.md).
# nf-test needs a writable java.io.tmpdir and a real user.home (users without a passwd entry have none).

NFT_JAVA_OPTS ?= -Djava.io.tmpdir=$(or $(TMPDIR),/tmp) -Duser.home=$(HOME)
NFT = JAVA_TOOL_OPTIONS="$(NFT_JAVA_OPTS)" nf-test test
LINT_PATHS = main.nf nextflow.config conf modules subworkflows workflows tests
# nextflow lint's default excludes, plus FragPipe work-dir fixtures (FragPipe's own *.config files)
LINT_EXCLUDES = $(foreach p,.git .lineage .nf-test .nextflow work fragpipe_plans ptmshepherd_parity,-exclude $(p))

.PHONY: all test test-integration test-parity test-python lint

all: lint test test-python

## Stub and parity tests, the command of CI's blocking nf-test job (.github/workflows/nf-test.yml): every
## stub- or parity-tagged test in one pass under the docker profile, as nf-core/modules runs them
test:
	$(NFT) --profile=+docker --tag stub,parity --ci

## Integration tests: need docker and the licensed FragPipe image passed at runtime
test-integration:
	$(NFT) --tag integration

## FragPipe parity tests
test-parity:
	$(NFT) --tag parity

## Python unit tests of bin/, tests/tools/ and the module scripts (deps: tests/python/requirements.txt)
PYTHON_TEST_PATHS = tests/python
test-python:
	python3 -m pytest $(PYTHON_TEST_PATHS)

## pre-commit hooks and nextflow lint (0 errors required)
lint:
	pre-commit run --all-files
	nextflow lint -o concise $(LINT_EXCLUDES) $(LINT_PATHS)
