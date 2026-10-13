# Developer entry points required by .claude/rules/general.md (make test / make lint / make all).
# nf-test needs a writable java.io.tmpdir and a real user.home (sandboxed uids have no passwd entry).

NFT_JAVA_OPTS ?= -Djava.io.tmpdir=$(or $(TMPDIR),/tmp) -Duser.home=$(HOME)
NFT = JAVA_TOOL_OPTIONS="$(NFT_JAVA_OPTS)" nf-test test
LINT_PATHS = main.nf nextflow.config conf modules subworkflows workflows tests

.PHONY: all test test-integration test-parity lint

all: lint test

## Stub tests, the command of CI's blocking nf-test job (.github/actions/nf-test/action.yml): every
## stub-tagged test in one pass under the docker profile, as nf-core/modules runs them
test:
	$(NFT) --profile=+docker --tag stub --ci

## Integration tests: need docker and the licensed FragPipe image passed at runtime
test-integration:
	$(NFT) --tag integration

## FragPipe parity tests
test-parity:
	$(NFT) --tag parity

## pre-commit hooks and nextflow lint (0 errors required)
lint:
	pre-commit run --all-files
	nextflow lint -o concise $(LINT_PATHS)
