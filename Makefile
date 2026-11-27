# Unified local commands. `make verify` is the mandatory local quality gate.
SHELL := /bin/bash
LLVM_BIN ?= $(shell if [ -d /opt/homebrew/opt/llvm/bin ]; then echo /opt/homebrew/opt/llvm/bin; elif [ -d /usr/local/opt/llvm/bin ]; then echo /usr/local/opt/llvm/bin; fi)
CLANG_FORMAT ?= $(if $(LLVM_BIN),$(LLVM_BIN)/clang-format,clang-format)
CLANG_TIDY ?= $(if $(LLVM_BIN),$(LLVM_BIN)/clang-tidy,clang-tidy)
# Homebrew clang-tidy does not know the macOS SDK; pass it explicitly on Darwin.
TIDY_EXTRA ?= $(shell if [ "$$(uname -s)" = Darwin ]; then echo "--extra-arg=-isysroot$$(xcrun --show-sdk-path)"; fi)
JOBS ?= $(shell sysctl -n hw.ncpu 2>/dev/null || nproc)
PRESET ?= debug
CXX_SOURCES = $(shell git ls-files --cached --others --exclude-standard -- '*.cpp' '*.hpp' '*.h')
BENCH_TIMEOUT ?= 3600

.PHONY: public-safety fetch-paper bootstrap configure build test test-unit test-integration test-differential test-property \
        test-regression test-harness test-python format format-check lint tidy sanitizer coverage \
        verify benchmark-smoke benchmark benchmark-full report clean

public-safety:
	uv run python scripts/check_public_safety.py

fetch-paper:
	./scripts/fetch_paper.sh

bootstrap:
	./scripts/bootstrap.sh

configure:
	cmake --preset $(PRESET)

build: configure
	cmake --build --preset $(PRESET) -j $(JOBS)

test: build
	ctest --preset $(PRESET)
	$(MAKE) test-harness

test-unit: build
	ctest --preset $(PRESET) -L unit
test-integration: build
	ctest --preset $(PRESET) -L integration
test-differential: build
	ctest --preset $(PRESET) -L differential
test-property: build
	ctest --preset $(PRESET) -L property
test-regression: build
	ctest --preset $(PRESET) -L regression

test-harness:
	uv run pytest
	uv run python scripts/tdd.py check

format:
	@[ -z "$(CXX_SOURCES)" ] || $(CLANG_FORMAT) -i $(CXX_SOURCES)
	uv run ruff format .
	uv run ruff check --fix .

format-check:
	@[ -z "$(CXX_SOURCES)" ] || $(CLANG_FORMAT) --dry-run --Werror $(CXX_SOURCES)
	uv run ruff format --check .

lint: format-check
	uv run ruff check .
	$(MAKE) tidy

tidy: configure
	@files="$$(git ls-files --cached --others --exclude-standard -- 'src/*.cpp' 'tools/*.cpp')"; \
	if [ -n "$$files" ]; then $(CLANG_TIDY) -p build/$(PRESET) --quiet $(TIDY_EXTRA) $$files; fi

sanitizer:
	cmake --preset sanitizer && cmake --build --preset sanitizer -j $(JOBS) && ctest --preset sanitizer

coverage:
	./scripts/coverage.sh

verify:
	uv run python scripts/check_public_safety.py
	$(MAKE) format-check
	uv run ruff check .
	$(MAKE) PRESET=ci build
	ctest --preset ci
	$(MAKE) tidy PRESET=ci
	uv run pytest
	uv run python scripts/tdd.py check

benchmark-smoke:
	./scripts/run_benchmarks.sh smoke

benchmark:
	./scripts/run_benchmarks.sh standard

benchmark-full:
	BENCH_TIMEOUT=$(BENCH_TIMEOUT) ./scripts/run_benchmarks.sh full

report:
	uv run python scripts/analyze_results.py

clean:
	rm -rf build
