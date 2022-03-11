.PHONY: format test

help:
	@echo "Choose targets."

format:
	@./bin/format_cpp.sh

test:
	rm -rf ./BuildTests
	cmake -BBuildTests -DBUILD_TESTS=ON .
	cmake --build ./BuildTests --target runner
	cmake --build ./BuildTests --target test
