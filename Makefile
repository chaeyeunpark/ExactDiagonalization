.PHONY: format test

help:
	@echo "Choose targets."

format:
	@./bin/format_cpp.sh

test:
	rm -rf ./BuildTests
	cmake -BBuildTests -DBUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Debug .
	cmake --build ./BuildTests --target runner
	cmake --build ./BuildTests --target test

check-tidy:
	rm -rf ./BuildTidy
	cmake -BBuildTidy -DBUILD_TESTS=ON -DENABLE_CLANG_TIDY=ON -DCMAKE_BUILD_TYPE=Debug .
	cmake --build ./BuildTidy --target runner
	cmake --build ./BuildTidy --target test
